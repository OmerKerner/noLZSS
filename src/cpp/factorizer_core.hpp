/**
 * @file factorizer_core.hpp
 * @brief Template implementations for core noLZSS factorization algorithms
 *
 * This header contains the template implementations of the sink-based factorization
 * algorithms. These are in the `detail` namespace to indicate they are internal
 * implementation details. Public API wrappers are in factorizer.hpp.
 *
 * By keeping templates in a header, they can be instantiated in any compilation unit
 * without duplication, solving the "declared using local type" lambda issue.
 */

#pragma once

#include "factorizer.hpp"
#include "factorizer_helpers.hpp"
#include <sdsl/suffix_trees.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <limits>

namespace noLZSS {
namespace detail {

// Forward declarations for template functions (default arguments specified here only)
template<class Sink>
size_t nolzss(cst_t& cst, Sink&& sink, size_t start_pos = 0);

template<class Sink>
size_t nolzss_dna_w_rc(const std::string& T, Sink&& sink);

template<class Sink>
size_t nolzss_multiple_dna_w_rc(const std::string& S, Sink&& sink, size_t start_pos = 0);

/**
 * @brief Core noLZSS factorization algorithm implementation.
 *
 * Implements the non-overlapping Lempel-Ziv-Storer-Szymanski factorization
 * using a compressed suffix tree. The algorithm finds the longest previous
 * factor for each position in the text and emits factors through a sink.
 *
 * @tparam Sink Callable type that accepts Factor objects (e.g., lambda, function)
 * @param cst The compressed suffix tree built from the input text
 * @param sink Callable that receives each computed factor
 * @param start_pos Starting position for factorization (default: 0)
 * @return Number of factors emitted
 *
 * @note This is the core algorithm that all public functions use
 * @note The sink pattern allows for memory-efficient processing
 * @note All factors are emitted, including the last one
 */
template<class Sink>
size_t nolzss(cst_t& cst, Sink&& sink, size_t start_pos) {
    sdsl::rmq_succinct_sct<> rmq(&cst.csa);
    const size_t str_len = cst.size() - 1; // the length of the string is the size of the CST minus the sentinel

    auto lambda = cst.select_leaf(cst.csa.isa[start_pos] + 1);
    size_t lambda_node_depth = cst.node_depth(lambda);
    size_t lambda_sufnum = start_pos;

    cst_t::node_type v;
    size_t v_min_leaf_sufnum = 0;
    size_t u_min_leaf_sufnum = 0;

    size_t count = 0;

    while (lambda_sufnum < str_len) {
        // Compute current factor
        size_t d = 1;
        size_t l = 1;
        while (true) {
            v = cst.bp_support.level_anc(lambda, lambda_node_depth - d);
            v_min_leaf_sufnum = cst.csa[rmq(cst.lb(v), cst.rb(v))];
            l = cst.depth(v);

            if (v_min_leaf_sufnum + l - 1 < lambda_sufnum) {
                u_min_leaf_sufnum = v_min_leaf_sufnum;
                ++d; continue;
            }
            auto u = cst.parent(v);
            auto u_depth = cst.depth(u);

            if (v_min_leaf_sufnum == lambda_sufnum) {
                if (u == cst.root()) {
                    l = 1;
                    Factor factor{static_cast<uint64_t>(lambda_sufnum), static_cast<uint64_t>(l), static_cast<uint64_t>(lambda_sufnum)};
                    sink(factor);
                    break;
                }
                else {
                    l = u_depth;
                    Factor factor{static_cast<uint64_t>(lambda_sufnum), static_cast<uint64_t>(l), static_cast<uint64_t>(u_min_leaf_sufnum)};
                    sink(factor);
                    break;
                }
            }
            l = std::min(lcp(cst, lambda_sufnum, v_min_leaf_sufnum),
                         (lambda_sufnum - v_min_leaf_sufnum));
            if (l <= u_depth) {
                l = u_depth;
                Factor factor{static_cast<uint64_t>(lambda_sufnum), static_cast<uint64_t>(l), static_cast<uint64_t>(u_min_leaf_sufnum)};
                sink(factor);
                break;
            }
            else {
                Factor factor{static_cast<uint64_t>(lambda_sufnum), static_cast<uint64_t>(l), static_cast<uint64_t>(v_min_leaf_sufnum)};
                sink(factor);
                break;
            }
        }

        ++count;
        // Advance to next position
        lambda = next_leaf(cst, lambda, l);
        lambda_node_depth = cst.node_depth(lambda);
        lambda_sufnum = cst.sn(lambda);
    }

    return count;
}

/**
 * @brief Core noLZSS factorization algorithm implementation with reverse complement awareness for DNA.
 *
 * Implements the non-overlapping Lempel-Ziv-Storer-Szymanski factorization
 * using a compressed suffix tree, extended to handle DNA sequences with reverse complement matches.
 * The algorithm constructs a combined string S = T '$' rc(T) '#' where rc(T) is the reverse complement,
 * builds a suffix tree over S, and finds the longest previous factor (either forward or reverse complement)
 * for each position in the original text T, emitting factors through a sink.
 *
 * @tparam Sink Callable type that accepts Factor objects (e.g., lambda, function)
 * @param T Input DNA text string
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This is the core algorithm for DNA-aware factorization that all DNA public functions use
 * @note The sink pattern allows for memory-efficient processing
 * @note All factors are emitted, including the last one
 * @note Reverse complement matches are encoded with the RC_MASK in the ref field
 */
template<class Sink>
size_t nolzss_dna_w_rc(const std::string& T, Sink&& sink) {
    const size_t n = T.size();
    if (n == 0) return 0;

    // Use prepare_multiple_dna_sequences_w_rc with a single sequence
    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc({T});
    std::string S = prep_result.prepared_string;

    // Use the multiple sequence algorithm
    return nolzss_multiple_dna_w_rc(S, std::forward<Sink>(sink));
}

/**
 * @brief Outcome of the longest-previous-match search at one position of the RC-aware algorithm.
 */
struct RcMatchResult {
    size_t   len = 1;          /**< Length of the factor to emit (1 for a literal) */
    uint64_t ref = 0;          /**< Reference to emit (RC_MASK set for reverse-complement matches; == position for a literal) */
    bool     is_literal = true;/**< True when no usable previous occurrence exists */
};

/**
 * @brief Finds the longest non-overlapping previous occurrence (forward or reverse complement) of the suffix at position i.
 *
 * Shared by nolzss_multiple_dna_w_rc() and ParallelFactorizer::factorize_dna_w_rc_thread() so that both
 * paths select factors identically.
 *
 * The walk goes from the shallowest non-root ancestor of the current leaf toward the leaf. At each node v
 * (string depth ell) it inspects
 *   - the earliest forward start jF in v's SA interval (via rmqF): usable if jF + ell - 1 < i, and
 *   - the earliest reverse-complement end endRC in v's SA interval (via rmqRcEnd): usable if endRC < i.
 * Both conditions are monotone along the path (min start / min end never decrease with depth), so the walk
 * stops as soon as neither holds.
 *
 * Forward "partial-edge" candidate: the first node c whose earliest start jF overlaps i (jF + ell - 1 >= i,
 * jF < i) still provides a non-overlapping match of length i - jF, because every leaf of c shares at least
 * ell > i - jF characters with the suffix at i. This is the case the base algorithm handles on the edge
 * between the last fitting ancestor and c, and it can be longer than the deepest fully fitting node.
 * Reverse-complement matches need no such treatment: their non-overlap test (endRC < i) does not depend on
 * the match length, so the deepest usable node already yields the longest RC match.
 *
 * Tie-breaking: at equal length a forward match wins over an RC match, and an RC match must be strictly
 * longer than a literal (length >= 2) to be used. Among forward candidates of equal length the earliest
 * start wins; among RC candidates the earliest end wins.
 *
 * @param cst Compressed suffix tree over the prepared string S = T sentinel rc(T) sentinel
 * @param rmqF RMQ over fwd_starts
 * @param rmqRcEnd RMQ over rc_ends
 * @param fwd_starts Per-SA-rank forward start in T (INF for suffixes outside T)
 * @param rc_ends Per-SA-rank mapped end in T of RC suffixes (INF for suffixes outside the RC part)
 * @param INF Marker used in fwd_starts / rc_ends
 * @param lambda Leaf of the suffix starting at S position i
 * @param lambda_node_depth Node depth of lambda
 * @param i Current position (0-based in T)
 * @return The factor length and reference to emit at position i
 */
inline RcMatchResult find_longest_rc_match(const cst_t& cst,
                                           const sdsl::rmq_succinct_sct<>& rmqF,
                                           const sdsl::rmq_succinct_sct<>& rmqRcEnd,
                                           const sdsl::int_vector<64>& fwd_starts,
                                           const sdsl::int_vector<64>& rc_ends,
                                           uint64_t INF,
                                           cst_t::node_type lambda,
                                           size_t lambda_node_depth,
                                           size_t i) {
    // Best node-aligned forward candidate (deepest node whose earliest start fits before i)
    bool   have_fwd       = false;
    size_t best_fwd_start = 0;
    size_t best_fwd_depth = 0;

    // Forward partial-edge candidate (first node whose earliest start overlaps i)
    bool   have_partial   = false;
    size_t partial_start  = 0;
    size_t partial_len    = 0;

    // Best RC candidate (deepest node whose earliest end lies before i)
    bool   have_rc        = false;
    size_t best_rc_end    = 0;
    size_t best_rc_posS   = 0;
    size_t best_rc_depth  = 0;

    for (size_t step = 1; step <= lambda_node_depth; ++step) {
        auto v = cst.bp_support.level_anc(lambda, lambda_node_depth - step);
        size_t ell = cst.depth(v);
        if (ell == 0) break; // root

        auto lb = cst.lb(v), rb = cst.rb(v);

        size_t   kF  = rmqF(lb, rb);
        uint64_t jF  = fwd_starts[kF];
        bool     okF = (jF != INF) && (jF + ell - 1 < i);

        size_t   kR    = rmqRcEnd(lb, rb);
        uint64_t endRC = rc_ends[kR];
        bool     okR   = (endRC != INF) && (endRC < i);

        if (!okF && !have_partial && jF != INF && jF < i) {
            // Earliest occurrence overlaps i: the prefix T[jF..i) still matches T[i..) for i - jF characters.
            have_partial  = true;
            partial_start = jF;
            partial_len   = i - jF;
        }

        if (!okF && !okR) {
            // Deeper nodes only increase the minimal forward start and minimal RC end.
            break;
        }

        if (okF && ell > best_fwd_depth) {
            have_fwd       = true;
            best_fwd_start = jF;
            best_fwd_depth = ell;
        }

        if (okR && ell > best_rc_depth) {
            have_rc       = true;
            best_rc_end   = endRC;
            best_rc_posS  = cst.csa[kR];
            best_rc_depth = ell;
        }
    }

    // True forward length: node-aligned candidate (possibly extending along the edge below it), then the
    // partial-edge candidate if it is strictly longer (at equal length the node-aligned one has the earlier start).
    size_t fwd_len   = 0;
    size_t fwd_start = 0;
    if (have_fwd) {
        size_t cap = i - best_fwd_start;
        size_t L   = lcp(cst, i, best_fwd_start);
        fwd_len    = std::min(L, cap);
        fwd_start  = best_fwd_start;
    }
    if (have_partial && partial_len > fwd_len) {
        fwd_len   = partial_len;
        fwd_start = partial_start;
    }

    size_t rc_len = 0;
    if (have_rc) {
        rc_len = lcp(cst, i, best_rc_posS);
    }

    RcMatchResult out;
    if (fwd_len >= 1 && !(have_rc && rc_len > fwd_len)) {
        out.len = fwd_len;
        out.ref = static_cast<uint64_t>(fwd_start);
        out.is_literal = false;
    } else if (have_rc && rc_len > 1 && rc_len > fwd_len) {
        out.len = rc_len;
        size_t rc_start = best_rc_end - rc_len + 1;
        out.ref = RC_MASK | static_cast<uint64_t>(rc_start);
        out.is_literal = false;
    } else {
        out.len = 1;
        out.ref = static_cast<uint64_t>(i);
        out.is_literal = true;
    }

    if (out.len == 0) {
        throw std::runtime_error("factor length must be positive to ensure factorization progress");
    }
    return out;
}

/**
 * @brief Core noLZSS factorization algorithm implementation with reverse complement awareness for multiple DNA sequences.
 *
 * Implements the non-overlapping Lempel-Ziv-Storer-Szymanski factorization
 * using a compressed suffix tree, extended to handle multiple DNA sequences with reverse complement matches.
 * The algorithm takes a concatenated string S of multiple sequences with sentinels and their reverse complements,
 * builds a suffix tree over S, and finds the longest previous factor (either forward or reverse complement)
 * for each position in the original sequences, emitting factors through a sink.
 *
 * @tparam Sink Callable type that accepts Factor objects (e.g., lambda, function)
 * @param S Input concatenated DNA text string with sentinels and reverse complements
 * @param sink Callable that receives each computed factor
 * @param start_pos Starting position for factorization (default: 0)
 * @return Number of factors emitted
 *
 * @note This is the core algorithm for multiple DNA sequences factorization that all multiple DNA public functions use
 * @note The sink pattern allows for memory-efficient processing
 * @note All factors are emitted, including the last one
 * @note Reverse complement matches are encoded with the RC_MASK in the ref field
 * @note start_pos allows factorization to begin from a specific position, useful for reference+target factorization
 * @note Candidate selection (walk, partial-edge forward candidate, tie-breaking) lives in find_longest_rc_match()
 *       and is shared with the parallel implementation.
 */
template<class Sink>
size_t nolzss_multiple_dna_w_rc(const std::string& S, Sink&& sink, size_t start_pos) {
    // Early return for empty input to avoid CST construction on empty string
    if (S.empty()) return 0;
    
    // Validate minimum string length to avoid underflow and SDSL errors
    // Format for w/ RC: T1!T2@...Tn$rt(Tn)%...rt(T1)&
    // Minimum valid input: "A$A&" (4 chars: sequence + sentinel + reverse_complement + sentinel)
    // - Position 0: Single nucleotide (e.g., 'A')
    // - Position 1: Sentinel separating original from RC (e.g., '$')
    // - Position 2: Reverse complement of the nucleotide (e.g., 'A' is revcomp of 'T')
    // - Position 3: Final sentinel (e.g., '&')
    if (S.size() < 4) {
        std::cerr << "Warning: Input string too short for factorization with reverse complement (size=" 
                  << S.size() << "). Returning 0 factors." << std::endl;
        return 0;
    }
    
    const size_t N = (S.size() / 2) - 1;
    if (N == 0) {
        std::cerr << "Warning: Computed N=0 from input size=" << S.size() 
                  << ". Returning 0 factors." << std::endl;
        return 0;
    }
    
    // Validate start_pos
    if (start_pos >= N) {
        throw std::invalid_argument("start_pos must be less than the original sequence length");
    }

    // Build CST over S
    cst_t cst; construct_im(cst, S, 1);

    // Build RMQ inputs aligned to SA: forward starts and RC ends (in T-coords)
    const uint64_t INF = std::numeric_limits<uint64_t>::max()/2ULL;
    sdsl::int_vector<64> fwd_starts(cst.csa.size(), INF);
    sdsl::int_vector<64> rc_ends   (cst.csa.size(), INF);

    const size_t T_end = N;           // end of original (exclusive)
    const size_t R_beg = N + 1;       // first char of rc (skip sentinel at position N)
    const size_t R_end = S.size() - 1;   // exclude final sentinel

    for (size_t k = 0; k < cst.csa.size(); ++k) {
        size_t posS = cst.csa[k];
        if (posS < T_end) {
            // suffix starts in T
            fwd_starts[k] = posS;     // 0-based start in T
        } else if (posS >= R_beg && posS < R_end) {
            // suffix starts in R
            size_t jR0   = posS - R_beg;         // 0-based start in R
            size_t endT0 = N - jR0 - 1;          // mapped end in T (0-based)
            rc_ends[k] = endT0;
        }
    }
    sdsl::rmq_succinct_sct<> rmqF(&fwd_starts);
    sdsl::rmq_succinct_sct<> rmqRcEnd(&rc_ends);

    // Initialize to the leaf of suffix starting at S position start_pos (i.e., T[start_pos])
    auto lambda = cst.select_leaf(cst.csa.isa[start_pos] + 1);
    size_t lambda_node_depth = cst.node_depth(lambda);
    size_t i = cst.sn(lambda); // suffix start in S, begins at start_pos

    size_t factors = 0;

    while (i < N) { // only factorize inside T
        RcMatchResult m = find_longest_rc_match(cst, rmqF, rmqRcEnd, fwd_starts, rc_ends, INF,
                                                lambda, lambda_node_depth, i);

        Factor f{static_cast<uint64_t>(i), static_cast<uint64_t>(m.len), m.ref};
        sink(f);
        ++factors;

        // Advance to next phrase start
        lambda = next_leaf(cst, lambda, m.len);
        lambda_node_depth = cst.node_depth(lambda);
        i = cst.sn(lambda);
    }

    return factors;
}

} // namespace detail
} // namespace noLZSS
