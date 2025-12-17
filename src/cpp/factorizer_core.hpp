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
 * @note Tie-breaking: Both forward and RC candidates are tracked independently during the tree walk,
 *       then their TRUE LCPs are computed and compared. When true lengths are equal, forward is preferred.
 *       Among candidates of the same type, the one with the earliest position (smallest start for forward, smallest end for RC) wins.
 */
template<class Sink>
size_t nolzss_multiple_dna_w_rc(const std::string& S, Sink&& sink, size_t start_pos) {
    const size_t N = (S.size() / 2) - 1;
    if (N == 0) return 0;
    
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
        // At factor start i (0-based in T), walk up ancestors and track best FWD and RC candidates independently
        
        // Track best FWD candidate
        bool   have_fwd       = false;
        size_t best_fwd_start = 0;  // start in T (for FWD)
        size_t best_fwd_depth = 0;  // tree depth of best FWD candidate
        
        // Track best RC candidate
        bool   have_rc        = false;
        size_t best_rc_end    = 0;  // end in T (for RC)
        size_t best_rc_posS   = 0;  // pos in S where RC candidate suffix starts (for LCP)
        size_t best_rc_depth  = 0;  // tree depth of best RC candidate

        // Walk from leaf to root via level_anc
        for (size_t step = 1; step <= lambda_node_depth; ++step) {
            auto v = cst.bp_support.level_anc(lambda, lambda_node_depth - step);
            size_t ell = cst.depth(v);
            if (ell == 0) break; // reached root

            auto lb = cst.lb(v), rb = cst.rb(v);

            // Forward candidate (min start in T within v's interval)
            size_t kF = rmqF(lb, rb);
            uint64_t jF = fwd_starts[kF];
            bool okF = (jF != INF) && (jF + ell - 1 < i); // non-overlap: endF <= i-1

            // RC candidate (min END in T within v's interval; monotone with depth)
            size_t kR = rmqRcEnd(lb, rb);
            uint64_t endRC = rc_ends[kR];
            bool okR = (endRC != INF) && (endRC < i); // endRC <= i-1

            if (!okF && !okR) {
                // deeper nodes can only increase jF and the minimal RC end
                // -> non-overlap won't become true again for either; stop
                break;
            }

            // Update best FWD candidate independently
            if (okF) {
                if (ell > best_fwd_depth ||
                    (ell == best_fwd_depth && (jF + ell - 1) < (best_fwd_start + best_fwd_depth - 1))) {
                    best_fwd_depth = ell;
                    best_fwd_start = jF;
                    have_fwd = true;
                }
            }
            
            // Update best RC candidate independently
            if (okR) {
                size_t posS_R = cst.csa[kR]; // suffix position in S for LCP
                if (ell > best_rc_depth ||
                    (ell == best_rc_depth && (endRC < best_rc_end))) {
                    best_rc_depth = ell;
                    best_rc_end   = endRC;
                    best_rc_posS  = posS_R;
                    have_rc = true;
                }
            }
        }

        size_t emit_len = 1;
        uint64_t emit_ref = i; // default for literal
        
        if (!have_fwd && !have_rc) {
            // No previous occurrence (FWD nor RC) — literal of length 1
            Factor f{static_cast<uint64_t>(i), static_cast<uint64_t>(emit_len), static_cast<uint64_t>(emit_ref)};
            sink(f);
            ++factors;

            // Advance
            lambda = next_leaf(cst, lambda, emit_len);
            lambda_node_depth = cst.node_depth(lambda);
            i = cst.sn(lambda);
            continue;
        }

        // Compute true lengths for both candidates
        size_t fwd_true_len = 0;
        size_t rc_true_len = 0;
        
        if (have_fwd) {
            size_t cap = i - best_fwd_start; // non-overlap cap
            size_t L = lcp(cst, i, best_fwd_start);
            fwd_true_len = std::min(L, cap);
        }
        
        if (have_rc) {
            rc_true_len = lcp(cst, i, best_rc_posS);
        }
        
        // Choose based on TRUE length, with FWD preference at equal length.
        // A literal (self-reference, length 1) is always available as a fallback.
        // Only use RC if it provides strictly longer match than both FWD and literal.
        bool use_fwd = false;
        bool use_literal = false;
        
        if (have_fwd && fwd_true_len >= 1) {
            // Have a real forward match
            if (have_rc && rc_true_len > fwd_true_len) {
                use_fwd = false;  // RC is strictly longer
            } else {
                use_fwd = true;   // FWD wins (including ties)
            }
        } else {
            // No real forward match; compare RC against literal (length 1)
            if (have_rc && rc_true_len > 1) {
                use_fwd = false;  // RC is strictly longer than literal
            } else {
                use_literal = true;  // Use literal (self-reference)
            }
        }
        
        if (use_literal) {
            // Emit literal
            emit_len = 1;
            emit_ref = static_cast<uint64_t>(i);
        } else if (use_fwd) {
            emit_len = fwd_true_len;
            emit_ref = static_cast<uint64_t>(best_fwd_start);
        } else {
            emit_len = rc_true_len;
            size_t start_pos_val = best_rc_end - emit_len + 1;  // start = end - len + 1
            emit_ref = RC_MASK | static_cast<uint64_t>(start_pos_val); // start-anchored + RC flag
        }
        
        // Safety: ensure progress
        if (emit_len <= 0) {
            throw std::runtime_error("emit_len must be positive to ensure factorization progress");
        }

        Factor f{static_cast<uint64_t>(i), static_cast<uint64_t>(emit_len), emit_ref};
        sink(f);
        ++factors;

        // Advance to next phrase start
        lambda = next_leaf(cst, lambda, emit_len);
        lambda_node_depth = cst.node_depth(lambda);
        i = cst.sn(lambda);
    }

    return factors;
}

} // namespace detail
} // namespace noLZSS
