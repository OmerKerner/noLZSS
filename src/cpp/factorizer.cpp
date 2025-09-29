#include "factorizer.hpp"
#include <sdsl/suffix_trees.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <optional>
#include <limits>

namespace noLZSS {
using cst_t = sdsl::cst_sada<>;

/**
 * @brief Computes the longest common prefix between two suffixes.
 *
 * Uses the suffix tree's LCA (Lowest Common Ancestor) to efficiently
 * compute the length of the longest common prefix between suffixes
 * starting at positions i and j.
 *
 * @param cst The compressed suffix tree
 * @param i Starting position of first suffix
 * @param j Starting position of second suffix
 * @return Length of the longest common prefix
 */
static size_t lcp(cst_t& cst, size_t i, size_t j) {
    if (i == j) return cst.csa.size() - cst.csa[i];
    auto lca = cst.lca(cst.select_leaf(cst.csa.isa[i]+1), cst.select_leaf(cst.csa.isa[j]+1));
    return cst.depth(lca);
}

/**
 * @brief Advances a leaf node by a specified number of positions.
 *
 * Moves from the current leaf node forward by 'iterations' positions
 * in the suffix array order. This is used to advance the current
 * factorization position.
 *
 * @param cst The compressed suffix tree
 * @param lambda Current leaf node
 * @param iterations Number of positions to advance (default: 1)
 * @return The leaf node at the new position
 */
static cst_t::node_type next_leaf(cst_t& cst, cst_t::node_type lambda, size_t iterations = 1) {
    assert(cst.is_leaf(lambda));
    auto lambda_rank = cst.lb(lambda);
    for (size_t i = 0; i < iterations; i++) lambda_rank = cst.csa.psi[lambda_rank];
    return cst.select_leaf(lambda_rank + 1);
}

// ---------- genomic utilities ----------
char complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:
            // Handle invalid input, e.g., throw an exception or return a sentinel value
            throw std::invalid_argument("Invalid nucleotide: " + std::string(1, c));
    }
}
static std::string revcomp(std::string_view s) {
    std::string r; r.resize(s.size());
    for (size_t i = 0, n = s.size(); i < n; ++i) r[i] = complement(s[n-1-i]);
    return r;
}


// ---------- factorization utilities ----------
/**
 * @brief Prepares multiple DNA sequences for factorization with reverse complement and tracks sentinel positions.
 *
 * Takes multiple DNA sequences, concatenates them with unique sentinels, appends
 * their reverse complements with unique sentinels, and tracks sentinel positions.
 * The output format is compatible with nolzss_multiple_dna_w_rc(): S = T1!T2@T3$rt(T3)%rt(T2)^rt(T1)&
 *
 * @param sequences Vector of DNA sequence strings (should contain only A, C, T, G)
 * @return PreparedSequenceResult containing:
 *         - prepared_string: The formatted string with sequences and reverse complements
 *         - original_length: Length of the original sequences part (before reverse complements)
 *         - sentinel_positions: Positions of all sentinels in the prepared string
 *
 * @throws std::invalid_argument If too many sequences (>125) or invalid nucleotides found
 * @throws std::runtime_error If sequences contain invalid characters
 *
 * @note Sentinels avoid 0, A(65), C(67), G(71), T(84) - lowercase nucleotides are safe as sentinels
 * @note The function validates that all sequences contain only valid DNA nucleotides
 * @note Input sequences can be lowercase or uppercase, output is always uppercase
 */
PreparedSequenceResult prepare_multiple_dna_sequences_w_rc(const std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return {"", 0, {}};
    }
    
    // Check if we have too many sequences
    if (sequences.size() > 125) {
        throw std::invalid_argument("Too many sequences: maximum 125 sequences supported (due to sentinel character limitations)");
    }
    
    // Validate sequences contain only valid DNA nucleotides
    for (size_t i = 0; i < sequences.size(); ++i) {
        for (char c : sequences[i]) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && 
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                throw std::runtime_error("Invalid nucleotide '" + std::string(1, c) + 
                                       "' found in sequence " + std::to_string(i));
            }
        }
    }
    
    PreparedSequenceResult result;
    
    // Calculate total size for reservation
    size_t total_size = 0;
    for (const auto& seq : sequences) {
        total_size += seq.length() * 2; // Original + reverse complement
    }
    total_size += sequences.size() * 2; // Add space for sentinels
    result.prepared_string.reserve(total_size);
    
    // Generate sentinel characters avoiding 0 and uppercase DNA nucleotides A(65), C(67), G(71), T(84)
    auto get_sentinel = [](size_t index) -> char {
        char sentinel = 1;
        size_t count = 0;
        
        while (true) {
            // Check if current sentinel is valid (not 0 and not uppercase DNA nucleotides)
            if (sentinel != 0 && sentinel != 'A' && sentinel != 'C' && sentinel != 'G' && sentinel != 'T') {
                if (count == index) {
                    return sentinel;  // Found the sentinel for this index
                }
                count++;
            }
            sentinel++;
            if (sentinel == 0) sentinel = 1; // wrap around, skip 0
        }
    };
    
    // First, add original sequences with sentinels
    for (size_t i = 0; i < sequences.size(); ++i) {
        // Convert to uppercase
        std::string seq = sequences[i];
        for (char& c : seq) {
            if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        }
        result.prepared_string += seq;
        
        // Add sentinel and track its position
        size_t sentinel_pos = result.prepared_string.length();
        char sentinel = get_sentinel(i);
        result.prepared_string += sentinel;
        result.sentinel_positions.push_back(sentinel_pos);
    }
    
    result.original_length = result.prepared_string.length();
    
    // Then, add reverse complements with different sentinels
    for (int i = static_cast<int>(sequences.size()) - 1; i >= 0; --i) {
        // Convert to uppercase first
        std::string seq = sequences[i];
        for (char& c : seq) {
            if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        }
        
        // Add reverse complement
        std::string rc = revcomp(seq);
        result.prepared_string += rc;
        
        // Add sentinel and track its position (offset by sequences.size() to make them unique)
        size_t sentinel_pos = result.prepared_string.length();
        char sentinel = get_sentinel(sequences.size() + (sequences.size() - 1 - i));
        result.prepared_string += sentinel;
        result.sentinel_positions.push_back(sentinel_pos);
    }
    
    return result;
}

/**
 * @brief Prepares multiple DNA sequences for factorization without reverse complement and tracks sentinel positions.
 *
 * Takes multiple DNA sequences, concatenates them with unique sentinels, and tracks sentinel positions.
 * Unlike prepare_multiple_dna_sequences_w_rc(), this function does not append
 * reverse complements. The output format is: S = T1!T2@T3$
 *
 * @param sequences Vector of DNA sequence strings (should contain only A, C, T, G)
 * @return PreparedSequenceResult containing:
 *         - prepared_string: The formatted string with sequences and sentinels
 *         - original_length: Total length of the concatenated string (same as prepared_string.length())
 *         - sentinel_positions: Positions of all sentinels in the prepared string
 *
 * @throws std::invalid_argument If too many sequences (>250) or invalid nucleotides found
 * @throws std::runtime_error If sequences contain invalid characters
 *
 * @note Sentinels range from 1-251, avoiding 0, A(65), C(67), G(71), T(84)
 * @note The function validates that all sequences contain only valid DNA nucleotides
 * @note Input sequences can be lowercase or uppercase, output is always uppercase
 */
PreparedSequenceResult prepare_multiple_dna_sequences_no_rc(const std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return {"", 0, {}};
    }
    
    // Check if we have too many sequences
    if (sequences.size() > 250) {
        throw std::invalid_argument("Too many sequences: maximum 250 sequences supported (due to sentinel character limitations)");
    }
    
    // Validate sequences contain only valid DNA nucleotides
    for (size_t i = 0; i < sequences.size(); ++i) {
        for (char c : sequences[i]) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && 
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                throw std::runtime_error("Invalid nucleotide '" + std::string(1, c) + 
                                       "' found in sequence " + std::to_string(i));
            }
        }
    }
    
    PreparedSequenceResult result;
    
    // Calculate total size for reservation
    size_t total_size = 0;
    for (const auto& seq : sequences) {
        total_size += seq.length();
    }
    total_size += sequences.size(); // Add space for sentinels
    result.prepared_string.reserve(total_size);
    
    // Generate sentinel characters avoiding 0 and uppercase DNA nucleotides A(65), C(67), G(71), T(84)
    auto get_sentinel = [](size_t index) -> char {
        char sentinel = 1;
        size_t count = 0;
        
        while (true) {
            // Check if current sentinel is valid (not 0 and not uppercase DNA nucleotides)
            if (sentinel != 0 && sentinel != 'A' && sentinel != 'C' && sentinel != 'G' && sentinel != 'T') {
                if (count == index) {
                    return sentinel;  // Found the sentinel for this index
                }
                count++;
            }
            sentinel++;
            if (sentinel == 0) sentinel = 1; // wrap around, skip 0
        }
    };
    
    // Add original sequences with sentinels
    for (size_t i = 0; i < sequences.size(); ++i) {
        // Convert to uppercase
        std::string seq = sequences[i];
        for (char& c : seq) {
            if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        }
        result.prepared_string += seq;
        
        // Add sentinel and track its position (only between sequences, not after the last one)
        if (i < sequences.size() - 1) {
            size_t sentinel_pos = result.prepared_string.length();
            char sentinel = get_sentinel(i);
            result.prepared_string += sentinel;
            result.sentinel_positions.push_back(sentinel_pos);
        }
    }
    
    result.original_length = result.prepared_string.length();
    
    return result;
}


// ---------- generic, sink-driven noLZSS ----------

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
 * @return Number of factors emitted
 *
 * @note This is the core algorithm that all public functions use
 * @note The sink pattern allows for memory-efficient processing
 * @note All factors are emitted, including the last one
 */
template<class Sink>
static size_t nolzss(cst_t& cst, Sink&& sink) {
    sdsl::rmq_succinct_sct<> rmq(&cst.csa);
    const size_t str_len = cst.size() - 1; // the length of the string is the size of the CST minus the sentinel

    auto lambda = cst.select_leaf(cst.csa.isa[0] + 1);
    size_t lambda_node_depth = cst.node_depth(lambda);
    size_t lambda_sufnum = 0;

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
static size_t nolzss_dna_w_rc(const std::string& T, Sink&& sink) {
    const size_t n = T.size();
    if (n == 0) return 0;

    // Use prepare_multiple_dna_sequences_w_rc with a single sequence
    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc({T});
    std::string S = prep_result.prepared_string;
    size_t original_length = prep_result.original_length;

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
 */
template<class Sink>
static size_t nolzss_multiple_dna_w_rc(const std::string& S, Sink&& sink, size_t start_pos = 0) {
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

    const size_t T_beg = 0;
    const size_t T_end = N;           // end of original
    const size_t R_beg = N;       // first char of rc
    const size_t R_end = S.size();   // end of S

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
        // At factor start i (0-based in T), walk up ancestors and pick best candidate
        size_t best_len_depth = 0;   // best candidate's depth (proxy for length)
        bool   best_is_rc      = false;
        size_t best_fwd_start  = 0;  // start in T (for FWD)
        size_t best_rc_end     = 0;  // end in T (for RC)
        size_t best_rc_posS    = 0;  // pos in S where RC candidate suffix starts (for LCP)

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

            // Choose the better of the valid candidates at this depth
            if (okF) {
                if (ell > best_len_depth ||
                    (ell == best_len_depth && !best_is_rc && (jF + ell - 1) < (best_fwd_start + best_len_depth - 1))) {
                    best_len_depth = ell;
                    best_is_rc     = false;
                    best_fwd_start = jF;
                }
            }
            if (okR) {
                size_t posS_R = cst.csa[kR]; // suffix position in S for LCP
                if (ell > best_len_depth ||
                    (ell == best_len_depth && (best_is_rc ? (endRC < best_rc_end) : true))) {
                    best_len_depth = ell;
                    best_is_rc     = true;
                    best_rc_end    = endRC;
                    best_rc_posS   = posS_R;
                }
            }
        }

        size_t emit_len = 1;
        uint64_t emit_ref = i; // default for literal
        if (best_len_depth == 0) {
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

        if (!best_is_rc) {
            // Finalize FWD with true LCP and non-overlap cap
            size_t cap = i - best_fwd_start; // i-1 - (best_fwd_start) + 1
            size_t L   = lcp(cst, i, best_fwd_start);
            emit_len   = std::min(L, cap);
            emit_ref   = static_cast<uint64_t>(best_fwd_start);
        } else {
            // Finalize RC with true LCP (against suffix in R) 
            size_t L   = lcp(cst, i, best_rc_posS);
            emit_len   = L;
            size_t start_pos = best_rc_end - L + 2;
            emit_ref   = RC_MASK | static_cast<uint64_t>(start_pos); // start-anchored + RC flag
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

// ------------- public wrappers -------------

/**
 * @brief Factorizes a text string using the noLZSS algorithm.
 *
 * This is a template function that provides the core factorization functionality
 * for in-memory text. It builds a compressed suffix tree and applies the noLZSS
 * algorithm to find all factors.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param text Input text string
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This function copies the input string for suffix tree construction
 * @note For large inputs, consider using factorize_file_stream() instead
 * @see factorize() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_stream(std::string_view text, Sink&& sink) {
    std::string tmp(text);
    cst_t cst; construct_im(cst, tmp, 1);
    return nolzss(cst, std::forward<Sink>(sink));
}

/**
 * @brief Factorizes text from a file using the noLZSS algorithm.
 *
 * This template function reads text directly from a file and performs factorization
 * without loading the entire file into memory. This is more memory-efficient for
 * large files.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param path Path to input file containing text
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This function builds the suffix tree directly from the file
 * @see factorize_file() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_file_stream(const std::string& path, Sink&& sink) {
    // sdsl-lite will automatically add the sentinel when needed
    cst_t cst; construct(cst, path, 1);
    return nolzss(cst, std::forward<Sink>(sink));
}

/**
 * @brief Counts noLZSS factors in a text string.
 *
 * This function provides a convenient way to count factors without storing them.
 * It uses the sink-based factorization internally with a counting lambda.
 *
 * @param text Input text string
 * @return Number of factors in the factorization
 *
 * @note This is more memory-efficient than factorize() when you only need the count
 * @see factorize() for getting the actual factors
 * @see count_factors_file() for file-based counting
 */
size_t count_factors(std::string_view text) {
    size_t n = 0;
    factorize_stream(text, [&](const Factor&){ ++n; });
    return n;
}

/**
 * @brief Counts noLZSS factors in a file.
 *
 * This function reads text from a file and counts factors without storing them
 * or loading the entire file into memory. It's the most memory-efficient way
 * to get factor counts for large files.
 *
 * @param path Path to input file containing text
 * @return Number of factors in the factorization
 *
 * @note This function builds the suffix tree directly from the file
 * @see count_factors() for in-memory counting
 * @see factorize_file() for getting the actual factors from a file
 */
size_t count_factors_file(const std::string& path) {
    size_t n = 0;
    factorize_file_stream(path, [&](const Factor&){ ++n; });
    return n;
}

/**
 * @brief Factorizes a text string and returns factors as a vector.
 *
 * This is the main user-facing function for in-memory factorization.
 * It performs noLZSS factorization and returns all factors in a vector.
 *
 * @param text Input text string
 * @return Vector containing all factors from the factorization
 *
 * @note Factors are returned in order of appearance in the text
 * @note The returned factors are non-overlapping and cover the entire input
 * @see factorize_file() for file-based factorization
 */
std::vector<Factor> factorize(std::string_view text) {
    std::vector<Factor> out;
    factorize_stream(text, [&](const Factor& f){ out.push_back(f); });
    return out;
}

/**
 * @brief Factorizes text from a file and returns factors as a vector.
 *
 * This function reads text from a file, performs factorization, and returns
 * all factors in a vector. The reserve_hint parameter can improve performance
 * when you have an estimate of the number of factors.
 *
 * @param path Path to input file containing text
 * @param reserve_hint Optional hint for reserving space in output vector (0 = no hint)
 * @return Vector containing all factors from the factorization
 *
 * @note Use reserve_hint for better performance when you know approximate factor count
 * @note This is more memory-efficient than factorize() for large files
 * @see factorize() for in-memory factorization
 */
std::vector<Factor> factorize_file(const std::string& path, size_t reserve_hint) {
    std::vector<Factor> out;
    if (reserve_hint) out.reserve(reserve_hint);
    factorize_file_stream(path, [&](const Factor& f){ out.push_back(f); });
    return out;
}

/**
 * @brief Writes noLZSS factors from a file to a binary output file.
 *
 * This function reads text from an input file, performs factorization, and
 * writes the resulting factors in binary format to an output file. Each factor
 * is written as two uint64_t values (start position, length).
 *
 * @param in_path Path to input file containing text
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref)
 * @note This function overwrites the output file if it exists
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_binary_file(const std::string& in_path, const std::string& out_path) {
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Collect factors in memory to write header
    std::vector<Factor> factors;
    size_t n = factorize_file_stream(in_path, [&](const Factor& f){
        factors.push_back(f);
    });
    
    // For non-FASTA files, create minimal header with no sequence names or sentinels
    FactorFileHeader header;
    header.num_factors = factors.size();
    header.num_sequences = 0;  // Unknown for raw text files
    header.num_sentinels = 0;  // No sentinels for raw text files
    header.header_size = sizeof(FactorFileHeader);
    
    // Write header
    os.write(reinterpret_cast<const char*>(&header), sizeof(header));
    
    // Write factors
    for (const Factor& f : factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    return n;
}

/**
 * @brief Factorizes a DNA text string with reverse complement awareness using the noLZSS algorithm.
 *
 * This is a template function that provides the core factorization functionality
 * for in-memory DNA text, considering both forward and reverse complement matches.
 * It builds a compressed suffix tree and applies the noLZSS algorithm to find all factors.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param text Input DNA text string
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This function copies the input string for suffix tree construction
 * @note For large inputs, consider using factorize_file_stream_dna_w_rc() instead
 * @see factorize_dna_w_rc() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_stream_dna_w_rc(std::string_view text, Sink&& sink) {
    std::string tmp(text);
    return nolzss_dna_w_rc(tmp, std::forward<Sink>(sink));
}

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness using the noLZSS algorithm.
 *
 * This template function reads DNA text directly from a file and performs factorization
 * without loading the entire file into memory, considering both forward and reverse complement matches.
 * This is more memory-efficient for large files.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param path Path to input file containing DNA text
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This function builds the suffix tree directly from the file
 * @see factorize_file_dna_w_rc() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_file_stream_dna_w_rc(const std::string& path, Sink&& sink) {
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("Cannot open input file: " + path);
    std::string data((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    return nolzss_dna_w_rc(data, std::forward<Sink>(sink));
}

/**
 * @brief Factorizes a DNA text string with reverse complement awareness and returns factors as a vector.
 *
 * This is the main user-facing function for in-memory DNA factorization with reverse complement.
 * It performs noLZSS factorization and returns all factors in a vector.
 *
 * @param text Input DNA text string
 * @return Vector containing all factors from the factorization
 *
 * @note Factors are returned in order of appearance in the text
 * @note The returned factors are non-overlapping and cover the entire input
 * @see factorize_file_dna_w_rc() for file-based factorization
 */
std::vector<Factor> factorize_dna_w_rc(std::string_view text) {
    std::vector<Factor> out;
    factorize_stream_dna_w_rc(text, [&](const Factor& f){ out.push_back(f); });
    return out;
}

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness and returns factors as a vector.
 *
 * This function reads DNA text from a file, performs factorization with reverse complement, and returns
 * all factors in a vector. The reserve_hint parameter can improve performance
 * when you have an estimate of the number of factors.
 *
 * @param path Path to input file containing DNA text
 * @param reserve_hint Optional hint for reserving space in output vector (0 = no hint)
 * @return Vector containing all factors from the factorization
 *
 * @note Use reserve_hint for better performance when you know approximate factor count
 * @note This is more memory-efficient than factorize_dna_w_rc() for large files
 * @see factorize_dna_w_rc() for in-memory factorization
 */
std::vector<Factor> factorize_file_dna_w_rc(const std::string& path, size_t reserve_hint) {
    std::vector<Factor> out; if (reserve_hint) out.reserve(reserve_hint);
    factorize_file_stream_dna_w_rc(path, [&](const Factor& f){ out.push_back(f); });
    return out;
}

/**
 * @brief Counts noLZSS factors in a DNA text string with reverse complement awareness.
 *
 * This function provides a convenient way to count factors in DNA text without storing them.
 * It uses the sink-based factorization internally with a counting lambda.
 *
 * @param text Input DNA text string
 * @return Number of factors in the factorization
 *
 * @note This is more memory-efficient than factorize_dna_w_rc() when you only need the count
 * @see factorize_dna_w_rc() for getting the actual factors
 * @see count_factors_file_dna_w_rc() for file-based counting
 */
size_t count_factors_dna_w_rc(std::string_view text) {
    size_t n = 0; factorize_stream_dna_w_rc(text, [&](const Factor&){ ++n; }); return n;
}

/**
 * @brief Counts noLZSS factors in a DNA file with reverse complement awareness.
 *
 * This function reads DNA text from a file and counts factors without storing them
 * or loading the entire file into memory. It's the most memory-efficient way
 * to get factor counts for large DNA files.
 *
 * @param path Path to input file containing DNA text
 * @return Number of factors in the factorization
 *
 * @note This function builds the suffix tree directly from the file
 * @see count_factors_dna_w_rc() for in-memory counting
 * @see factorize_file_dna_w_rc() for getting the actual factors from a file
 */
size_t count_factors_file_dna_w_rc(const std::string& path) {
    size_t n = 0; factorize_file_stream_dna_w_rc(path, [&](const Factor&){ ++n; }); return n;
}

/**
 * @brief Writes noLZSS factors from a DNA file with reverse complement awareness to a binary output file.
 *
 * This function reads DNA text from an input file, performs factorization with reverse complement, and
 * writes the resulting factors in binary format to an output file. Each factor
 * is written as three uint64_t values (start position, length, ref).
 *
 * @param in_path Path to input file containing DNA text
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref)
 * @note This function overwrites the output file if it exists
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_binary_file_dna_w_rc(const std::string& in_path, const std::string& out_path) {
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Collect factors in memory to write header
    std::vector<Factor> factors;
    size_t n = factorize_file_stream_dna_w_rc(in_path, [&](const Factor& f){
        factors.push_back(f);
    });
    
    // For single DNA sequence files, create minimal header
    FactorFileHeader header;
    header.num_factors = factors.size();
    header.num_sequences = 1;  // Single DNA sequence
    header.num_sentinels = 0;  // No sentinels for single sequence
    header.header_size = sizeof(FactorFileHeader) + 1; // Empty sequence name
    
    // Write header
    os.write(reinterpret_cast<const char*>(&header), sizeof(header));
    
    // Write empty sequence name (single null terminator)
    os.write("\0", 1);
    
    // Write factors
    for (const Factor& f : factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    return n;
}

/**
 * @brief Factorizes a DNA text string with reverse complement awareness for multiple sequences using the noLZSS algorithm.
 *
 * This is a template function that provides the core factorization functionality
 * for in-memory DNA text with multiple sequences, considering both forward and reverse complement matches.
 * It builds a compressed suffix tree and applies the noLZSS algorithm to find all factors.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param text Input DNA text string with multiple sequences and sentinels
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This function copies the input string for suffix tree construction
 * @note For large inputs, consider using factorize_file_stream_multiple_dna_w_rc() instead
 * @see factorize_multiple_dna_w_rc() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_stream_multiple_dna_w_rc(std::string_view text, Sink&& sink) {
    std::string tmp(text);
    return nolzss_multiple_dna_w_rc(tmp, std::forward<Sink>(sink));
}

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness for multiple sequences using the noLZSS algorithm.
 *
 * This template function reads DNA text directly from a file and performs factorization
 * without loading the entire file into memory, considering both forward and reverse complement matches.
 * This is more memory-efficient for large files.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param path Path to input file containing DNA text with multiple sequences
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This function builds the suffix tree directly from the file
 * @see factorize_file_multiple_dna_w_rc() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_file_stream_multiple_dna_w_rc(const std::string& path, Sink&& sink) {
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("Cannot open input file: " + path);
    std::string data((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    return nolzss_multiple_dna_w_rc(data, std::forward<Sink>(sink));
}

/**
 * @brief Factorizes a DNA text string with reverse complement awareness for multiple sequences and returns factors as a vector.
 *
 * This is the main user-facing function for in-memory DNA factorization with multiple sequences and reverse complement.
 * It performs noLZSS factorization and returns all factors in a vector.
 *
 * @param text Input DNA text string with multiple sequences and sentinels
 * @return Vector containing all factors from the factorization
 *
 * @note Factors are returned in order of appearance in the text
 * @note The returned factors are non-overlapping and cover the entire input
 * @see factorize_file_multiple_dna_w_rc() for file-based factorization
 */
std::vector<Factor> factorize_multiple_dna_w_rc(std::string_view text) {
    std::vector<Factor> out;
    factorize_stream_multiple_dna_w_rc(text, [&](const Factor& f){ out.push_back(f); });
    return out;
}

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness for multiple sequences and returns factors as a vector.
 *
 * This function reads DNA text from a file, performs factorization with reverse complement for multiple sequences, and returns
 * all factors in a vector. The reserve_hint parameter can improve performance
 * when you have an estimate of the number of factors.
 *
 * @param path Path to input file containing DNA text with multiple sequences
 * @param reserve_hint Optional hint for reserving space in output vector (0 = no hint)
 * @return Vector containing all factors from the factorization
 *
 * @note Use reserve_hint for better performance when you know approximate factor count
 * @note This is more memory-efficient than factorize_multiple_dna_w_rc() for large files
 * @see factorize_multiple_dna_w_rc() for in-memory factorization
 */
std::vector<Factor> factorize_file_multiple_dna_w_rc(const std::string& path, size_t reserve_hint) {
    std::vector<Factor> out; if (reserve_hint) out.reserve(reserve_hint);
    factorize_file_stream_multiple_dna_w_rc(path, [&](const Factor& f){ out.push_back(f); });
    return out;
}

/**
 * @brief Counts noLZSS factors in a DNA text string with reverse complement awareness for multiple sequences.
 *
 * This function provides a convenient way to count factors in DNA text with multiple sequences without storing them.
 * It uses the sink-based factorization internally with a counting lambda.
 *
 * @param text Input DNA text string with multiple sequences
 * @return Number of factors in the factorization
 *
 * @note This is more memory-efficient than factorize_multiple_dna_w_rc() when you only need the count
 * @see factorize_multiple_dna_w_rc() for getting the actual factors
 * @see count_factors_file_multiple_dna_w_rc() for file-based counting
 */
size_t count_factors_multiple_dna_w_rc(std::string_view text) {
    size_t n = 0; factorize_stream_multiple_dna_w_rc(text, [&](const Factor&){ ++n; }); return n;
}

/**
 * @brief Counts noLZSS factors in a DNA file with reverse complement awareness for multiple sequences.
 *
 * This function reads DNA text from a file and counts factors without storing them
 * or loading the entire file into memory. It's the most memory-efficient way
 * to get factor counts for large DNA files with multiple sequences.
 *
 * @param path Path to input file containing DNA text with multiple sequences
 * @return Number of factors in the factorization
 *
 * @note This function builds the suffix tree directly from the file
 * @see count_factors_multiple_dna_w_rc() for in-memory counting
 * @see factorize_file_multiple_dna_w_rc() for getting the actual factors from a file
 */
size_t count_factors_file_multiple_dna_w_rc(const std::string& path) {
    size_t n = 0; factorize_file_stream_multiple_dna_w_rc(path, [&](const Factor&){ ++n; }); return n;
}

/**
 * @brief Writes noLZSS factors from a DNA file with reverse complement awareness for multiple sequences to a binary output file.
 *
 * This function reads DNA text from an input file, performs factorization with reverse complement for multiple sequences, and
 * writes the resulting factors in binary format to an output file. Each factor
 * is written as three uint64_t values (start position, length, ref).
 *
 * @param in_path Path to input file containing DNA text with multiple sequences
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref)
 * @note This function overwrites the output file if it exists
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_binary_file_multiple_dna_w_rc(const std::string& in_path, const std::string& out_path) {
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Collect factors in memory to write header
    std::vector<Factor> factors;
    size_t n = factorize_file_stream_multiple_dna_w_rc(in_path, [&](const Factor& f){
        factors.push_back(f);
    });
    
    // For multiple DNA sequences from text file, we don't know sequence names or sentinels
    FactorFileHeader header;
    header.num_factors = factors.size();
    header.num_sequences = 0;  // Unknown for raw text files with multiple sequences
    header.num_sentinels = 0;  // Cannot identify sentinels without preparation function
    header.header_size = sizeof(FactorFileHeader);
    
    // Write header
    os.write(reinterpret_cast<const char*>(&header), sizeof(header));
    
    // Write factors
    for (const Factor& f : factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    return n;
}

// Reference sequence factorization functions

std::vector<Factor> factorize_w_reference_seq(const std::string& reference_seq, const std::string& target_seq) {
    // Prepare reference and target sequences together
    std::vector<std::string> sequences = {reference_seq, target_seq};
    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc(sequences);
    
    // Calculate the starting position of the target sequence in the prepared string
    // The format is: REF[sentinel]TARGET[sentinel]RC(TARGET)[sentinel]RC(REF)[sentinel]
    // We want to start factorization from where TARGET begins
    size_t target_start_pos = reference_seq.length() + 1; // +1 for the sentinel between ref and target
    
    // Perform factorization starting from target sequence
    std::vector<Factor> factors;
    nolzss_multiple_dna_w_rc(prep_result.prepared_string, [&](const Factor& f) {
        factors.push_back(f);
    }, target_start_pos);
    
    return factors;
}

size_t factorize_w_reference_seq_file(const std::string& reference_seq, const std::string& target_seq, const std::string& out_path) {
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Get factors using the in-memory function
    std::vector<Factor> factors = factorize_w_reference_seq(reference_seq, target_seq);
    
    // Create header for reference+target factorization
    FactorFileHeader header;
    header.num_factors = factors.size();
    header.num_sequences = 2;  // Reference + target sequences
    header.num_sentinels = 1;  // One sentinel between ref and target (in factorized region)
    header.header_size = sizeof(FactorFileHeader);
    
    // Write header
    os.write(reinterpret_cast<const char*>(&header), sizeof(header));
    
    // Write factors
    for (const Factor& f : factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    return factors.size();
}

} // namespace noLZSS
