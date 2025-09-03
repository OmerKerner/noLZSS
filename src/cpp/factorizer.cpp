#include "factorizer.hpp"
#include <sdsl/suffix_trees.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <cassert>
#include <fstream>
#include <optional>

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

// ---------- generic, sink-driven LZSS ----------

/**
 * @brief Core LZSS factorization algorithm implementation.
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
static size_t lzss(cst_t& cst, Sink&& sink) {
    sdsl::rmq_succinct_sct<> rmq(&cst.csa);
    const size_t str_len = cst.size(); // Use the full string length including sentinel

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
                if (u == cst.root()) { l = 1; break; }
                else { l = u_depth; break; }
            }
            l = std::min(lcp(cst, lambda_sufnum, v_min_leaf_sufnum),
                         (lambda_sufnum - v_min_leaf_sufnum));
            if (l <= u_depth) { l = u_depth; break; }
            else { break; }
        }

        // Emit factor immediately
        Factor factor{static_cast<uint64_t>(lambda_sufnum), static_cast<uint64_t>(l)};
        sink(factor);
        ++count;

        // Advance to next position
        lambda = next_leaf(cst, lambda, l);
        lambda_node_depth = cst.node_depth(lambda);
        lambda_sufnum = cst.sn(lambda);
    }

    return count;
}

// ------------- public wrappers -------------

/**
 * @brief Factorizes a text string using the LZSS algorithm.
 *
 * This is a template function that provides the core factorization functionality
 * for in-memory text. It builds a compressed suffix tree and applies the LZSS
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
    // sdsl-lite will automatically add the sentinel when needed
    std::string tmp(text);
    cst_t cst; construct_im(cst, tmp, 1);
    return lzss(cst, std::forward<Sink>(sink));
}

/**
 * @brief Factorizes text from a file using the LZSS algorithm.
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
    return lzss(cst, std::forward<Sink>(sink));
}

/**
 * @brief Counts LZSS factors in a text string.
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
 * @brief Counts LZSS factors in a file.
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
 * It performs LZSS factorization and returns all factors in a vector.
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
 * @brief Writes LZSS factors from a file to a binary output file.
 *
 * This function reads text from an input file, performs factorization, and
 * writes the resulting factors in binary format to an output file. Each factor
 * is written as two uint64_t values (start position, length).
 *
 * @param in_path Path to input file containing text
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is 16 bytes (2 × uint64_t)
 * @note This function overwrites the output file if it exists
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_binary_file(const std::string& in_path, const std::string& out_path) {
    std::ofstream os(out_path, std::ios::binary);
    std::vector<char> buf(1<<20);
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    size_t n = factorize_file_stream(in_path, [&](const Factor& f){
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    });
    return n;
}

} // namespace noLZSS
