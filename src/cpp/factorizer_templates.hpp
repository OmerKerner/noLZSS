#ifndef FACTORIZER_TEMPLATES_HPP
#define FACTORIZER_TEMPLATES_HPP

#include "factorizer.hpp"
#include "factorizer_helpers.hpp"

namespace noLZSS {

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
 * @param start_pos Starting position for factorization (default 0)
 * @return Number of factors emitted
 *
 * @note This function copies the input string for suffix tree construction
 * @note For large inputs, consider using factorize_file_stream() instead
 * @see factorize() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_stream(std::string_view text, Sink&& sink, size_t start_pos) {
    std::string tmp(text);
    cst_t cst; construct_im(cst, tmp, 1);
    return nolzss(cst, std::forward<Sink>(sink), start_pos);
}

/**
 * @brief Factorize DNA sequences with reverse complement processing using streaming
 *
 * Template version that accepts a sink for processing factors without storing in memory.
 * Processes multiple DNA sequences concatenated with sentinels, handling both forward
 * and reverse complement strands.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param text Concatenated DNA sequences with sentinels
 * @param sentinel_positions Positions where sentinels were inserted
 * @param sink Callable that receives each computed factor
 * @param start_pos Starting position for factor positions (default 0)
 * @return Total number of factors emitted
 *
 * @note Memory-efficient for large genomic datasets
 * @see factorize_multiple_dna_w_rc() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_stream_multiple_dna_w_rc(
    std::string_view text,
    Sink&& sink,
    size_t start_pos
) {
    std::string tmp(text);
    return nolzss_multiple_dna_w_rc(tmp, std::forward<Sink>(sink), start_pos);
}

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness using streaming.
 *
 * Template version that accepts a sink for processing factors without storing in memory.
 * Reads DNA text from a file and processes it efficiently with reverse complement handling.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param filepath Path to the input file
 * @param sentinel_positions Positions where sentinels were inserted
 * @param sink Callable that receives each computed factor
 * @param start_pos Starting position for factor positions (default 0)
 * @return Total number of factors emitted
 *
 * @note Uses cache-efficient file reading for large files
 * @see factorize_file_multiple_dna_w_rc() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_file_stream_multiple_dna_w_rc(
    const std::string& filepath,
    Sink&& sink,
    size_t start_pos
) {
    // Read file into string first
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    std::string text((std::istreambuf_iterator<char>(file)),
                      std::istreambuf_iterator<char>());
    return nolzss_multiple_dna_w_rc(text, std::forward<Sink>(sink), start_pos);
}

} // namespace noLZSS

#endif // FACTORIZER_TEMPLATES_HPP
