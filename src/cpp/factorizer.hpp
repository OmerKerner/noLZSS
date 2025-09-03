#pragma once
#include <vector>
#include <utility>
#include <string_view>
#include <string>
#include <cstdint>

namespace noLZSS {

/**
 * @brief Represents a factorization factor with start position, length, and reference position.
 *
 * A factor represents a substring in the original text that was identified
 * during LZSS factorization. The factor covers text from position 'start'
 * with the specified 'length', and 'ref' indicates the position of the
 * previous occurrence (reference) for compression purposes.
 */
struct Factor {
    uint64_t start;   /**< Starting position of the factor in the original text */
    uint64_t length;  /**< Length of the factor substring */
    uint64_t ref;     /**< Reference position of the previous occurrence */
};

// Core factorization functions

/**
 * @brief Factorizes a text string into LZSS factors.
 *
 * Performs non-overlapping Lempel-Ziv-Storer-Szymanski factorization on the input text.
 * The algorithm uses a suffix tree to find the longest previous factors for each position.
 *
 * @param text Input text string
 * @return Vector of Factor objects representing the factorization
 *
 * @note Factors are non-overlapping and cover the entire input
 * @see factorize_file() for file-based factorization
 */
std::vector<Factor> factorize(std::string_view text);

/**
 * @brief Factorizes text from a file into LZSS factors.
 *
 * Reads text from a file and performs LZSS factorization. This is more memory-efficient
 * for large files as it avoids loading the entire file into memory.
 *
 * @param path Path to the input file containing text
 * @param reserve_hint Optional hint for reserving space in the output vector (0 = no hint)
 * @return Vector of Factor objects representing the factorization
 *
 * @note Use reserve_hint for better performance when you know approximate factor count
 * @see factorize() for in-memory factorization
 */
std::vector<Factor> factorize_file(const std::string& path, size_t reserve_hint = 0);

// Counting functions

/**
 * @brief Counts the number of LZSS factors in a text string.
 *
 * This is a memory-efficient alternative to factorize() when you only need
 * the count of factors rather than the factors themselves.
 *
 * @param text Input text string
 * @return Number of factors in the factorization
 *
 * @see count_factors_file() for file-based counting
 */
size_t count_factors(std::string_view text);

/**
 * @brief Counts the number of LZSS factors in a file.
 *
 * Reads text from a file and counts LZSS factors without storing them.
 * This is the most memory-efficient way to get factor counts for large files.
 *
 * @param path Path to the input file containing text
 * @return Number of factors in the factorization
 *
 * @see count_factors() for in-memory counting
 */
size_t count_factors_file(const std::string& path);

// Binary output

/**
 * @brief Writes LZSS factors from a file to a binary output file.
 *
 * Reads text from an input file, performs factorization, and writes the factors
 * in binary format to an output file. This is useful for storing factorizations
 * efficiently or for further processing.
 *
 * @param in_path Path to input file containing text
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is written as two uint64_t values (start, length)
 * @warning This function overwrites the output file if it exists
 */
size_t write_factors_binary_file(const std::string& in_path, const std::string& out_path);

}
