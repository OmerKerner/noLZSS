#pragma once
#include <vector>
#include <utility>
#include <string_view>
#include <string>
#include <cstdint>

namespace noLZSS {

// Constants and utility functions for reverse complement handling
/** @brief Mask for reverse complement flag in ref field (MSB of uint64_t) */
constexpr uint64_t RC_MASK = (1ULL << 63);

/**
 * @brief Result of preparing sequences with sentinel position tracking.
 */
struct PreparedSequenceResult {
    std::string prepared_string;              /**< Prepared string with sequences and sentinels */
    size_t original_length;                   /**< Length of original sequences (before RC if applicable) */
    std::vector<size_t> sentinel_positions;  /**< Positions of sentinels in prepared_string */
};

/**
 * @brief Binary file header for factor files with metadata.
 */
struct FactorFileHeader {
    char magic[8] = {'n', 'o', 'L', 'Z', 'S', 'S', 'v', '1'};  /**< Format identifier */
    uint64_t num_factors;     /**< Number of factors in file */
    uint64_t num_sequences;   /**< Number of sequences processed */
    uint64_t num_sentinels;   /**< Number of sentinel factors */
    uint64_t header_size;     /**< Total header size for extensibility */
};

/**
 * @brief Check if a reference value indicates a reverse complement match
 * @param ref The reference value to check
 * @return true if the MSB is set (reverse complement), false otherwise
 */
inline bool is_rc(uint64_t ref) { return (ref & RC_MASK) != 0; }

/**
 * @brief Extract the clean reference position by stripping the RC_MASK
 * @param ref The reference value (potentially with RC_MASK set)
 * @return The reference position with RC_MASK stripped
 */
inline uint64_t rc_end(uint64_t ref) { return (ref & ~RC_MASK); }

// Utility functions for DNA sequence preparation

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
PreparedSequenceResult prepare_multiple_dna_sequences_w_rc(const std::vector<std::string>& sequences);

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
PreparedSequenceResult prepare_multiple_dna_sequences_no_rc(const std::vector<std::string>& sequences);

/**
 * @brief Represents a factorization factor with start position, length, and reference position.
 *
 * A factor represents a substring in the original text that was identified
 * during noLZSS factorization. The factor covers text from position 'start'
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
 * @brief Factorizes a text string into noLZSS factors.
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
 * @brief Factorizes text from a file into noLZSS factors.
 *
 * Reads text from a file and performs noLZSS factorization. This is more memory-efficient
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
 * @brief Counts the number of noLZSS factors in a text string.
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
 * @brief Counts the number of noLZSS factors in a file.
 *
 * Reads text from a file and counts noLZSS factors without storing them.
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
 * @brief Writes noLZSS factors from a file to a binary output file.
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

// DNA-aware factorization functions with reverse complement support

/**
 * @brief Factorizes a DNA text string with reverse complement awareness into noLZSS factors.
 *
 * Performs non-overlapping Lempel-Ziv-Storer-Szymanski factorization on DNA sequences,
 * considering both forward and reverse complement matches. This is particularly useful
 * for genomic data where reverse complement patterns are biologically significant.
 *
 * @param text Input DNA text string
 * @return Vector of Factor objects representing the factorization
 *
 * @note Reverse complement matches are encoded with RC_MASK in the ref field
 * @note Factors are non-overlapping and cover the entire input
 * @see factorize_file_dna_w_rc() for file-based factorization
 */
std::vector<Factor> factorize_dna_w_rc(std::string_view text);

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness into noLZSS factors.
 *
 * Reads DNA text from a file and performs noLZSS factorization considering both forward
 * and reverse complement matches. This is more memory-efficient for large genomic files.
 *
 * @param path Path to the input file containing DNA text
 * @param reserve_hint Optional hint for reserving space in the output vector (0 = no hint)
 * @return Vector of Factor objects representing the factorization
 *
 * @note Use reserve_hint for better performance when you know approximate factor count
 * @see factorize_dna_w_rc() for in-memory factorization
 */
std::vector<Factor> factorize_file_dna_w_rc(const std::string& path, size_t reserve_hint = 0);

/**
 * @brief Counts the number of noLZSS factors in a DNA text string with reverse complement awareness.
 *
 * This is a memory-efficient alternative to factorize_dna_w_rc() when you only need
 * the count of factors rather than the factors themselves.
 *
 * @param text Input DNA text string
 * @return Number of factors in the factorization
 *
 * @see count_factors_file_dna_w_rc() for file-based counting
 */
size_t count_factors_dna_w_rc(std::string_view text);

/**
 * @brief Counts the number of noLZSS factors in a DNA file with reverse complement awareness.
 *
 * Reads DNA text from a file and counts noLZSS factors without storing them.
 * This is the most memory-efficient way to get factor counts for large genomic files.
 *
 * @param path Path to the input file containing DNA text
 * @return Number of factors in the factorization
 *
 * @see count_factors_dna_w_rc() for in-memory counting
 */
size_t count_factors_file_dna_w_rc(const std::string& path);

/**
 * @brief Writes noLZSS factors from a DNA file with reverse complement awareness to a binary output file.
 *
 * Reads DNA text from an input file, performs factorization with reverse complement support,
 * and writes the factors in binary format to an output file.
 *
 * @param in_path Path to input file containing DNA text
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is written as three uint64_t values (start, length, ref)
 * @note Reverse complement factors have RC_MASK set in the ref field
 * @warning This function overwrites the output file if it exists
 */
size_t write_factors_binary_file_dna_w_rc(const std::string& in_path, const std::string& out_path);

// Template functions for advanced usage

/**
 * @brief Advanced factorization function for DNA text with reverse complement awareness.
 *
 * This template function provides low-level access to the factorization process,
 * allowing custom handling of factors through a sink callable.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param text Input DNA text string
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 */
template<class Sink>
size_t factorize_stream_dna_w_rc(std::string_view text, Sink&& sink);

/**
 * @brief Advanced factorization function for DNA files with reverse complement awareness.
 *
 * This template function reads DNA text from a file and provides low-level access
 * to the factorization process through a sink callable.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param path Path to input file containing DNA text
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 */
template<class Sink>
size_t factorize_file_stream_dna_w_rc(const std::string& path, Sink&& sink);

// Multiple DNA sequences factorization functions with reverse complement support

/**
 * @brief Factorizes a DNA text string with reverse complement awareness for multiple sequences into noLZSS factors.
 *
 * Performs non-overlapping Lempel-Ziv-Storer-Szymanski factorization on DNA sequences
 * containing multiple sequences, considering both forward and reverse complement matches.
 * This is particularly useful for genomic data where reverse complement patterns are
 * biologically significant across multiple sequences.
 *
 * @param text Input DNA text string with multiple sequences and sentinels
 * @return Vector of Factor objects representing the factorization
 *
 * @note Reverse complement matches are encoded with RC_MASK in the ref field
 * @note Factors are non-overlapping and cover the entire input
 * @see factorize_file_multiple_dna_w_rc() for file-based factorization
 */
std::vector<Factor> factorize_multiple_dna_w_rc(std::string_view text);

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness for multiple sequences into noLZSS factors.
 *
 * Reads DNA text from a file containing multiple sequences and performs noLZSS factorization
 * considering both forward and reverse complement matches. This is more memory-efficient
 * for large genomic files with multiple sequences.
 *
 * @param path Path to the input file containing DNA text with multiple sequences
 * @param reserve_hint Optional hint for reserving space in the output vector (0 = no hint)
 * @return Vector of Factor objects representing the factorization
 *
 * @note Use reserve_hint for better performance when you know approximate factor count
 * @see factorize_multiple_dna_w_rc() for in-memory factorization
 */
std::vector<Factor> factorize_file_multiple_dna_w_rc(const std::string& path, size_t reserve_hint = 0);

/**
 * @brief Counts the number of noLZSS factors in a DNA text string with reverse complement awareness for multiple sequences.
 *
 * This is a memory-efficient alternative to factorize_multiple_dna_w_rc() when you only need
 * the count of factors rather than the factors themselves.
 *
 * @param text Input DNA text string with multiple sequences
 * @return Number of factors in the factorization
 *
 * @see count_factors_file_multiple_dna_w_rc() for file-based counting
 */
size_t count_factors_multiple_dna_w_rc(std::string_view text);

/**
 * @brief Counts the number of noLZSS factors in a DNA file with reverse complement awareness for multiple sequences.
 *
 * Reads DNA text from a file and counts noLZSS factors without storing them.
 * This is the most memory-efficient way to get factor counts for large genomic files
 * with multiple sequences.
 *
 * @param path Path to the input file containing DNA text with multiple sequences
 * @return Number of factors in the factorization
 *
 * @see count_factors_multiple_dna_w_rc() for in-memory counting
 */
size_t count_factors_file_multiple_dna_w_rc(const std::string& path);

/**
 * @brief Writes noLZSS factors from a DNA file with reverse complement awareness for multiple sequences to a binary output file.
 *
 * Reads DNA text from an input file containing multiple sequences, performs factorization
 * with reverse complement support, and writes the factors in binary format to an output file.
 *
 * @param in_path Path to input file containing DNA text with multiple sequences
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is written as three uint64_t values (start, length, ref)
 * @note Reverse complement factors have RC_MASK set in the ref field
 * @warning This function overwrites the output file if it exists
 */
size_t write_factors_binary_file_multiple_dna_w_rc(const std::string& in_path, const std::string& out_path);

// Template functions for advanced usage with multiple sequences

/**
 * @brief Advanced factorization function for DNA text with reverse complement awareness for multiple sequences.
 *
 * This template function provides low-level access to the factorization process
 * for multiple DNA sequences, allowing custom handling of factors through a sink callable.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param text Input DNA text string with multiple sequences and sentinels
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 */
template<class Sink>
size_t factorize_stream_multiple_dna_w_rc(std::string_view text, Sink&& sink);

/**
 * @brief Advanced factorization function for DNA files with reverse complement awareness for multiple sequences.
 *
 * This template function reads DNA text from a file containing multiple sequences
 * and provides low-level access to the factorization process through a sink callable.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param path Path to input file containing DNA text with multiple sequences
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 */
template<class Sink>
size_t factorize_file_stream_multiple_dna_w_rc(const std::string& path, Sink&& sink);

// Reference sequence factorization functions

/**
 * @brief Factorizes a target DNA sequence with reverse complement awareness using a reference sequence.
 *
 * Concatenates a reference sequence and target sequence (ref@target), then performs
 * noLZSS factorization with reverse complement awareness starting from where the target
 * sequence begins. This allows the target sequence to reference patterns in the reference
 * sequence without factorizing the reference itself.
 *
 * @param reference_seq Reference DNA sequence string
 * @param target_seq Target DNA sequence string to be factorized
 * @return Vector of Factor objects representing the factorization of the target sequence
 *
 * @note Factors start positions are absolute positions in the combined reference+target string
 * @note Both sequences should contain only A, C, T, G nucleotides (case insensitive)
 * @note Reverse complement matches are encoded with RC_MASK in the ref field
 * @throws std::invalid_argument If sequences contain invalid nucleotides
 */
std::vector<Factor> factorize_w_reference_seq(const std::string& reference_seq, const std::string& target_seq);

/**
 * @brief Factorizes a target DNA sequence using a reference sequence and writes factors to a binary file.
 *
 * Concatenates a reference sequence and target sequence (ref@target), then performs
 * noLZSS factorization with reverse complement awareness starting from where the target
 * sequence begins, and writes the resulting factors to a binary file.
 *
 * @param reference_seq Reference DNA sequence string
 * @param target_seq Target DNA sequence string to be factorized
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Factors start positions are absolute positions in the combined reference+target string
 * @note Both sequences should contain only A, C, T, G nucleotides (case insensitive)
 * @note Binary format follows the same structure as other DNA factorization binary outputs
 * @throws std::invalid_argument If sequences contain invalid nucleotides
 * @warning This function overwrites the output file if it exists
 */
size_t factorize_w_reference_seq_file(const std::string& reference_seq, const std::string& target_seq, const std::string& out_path);


} // namespace noLZSS
