#pragma once
#include "fasta_processor.hpp"
#include <string>

namespace noLZSS {

/**
 * @brief Parallel version of write_factors_binary_file_fasta_multiple_dna_w_rc
 * 
 * Reads a FASTA file containing DNA sequences, prepares them for factorization with
 * reverse complement awareness, and performs parallel factorization writing results
 * to a binary output file with metadata.
 * 
 * @param fasta_path Path to input FASTA file containing DNA sequences
 * @param out_path Path to output file where binary factors will be written
 * @param num_threads Number of threads to use (0 = auto-detect based on input size)
 * @return Number of factors written to the output file
 * 
 * @throws std::runtime_error If FASTA file cannot be opened or contains no valid sequences
 * @throws std::invalid_argument If too many sequences (>125) in the FASTA file or invalid nucleotides found
 * 
 * @note Binary format includes factors, sequence IDs, sentinel indices, and footer metadata
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note This function overwrites the output file if it exists
 * @note Reverse complement matches are supported during factorization
 * @note For single-threaded execution (num_threads=1), no temporary files are created
 * @warning Ensure sufficient disk space for the output file and temporary files
 */
size_t parallel_write_factors_binary_file_fasta_multiple_dna_w_rc(
    const std::string& fasta_path, 
    const std::string& out_path,
    size_t num_threads = 0
);

/**
 * @brief Parallel version of write_factors_binary_file_fasta_multiple_dna_no_rc
 * 
 * Reads a FASTA file containing DNA sequences, prepares them for factorization without
 * reverse complement awareness, and performs parallel factorization writing results
 * to a binary output file with metadata.
 * 
 * @param fasta_path Path to input FASTA file containing DNA sequences
 * @param out_path Path to output file where binary factors will be written
 * @param num_threads Number of threads to use (0 = auto-detect based on input size)
 * @return Number of factors written to the output file
 * 
 * @throws std::runtime_error If FASTA file cannot be opened or contains no valid sequences
 * @throws std::invalid_argument If too many sequences (>250) in the FASTA file or invalid nucleotides found
 * 
 * @note Binary format includes factors, sequence IDs, sentinel indices, and footer metadata
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note This function overwrites the output file if it exists
 * @note Reverse complement matches are NOT supported during factorization
 * @note For single-threaded execution (num_threads=1), no temporary files are created
 * @warning Ensure sufficient disk space for the output file and temporary files
 */
size_t parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
    const std::string& fasta_path,
    const std::string& out_path,
    size_t num_threads = 0
);

/**
 * @brief Parallel version of write_factors_dna_w_reference_fasta_files_to_binary
 * 
 * Reads DNA sequences from reference and target FASTA files, concatenates them with
 * sentinels, and performs parallel factorization starting from target sequences,
 * writing results to a binary output file with metadata.
 * 
 * @param reference_fasta_path Path to FASTA file containing reference DNA sequences
 * @param target_fasta_path Path to FASTA file containing target DNA sequences
 * @param out_path Path to output file where binary factors will be written
 * @param num_threads Number of threads to use (0 = auto-detect based on input size)
 * @return Number of factors written to the output file
 * 
 * @throws std::runtime_error If FASTA files cannot be opened or contain no valid sequences
 * @throws std::invalid_argument If too many sequences total or invalid nucleotides found
 * 
 * @note Binary format includes factors, sequence IDs, sentinel indices, and footer metadata
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note This function overwrites the output file if it exists
 * @note Reverse complement matches are supported during factorization
 * @note Factorization starts from target sequence positions only
 * @note For single-threaded execution (num_threads=1), no temporary files are created
 * @warning Ensure sufficient disk space for the output file and temporary files
 */
size_t parallel_write_factors_dna_w_reference_fasta_files_to_binary(
    const std::string& reference_fasta_path,
    const std::string& target_fasta_path,
    const std::string& out_path,
    size_t num_threads = 0
);

/**
 * @brief Parallel version of write_factors_binary_file_fasta_dna_w_rc_per_sequence
 * 
 * Reads a FASTA file, factorizes each sequence independently with reverse complement
 * awareness using parallel processing, and writes each sequence to a separate binary file.
 * 
 * @param fasta_path Path to input FASTA file containing DNA sequences
 * @param out_dir Path to output directory where binary factor files will be written
 * @param num_threads Number of threads to use (0 = auto-detect based on sequence count)
 * @return Total number of factors written across all sequences
 * 
 * @throws std::runtime_error If FASTA file cannot be opened or contains no valid sequences
 * @throws std::invalid_argument If invalid nucleotides found
 * 
 * @note Each sequence is factorized independently in parallel
 * @note Creates separate binary file for each sequence: <out_dir>/<seq_id>.bin
 * @note Binary format per file: factors + metadata footer
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note Reverse complement matches are supported during factorization
 * @warning Ensure sufficient disk space for the output files
 */
size_t parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence(
    const std::string& fasta_path,
    const std::string& out_dir,
    size_t num_threads = 0
);

/**
 * @brief Parallel version of write_factors_binary_file_fasta_dna_no_rc_per_sequence
 * 
 * Reads a FASTA file, factorizes each sequence independently without reverse complement
 * awareness using parallel processing, and writes each sequence to a separate binary file.
 * 
 * @param fasta_path Path to input FASTA file containing DNA sequences
 * @param out_dir Path to output directory where binary factor files will be written
 * @param num_threads Number of threads to use (0 = auto-detect based on sequence count)
 * @return Total number of factors written across all sequences
 * 
 * @throws std::runtime_error If FASTA file cannot be opened or contains no valid sequences
 * @throws std::invalid_argument If invalid nucleotides found
 * 
 * @note Each sequence is factorized independently in parallel
 * @note Creates separate binary file for each sequence: <out_dir>/<seq_id>.bin
 * @note Binary format per file: factors + metadata footer
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note Reverse complement matches are NOT supported during factorization
 * @warning Ensure sufficient disk space for the output files
 */
size_t parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence(
    const std::string& fasta_path,
    const std::string& out_dir,
    size_t num_threads = 0
);

} // namespace noLZSS
