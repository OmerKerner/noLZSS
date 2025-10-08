#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <cctype>
#include "factorizer.hpp"

namespace noLZSS {

/**
 * @brief Result of FASTA file processing containing concatenated sequences and metadata.
 */
struct FastaProcessResult {
    std::string sequence;      /**< Concatenated sequences with sentinels */
    size_t num_sequences;      /**< Number of sequences processed */
    std::vector<std::string> sequence_ids;  /**< IDs of processed sequences */
    std::vector<size_t> sequence_lengths;   /**< Lengths of each sequence (excluding sentinels) */
    std::vector<size_t> sequence_positions; /**< Start positions of each sequence in concatenated string */
    std::vector<size_t> file_boundaries;    /**< End positions of each file in concatenated string */
};

/**
 * @brief Result of parsing FASTA file containing both sequences and IDs.
 */
struct FastaParseResult {
    std::vector<std::string> sequences;    /**< Individual sequences from FASTA file */
    std::vector<std::string> sequence_ids; /**< Sequence identifiers from headers */
};

/**
 * @brief Result of FASTA factorization containing factors and sentinel information.
 */
struct FastaFactorizationResult {
    std::vector<Factor> factors;                    /**< Factorization result */
    std::vector<uint64_t> sentinel_factor_indices;  /**< Indices of factors that are sentinels */
    std::vector<std::string> sequence_ids;          /**< Sequence identifiers from FASTA headers */
};

/**
 * @brief Combined result of FASTA preparation containing prepared sequence and FASTA metadata.
 */
struct FastaReferenceTargetResult {
    PreparedSequenceResult concatinated_sequences;  /**< Prepared sequence with sentinels and metadata */
    std::vector<std::string> sequence_ids;          /**< Sequence identifiers from FASTA headers */
    size_t num_ref_sequences;                       /**< Number of reference sequences */
    size_t num_target_sequences;                    /**< Number of target sequences */
    size_t target_start_index;                      /**< Start index of target sequences in prepared string */
};

/**
 * @brief Processes a nucleotide FASTA file into a concatenated string with sentinels.
 *
 * Reads a FASTA file containing nucleotide sequences and creates a single concatenated
 * string with sentinel characters separating sequences. Only A, C, T, G nucleotides
 * are allowed (case insensitive, converted to uppercase).
 */
FastaProcessResult process_nucleotide_fasta(const std::string& fasta_path);

/**
 * @brief Prepares reference and target DNA sequences from FASTA files without reverse complement.
 *
 * Reads reference and target FASTA files, parses DNA sequences from each, and concatenates them with
 * reference sequences first, followed by target sequences, using sentinel characters between sequences.
 * Only nucleotides A, C, T, G are allowed (case insensitive, converted to uppercase).
 * This version does not append reverse complements.
 *
 * @param reference_fasta_path Path to the reference FASTA file
 * @param target_fasta_path Path to the target FASTA file
 * @return FastaReferenceTargetResult containing the prepared sequence data and sequence IDs
 */
FastaReferenceTargetResult prepare_ref_target_dna_no_rc_from_fasta(const std::string& reference_fasta_path,
                                                           const std::string& target_fasta_path);

/**
 * @brief Prepares reference and target DNA sequences from FASTA files with reverse complement.
 *
 * Reads reference and target FASTA files, parses DNA sequences from each, and prepares them using
 * prepare_multiple_dna_sequences_w_rc which concatenates sequences with sentinels
 * and appends reverse complements. Reference sequences come first, followed by target sequences.
 * Only nucleotides A, C, T, G are allowed.
 *
 * @param reference_fasta_path Path to the reference FASTA file
 * @param target_fasta_path Path to the target FASTA file
 * @return FastaReferenceTargetResult containing the prepared sequence data and sequence IDs
 */
FastaReferenceTargetResult prepare_ref_target_dna_w_rc_from_fasta(const std::string& reference_fasta_path,
                                                          const std::string& target_fasta_path);

/**
 * @brief Processes an amino acid FASTA file into a concatenated string with sentinels.
 *
 * Reads a FASTA file containing amino acid sequences and creates a single concatenated
 * string with sentinel characters separating sequences. Only canonical amino acids
 * are allowed (case insensitive, converted to uppercase).
 */
FastaProcessResult process_amino_acid_fasta(const std::string& fasta_path);

/**
 * @brief Factorizes multiple DNA sequences from a FASTA file with reverse complement awareness.
 *
 * Reads a FASTA file containing DNA sequences, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_w_rc(), and then
 * performs noLZSS factorization with reverse complement awareness.
 *
 * @param fasta_path Path to the FASTA file containing DNA sequences
 * @return FastaFactorizationResult containing factors and sentinel factor indices
 *
 * @throws std::runtime_error If FASTA file cannot be opened or contains no valid sequences
 * @throws std::invalid_argument If too many sequences (>125) in the FASTA file or invalid nucleotides found
 *
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note Sequences are converted to uppercase before factorization
 * @note Reverse complement matches are supported during factorization
 * @note Nucleotide validation is performed by prepare_multiple_dna_sequences_w_rc()
 */
FastaFactorizationResult factorize_fasta_multiple_dna_w_rc(const std::string& fasta_path);

/**
 * @brief Factorizes multiple DNA sequences from a FASTA file without reverse complement awareness.
 *
 * Reads a FASTA file containing DNA sequences, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_no_rc(), and then
 * performs noLZSS factorization without reverse complement awareness.
 *
 * @param fasta_path Path to the FASTA file containing DNA sequences
 * @return FastaFactorizationResult containing factors and sentinel factor indices
 *
 * @throws std::runtime_error If FASTA file cannot be opened or contains no valid sequences
 * @throws std::invalid_argument If too many sequences (>250) in the FASTA file or invalid nucleotides found
 *
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note Sequences are converted to uppercase before factorization
 * @note Reverse complement matches are NOT supported during factorization
 * @note Nucleotide validation is performed by prepare_multiple_dna_sequences_no_rc()
 */
FastaFactorizationResult factorize_fasta_multiple_dna_no_rc(const std::string& fasta_path);

/**
 * @brief Factorizes DNA sequences from reference and target FASTA files with reverse complement awareness.
 *
 * Reads two FASTA files (reference and target), concatenates their sequences with sentinels,
 * prepares them using prepare_multiple_dna_sequences_w_rc(), and performs factorization starting
 * from the target sequence region. Returns factors alongside sentinel metadata and sequence IDs.
 *
 * @param reference_fasta_path Path to the reference FASTA file
 * @param target_fasta_path Path to the target FASTA file
 * @return FastaFactorizationResult containing factors, sentinel indices, and sequence IDs
 */
FastaFactorizationResult factorize_dna_rc_w_ref_fasta_files(const std::string& reference_fasta_path,
                                                           const std::string& target_fasta_path);

/**
 * @brief Writes noLZSS factors from multiple DNA sequences in a FASTA file with reverse complement awareness to a binary output file.
 *
 * This function reads DNA sequences from a FASTA file, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_w_rc(), performs 
 * factorization with reverse complement awareness, and writes the resulting factors in 
 * binary format to an output file. Each factor is written as three uint64_t values.
 *
 * @param fasta_path Path to input FASTA file containing DNA sequences
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @throws std::runtime_error If FASTA file cannot be opened or contains no valid sequences
 * @throws std::invalid_argument If too many sequences (>125) in the FASTA file or invalid nucleotides found
 *
 * @note Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref)
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note This function overwrites the output file if it exists
 * @note Reverse complement matches are supported during factorization
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_binary_file_fasta_multiple_dna_w_rc(const std::string& fasta_path, const std::string& out_path);

/**
 * @brief Writes noLZSS factors from multiple DNA sequences in a FASTA file without reverse complement awareness to a binary output file.
 *
 * This function reads DNA sequences from a FASTA file, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_no_rc(), performs 
 * factorization without reverse complement awareness, and writes the resulting factors in 
 * binary format to an output file. Each factor is written as three uint64_t values.
 *
 * @param fasta_path Path to input FASTA file containing DNA sequences
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @throws std::runtime_error If FASTA file cannot be opened or contains no valid sequences
 * @throws std::invalid_argument If too many sequences (>250) in the FASTA file or invalid nucleotides found
 *
 * @note Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref)
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note This function overwrites the output file if it exists
 * @note Reverse complement matches are NOT supported during factorization
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_binary_file_fasta_multiple_dna_no_rc(const std::string& fasta_path, const std::string& out_path);

/**
 * @brief Writes noLZSS factors from DNA sequences in reference and target FASTA files to a binary output file.
 *
 * This function reads DNA sequences from two FASTA files (reference and target), concatenates them
 * with a sentinel separator, performs general factorization starting from the target sequences,
 * and writes the resulting factors in binary format to an output file.
 *
 * @param reference_fasta_path Path to FASTA file containing reference DNA sequences
 * @param target_fasta_path Path to FASTA file containing target DNA sequences
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @throws std::runtime_error If FASTA files cannot be opened or contain no valid sequences
 * @throws std::invalid_argument If too many sequences total or invalid nucleotides found
 *
 * @note Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref)
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note This function overwrites the output file if it exists
 * @note Uses general factorization (no reverse complement awareness)
 * @note Factorization starts from target sequence positions only
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_dna_w_reference_fasta_files_to_binary(const std::string& reference_fasta_path, 
                                                          const std::string& target_fasta_path, 
                                                          const std::string& out_path);

} // namespace noLZSS
