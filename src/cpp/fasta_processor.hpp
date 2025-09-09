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
 * prepares them for factorization using prepare_multiple_dna_sequences(), and then
 * performs noLZSS factorization with reverse complement awareness.
 *
 * @param fasta_path Path to the FASTA file containing DNA sequences
 * @return Vector containing all factors from the factorization
 *
 * @throws std::runtime_error If FASTA file cannot be opened or contains no valid sequences
 * @throws std::invalid_argument If too many sequences (>125) in the FASTA file or invalid nucleotides found
 *
 * @note Only A, C, T, G nucleotides are allowed (case insensitive)
 * @note Sequences are converted to uppercase before factorization
 * @note Reverse complement matches are supported during factorization
 * @note Nucleotide validation is performed by prepare_multiple_dna_sequences()
 */
std::vector<Factor> factorize_fasta_multiple_dna_w_rc(const std::string& fasta_path);

} // namespace noLZSS
