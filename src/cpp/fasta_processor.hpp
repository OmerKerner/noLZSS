#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <cctype>

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

} // namespace noLZSS
