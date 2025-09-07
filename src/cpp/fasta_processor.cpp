#include "fasta_processor.hpp"
#include <iostream>
#include <algorithm>

namespace noLZSS {

/**
 * @brief Processes a nucleotide FASTA file into a concatenated string with sentinels.
 *
 * Reads a FASTA file containing nucleotide sequences and creates a single concatenated
 * string with sentinel characters separating sequences. Only A, C, T, G nucleotides
 * are allowed (case insensitive, converted to uppercase).
 */
FastaProcessResult process_nucleotide_fasta(const std::string& fasta_path) {
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + fasta_path);
    }

    FastaProcessResult result;
    result.num_sequences = 0;  // Initialize to 0 to avoid garbage values
    result.sequence.reserve(1024 * 1024); // Start with 1MB reservation

    std::string line;
    std::string current_id;
    size_t current_seq_length = 0;
    size_t current_seq_start = 0;
    bool in_sequence = false;

    // Generate sentinel characters (1-251, avoiding 0, A=65, C=67, G=71, T=84)
    auto get_sentinel = [](size_t seq_index) -> char {
        if (seq_index >= 251) {
            throw std::runtime_error("Too many sequences in FASTA file (max 251)");
        }
        char sentinel = static_cast<char>(seq_index + 1);
        // Skip nucleotide characters
        while (sentinel == 65 || sentinel == 67 || sentinel == 71 || sentinel == 84) { // A, C, G, T
            sentinel++;
            if (sentinel > 251) {
                throw std::runtime_error("Sentinel generation failed - too many sequences");
            }
        }
        return sentinel;
    };

    while (std::getline(file, line)) {
        // Remove trailing whitespace
        while (!line.empty() && std::isspace(line.back())) {
            line.pop_back();
        }

        if (line.empty()) {
            continue; // Skip empty lines
        }

        if (line[0] == '>') {
            // Header line - finish previous sequence if exists
            if (in_sequence && !current_id.empty()) {
                if (current_seq_length > 0) {
                    result.sequence_ids.push_back(current_id);
                    result.sequence_lengths.push_back(current_seq_length);
                    result.sequence_positions.push_back(current_seq_start);
                    result.num_sequences++;

                    // Add sentinel after sequence (except for last sequence)
                    char sentinel = get_sentinel(result.num_sequences - 1);
                    result.sequence.push_back(sentinel);
                }
            }

            // Parse new header
            size_t start = 1; // Skip '>'
            while (start < line.size() && std::isspace(line[start])) {
                start++;
            }
            size_t end = start;
            while (end < line.size() && !std::isspace(line[end])) {
                end++;
            }

            if (start < line.size()) {
                current_id = line.substr(start, end - start);
                current_seq_length = 0;
                current_seq_start = result.sequence.size();
                in_sequence = true;
            } else {
                throw std::runtime_error("Empty sequence header in FASTA file");
            }
        } else if (in_sequence) {
            // Sequence line - process each character
            for (char c : line) {
                if (std::isspace(c)) {
                    continue; // Skip whitespace
                }

                char upper_c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                if (upper_c == 'A' || upper_c == 'C' || upper_c == 'G' || upper_c == 'T') {
                    result.sequence.push_back(upper_c);
                    current_seq_length++;
                } else {
                    throw std::runtime_error("Invalid nucleotide character '" + std::string(1, c) +
                                           "' in sequence " + current_id);
                }
            }
        }
    }

    // Finish last sequence
    if (in_sequence && !current_id.empty() && current_seq_length > 0) {
        result.sequence_ids.push_back(current_id);
        result.sequence_lengths.push_back(current_seq_length);
        result.sequence_positions.push_back(current_seq_start);
        result.num_sequences++;
    }

    if (result.num_sequences == 0) {
        throw std::runtime_error("No valid sequences found in FASTA file");
    }

    file.close();
    return result;
}

/**
 * @brief Processes an amino acid FASTA file into a concatenated string with sentinels.
 *
 * Reads a FASTA file containing amino acid sequences and creates a single concatenated
 * string with sentinel characters separating sequences. Only canonical amino acids
 * are allowed (case insensitive, converted to uppercase).
 */
FastaProcessResult process_amino_acid_fasta(const std::string& fasta_path) {
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + fasta_path);
    }

    FastaProcessResult result;
    result.num_sequences = 0;  // Initialize to 0 to avoid garbage values
    result.sequence.reserve(1024 * 1024); // Start with 1MB reservation

    std::string line;
    std::string current_id;
    size_t current_seq_length = 0;
    size_t current_seq_start = 0;
    bool in_sequence = false;

    // Generate sentinel characters (1-251, avoiding 0 and amino acids)
    auto get_sentinel = [](size_t seq_index) -> char {
        if (seq_index >= 235) {  // 256 - 20 amino acids - 1 null = 235 available sentinels
            throw std::runtime_error("Too many sequences in FASTA file (max 235)");
        }

        // Amino acid ASCII values to avoid
        const std::string amino_acids = "ACDEFGHIKLMNPQRSTVWY";
        char sentinel = static_cast<char>(seq_index + 1);

        // Find next available character that's not an amino acid
        while (amino_acids.find(sentinel) != std::string::npos || sentinel == 0) {
            sentinel++;
            if (sentinel > 235) {
                throw std::runtime_error("Sentinel generation failed - no available characters");
            }
        }

        return sentinel;
    };

    // Canonical amino acids (20 standard)
    const std::string valid_aa = "ACDEFGHIKLMNPQRSTVWY";

    while (std::getline(file, line)) {
        // Remove trailing whitespace
        while (!line.empty() && std::isspace(line.back())) {
            line.pop_back();
        }

        if (line.empty()) {
            continue; // Skip empty lines
        }

        if (line[0] == '>') {
            // Header line - finish previous sequence if exists
            if (in_sequence && !current_id.empty()) {
                if (current_seq_length > 0) {
                    result.sequence_ids.push_back(current_id);
                    result.sequence_lengths.push_back(current_seq_length);
                    result.sequence_positions.push_back(current_seq_start);
                    result.num_sequences++;

                    // Add sentinel after sequence (except for last sequence)
                    char sentinel = get_sentinel(result.num_sequences - 1);
                    result.sequence.push_back(sentinel);
                }
            }

            // Parse new header
            size_t start = 1; // Skip '>'
            while (start < line.size() && std::isspace(line[start])) {
                start++;
            }
            size_t end = start;
            while (end < line.size() && !std::isspace(line[end])) {
                end++;
            }

            if (start < line.size()) {
                current_id = line.substr(start, end - start);
                current_seq_length = 0;
                current_seq_start = result.sequence.size();
                in_sequence = true;
            } else {
                throw std::runtime_error("Empty sequence header in FASTA file");
            }
        } else if (in_sequence) {
            // Sequence line - process each character
            for (char c : line) {
                if (std::isspace(c)) {
                    continue; // Skip whitespace
                }

                char upper_c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                if (valid_aa.find(upper_c) != std::string::npos) {
                    result.sequence.push_back(upper_c);
                    current_seq_length++;
                } else {
                    throw std::runtime_error("Invalid amino acid character '" + std::string(1, c) +
                                           "' in sequence " + current_id +
                                           ". Only canonical amino acids (ACDEFGHIKLMNPQRSTVWY) are allowed.");
                }
            }
        }
    }

    // Finish last sequence
    if (in_sequence && !current_id.empty() && current_seq_length > 0) {
        result.sequence_ids.push_back(current_id);
        result.sequence_lengths.push_back(current_seq_length);
        result.sequence_positions.push_back(current_seq_start);
        result.num_sequences++;
    }

    if (result.num_sequences == 0) {
        throw std::runtime_error("No valid sequences found in FASTA file");
    }

    file.close();
    return result;
}

} // namespace noLZSS
