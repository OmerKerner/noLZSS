#include "fasta_processor.hpp"
#include "factorizer.hpp"
#include <iostream>
#include <algorithm>

namespace noLZSS {

// Helper function to parse FASTA file into individual sequences
static std::vector<std::string> parse_fasta_sequences(const std::string& fasta_path) {
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + fasta_path);
    }

    std::vector<std::string> sequences;
    std::string line;
    std::string current_sequence;

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
            if (!current_sequence.empty()) {
                sequences.push_back(current_sequence);
                current_sequence.clear();
            }
            // Skip header, continue to next line
        } else {
            // Sequence line - append to current sequence
            for (char c : line) {
                if (!std::isspace(c)) {
                    current_sequence += c;
                }
            }
        }
    }

    // Add the last sequence if it exists
    if (!current_sequence.empty()) {
        sequences.push_back(current_sequence);
    }

    file.close();

    if (sequences.empty()) {
        throw std::runtime_error("No valid sequences found in FASTA file");
    }

    return sequences;
}

// Helper function to extract sequence IDs from FASTA headers
std::vector<std::string> parse_fasta_sequence_ids(const std::string& fasta_path) {
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + fasta_path);
    }
    
    std::vector<std::string> sequence_ids;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Parse header to extract ID
            size_t start = 1; // Skip '>'
            while (start < line.size() && std::isspace(line[start])) {
                start++;
            }
            size_t end = start;
            while (end < line.size() && !std::isspace(line[end])) {
                end++;
            }
            
            if (start < line.size()) {
                sequence_ids.push_back(line.substr(start, end - start));
            }
        }
    }
    
    return sequence_ids;
}

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

/**
 * @brief Factorizes multiple DNA sequences from a FASTA file with reverse complement awareness.
 *
 * Reads a FASTA file containing DNA sequences, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_w_rc(), and then
 * performs noLZSS factorization with reverse complement awareness.
 */
std::vector<Factor> factorize_fasta_multiple_dna_w_rc(const std::string& fasta_path) {
    // Parse FASTA file into individual sequences
    std::vector<std::string> sequences = parse_fasta_sequences(fasta_path);

    // Prepare sequences for factorization (this will validate nucleotides)
    auto [prepared_string, original_length] = prepare_multiple_dna_sequences_w_rc(sequences);
    
    // Perform factorization
    return factorize_multiple_dna_w_rc(prepared_string);
}

/**
 * @brief Factorizes multiple DNA sequences from a FASTA file without reverse complement awareness.
 *
 * Reads a FASTA file containing DNA sequences, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_no_rc(), and then
 * performs noLZSS factorization without reverse complement awareness.
 */
std::vector<Factor> factorize_fasta_multiple_dna_no_rc(const std::string& fasta_path) {
    // Parse FASTA file into individual sequences
    std::vector<std::string> sequences = parse_fasta_sequences(fasta_path);

    // Prepare sequences for factorization (this will validate nucleotides)
    auto [prepared_string, total_length] = prepare_multiple_dna_sequences_no_rc(sequences);
    
    // Perform factorization using regular factorize function
    return factorize(prepared_string);
}

/**
 * @brief Writes noLZSS factors from multiple DNA sequences in a FASTA file with reverse complement awareness to a binary output file.
 *
 * This function reads DNA sequences from a FASTA file, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_w_rc(), performs 
 * factorization with reverse complement awareness, and writes the resulting factors in 
 * binary format to an output file. Each factor is written as three uint64_t values.
 */
size_t write_factors_binary_file_fasta_multiple_dna_w_rc(const std::string& fasta_path, const std::string& out_path) {
    // Parse FASTA file into individual sequences
    std::vector<std::string> sequences = parse_fasta_sequences(fasta_path);

    // Prepare sequences for factorization (this will validate nucleotides)
    auto [prepared_string, original_length] = prepare_multiple_dna_sequences_w_rc(sequences);
    
    // Perform factorization with reverse complement awareness
    std::vector<Factor> factors = factorize_multiple_dna_w_rc(prepared_string);
    
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Write factors to file
    for (const Factor& f : factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    return factors.size();
}

/**
 * @brief Writes noLZSS factors from multiple DNA sequences in a FASTA file without reverse complement awareness to a binary output file.
 *
 * This function reads DNA sequences from a FASTA file, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_no_rc(), performs 
 * factorization without reverse complement awareness, and writes the resulting factors in 
 * binary format to an output file. Each factor is written as three uint64_t values.
 */
size_t write_factors_binary_file_fasta_multiple_dna_no_rc(const std::string& fasta_path, const std::string& out_path) {
    // Parse FASTA file into individual sequences
    std::vector<std::string> sequences = parse_fasta_sequences(fasta_path);

    // Prepare sequences for factorization (this will validate nucleotides)
    auto [prepared_string, total_length] = prepare_multiple_dna_sequences_no_rc(sequences);
    
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Perform factorization and write factors to file using the regular factorize
    std::vector<Factor> factors = factorize(prepared_string);
    
    // Write factors to file
    for (const Factor& f : factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    return factors.size();
}

/**
 * @brief Factorizes multiple DNA sequences from a FASTA file with reverse complement awareness and returns metadata.
 */
FactorizationResult factorize_fasta_multiple_dna_w_rc_with_metadata(const std::string& fasta_path) {
    // Parse FASTA file into individual sequences and IDs
    std::vector<std::string> sequences = parse_fasta_sequences(fasta_path);
    std::vector<std::string> sequence_ids = parse_fasta_sequence_ids(fasta_path);
    
    // Prepare sequences for factorization (this will validate nucleotides)
    auto [prepared_string, original_length] = prepare_multiple_dna_sequences_w_rc(sequences);
    
    // Perform factorization
    std::vector<Factor> factors = factorize_multiple_dna_w_rc(prepared_string);
    
    FactorizationResult result;
    result.factors = std::move(factors);
    result.sequence_ids = std::move(sequence_ids);
    
    // Calculate sequence lengths and positions, identify sentinel factors
    std::vector<uint64_t> sentinel_indices;
    std::vector<size_t> sequence_lengths;
    std::vector<size_t> sequence_positions;
    
    size_t position = 0;
    for (size_t seq_idx = 0; seq_idx < sequences.size(); ++seq_idx) {
        sequence_lengths.push_back(sequences[seq_idx].length());
        sequence_positions.push_back(position);
        position += sequences[seq_idx].length();
        
        if (seq_idx < sequences.size() - 1) {  // Not the last sequence
            // Find factor at this position with length 1 (sentinel)
            for (size_t factor_idx = 0; factor_idx < result.factors.size(); ++factor_idx) {
                const Factor& f = result.factors[factor_idx];
                if (f.start == position && f.length == 1) {
                    // This is likely a sentinel factor
                    sentinel_indices.push_back(factor_idx);
                    break;
                }
            }
            position += 1;  // Skip sentinel character
        }
    }
    
    result.sentinel_factor_indices = std::move(sentinel_indices);
    result.sequence_lengths = std::move(sequence_lengths);
    result.sequence_positions = std::move(sequence_positions);
    
    return result;
}

/**
 * @brief Factorizes multiple DNA sequences from a FASTA file without reverse complement awareness and returns metadata.
 */
FactorizationResult factorize_fasta_multiple_dna_no_rc_with_metadata(const std::string& fasta_path) {
    // Parse FASTA file into individual sequences and IDs
    std::vector<std::string> sequences = parse_fasta_sequences(fasta_path);
    std::vector<std::string> sequence_ids = parse_fasta_sequence_ids(fasta_path);
    
    // Prepare sequences for factorization (this will validate nucleotides)
    auto [prepared_string, total_length] = prepare_multiple_dna_sequences_no_rc(sequences);
    
    // Perform factorization using regular factorize function
    std::vector<Factor> factors = factorize(prepared_string);
    
    FactorizationResult result;
    result.factors = std::move(factors);
    result.sequence_ids = std::move(sequence_ids);
    
    // Calculate sequence lengths and positions, identify sentinel factors
    std::vector<uint64_t> sentinel_indices;
    std::vector<size_t> sequence_lengths;
    std::vector<size_t> sequence_positions;
    
    size_t position = 0;
    for (size_t seq_idx = 0; seq_idx < sequences.size(); ++seq_idx) {
        sequence_lengths.push_back(sequences[seq_idx].length());
        sequence_positions.push_back(position);
        position += sequences[seq_idx].length();
        
        if (seq_idx < sequences.size() - 1) {  // Not the last sequence
            // Find factor at this position with length 1 (sentinel)
            for (size_t factor_idx = 0; factor_idx < result.factors.size(); ++factor_idx) {
                const Factor& f = result.factors[factor_idx];
                if (f.start == position && f.length == 1) {
                    // This is likely a sentinel factor
                    sentinel_indices.push_back(factor_idx);
                    break;
                }
            }
            position += 1;  // Skip sentinel character
        }
    }
    
    result.sentinel_factor_indices = std::move(sentinel_indices);
    result.sequence_lengths = std::move(sequence_lengths);
    result.sequence_positions = std::move(sequence_positions);
    
    return result;
}

/**
 * @brief Writes factorization result with metadata to extended binary format.
 */
size_t write_factorization_result_binary(const FactorizationResult& result, const std::string& out_path) {
    // Calculate header size
    size_t names_size = 0;
    for (const auto& name : result.sequence_ids) {
        names_size += name.length() + 1;  // +1 for null terminator
    }
    
    size_t header_size = sizeof(FactorFileHeader) + names_size + 
                        result.sentinel_factor_indices.size() * sizeof(uint64_t) +
                        result.sequence_lengths.size() * sizeof(size_t) +
                        result.sequence_positions.size() * sizeof(size_t);
    
    // Write to file
    std::ofstream os(out_path, std::ios::binary);
    std::vector<char> buf(1<<20); // 1 MB buffer
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Write header
    FactorFileHeader header;
    header.num_factors = result.factors.size();
    header.num_sequences = result.sequence_ids.size();
    header.num_sentinels = result.sentinel_factor_indices.size();
    header.header_size = header_size;
    header.reserved[0] = header.reserved[1] = header.reserved[2] = 0;
    
    os.write(reinterpret_cast<const char*>(&header), sizeof(header));
    
    // Write sequence names
    for (const auto& name : result.sequence_ids) {
        os.write(name.c_str(), name.length() + 1);  // Include null terminator
    }
    
    // Write sentinel indices
    for (uint64_t idx : result.sentinel_factor_indices) {
        os.write(reinterpret_cast<const char*>(&idx), sizeof(idx));
    }
    
    // Write sequence lengths
    for (size_t len : result.sequence_lengths) {
        os.write(reinterpret_cast<const char*>(&len), sizeof(len));
    }
    
    // Write sequence positions
    for (size_t pos : result.sequence_positions) {
        os.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
    }
    
    // Write factors (existing format)
    for (const Factor& f : result.factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    return result.factors.size();
}

} // namespace noLZSS
