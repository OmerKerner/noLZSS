#include "fasta_processor.hpp"
#include "factorizer.hpp"
#include <iostream>
#include <algorithm>
#include <set>

namespace noLZSS {

// Helper function to parse FASTA file into individual sequences and IDs
static FastaParseResult parse_fasta_sequences_and_ids(const std::string& fasta_path) {
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + fasta_path);
    }

    FastaParseResult result;
    std::string line;
    std::string current_sequence;
    std::string current_id;

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
                result.sequences.push_back(current_sequence);
                result.sequence_ids.push_back(current_id);
                current_sequence.clear();
            }
            
            // Parse new header to extract ID
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
            } else {
                throw std::runtime_error("Empty sequence header in FASTA file");
            }
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
        result.sequences.push_back(current_sequence);
        result.sequence_ids.push_back(current_id);
    }

    file.close();

    if (result.sequences.empty()) {
        throw std::runtime_error("No valid sequences found in FASTA file");
    }

    return result;
}

// Helper function to identify sentinel factors from factorization results
static std::vector<uint64_t> identify_sentinel_factors(const std::vector<Factor>& factors, 
                                                      const std::vector<size_t>& sentinel_positions) {
    std::vector<uint64_t> sentinel_factor_indices;
    size_t sentinel_idx = 0;  // Current index in sentinel_positions
    
    for (size_t i = 0; i < factors.size(); ++i) {
        const Factor& f = factors[i];

        // Advance sentinel index past positions that occur before the current factor start
        while (sentinel_idx < sentinel_positions.size() &&
               sentinel_positions[sentinel_idx] < f.start) {
            sentinel_idx++;
        }

        // Check if this factor's start position matches current sentinel position
        if (sentinel_idx < sentinel_positions.size() &&
            f.start == sentinel_positions[sentinel_idx]) {

            // Sanity checks for sentinel factors
            if (f.length != 1) {
                throw std::runtime_error("Sentinel factor has unexpected length: " + std::to_string(f.length));
            }
            if (f.ref != f.start) {
                throw std::runtime_error("Sentinel factor reference mismatch: ref=" +
                                       std::to_string(f.ref) + ", pos=" + std::to_string(f.start));
            }
            sentinel_factor_indices.push_back(i);
            sentinel_idx++;  // Move to next sentinel position
        }
    }
    
    return sentinel_factor_indices;
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
 * @brief Prepares reference and target DNA sequences from FASTA files without reverse complement.
 *
 * Reads reference and target FASTA files, parses DNA sequences from each, and concatenates them with
 * reference sequences first, followed by target sequences, using sentinel characters between sequences.
 * Only nucleotides A, C, T, G are allowed (case insensitive, converted to uppercase).
 * This version does not append reverse complements.
 *
 * @param reference_fasta_path Path to the reference FASTA file
 * @param target_fasta_path Path to the target FASTA file
 * @return FastaReferenceTargetResult containing the prepared sequence data, sequence IDs, and counts of reference and target sequences
 */
FastaReferenceTargetResult prepare_ref_target_dna_no_rc_from_fasta(const std::string& reference_fasta_path,
                                                           const std::string& target_fasta_path) {
    // Process reference FASTA file first
    FastaParseResult ref_parse_result = parse_fasta_sequences_and_ids(reference_fasta_path);
    
    // Calculate target start index BEFORE moving sequences
    size_t target_start_index = 0;
    for (const auto& seq : ref_parse_result.sequences) {
        target_start_index += seq.length() + 1; // +1 for sentinel
    }
    
    size_t num_ref_sequences = ref_parse_result.sequences.size();
    
    // Process target FASTA file second
    FastaParseResult target_parse_result = parse_fasta_sequences_and_ids(target_fasta_path);
    
    size_t num_target_sequences = target_parse_result.sequences.size();
    
    // Reserve and move sequences (avoid copying)
    std::vector<std::string> all_original_sequences;
    all_original_sequences.reserve(num_ref_sequences + num_target_sequences);
    all_original_sequences.insert(all_original_sequences.end(),
                                 std::make_move_iterator(ref_parse_result.sequences.begin()),
                                 std::make_move_iterator(ref_parse_result.sequences.end()));
    all_original_sequences.insert(all_original_sequences.end(),
                                 std::make_move_iterator(target_parse_result.sequences.begin()),
                                 std::make_move_iterator(target_parse_result.sequences.end()));
    
    // Reserve and move IDs
    std::vector<std::string> all_sequence_ids;
    all_sequence_ids.reserve(num_ref_sequences + num_target_sequences);
    all_sequence_ids.insert(all_sequence_ids.end(),
                           std::make_move_iterator(ref_parse_result.sequence_ids.begin()),
                           std::make_move_iterator(ref_parse_result.sequence_ids.end()));
    all_sequence_ids.insert(all_sequence_ids.end(),
                           std::make_move_iterator(target_parse_result.sequence_ids.begin()),
                           std::make_move_iterator(target_parse_result.sequence_ids.end()));

    if (all_original_sequences.empty()) {
        throw std::runtime_error("No valid sequences found in FASTA files");
    }

    // Use prepare_multiple_dna_sequences_no_rc directly with collected sequences
    return {prepare_multiple_dna_sequences_no_rc(all_original_sequences), all_sequence_ids, num_ref_sequences, num_target_sequences, target_start_index};
}

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
 * @return FastaReferenceTargetResult containing the prepared sequence data, sequence IDs, and counts of reference and target sequences
 */
FastaReferenceTargetResult prepare_ref_target_dna_w_rc_from_fasta(const std::string& reference_fasta_path,
                                                          const std::string& target_fasta_path) {
    // Process reference FASTA file first
    FastaParseResult ref_parse_result = parse_fasta_sequences_and_ids(reference_fasta_path);
    
    // Calculate target start index BEFORE moving sequences
    size_t target_start_index = 0;
    for (const auto& seq : ref_parse_result.sequences) {
        target_start_index += seq.length() + 1; // +1 for sentinel
    }
    
    size_t num_ref_sequences = ref_parse_result.sequences.size();
    
    // Process target FASTA file second
    FastaParseResult target_parse_result = parse_fasta_sequences_and_ids(target_fasta_path);
    
    size_t num_target_sequences = target_parse_result.sequences.size();
    
    // Reserve and move sequences (avoid copying)
    std::vector<std::string> all_original_sequences;
    all_original_sequences.reserve(num_ref_sequences + num_target_sequences);
    all_original_sequences.insert(all_original_sequences.end(),
                                 std::make_move_iterator(ref_parse_result.sequences.begin()),
                                 std::make_move_iterator(ref_parse_result.sequences.end()));
    all_original_sequences.insert(all_original_sequences.end(),
                                 std::make_move_iterator(target_parse_result.sequences.begin()),
                                 std::make_move_iterator(target_parse_result.sequences.end()));
    
    // Reserve and move IDs
    std::vector<std::string> all_sequence_ids;
    all_sequence_ids.reserve(num_ref_sequences + num_target_sequences);
    all_sequence_ids.insert(all_sequence_ids.end(),
                           std::make_move_iterator(ref_parse_result.sequence_ids.begin()),
                           std::make_move_iterator(ref_parse_result.sequence_ids.end()));
    all_sequence_ids.insert(all_sequence_ids.end(),
                           std::make_move_iterator(target_parse_result.sequence_ids.begin()),
                           std::make_move_iterator(target_parse_result.sequence_ids.end()));

    if (all_original_sequences.empty()) {
        throw std::runtime_error("No valid sequences found in FASTA files");
    }

    // Use prepare_multiple_dna_sequences_w_rc directly with collected sequences
    return {prepare_multiple_dna_sequences_w_rc(all_original_sequences), all_sequence_ids, num_ref_sequences, num_target_sequences, target_start_index};
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
FastaFactorizationResult factorize_fasta_multiple_dna_w_rc(const std::string& fasta_path) {
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);

    // Prepare sequences for factorization (this will validate nucleotides)
    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc(parse_result.sequences);
    
    // Perform factorization
    std::vector<Factor> factors = factorize_multiple_dna_w_rc(prep_result.prepared_string);
    
    // Identify sentinel factors using helper function
    std::vector<uint64_t> sentinel_factor_indices = identify_sentinel_factors(factors, prep_result.sentinel_positions);
    
    return {factors, sentinel_factor_indices, parse_result.sequence_ids};
}

/**
 * @brief Factorizes multiple DNA sequences from a FASTA file without reverse complement awareness.
 *
 * Reads a FASTA file containing DNA sequences, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_no_rc(), and then
 * performs noLZSS factorization without reverse complement awareness.
 */
FastaFactorizationResult factorize_fasta_multiple_dna_no_rc(const std::string& fasta_path) {
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);

    // Prepare sequences for factorization (this will validate nucleotides)
    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_no_rc(parse_result.sequences);
    
    // Perform factorization using regular factorize function
    std::vector<Factor> factors = factorize(prep_result.prepared_string);
    
    // Identify sentinel factors using helper function
    std::vector<uint64_t> sentinel_factor_indices = identify_sentinel_factors(factors, prep_result.sentinel_positions);
    
    return {factors, sentinel_factor_indices, parse_result.sequence_ids};
}

/**
 * @brief Writes noLZSS factors from multiple DNA sequences in a FASTA file with reverse complement awareness to a binary output file.
 *
 * This function reads DNA sequences from a FASTA file, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_w_rc(), performs 
 * factorization with reverse complement awareness, and writes the resulting factors in 
 * binary format to an output file with metadata including sequence IDs and sentinel factor indices.
 */
size_t write_factors_binary_file_fasta_multiple_dna_w_rc(const std::string& fasta_path, const std::string& out_path) {
    // Get factorization result with sentinel information
    FastaFactorizationResult factorization_result = factorize_fasta_multiple_dna_w_rc(fasta_path);
    
    // Calculate footer size
    size_t names_size = 0;
    for (const auto& name : factorization_result.sequence_ids) {
        names_size += name.length() + 1;  // +1 for null terminator
    }
    
    size_t footer_size = sizeof(FactorFileFooter) + names_size + 
                        factorization_result.sentinel_factor_indices.size() * sizeof(uint64_t);
    
    // Write to file
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Write factors first
    for (const Factor& f : factorization_result.factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    // Write sequence names
    for (const auto& name : factorization_result.sequence_ids) {
        os.write(name.c_str(), name.length() + 1);  // Include null terminator
    }
    
    // Write sentinel factor indices
    for (uint64_t idx : factorization_result.sentinel_factor_indices) {
        os.write(reinterpret_cast<const char*>(&idx), sizeof(idx));
    }
    
    // Write footer at the end
    FactorFileFooter footer;
    footer.num_factors = factorization_result.factors.size();
    footer.num_sequences = factorization_result.sequence_ids.size();
    footer.num_sentinels = factorization_result.sentinel_factor_indices.size();
    footer.footer_size = footer_size;
    
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
    
    return factorization_result.factors.size();
}

/**
 * @brief Writes noLZSS factors from multiple DNA sequences in a FASTA file without reverse complement awareness to a binary output file.
 *
 * This function reads DNA sequences from a FASTA file, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_no_rc_with_sentinels(), performs 
 * factorization without reverse complement awareness, and writes the resulting factors in 
 * binary format to an output file with metadata including sequence IDs and sentinel factor indices.
 */
size_t write_factors_binary_file_fasta_multiple_dna_no_rc(const std::string& fasta_path, const std::string& out_path) {
    // Get factorization result with sentinel information
    FastaFactorizationResult factorization_result = factorize_fasta_multiple_dna_no_rc(fasta_path);
    
    // Calculate footer size
    size_t names_size = 0;
    for (const auto& name : factorization_result.sequence_ids) {
        names_size += name.length() + 1;  // +1 for null terminator
    }
    
    size_t footer_size = sizeof(FactorFileFooter) + names_size + 
                        factorization_result.sentinel_factor_indices.size() * sizeof(uint64_t);
    
    // Write to file
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Write factors first
    for (const Factor& f : factorization_result.factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    // Write sequence names
    for (const auto& name : factorization_result.sequence_ids) {
        os.write(name.c_str(), name.length() + 1);  // Include null terminator
    }
    
    // Write sentinel factor indices
    for (uint64_t idx : factorization_result.sentinel_factor_indices) {
        os.write(reinterpret_cast<const char*>(&idx), sizeof(idx));
    }
    
    // Write footer at the end
    FactorFileFooter footer;
    footer.num_factors = factorization_result.factors.size();
    footer.num_sequences = factorization_result.sequence_ids.size();
    footer.num_sentinels = factorization_result.sentinel_factor_indices.size();
    footer.footer_size = footer_size;
    
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
    
    return factorization_result.factors.size();
}

/**
 * @brief Factorizes DNA sequences from reference and target FASTA files with reverse complement awareness.
 *
 * Reads two FASTA files (reference and target), concatenates their sequences with sentinels,
 * prepares them for factorization using prepare_multiple_dna_sequences_w_rc(), and then
 * performs noLZSS factorization starting from the target sequences.
 *
 * @param reference_fasta_path Path to the reference FASTA file
 * @param target_fasta_path Path to the target FASTA file
 * @return FastaFactorizationResult containing the factors, sentinel factor indices, and sequence IDs
 */
FastaFactorizationResult factorize_dna_rc_w_ref_fasta_files(const std::string& reference_fasta_path, 
                                               const std::string& target_fasta_path) {
    // Process both FASTA files and get prepared sequence with reverse complement
    FastaReferenceTargetResult ref_target_concat_w_rc = prepare_ref_target_dna_w_rc_from_fasta(reference_fasta_path, target_fasta_path);
    
    // Perform factorization starting from the target start index
    std::vector<Factor> factors = factorize_multiple_dna_w_rc(ref_target_concat_w_rc.concatinated_sequences.prepared_string, 
                                                             ref_target_concat_w_rc.target_start_index);
    
    // Identify sentinel factors using helper function
    std::vector<uint64_t> sentinel_factor_indices = identify_sentinel_factors(factors, 
                                                                             ref_target_concat_w_rc.concatinated_sequences.sentinel_positions);
    
    return {factors, sentinel_factor_indices, ref_target_concat_w_rc.sequence_ids};
}


size_t write_factors_dna_w_reference_fasta_files_to_binary(const std::string& reference_fasta_path, 
                                                          const std::string& target_fasta_path, 
                                                          const std::string& out_path) {
    // Get factorization result with sentinel information
    FastaFactorizationResult factorization_result = factorize_dna_rc_w_ref_fasta_files(reference_fasta_path, target_fasta_path);

    // Calculate footer size
    size_t names_size = 0;
    for (const auto& name : factorization_result.sequence_ids) {
        names_size += name.length() + 1;  // +1 for null terminator
    }

    size_t footer_size = sizeof(FactorFileFooter) + names_size + 
                        factorization_result.sentinel_factor_indices.size() * sizeof(uint64_t);

    // Write to file
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }

    // Write factors first
    for (const Factor& f : factorization_result.factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }

    // Write sequence IDs
    for (const auto& name : factorization_result.sequence_ids) {
        os.write(name.c_str(), name.length() + 1);  // +1 for null terminator
    }

    // Write sentinel factor indices
    os.write(reinterpret_cast<const char*>(factorization_result.sentinel_factor_indices.data()),
             factorization_result.sentinel_factor_indices.size() * sizeof(uint64_t));

    // Write footer at the end
    FactorFileFooter footer;
    footer.num_factors = factorization_result.factors.size();
    footer.num_sequences = factorization_result.sequence_ids.size();
    footer.num_sentinels = factorization_result.sentinel_factor_indices.size();
    footer.footer_size = footer_size;

    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));

    return factorization_result.factors.size();
}

} // namespace noLZSS
