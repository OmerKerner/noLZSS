#include "fasta_processor.hpp"
#include "parallel_fasta_processor.hpp"  // For parallel implementations
#include "factorizer.hpp"
#include "factorizer_core.hpp"      // Template implementations in detail:: namespace
#include "factorizer_helpers.hpp"  // For cst_t and lcp
#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/construct.hpp>
#include <iostream>
#include <algorithm>
#include <set>
#include <fstream>

namespace noLZSS {

// Helper function to parse FASTA file into individual sequences and IDs
FastaParseResult parse_fasta_sequences_and_ids(const std::string& fasta_path) {
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
std::vector<uint64_t> identify_sentinel_factors(const std::vector<Factor>& factors, 
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

// Note: Template implementations now in factorizer_core.hpp in the detail:: namespace
// No more local duplicates needed!

/**
 * @brief Writes noLZSS factors from multiple DNA sequences in a FASTA file with reverse complement awareness to a binary output file.
 *
 * This function reads DNA sequences from a FASTA file, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_w_rc(), performs 
 * factorization with reverse complement awareness, and writes the resulting factors in 
 * binary format to an output file with metadata including sequence IDs and sentinel factor indices.
 * Uses streaming to avoid storing all factors in memory.
 */
size_t write_factors_binary_file_fasta_multiple_dna_w_rc(const std::string& fasta_path, const std::string& out_path) {
    // Delegate to parallel version with 1 thread
    return parallel_write_factors_binary_file_fasta_multiple_dna_w_rc(fasta_path, out_path, 1);
}

/**
 * @brief Writes noLZSS factors from multiple DNA sequences in a FASTA file without reverse complement awareness to a binary output file.
 *
 * This function reads DNA sequences from a FASTA file, parses them into individual sequences,
 * prepares them for factorization using prepare_multiple_dna_sequences_no_rc(), performs 
 * factorization without reverse complement awareness, and writes the resulting factors in 
 * binary format to an output file with metadata including sequence IDs and sentinel factor indices.
 * Uses streaming to avoid storing all factors in memory.
 */
size_t write_factors_binary_file_fasta_multiple_dna_no_rc(const std::string& fasta_path, const std::string& out_path) {
    // Delegate to parallel version with 1 thread
    return parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(fasta_path, out_path, 1);
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
    // Delegate to parallel version with 1 thread
    return parallel_write_factors_dna_w_reference_fasta_files_to_binary(reference_fasta_path, target_fasta_path, out_path, 1);
}

/**
 * @brief Factorizes each DNA sequence in a FASTA file separately with reverse complement awareness.
 */
FastaPerSequenceFactorizationResult factorize_fasta_dna_w_rc_per_sequence(const std::string& fasta_path) {
    // Delegate to parallel version with 1 thread
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);
    
    FastaPerSequenceFactorizationResult result;
    result.sequence_ids = std::move(parse_result.sequence_ids);
    result.per_sequence_factors.reserve(parse_result.sequences.size());
    
    // Process sequentially (equivalent to parallel with 1 thread, but no threading overhead)
    for (const auto& sequence : parse_result.sequences) {
        std::vector<std::string> single_seq = {sequence};
        PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc(single_seq);
        std::vector<Factor> factors = factorize_multiple_dna_w_rc(prep_result.prepared_string);
        result.per_sequence_factors.push_back(std::move(factors));
    }
    
    return result;
}

/**
 * @brief Factorizes each DNA sequence in a FASTA file separately without reverse complement awareness.
 */
FastaPerSequenceFactorizationResult factorize_fasta_dna_no_rc_per_sequence(const std::string& fasta_path) {
    // Delegate to parallel version with 1 thread
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);
    
    FastaPerSequenceFactorizationResult result;
    result.sequence_ids = std::move(parse_result.sequence_ids);
    result.per_sequence_factors.reserve(parse_result.sequences.size());
    
    // Process sequentially (equivalent to parallel with 1 thread, but no threading overhead)
    for (const auto& sequence : parse_result.sequences) {
        std::vector<std::string> single_seq = {sequence};
        PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_no_rc(single_seq);
        std::string seq_without_sentinel = prep_result.prepared_string.substr(0, prep_result.prepared_string.length() - 1);
        std::vector<Factor> factors = factorize(seq_without_sentinel);
        result.per_sequence_factors.push_back(std::move(factors));
    }
    
    return result;
}

/**
 * @brief Writes factors from per-sequence DNA factorization with reverse complement to separate binary files.
 */
size_t write_factors_binary_file_fasta_dna_w_rc_per_sequence(const std::string& fasta_path, const std::string& out_dir) {
    // Delegate to parallel version with 1 thread
    return parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence(fasta_path, out_dir, 1);
}

/**
 * @brief Writes factors from per-sequence DNA factorization without reverse complement to separate binary files.
 */
size_t write_factors_binary_file_fasta_dna_no_rc_per_sequence(const std::string& fasta_path, const std::string& out_dir) {
    // Delegate to parallel version with 1 thread
    return parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence(fasta_path, out_dir, 1);
}

/**
 * @brief Counts per-sequence factors from DNA factorization with reverse complement.
 */
FastaPerSequenceCountResult count_factors_fasta_dna_w_rc_per_sequence(const std::string& fasta_path) {
    // Parse FASTA file into individual sequences
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);
    
    FastaPerSequenceCountResult result;
    result.sequence_ids = parse_result.sequence_ids;
    result.factor_counts.reserve(parse_result.sequences.size());
    
    // Count factors for each sequence independently
    for (const auto& sequence : parse_result.sequences) {
        // Prepare single DNA sequence with reverse complement
        std::vector<std::string> single_seq = {sequence};
        PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc(single_seq);
        
        // Count factors using count_factors_multiple_dna_w_rc
        size_t count = count_factors_multiple_dna_w_rc(prep_result.prepared_string);
        result.factor_counts.push_back(count);
        result.total_factors += count;
    }
    
    return result;
}

/**
 * @brief Counts per-sequence factors from DNA factorization without reverse complement.
 */
FastaPerSequenceCountResult count_factors_fasta_dna_no_rc_per_sequence(const std::string& fasta_path) {
    // Parse FASTA file into individual sequences
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);
    
    FastaPerSequenceCountResult result;
    result.sequence_ids = parse_result.sequence_ids;
    result.factor_counts.reserve(parse_result.sequences.size());
    
    // Count factors for each sequence independently
    for (const auto& sequence : parse_result.sequences) {
        // Prepare single DNA sequence (no reverse complement)
        std::vector<std::string> single_seq = {sequence};
        PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_no_rc(single_seq);
        
        // Remove the sentinel at the end
        std::string seq_without_sentinel = prep_result.prepared_string.substr(0, prep_result.prepared_string.length() - 1);
        
        // Count factors
        size_t count = count_factors(seq_without_sentinel);
        result.factor_counts.push_back(count);
        result.total_factors += count;
    }
    
    return result;
}

} // namespace noLZSS
