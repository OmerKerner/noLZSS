#include "parallel_fasta_processor.hpp"
#include "parallel_factorizer.hpp"
#include "factorizer_core.hpp"
#include "factorizer_helpers.hpp"
#include <fstream>
#include <atomic>
#include <sdsl/rmq_succinct_sct.hpp>
#include <sdsl/construct.hpp>
#include <thread>
#include <mutex>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;
namespace noLZSS {

/**
 * @brief Helper function to write metadata (sequence IDs and sentinel indices) to binary file
 * 
 * This is called after all factors have been written to the file.
 * It writes the sequence names, sentinel indices, and footer.
 * 
 * @param os Output stream positioned after all factors
 * @param sequence_ids Vector of sequence ID strings
 * @param sentinel_factor_indices Vector of sentinel factor indices
 * @param factor_count Total number of factors written
 * @param total_length Sum of all factor lengths
 */
static void write_fasta_metadata(std::ofstream& os,
                                 const std::vector<std::string>& sequence_ids,
                                 const std::vector<uint64_t>& sentinel_factor_indices,
                                 size_t factor_count,
                                 uint64_t total_length) {
    // Calculate footer size
    size_t names_size = 0;
    for (const auto& name : sequence_ids) {
        names_size += name.length() + 1;  // +1 for null terminator
    }
    
    size_t footer_size = sizeof(FactorFileFooter) + names_size + 
                        sentinel_factor_indices.size() * sizeof(uint64_t);
    
    // Write sequence names
    for (const auto& name : sequence_ids) {
        os.write(name.c_str(), name.length() + 1);  // Include null terminator
    }
    
    // Write sentinel factor indices
    for (uint64_t idx : sentinel_factor_indices) {
        os.write(reinterpret_cast<const char*>(&idx), sizeof(idx));
    }
    
    // Write footer at the end
    FactorFileFooter footer;
    footer.num_factors = factor_count;
    footer.num_sequences = sequence_ids.size();
    footer.num_sentinels = sentinel_factor_indices.size();
    footer.footer_size = footer_size;
    footer.total_length = total_length;
    
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
}

size_t parallel_write_factors_binary_file_fasta_multiple_dna_w_rc(
    const std::string& fasta_path, 
    const std::string& out_path,
    size_t num_threads) {
    
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);

    // Prepare sequences for factorization (this will validate nucleotides)
    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc(parse_result.sequences);
    
    // Use the core parallel function with the prepared string
    ParallelFactorizer factorizer;
    size_t factor_count = factorizer.parallel_factorize_multiple_dna_w_rc(
        prep_result.prepared_string,
        prep_result.original_length,
        out_path,
        num_threads);
    
    // Now identify sentinel factors by reading back the factors
    std::vector<uint64_t> sentinel_factor_indices;
    uint64_t total_length = 0;
    
    std::ifstream temp_is(out_path, std::ios::binary);
    if (!temp_is) {
        throw std::runtime_error("Cannot read output file to identify sentinels: " + out_path);
    }
    
    size_t sentinel_idx = 0;
    for (size_t i = 0; i < factor_count; ++i) {
        Factor f;
        temp_is.read(reinterpret_cast<char*>(&f), sizeof(Factor));
        total_length += f.length;
        
        // Check if this is a sentinel factor
        while (sentinel_idx < prep_result.sentinel_positions.size() &&
               prep_result.sentinel_positions[sentinel_idx] < f.start) {
            sentinel_idx++;
        }
        
        if (sentinel_idx < prep_result.sentinel_positions.size() &&
            f.start == prep_result.sentinel_positions[sentinel_idx]) {
            sentinel_factor_indices.push_back(i);
            sentinel_idx++;
        }
    }
    temp_is.close();
    
    // Truncate file to remove the basic footer written by merge_temp_files
    // We need to replace it with the metadata-enriched footer
    fs::resize_file(out_path, factor_count * sizeof(Factor));
    
    // Append metadata footer to the file
    std::ofstream os(out_path, std::ios::binary | std::ios::app);
    if (!os) {
        throw std::runtime_error("Cannot append metadata to output file: " + out_path);
    }
    
    write_fasta_metadata(os, parse_result.sequence_ids, sentinel_factor_indices, factor_count, total_length);
    
    return factor_count;
}

size_t parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
    const std::string& fasta_path,
    const std::string& out_path,
    size_t num_threads) {
    
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);

    // Prepare sequences for factorization (this will validate nucleotides)
    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_no_rc(parse_result.sequences);
    
    // Use ParallelFactorizer for the actual work
    ParallelFactorizer factorizer;
    size_t factor_count = factorizer.parallel_factorize(prep_result.prepared_string, out_path, num_threads);
    
    // Now we need to identify sentinel factors by reading back the factors
    std::vector<uint64_t> sentinel_factor_indices;
    uint64_t total_length = 0;
    
    // Read factors back to identify sentinels
    std::ifstream temp_is(out_path, std::ios::binary);
    if (!temp_is) {
        throw std::runtime_error("Cannot read output file to identify sentinels: " + out_path);
    }
    
    size_t sentinel_idx = 0;
    for (size_t i = 0; i < factor_count; ++i) {
        Factor f;
        temp_is.read(reinterpret_cast<char*>(&f), sizeof(Factor));
        total_length += f.length;
        
        // Check if this is a sentinel factor
        while (sentinel_idx < prep_result.sentinel_positions.size() &&
               prep_result.sentinel_positions[sentinel_idx] < f.start) {
            sentinel_idx++;
        }
        
        if (sentinel_idx < prep_result.sentinel_positions.size() &&
            f.start == prep_result.sentinel_positions[sentinel_idx]) {
            sentinel_factor_indices.push_back(i);
            sentinel_idx++;
        }
    }
    temp_is.close();
    
    // Truncate file to remove the basic footer written by merge_temp_files
    // We need to replace it with the metadata-enriched footer
    fs::resize_file(out_path, factor_count * sizeof(Factor));
    
    // Append metadata footer to the file
    std::ofstream os(out_path, std::ios::binary | std::ios::app);
    if (!os) {
        throw std::runtime_error("Cannot append metadata to output file: " + out_path);
    }
    
    write_fasta_metadata(os, parse_result.sequence_ids, sentinel_factor_indices, factor_count, total_length);
    
    return factor_count;
}

size_t parallel_write_factors_dna_w_reference_fasta_files_to_binary(
    const std::string& reference_fasta_path,
    const std::string& target_fasta_path,
    const std::string& out_path,
    size_t num_threads) {
    
    // Process both FASTA files and get prepared sequence with reverse complement
    FastaReferenceTargetResult ref_target_concat_w_rc = 
        prepare_ref_target_dna_w_rc_from_fasta(reference_fasta_path, target_fasta_path);
    
    // Use the core parallel function with the prepared string and start position
    ParallelFactorizer factorizer;
    size_t factor_count = factorizer.parallel_factorize_multiple_dna_w_rc(
        ref_target_concat_w_rc.concatinated_sequences.prepared_string,
        ref_target_concat_w_rc.concatinated_sequences.original_length,
        out_path,
        num_threads,
        ref_target_concat_w_rc.target_start_index);
    
    // Now identify sentinel factors by reading back the factors
    std::vector<uint64_t> sentinel_factor_indices;
    uint64_t total_length = 0;
    
    std::ifstream temp_is(out_path, std::ios::binary);
    if (!temp_is) {
        throw std::runtime_error("Cannot read output file to identify sentinels: " + out_path);
    }
    
    size_t sentinel_idx = 0;
    for (size_t i = 0; i < factor_count; ++i) {
        Factor f;
        temp_is.read(reinterpret_cast<char*>(&f), sizeof(Factor));
        total_length += f.length;
        
        // Check if this is a sentinel factor
        while (sentinel_idx < ref_target_concat_w_rc.concatinated_sequences.sentinel_positions.size() &&
               ref_target_concat_w_rc.concatinated_sequences.sentinel_positions[sentinel_idx] < f.start) {
            sentinel_idx++;
        }
        
        if (sentinel_idx < ref_target_concat_w_rc.concatinated_sequences.sentinel_positions.size() &&
            f.start == ref_target_concat_w_rc.concatinated_sequences.sentinel_positions[sentinel_idx]) {
            sentinel_factor_indices.push_back(i);
            sentinel_idx++;
        }
    }
    temp_is.close();
    
    // Truncate file to remove the basic footer written by merge_temp_files
    // We need to replace it with the metadata-enriched footer
    fs::resize_file(out_path, factor_count * sizeof(Factor));
    
    // Append metadata footer to the file
    std::ofstream os(out_path, std::ios::binary | std::ios::app);
    if (!os) {
        throw std::runtime_error("Cannot append metadata to output file: " + out_path);
    }
    
    write_fasta_metadata(os, ref_target_concat_w_rc.sequence_ids, sentinel_factor_indices, factor_count, total_length);
    
    return factor_count;
}

/**
 * @brief Helper function to write single sequence factors to binary file
 * 
 * Writes factors for a single sequence to a binary file with metadata footer.
 * 
 * @param factors Vector of factors for the sequence
 * @param sequence_id Sequence identifier
 * @param output_path Path to output file
 */
static void write_single_sequence_factors(const std::vector<Factor>& factors,
                                          const std::string& sequence_id,
                                          const std::string& output_path) {
    std::ofstream os(output_path, std::ios::binary | std::ios::trunc);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + output_path);
    }
    
    // Write factors
    uint64_t total_length = 0;
    for (const auto& f : factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
        total_length += f.length;
    }
    
    // Write sequence ID
    os.write(sequence_id.c_str(), sequence_id.length() + 1);  // Include null terminator
    
    // Write footer
    size_t footer_size = sizeof(FactorFileFooter) + sequence_id.length() + 1;
    
    FactorFileFooter footer;
    footer.num_factors = factors.size();
    footer.num_sequences = 1;  // Single sequence per file
    footer.num_sentinels = 0;  // No sentinels in per-sequence factorization
    footer.footer_size = footer_size;
    footer.total_length = total_length;
    
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
}

/**
 * @brief Helper function to sanitize sequence ID for use in filename
 * 
 * Replaces characters that are problematic in filenames with underscores.
 * 
 * @param seq_id Original sequence ID
 * @return Sanitized sequence ID safe for filenames
 */
static std::string sanitize_filename(const std::string& seq_id) {
    std::string safe_name = seq_id;
    for (char& c : safe_name) {
        // Replace problematic characters with underscore
        if (c == '/' || c == '\\' || c == ':' || c == '*' || c == '?' || 
            c == '"' || c == '<' || c == '>' || c == '|' || c == ' ') {
            c = '_';
        }
    }
    return safe_name;
}

size_t parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence(
    const std::string& fasta_path,
    const std::string& out_dir,
    size_t num_threads) {
    
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);
    
    size_t num_sequences = parse_result.sequences.size();
    
    // Determine actual number of threads to use
    if (num_threads == 0) {
        num_threads = std::min(num_sequences, static_cast<size_t>(std::thread::hardware_concurrency()));
    }
    num_threads = std::min(num_threads, num_sequences);
    
    // Create output directory if it doesn't exist
    fs::create_directories(out_dir);
    
    // Storage for results
    std::vector<std::vector<Factor>> all_factors(num_sequences);
    std::atomic<size_t> total_factor_count(0);
    
    // Parallel processing of sequences
    if (num_threads == 1) {
        // Sequential processing
        for (size_t i = 0; i < num_sequences; ++i) {
            std::vector<std::string> single_seq = {parse_result.sequences[i]};
            PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc(single_seq);
            all_factors[i] = factorize_multiple_dna_w_rc(prep_result.prepared_string);
            
            // Write to separate file
            std::string safe_id = sanitize_filename(parse_result.sequence_ids[i]);
            std::string output_path = out_dir + "/" + safe_id + ".bin";
            write_single_sequence_factors(all_factors[i], parse_result.sequence_ids[i], output_path);
            
            total_factor_count.fetch_add(all_factors[i].size());
        }
    } else {
        // Parallel processing
        std::vector<std::thread> threads;
        std::atomic<size_t> next_sequence(0);
        
        for (size_t t = 0; t < num_threads; ++t) {
            threads.emplace_back([&]() {
                while (true) {
                    size_t seq_idx = next_sequence.fetch_add(1);
                    if (seq_idx >= num_sequences) break;
                    
                    std::vector<std::string> single_seq = {parse_result.sequences[seq_idx]};
                    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc(single_seq);
                    std::vector<Factor> factors = factorize_multiple_dna_w_rc(prep_result.prepared_string);
                    
                    // Write to separate file
                    std::string safe_id = sanitize_filename(parse_result.sequence_ids[seq_idx]);
                    std::string output_path = out_dir + "/" + safe_id + ".bin";
                    write_single_sequence_factors(factors, parse_result.sequence_ids[seq_idx], output_path);
                    
                    total_factor_count.fetch_add(factors.size());
                }
            });
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
    }
    
    return total_factor_count.load();
}

size_t parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence(
    const std::string& fasta_path,
    const std::string& out_dir,
    size_t num_threads) {
    
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);
    
    size_t num_sequences = parse_result.sequences.size();
    
    // Determine actual number of threads to use
    if (num_threads == 0) {
        num_threads = std::min(num_sequences, static_cast<size_t>(std::thread::hardware_concurrency()));
    }
    num_threads = std::min(num_threads, num_sequences);
    
    // Create output directory if it doesn't exist
    fs::create_directories(out_dir);
    
    // Storage for results
    std::vector<std::vector<Factor>> all_factors(num_sequences);
    std::atomic<size_t> total_factor_count(0);
    
    // Parallel processing of sequences
    if (num_threads == 1) {
        // Sequential processing
        for (size_t i = 0; i < num_sequences; ++i) {
            std::vector<std::string> single_seq = {parse_result.sequences[i]};
            PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_no_rc(single_seq);
            // Remove the sentinel at the end
            std::string seq_without_sentinel = prep_result.prepared_string.substr(0, prep_result.prepared_string.length() - 1);
            all_factors[i] = factorize(seq_without_sentinel);
            
            // Write to separate file
            std::string safe_id = sanitize_filename(parse_result.sequence_ids[i]);
            std::string output_path = out_dir + "/" + safe_id + ".bin";
            write_single_sequence_factors(all_factors[i], parse_result.sequence_ids[i], output_path);
            
            total_factor_count.fetch_add(all_factors[i].size());
        }
    } else {
        // Parallel processing
        std::vector<std::thread> threads;
        std::atomic<size_t> next_sequence(0);
        
        for (size_t t = 0; t < num_threads; ++t) {
            threads.emplace_back([&]() {
                while (true) {
                    size_t seq_idx = next_sequence.fetch_add(1);
                    if (seq_idx >= num_sequences) break;
                    
                    std::vector<std::string> single_seq = {parse_result.sequences[seq_idx]};
                    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_no_rc(single_seq);
                    // Remove the sentinel at the end
                    std::string seq_without_sentinel = prep_result.prepared_string.substr(0, prep_result.prepared_string.length() - 1);
                    std::vector<Factor> factors = factorize(seq_without_sentinel);
                    
                    // Write to separate file
                    std::string safe_id = sanitize_filename(parse_result.sequence_ids[seq_idx]);
                    std::string output_path = out_dir + "/" + safe_id + ".bin";
                    write_single_sequence_factors(factors, parse_result.sequence_ids[seq_idx], output_path);
                    
                    total_factor_count.fetch_add(factors.size());
                }
            });
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
    }
    
    return total_factor_count.load();
}

} // namespace noLZSS
