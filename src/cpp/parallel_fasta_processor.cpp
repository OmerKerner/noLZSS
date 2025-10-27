#include "parallel_fasta_processor.hpp"
#include "parallel_factorizer.hpp"
#include "factorizer_core.hpp"
#include "factorizer_helpers.hpp"
#include <fstream>
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
 * @brief Helper function to write per-sequence factors metadata to binary file
 * 
 * This writes metadata for per-sequence factorization results, including
 * sequence boundaries (factor counts per sequence) and sequence IDs.
 * 
 * @param os Output stream positioned after all factors
 * @param sequence_ids Vector of sequence ID strings
 * @param factors_per_sequence Vector containing factor count for each sequence
 * @param total_factor_count Total number of factors across all sequences
 * @param total_length Sum of all factor lengths
 */
static void write_per_sequence_metadata(std::ofstream& os,
                                        const std::vector<std::string>& sequence_ids,
                                        const std::vector<size_t>& factors_per_sequence,
                                        size_t total_factor_count,
                                        uint64_t total_length) {
    // Calculate footer size
    size_t names_size = 0;
    for (const auto& name : sequence_ids) {
        names_size += name.length() + 1;  // +1 for null terminator
    }
    
    size_t footer_size = sizeof(FactorFileFooter) + names_size + 
                        factors_per_sequence.size() * sizeof(uint64_t);
    
    // Write sequence names
    for (const auto& name : sequence_ids) {
        os.write(name.c_str(), name.length() + 1);  // Include null terminator
    }
    
    // Write factors per sequence (as boundary markers)
    for (size_t count : factors_per_sequence) {
        uint64_t count_u64 = static_cast<uint64_t>(count);
        os.write(reinterpret_cast<const char*>(&count_u64), sizeof(count_u64));
    }
    
    // Write footer at the end
    FactorFileFooter footer;
    footer.num_factors = total_factor_count;
    footer.num_sequences = sequence_ids.size();
    footer.num_sentinels = 0;  // No sentinels in per-sequence factorization
    footer.footer_size = footer_size;
    footer.total_length = total_length;
    
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
}

size_t parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence(
    const std::string& fasta_path,
    const std::string& out_path,
    size_t num_threads) {
    
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);
    
    size_t num_sequences = parse_result.sequences.size();
    
    // Determine actual number of threads to use
    if (num_threads == 0) {
        num_threads = std::min(num_sequences, static_cast<size_t>(std::thread::hardware_concurrency()));
    }
    num_threads = std::min(num_threads, num_sequences);
    
    // Prepare output file
    std::ofstream os(out_path, std::ios::binary | std::ios::trunc);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    // Storage for results
    std::vector<std::vector<Factor>> all_factors(num_sequences);
    std::vector<size_t> factors_per_sequence(num_sequences);
    std::mutex write_mutex;
    
    // Parallel processing of sequences
    if (num_threads == 1) {
        // Sequential processing
        for (size_t i = 0; i < num_sequences; ++i) {
            std::vector<std::string> single_seq = {parse_result.sequences[i]};
            PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc(single_seq);
            all_factors[i] = factorize_dna_w_rc(prep_result.prepared_string);
            factors_per_sequence[i] = all_factors[i].size();
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
                    std::vector<Factor> factors = factorize_dna_w_rc(prep_result.prepared_string);
                    
                    std::lock_guard<std::mutex> lock(write_mutex);
                    all_factors[seq_idx] = std::move(factors);
                    factors_per_sequence[seq_idx] = all_factors[seq_idx].size();
                }
            });
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
    }
    
    // Write all factors sequentially and compute total length
    size_t total_factor_count = 0;
    uint64_t total_length = 0;
    
    for (const auto& factors : all_factors) {
        for (const auto& f : factors) {
            os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
            total_length += f.length;
        }
        total_factor_count += factors.size();
    }
    
    // Write metadata footer
    write_per_sequence_metadata(os, parse_result.sequence_ids, factors_per_sequence, 
                                total_factor_count, total_length);
    
    return total_factor_count;
}

size_t parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence(
    const std::string& fasta_path,
    const std::string& out_path,
    size_t num_threads) {
    
    // Parse FASTA file into individual sequences with IDs
    FastaParseResult parse_result = parse_fasta_sequences_and_ids(fasta_path);
    
    size_t num_sequences = parse_result.sequences.size();
    
    // Determine actual number of threads to use
    if (num_threads == 0) {
        num_threads = std::min(num_sequences, static_cast<size_t>(std::thread::hardware_concurrency()));
    }
    num_threads = std::min(num_threads, num_sequences);
    
    // Prepare output file
    std::ofstream os(out_path, std::ios::binary | std::ios::trunc);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    // Storage for results
    std::vector<std::vector<Factor>> all_factors(num_sequences);
    std::vector<size_t> factors_per_sequence(num_sequences);
    std::mutex write_mutex;
    
    // Parallel processing of sequences
    if (num_threads == 1) {
        // Sequential processing
        for (size_t i = 0; i < num_sequences; ++i) {
            std::vector<std::string> single_seq = {parse_result.sequences[i]};
            PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_no_rc(single_seq);
            // Remove the sentinel at the end
            std::string seq_without_sentinel = prep_result.prepared_string.substr(0, prep_result.prepared_string.length() - 1);
            all_factors[i] = factorize(seq_without_sentinel);
            factors_per_sequence[i] = all_factors[i].size();
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
                    
                    std::lock_guard<std::mutex> lock(write_mutex);
                    all_factors[seq_idx] = std::move(factors);
                    factors_per_sequence[seq_idx] = all_factors[seq_idx].size();
                }
            });
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
    }
    
    // Write all factors sequentially and compute total length
    size_t total_factor_count = 0;
    uint64_t total_length = 0;
    
    for (const auto& factors : all_factors) {
        for (const auto& f : factors) {
            os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
            total_length += f.length;
        }
        total_factor_count += factors.size();
    }
    
    // Write metadata footer
    write_per_sequence_metadata(os, parse_result.sequence_ids, factors_per_sequence, 
                                total_factor_count, total_length);
    
    return total_factor_count;
}

} // namespace noLZSS
