#include "parallel_factorizer.hpp"
#include "factorizer_helpers.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <chrono>

namespace fs = std::filesystem;
namespace noLZSS {
// Helper functions lcp() and next_leaf() are now in factorizer_helpers.hpp

/**
 * @brief Creates a unique temporary file path for a worker thread
 * 
 * Generates a unique temporary file path in the system's temp directory using
 * a timestamp and thread ID to ensure uniqueness across concurrent operations.
 * Thread 0 doesn't use this (writes directly to output), but threads 1+ use it.
 * 
 * @param thread_id The ID of the thread requesting a temp file path
 * @return std::string Absolute path to a unique temporary file
 */
std::string ParallelFactorizer::create_temp_file_path(size_t thread_id) {
    auto timestamp = std::chrono::system_clock::now().time_since_epoch().count();
    std::string temp_dir = fs::temp_directory_path().string();
    return temp_dir + "/noLZSS_temp_" + std::to_string(timestamp) + "_" + 
           std::to_string(thread_id) + ".bin";
}

/**
 * @brief Main parallel factorization function
 * 
 * Performs non-overlapping LZSS factorization using multiple threads with convergence
 * detection. The algorithm:
 * 1. Auto-detects optimal thread count based on input size (if num_threads=0)
 * 2. Builds a single compressed suffix tree (CST) shared by all threads
 * 3. Divides input text into chunks, one per thread
 * 4. Each thread factorizes its chunk independently, with thread 0 writing directly
 *    to the final output file for efficiency
 * 5. Threads detect convergence by monitoring when their factors align with the
 *    next thread's factors, allowing early termination
 * 6. Merges results from thread 1+ into the output file (thread 0 already there)
 * 7. Appends binary footer with metadata
 * 
 * OPTIMIZATION: Thread 0 writes directly to output_path, avoiding a copy operation.
 * For single-threaded execution, this means no temp files and no merge step.
 * 
 * @param text Input text to factorize (string_view for zero-copy efficiency)
 * @param output_path Path where the binary factor file will be written
 * @param num_threads Number of threads to use (0 = auto-detect based on input size)
 * @return size_t Total number of factors produced
 * 
 * @throws std::runtime_error If file I/O operations fail
 */
size_t ParallelFactorizer::parallel_factorize(std::string_view text, const std::string& output_path, 
                                           size_t num_threads) {
    if (text.empty()) return 0;
    
    // Determine optimal thread count if not specified
    if (num_threads == 0) {
        // Minimum characters per thread to make parallelization worthwhile
        constexpr size_t MIN_CHARS_PER_THREAD = 10000;
        
        // Calculate maximum useful threads based on text size
        size_t max_useful_threads = text.length() / MIN_CHARS_PER_THREAD;
        
        // Use available hardware threads, but don't exceed what's useful for this input size
        size_t hardware_threads = std::thread::hardware_concurrency();
        num_threads = std::min(hardware_threads, max_useful_threads);
        
        // Ensure at least one thread
        num_threads = std::max(1UL, num_threads);
    }
    
    // Create suffix tree for all threads to use
    std::string tmp(text);
    cst_t cst;
    construct_im(cst, tmp, 1);
    
    // Create RMQ support once for all threads
    sdsl::rmq_succinct_sct<> rmq(&cst.csa);
    
    // Create contexts for each thread
    std::vector<ThreadContext> contexts(num_threads);
    
    // Divide work among threads
    const size_t chunk_size = text.length() / num_threads;
    for (size_t i = 0; i < num_threads; ++i) {
        contexts[i].thread_id = i;
        contexts[i].start_pos = i * chunk_size;
        contexts[i].end_pos = (i + 1 < num_threads) ? (i + 1) * chunk_size : text.length();
        contexts[i].text_length = text.length();
        // Thread 0 writes directly to output file; others use temp files
        contexts[i].temp_file_path = (i == 0) ? output_path : create_temp_file_path(i);
        contexts[i].is_last_thread = (i == num_threads - 1);
    }
    
    // Create mutexes for file access
    std::vector<std::mutex> file_mutexes(num_threads);
    
    // Create and start worker threads
    std::vector<std::thread> threads;
    for (size_t i = 0; i < num_threads; ++i) {
        threads.emplace_back(&ParallelFactorizer::factorize_thread, this,
                            std::ref(cst), std::ref(rmq),
                            std::ref(contexts[i]), std::ref(contexts),
                            std::ref(file_mutexes));
    }
    
    // Wait for all threads to complete
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
    
    // Merge temporary files and create final output
    size_t total_factors = merge_temp_files(output_path, contexts);
    
    // Cleanup temporary files
    cleanup_temp_files(contexts);
    
    return total_factors;
}

/**
 * @brief Worker thread function that performs factorization on a text chunk
 * 
 * This is the core worker function executed by each thread. It performs non-overlapping
 * LZSS factorization on a specific chunk of the input text using the shared CST.
 * 
 * Algorithm details:
 * 1. Initializes at the thread's start position in the text
 * 2. Creates/truncates its output file (thread 0 uses final output, others use temp files)
 * 3. Iteratively computes factors using CST traversal:
 *    - Finds the longest previous match using level ancestors and RMQ
 *    - Ensures non-overlapping constraint (reference position + length <= current position)
 *    - Writes each factor to its designated file
 * 4. Checks for convergence with the next thread after each factor
 * 5. Stops when reaching text end or detecting convergence
 * 
 * Convergence detection: When a thread's factor extends into or past the next thread's
 * region, it checks if the next thread has reached the same factorization state. If so,
 * both threads will produce identical factors from that point, so this thread can stop.
 * 
 * @param cst Compressed suffix tree built from the entire input text (shared, read-only)
 * @param rmq Range minimum query structure for finding leftmost occurrences (shared, read-only)
 * @param ctx Thread's context containing start/end positions and file path
 * @param all_contexts Vector of all thread contexts (for accessing next thread's state)
 * @param file_mutexes Mutexes protecting file access for each thread
 */
void ParallelFactorizer::factorize_thread(const cst_t& cst, const sdsl::rmq_succinct_sct<>& rmq,
                                        ThreadContext& ctx,
                                        std::vector<ThreadContext>& all_contexts,
                                        std::vector<std::mutex>& file_mutexes) {
    // Initialize factorization at our starting position
    auto lambda = cst.select_leaf(cst.csa.isa[ctx.start_pos] + 1);
    size_t lambda_node_depth = cst.node_depth(lambda);
    size_t lambda_sufnum = ctx.start_pos;
    
    // Create or truncate output file
    {
        std::lock_guard<std::mutex> lock(file_mutexes[ctx.thread_id]);
        std::ofstream ofs(ctx.temp_file_path, std::ios::binary | std::ios::trunc);
    }
    
    // Track the next thread for convergence checking
    ThreadContext* next_ctx = nullptr;
    if (!ctx.is_last_thread && ctx.thread_id + 1 < all_contexts.size()) {
        next_ctx = &all_contexts[ctx.thread_id + 1];
    }
    
    // Main factorization loop
    while (lambda_sufnum < ctx.text_length) {
        // Compute current factor
        size_t d = 1;
        size_t l = 1;
        cst_t::node_type v;
        size_t v_min_leaf_sufnum = 0;
        size_t u_min_leaf_sufnum = 0;
        Factor current_factor;
        
        // Factor computation logic (similar to original nolzss)
        while (true) {
            v = cst.bp_support.level_anc(lambda, lambda_node_depth - d);
            v_min_leaf_sufnum = cst.csa[rmq(cst.lb(v), cst.rb(v))];
            l = cst.depth(v);
            
            if (v_min_leaf_sufnum + l - 1 < lambda_sufnum) {
                u_min_leaf_sufnum = v_min_leaf_sufnum;
                ++d;
                continue;
            }
            
            auto u = cst.parent(v);
            auto u_depth = cst.depth(u);
            
            if (v_min_leaf_sufnum == lambda_sufnum) {
                if (u == cst.root()) {
                    l = 1;
                    current_factor = Factor{static_cast<uint64_t>(lambda_sufnum), 
                                         static_cast<uint64_t>(l), 
                                         static_cast<uint64_t>(lambda_sufnum)};
                } else {
                    l = u_depth;
                    current_factor = Factor{static_cast<uint64_t>(lambda_sufnum), 
                                         static_cast<uint64_t>(l), 
                                         static_cast<uint64_t>(u_min_leaf_sufnum)};
                }
            } else {
                l = std::min(lcp(cst, lambda_sufnum, v_min_leaf_sufnum),
                           (lambda_sufnum - v_min_leaf_sufnum));
                if (l <= u_depth) {
                    l = u_depth;
                    current_factor = Factor{static_cast<uint64_t>(lambda_sufnum), 
                                         static_cast<uint64_t>(l), 
                                         static_cast<uint64_t>(u_min_leaf_sufnum)};
                } else {
                    current_factor = Factor{static_cast<uint64_t>(lambda_sufnum), 
                                         static_cast<uint64_t>(l), 
                                         static_cast<uint64_t>(v_min_leaf_sufnum)};
                }
            }
            break;
        }
        
        // Write the factor to temporary file
        write_factor(current_factor, ctx.temp_file_path, file_mutexes[ctx.thread_id]);
        
        // Check for convergence when factor extends into next thread's region
        // Use lambda_sufnum + l (exclusive end) to check against next thread's start
        if (next_ctx && lambda_sufnum + l >= next_ctx->start_pos) {
            size_t current_end = lambda_sufnum + l;  // Exclusive end position
            if (check_convergence(current_end, *next_ctx, file_mutexes[next_ctx->thread_id])) {
                break;  // Convergence detected
            }
        }
        
        // Advance to next position
        lambda = next_leaf(cst, lambda, l);
        lambda_node_depth = cst.node_depth(lambda);
        lambda_sufnum = cst.sn(lambda);
    }
}

/**
 * @brief Reads a single factor from a binary file at a specific index
 * 
 * Seeks to the byte position corresponding to the given factor index and reads
 * one Factor struct. Used during convergence checking to read factors from the
 * next thread's file without holding a lock for the entire read operation.
 * 
 * @param file_path Path to the binary factor file
 * @param factor_index Zero-based index of the factor to read
 * @return std::optional<Factor> The factor if successfully read, std::nullopt otherwise
 *         (file doesn't exist, index out of bounds, or incomplete read)
 */
std::optional<Factor> ParallelFactorizer::read_factor_at_index(const std::string& file_path, 
                                                                size_t factor_index) {
    std::ifstream ifs(file_path, std::ios::binary);
    if (!ifs) {
        return std::nullopt;
    }
    
    // Seek to the position of the requested factor
    ifs.seekg(factor_index * sizeof(Factor), std::ios::beg);
    if (!ifs) {
        return std::nullopt; // Index out of bounds
    }
    
    // Read the factor
    Factor factor;
    ifs.read(reinterpret_cast<char*>(&factor), sizeof(Factor));
    if (!ifs || ifs.gcount() != sizeof(Factor)) {
        return std::nullopt; // Failed to read complete factor
    }
    
    return factor;
}

/**
 * @brief Checks if this thread has converged with the next thread's factorization
 * 
 * Convergence occurs when the current thread's position aligns with the next thread's
 * factorization such that both would produce identical factors from that point forward.
 * This allows the current thread to stop early, avoiding redundant computation.
 * 
 * The algorithm:
 * 1. First checks the cached factor from the next thread (if available)
 * 2. If cache is insufficient, reads factors from the next thread's file sequentially
 * 3. Looks for exact alignment: next_factor.start == current_end
 * 4. Caches the last read factor for future checks
 * 
 * This function is called after each factor is written, when that factor extends into
 * or past the next thread's starting region.
 * 
 * @param current_end The exclusive end position of the current thread's last factor
 * @param next_ctx Reference to the next thread's context (for reading its factors)
 * @param next_file_mutex Mutex protecting access to the next thread's file
 * @return true if convergence is detected (current thread should stop)
 * @return false if convergence not yet reached (current thread should continue)
 */
bool ParallelFactorizer::check_convergence(size_t current_end, ThreadContext& next_ctx,
                                        std::mutex& next_file_mutex) {
    // If we have a cached factor, check it first
    if (next_ctx.last_read_factor.has_value()) {
        const auto& cached_factor = next_ctx.last_read_factor.value();
        size_t cached_end = cached_factor.start + cached_factor.length;
        
        if (cached_end > current_end) {
            // The cached factor ends beyond our current position,
            // so we haven't reached convergence yet
            return false;
        }
        
        if (cached_factor.start == current_end) {
            // Found convergence: next thread's factor starts exactly where we end
            return true;
        }
        
        // cached_end <= current_end but cached_factor.start < current_end
        // Need to read the next factor
    }
    
    // Read factors one at a time until we find convergence or go past current_end
    while (true) {
        std::lock_guard<std::mutex> lock(next_file_mutex);
        
        auto factor_opt = read_factor_at_index(next_ctx.temp_file_path, 
                                               next_ctx.next_thread_factor_index);
        
        if (!factor_opt.has_value()) {
            // No more factors in next thread's file - convergence not found yet
            return false;
        }
        
        Factor factor = factor_opt.value();
        next_ctx.last_read_factor = factor;
        next_ctx.next_thread_factor_index++;
        
        if (factor.start == current_end) {
            // Found convergence
            return true;
        }
        
        if (factor.start > current_end) {
            // Next thread's factor starts beyond our current position
            // We haven't reached convergence yet
            return false;
        }
        
        // factor.start < current_end, keep reading next factors
    }
}

/**
 * @brief Writes a single factor to a binary file (thread-safe)
 * 
 * Appends a Factor struct to the specified file in binary format. Uses a mutex
 * to ensure thread-safe access when multiple threads might write to the same file
 * (though in practice, each thread writes to its own file).
 * 
 * The factor is written as a raw binary struct with no serialization overhead.
 * 
 * @param factor The Factor struct to write
 * @param file_path Path to the output file (opened in append mode)
 * @param file_mutex Mutex protecting access to this file
 * @throws std::runtime_error If the file cannot be opened for writing
 */
void ParallelFactorizer::write_factor(const Factor& factor, const std::string& file_path, 
                                   std::mutex& file_mutex) {
    std::lock_guard<std::mutex> lock(file_mutex);
    
    std::ofstream ofs(file_path, std::ios::binary | std::ios::app);
    if (!ofs) {
        throw std::runtime_error("Cannot open temporary file for writing: " + file_path);
    }
    
    ofs.write(reinterpret_cast<const char*>(&factor), sizeof(Factor));
}

/**
 * @brief Reads a factor from a file at a specific index (thread-safe version)
 * 
 * Similar to read_factor_at_index(), but acquires a mutex before reading.
 * This ensures thread-safe access when multiple threads might read from the same file.
 * 
 * NOTE: Currently unused in favor of read_factor_at_index() which doesn't hold
 * the lock during the entire read operation for better concurrency.
 * 
 * @param file_path Path to the binary factor file
 * @param index Zero-based index of the factor to read
 * @param file_mutex Mutex protecting access to this file
 * @return std::optional<Factor> The factor if successfully read, std::nullopt otherwise
 */
std::optional<Factor> ParallelFactorizer::read_factor_at(const std::string& file_path, 
                                                     size_t index, 
                                                     std::mutex& file_mutex) {
    std::lock_guard<std::mutex> lock(file_mutex);
    
    std::ifstream ifs(file_path, std::ios::binary);
    if (!ifs) {
        return std::nullopt;
    }
    
    // Seek to the position of the requested factor
    ifs.seekg(index * sizeof(Factor));
    
    Factor factor;
    if (ifs.read(reinterpret_cast<char*>(&factor), sizeof(Factor))) {
        return factor;
    }
    
    return std::nullopt;
}

/**
 * @brief Merges temporary factor files and appends footer to create final output
 * 
 * This function handles the final stage of parallel factorization by:
 * 1. Reading thread 0's factors (already in output_path) to count them and find the end position
 * 2. If multiple threads and sequence not complete: appending factors from threads 1+ after
 *    finding the convergence point where their factorization aligns with prior threads
 * 3. Writing the binary footer with metadata (factor count, sequences, sentinels)
 * 
 * OPTIMIZATION: Thread 0 wrote directly to output_path during factorization, so we only
 * need to append from other threads. For single-threaded execution, we just count factors
 * and add the footer.
 * 
 * Convergence-based merging: For each subsequent thread, we skip factors until we find
 * one that starts exactly where the previous thread's last factor ended. From that point,
 * we copy all remaining factors (they represent new coverage of the text).
 * 
 * @param output_path Path to the output file (already contains thread 0's factors)
 * @param contexts Vector of all thread contexts (for accessing temp file paths)
 * @return size_t Total number of factors in the final output
 * @throws std::runtime_error If file I/O operations fail
 */
size_t ParallelFactorizer::merge_temp_files(const std::string& output_path,
                                         std::vector<ThreadContext>& contexts) {
    // Thread 0 already wrote to output_path, so we open in append mode
    // to add factors from other threads
    
    size_t total_factors = 0;
    size_t current_position = 0;  // Track the current end position in the merged output
    size_t text_length = contexts.empty() ? 0 : contexts[0].text_length;
    std::optional<Factor> last_written_factor;
    
    // First, read thread 0's factors to find the last one and count them
    {
        std::ifstream ifs(output_path, std::ios::binary);
        if (!ifs) {
            throw std::runtime_error("Cannot open output file for reading: " + output_path);
        }
        
        Factor factor;
        while (ifs.read(reinterpret_cast<char*>(&factor), sizeof(Factor))) {
            total_factors++;
            last_written_factor = factor;
            
            size_t factor_end = factor.start + factor.length;
            if (factor_end > current_position) {
                current_position = factor_end;
            }
        }
    }
    
    // Process remaining threads if there are multiple threads and sequence isn't complete
    if (contexts.size() > 1 && current_position < text_length) {
        // Open output file in append mode for adding factors from other threads
        std::ofstream ofs(output_path, std::ios::binary | std::ios::app);
        if (!ofs) {
            throw std::runtime_error("Cannot open output file for appending: " + output_path);
        }
        
        // Process remaining threads (i >= 1)
        for (size_t i = 1; i < contexts.size(); i++) {
            // Check if we've already covered the entire sequence
            if (current_position >= text_length) {
                break;  // Stop merging - sequence is complete
            }
            
            std::ifstream ifs(contexts[i].temp_file_path, std::ios::binary);
            if (!ifs) continue;
            
            Factor factor;
            bool found_convergence = false;
            
            // For subsequent threads, skip factors until we find convergence point
            // Convergence: last_written_factor.end == current_factor.start
            
            if (!last_written_factor.has_value()) {
                // No last written factor - shouldn't happen, but handle gracefully
                continue;
            }
            
            size_t last_end = last_written_factor->start + last_written_factor->length;
            
            // Read factors and look for convergence
            while (ifs.read(reinterpret_cast<char*>(&factor), sizeof(Factor))) {
                if (!found_convergence) {
                    // Still looking for convergence point
                    if (factor.start == last_end) {
                        // Found convergence! Start copying from this factor onwards
                        found_convergence = true;
                    } else {
                        // Skip this factor - it's before convergence
                        continue;
                    }
                }
                
                // Write factor to output (either we found convergence or we're continuing)
                ofs.write(reinterpret_cast<const char*>(&factor), sizeof(Factor));
                total_factors++;
                last_written_factor = factor;
                
                size_t factor_end = factor.start + factor.length;
                if (factor_end > current_position) {
                    current_position = factor_end;
                }
                last_end = factor_end;  // Update for next iteration
                
                // Stop if we've covered the entire sequence
                if (current_position >= text_length) {
                    break;
                }
            }
            
            // If we've covered the entire sequence, stop processing more threads
            if (current_position >= text_length) {
                break;
            }
        }
    }
    
    // Write footer at the end (v2 format) - single location for all cases
    std::ofstream ofs(output_path, std::ios::binary | std::ios::app);
    if (!ofs) {
        throw std::runtime_error("Cannot open output file for appending footer: " + output_path);
    }
    
    FactorFileFooter footer;
    footer.num_factors = total_factors;
    footer.num_sequences = 0;  // Unknown for general (non-FASTA) factorization
    footer.num_sentinels = 0;  // No sentinels for general factorization
    footer.footer_size = sizeof(FactorFileFooter);
    
    ofs.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
    
    return total_factors;
}

/**
 * @brief Removes temporary files created by worker threads
 * 
 * Deletes temporary factor files created by threads 1 and higher. Thread 0's "file"
 * is actually the final output file, so it's skipped during cleanup.
 * 
 * Errors during cleanup are logged as warnings but don't cause failure, since the
 * factorization has already completed successfully.
 * 
 * @param contexts Vector of all thread contexts (containing temp file paths)
 */
void ParallelFactorizer::cleanup_temp_files(const std::vector<ThreadContext>& contexts) {
    for (const auto& ctx : contexts) {
        // Skip thread 0 - it wrote directly to the output file
        if (ctx.thread_id == 0) continue;
        
        try {
            fs::remove(ctx.temp_file_path);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to remove temporary file: " << ctx.temp_file_path
                     << " - " << e.what() << std::endl;
        }
    }
}

/**
 * @brief Performs parallel factorization on a file
 * 
 * Convenience wrapper that reads an entire file into memory and then factorizes it
 * using the in-memory parallel_factorize() function.
 * 
 * NOTE: For very large files, this may consume significant memory. Consider using
 * streaming approaches for files larger than available RAM.
 * 
 * @param input_path Path to the input text file
 * @param output_path Path where the binary factor file will be written
 * @param num_threads Number of threads to use (0 = auto-detect)
 * @return size_t Total number of factors produced
 * @throws std::runtime_error If the input file cannot be opened
 */
size_t ParallelFactorizer::parallel_factorize_file(const std::string& input_path, 
                                                const std::string& output_path,
                                                size_t num_threads) {
    // Read the file content
    std::ifstream is(input_path, std::ios::binary);
    if (!is) {
        throw std::runtime_error("Cannot open input file: " + input_path);
    }
    
    std::string data((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    return parallel_factorize(data, output_path, num_threads);
}

/**
 * @brief Performs parallel DNA factorization with reverse complement support
 * 
 * PLACEHOLDER: This function is not yet implemented.
 * 
 * When implemented, this will perform parallel factorization on DNA sequences with
 * reverse complement awareness, allowing factors to reference either the forward
 * strand or the reverse complement. The MSB of the reference field will indicate
 * whether the match is to the reverse complement.
 * 
 * @param text Input DNA text to factorize (should contain only A, C, G, T)
 * @param output_path Path where the binary factor file will be written
 * @param num_threads Number of threads to use (0 = auto-detect)
 * @return size_t Total number of factors produced (currently returns 0)
 */
size_t ParallelFactorizer::parallel_factorize_dna_w_rc(std::string_view text, 
                                                    const std::string& output_path,
                                                    size_t num_threads) {
    // Implementation would be similar to parallel_factorize but using DNA-specific
    // algorithms from the original code - this is a placeholder
    // In a complete implementation, we'd use the nolzss_dna_w_rc algorithm
    std::cerr << "DNA parallel factorization not yet implemented" << std::endl;
    return 0;
}

/**
 * @brief Performs parallel DNA factorization with reverse complement on a file
 * 
 * PLACEHOLDER: This function is not yet implemented.
 * 
 * File-based wrapper for parallel_factorize_dna_w_rc(). Reads the entire file
 * into memory and processes it.
 * 
 * @param input_path Path to the input DNA file
 * @param output_path Path where the binary factor file will be written
 * @param num_threads Number of threads to use (0 = auto-detect)
 * @return size_t Total number of factors produced (currently returns 0)
 * @throws std::runtime_error If the input file cannot be opened
 */
size_t ParallelFactorizer::parallel_factorize_file_dna_w_rc(const std::string& input_path, 
                                                         const std::string& output_path,
                                                         size_t num_threads) {
    std::ifstream is(input_path, std::ios::binary);
    if (!is) {
        throw std::runtime_error("Cannot open input file: " + input_path);
    }
    
    std::string data((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    return parallel_factorize_dna_w_rc(data, output_path, num_threads);
}

} // namespace noLZSS
