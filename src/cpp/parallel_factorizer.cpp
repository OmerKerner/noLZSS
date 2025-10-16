#include "parallel_factorizer.hpp"
#include "factorizer_helpers.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <chrono>

namespace fs = std::filesystem;
namespace noLZSS {
// Helper functions lcp() and next_leaf() are now in factorizer_helpers.hpp

std::string ParallelFactorizer::create_temp_file_path(size_t thread_id) {
    auto timestamp = std::chrono::system_clock::now().time_since_epoch().count();
    std::string temp_dir = fs::temp_directory_path().string();
    return temp_dir + "/noLZSS_temp_" + std::to_string(timestamp) + "_" + 
           std::to_string(thread_id) + ".bin";
}

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
        contexts[i].temp_file_path = create_temp_file_path(i);
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

void ParallelFactorizer::write_factor(const Factor& factor, const std::string& file_path, 
                                   std::mutex& file_mutex) {
    std::lock_guard<std::mutex> lock(file_mutex);
    
    std::ofstream ofs(file_path, std::ios::binary | std::ios::app);
    if (!ofs) {
        throw std::runtime_error("Cannot open temporary file for writing: " + file_path);
    }
    
    ofs.write(reinterpret_cast<const char*>(&factor), sizeof(Factor));
}

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

size_t ParallelFactorizer::merge_temp_files(const std::string& output_path,
                                         std::vector<ThreadContext>& contexts) {
    // Open output file
    std::ofstream ofs(output_path, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Cannot open output file for writing: " + output_path);
    }
    
    // Reserve space for the header, we'll write it at the end
    FactorFileHeader header;
    ofs.seekp(sizeof(FactorFileHeader));
    
    size_t total_factors = 0;
    size_t current_position = 0;  // Track the current end position in the merged output
    size_t text_length = contexts.empty() ? 0 : contexts[0].text_length;
    std::optional<Factor> last_written_factor;
    
    // Process each thread's temp file
    for (size_t i = 0; i < contexts.size(); i++) {
        // Check if we've already covered the entire sequence
        if (current_position >= text_length) {
            break;  // Stop merging - sequence is complete
        }
        
        std::ifstream ifs(contexts[i].temp_file_path, std::ios::binary);
        if (!ifs) continue;
        
        Factor factor;
        bool found_convergence = false;
        
        // For first thread (i == 0), copy all factors until convergence or end
        if (i == 0) {
            // Copy all factors directly to output
            while (ifs.read(reinterpret_cast<char*>(&factor), sizeof(Factor))) {
                ofs.write(reinterpret_cast<const char*>(&factor), sizeof(Factor));
                total_factors++;
                last_written_factor = factor;
                
                size_t factor_end = factor.start + factor.length;
                if (factor_end > current_position) {
                    current_position = factor_end;
                }
                
                // Stop if we've covered the entire sequence
                if (current_position >= text_length) {
                    break;
                }
            }
        } else {
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
        }
        
        // If we've covered the entire sequence, stop processing more threads
        if (current_position >= text_length) {
            break;
        }
    }
    
    // Update and write header
    header.num_factors = total_factors;
    header.num_sequences = 1;  // Default for general factorization
    header.num_sentinels = 0;
    header.header_size = sizeof(FactorFileHeader);
    
    // Go back and write the header
    ofs.seekp(0);
    ofs.write(reinterpret_cast<const char*>(&header), sizeof(header));
    
    return total_factors;
}

void ParallelFactorizer::cleanup_temp_files(const std::vector<ThreadContext>& contexts) {
    for (const auto& ctx : contexts) {
        try {
            fs::remove(ctx.temp_file_path);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to remove temporary file: " << ctx.temp_file_path
                     << " - " << e.what() << std::endl;
        }
    }
}

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

size_t ParallelFactorizer::parallel_factorize_dna_w_rc(std::string_view text, 
                                                    const std::string& output_path,
                                                    size_t num_threads) {
    // Implementation would be similar to parallel_factorize but using DNA-specific
    // algorithms from the original code - this is a placeholder
    // In a complete implementation, we'd use the nolzss_dna_w_rc algorithm
    std::cerr << "DNA parallel factorization not yet implemented" << std::endl;
    return 0;
}

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
