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
        num_threads = std::thread::hardware_concurrency();
        // Limit threads for small inputs
        if (text.length() < 100000) {
            num_threads = std::max(1UL, num_threads / 2);
        }
    }
    
    // Ensure at least one thread
    num_threads = std::max(1UL, num_threads);
    // Limit threads if text is too small
    num_threads = std::min(num_threads, text.length() / 10000);
    if (num_threads == 0) num_threads = 1;
    
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
        
        // Write the factor to our thread's temporary file
        write_factor(current_factor, ctx.temp_file_path, file_mutexes[ctx.thread_id]);
        
        // Check for convergence if we're in the next thread's region
        if (next_ctx && lambda_sufnum >= next_ctx->start_pos) {
            size_t current_end = lambda_sufnum + l;
            if (check_convergence(current_end, *next_ctx, file_mutexes[next_ctx->thread_id])) {
                // We've converged, stop this thread's factorization
                break;
            }
        }
        
        // Advance to next position
        lambda = next_leaf(cst, lambda, l);
        lambda_node_depth = cst.node_depth(lambda);
        lambda_sufnum = cst.sn(lambda);
    }
}

bool ParallelFactorizer::check_convergence(size_t current_end, ThreadContext& next_ctx,
                                        std::mutex& next_file_mutex) {
    // Read all factors from the next thread's file
    auto next_factors = read_factors(next_ctx.temp_file_path, next_file_mutex);
    
    // Look for a factor that starts exactly where our current factorization ends
    for (const auto& factor : next_factors) {
        if (factor.start == current_end) {
            return true; // Found convergence
        }
    }
    
    return false;
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

std::vector<Factor> ParallelFactorizer::read_factors(const std::string& file_path, 
                                                 std::mutex& file_mutex) {
    std::lock_guard<std::mutex> lock(file_mutex);
    
    std::vector<Factor> factors;
    std::ifstream ifs(file_path, std::ios::binary);
    
    if (!ifs) {
        return factors; // Return empty vector if file can't be opened
    }
    
    Factor factor;
    while (ifs.read(reinterpret_cast<char*>(&factor), sizeof(Factor))) {
        factors.push_back(factor);
    }
    
    return factors;
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
    std::vector<size_t> convergence_points(contexts.size());
    
    // Find the convergence points between threads
    for (size_t i = 0; i < contexts.size() - 1; i++) {
        std::vector<Factor> current_factors;
        std::vector<Factor> next_factors;
        
        // Read all factors from current and next thread
        {
            std::mutex dummy_mutex;
            current_factors = read_factors(contexts[i].temp_file_path, dummy_mutex);
            next_factors = read_factors(contexts[i+1].temp_file_path, dummy_mutex);
        }
        
        if (current_factors.empty()) {
            convergence_points[i] = 0;
            continue;
        }
        
        // Default convergence point is the start of next thread
        convergence_points[i] = contexts[i+1].start_pos;
        
        // Look for the actual convergence point
        for (size_t j = 0; j < current_factors.size(); j++) {
            const auto& current_factor = current_factors[j];
            size_t current_end = current_factor.start + current_factor.length;
            
            // If this factor ends beyond the next thread's start, check for convergence
            if (current_factor.start >= contexts[i+1].start_pos) {
                for (const auto& next_factor : next_factors) {
                    if (next_factor.start == current_end) {
                        // Found convergence point - use the position of this factor
                        convergence_points[i] = current_factor.start;
                        break;
                    }
                }
                
                // If we found a convergence point, break out of the loop
                if (convergence_points[i] != contexts[i+1].start_pos) {
                    break;
                }
            }
        }
    }
    
    // Set the convergence point for the last thread to the end of text
    if (!contexts.empty()) {
        convergence_points[contexts.size() - 1] = contexts.back().text_length;
    }
    
    // Now merge the files, only including factors up to the convergence point
    for (size_t i = 0; i < contexts.size(); i++) {
        std::ifstream ifs(contexts[i].temp_file_path, std::ios::binary);
        if (!ifs) continue;
        
        Factor factor;
        while (ifs.read(reinterpret_cast<char*>(&factor), sizeof(Factor))) {
            // Only include factors that start before the convergence point
            if (factor.start < convergence_points[i]) {
                ofs.write(reinterpret_cast<const char*>(&factor), sizeof(Factor));
                total_factors++;
            }
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
