#pragma once
#include "factorizer.hpp"
#include <sdsl/suffix_trees.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <optional>

namespace noLZSS {

// Define cst_t for use in parallel factorization
using cst_t = sdsl::cst_sada<>;

/**
 * @brief Thread context for parallel factorization
 * 
 * Contains all the information needed for a single thread
 * in the parallel factorization process.
 */
struct ThreadContext {
    size_t thread_id;               // Thread identifier
    size_t start_pos;               // Starting position in text
    size_t end_pos;                 // Position to stop (start of next thread)
    size_t text_length;             // Total text length
    std::string temp_file_path;     // Path to temporary file for this thread
    bool is_last_thread;            // Flag indicating if this is the last thread
    
    // Convergence tracking state
    size_t next_thread_factor_index = 0;  // Index of next factor to read from next thread
    std::optional<Factor> last_read_factor;  // Last factor read from next thread's file
};

/**
 * @brief Parallel implementation of noLZSS factorization
 * 
 * This class implements parallel factorization for noLZSS by dividing
 * the input text among multiple threads and using convergence detection
 * to efficiently merge results.
 */
class ParallelFactorizer {
public:
    /**
     * @brief Minimum number of characters per thread to make parallelization worthwhile
     * 
     * If the input text is smaller than this value multiplied by the requested thread count,
     * fewer threads will be used. This prevents excessive overhead from thread management
     * when the per-thread workload would be too small.
     * 
     * Can be adjusted based on hardware characteristics and typical input sizes.
     */
    static constexpr size_t MIN_CHARS_PER_THREAD = 100000;

    /**
     * @brief Constructor
     */
    ParallelFactorizer() = default;
    
    /**
     * @brief Destructor
     */
    ~ParallelFactorizer() = default;

    /**
     * @brief Main parallel factorization function
     * 
     * Factorizes the input text using multiple threads and writes results to output file
     * 
     * @param text Input text to factorize
     * @param output_path Path to output binary factor file
     * @param num_threads Number of threads to use (0 for auto-detection)
     * @param start_pos Starting position in the text for factorization (default: 0)
     * @return Number of factors produced
     */
    size_t parallel_factorize(std::string_view text, const std::string& output_path, 
                             size_t num_threads = 0, size_t start_pos = 0);
    
    /**
     * @brief File-based parallel factorization
     * 
     * Reads input from a file, factorizes using multiple threads, writes results to output file
     * 
     * @param input_path Path to input text file
     * @param output_path Path to output binary factor file
     * @param num_threads Number of threads to use (0 for auto-detection)
     * @param start_pos Starting position in the text for factorization (default: 0)
     * @return Number of factors produced
     */
    size_t parallel_factorize_file(const std::string& input_path, const std::string& output_path,
                                  size_t num_threads = 0, size_t start_pos = 0);
    
    /**
     * @brief DNA-specific parallel factorization with reverse complement support
     * 
     * Takes raw DNA text, prepares it with reverse complement, and calls the core
     * parallel_factorize_multiple_dna_w_rc() function.
     * 
     * @param text Input DNA text (raw nucleotides)
     * @param output_path Path to output binary factor file
     * @param num_threads Number of threads to use (0 for auto-detection)
     * @param start_pos Starting position in the text for factorization (default: 0)
     * @return Number of factors produced
     */
    size_t parallel_factorize_dna_w_rc(std::string_view text, const std::string& output_path,
                                      size_t num_threads = 0, size_t start_pos = 0);
    
    /**
     * @brief Core parallel DNA factorization with reverse complement for prepared strings
     * 
     * Performs parallel DNA factorization on already-prepared strings that include
     * reverse complement. This is the core function that parallel_factorize_dna_w_rc()
     * calls after preparing the input, and can be used directly by FASTA processors
     * that have already prepared their strings.
     * 
     * The input string S should be in the format: T + sentinel + RC(T) + sentinel
     * where T is the original sequence(s) concatenated with sentinels.
     * 
     * @param prepared_string The prepared string with reverse complement
     * @param original_length Length of the original sequence (before RC, excluding final sentinel)
     * @param output_path Path to output binary factor file
     * @param num_threads Number of threads to use (0 for auto-detection)
     * @param start_pos Starting position in the original sequence for factorization (default: 0)
     * @return Number of factors produced
     */
    size_t parallel_factorize_multiple_dna_w_rc(const std::string& prepared_string,
                                                size_t original_length,
                                                const std::string& output_path,
                                                size_t num_threads = 0,
                                                size_t start_pos = 0);
    
    /**
     * @brief File-based DNA-specific parallel factorization with reverse complement support
     * 
     * @param input_path Path to input DNA text file
     * @param output_path Path to output binary factor file
     * @param num_threads Number of threads to use (0 for auto-detection)
     * @return Number of factors produced
     */
    size_t parallel_factorize_file_dna_w_rc(const std::string& input_path, const std::string& output_path,
                                           size_t num_threads = 0);
    
    /**
     * @brief Create temporary file path (public for use by fasta processor)
     * 
     * @param thread_id Thread identifier
     * @return Unique temporary file path
     */
    std::string create_temp_file_path(size_t thread_id);
    
    /**
     * @brief Merge temporary files into final output (public for use by fasta processor)
     * 
     * @param output_path Path to final output file
     * @param contexts Thread contexts containing temporary file info
     * @return Number of factors in final output
     */
    size_t merge_temp_files(const std::string& output_path,
                           std::vector<ThreadContext>& contexts);
    
    /**
     * @brief Clean up temporary files (public for use by fasta processor)
     * 
     * @param contexts Thread contexts containing temporary file paths
     */
    void cleanup_temp_files(const std::vector<ThreadContext>& contexts);
    
    /**
     * @brief Thread worker function for DNA factorization with reverse complement (public for use by fasta processor)
     * 
     * Worker thread that performs DNA-specific factorization with reverse complement
     * awareness. Similar to factorize_thread() but uses the DNA w/ RC algorithm.
     * 
     * @param cst The compressed suffix tree (shared by all threads)
     * @param rmqF RMQ for forward starts
     * @param rmqRcEnd RMQ for reverse complement ends
     * @param fwd_starts Forward start positions vector
     * @param rc_ends Reverse complement end positions vector
     * @param INF Infinity value for invalid positions
     * @param N Original sequence length
     * @param ctx Thread context
     * @param all_contexts All thread contexts (for convergence checking)
     * @param file_mutexes Mutexes for protecting file access
     */
    void factorize_dna_w_rc_thread(const cst_t& cst,
                                   const sdsl::rmq_succinct_sct<>& rmqF,
                                   const sdsl::rmq_succinct_sct<>& rmqRcEnd,
                                   const sdsl::int_vector<64>& fwd_starts,
                                   const sdsl::int_vector<64>& rc_ends,
                                   uint64_t INF,
                                   size_t N,
                                   ThreadContext& ctx,
                                   std::vector<ThreadContext>& all_contexts,
                                   std::vector<std::mutex>& file_mutexes);

private:
    /**
     * @brief Thread worker function that performs factorization
     * 
     * @param cst The compressed suffix tree (shared by all threads)
     * @param rmq RMQ support for the suffix tree (shared by all threads)
     * @param ctx Thread context
     * @param all_contexts All thread contexts (for accessing next thread's info)
     * @param file_mutexes Mutexes for protecting file access
     */
    void factorize_thread(const cst_t& cst, const sdsl::rmq_succinct_sct<>& rmq,
                         ThreadContext& ctx, 
                         std::vector<ThreadContext>& all_contexts,
                         std::vector<std::mutex>& file_mutexes);
    
    /**
     * @brief Read a single factor from a file at a specific index
     * 
     * @param file_path Path to the factor file
     * @param factor_index Index of the factor to read (0-based)
     * @return Factor if successful, std::nullopt if index is out of bounds
     */
    std::optional<Factor> read_factor_at_index(const std::string& file_path, size_t factor_index);
    
    /**
     * @brief Check if current thread has converged with next thread
     * 
     * Checks if the current factor's end position matches the start of any
     * factor in the next thread, indicating convergence. Uses stateful reading
     * to efficiently read factors one at a time.
     * 
     * @param current_end Current factorization end position
     * @param next_ctx Next thread context (contains state tracking)
     * @param next_file_mutex Mutex for next thread's file
     * @return true if convergence detected, false otherwise
     */
    bool check_convergence(size_t current_end, ThreadContext& next_ctx,
                          std::mutex& next_file_mutex);
    
    /**
     * @brief Thread-safe factor writing to temporary file
     * 
     * @param factor Factor to write
     * @param file_path Path to temporary file
     * @param file_mutex Mutex for protecting file access
     */
    void write_factor(const Factor& factor, const std::string& file_path, 
                      std::mutex& file_mutex);

    /**
     * @brief Read a specific factor from a temporary file
     * 
     * @param file_path Path to temporary file
     * @param index Index of the factor to read
     * @param file_mutex Mutex for protecting file access
     * @return Factor at the specified index or nullopt if not found
     */
    std::optional<Factor> read_factor_at(const std::string& file_path, 
                                      size_t index, 
                                      std::mutex& file_mutex);
};

} // namespace noLZSS
