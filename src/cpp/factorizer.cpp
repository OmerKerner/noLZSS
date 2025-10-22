#include "factorizer.hpp"
#include "factorizer_core.hpp"  // Template implementations in detail:: namespace
#include "factorizer_helpers.hpp"
#include "parallel_factorizer.hpp"
#include <sdsl/suffix_trees.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <optional>
#include <limits>

namespace noLZSS {
// Helper functions lcp() and next_leaf() are now in factorizer_helpers.hpp

// ---------- genomic utilities ----------
char complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:
            // Handle invalid input, e.g., throw an exception or return a sentinel value
            throw std::invalid_argument("Invalid nucleotide: " + std::string(1, c));
    }
}
static std::string revcomp(std::string_view s) {
    std::string r; r.resize(s.size());
    for (size_t i = 0, n = s.size(); i < n; ++i) r[i] = complement(s[n-1-i]);
    return r;
}

/**
 * @brief Prepares multiple DNA sequences for factorization with reverse complement and tracks sentinel positions.
 *
 * Takes multiple DNA sequences, concatenates them with unique sentinels, appends
 * their reverse complements with unique sentinels, and tracks sentinel positions.
 * The output format is compatible with nolzss_multiple_dna_w_rc(): S = T1!T2@T3$rt(T3)%rt(T2)^rt(T1)&
 *
 * @param sequences Vector of DNA sequence strings (should contain only A, C, T, G)
 * @return PreparedSequenceResult containing:
 *         - prepared_string: The formatted string with sequences and reverse complements
 *         - original_length: Length of the original sequences part (before reverse complements)
 *         - sentinel_positions: Positions of all sentinels in the prepared string
 *
 * @throws std::invalid_argument If too many sequences (>125) or invalid nucleotides found
 * @throws std::runtime_error If sequences contain invalid characters
 *
 * @note Sentinels avoid 0, A(65), C(67), G(71), T(84) - lowercase nucleotides are safe as sentinels
 * @note The function validates that all sequences contain only valid DNA nucleotides
 * @note Input sequences can be lowercase or uppercase, output is always uppercase
 */
PreparedSequenceResult prepare_multiple_dna_sequences_w_rc(const std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return {"", 0, {}};
    }
    
    // Check if we have too many sequences
    if (sequences.size() > 125) {
        throw std::invalid_argument("Too many sequences: maximum 125 sequences supported (due to sentinel character limitations)");
    }
    
    // Validate sequences contain only valid DNA nucleotides
    for (size_t i = 0; i < sequences.size(); ++i) {
        for (char c : sequences[i]) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && 
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                throw std::runtime_error("Invalid nucleotide '" + std::string(1, c) + 
                                       "' found in sequence " + std::to_string(i));
            }
        }
    }
    
    PreparedSequenceResult result;
    
    // Calculate total size for reservation
    size_t total_size = 0;
    for (const auto& seq : sequences) {
        total_size += seq.length() * 2; // Original + reverse complement
    }
    total_size += sequences.size() * 2; // Add space for sentinels
    result.prepared_string.reserve(total_size);
    
    // Generate sentinel characters avoiding 0 and uppercase DNA nucleotides A(65), C(67), G(71), T(84)
    auto get_sentinel = [](size_t index) -> char {
        char sentinel = 1;
        size_t count = 0;
        
        while (true) {
            // Check if current sentinel is valid (not 0 and not uppercase DNA nucleotides)
            if (sentinel != 0 && sentinel != 'A' && sentinel != 'C' && sentinel != 'G' && sentinel != 'T') {
                if (count == index) {
                    return sentinel;  // Found the sentinel for this index
                }
                count++;
            }
            sentinel++;
            if (sentinel == 0) sentinel = 1; // wrap around, skip 0
        }
    };
    
    // First, add original sequences with sentinels
    for (size_t i = 0; i < sequences.size(); ++i) {
        // Convert to uppercase
        std::string seq = sequences[i];
        for (char& c : seq) {
            if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        }
        result.prepared_string += seq;
        
        // Add sentinel and track its position
        size_t sentinel_pos = result.prepared_string.length();
        char sentinel = get_sentinel(i);
        result.prepared_string += sentinel;
        result.sentinel_positions.push_back(sentinel_pos);
    }
    
    result.original_length = result.prepared_string.length();
    
    // Then, add reverse complements with different sentinels
    for (int i = static_cast<int>(sequences.size()) - 1; i >= 0; --i) {
        // Convert to uppercase first
        std::string seq = sequences[i];
        for (char& c : seq) {
            if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        }
        
        // Add reverse complement
        std::string rc = revcomp(seq);
        result.prepared_string += rc;
        
        // Add sentinel and track its position (offset by sequences.size() to make them unique)
        size_t sentinel_pos = result.prepared_string.length();
        char sentinel = get_sentinel(sequences.size() + (sequences.size() - 1 - i));
        result.prepared_string += sentinel;
        result.sentinel_positions.push_back(sentinel_pos);
    }
    
    return result;
}

/**
 * @brief Prepares multiple DNA sequences for factorization without reverse complement and tracks sentinel positions.
 *
 * Takes multiple DNA sequences, concatenates them with unique sentinels, and tracks sentinel positions.
 * Unlike prepare_multiple_dna_sequences_w_rc(), this function does not append
 * reverse complements. The output format is: S = T1!T2@T3$
 *
 * @param sequences Vector of DNA sequence strings (should contain only A, C, T, G)
 * @return PreparedSequenceResult containing:
 *         - prepared_string: The formatted string with sequences and sentinels
 *         - original_length: Total length of the concatenated string (same as prepared_string.length())
 *         - sentinel_positions: Positions of all sentinels in the prepared string
 *
 * @throws std::invalid_argument If too many sequences (>250) or invalid nucleotides found
 * @throws std::runtime_error If sequences contain invalid characters
 *
 * @note Sentinels range from 1-251, avoiding 0, A(65), C(67), G(71), T(84)
 * @note The function validates that all sequences contain only valid DNA nucleotides
 * @note Input sequences can be lowercase or uppercase, output is always uppercase
 */
PreparedSequenceResult prepare_multiple_dna_sequences_no_rc(const std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return {"", 0, {}};
    }
    
    // Check if we have too many sequences
    if (sequences.size() > 250) {
        throw std::invalid_argument("Too many sequences: maximum 250 sequences supported (due to sentinel character limitations)");
    }
    
    // Validate sequences contain only valid DNA nucleotides
    for (size_t i = 0; i < sequences.size(); ++i) {
        for (char c : sequences[i]) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && 
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                throw std::runtime_error("Invalid nucleotide '" + std::string(1, c) + 
                                       "' found in sequence " + std::to_string(i));
            }
        }
    }
    
    PreparedSequenceResult result;
    
    // Calculate total size for reservation
    size_t total_size = 0;
    for (const auto& seq : sequences) {
        total_size += seq.length();
    }
    total_size += sequences.size(); // Add space for sentinels
    result.prepared_string.reserve(total_size);
    
    // Generate sentinel characters avoiding 0 and uppercase DNA nucleotides A(65), C(67), G(71), T(84)
    auto get_sentinel = [](size_t index) -> char {
        char sentinel = 1;
        size_t count = 0;
        
        while (true) {
            // Check if current sentinel is valid (not 0 and not uppercase DNA nucleotides)
            if (sentinel != 0 && sentinel != 'A' && sentinel != 'C' && sentinel != 'G' && sentinel != 'T') {
                if (count == index) {
                    return sentinel;  // Found the sentinel for this index
                }
                count++;
            }
            sentinel++;
            if (sentinel == 0) sentinel = 1; // wrap around, skip 0
        }
    };
    
    // Add original sequences with sentinels
    for (size_t i = 0; i < sequences.size(); ++i) {
        // Convert to uppercase
        std::string seq = sequences[i];
        for (char& c : seq) {
            if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        }
        result.prepared_string += seq;
        
        // Add sentinel and track its position (only between sequences, not after the last one)
        if (i < sequences.size() - 1) {
            size_t sentinel_pos = result.prepared_string.length();
            char sentinel = get_sentinel(i);
            result.prepared_string += sentinel;
            result.sentinel_positions.push_back(sentinel_pos);
        }
    }
    
    result.original_length = result.prepared_string.length();
    
    return result;
}

// ---------- Core template algorithms ----------
// Template implementations are in factorizer_core.hpp in the detail:: namespace
// This allows them to be used from any compilation unit without code duplication

// ------------- public wrappers -------------

/**
 * @brief Factorizes text from a file using the noLZSS algorithm.
 *
 * This template function reads text directly from a file and performs factorization
 * without loading the entire file into memory. This is more memory-efficient for
 * large files.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param path Path to input file containing text
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This function builds the suffix tree directly from the file
 * @see factorize_file() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_file_stream(const std::string& path, Sink&& sink, size_t start_pos) {
    // sdsl-lite will automatically add the sentinel when needed
    cst_t cst; construct(cst, path, 1);
    return detail::nolzss(cst, std::forward<Sink>(sink), start_pos);
}

/**
 * @brief Counts noLZSS factors in a text string.
 *
 * This function provides a convenient way to count factors without storing them.
 * It uses the sink-based factorization internally with a counting lambda.
 *
 * @param text Input text string
 * @return Number of factors in the factorization
 *
 * @note This is more memory-efficient than factorize() when you only need the count
 * @see factorize() for getting the actual factors
 * @see count_factors_file() for file-based counting
 */
size_t count_factors(std::string_view text, size_t start_pos) {
    size_t n = 0;
    std::string tmp(text);
    cst_t cst; construct_im(cst, tmp, 1);
    detail::nolzss(cst, [&](const Factor&){ ++n; }, start_pos);
    return n;
}

/**
 * @brief Counts noLZSS factors in a file.
 *
 * This function reads text from a file and counts factors without storing them
 * or loading the entire file into memory. It's the most memory-efficient way
 * to get factor counts for large files.
 *
 * @param path Path to input file containing text
 * @return Number of factors in the factorization
 *
 * @note This function builds the suffix tree directly from the file
 * @see count_factors() for in-memory counting
 * @see factorize_file() for getting the actual factors from a file
 */
size_t count_factors_file(const std::string& path, size_t start_pos) {
    size_t n = 0;
    factorize_file_stream(path, [&](const Factor&){ ++n; }, start_pos);
    return n;
}

/**
 * @brief Factorizes a text string and returns factors as a vector.
 *
 * This is the main user-facing function for in-memory factorization.
 * It performs noLZSS factorization and returns all factors in a vector.
 *
 * @param text Input text string
 * @return Vector containing all factors from the factorization
 *
 * @note Factors are returned in order of appearance in the text
 * @note The returned factors are non-overlapping and cover the entire input
 * @see factorize_file() for file-based factorization
 */
std::vector<Factor> factorize(std::string_view text, size_t start_pos) {
    std::vector<Factor> out;
    std::string tmp(text);
    cst_t cst; construct_im(cst, tmp, 1);
    detail::nolzss(cst, [&](const Factor& f){ out.push_back(f); }, start_pos);
    return out;
}

/**
 * @brief Factorizes text from a file and returns factors as a vector.
 *
 * This function reads text from a file, performs factorization, and returns
 * all factors in a vector. The reserve_hint parameter can improve performance
 * when you have an estimate of the number of factors.
 *
 * @param path Path to input file containing text
 * @param reserve_hint Optional hint for reserving space in output vector (0 = no hint)
 * @return Vector containing all factors from the factorization
 *
 * @note Use reserve_hint for better performance when you know approximate factor count
 * @note This is more memory-efficient than factorize() for large files
 * @see factorize() for in-memory factorization
 */
std::vector<Factor> factorize_file(const std::string& path, size_t reserve_hint, size_t start_pos) {
    std::vector<Factor> out;
    if (reserve_hint) out.reserve(reserve_hint);
    factorize_file_stream(path, [&](const Factor& f){ out.push_back(f); }, start_pos);
    return out;
}

/**
 * @brief Writes noLZSS factors from a file to a binary output file.
 *
 * This function reads text from an input file, performs factorization, and
 * writes the resulting factors in binary format to an output file. Each factor
 * is written as two uint64_t values (start position, length).
 *
 * @param in_path Path to input file containing text
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref)
 * @note Footer is written at the END of the file
 * @note This function overwrites the output file if it exists
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_binary_file(const std::string& in_path, const std::string& out_path) {
    // Get input file size for total_length
    std::ifstream infile(in_path, std::ios::binary | std::ios::ate);
    if (!infile) {
        throw std::runtime_error("Cannot open input file: " + in_path);
    }
    uint64_t total_length = infile.tellg();
    infile.close();
    
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Stream factors directly to file without collecting in memory
    size_t n = factorize_file_stream(in_path, [&](const Factor& f){
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    });
    
    // For non-FASTA files, create minimal footer with no sequence names or sentinels
    FactorFileFooter footer;
    footer.num_factors = n;
    footer.num_sequences = 0;  // Unknown for raw text files
    footer.num_sentinels = 0;  // No sentinels for raw text files
    footer.footer_size = sizeof(FactorFileFooter);
    footer.total_length = total_length;
    
    // Write footer at the end
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
    
    return n;
}

/**
 * @brief Factorizes a DNA text string with reverse complement awareness using the noLZSS algorithm.
 *
 * This is a template function that provides the core factorization functionality
 * for in-memory DNA text, considering both forward and reverse complement matches.
 * It builds a compressed suffix tree and applies the noLZSS algorithm to find all factors.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param text Input DNA text string
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This function copies the input string for suffix tree construction
 * @note For large inputs, consider using factorize_file_stream_dna_w_rc() instead
 * @see factorize_dna_w_rc() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_stream_dna_w_rc(std::string_view text, Sink&& sink) {
    std::string tmp(text);
    return detail::nolzss_dna_w_rc(tmp, std::forward<Sink>(sink));
}

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness using the noLZSS algorithm.
 *
 * This template function reads DNA text directly from a file and performs factorization
 * without loading the entire file into memory, considering both forward and reverse complement matches.
 * This is more memory-efficient for large files.
 *
 * @tparam Sink Callable type that accepts Factor objects
 * @param path Path to input file containing DNA text
 * @param sink Callable that receives each computed factor
 * @return Number of factors emitted
 *
 * @note This function builds the suffix tree directly from the file
 * @see factorize_file_dna_w_rc() for the non-template version that returns a vector
 */
template<class Sink>
size_t factorize_file_stream_dna_w_rc(const std::string& path, Sink&& sink) {
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("Cannot open input file: " + path);
    std::string data((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    return detail::nolzss_dna_w_rc(data, std::forward<Sink>(sink));
}

/**
 * @brief Factorizes a DNA text string with reverse complement awareness and returns factors as a vector.
 *
 * This is the main user-facing function for in-memory DNA factorization with reverse complement.
 * It performs noLZSS factorization and returns all factors in a vector.
 *
 * @param text Input DNA text string
 * @return Vector containing all factors from the factorization
 *
 * @note Factors are returned in order of appearance in the text
 * @note The returned factors are non-overlapping and cover the entire input
 * @see factorize_file_dna_w_rc() for file-based factorization
 */
std::vector<Factor> factorize_dna_w_rc(std::string_view text) {
    std::vector<Factor> out;
    factorize_stream_dna_w_rc(text, [&](const Factor& f){ out.push_back(f); });
    return out;
}

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness and returns factors as a vector.
 *
 * This function reads DNA text from a file, performs factorization with reverse complement, and returns
 * all factors in a vector. The reserve_hint parameter can improve performance
 * when you have an estimate of the number of factors.
 *
 * @param path Path to input file containing DNA text
 * @param reserve_hint Optional hint for reserving space in output vector (0 = no hint)
 * @return Vector containing all factors from the factorization
 *
 * @note Use reserve_hint for better performance when you know approximate factor count
 * @note This is more memory-efficient than factorize_dna_w_rc() for large files
 * @see factorize_dna_w_rc() for in-memory factorization
 */
std::vector<Factor> factorize_file_dna_w_rc(const std::string& path, size_t reserve_hint) {
    std::vector<Factor> out; if (reserve_hint) out.reserve(reserve_hint);
    factorize_file_stream_dna_w_rc(path, [&](const Factor& f){ out.push_back(f); });
    return out;
}

/**
 * @brief Counts noLZSS factors in a DNA text string with reverse complement awareness.
 *
 * This function provides a convenient way to count factors in DNA text without storing them.
 * It uses the sink-based factorization internally with a counting lambda.
 *
 * @param text Input DNA text string
 * @return Number of factors in the factorization
 *
 * @note This is more memory-efficient than factorize_dna_w_rc() when you only need the count
 * @see factorize_dna_w_rc() for getting the actual factors
 * @see count_factors_file_dna_w_rc() for file-based counting
 */
size_t count_factors_dna_w_rc(std::string_view text) {
    size_t n = 0; factorize_stream_dna_w_rc(text, [&](const Factor&){ ++n; }); return n;
}

/**
 * @brief Counts noLZSS factors in a DNA file with reverse complement awareness.
 *
 * This function reads DNA text from a file and counts factors without storing them
 * or loading the entire file into memory. It's the most memory-efficient way
 * to get factor counts for large DNA files.
 *
 * @param path Path to input file containing DNA text
 * @return Number of factors in the factorization
 *
 * @note This function builds the suffix tree directly from the file
 * @see count_factors_dna_w_rc() for in-memory counting
 * @see factorize_file_dna_w_rc() for getting the actual factors from a file
 */
size_t count_factors_file_dna_w_rc(const std::string& path) {
    size_t n = 0; factorize_file_stream_dna_w_rc(path, [&](const Factor&){ ++n; }); return n;
}

/**
 * @brief Writes noLZSS factors from a DNA file with reverse complement awareness to a binary output file.
 *
 * This function reads DNA text from an input file, performs factorization with reverse complement, and
 * writes the resulting factors in binary format to an output file. Each factor
 * is written as three uint64_t values (start position, length, ref).
 *
 * @param in_path Path to input file containing DNA text
 * @param out_path Path to output file where binary factors will be written
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref)
 * @note Footer is written at the END of the file
 * @note This function overwrites the output file if it exists
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_binary_file_dna_w_rc(const std::string& in_path, const std::string& out_path) {
    // Get input file size for total_length
    std::ifstream infile(in_path, std::ios::binary | std::ios::ate);
    if (!infile) {
        throw std::runtime_error("Cannot open input file: " + in_path);
    }
    uint64_t total_length = infile.tellg();
    infile.close();
    
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Stream factors directly to file without collecting in memory
    size_t n = factorize_file_stream_dna_w_rc(in_path, [&](const Factor& f){
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    });
    
    // Write empty sequence name (single null terminator)
    os.write("\0", 1);
    
    // For single DNA sequence files, create minimal footer
    FactorFileFooter footer;
    footer.num_factors = n;
    footer.num_sequences = 1;  // Single DNA sequence
    footer.num_sentinels = 0;  // No sentinels for single sequence
    footer.footer_size = sizeof(FactorFileFooter) + 1; // Empty sequence name
    footer.total_length = total_length;
    
    // Write footer at the end
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
    
    return n;
}

/**
 * @brief Factorizes a DNA text string with reverse complement awareness for multiple sequences and returns factors as a vector.
 *
 * This is the main user-facing function for in-memory DNA factorization with multiple sequences and reverse complement.
 * It performs noLZSS factorization and returns all factors in a vector.
 *
 * @param text Input DNA text string with multiple sequences and sentinels
 * @param start_pos Starting position for factorization (default: 0)
 * @return Vector containing all factors from the factorization
 *
 * @note Factors are returned in order of appearance in the text
 * @note The returned factors are non-overlapping and cover the entire input
 * @see factorize_file_multiple_dna_w_rc() for file-based factorization
 */
std::vector<Factor> factorize_multiple_dna_w_rc(std::string_view text, size_t start_pos) {
    std::vector<Factor> out;
    std::string tmp(text);
    detail::nolzss_multiple_dna_w_rc(tmp, [&](const Factor& f){ out.push_back(f); }, start_pos);
    return out;
}

/**
 * @brief Factorizes DNA text from a file with reverse complement awareness for multiple sequences and returns factors as a vector.
 *
 * This function reads DNA text from a file, performs factorization with reverse complement for multiple sequences, and returns
 * all factors in a vector. The reserve_hint parameter can improve performance
 * when you have an estimate of the number of factors.
 *
 * @param path Path to input file containing DNA text with multiple sequences
 * @param reserve_hint Optional hint for reserving space in output vector (0 = no hint)
 * @param start_pos Starting position for factorization (default: 0)
 * @return Vector containing all factors from the factorization
 *
 * @note Use reserve_hint for better performance when you know approximate factor count
 * @note This is more memory-efficient than factorize_multiple_dna_w_rc() for large files
 * @see factorize_multiple_dna_w_rc() for in-memory factorization
 */
std::vector<Factor> factorize_file_multiple_dna_w_rc(const std::string& path, size_t reserve_hint, size_t start_pos) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + path);
    }
    std::string text((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    
    std::vector<Factor> out; if (reserve_hint) out.reserve(reserve_hint);
    detail::nolzss_multiple_dna_w_rc(text, [&](const Factor& f){ out.push_back(f); }, start_pos);
    return out;
}

/**
 * @brief Counts noLZSS factors in a DNA text string with reverse complement awareness for multiple sequences.
 *
 * This function provides a convenient way to count factors in DNA text with multiple sequences without storing them.
 * It uses the sink-based factorization internally with a counting lambda.
 *
 * @param text Input DNA text string with multiple sequences
 * @param start_pos Starting position for factorization (default: 0)
 * @return Number of factors in the factorization
 *
 * @note This is more memory-efficient than factorize_multiple_dna_w_rc() when you only need the count
 * @see factorize_multiple_dna_w_rc() for getting the actual factors
 * @see count_factors_file_multiple_dna_w_rc() for file-based counting
 */
size_t count_factors_multiple_dna_w_rc(std::string_view text, size_t start_pos) {
    size_t n = 0; 
    std::string tmp(text);
    detail::nolzss_multiple_dna_w_rc(tmp, [&](const Factor&){ ++n; }, start_pos); 
    return n;
}

/**
 * @brief Counts noLZSS factors in a DNA file with reverse complement awareness for multiple sequences.
 *
 * This function reads DNA text from a file and counts factors without storing them
 * or loading the entire file into memory. It's the most memory-efficient way
 * to get factor counts for large DNA files with multiple sequences.
 *
 * @param path Path to input file containing DNA text with multiple sequences
 * @param start_pos Starting position for factorization (default: 0)
 * @return Number of factors in the factorization
 *
 * @note This function builds the suffix tree directly from the file
 * @see count_factors_multiple_dna_w_rc() for in-memory counting
 * @see factorize_file_multiple_dna_w_rc() for getting the actual factors from a file
 */
size_t count_factors_file_multiple_dna_w_rc(const std::string& path, size_t start_pos) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + path);
    }
    std::string text((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    
    size_t n = 0; 
    detail::nolzss_multiple_dna_w_rc(text, [&](const Factor&){ ++n; }, start_pos); 
    return n;
}

/**
 * @brief Writes noLZSS factors from a DNA file with reverse complement awareness for multiple sequences to a binary output file.
 *
 * This function reads DNA text from an input file, performs factorization with reverse complement for multiple sequences, and
 * writes the resulting factors in binary format to an output file. Each factor
 * is written as three uint64_t values (start position, length, ref).
 *
 * @param in_path Path to input file containing DNA text with multiple sequences
 * @param out_path Path to output file where binary factors will be written
 * @param start_pos Starting position for factorization (default: 0)
 * @return Number of factors written to the output file
 *
 * @note Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref)
 * @note Footer is written at the END of the file
 * @note This function overwrites the output file if it exists
 * @warning Ensure sufficient disk space for the output file
 */
size_t write_factors_binary_file_multiple_dna_w_rc(const std::string& in_path, const std::string& out_path, size_t start_pos) {
    // Read input file
    std::ifstream infile(in_path, std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Cannot open input file: " + in_path);
    }
    std::string text((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
    
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Calculate total_length from input text length minus start position
    uint64_t total_length = text.length() - start_pos;
    
    // Stream factors directly to file without collecting in memory
    size_t n = 0;
    detail::nolzss_multiple_dna_w_rc(text, [&](const Factor& f){
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
        ++n;
    }, start_pos);
    
    // For multiple DNA sequences from text file, we don't know sequence names or sentinels
    FactorFileFooter footer;
    footer.num_factors = n;
    footer.num_sequences = 0;  // Unknown for raw text files with multiple sequences
    footer.num_sentinels = 0;  // Cannot identify sentinels without preparation function
    footer.footer_size = sizeof(FactorFileFooter);
    footer.total_length = total_length;
    
    // Write footer at the end
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
    
    return n;
}

// Reference sequence factorization functions

/**
 * @brief Factorizes a target DNA sequence using a reference DNA sequence with reverse complement awareness.
 *
 * This function allows factorization of a target DNA sequence where factors can reference
 * positions in a reference sequence (or its reverse complement). This is useful for:
 * - Comparing related genomes (e.g., different strains of the same organism)
 * - Identifying similarities and differences between sequences
 * - Compression where a reference genome is available
 *
 * Algorithm:
 * 1. Prepares both sequences with reverse complements: REF[s1]TARGET[s2]RC(TARGET)[s3]RC(REF)[s4]
 * 2. Builds a single suffix tree containing all sequences
 * 3. Factorizes starting from the TARGET sequence position
 * 4. Factors can reference any position in REF, TARGET, or their reverse complements
 *
 * The MSB of the reference field indicates reverse complement matches.
 *
 * @param reference_seq Reference DNA sequence (should contain only A, C, T, G)
 * @param target_seq Target DNA sequence to factorize (should contain only A, C, T, G)
 * @return Vector of factors representing the target sequence
 *
 * @note Factors cover only the target sequence, but can reference the reference sequence
 * @note The reference field in factors points to positions in the combined prepared string
 * @note Both sequences are converted to uppercase
 *
 * @throws std::invalid_argument If sequences contain invalid nucleotides
 * @throws std::runtime_error If sequence preparation fails
 *
 * @see factorize_dna_w_reference_seq_file() for file output version
 * @see factorize_w_reference() for general (non-DNA) reference factorization
 */
std::vector<Factor> factorize_dna_w_reference_seq(const std::string& reference_seq, const std::string& target_seq) {
    // Prepare reference and target sequences together
    std::vector<std::string> sequences = {reference_seq, target_seq};
    PreparedSequenceResult prep_result = prepare_multiple_dna_sequences_w_rc(sequences);
    
    // Calculate the starting position of the target sequence in the prepared string
    // The format is: REF[sentinel]TARGET[sentinel]RC(TARGET)[sentinel]RC(REF)[sentinel]
    // We want to start factorization from where TARGET begins
    size_t target_start_pos = reference_seq.length() + 1; // +1 for the sentinel between ref and target
    
    // Perform factorization starting from target sequence
    std::vector<Factor> factors;
    detail::nolzss_multiple_dna_w_rc(prep_result.prepared_string, [&](const Factor& f) {
        factors.push_back(f);
    }, target_start_pos);
    
    return factors;
}

/**
 * @brief Factorizes a target DNA sequence using a reference sequence and writes factors to a binary file.
 *
 * This is the file output version of factorize_dna_w_reference_seq(). It performs the same
 * reference-based DNA factorization but writes the results directly to a binary file in the
 * noLZSS factor format with metadata footer.
 *
 * The output file format:
 * - Factors: Binary array of Factor structs (24 bytes each: start, length, reference)
 * - Footer: Metadata including factor count, sequence count (2), sentinel count (1)
 *
 * This is useful for:
 * - Processing large sequences without storing all factors in memory
 * - Saving factorization results for later analysis
 * - Feeding results into other tools that read noLZSS binary format
 *
 * @param reference_seq Reference DNA sequence (should contain only A, C, T, G)
 * @param target_seq Target DNA sequence to factorize (should contain only A, C, T, G)
 * @param out_path Path to output binary file (will be overwritten if exists)
 * @return Number of factors written to the file
 *
 * @note Output file includes footer with num_sequences=2, num_sentinels=1
 * @note File uses buffered I/O (1MB buffer) for performance
 * @note The reference field in factors points to positions in the combined prepared string
 *
 * @throws std::runtime_error If output file cannot be created
 * @throws std::invalid_argument If sequences contain invalid nucleotides
 *
 * @see factorize_dna_w_reference_seq() for in-memory version
 * @see factorize_w_reference_file() for general (non-DNA) reference factorization to file
 */
size_t factorize_dna_w_reference_seq_file(const std::string& reference_seq, const std::string& target_seq, const std::string& out_path) {
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Get factors using the in-memory function
    std::vector<Factor> factors = factorize_dna_w_reference_seq(reference_seq, target_seq);
    
    // Calculate total_length from target sequence length (what we're factorizing)
    uint64_t total_length = target_seq.length();
    
    // Write factors first
    for (const Factor& f : factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    // Create footer for reference+target factorization
    FactorFileFooter footer;
    footer.num_factors = factors.size();
    footer.num_sequences = 2;  // Reference + target sequences
    footer.num_sentinels = 1;  // One sentinel between ref and target (in factorized region)
    footer.footer_size = sizeof(FactorFileFooter);
    footer.total_length = total_length;
    
    // Write footer at the end
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
    
    return factors.size();
}

// General reference sequence factorization functions (no reverse complement)

/**
 * @brief Factorizes a target sequence using a reference sequence (general, non-DNA version).
 *
 * This is the general-purpose version of reference-based factorization that works with
 * any text, not just DNA. Unlike the DNA version, it does NOT consider reverse complements.
 * This is useful for:
 * - Non-genomic text compression with a reference document
 * - Finding similarities between general text documents
 * - Analyzing protein sequences or other non-DNA biological data
 *
 * Algorithm:
 * 1. Concatenates reference and target with a sentinel (ASCII 1) between them
 * 2. Builds a compressed suffix tree on the combined string
 * 3. Factorizes starting from the target sequence position
 * 4. Factors can reference any position in the reference or target
 *
 * @param reference_seq Reference sequence (any text)
 * @param target_seq Target sequence to factorize (any text)
 * @return Vector of factors representing the target sequence
 *
 * @note Factors cover only the target sequence, but can reference the reference sequence
 * @note The reference field in factors points to positions in the combined string
 * @note Sentinel character (ASCII 1) separates reference from target
 * @note No reverse complement awareness - use factorize_dna_w_reference_seq() for DNA
 *
 * @see factorize_w_reference_file() for file output version
 * @see factorize_dna_w_reference_seq() for DNA-specific version with reverse complement
 */
std::vector<Factor> factorize_w_reference(const std::string& reference_seq, const std::string& target_seq) {
    // Concatenate reference and target with a sentinel (use ASCII character 1)
    std::string combined = reference_seq + '\x01' + target_seq;
    
    // Calculate the starting position of the target sequence in the combined string
    size_t target_start_pos = reference_seq.length() + 1; // +1 for the sentinel
    
    // Perform general factorization starting from target sequence position
    std::vector<Factor> factors;
    cst_t cst; construct_im(cst, combined, 1);
    detail::nolzss(cst, [&](const Factor& f) {
        factors.push_back(f);
    }, target_start_pos);
    
    return factors;
}

/**
 * @brief Factorizes a target sequence using a reference sequence and writes factors to a binary file (general version).
 *
 * This is the file output version of factorize_w_reference(). It performs general
 * reference-based factorization (no reverse complement) and writes results directly
 * to a binary file in the noLZSS factor format with metadata footer.
 *
 * The output file format:
 * - Factors: Binary array of Factor structs (24 bytes each: start, length, reference)
 * - Footer: Metadata including factor count, sequence count (2), sentinel count (1)
 *
 * Use cases:
 * - Processing large non-DNA sequences without storing all factors in memory
 * - Saving factorization results for later analysis
 * - Comparing general text documents with a reference
 *
 * @param reference_seq Reference sequence (any text)
 * @param target_seq Target sequence to factorize (any text)
 * @param out_path Path to output binary file (will be overwritten if exists)
 * @return Number of factors written to the file
 *
 * @note Output file includes footer with num_sequences=2, num_sentinels=1
 * @note File uses buffered I/O (1MB buffer) for performance
 * @note The reference field in factors points to positions in the combined string
 * @note No reverse complement awareness - this is for general text, not DNA
 *
 * @throws std::runtime_error If output file cannot be created
 *
 * @see factorize_w_reference() for in-memory version
 * @see factorize_dna_w_reference_seq_file() for DNA-specific version with reverse complement
 */
size_t factorize_w_reference_file(const std::string& reference_seq, const std::string& target_seq, const std::string& out_path) {
    // Set up binary output file with buffering
    std::ofstream os(out_path, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Cannot create output file: " + out_path);
    }
    
    std::vector<char> buf(1<<20); // 1 MB buffer for performance
    os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    
    // Get factors using the in-memory function
    std::vector<Factor> factors = factorize_w_reference(reference_seq, target_seq);
    
    // Calculate total_length from target sequence length (what we're factorizing)
    uint64_t total_length = target_seq.length();
    
    // Write factors first
    for (const Factor& f : factors) {
        os.write(reinterpret_cast<const char*>(&f), sizeof(Factor));
    }
    
    // Create footer for reference+target factorization
    FactorFileFooter footer;
    footer.num_factors = factors.size();
    footer.num_sequences = 2;  // Reference + target sequences
    footer.num_sentinels = 1;  // One sentinel between ref and target
    footer.footer_size = sizeof(FactorFileFooter);
    footer.total_length = total_length;
    
    // Write footer at the end
    os.write(reinterpret_cast<const char*>(&footer), sizeof(footer));
    
    return factors.size();
}

// Parallel factorization implementations

size_t parallel_factorize_to_file(::std::string_view text, const ::std::string& output_path, size_t num_threads, size_t start_pos) {
    ParallelFactorizer factorizer;
    return factorizer.parallel_factorize(text, output_path, num_threads, start_pos);
}

size_t parallel_factorize_file_to_file(const ::std::string& input_path, const ::std::string& output_path, size_t num_threads, size_t start_pos) {
    ParallelFactorizer factorizer;
    return factorizer.parallel_factorize_file(input_path, output_path, num_threads, start_pos);
}

size_t parallel_factorize_dna_w_rc_to_file(::std::string_view text, const ::std::string& output_path, size_t num_threads) {
    ParallelFactorizer factorizer;
    return factorizer.parallel_factorize_dna_w_rc(text, output_path, num_threads);
}

size_t parallel_factorize_file_dna_w_rc_to_file(const ::std::string& input_path, const ::std::string& output_path, size_t num_threads) {
    ParallelFactorizer factorizer;
    return factorizer.parallel_factorize_file_dna_w_rc(input_path, output_path, num_threads);
}

} // namespace noLZSS
