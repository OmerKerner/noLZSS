/**
 * @file bindings.cpp
 * @brief Python bindings for the noLZSS factorization library.
 *
 * This file contains the Python bindings for the non-overlapping Lempel-Ziv-Storer-Szymanski
 * factorization algorithm. The bindings provide both in-memory and file-based factorization
 * capabilities with proper GIL management for performance.
 *
 * The module exposes the following functions:
 * - factorize(): Factorize in-memory text
 * - factorize_file(): Factorize text from file
 * - count_factors(): Count factors in text
 * - count_factors_file(): Count factors in file
 * - write_factors_binary_file(): Write factors to binary file
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/pytypes.h>
#include <string>
#include <string_view>
#include <stdexcept>
#include "factorizer.hpp"
#include "fasta_processor.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_noLZSS, m) {
    m.doc() = "Non-overlapping Lempel-Ziv-Storer-Szymanski factorization\n\n"
              "This module provides efficient text factorization using compressed suffix trees.";

    // Factor class documentation
    py::class_<noLZSS::Factor>(m, "Factor", "Represents a single factorization factor with start position, length, and reference position")
        .def_readonly("start", &noLZSS::Factor::start, "Starting position of the factor in the original text")
        .def_readonly("length", &noLZSS::Factor::length, "Length of the factor substring")
        .def_property_readonly("ref", [](const noLZSS::Factor& f) { return f.ref & ~noLZSS::RC_MASK; }, "Reference position with RC_MASK stripped if it's a reverse complement match")
        .def_property_readonly("is_rc", [](const noLZSS::Factor& f) { return noLZSS::is_rc(f.ref); }, "Whether this factor is a reverse complement match");

    // FastaFactorizationResult class documentation
    py::class_<noLZSS::FastaFactorizationResult>(m, "FastaFactorizationResult", "Result of FASTA factorization containing factors and sentinel information")
        .def_readonly("factors", &noLZSS::FastaFactorizationResult::factors, "List of factorization factors")
        .def_readonly("sentinel_factor_indices", &noLZSS::FastaFactorizationResult::sentinel_factor_indices, "Indices of factors that are sentinels (sequence separators)");

    // factorize function documentation
    m.def("factorize", [](py::buffer b) {
        // Accept any bytes-like 1-byte-per-item contiguous buffer (e.g. bytes, bytearray, memoryview)
        py::buffer_info info = b.request();
        if (info.itemsize != 1) {
            throw std::invalid_argument("factorize: buffer must be a bytes-like object with itemsize==1");
        }
        if (info.ndim != 1) {
            throw std::invalid_argument("factorize: buffer must be a 1-dimensional bytes-like object");
        }

        const char* data = static_cast<const char*>(info.ptr);
        std::string_view sv(data, static_cast<size_t>(info.size));

        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize(sv);
        py::gil_scoped_acquire acquire;

        py::list out;
        for (auto &f : factors) out.append(py::make_tuple(f.start, f.length, f.ref));
        return out;
    }, py::arg("data"), R"doc(Factorize a text string into noLZSS factors.

This is the main factorization function for in-memory text processing.
It accepts any Python bytes-like object and returns a list of (start, length) tuples.

Args:
    data: Python bytes-like object containing text

Returns:
    List of (start, length, ref) tuples representing the factorization

Raises:
    ValueError: if data is not a valid bytes-like object

Note:
    GIL is released during computation for better performance with large data.
)doc");

    // factorize_file function documentation
    m.def("factorize_file", [](const std::string& path, size_t reserve_hint) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize_file(path, reserve_hint);
        py::gil_scoped_acquire acquire;

        py::list out;
        for (auto &f : factors) out.append(py::make_tuple(f.start, f.length, f.ref));
        return out;
    }, py::arg("path"), py::arg("reserve_hint") = 0, R"doc(Factorize text from file into noLZSS factors.

Reads text from a file and performs factorization. This is more memory-efficient
for large files as it avoids loading the entire file into memory.

Args:
    path: Path to input file containing text
    reserve_hint: Optional hint for reserving space in output vector (0 = no hint)

Returns:
    List of (start, length, ref) tuples representing the factorization

Note:
    Use reserve_hint for better performance when you know approximate factor count.
)doc");

    // count_factors function documentation
    m.def("count_factors", [](py::buffer b) {
        // Accept any bytes-like 1-byte-per-item contiguous buffer
        py::buffer_info info = b.request();
        if (info.itemsize != 1) {
            throw std::invalid_argument("count_factors: buffer must be a bytes-like object with itemsize==1");
        }
        if (info.ndim != 1) {
            throw std::invalid_argument("count_factors: buffer must be a 1-dimensional bytes-like object");
        }

        const char* data = static_cast<const char*>(info.ptr);
        std::string_view sv(data, static_cast<size_t>(info.size));

        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        size_t count = noLZSS::count_factors(sv);
        py::gil_scoped_acquire acquire;

        return count;
    }, py::arg("data"), R"doc(Count number of LZSS factors in text.

This is a memory-efficient alternative to factorize() when you only need
the count of factors rather than the factors themselves.

Args:
    data: Python bytes-like object containing text

Returns:
    Number of factors in the factorization

Note:
    GIL is released during computation for better performance with large data.
)doc");

    // count_factors_file function documentation
    m.def("count_factors_file", [](const std::string& path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        size_t count = noLZSS::count_factors_file(path);
        py::gil_scoped_acquire acquire;

        return count;
    }, py::arg("path"), R"doc(Count number of noLZSS factors in a file.

Reads text from a file and counts factors without storing them.
This is the most memory-efficient way to get factor counts for large files.

Args:
    path: Path to input file containing text

Returns:
    Number of factors in the factorization

Note:
    GIL is released during computation for better performance.
)doc");

    // write_factors_binary_file function documentation
    m.def("write_factors_binary_file", [](const std::string& in_path, const std::string& out_path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        size_t count = noLZSS::write_factors_binary_file(in_path, out_path);
        py::gil_scoped_acquire acquire;

        return count;
    }, py::arg("in_path"), py::arg("out_path"), R"doc(Write noLZSS factors from file to binary output file.

Reads text from an input file, performs factorization, and writes the factors
in binary format with metadata header to an output file.

Args:
    in_path: Path to input file containing text
    out_path: Path to output file where binary factors will be written

Returns:
    Number of factors written to the output file

Note:
    Binary format: header with metadata + factors as 24 bytes each (3 × uint64_t: start, length, ref).
    This function overwrites the output file if it exists.
)doc");

    // DNA-aware factorization functions with reverse complement support

    // factorize_dna_w_rc function documentation
    m.def("factorize_dna_w_rc", [](py::buffer b) {
        // Accept any bytes-like 1-byte-per-item contiguous buffer (e.g. bytes, bytearray, memoryview)
        py::buffer_info info = b.request();
        if (info.itemsize != 1) {
            throw std::invalid_argument("factorize_dna_w_rc: buffer must be a bytes-like object with itemsize==1");
        }
        if (info.ndim != 1) {
            throw std::invalid_argument("factorize_dna_w_rc: buffer must be a 1-dimensional bytes-like object");
        }

        const char* data = static_cast<const char*>(info.ptr);
        std::string_view sv(data, static_cast<size_t>(info.size));

        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize_dna_w_rc(sv);
        py::gil_scoped_acquire acquire;

        py::list out;
        for (auto &f : factors) out.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
        return out;
    }, py::arg("data"), R"doc(Factorize DNA text with reverse complement awareness into noLZSS factors.

Performs non-overlapping Lempel-Ziv-Storer-Szymanski factorization on DNA sequences,
considering both forward and reverse complement matches. This is particularly useful
for genomic data where reverse complement patterns are biologically significant.

Args:
    data: Python bytes-like object containing DNA text

Returns:
    List of (start, length, ref, is_rc) tuples representing the factorization

Raises:
    ValueError: if data is not a valid bytes-like object

Note:
    ref field has RC_MASK cleared. is_rc boolean indicates if this was a reverse complement match.
    GIL is released during computation for better performance with large data.
)doc");

    // factorize_file_dna_w_rc function documentation
    m.def("factorize_file_dna_w_rc", [](const std::string& path, size_t reserve_hint) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize_file_dna_w_rc(path, reserve_hint);
        py::gil_scoped_acquire acquire;

        py::list out;
        for (auto &f : factors) out.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
        return out;
    }, py::arg("path"), py::arg("reserve_hint") = 0, R"doc(Factorize DNA text from file with reverse complement awareness into noLZSS factors.

Reads DNA text from a file and performs factorization considering both forward
and reverse complement matches. This is more memory-efficient for large genomic files.

Args:
    path: Path to input file containing DNA text
    reserve_hint: Optional hint for reserving space in output vector (0 = no hint)

Returns:
    List of (start, length, ref, is_rc) tuples representing the factorization

Note:
    Use reserve_hint for better performance when you know approximate factor count.
    ref field has RC_MASK cleared. is_rc boolean indicates if this was a reverse complement match.
)doc");

    // count_factors_dna_w_rc function documentation
    m.def("count_factors_dna_w_rc", [](py::buffer b) {
        // Accept any bytes-like 1-byte-per-item contiguous buffer
        py::buffer_info info = b.request();
        if (info.itemsize != 1) {
            throw std::invalid_argument("count_factors_dna_w_rc: buffer must be a bytes-like object with itemsize==1");
        }
        if (info.ndim != 1) {
            throw std::invalid_argument("count_factors_dna_w_rc: buffer must be a 1-dimensional bytes-like object");
        }

        const char* data = static_cast<const char*>(info.ptr);
        std::string_view sv(data, static_cast<size_t>(info.size));

        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        size_t count = noLZSS::count_factors_dna_w_rc(sv);
        py::gil_scoped_acquire acquire;

        return count;
    }, py::arg("data"), R"doc(Count number of noLZSS factors in DNA text with reverse complement awareness.

This is a memory-efficient alternative to factorize_dna_w_rc() when you only need
the count of factors rather than the factors themselves.

Args:
    data: Python bytes-like object containing DNA text

Returns:
    Number of factors in the factorization

Note:
    GIL is released during computation for better performance with large data.
)doc");

    // count_factors_file_dna_w_rc function documentation
    m.def("count_factors_file_dna_w_rc", [](const std::string& path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        size_t count = noLZSS::count_factors_file_dna_w_rc(path);
        py::gil_scoped_acquire acquire;

        return count;
    }, py::arg("path"), R"doc(Count number of noLZSS factors in a DNA file with reverse complement awareness.

Reads DNA text from a file and counts factors without storing them.
This is the most memory-efficient way to get factor counts for large genomic files.

Args:
    path: Path to input file containing DNA text

Returns:
    Number of factors in the factorization

Note:
    GIL is released during computation for better performance.
)doc");

    // write_factors_binary_file_dna_w_rc function documentation
    m.def("write_factors_binary_file_dna_w_rc", [](const std::string& in_path, const std::string& out_path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        size_t count = noLZSS::write_factors_binary_file_dna_w_rc(in_path, out_path);
        py::gil_scoped_acquire acquire;

        return count;
    }, py::arg("in_path"), py::arg("out_path"), R"doc(Write noLZSS factors from DNA file with reverse complement awareness to binary output file.

Reads DNA text from an input file, performs factorization with reverse complement support,
and writes the factors in binary format with metadata header to an output file.

Args:
    in_path: Path to input file containing DNA text
    out_path: Path to output file where binary factors will be written

Returns:
    Number of factors written to the output file

Note:
    Binary format: header with metadata + factors as 24 bytes each (3 × uint64_t: start, length, ref).
    Reverse complement factors have RC_MASK set in the ref field.
    This function overwrites the output file if it exists.
)doc");


    // factorize_multiple_dna_w_rc function documentation
    m.def("factorize_multiple_dna_w_rc", [](py::buffer b) {
        // Accept any bytes-like 1-byte-per-item contiguous buffer
        py::buffer_info info = b.request();
        if (info.itemsize != 1) {
            throw std::invalid_argument("factorize_multiple_dna_w_rc: buffer must be a bytes-like object with itemsize==1");
        }
        if (info.ndim != 1) {
            throw std::invalid_argument("factorize_multiple_dna_w_rc: buffer must be a 1-dimensional bytes-like object");
        }

        const char* data = static_cast<const char*>(info.ptr);
        std::string_view sv(data, static_cast<size_t>(info.size));

        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize_multiple_dna_w_rc(sv);
        py::gil_scoped_acquire acquire;

        py::list out;
        for (auto &f : factors) out.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
        return out;
    }, py::arg("data"), R"doc(Factorize DNA text with multiple sequences and reverse complement awareness.

Performs non-overlapping Lempel-Ziv-Storer-Szymanski factorization on DNA text
containing multiple sequences separated by sentinels, considering both forward
and reverse complement matches.

Args:
    data: Python bytes-like object containing DNA text with multiple sequences and sentinels

Returns:
    List of (start, length, ref, is_rc) tuples representing the factorization

Note:
    ref field has RC_MASK cleared. is_rc boolean indicates if this was a reverse complement match.
    GIL is released during computation for better performance with large data.
)doc");

    // factorize_file_multiple_dna_w_rc function documentation
    m.def("factorize_file_multiple_dna_w_rc", [](const std::string& path, size_t reserve_hint) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize_file_multiple_dna_w_rc(path, reserve_hint);
        py::gil_scoped_acquire acquire;

        py::list out;
        for (auto &f : factors) out.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
        return out;
    }, py::arg("path"), py::arg("reserve_hint") = 0, R"doc(Factorize DNA text from file with multiple sequences and reverse complement awareness.

Reads DNA text from a file and performs factorization with multiple sequences
and reverse complement matches.

Args:
    path: Path to input file containing DNA text with multiple sequences
    reserve_hint: Optional hint for reserving space in output vector (0 = no hint)

Returns:
    List of (start, length, ref, is_rc) tuples representing the factorization

Note:
    Use reserve_hint for better performance when you know approximate factor count.
    ref field has RC_MASK cleared. is_rc boolean indicates if this was a reverse complement match.
)doc");

    // count_factors_multiple_dna_w_rc function documentation
    m.def("count_factors_multiple_dna_w_rc", [](py::buffer b) {
        // Accept any bytes-like 1-byte-per-item contiguous buffer
        py::buffer_info info = b.request();
        if (info.itemsize != 1) {
            throw std::invalid_argument("count_factors_multiple_dna_w_rc: buffer must be a bytes-like object with itemsize==1");
        }
        if (info.ndim != 1) {
            throw std::invalid_argument("count_factors_multiple_dna_w_rc: buffer must be a 1-dimensional bytes-like object");
        }

        const char* data = static_cast<const char*>(info.ptr);
        std::string_view sv(data, static_cast<size_t>(info.size));

        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        size_t count = noLZSS::count_factors_multiple_dna_w_rc(sv);
        py::gil_scoped_acquire acquire;

        return count;
    }, py::arg("data"), R"doc(Count number of LZSS factors in DNA text with multiple sequences and reverse complement awareness.

This is a memory-efficient alternative to factorize_multiple_dna_w_rc() when you only need
the count of factors rather than the factors themselves.

Args:
    data: Python bytes-like object containing DNA text with multiple sequences

Returns:
    Number of factors in the factorization

Note:
    GIL is released during computation for better performance with large data.
)doc");

    // count_factors_file_multiple_dna_w_rc function documentation
    m.def("count_factors_file_multiple_dna_w_rc", [](const std::string& path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        size_t count = noLZSS::count_factors_file_multiple_dna_w_rc(path);
        py::gil_scoped_acquire acquire;

        return count;
    }, py::arg("path"), R"doc(Count number of noLZSS factors in a DNA file with multiple sequences and reverse complement awareness.

Reads DNA text from a file and counts factors with multiple sequences and
reverse complement matches.

Args:
    path: Path to input file containing DNA text with multiple sequences

Returns:
    Number of factors in the factorization

Note:
    GIL is released during computation for better performance with large data.
)doc");

    // write_factors_binary_file_multiple_dna_w_rc function documentation
    m.def("write_factors_binary_file_multiple_dna_w_rc", [](const std::string& in_path, const std::string& out_path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        size_t count = noLZSS::write_factors_binary_file_multiple_dna_w_rc(in_path, out_path);
        py::gil_scoped_acquire acquire;

        return count;
    }, py::arg("in_path"), py::arg("out_path"), R"doc(Write noLZSS factors from DNA file with multiple sequences and reverse complement awareness to binary output file.

Reads DNA text from an input file, performs factorization with multiple sequences and reverse complement support,
and writes the factors in binary format with metadata header to an output file.

Args:
    in_path: Path to input file containing DNA text with multiple sequences
    out_path: Path to output file where binary factors will be written

Returns:
    Number of factors written to the output file

Note:
    Binary format: header with metadata + factors as 24 bytes each (3 × uint64_t: start, length, ref).
    Reverse complement factors have RC_MASK set in the ref field.
    This function overwrites the output file if it exists.
)doc");

    // FASTA processing function
    m.def("process_nucleotide_fasta", [](const std::string& fasta_path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto result = noLZSS::process_nucleotide_fasta(fasta_path);
        py::gil_scoped_acquire acquire;

        // Return as Python dictionary
        py::dict py_result;
        py_result["sequence"] = result.sequence;
        py_result["num_sequences"] = result.num_sequences;
        py_result["sequence_ids"] = result.sequence_ids;
        py_result["sequence_lengths"] = result.sequence_lengths;
        py_result["sequence_positions"] = result.sequence_positions;
        return py_result;
    }, py::arg("fasta_path"), R"doc(Process a nucleotide FASTA file into concatenated string with sentinels.

Reads a FASTA file containing nucleotide sequences and creates a single concatenated
string with sentinel characters separating sequences. Only A, C, T, G nucleotides
are allowed (case insensitive, converted to uppercase).

Args:
    fasta_path: Path to the FASTA file

Returns:
    Dictionary containing:
    - 'sequence': Concatenated sequences with sentinels
    - 'num_sequences': Number of sequences processed
    - 'sequence_ids': List of sequence IDs
    - 'sequence_lengths': List of sequence lengths (excluding sentinels)  
    - 'sequence_positions': List of start positions in concatenated string

Raises:
    RuntimeError: If file cannot be read, contains invalid nucleotides,
                 or has more than 251 sequences (sentinel limit)

    Note:
    Sentinels are characters 1-251 (avoiding 0, A=65, C=67, G=71, T=84).
    Empty sequences are skipped. Only whitespace is ignored in sequences.
)doc");

    // Amino acid FASTA processing function
    m.def("process_amino_acid_fasta", [](const std::string& fasta_path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto result = noLZSS::process_amino_acid_fasta(fasta_path);
        py::gil_scoped_acquire acquire;

        // Return as Python dictionary
        py::dict py_result;
        py_result["sequence"] = result.sequence;
        py_result["num_sequences"] = result.num_sequences;
        py_result["sequence_ids"] = result.sequence_ids;
        py_result["sequence_lengths"] = result.sequence_lengths;
        py_result["sequence_positions"] = result.sequence_positions;
        return py_result;
    }, py::arg("fasta_path"), R"doc(Process an amino acid FASTA file into concatenated string with sentinels.

Reads a FASTA file containing amino acid sequences and creates a single concatenated
string with sentinel characters separating sequences. Only canonical amino acids
are allowed (case insensitive, converted to uppercase).

Args:
    fasta_path: Path to the FASTA file

Returns:
    Dictionary containing:
    - 'sequence': Concatenated sequences with sentinels
    - 'num_sequences': Number of sequences processed
    - 'sequence_ids': List of sequence IDs
    - 'sequence_lengths': List of sequence lengths (excluding sentinels)  
    - 'sequence_positions': List of start positions in concatenated string

Raises:
    RuntimeError: If file cannot be read, contains invalid amino acids,
                 or has more than 235 sequences (sentinel limit)

    Note:
    Canonical amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
    Sentinels are characters 1-251 (avoiding 0 and common amino acids).
    Empty sequences are skipped. Only whitespace is ignored in sequences.
)doc");
// FASTA factorization function
m.def("factorize_fasta_multiple_dna_w_rc", [](const std::string& fasta_path) {
    // Release GIL while doing heavy C++ work
    py::gil_scoped_release release;
    auto result = noLZSS::factorize_fasta_multiple_dna_w_rc(fasta_path);
    py::gil_scoped_acquire acquire;

    // Convert factors to Python tuples
    py::list factors_list;
    for (auto &f : result.factors) {
        factors_list.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
    }
    
    // Convert sentinel indices to Python list
    py::list sentinel_indices_list;
    for (auto idx : result.sentinel_factor_indices) {
        sentinel_indices_list.append(idx);
    }
    
    // Convert sequence IDs to Python list
    py::list sequence_ids_list;
    for (const auto& seq_id : result.sequence_ids) {
        sequence_ids_list.append(seq_id);
    }
    
    return py::make_tuple(factors_list, sentinel_indices_list, sequence_ids_list);
}, py::arg("fasta_path"), R"doc(Factorize multiple DNA sequences from a FASTA file with reverse complement awareness.

Reads a FASTA file containing DNA sequences, parses them into individual sequences,
prepares them for factorization using prepare_multiple_dna_sequences_w_rc(), and then
performs noLZSS factorization with reverse complement awareness.

Args:
fasta_path: Path to the FASTA file containing DNA sequences

Returns:
Tuple of (factors, sentinel_factor_indices, sequence_ids) where:
- factors: List of (start, length, ref, is_rc) tuples representing the factorization
- sentinel_factor_indices: List of factor indices that represent sequence separators
- sequence_ids: List of sequence identifiers from FASTA headers

Raises:
RuntimeError: If FASTA file cannot be opened or contains no valid sequences
ValueError: If too many sequences (>125) in the FASTA file or invalid nucleotides found

Note:
Only A, C, T, G nucleotides are allowed (case insensitive)
Sequences are converted to uppercase before factorization
Reverse complement matches are supported during factorization
Nucleotide validation is performed by prepare_multiple_dna_sequences_w_rc()
ref field has RC_MASK cleared. is_rc boolean indicates if this was a reverse complement match.
)doc");

// FASTA factorization function (no reverse complement)
m.def("factorize_fasta_multiple_dna_no_rc", [](const std::string& fasta_path) {
    // Release GIL while doing heavy C++ work
    py::gil_scoped_release release;
    auto result = noLZSS::factorize_fasta_multiple_dna_no_rc(fasta_path);
    py::gil_scoped_acquire acquire;

    // Convert factors to Python tuples
    py::list factors_list;
    for (auto &f : result.factors) {
        factors_list.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
    }
    
    // Convert sentinel indices to Python list
    py::list sentinel_indices_list;
    for (auto idx : result.sentinel_factor_indices) {
        sentinel_indices_list.append(idx);
    }
    
    // Convert sequence IDs to Python list
    py::list sequence_ids_list;
    for (const auto& seq_id : result.sequence_ids) {
        sequence_ids_list.append(seq_id);
    }
    
    return py::make_tuple(factors_list, sentinel_indices_list, sequence_ids_list);
}, py::arg("fasta_path"), R"doc(Factorize multiple DNA sequences from a FASTA file without reverse complement awareness.

Reads a FASTA file containing DNA sequences, parses them into individual sequences,
prepares them for factorization using prepare_multiple_dna_sequences_no_rc(), and then
performs noLZSS factorization without reverse complement awareness.

Args:
fasta_path: Path to the FASTA file containing DNA sequences

Returns:
Tuple of (factors, sentinel_factor_indices, sequence_ids) where:
- factors: List of (start, length, ref, is_rc) tuples representing the factorization
- sentinel_factor_indices: List of factor indices that represent sequence separators
- sequence_ids: List of sequence identifiers from FASTA headers

Raises:
RuntimeError: If FASTA file cannot be opened or contains no valid sequences
ValueError: If too many sequences (>250) in the FASTA file or invalid nucleotides found

Note:
Only A, C, T, G nucleotides are allowed (case insensitive)
Sequences are converted to uppercase before factorization
Reverse complement matches are NOT supported during factorization
Nucleotide validation is performed by prepare_multiple_dna_sequences_no_rc()
ref field has RC_MASK cleared. is_rc boolean indicates if this was a reverse complement match.
)doc");

// FASTA binary factorization function (with reverse complement)
m.def("write_factors_binary_file_fasta_multiple_dna_w_rc", [](const std::string& fasta_path, const std::string& out_path) {
    // Release GIL while doing heavy C++ work
    py::gil_scoped_release release;
    size_t count = noLZSS::write_factors_binary_file_fasta_multiple_dna_w_rc(fasta_path, out_path);
    py::gil_scoped_acquire acquire;
    return count;
}, py::arg("fasta_path"), py::arg("out_path"), R"doc(Write noLZSS factors from multiple DNA sequences in a FASTA file with reverse complement awareness to a binary output file.

This function reads DNA sequences from a FASTA file, parses them into individual sequences,
prepares them for factorization, performs factorization with reverse complement awareness, 
and writes the resulting factors in binary format with metadata including sequence IDs and 
sentinel factor indices to an output file.

Args:
    fasta_path: Path to input FASTA file containing DNA sequences
    out_path: Path to output file where binary factors will be written

Returns:
    int: Number of factors written to the output file

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If too many sequences (>125) in the FASTA file or invalid nucleotides found

Note:
    Binary format: header with sequence metadata + factors as 24 bytes each (3 × uint64_t: start, length, ref)
    Header includes sequence IDs, sentinel factor indices, and other metadata
    Only A, C, T, G nucleotides are allowed (case insensitive)
    This function overwrites the output file if it exists
    Reverse complement matches are supported during factorization
)doc");

// FASTA binary factorization function (no reverse complement)
m.def("write_factors_binary_file_fasta_multiple_dna_no_rc", [](const std::string& fasta_path, const std::string& out_path) {
    // Release GIL while doing heavy C++ work
    py::gil_scoped_release release;
    size_t count = noLZSS::write_factors_binary_file_fasta_multiple_dna_no_rc(fasta_path, out_path);
    py::gil_scoped_acquire acquire;
    return count;
}, py::arg("fasta_path"), py::arg("out_path"), R"doc(Write noLZSS factors from multiple DNA sequences in a FASTA file without reverse complement awareness to a binary output file.

This function reads DNA sequences from a FASTA file, parses them into individual sequences,
prepares them for factorization, performs factorization without reverse complement awareness, 
and writes the resulting factors in binary format with metadata including sequence IDs and 
sentinel factor indices to an output file.

Args:
    fasta_path: Path to input FASTA file containing DNA sequences
    out_path: Path to output file where binary factors will be written

Returns:
    int: Number of factors written to the output file

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If too many sequences (>250) in the FASTA file or invalid nucleotides found

Note:
    Binary format: header with sequence metadata + factors as 24 bytes each (3 × uint64_t: start, length, ref)
    Header includes sequence IDs, sentinel factor indices, and other metadata
    Only A, C, T, G nucleotides are allowed (case insensitive)
    This function overwrites the output file if it exists
    Reverse complement matches are NOT supported during factorization
)doc");

    // DNA sequence preparation utility
    m.def("prepare_multiple_dna_sequences_w_rc", [](const std::vector<std::string>& sequences) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto result = noLZSS::prepare_multiple_dna_sequences_w_rc(sequences);
        py::gil_scoped_acquire acquire;

        // Return as Python tuple (concatenated_string, original_length, sentinel_positions)
        return py::make_tuple(result.prepared_string, result.original_length, result.sentinel_positions);
    }, py::arg("sequences"), R"doc(Prepare multiple DNA sequences for factorization with reverse complement awareness.

Takes multiple DNA sequences, concatenates them with unique sentinels, and appends
their reverse complements with unique sentinels. The output format is compatible
with nolzss_multiple_dna_w_rc(): S = T1!T2@T3$rt(T3)%rt(T2)^rt(T1)&

Args:
    sequences: List of DNA sequence strings (should contain only A, C, T, G)

Returns:
    Tuple containing:
    - concatenated_string: The formatted string with sequences and reverse complements
    - original_length: Length of the original sequences part (before reverse complements)
    - sentinel_positions: List of positions where sentinels are located

Raises:
    ValueError: If too many sequences (>125) or invalid nucleotides found
    RuntimeError: If sequences contain invalid characters

Note:
    Sentinels range from 1-251, avoiding 0, A(65), C(67), G(71), T(84).
    Input sequences can be lowercase or uppercase, output is always uppercase.
    The function validates that all sequences contain only valid DNA nucleotides.
)doc");

    // DNA sequence preparation utility (no reverse complement)
    m.def("prepare_multiple_dna_sequences_no_rc", [](const std::vector<std::string>& sequences) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto result = noLZSS::prepare_multiple_dna_sequences_no_rc(sequences);
        py::gil_scoped_acquire acquire;

        // Return as Python tuple (concatenated_string, total_length, sentinel_positions)
        return py::make_tuple(result.prepared_string, result.original_length, result.sentinel_positions);
    }, py::arg("sequences"), R"doc(Prepare multiple DNA sequences for factorization without reverse complement.

Takes multiple DNA sequences and concatenates them with unique sentinels.
Unlike prepare_multiple_dna_sequences_w_rc(), this function does not append
reverse complements. The output format is: S = T1!T2@T3$

Args:
    sequences: List of DNA sequence strings (should contain only A, C, T, G)

Returns:
    Tuple containing:
    - concatenated_string: The formatted string with sequences and sentinels
    - total_length: Total length of the concatenated string
    - sentinel_positions: List of positions where sentinels are located

Raises:
    ValueError: If too many sequences (>250) or invalid nucleotides found
    RuntimeError: If sequences contain invalid characters

Note:
    Sentinels range from 1-251, avoiding 0, A(65), C(67), G(71), T(84).
    Input sequences can be lowercase or uppercase, output is always uppercase.
    The function validates that all sequences contain only valid DNA nucleotides.
)doc");

    // Reference sequence factorization functions
    m.def("factorize_w_reference_seq", [](const std::string& reference_seq, const std::string& target_seq) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize_w_reference_seq(reference_seq, target_seq);
        py::gil_scoped_acquire acquire;

        // Convert factors to Python tuples
        py::list out;
        for (auto &f : factors) {
            out.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
        }
        return out;
    }, py::arg("reference_seq"), py::arg("target_seq"), R"doc(Factorize target DNA sequence using a reference sequence with reverse complement awareness.

Concatenates a reference sequence and target sequence, then performs noLZSS factorization
with reverse complement awareness starting from where the target sequence begins. This allows
the target sequence to reference patterns in the reference sequence without factorizing the
reference itself.

Args:
    reference_seq: Reference DNA sequence string (A, C, T, G - case insensitive)
    target_seq: Target DNA sequence string to be factorized (A, C, T, G - case insensitive)

Returns:
    List of (start, length, ref, is_rc) tuples representing the factorization of target sequence

Raises:
    ValueError: If sequences contain invalid nucleotides or are empty
    RuntimeError: If too many sequences or other processing errors

Note:
    Factor start positions are relative to the beginning of the target sequence.
    Both sequences are converted to uppercase before factorization.
    Reverse complement matches are supported during factorization.
    ref field has RC_MASK cleared. is_rc boolean indicates if this was a reverse complement match.
)doc");

    m.def("factorize_w_reference_seq_file", [](const std::string& reference_seq, const std::string& target_seq, const std::string& out_path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto num_factors = noLZSS::factorize_w_reference_seq_file(reference_seq, target_seq, out_path);
        py::gil_scoped_acquire acquire;
        return num_factors;
    }, py::arg("reference_seq"), py::arg("target_seq"), py::arg("out_path"), R"doc(Factorize target DNA sequence using a reference sequence and write factors to binary file.

Concatenates a reference sequence and target sequence, then performs noLZSS factorization
with reverse complement awareness starting from where the target sequence begins, and writes
the resulting factors to a binary file.

Args:
    reference_seq: Reference DNA sequence string (A, C, T, G - case insensitive)
    target_seq: Target DNA sequence string to be factorized (A, C, T, G - case insensitive)
    out_path: Path to output file where binary factors will be written

Returns:
    Number of factors written to the output file

Raises:
    ValueError: If sequences contain invalid nucleotides or are empty
    RuntimeError: If unable to create output file or other processing errors

Note:
    Factor start positions are relative to the beginning of the target sequence.
    Binary format follows the same structure as other DNA factorization binary outputs.
    This function overwrites the output file if it exists.
)doc");

    // Version information
#ifdef NOLZSS_VERSION
    m.attr("__version__") = NOLZSS_VERSION;
#else
    m.attr("__version__") = "0.0.0";
#endif
}
