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
#include "parallel_fasta_processor.hpp"

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

m.def("factorize_dna_rc_w_ref_fasta_files", [](const std::string& reference_fasta_path, const std::string& target_fasta_path) {
    py::gil_scoped_release release;
    auto result = noLZSS::factorize_dna_rc_w_ref_fasta_files(reference_fasta_path, target_fasta_path);
    py::gil_scoped_acquire acquire;

    py::list factors_list;
    for (auto &f : result.factors) {
        factors_list.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
    }

    py::list sentinel_indices_list;
    for (auto idx : result.sentinel_factor_indices) {
        sentinel_indices_list.append(idx);
    }

    py::list sequence_ids_list;
    for (const auto& seq_id : result.sequence_ids) {
        sequence_ids_list.append(seq_id);
    }

    return py::make_tuple(factors_list, sentinel_indices_list, sequence_ids_list);
}, py::arg("reference_fasta_path"), py::arg("target_fasta_path"), R"doc(Factorize DNA sequences from reference and target FASTA files with reverse complement awareness.

Reads reference and target FASTA files, concatenates their sequences with sentinels, and performs
reverse-complement-aware factorization starting from the target region.

Args:
reference_fasta_path: Path to reference FASTA file containing DNA sequences
target_fasta_path: Path to target FASTA file containing DNA sequences to factorize

Returns:
Tuple of (factors, sentinel_factor_indices, sequence_ids) where:
- factors: List of (start, length, ref, is_rc) tuples representing the factorization
- sentinel_factor_indices: List of factor indices that correspond to sentinel separators
- sequence_ids: List of sequence identifiers from both FASTA files (reference first)

Raises:
RuntimeError: If FASTA files cannot be read or contain no valid sequences
ValueError: If invalid nucleotides or too many sequences are encountered

Note:
- Factor start positions are absolute indices in the concatenated reference+target string
- Reverse complement matches are indicated by is_rc=True
- Sentinels mark sequence boundaries for both original and reverse-complement sections
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

    // DNA-specific reference sequence factorization functions
    m.def("factorize_dna_w_reference_seq", [](const std::string& reference_seq, const std::string& target_seq) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize_dna_w_reference_seq(reference_seq, target_seq);
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

    m.def("factorize_dna_w_reference_seq_file", [](const std::string& reference_seq, const std::string& target_seq, const std::string& out_path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto num_factors = noLZSS::factorize_dna_w_reference_seq_file(reference_seq, target_seq, out_path);
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

    // General reference sequence factorization functions (no reverse complement)
    m.def("factorize_w_reference", [](const std::string& reference_seq, const std::string& target_seq) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize_w_reference(reference_seq, target_seq);
        py::gil_scoped_acquire acquire;

        // Convert factors to Python tuples (no is_rc field for general factorization)
        py::list out;
        for (auto &f : factors) {
            out.append(py::make_tuple(f.start, f.length, f.ref));
        }
        return out;
    }, py::arg("reference_seq"), py::arg("target_seq"), R"doc(Factorize target sequence using a reference sequence without reverse complement.

Concatenates a reference sequence and target sequence, then performs noLZSS factorization
starting from where the target sequence begins. This allows the target sequence to reference
patterns in the reference sequence without factorizing the reference itself. Suitable for
general text or amino acid sequences.

Args:
    reference_seq: Reference sequence string (any text)
    target_seq: Target sequence string to be factorized (any text)

Returns:
    List of (start, length, ref) tuples representing the factorization of target sequence

Raises:
    ValueError: If sequences are empty
    RuntimeError: If processing errors occur

Note:
    Factor start positions are absolute positions in the combined reference+target string.
    No reverse complement matching is performed - suitable for text or amino acid sequences.

Warning:
    The sentinel character '\x01' (ASCII 1) must not appear in either input sequence,
    as it is used internally to separate the reference and target sequences.
)doc");

    m.def("factorize_w_reference_file", [](const std::string& reference_seq, const std::string& target_seq, const std::string& out_path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto num_factors = noLZSS::factorize_w_reference_file(reference_seq, target_seq, out_path);
        py::gil_scoped_acquire acquire;
        return num_factors;
    }, py::arg("reference_seq"), py::arg("target_seq"), py::arg("out_path"), R"doc(Factorize target sequence using a reference sequence and write factors to binary file.

Concatenates a reference sequence and target sequence, then performs noLZSS factorization
starting from where the target sequence begins, and writes the resulting factors to a binary file.
Suitable for general text or amino acid sequences.

Args:
    reference_seq: Reference sequence string (any text)
    target_seq: Target sequence string to be factorized (any text)
    out_path: Path to output file where binary factors will be written

Returns:
    Number of factors written to the output file

Raises:
    ValueError: If sequences are empty
    RuntimeError: If unable to create output file or processing errors

Note:
    Factor start positions are absolute positions in the combined reference+target string.
    No reverse complement matching is performed - suitable for text or amino acid sequences.
    This function overwrites the output file if it exists.

Warning:
    The sentinel character '\x01' (ASCII 1) must not appear in either input sequence,
    as it is used internally to separate the reference and target sequences.
)doc");

    m.def("write_factors_dna_w_reference_fasta_files_to_binary", [](const std::string& reference_fasta_path, const std::string& target_fasta_path, const std::string& out_path) {
        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto num_factors = noLZSS::write_factors_dna_w_reference_fasta_files_to_binary(reference_fasta_path, target_fasta_path, out_path);
        py::gil_scoped_acquire acquire;
        return num_factors;
    }, py::arg("reference_fasta_path"), py::arg("target_fasta_path"), py::arg("out_path"), R"doc(Factorize DNA sequences from FASTA files with reference and write factors to binary file.

Reads DNA sequences from reference and target FASTA files, performs noLZSS factorization
of the target using the reference, and writes the resulting factors to a binary file.
Specialized for nucleotide sequences (A, C, T, G) with reverse complement matching capability.

Args:
    reference_fasta_path: Path to reference FASTA file containing DNA sequences
    target_fasta_path: Path to target FASTA file containing DNA sequences to factorize
    out_path: Path to output file where binary factors will be written

Returns:
    Number of factors written to the output file

Raises:
    ValueError: If files contain empty sequences or invalid nucleotides
    RuntimeError: If unable to read FASTA files, create output file, or processing errors

Note:
    - Factor start positions are absolute positions in the combined reference+target string
    - Supports reverse complement matching for DNA sequences (indicated by MSB in reference field)
    - All sequences from both FASTA files are concatenated with sentinel separators
    - Only nucleotides A, C, T, G are allowed (case-insensitive)
    - This function overwrites the output file if it exists

Warning:
    Characters 1-251 are used as sentinel separators and must not appear in sequences.
)doc");

    // Parallel factorization bindings
    m.def("parallel_factorize_to_file", &noLZSS::parallel_factorize_to_file,
          py::arg("text"), py::arg("output_path"), py::arg("num_threads") = 0, py::arg("start_pos") = 0,
          R"doc(
Factorizes text in parallel and writes results to a binary file.

This function uses multiple threads to factorize the input text in parallel,
significantly improving performance on large inputs. The results are automatically
merged and written to a binary output file.

Args:
    text: Input text to factorize (str or bytes)
    output_path: Path to output binary file
    num_threads: Number of threads (0 for auto-detection based on CPU cores)
    start_pos: Starting position in the text for factorization (default: 0)

Returns:
    Number of factors produced

Raises:
    RuntimeError: If factorization fails or output file cannot be written

Note:
    - Thread count is automatically adjusted based on input size
    - Small inputs may use fewer threads than specified for efficiency
    - Temporary files are created and cleaned up automatically
)doc");

    m.def("parallel_factorize_file_to_file", &noLZSS::parallel_factorize_file_to_file,
          py::arg("input_path"), py::arg("output_path"), py::arg("num_threads") = 0, py::arg("start_pos") = 0,
          R"doc(
Factorizes text from file in parallel and writes results to a binary file.

Reads text from an input file, factorizes it using multiple threads, and writes
the results to a binary output file. This is more memory-efficient for large files.

Args:
    input_path: Path to input text file
    output_path: Path to output binary file
    num_threads: Number of threads (0 for auto-detection)
    start_pos: Starting position in the text for factorization (default: 0)

Returns:
    Number of factors produced

Raises:
    FileNotFoundError: If input file doesn't exist
    RuntimeError: If factorization fails

Note:
    - Suitable for processing very large text files
    - Memory usage scales with input size but is optimized for large files
)doc");

    m.def("parallel_factorize_dna_w_rc_to_file", &noLZSS::parallel_factorize_dna_w_rc_to_file,
          py::arg("text"), py::arg("output_path"), py::arg("num_threads") = 0,
          R"doc(
Factorizes DNA text in parallel with reverse complement and writes results to a binary file.

Performs parallel factorization on DNA sequences with reverse complement awareness.
This is particularly useful for genomic data where reverse complement patterns
are biologically significant.

Args:
    text: Input DNA text (should contain only A, C, T, G)
    output_path: Path to output binary file
    num_threads: Number of threads (0 for auto-detection)

Returns:
    Number of factors produced

Raises:
    ValueError: If input contains invalid nucleotides
    RuntimeError: If factorization fails

Note:
    - Reverse complement matches are encoded with RC_MASK in the ref field
    - Currently under development - may not be fully implemented
)doc");

    m.def("parallel_factorize_file_dna_w_rc_to_file", &noLZSS::parallel_factorize_file_dna_w_rc_to_file,
          py::arg("input_path"), py::arg("output_path"), py::arg("num_threads") = 0,
          R"doc(
Factorizes DNA text from file in parallel with reverse complement and writes results to a binary file.

Reads DNA text from a file, factorizes it in parallel with reverse complement awareness,
and writes the results to a binary output file.

Args:
    input_path: Path to input DNA text file
    output_path: Path to output binary file
    num_threads: Number of threads (0 for auto-detection)

Returns:
    Number of factors produced

Raises:
    FileNotFoundError: If input file doesn't exist
    ValueError: If file contains invalid nucleotides
    RuntimeError: If factorization fails

Note:
    - Optimized for large genomic datasets
    - Currently under development - may not be fully implemented
)doc");

    // =========================================================================
    // Parallel FASTA Processing Bindings
    // =========================================================================

    m.def("parallel_write_factors_binary_file_fasta_multiple_dna_w_rc",
          &noLZSS::parallel_write_factors_binary_file_fasta_multiple_dna_w_rc,
          py::arg("fasta_path"),
          py::arg("out_path"),
          py::arg("num_threads") = 0,
          R"doc(Parallel version of write_factors_binary_file_fasta_multiple_dna_w_rc.

Reads a FASTA file containing DNA sequences, prepares them for factorization with
reverse complement awareness, and performs parallel factorization writing results
to a binary output file with metadata.

Args:
    fasta_path: Path to input FASTA file containing DNA sequences
    out_path: Path to output file where binary factors will be written
    num_threads: Number of threads to use (0 = auto-detect based on input size, default)

Returns:
    Number of factors written to the output file

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If too many sequences (>125) or invalid nucleotides found

Note:
    - Binary format includes factors, sequence IDs, sentinel indices, and footer metadata
    - Only A, C, T, G nucleotides are allowed (case insensitive)
    - This function overwrites the output file if it exists
    - Reverse complement matches are supported during factorization
    - For single-threaded execution (num_threads=1), no temporary files are created
    - GIL is released during computation for better performance
)doc");

    m.def("parallel_write_factors_binary_file_fasta_multiple_dna_no_rc",
          &noLZSS::parallel_write_factors_binary_file_fasta_multiple_dna_no_rc,
          py::arg("fasta_path"),
          py::arg("out_path"),
          py::arg("num_threads") = 0,
          R"doc(Parallel version of write_factors_binary_file_fasta_multiple_dna_no_rc.

Reads a FASTA file containing DNA sequences, prepares them for factorization without
reverse complement awareness, and performs parallel factorization writing results
to a binary output file with metadata.

Args:
    fasta_path: Path to input FASTA file containing DNA sequences
    out_path: Path to output file where binary factors will be written
    num_threads: Number of threads to use (0 = auto-detect based on input size, default)

Returns:
    Number of factors written to the output file

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If too many sequences (>250) or invalid nucleotides found

Note:
    - Binary format includes factors, sequence IDs, sentinel indices, and footer metadata
    - Only A, C, T, G nucleotides are allowed (case insensitive)
    - This function overwrites the output file if it exists
    - Reverse complement matches are NOT supported during factorization
    - For single-threaded execution (num_threads=1), no temporary files are created
    - GIL is released during computation for better performance
)doc");

    m.def("parallel_write_factors_dna_w_reference_fasta_files_to_binary",
          &noLZSS::parallel_write_factors_dna_w_reference_fasta_files_to_binary,
          py::arg("reference_fasta_path"),
          py::arg("target_fasta_path"),
          py::arg("out_path"),
          py::arg("num_threads") = 0,
          R"doc(Parallel version of write_factors_dna_w_reference_fasta_files_to_binary.

Reads DNA sequences from reference and target FASTA files, concatenates them with
sentinels, and performs parallel factorization starting from target sequences,
writing results to a binary output file with metadata.

Args:
    reference_fasta_path: Path to FASTA file containing reference DNA sequences
    target_fasta_path: Path to FASTA file containing target DNA sequences
    out_path: Path to output file where binary factors will be written
    num_threads: Number of threads to use (0 = auto-detect based on input size, default)

Returns:
    Number of factors written to the output file

Raises:
    RuntimeError: If FASTA files cannot be opened or contain no valid sequences
    ValueError: If too many sequences total or invalid nucleotides found

Note:
    - Binary format includes factors, sequence IDs, sentinel indices, and footer metadata
    - Only A, C, T, G nucleotides are allowed (case insensitive)
    - This function overwrites the output file if it exists
    - Reverse complement matches are supported during factorization
    - Factorization starts from target sequence positions only
    - For single-threaded execution (num_threads=1), no temporary files are created
    - GIL is released during computation for better performance
)doc");

    // Per-sequence FASTA factorization (with reverse complement)
    py::class_<noLZSS::FastaPerSequenceFactorizationResult>(m, "FastaPerSequenceFactorizationResult", 
        "Result of per-sequence FASTA factorization containing factors for each sequence separately")
        .def_readonly("per_sequence_factors", &noLZSS::FastaPerSequenceFactorizationResult::per_sequence_factors, 
            "List of factor lists, one for each sequence")
        .def_readonly("sequence_ids", &noLZSS::FastaPerSequenceFactorizationResult::sequence_ids, 
            "List of sequence identifiers from FASTA headers");

    m.def("factorize_fasta_dna_w_rc_per_sequence", [](const std::string& fasta_path) {
        py::gil_scoped_release release;
        auto result = noLZSS::factorize_fasta_dna_w_rc_per_sequence(fasta_path);
        py::gil_scoped_acquire acquire;

        // Convert per-sequence factors to Python lists
        py::list per_seq_factors_list;
        for (const auto& seq_factors : result.per_sequence_factors) {
            py::list factors_list;
            for (const auto& f : seq_factors) {
                factors_list.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
            }
            per_seq_factors_list.append(factors_list);
        }

        // Convert sequence IDs to Python list
        py::list sequence_ids_list;
        for (const auto& seq_id : result.sequence_ids) {
            sequence_ids_list.append(seq_id);
        }

        return py::make_tuple(per_seq_factors_list, sequence_ids_list);
    }, py::arg("fasta_path"), R"doc(Factorize each DNA sequence in a FASTA file separately with reverse complement awareness.

Unlike factorize_fasta_multiple_dna_w_rc which concatenates sequences with sentinels,
this function factorizes each sequence independently. Each sequence gets its own
compressed suffix tree and factorization, which avoids sentinel limitations and
produces cleaner per-sequence results.

Args:
    fasta_path: Path to the FASTA file containing DNA sequences

Returns:
    Tuple of (per_sequence_factors, sequence_ids) where:
    - per_sequence_factors: List of factor lists, one for each sequence. Each factor is (start, length, ref, is_rc)
    - sequence_ids: List of sequence identifiers from FASTA headers

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If invalid nucleotides found in sequences

Note:
    - Only A, C, T, G nucleotides are allowed (case insensitive)
    - Sequences are converted to uppercase before factorization
    - Reverse complement matches are supported during factorization
    - Each sequence is factorized independently - no cross-sequence matches
    - ref field has RC_MASK cleared. is_rc boolean indicates if this was a reverse complement match
)doc");

    m.def("factorize_fasta_dna_no_rc_per_sequence", [](const std::string& fasta_path) {
        py::gil_scoped_release release;
        auto result = noLZSS::factorize_fasta_dna_no_rc_per_sequence(fasta_path);
        py::gil_scoped_acquire acquire;

        // Convert per-sequence factors to Python lists
        py::list per_seq_factors_list;
        for (const auto& seq_factors : result.per_sequence_factors) {
            py::list factors_list;
            for (const auto& f : seq_factors) {
                factors_list.append(py::make_tuple(f.start, f.length, f.ref & ~noLZSS::RC_MASK, noLZSS::is_rc(f.ref)));
            }
            per_seq_factors_list.append(factors_list);
        }

        // Convert sequence IDs to Python list
        py::list sequence_ids_list;
        for (const auto& seq_id : result.sequence_ids) {
            sequence_ids_list.append(seq_id);
        }

        return py::make_tuple(per_seq_factors_list, sequence_ids_list);
    }, py::arg("fasta_path"), R"doc(Factorize each DNA sequence in a FASTA file separately without reverse complement awareness.

Unlike factorize_fasta_multiple_dna_no_rc which concatenates sequences with sentinels,
this function factorizes each sequence independently. Each sequence gets its own
compressed suffix tree and factorization, which avoids sentinel limitations and
produces cleaner per-sequence results.

Args:
    fasta_path: Path to the FASTA file containing DNA sequences

Returns:
    Tuple of (per_sequence_factors, sequence_ids) where:
    - per_sequence_factors: List of factor lists, one for each sequence. Each factor is (start, length, ref, is_rc)
    - sequence_ids: List of sequence identifiers from FASTA headers

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If invalid nucleotides found in sequences

Note:
    - Only A, C, T, G nucleotides are allowed (case insensitive)
    - Sequences are converted to uppercase before factorization
    - Reverse complement matches are NOT supported during factorization
    - Each sequence is factorized independently - no cross-sequence matches
)doc");

    m.def("write_factors_binary_file_fasta_dna_w_rc_per_sequence", 
        [](const std::string& fasta_path, const std::string& out_dir) {
        py::gil_scoped_release release;
        size_t count = noLZSS::write_factors_binary_file_fasta_dna_w_rc_per_sequence(fasta_path, out_dir);
        py::gil_scoped_acquire acquire;
        return count;
    }, py::arg("fasta_path"), py::arg("out_dir"), R"doc(Write factors from per-sequence DNA factorization with reverse complement to separate binary files.

Reads a FASTA file, factorizes each sequence independently with reverse complement awareness,
and writes each sequence's factors to a separate binary file. File names include the sequence ID.

Args:
    fasta_path: Path to input FASTA file containing DNA sequences
    out_dir: Path to output directory where binary factor files will be written

Returns:
    int: Total number of factors written across all sequences

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If invalid nucleotides found

Note:
    - Creates separate binary file for each sequence: <out_dir>/<seq_id>.bin
    - Binary format per file: factors + metadata footer
    - Only A, C, T, G nucleotides are allowed (case insensitive)
    - Reverse complement matches are supported during factorization
)doc");

    m.def("write_factors_binary_file_fasta_dna_no_rc_per_sequence", 
        [](const std::string& fasta_path, const std::string& out_dir) {
        py::gil_scoped_release release;
        size_t count = noLZSS::write_factors_binary_file_fasta_dna_no_rc_per_sequence(fasta_path, out_dir);
        py::gil_scoped_acquire acquire;
        return count;
    }, py::arg("fasta_path"), py::arg("out_dir"), R"doc(Write factors from per-sequence DNA factorization without reverse complement to separate binary files.

Reads a FASTA file, factorizes each sequence independently without reverse complement awareness,
and writes each sequence's factors to a separate binary file. File names include the sequence ID.

Args:
    fasta_path: Path to input FASTA file containing DNA sequences
    out_dir: Path to output directory where binary factor files will be written

Returns:
    int: Total number of factors written across all sequences

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If invalid nucleotides found

Note:
    - Creates separate binary file for each sequence: <out_dir>/<seq_id>.bin
    - Binary format per file: factors + metadata footer
    - Only A, C, T, G nucleotides are allowed (case insensitive)
    - Reverse complement matches are NOT supported during factorization
)doc");

    m.def("count_factors_fasta_dna_w_rc_per_sequence", [](const std::string& fasta_path) {
        py::gil_scoped_release release;
        auto cpp_result = noLZSS::count_factors_fasta_dna_w_rc_per_sequence(fasta_path);
        py::gil_scoped_acquire acquire;

        py::list counts_list;
        for (size_t count : cpp_result.factor_counts) {
            counts_list.append(py::int_(count));
        }

        py::list sequence_ids_list;
        for (const auto& seq_id : cpp_result.sequence_ids) {
            sequence_ids_list.append(py::str(seq_id));
        }

        return py::make_tuple(counts_list, sequence_ids_list, cpp_result.total_factors);
    }, py::arg("fasta_path"), R"doc(Count per-sequence factors from DNA FASTA (with reverse complement).

Reads a FASTA file and factorizes each sequence independently with reverse complement awareness,
returning per-sequence counts, sequence IDs, and the aggregate total without storing factors.

Args:
    fasta_path: Path to the FASTA file containing DNA sequences

Returns:
    Tuple[List[int], List[str], int]: (per-sequence counts, sequence IDs, total factors)

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If invalid nucleotides found

Note:
    - Memory-efficient - only counts factors without storing them
    - Only A, C, T, G nucleotides are allowed (case insensitive)
)doc");

    m.def("count_factors_fasta_dna_no_rc_per_sequence", [](const std::string& fasta_path) {
        py::gil_scoped_release release;
        auto cpp_result = noLZSS::count_factors_fasta_dna_no_rc_per_sequence(fasta_path);
        py::gil_scoped_acquire acquire;

        py::list counts_list;
        for (size_t count : cpp_result.factor_counts) {
            counts_list.append(py::int_(count));
        }

        py::list sequence_ids_list;
        for (const auto& seq_id : cpp_result.sequence_ids) {
            sequence_ids_list.append(py::str(seq_id));
        }

        return py::make_tuple(counts_list, sequence_ids_list, cpp_result.total_factors);
    }, py::arg("fasta_path"), R"doc(Count per-sequence factors from DNA FASTA (without reverse complement).

Reads a FASTA file and factorizes each sequence independently without reverse complement awareness,
returning per-sequence counts, sequence IDs, and the aggregate total without storing factors.

Args:
    fasta_path: Path to the FASTA file containing DNA sequences

Returns:
    Tuple[List[int], List[str], int]: (per-sequence counts, sequence IDs, total factors)

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If invalid nucleotides found

Note:
    - Memory-efficient - only counts factors without storing them
    - Only A, C, T, G nucleotides are allowed (case insensitive)
)doc");

    m.def("parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence", 
        [](const std::string& fasta_path, const std::string& out_dir, size_t num_threads) {
        py::gil_scoped_release release;
        size_t count = noLZSS::parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence(fasta_path, out_dir, num_threads);
        py::gil_scoped_acquire acquire;
        return count;
    }, py::arg("fasta_path"), py::arg("out_dir"), py::arg("num_threads") = 0, R"doc(Parallel version of write_factors_binary_file_fasta_dna_w_rc_per_sequence.

Reads a FASTA file, factorizes each sequence independently with reverse complement
awareness using parallel processing, and writes each sequence to a separate binary file.

Args:
    fasta_path: Path to input FASTA file containing DNA sequences
    out_dir: Path to output directory where binary factor files will be written
    num_threads: Number of threads to use (0 = auto-detect based on sequence count)

Returns:
    int: Total number of factors written across all sequences

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If invalid nucleotides found

Note:
    - Each sequence is factorized independently in parallel
    - Creates separate binary file for each sequence: <out_dir>/<seq_id>.bin
    - Binary format per file: factors + metadata footer
    - Only A, C, T, G nucleotides are allowed (case insensitive)
    - Reverse complement matches are supported during factorization
    - GIL is released during computation for better performance
)doc");

    m.def("parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence", 
        [](const std::string& fasta_path, const std::string& out_dir, size_t num_threads) {
        py::gil_scoped_release release;
        size_t count = noLZSS::parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence(fasta_path, out_dir, num_threads);
        py::gil_scoped_acquire acquire;
        return count;
    }, py::arg("fasta_path"), py::arg("out_dir"), py::arg("num_threads") = 0, R"doc(Parallel version of write_factors_binary_file_fasta_dna_no_rc_per_sequence.

Reads a FASTA file, factorizes each sequence independently without reverse complement
awareness using parallel processing, and writes each sequence to a separate binary file.

Args:
    fasta_path: Path to input FASTA file containing DNA sequences
    out_dir: Path to output directory where binary factor files will be written
    num_threads: Number of threads to use (0 = auto-detect based on sequence count)

Returns:
    int: Total number of factors written across all sequences

Raises:
    RuntimeError: If FASTA file cannot be opened or contains no valid sequences
    ValueError: If invalid nucleotides found

Note:
    - Each sequence is factorized independently in parallel
    - Creates separate binary file for each sequence: <out_dir>/<seq_id>.bin
    - Binary format per file: factors + metadata footer
    - Only A, C, T, G nucleotides are allowed (case insensitive)
    - Reverse complement matches are NOT supported during factorization
    - GIL is released during computation for better performance
)doc");

    // Version information
#ifdef NOLZSS_VERSION
    m.attr("__version__") = NOLZSS_VERSION;
#else
    m.attr("__version__") = "0.0.0";
#endif
}
