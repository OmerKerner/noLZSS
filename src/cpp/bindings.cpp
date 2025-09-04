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
#include "version.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_noLZSS, m) {
    m.doc() = "Non-overlapping Lempel-Ziv-Storer-Szymanski factorization\n\n"
              "This module provides efficient text factorization using compressed suffix trees.";

    // Factor class documentation
    py::class_<noLZSS::Factor>(m, "Factor", "Represents a single factorization factor with start position, length, and reference position")
        .def_readonly("start", &noLZSS::Factor::start, "Starting position of the factor in the original text")
        .def_readonly("length", &noLZSS::Factor::length, "Length of the factor substring")
        .def_readonly("ref", &noLZSS::Factor::ref, "Reference position of the previous occurrence");

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
    }, py::arg("data"), R"doc(Factorize a text string into LZSS factors.

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
    }, py::arg("path"), py::arg("reserve_hint") = 0, R"doc(Factorize text from file into LZSS factors.

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
    }, py::arg("path"), R"doc(Count number of LZSS factors in a file.

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
    }, py::arg("in_path"), py::arg("out_path"), R"doc(Write LZSS factors from file to binary output file.

Reads text from an input file, performs factorization, and writes the factors
in binary format to an output file. Each factor is written as two uint64_t values.

Args:
    in_path: Path to input file containing text
    out_path: Path to output file where binary factors will be written

Returns:
    Number of factors written to the output file

Note:
    Binary format: each factor is 24 bytes (3 × uint64_t: start, length, ref).
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

    // Version information
    m.attr("__version__") = std::to_string(noLZSS::VERSION_MAJOR) + "." + std::to_string(noLZSS::VERSION_MINOR) + "." + std::to_string(noLZSS::VERSION_PATCH);
}
