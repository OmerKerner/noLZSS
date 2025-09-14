# noLZSS

[![Build Wheels](https://github.com/OmerKerner/noLZSS/actions/workflows/wheels.yml/badge.svg)](https://github.com/OmerKerner/noLZSS/actions/workflows/wheels.yml)
[![Documentation](https://github.com/OmerKerner/noLZSS/actions/workflows/docs.yml/badge.svg)](https://omerkerner.github.io/noLZSS/)
<img align="right" src="assets/logo.png" alt="noLZSS Logo" width=200px/>

**Non-overlapping Lempel–Ziv–Storer–Szymanski factorization**

High-performance Python library for text factorization using compressed suffix trees. The library provides efficient algorithms for finding non-overlapping factors in text data, with both in-memory and file-based processing capabilities. Based on a paper by Dominik Köppl - [Non-Overlapping LZ77 Factorization and LZ78 Substring Compression Queries with Suffix Trees](https://doi.org/10.3390/a14020044)

## Features

- 🚀 **High Performance**: Uses compressed suffix trees (SDSL) for optimal factorization speed
- 💾 **Memory Efficient**: File-based processing for large datasets without loading everything into memory
- 🐍 **Python Bindings**: Easy-to-use Python interface with proper GIL management
- 📊 **Multiple Output Formats**: Get factors as lists, counts, or binary files
- 🔧 **Flexible API**: Support for both strings and files with optional performance hints
- 🧬 **Genomics Support**: Specialized functions for FASTA file processing of nucleotide and protein sequences
- ⚡ **C++ Extensions**: High-performance C++ implementations for memory-intensive operations

## Installation

### From Source (Development)

```bash
git clone https://github.com/OmerKerner/noLZSS.git
cd noLZSS
pip install -e .
```

### Requirements

- Python 3.8+
- C++17 compatible compiler
- CMake 3.20+

## Quick Start

### Basic Usage

```python
import noLZSS

# Factorize a text string
text = b"abracadabra"
factors = noLZSS.factorize(text)
print(factors)  # [(0, 1, 0), (1, 1, 1), (2, 1, 2), (3, 1, 0), (4, 1, 4), (5, 1, 0), (6, 1, 6), (7, 4, 0)]
```

### Working with Files

```python
# Factorize text from a file
factors = noLZSS.factorize_file("large_text.txt")
print(f"Found {len(factors)} factors")

# Just count factors without storing them (memory efficient)
count = noLZSS.count_factors_file("large_text.txt")
print(f"Total factors: {count}")

# Write factors to binary file for later processing
noLZSS.write_factors_binary_file("input.txt", "factors.bin")
```

### Advanced Usage

```python
# Use reserve hint for better performance - An estimate of the number of compressed factors
factors = noLZSS.factorize_file("data.txt", reserve_hint=1000)

# Process factors efficiently
for start, length, ref in factors:
    substring = text[start:start+length]
    print(f"Factor at {start}: '{substring}' (length {length}, ref {ref})")
```

## API Reference

### Core Functions

#### `factorize(data)`
Factorize in-memory text into LZSS factors.

**Parameters:**
- `data` (bytes-like): Input text to factorize

**Returns:**
- `List[Tuple[int, int, int]]`: List of (start, length, reference) tuples representing factors

**Example:**
```python
factors = noLZSS.factorize(b"hello")
# Returns: [(0, 1, 0), (1, 1, 1), (2, 1, 2), (3, 1, 2), (4, 1, 4)]
```

#### `factorize_file(path, reserve_hint=0)`
Factorize text from a file into LZSS factors.

**Parameters:**
- `path` (str): Path to input file containing text
- `reserve_hint` (int, optional): Hint for reserving space in output vector

**Returns:**
- `List[Tuple[int, int, int]]`: List of (start, length, reference) tuples representing factors

#### `count_factors(data)`
Count LZSS factors in text without storing them.

**Parameters:**
- `data` (bytes-like): Input text to analyze

**Returns:**
- `int`: Number of factors in the factorization

#### `count_factors_file(path)`
Count LZSS factors in a file without storing them.

**Parameters:**
- `path` (str): Path to input file containing text

**Returns:**
- `int`: Number of factors in the factorization

#### `write_factors_binary_file(in_path, out_path)`
Write LZSS factors from input file to binary output file.

**Parameters:**
- `in_path` (str): Path to input file containing text
- `out_path` (str): Path to output file for binary factors

**Returns:**
- `int`: Number of factors written

### Genomics Functions

#### `read_nucleotide_fasta(filepath)`
Read and factorize nucleotide sequences from a FASTA file.

Only accepts sequences containing A, C, T, G (case insensitive, converted to uppercase).
Sequences are validated and factorized using LZSS.

**Parameters:**
- `filepath` (str or Path): Path to FASTA file

**Returns:**
- `List[Tuple[str, List[Tuple[int, int, int]]]]`: List of (sequence_id, factors) tuples

**Raises:**
- `FASTAError`: If file format is invalid or contains invalid nucleotides

**Example:**
```python
results = noLZSS.read_nucleotide_fasta("sequences.fasta")
for seq_id, factors in results:
    print(f"Sequence {seq_id}: {len(factors)} factors")
```

#### `read_protein_fasta(filepath)`
Read amino acid sequences from a FASTA file.

Only accepts sequences containing canonical amino acids (case insensitive, converted to uppercase).

**Parameters:**
- `filepath` (str or Path): Path to FASTA file

**Returns:**
- `List[Tuple[str, str]]`: List of (sequence_id, sequence) tuples

**Raises:**
- `FASTAError`: If file format is invalid or contains invalid amino acids

#### `read_fasta_auto(filepath)`
Automatically detect sequence type and read FASTA file accordingly.

For nucleotide sequences: validates A,C,T,G and returns factorized results
For amino acid sequences: validates canonical amino acids and returns sequences

**Parameters:**
- `filepath` (str or Path): Path to FASTA file

**Returns:**
- For nucleotides: `List[Tuple[str, List[Tuple[int, int, int]]]]`
- For proteins: `List[Tuple[str, str]]`

**Raises:**
- `FASTAError`: If file format is invalid or sequence type cannot be determined

**Example:**
```python
# Automatically handles both nucleotide and protein FASTA files
results = noLZSS.read_fasta_auto("mixed_sequences.fasta")
```

#### `process_fasta_with_plots(filepath, output_dir, max_sequences=None)`
Process a FASTA file and generate factors and plots for all sequences.

For each sequence, this function:
- Factorizes the sequence using LZSS
- Saves factor data as text files
- Optionally saves factors as binary files
- Creates and saves plots of factor length accumulation

**Parameters:**
- `filepath` (str or Path): Path to FASTA file
- `output_dir` (str or Path): Directory to save all output files
- `max_sequences` (int, optional): Maximum number of sequences to process
- `save_factors_text` (bool): Whether to save factors as text files (default: True)
- `save_factors_binary` (bool): Whether to save factors as binary files (default: False)

**Returns:**
- `dict`: Processing results for each sequence

**Example:**
```python
# Process all sequences in a FASTA file
results = noLZSS.process_fasta_with_plots(
    "genome.fasta", 
    "output_factors/",
    max_sequences=10  # Process only first 10 sequences
)
```

#### `process_nucleotide_fasta(filepath)`
High-performance C++ implementation for processing nucleotide FASTA files.

This function provides memory-efficient processing of large FASTA files by directly reading
and concatenating sequences in C++ without creating intermediate Python objects. It's ideal
for processing large genome files that would otherwise cause memory issues.

**Parameters:**
- `filepath` (str): Path to FASTA file containing nucleotide sequences

**Returns:**
- `dict`: Dictionary containing:
  - `"sequence"`: Concatenated sequences with sentinels
  - `"num_sequences"`: Number of sequences processed
  - `"sequence_ids"`: List of sequence IDs
  - `"sequence_lengths"`: List of sequence lengths
  - `"sequence_positions"`: List of sequence start positions

**Raises:**
- `RuntimeError`: If file cannot be read, contains invalid nucleotides, or has >251 sequences

**Example:**
```python
from noLZSS._noLZSS import process_nucleotide_fasta

# Process large FASTA file efficiently
result = process_nucleotide_fasta("large_genome.fasta")
print(f"Processed {result['num_sequences']} sequences")
print(f"Total length: {len(result['sequence']):,} characters")

# Use the concatenated sequence for factorization
factors = noLZSS.factorize(result["sequence"])
```

**Note:** This function uses sentinels (characters 1-251) to separate sequences, avoiding conflicts with nucleotides A, C, G, T.

#### `process_amino_acid_fasta(filepath)`
High-performance C++ implementation for processing amino acid FASTA files.

This function provides memory-efficient processing of large FASTA files containing amino acid sequences by directly reading and concatenating sequences in C++ without creating intermediate Python objects. It's ideal for processing large protein databases that would otherwise cause memory issues.

**Parameters:**
- `filepath` (str): Path to FASTA file containing amino acid sequences

**Returns:**
- `dict`: Dictionary containing:
  - `"sequence"`: Concatenated sequences with sentinels
  - `"num_sequences"`: Number of sequences processed
  - `"sequence_ids"`: List of sequence IDs
  - `"sequence_lengths"`: List of sequence lengths
  - `"sequence_positions"`: List of sequence start positions

**Raises:**
- `RuntimeError`: If file cannot be read, contains invalid amino acids, or has more than 235 sequences

**Example:**
```python
from noLZSS._noLZSS import process_amino_acid_fasta

# Process large protein FASTA file efficiently
result = process_amino_acid_fasta("proteins.fasta")
print(f"Processed {result['num_sequences']} protein sequences")
print(f"Total length: {len(result['sequence']):,} characters")

# Use the concatenated sequence for factorization
factors = noLZSS.factorize(result["sequence"])
```

**Note:** Canonical amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y. Sentinels avoid these characters and null (0).

### Genomics Usage Guide

The library provides specialized functions for biological sequence analysis:

#### Choosing the Right Function

- **For small FASTA files (< 100MB)**: Use `read_nucleotide_fasta()` or `read_protein_fasta()`
- **For large FASTA files (> 100MB)**: Use `process_nucleotide_fasta()` (C++ implementation)
- **For mixed sequence types**: Use `read_fasta_auto()` for automatic detection
- **For batch processing with visualization**: Use `process_fasta_with_plots()`

#### Memory Considerations

```python
# Memory-efficient processing of large genomes
from noLZSS._noLZSS import process_nucleotide_fasta

# This avoids loading the entire file into Python memory
result = process_nucleotide_fasta("human_genome.fasta")
factors = noLZSS.factorize(result["sequence"])
```

#### Sequence Validation

All genomics functions perform strict validation:
- **Nucleotides**: Only A, C, T, G allowed (case insensitive)
- **Proteins**: Only canonical amino acids allowed (case insensitive)
- **Auto-detection**: Based on character composition analysis

### Usage Notes

```python
# Simple usage
text = b"hello"
factors = noLZSS.factorize(text)

# With validation and analysis
import noLZSS
result = noLZSS.factorize_with_info("hello world")
print(f"Found {result['num_factors']} factors")
print(f"Text entropy: {result['alphabet_info']['entropy']:.2f} bits")

# Read factors from binary file
factors = noLZSS.read_factors_binary_file("factors.bin")

# Plot factor length accumulation
noLZSS.plot_factor_lengths(factors, save_path="plot.png")
```

### Utility Functions

#### `analyze_alphabet(data)`
Analyze the alphabet and entropy of input text.

**Parameters:**
- `data` (str or bytes): Input text to analyze

**Returns:**
- `dict`: Dictionary with alphabet statistics including size, entropy, and character frequencies

#### `read_factors_binary_file(path)`
Read LZSS factors from a binary file created by `write_factors_binary_file`.

**Parameters:**
- `path` (str): Path to binary factors file

**Returns:**
- `List[Tuple[int, int, int]]`: List of (start, length, reference) tuples

#### `plot_factor_lengths(factors_or_file, save_path=None, show_plot=True)`
Plot the cumulative factor lengths vs factor index.

Creates a scatter plot showing how factor lengths accumulate:
- X-axis: Cumulative sum of factor lengths
- Y-axis: Factor index (0-based)

**Parameters:**
- `factors_or_file`: Either a list of factors or path to binary factors file
- `save_path` (str, optional): Path to save the plot image
- `show_plot` (bool): Whether to display the plot (default: True)

**Requires:** `pip install matplotlib` or `pip install noLZSS[plotting]`

## Algorithm Details

The library implements the **Non-overlapping Lempel-Ziv-Storer-Szymanski (LZSS)** factorization algorithm using:

- **Compressed Suffix Trees**: Built using the SDSL (Succinct Data Structure Library)
- **Range Minimum Queries**: For efficient lowest common ancestor computations
- **Sink-based Processing**: Memory-efficient processing using callback functions


## Performance

- **Time Complexity**: 𝒪(𝑛 lg<sup>ϵ</sup> 𝑛) for factorization, where n is input length, and 𝜖 ∈ (0,1]
- **Space Complexity**: 𝒪(𝑛lg𝜎) for suffix tree construction, where 𝜎 is the alphabet size
- **Memory Usage**: File-based processing uses minimal memory for large files
- **C++ Extensions**: Specialized high-performance functions for memory-intensive genomics operations

### Performance Tips

```python
# For large files, use file-based functions
factors = noLZSS.factorize_file("large_file.txt", reserve_hint=1000000)

# For genomics, use C++ implementation for large FASTA files
from noLZSS._noLZSS import process_nucleotide_fasta
result = process_nucleotide_fasta("genome.fasta")  # Memory efficient

# Use reserve_hint for better performance when you know factor count
factors = noLZSS.factorize_file("data.txt", reserve_hint=50000)
```

## Documentation

Complete documentation is available at **[omerkerner.github.io/noLZSS](https://omerkerner.github.io/noLZSS/)**

The documentation includes:
- **Python API Reference**: Complete Python API with examples and parameter descriptions
- **C++ API Reference**: Auto-generated C++ API documentation from source code
- **Genomics Module**: Specialized functions for biological sequence analysis
- **Examples and Tutorials**: Comprehensive usage examples and best practices

## Development

### Building from Source

```bash
# Clone repository
git clone https://github.com/OmerKerner/noLZSS.git
cd noLZSS

# Install in development mode
pip install -e .

# Run tests
python -m pytest tests/
```

## License

This project is licensed under the BSD 3-Clause License (see `LICENSE`).

The repository vendors third-party components (notably SDSL v3). Third-party license texts and attribution are provided in `THIRD_PARTY_LICENSES.txt`.

## Citation

If you use this library in your research, please cite:

```bibtex
@software{noLZSS,
  title = {noLZSS: Non-overlapping Lempel-Ziv-Storer-Szymanski factorization},
  author = {Kerner, Omer},
  url = {https://github.com/OmerKerner/noLZSS},
  year = {2024}
}
```
