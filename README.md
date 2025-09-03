# noLZSS

[![Build Wheels](https://github.com/OmerKerner/noLZSS/actions/workflows/wheels.yml/badge.svg)](https://github.com/OmerKerner/noLZSS/actions/workflows/wheels.yml)
<img align="right" src="assets/logo.png" alt="noLZSS Logo" width=200px/>

**Non-overlapping Lempel–Ziv–Storer–Szymanski factorization**

High-performance Python library for text factorization using compressed suffix trees. The library provides efficient algorithms for finding non-overlapping factors in text data, with both in-memory and file-based processing capabilities. Based on a paper by Dominik Köppl - [Non-Overlapping LZ77 Factorization and LZ78 Substring Compression Queries with Suffix Trees](https://doi.org/10.3390/a14020044)

## Features

- 🚀 **High Performance**: Uses compressed suffix trees (SDSL) for optimal factorization speed
- 💾 **Memory Efficient**: File-based processing for large datasets without loading everything into memory
- 🐍 **Python Bindings**: Easy-to-use Python interface with proper GIL management
- 📊 **Multiple Output Formats**: Get factors as lists, counts, or binary files
- 🔧 **Flexible API**: Support for both strings and files with optional performance hints

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
```

## Algorithm Details

The library implements the **Non-overlapping Lempel-Ziv-Storer-Szymanski (LZSS)** factorization algorithm using:

- **Compressed Suffix Trees**: Built using the SDSL (Succinct Data Structure Library)
- **Range Minimum Queries**: For efficient lowest common ancestor computations
- **Sink-based Processing**: Memory-efficient processing using callback functions


## Performance

- **Time Complexity**: 𝒪(𝑛 lg<sup>ϵ</sup> 𝑛) for factorization, where n is input length, and 𝜖 ∈ (0,1]
- **Space Complexity**: 𝒪(𝑛lg𝜎) for suffix tree construction, where 𝜎 is the alphabet size
- **Memory Usage**: File-based processing uses minimal memory for large files

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
