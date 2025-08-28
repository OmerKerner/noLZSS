# noLZSS

[![Build Wheels](https://github.com/OmerKerner/noLZSS/actions/workflows/wheels.yml/badge.svg)](https://github.com/OmerKerner/noLZSS/actions/workflows/wheels.yml)

**Non-overlapping Lempel–Ziv–Storer–Szymanski factorization**

<img alighn="right" src="assets/logo.png" alt="noLZSS Logo" width="200" height="200"/>

High-performance Python library for text factorization using compressed suffix trees. The library provides efficient algorithms for finding non-overlapping factors in text data, with both in-memory and file-based processing capabilities.

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

# Factorize a text string (must end with '$' sentinel)
text = b"abracadabra$"
factors = noLZSS.factorize(text)
print(factors)  # [(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (7, 1), (8, 1), (9, 1), (10, 1)]
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
# Use reserve hint for better performance
factors = noLZSS.factorize_file("data.txt", reserve_hint=1000)

# Process factors efficiently
for start, length in factors:
    substring = text[start:start+length]
    print(f"Factor at {start}: '{substring}' (length {length})")
```

## API Reference

### Core Functions

#### `factorize(data)`
Factorize in-memory text into LZSS factors.

**Parameters:**
- `data` (bytes-like): Input text that must end with '$' sentinel character

**Returns:**
- `List[Tuple[int, int]]`: List of (start, length) tuples representing factors

**Example:**
```python
factors = noLZSS.factorize(b"hello$")
# Returns: [(0, 1), (1, 1), (2, 1), (3, 1), (4, 1)]
```

#### `factorize_file(path, reserve_hint=0)`
Factorize text from a file into LZSS factors.

**Parameters:**
- `path` (str): Path to input file containing text that must end with '$' sentinel
- `reserve_hint` (int, optional): Hint for reserving space in output vector

**Returns:**
- `List[Tuple[int, int]]`: List of (start, length) tuples representing factors

#### `count_factors(data)`
Count LZSS factors in text without storing them.

**Parameters:**
- `data` (bytes-like): Input text that must end with '$' sentinel character

**Returns:**
- `int`: Number of factors in the factorization

#### `count_factors_file(path)`
Count LZSS factors in a file without storing them.

**Parameters:**
- `path` (str): Path to input file containing text that must end with '$' sentinel

**Returns:**
- `int`: Number of factors in the factorization

#### `write_factors_binary_file(in_path, out_path, assume_has_sentinel=False)`
Write LZSS factors from input file to binary output file.

**Parameters:**
- `in_path` (str): Path to input file containing text that must end with '$' sentinel
- `out_path` (str): Path to output file for binary factors
- `assume_has_sentinel` (bool, optional): Unused parameter for API consistency

**Returns:**
- `int`: Number of factors written

### Important Notes

⚠️ **Sentinel Requirement**: All input strings and files **must** end with the '$' character. This is essential for correct factorization results.

```python
# Correct usage
text = b"hello$"
factors = noLZSS.factorize(text)

# Incorrect usage (will produce wrong results)
text = b"hello"  # Missing '$' sentinel
```

## Algorithm Details

The library implements the **Non-overlapping Lempel-Ziv-Storer-Szymanski (LZSS)** factorization algorithm using:

- **Compressed Suffix Trees**: Built using the SDSL (Succinct Data Structure Library)
- **Range Minimum Queries**: For efficient lowest common ancestor computations
- **Sink-based Processing**: Memory-efficient processing using callback functions

### Factor Properties

- **Non-overlapping**: Factors don't overlap in the original text
- **Complete Coverage**: Factors cover the entire input text (excluding sentinel)
- **Greedy**: Each factor is the longest possible match at its position

## Performance

- **Time Complexity**: O(n) for factorization, where n is input length
- **Space Complexity**: O(n) for suffix tree construction
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

### Project Structure

```
src/noLZSS/
├── factorizer.hpp      # C++ header with function declarations
├── factorizer.cpp      # C++ implementation
├── bindings.cpp        # Python bindings
├── version.hpp         # Version information
└── __init__.py         # Python package initialization

tests/
└── test_basic.py       # Unit tests
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

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
