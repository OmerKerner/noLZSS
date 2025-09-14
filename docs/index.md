# noLZSS Documentation

Welcome to the noLZSS documentation! This package provides high-performance **Non-overlapping Lempel-Ziv-Storer-Szymanski (LZSS) factorization** with a C++ core and Python bindings.

## What is noLZSS?

noLZSS computes non-overlapping LZ factorizations of strings and files, particularly optimized for genomics applications. It uses compressed suffix trees (SDSL v3) for efficient computation and provides both Python and C++ APIs.

## Quick Start

```python
import noLZSS

# Factorize a string
factors = noLZSS.factorize("abcabcabc")
print(factors)  # [(0, 1, 0), (1, 1, 1), (2, 1, 2), (3, 3, 0), (6, 3, 0)]

# Factorize a file
factors = noLZSS.factorize_file("input.txt")

# Count factors without storing them
count = noLZSS.count_factors("large_string")
```

## Installation

```bash
pip install noLZSS
```

For development:
```bash
git clone https://github.com/OmerKerner/noLZSS.git
cd noLZSS
pip install -e .
```

## Key Features

- **High Performance**: C++ core with compressed suffix trees
- **Python Integration**: Clean Python API with comprehensive error handling
- **Genomics Support**: Specialized functions for DNA/protein sequences
- **Memory Efficient**: File-based processing for large datasets
- **Cross-platform**: Works on Linux and macOS

## API Documentation

```{toctree}
:maxdepth: 2

python_api
cpp_api
genomics
examples
```

## Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`