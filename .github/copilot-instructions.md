# noLZSS Copilot Instructions

## Project Overview

This is a high-performance C++/Python hybrid library for **Non-overlapping Lempel-Ziv-Storer-Szymanski (LZSS) factorization**. The library is specialized for genomics applications with compressed suffix trees (SDSL v3) as the core algorithm.

## Architecture Patterns

### Hybrid C++/Python Design
- **C++ Core**: High-performance algorithms in `src/cpp/` using SDSL v3 and pybind11
- **Python Wrappers**: Clean APIs in `src/noLZSS/` with validation and error handling
- **Build System**: scikit-build-core + CMake for seamless Python package building

### Key Components
```
src/cpp/                    # C++ implementation
├── factorizer.cpp/.hpp     # Core LZSS algorithm
├── fasta_processor.cpp/.hpp # Genomics-specific processing
└── bindings.cpp            # pybind11 Python bindings

src/noLZSS/                 # Python package
├── core.py                 # Main API wrappers
├── utils.py                # Input validation & utilities
└── genomics/               # Specialized genomics functions
```

## Development Workflows

### Building from Source
```bash
pip install -e .  # Triggers CMake build via scikit-build-core
```

### Testing Strategy
- **Modular tests**: `tests/test_*.py` files for different components
- **C++ dependency awareness**: Tests gracefully handle missing C++ extension
- **Run all tests**: `python tests/run_all_tests.py` or `python -m pytest tests/`

### Performance Testing
- Benchmarks in `benchmarks/fasta_benchmark.py` with memory tracking
- Benchmark results saved in `benchmarks/fasta_results/`

## Project-Specific Conventions

### Error Handling Pattern
Python wrappers add validation while preserving C++ performance:
```python
def factorize(data, validate=True):
    if validate:
        data = validate_input(data)  # Python validation
    return _factorize(data)  # Direct C++ call
```

### Memory Management
- **File-based processing**: Use `*_file()` functions for large datasets
- **Reserve hints**: Pass `reserve_hint` parameter for known factor counts
- **C++ implementations**: Direct memory-efficient processing (`process_nucleotide_fasta`)

### Genomics Specialization
- **Sequence validation**: Strict nucleotide (A,C,T,G) and amino acid validation
- **Sentinel handling**: Characters 1-251 as separators (avoiding nucleotide conflicts)
- **Reverse complement**: MSB flag in reference field for DNA analysis

## External Dependencies

### Critical Integration Points
- **SDSL v3 (xxsds fork)**: Vendored via FetchContent, provides compressed suffix trees
- **pybind11**: Python/C++ binding layer, also vendored if needed
- **matplotlib**: Optional dependency for plotting (`[plotting]` extra)

### Build Configuration
- **CMake 3.20+**: Required for FetchContent and C++17 support  
- **C++17**: Standard requirement across codebase
- **Position Independent Code**: Required for Python extensions

## Code Patterns to Follow

### Input Validation
Always validate in Python layer before calling C++:
```python
data = validate_input(data)  # Handles str->bytes, empty checks
```

### File Path Handling
Use `pathlib.Path` consistently and convert to string for C++:
```python
filepath = Path(filepath)
if not filepath.exists():
    raise FileNotFoundError(f"File not found: {filepath}")
return _cpp_function(str(filepath))
```

### Performance-Critical Paths
For large datasets, prefer C++ implementations:
- Use `process_nucleotide_fasta()` over Python `read_nucleotide_fasta()`
- Use `*_file()` functions instead of loading into memory
- Pass `reserve_hint` when factor count is known

## Testing Guidelines

### C++ Extension Awareness
Tests should handle missing C++ extension gracefully:
```python
try:
    from ._noLZSS import factorize
except ImportError:
    pytest.skip("C++ extension not available")
```

### Test Organization
- `test_utils.py`: Input validation, utilities (no C++ required)
- `test_core.py`: Python wrappers (requires C++ extension)  
- `test_cpp_bindings.py`: Direct C++ binding tests
- `test_genomics.py`: FASTA processing and genomics features

## Key Files for Understanding

- `pyproject.toml`: Build configuration and dependencies
- `CMakeLists.txt`: C++ build setup and SDSL integration
- `src/cpp/factorizer.hpp`: Core algorithm interface
- `src/noLZSS/core.py`: Main Python API patterns
- `tests/run_all_tests.py`: Complete test execution workflow