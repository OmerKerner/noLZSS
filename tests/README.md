# noLZSS Test Suite

This directory contains comprehensive tests for the noLZSS package.

## Test Files

### Core Functionality Tests

- **`test_cpp_bindings.py`**
  - Tests the C++ extension bindings directly
  - Includes invariant checking and consistency tests
  - Requires the C++ extension to be built

### Python Layer Tests

- **`test_utils.py`**
  - Comprehensive tests for utility functions
  - Input validation, alphabet analysis, sequence detection
  - File I/O helpers and error handling

- **`test_core.py`**
  - Tests the Python wrapper functions
  - Input validation integration with C++ calls
  - Error handling and type checking

### Module Structure Tests

- **`test_structure.py`**
  - Basic import and structure tests
  - Verifies the modular package organization

- **`test_genomics.py`**
  - Tests the genomics subpackage structure
  - Placeholder tests for future FASTA/sequence functionality

### Integration Tests

- **`test_integration.py`**
  - End-to-end workflow tests
  - Cross-module integration verification
  - File operations and error handling

## Test Runner

- **`run_all_tests.py`**
  - Executes all test files in the correct order
  - Provides summary of test results

## Running Tests

### Run All Tests
```bash
python tests/run_all_tests.py
```

### Run Individual Test Files
```bash
# Utils functionality (works without C++ build)
python tests/test_utils.py

# Package structure
python tests/test_structure.py

# Integration tests
python tests/test_integration.py

# Genomics subpackage
python tests/test_genomics.py

# Core wrappers (may fail without C++ build)
python tests/test_core.py

# C++ bindings (requires build)
python tests/test_cpp_bindings.py
```

### Using pytest (if available)
```bash
pytest tests/
```
