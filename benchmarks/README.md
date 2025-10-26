# noLZSS Benchmarking Suite

This directory contains comprehensive benchmarking tools for all noLZSS factorization functions, including core functions, DNA-specific functions, FASTA processing, and parallel computing implementations.

## Prerequisites

**Important**: Before running any benchmarks, you must build and install the noLZSS package:

```bash
# From the repository root
pip install -e .
```

This builds the C++ extension which is required for all benchmarking functions. If you encounter errors about missing functions, rebuild with:

```bash
pip install -e . --no-build-isolation --force-reinstall
```

## Quick Start

### Run All Benchmarks

```bash
# Run all benchmarks with default settings
python benchmarks/run_all_benchmarks.py --all

# Run quick test (smaller sizes, fewer runs)
python benchmarks/run_all_benchmarks.py --all --quick

# Run specific benchmarks
python benchmarks/run_all_benchmarks.py --core --dna --parallel
```

### Run Individual Benchmarks

```bash
# Core factorization functions
python benchmarks/core_benchmark.py

# DNA factorization with reverse complement
python benchmarks/dna_benchmark.py

# FASTA file processing
python benchmarks/fasta_benchmark.py

# Parallel factorization
python benchmarks/parallel_benchmark.py
```

### Use Predictions for Resource Estimation

```bash
# Predict resources for any FASTA input size
python benchmarks/fasta_predictor.py benchmarks/fasta_results/trend_parameters.pkl --size 5000000
```

## Benchmark Suites

### 1. Core Factorization Benchmark (`core_benchmark.py`)

Benchmarks the core factorization functions for general text processing.

**Functions benchmarked:**
- `factorize()` - Basic string factorization
- `factorize_file()` - File-based factorization
- `count_factors()` - Factor counting only
- `count_factors_file()` - File-based factor counting
- `write_factors_binary_file()` - Binary file output

**Usage:**
```bash
# Default: 1kb to 1Mb, 10 sizes
python benchmarks/core_benchmark.py

# Custom size range
python benchmarks/core_benchmark.py --min-size 1000 --max-size 1000000 --num-sizes 10

# Specific sizes
python benchmarks/core_benchmark.py --custom-sizes 1000 10000 100000 1000000

# More runs for better statistics
python benchmarks/core_benchmark.py --runs 5
```

**Output:**
- Time and memory benchmarks
- Throughput analysis
- Binary file size analysis
- Log-log plots with trend lines
- Trend parameters (JSON and pickle)

### 2. DNA Factorization Benchmark (`dna_benchmark.py`)

Benchmarks DNA-specific factorization functions with reverse complement support.

**Functions benchmarked:**
- `factorize_dna_w_rc()` - DNA factorization with reverse complement
- `factorize_file_dna_w_rc()` - File-based DNA factorization with RC
- `count_factors_dna_w_rc()` - DNA factor counting with RC
- `count_factors_file_dna_w_rc()` - File-based DNA factor counting with RC
- `write_factors_binary_file_dna_w_rc()` - Binary output with RC

**Usage:**
```bash
# Default settings
python benchmarks/dna_benchmark.py

# Custom nucleotide sizes
python benchmarks/dna_benchmark.py --custom-sizes 1000 10000 100000

# Custom output directory
python benchmarks/dna_benchmark.py --output-dir my_dna_results
```

**Output:**
- DNA-specific performance metrics
- Reverse complement overhead analysis
- Memory and time scaling
- Comprehensive plots

### 3. FASTA File Benchmark (`fasta_benchmark.py`)

Benchmarks FASTA file processing functions with both reverse complement and no-RC modes.

**Functions benchmarked:**
- `factorize_fasta_multiple_dna_w_rc()` - FASTA with reverse complement
- `factorize_fasta_multiple_dna_no_rc()` - FASTA without reverse complement
- `write_factors_binary_file_fasta_multiple_dna_w_rc()` - Binary output with RC
- `write_factors_binary_file_fasta_multiple_dna_no_rc()` - Binary output without RC

**Usage:**
```bash
# Default: 1kbp to 1Mbp
python benchmarks/fasta_benchmark.py

# Custom range
python benchmarks/fasta_benchmark.py --min-size 1000 --max-size 10000000 --num-sizes 15

# Specific sizes
python benchmarks/fasta_benchmark.py --custom-sizes 10000 100000 1000000
```

**Output:**
- Multi-sequence FASTA performance
- Comparison of RC vs no-RC modes
- Disk space analysis
- Resource prediction parameters

### 4. Parallel Factorization Benchmark (`parallel_benchmark.py`)

Benchmarks parallel factorization implementations with different thread counts.

**Functions benchmarked:**
- `parallel_factorize_to_file()` - Parallel text factorization
- `parallel_factorize_file_to_file()` - Parallel file-based factorization
- `parallel_factorize_dna_w_rc_to_file()` - Parallel DNA with RC
- `parallel_factorize_file_dna_w_rc_to_file()` - Parallel DNA file with RC

**Usage:**
```bash
# Default settings (auto-detect CPU count)
python benchmarks/parallel_benchmark.py

# Custom thread counts
python benchmarks/parallel_benchmark.py --custom-threads 1 2 4 8

# Specific sizes and threads
python benchmarks/parallel_benchmark.py --custom-sizes 10000 100000 --custom-threads 1 2 4
```

**Output:**
- Speedup analysis (vs single-threaded)
- Parallel efficiency metrics
- Thread scaling plots
- Performance vs input size

### 5. Unified Benchmark Runner (`run_all_benchmarks.py`)

Convenient script to run multiple benchmark suites with consistent parameters.

**Usage:**
```bash
# Run all benchmarks
python benchmarks/run_all_benchmarks.py --all

# Quick test (1k, 10k, 100k)
python benchmarks/run_all_benchmarks.py --all --quick

# Specific benchmarks with custom configuration
python benchmarks/run_all_benchmarks.py --core --parallel --custom-sizes 1000 50000 100000

# Large comprehensive test
python benchmarks/run_all_benchmarks.py --all --large --runs 5
```

**Size presets:**
- `--quick`: 1k, 10k, 100k (fast, for testing)
- `--small`: 1k-100k, 7 sizes
- `--medium`: 1k-1M, 10 sizes (default)
- `--large`: 1k-10M, 13 sizes

## Files

### Individual Benchmark Scripts

#### `core_benchmark.py`
Core factorization functions for general text.

**Output files:**
- `core_results/benchmark_results.json` - Raw data
- `core_results/trend_parameters.json/pkl` - Trend coefficients
- `core_results/core_benchmark_plots.png/pdf` - Visualizations

#### `dna_benchmark.py`
DNA-specific factorization with reverse complement.

**Output files:**
- `dna_results/benchmark_results.json` - Raw data
- `dna_results/trend_parameters.json/pkl` - Trend coefficients
- `dna_results/dna_benchmark_plots.png/pdf` - Visualizations

#### `fasta_benchmark.py`
FASTA file processing benchmark.

**Output files:**
- `fasta_results/benchmark_results.json` - Raw data
- `fasta_results/trend_parameters.json/pkl` - Trend coefficients
- `fasta_results/fasta_benchmark_plots.png/pdf` - Visualizations

#### `parallel_benchmark.py`
Parallel factorization with thread scaling.

**Output files:**
- `parallel_results/benchmark_results.json/pkl` - Raw data
- `parallel_results/parallel_benchmark_plots.png/pdf` - Speedup/efficiency plots

#### `run_all_benchmarks.py`
Unified runner for all benchmark suites.

**Output files:**
- `all_results/benchmark_summary.json` - Overall summary
- Individual results in subdirectories

#### `fasta_predictor.py`
Resource prediction utility for cluster job planning.

**Usage:**
```bash
# Predict resources for a specific size
python benchmarks/fasta_predictor.py benchmarks/fasta_results/trend_parameters.pkl --size 500000

# Generate resource table
python benchmarks/fasta_predictor.py benchmarks/fasta_results/trend_parameters.pkl \
    --function factorize_fasta_multiple_dna_w_rc \
    --table 10000 100000 1000000 10000000
```

## Benchmark Results

The benchmark generates log-log plots showing:
1. **Execution Time vs Input Size** - Nearly linear relationships (slope ≈ 1.0)
2. **Memory Usage vs Input Size** - Linear to slightly super-linear scaling
3. **Disk Space vs Input Size** - For binary functions, shows compression ratio
4. **Throughput vs Input Size** - MB/s processing rate
5. **Compression Ratio vs Input Size** - Factors per nucleotide
6. **Processing Efficiency vs Input Size** - Time per factor

## Trend Analysis

All functions show excellent power-law scaling relationships with R² > 0.99:

- **Time scaling**: ~O(n) where n is input size
- **Memory scaling**: ~O(n) for factorization functions
- **Disk space scaling**: ~O(n^0.88) for binary output (slight compression)

## Usage for Cluster Resource Estimation

1. Run benchmark once on your target machine:
   ```bash
   python benchmarks/fasta_benchmark.py --output-dir my_results
   ```

2. Use predictor for any future job sizes:
   ```bash
   python benchmarks/fasta_predictor.py my_results/trend_parameters.pkl --size 10000000
   ```

The predictor includes safety factors and rounds memory allocations to common cluster values (powers of 2).

## Example Predictions

For a 1 Mbp input (1,000,000 nucleotides):
- **factorize_fasta_multiple_dna_w_rc**: ~27 seconds, ~14 MB memory
- **write_factors_binary_file_fasta_multiple_dna_w_rc**: ~27 seconds, ~2.4 MB disk space
- **Cluster allocation**: 1 GB memory (with safety factor), 30-45 seconds time limit

## Dependencies

- numpy
- scipy  
- matplotlib
- noLZSS (with C++ extension built)

## Performance Characteristics

Based on the benchmark results, all FASTA functions show excellent power-law scaling:

- **Time complexity**: O(n) where n = input size
- **Memory complexity**: O(n) for factorization, O(1) for binary output
- **Disk space**: O(n^0.88) for binary files (slight compression)
- **R² values**: > 0.99 for all trend fits

This makes resource prediction highly accurate for any input size from 1 kbp to 10+ Mbp.