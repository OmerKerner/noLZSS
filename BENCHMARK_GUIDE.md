# noLZSS Benchmark Guide

This guide provides comprehensive information about the benchmarking suite for noLZSS factorization functions.

## Overview

The noLZSS benchmarking suite now includes comprehensive benchmarks for all factorization functions:

1. **Core Factorization** - General text factorization functions
2. **DNA Factorization** - DNA-specific functions with reverse complement support
3. **FASTA Processing** - Multi-sequence FASTA file processing
4. **Parallel Factorization** - Thread-parallel implementations with speedup analysis

## Quick Reference

### Run All Benchmarks

```bash
# Quick test (1k, 10k, 100k - takes ~2 minutes)
python benchmarks/run_all_benchmarks.py --all --quick

# Default test (1k to 1M, 10 sizes - takes ~15 minutes)
python benchmarks/run_all_benchmarks.py --all

# Comprehensive test (1k to 10M, 13 sizes - takes ~45 minutes)
python benchmarks/run_all_benchmarks.py --all --large --runs 5
```

### Run Individual Benchmarks

```bash
# Core factorization
python benchmarks/core_benchmark.py --custom-sizes 1000 10000 100000

# DNA with reverse complement
python benchmarks/dna_benchmark.py --custom-sizes 1000 10000 100000

# FASTA files
python benchmarks/fasta_benchmark.py --custom-sizes 1000 10000 100000

# Parallel factorization
python benchmarks/parallel_benchmark.py --custom-sizes 10000 --custom-threads 1 2 4
```

## Benchmark Suites

### 1. Core Factorization (`core_benchmark.py`)

Measures performance of general-purpose factorization functions.

**Functions:**
- `factorize()` - In-memory string factorization
- `factorize_file()` - File-based factorization
- `count_factors()` - Fast factor counting
- `count_factors_file()` - File-based factor counting
- `write_factors_binary_file()` - Direct binary output

**Key Metrics:**
- Execution time vs input size
- Memory usage vs input size
- Throughput (MB/s)
- Binary file compression ratio

**Typical Results:**
- Time complexity: O(n) where n = input size
- Memory complexity: O(n) for factorization, O(1) for binary output
- Throughput: 5-15 MB/s (hardware dependent)

### 2. DNA Factorization (`dna_benchmark.py`)

Specialized benchmarks for DNA sequences with reverse complement matching.

**Functions:**
- `factorize_dna_w_rc()` - DNA factorization with RC
- `factorize_file_dna_w_rc()` - File-based DNA with RC
- `count_factors_dna_w_rc()` - Fast counting with RC
- `count_factors_file_dna_w_rc()` - File-based counting with RC
- `write_factors_binary_file_dna_w_rc()` - Binary output with RC

**Key Metrics:**
- Performance with reverse complement enabled
- Memory overhead for RC support
- Compression efficiency on DNA sequences

**Typical Results:**
- Time complexity: O(n) similar to core functions
- RC overhead: ~10-20% additional time
- Better compression on DNA vs random text

### 3. FASTA Processing (`fasta_benchmark.py`)

Benchmarks for multi-sequence FASTA file processing.

**Functions:**
- `factorize_fasta_multiple_dna_w_rc()` - FASTA with RC
- `factorize_fasta_multiple_dna_no_rc()` - FASTA without RC
- `write_factors_binary_file_fasta_multiple_dna_w_rc()` - Binary with RC
- `write_factors_binary_file_fasta_multiple_dna_no_rc()` - Binary without RC

**Key Metrics:**
- Multi-sequence processing efficiency
- RC vs no-RC performance comparison
- Disk space requirements

**Typical Results:**
- Time scaling: O(n) for total sequence length
- R² > 0.99 for trend line fits
- Disk space: ~O(n^0.88) (compression)

### 4. Parallel Factorization (`parallel_benchmark.py`)

Measures parallel speedup and efficiency with multiple threads.

**Functions:**
- `parallel_factorize_to_file()` - Parallel text factorization
- `parallel_factorize_file_to_file()` - Parallel file processing
- `parallel_factorize_dna_w_rc_to_file()` - Parallel DNA with RC
- `parallel_factorize_file_dna_w_rc_to_file()` - Parallel DNA file with RC

**Key Metrics:**
- Speedup vs single-threaded
- Parallel efficiency percentage
- Scaling with thread count
- Best thread count for different input sizes

**Typical Results:**
- Speedup: 1.3-1.8x with 2 threads
- Efficiency: 65-90% at 2 threads
- Diminishing returns beyond 4-8 threads for smaller inputs

## Output Files

Each benchmark suite generates:

1. **`benchmark_results.json`** - Raw benchmark data with all measurements
2. **`trend_parameters.json`** - Fitted trend line parameters (JSON)
3. **`trend_parameters.pkl`** - Fitted trend line parameters (Python pickle)
4. **Plots** - Publication-quality PNG and PDF visualizations

### Plot Contents

**Core/DNA/FASTA Benchmarks:**
- Time vs Size (log-log scale)
- Memory vs Size (log-log scale)
- Throughput vs Size
- Disk Space vs Size (binary functions)
- Additional efficiency metrics

**Parallel Benchmark:**
- Speedup vs Thread Count
- Efficiency vs Thread Count
- Time vs Thread Count
- Time vs Size (max threads)

## Resource Prediction

Use trend parameters to predict resources for any input size:

```python
import pickle
import numpy as np

# Load trend parameters
with open('benchmarks/fasta_results/trend_parameters.pkl', 'rb') as f:
    trends = pickle.load(f)

# Predict time for 5 Mbp input
size = 5_000_000
time_trend = trends['factorize_fasta_multiple_dna_w_rc_time']
predicted_time = 10 ** (time_trend['slope'] * np.log10(size) + time_trend['intercept'])
print(f"Predicted time: {predicted_time:.2f} ms")
```

Or use the predictor script:

```bash
python benchmarks/fasta_predictor.py \
    benchmarks/fasta_results/trend_parameters.pkl \
    --size 5000000
```

## Size Presets

The unified runner provides convenient size presets:

- **`--quick`**: [1k, 10k, 100k] - Fast testing (~2 min)
- **`--small`**: [1k-100k, 7 sizes] - Basic profiling (~8 min)
- **`--medium`**: [1k-1M, 10 sizes] - Default, comprehensive (~15 min)
- **`--large`**: [1k-10M, 13 sizes] - Production profiling (~45 min)

## Best Practices

### For Development Testing

```bash
# Quick sanity check
python benchmarks/run_all_benchmarks.py --all --quick --runs 1

# Verify specific changes
python benchmarks/core_benchmark.py --custom-sizes 10000 --runs 3
```

### For Performance Profiling

```bash
# Get reliable measurements
python benchmarks/run_all_benchmarks.py --all --medium --runs 5

# Focus on specific area
python benchmarks/parallel_benchmark.py --custom-sizes 100000 1000000 --runs 5
```

### For Publication/Documentation

```bash
# Comprehensive, high-quality data
python benchmarks/run_all_benchmarks.py --all --large --runs 10 --output-dir paper_results

# Generate high-resolution plots
# (modify dpi=600 in the plotting functions)
```

## Interpreting Results

### Time Complexity

Look for the slope in log-log plots:
- Slope ≈ 1.0: Linear O(n) - Expected for most functions
- Slope > 1.0: Super-linear - May indicate algorithmic issue
- Slope < 1.0: Sub-linear - Indicates caching or overhead dominance

### Memory Scaling

Check memory trend:
- Linear with input size: Expected for in-memory factorization
- Constant: Expected for streaming/file-based functions
- Super-linear: Potential memory leak or inefficiency

### Parallel Efficiency

For parallel benchmarks:
- 90-100% efficiency: Excellent scaling
- 70-90% efficiency: Good scaling
- 50-70% efficiency: Acceptable for some workloads
- <50% efficiency: Consider reducing thread count

### R² Values

Trend fit quality:
- R² > 0.99: Excellent fit, reliable predictions
- R² > 0.95: Good fit, reasonable predictions
- R² < 0.95: Poor fit, investigate anomalies

## Troubleshooting

### Benchmark Fails

```bash
# Verify C++ extension is built
python -c "import noLZSS._noLZSS; print('OK')"

# Rebuild if needed
pip install -e .
```

### Memory Errors

```bash
# Use smaller sizes
python benchmarks/run_all_benchmarks.py --all --quick

# Or test specific sizes
python benchmarks/core_benchmark.py --custom-sizes 1000 10000
```

### Plots Not Generated

```bash
# Ensure matplotlib is installed
pip install matplotlib numpy scipy

# Check for plot files
ls benchmarks/*/results/*.png
```

## Dependencies

All benchmarks require:
- numpy
- scipy
- matplotlib
- noLZSS (with C++ extension built)

Install with:
```bash
pip install numpy scipy matplotlib
pip install -e .
```

## Contributing

When adding new benchmarks:

1. Follow existing patterns (time, memory, trend analysis)
2. Use consistent output structure (JSON + pickle + plots)
3. Add to unified runner (`run_all_benchmarks.py`)
4. Update this documentation
5. Test with `--quick` preset before committing

## Examples

### Example 1: Quick Development Check

```bash
# After making changes to core factorization
python benchmarks/core_benchmark.py --custom-sizes 10000 --runs 3
```

### Example 2: Parallel Performance Analysis

```bash
# Test parallel scaling on your machine
python benchmarks/parallel_benchmark.py \
    --custom-sizes 100000 \
    --custom-threads 1 2 4 8 16 \
    --runs 5
```

### Example 3: Production Profiling

```bash
# Comprehensive benchmark for cluster resource planning
python benchmarks/run_all_benchmarks.py \
    --all \
    --large \
    --runs 10 \
    --output-dir production_benchmarks
```

### Example 4: Compare FASTA Modes

```bash
# Benchmark both RC and no-RC modes
python benchmarks/fasta_benchmark.py \
    --custom-sizes 10000 100000 1000000 \
    --runs 5
    
# Check results
python benchmarks/fasta_predictor.py \
    benchmarks/fasta_results/trend_parameters.pkl \
    --size 5000000
```

## Conclusion

The noLZSS benchmarking suite provides comprehensive performance analysis for all factorization functions. Use it to:

- Validate performance characteristics
- Plan computational resources
- Compare implementations
- Identify optimization opportunities
- Generate publication-quality plots

For questions or issues, refer to the main README or open an issue on GitHub.
