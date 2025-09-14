# Examples and Usage

## Basic Usage

### String Factorization

```python
import noLZSS

# Simple string factorization
text = "abcabcabc"
factors = noLZSS.factorize(text)
print(factors)
# Output: [(0, 3, 'abc'), (0, 3, 'abc'), (0, 3, 'abc')]

# Get factorization with metadata
factors, info = noLZSS.factorize_with_info(text)
print(f"Number of factors: {info['num_factors']}")
print(f"Compression ratio: {info['compression_ratio']:.2f}")
```

### File Processing

```python
# Process large files efficiently
factors = noLZSS.factorize_file("large_input.txt")

# Count factors without storing them (memory efficient)
count = noLZSS.count_factors_file("large_input.txt")
print(f"File contains {count} factors")
```

## Genomics Applications

### DNA Sequence Analysis

```python
import noLZSS.genomics as genomics

# Process FASTA files
sequences = genomics.read_nucleotide_fasta("genome.fasta")
factors = genomics.factorize_sequences(sequences)

# Analyze factorization patterns
genomics.plot_factor_distribution(factors)
```

### Working with Biological Data

```python
# Validate nucleotide sequences
dna = "ATCGATCGATCG"
if genomics.is_valid_nucleotide_sequence(dna):
    factors = noLZSS.factorize(dna)
    
# Handle reverse complements
reverse_comp = genomics.reverse_complement(dna)
factors_rc = noLZSS.factorize(reverse_comp)
```

## Performance Optimization

### Memory-Efficient Processing

```python
# For very large files, use streaming approach
def process_large_genome(filepath):
    # Count first to get size estimates
    factor_count = noLZSS.count_factors_file(filepath)
    
    # Process with memory hint for better performance
    factors = noLZSS.factorize_file(filepath, reserve_hint=factor_count)
    return factors

# Use C++ implementations for best performance
factors = genomics.process_nucleotide_fasta("large_genome.fasta")
```

### Batch Processing

```python
import os
from pathlib import Path

def analyze_genome_directory(directory):
    """Analyze all FASTA files in a directory."""
    results = {}
    
    for fasta_file in Path(directory).glob("*.fasta"):
        print(f"Processing {fasta_file.name}...")
        
        # Use efficient file-based processing
        factor_count = noLZSS.count_factors_file(str(fasta_file))
        factors = noLZSS.factorize_file(str(fasta_file))
        
        results[fasta_file.name] = {
            'factor_count': factor_count,
            'factors': factors
        }
    
    return results
```

## Advanced Features

### Custom Input Validation

```python
# Disable validation for trusted input (performance gain)
factors = noLZSS.factorize(trusted_data, validate=False)

# Custom validation with error handling
try:
    data = noLZSS.validate_input(user_input)
    factors = noLZSS.factorize(data)
except noLZSS.InvalidInputError as e:
    print(f"Invalid input: {e}")
```

### Binary Factor Storage

```python
# Save factors in efficient binary format
factors = noLZSS.factorize_file("input.txt")
noLZSS.write_factors_binary_file(factors, "factors.bin")

# Load factors back
loaded_factors = noLZSS.read_factors_binary_file("factors.bin")
```

## Benchmarking and Analysis

```python
import time
from noLZSS.genomics import plot_factor_lengths

# Benchmark different approaches
def benchmark_factorization(data):
    start_time = time.time()
    factors = noLZSS.factorize(data)
    end_time = time.time()
    
    print(f"Factorization took {end_time - start_time:.2f} seconds")
    print(f"Found {len(factors)} factors")
    
    # Analyze factor length distribution
    plot_factor_lengths(factors, title="Factor Length Distribution")
    
    return factors

# Run benchmark
large_text = "your_large_text_here" * 1000
benchmark_factorization(large_text)
```