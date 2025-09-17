# Examples and Usage

This page provides comprehensive examples demonstrating the noLZSS library's capabilities, from basic string factorization to advanced genomics applications.

## Basic Usage

### String Factorization

```python
import noLZSS

# Simple string factorization
text = "abcabcabc"
factors = noLZSS.factorize(text)
print(factors)
# Output: [(0, 1, 0), (1, 1, 1), (2, 1, 2), (3, 3, 0), (6, 3, 0)]

# Factorize bytes (recommended for non-ASCII text)
data = b"hello world hello"
factors = noLZSS.factorize(data)
print(f"Found {len(factors)} factors")

# Understanding factor format: (position, length, reference)
text = "abracadabra"
factors = noLZSS.factorize(text)
for i, (pos, length, ref) in enumerate(factors):
    if ref == 0:
        print(f"Factor {i}: New character '{text[pos]}' at position {pos}")
    else:
        substring = text[pos:pos+length]
        ref_substring = text[ref:ref+length]
        print(f"Factor {i}: '{substring}' at pos {pos}, references pos {ref} ('{ref_substring}')")
```

### Enhanced Factorization with Metadata

```python
# Get detailed analysis with factorization
result = noLZSS.factorize_with_info("the quick brown fox jumps over the lazy dog")
factors = result['factors']

print(f"Input text: '{result['input_text']}'")
print(f"Number of factors: {result['num_factors']}")
print(f"Input size: {result['input_size']} characters")
print(f"Compression ratio: {result['num_factors'] / result['input_size']:.3f}")

# Alphabet analysis
alphabet_info = result['alphabet_info']
print(f"\nAlphabet analysis:")
print(f"  Size: {alphabet_info['size']} unique characters")
print(f"  Characters: {alphabet_info['characters']}")
print(f"  Entropy: {alphabet_info['entropy']:.3f} bits")
print(f"  Most common: {alphabet_info['most_common'][:3]}")  # Top 3
```

### File Processing

```python
# Create a sample file for demonstration
sample_text = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. " * 100
with open("sample.txt", "w") as f:
    f.write(sample_text)

# Process files efficiently
factors = noLZSS.factorize_file("sample.txt")
print(f"File factorization: {len(factors)} factors")

# Count factors without storing them (memory efficient for large files)
count = noLZSS.count_factors_file("sample.txt")
print(f"Factor count: {count}")

# Performance optimization with reserve hint
# If you know approximately how many factors to expect:
factors = noLZSS.factorize_file("sample.txt", reserve_hint=1000)
print(f"Optimized factorization: {len(factors)} factors")

# Clean up
import os
os.remove("sample.txt")
```

### Input Validation and Error Handling

```python
# Input validation examples
try:
    # Empty input
    factors = noLZSS.factorize("")
except ValueError as e:
    print(f"Empty input error: {e}")

try:
    # Invalid file path
    factors = noLZSS.factorize_file("nonexistent.txt")
except FileNotFoundError as e:
    print(f"File not found: {e}")

# Disable validation for performance (use with caution)
text = "valid input"
factors = noLZSS.factorize(text, validate=False)
print(f"Fast factorization: {len(factors)} factors")
```

## Genomics Applications

### DNA Sequence Analysis

```python
import noLZSS.genomics
import os
from pathlib import Path

# Create sample DNA FASTA file
fasta_content = """>sequence1
ATCGATCGATCGATCG
>sequence2  
GCTAGCTAGCTAGCTA
>sequence3
AAATTTCCCGGG
"""

with open("sample_dna.fasta", "w") as f:
    f.write(fasta_content)

# Read and factorize nucleotide sequences
try:
    results = noLZSS.genomics.read_nucleotide_fasta("sample_dna.fasta")
    
    for seq_id, factors in results:
        print(f"\nSequence: {seq_id}")
        print(f"  Factors: {len(factors)}")
        print(f"  First few factors: {factors[:3]}")
        
except Exception as e:
    print(f"Error processing FASTA: {e}")

# Automatic sequence type detection
results = noLZSS.genomics.read_fasta_auto("sample_dna.fasta")
print(f"Auto-detected {len(results)} DNA sequences")

# Clean up
os.remove("sample_dna.fasta")
```

### Protein Sequence Analysis

```python
# Create sample protein FASTA file
protein_fasta = """>protein1
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
>protein2
ARNDCQEGHILKMFPSTWYV
"""

with open("sample_proteins.fasta", "w") as f:
    f.write(protein_fasta)

# Read protein sequences (returns sequences, not factors)
try:
    results = noLZSS.genomics.read_protein_fasta("sample_proteins.fasta")
    
    for seq_id, sequence in results:
        print(f"Protein: {seq_id}")
        print(f"  Length: {len(sequence)} amino acids")
        print(f"  First 20 AA: {sequence[:20]}")
        
        # Factorize the protein sequence
        factors = noLZSS.factorize(sequence)
        print(f"  Factors: {len(factors)}")
        
except Exception as e:
    print(f"Error processing protein FASTA: {e}")

# Clean up
os.remove("sample_proteins.fasta")
```

### High-Performance Genomics Processing

```python
# For large genomic datasets, use C++ implementations
import noLZSS.genomics

# Create larger sample file
large_fasta = ""
for i in range(10):
    large_fasta += f">sequence_{i}\n"
    large_fasta += "ATCGATCGATCG" * 100 + "\n"

with open("large_genome.fasta", "w") as f:
    f.write(large_fasta)

# Memory-efficient processing with C++ implementation
try:
    result = noLZSS.genomics.process_nucleotide_fasta("large_genome.fasta")
    
    print(f"Processed {result['num_sequences']} sequences")
    print(f"Total concatenated length: {len(result['sequence']):,} characters")
    print(f"Sequence IDs: {result['sequence_ids'][:3]}...")  # First 3
    print(f"Sequence lengths: {result['sequence_lengths'][:3]}...")
    
    # Factorize the concatenated sequences
    factors = noLZSS.factorize(result["sequence"])
    print(f"Total factors in concatenated sequences: {len(factors):,}")
    
except RuntimeError as e:
    print(f"C++ processing error: {e}")

# Clean up
os.remove("large_genome.fasta")
```

## Performance Optimization

### Memory-Efficient Processing

```python
# For very large files, use counting instead of storing factors
def analyze_large_file(filepath):
    """Analyze a large file without loading all factors into memory."""
    # Count factors first
    factor_count = noLZSS.count_factors_file(filepath)
    print(f"File contains {factor_count:,} factors")
    
    # Calculate approximate memory usage
    approx_memory_mb = (factor_count * 3 * 8) / (1024 * 1024)  # 3 ints per factor
    print(f"Approximate memory for factors: {approx_memory_mb:.1f} MB")
    
    if approx_memory_mb < 100:  # Only load if reasonable size
        factors = noLZSS.factorize_file(filepath, reserve_hint=factor_count)
        return factors
    else:
        print("File too large for in-memory processing")
        return None

# Example usage
with open("medium_file.txt", "w") as f:
    f.write("repeated pattern " * 10000)

factors = analyze_large_file("medium_file.txt")
if factors:
    print(f"Loaded {len(factors)} factors successfully")

os.remove("medium_file.txt")
```

### Batch Processing

```python
# Process multiple files efficiently
import glob
from pathlib import Path

# Create test files
test_files = ["file1.txt", "file2.txt", "file3.txt"]
for i, filename in enumerate(test_files):
    with open(filename, "w") as f:
        f.write(f"test content for file {i} " * (100 + i * 50))

def batch_process_files(file_pattern="*.txt"):
    """Process multiple files and return summary statistics."""
    results = []
    
    for filepath in glob.glob(file_pattern):
        try:
            # Use counting for initial analysis
            factor_count = noLZSS.count_factors_file(filepath)
            file_size = Path(filepath).stat().st_size
            
            results.append({
                'file': filepath,
                'size_bytes': file_size,
                'factor_count': factor_count,
                'compression_ratio': factor_count / file_size if file_size > 0 else 0
            })
            
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
    
    return results

# Process all text files
results = batch_process_files("file*.txt")

print("Batch processing results:")
for result in results:
    print(f"  {result['file']}: {result['factor_count']} factors, "
          f"ratio: {result['compression_ratio']:.4f}")

# Clean up
for filename in test_files:
    os.remove(filename)
```

## Advanced Features

### Binary Factor Storage

```python
# Create test data
test_text = "the quick brown fox jumps over the lazy dog" * 50
with open("input.txt", "w") as f:
    f.write(test_text)

# Write factors directly to binary file (memory efficient)
num_factors = noLZSS.write_factors_binary_file("input.txt", "factors.bin")
print(f"Wrote {num_factors} factors to binary file")

# Read factors back from binary file
factors = noLZSS.read_factors_binary_file("factors.bin")
print(f"Read {len(factors)} factors from binary file")

# Verify integrity
factors_direct = noLZSS.factorize_file("input.txt")
assert factors == factors_direct, "Binary storage integrity check failed"
print("Binary storage integrity verified!")

# Check file sizes
import os
text_size = os.path.getsize("input.txt")
binary_size = os.path.getsize("factors.bin")
print(f"Original text: {text_size} bytes")
print(f"Binary factors: {binary_size} bytes")
print(f"Storage ratio: {binary_size / text_size:.3f}")

# Clean up
os.remove("input.txt")
os.remove("factors.bin")
```

### Alphabet Analysis

```python
# Detailed alphabet analysis
texts = [
    "hello world",
    "ATCGATCGATCG",
    "The quick brown fox jumps over the lazy dog",
    "αβγδεζηθικλμνξοπρστυφχψω",  # Greek
    "🚀🌟💻🔬📊"  # Emojis
]

for text in texts:
    alphabet_info = noLZSS.analyze_alphabet(text)
    
    print(f"\nText: '{text[:30]}{'...' if len(text) > 30 else ''}'")
    print(f"  Alphabet size: {alphabet_info['size']}")
    print(f"  Entropy: {alphabet_info['entropy']:.3f} bits")
    print(f"  Most common: {alphabet_info['most_common'][:3]}")
    
    # Factorize and analyze compression
    factors = noLZSS.factorize(text)
    compression_ratio = len(factors) / len(text)
    print(f"  Factors: {len(factors)}, Compression ratio: {compression_ratio:.3f}")
```

## Benchmarking and Analysis

### Plotting and Visualization

```python
# Install plotting dependencies if needed:
# pip install noLZSS[plotting]

try:
    # Generate sample data for plotting
    sample_data = "abcdefghijk" * 1000  # Repetitive pattern
    factors = noLZSS.factorize(sample_data)
    
    # Plot factor length accumulation
    noLZSS.plot_factor_lengths(factors, save_path="factor_plot.png", show_plot=False)
    print("Factor length plot saved to factor_plot.png")
    
    # Plot from binary file
    with open("plot_input.txt", "w") as f:
        f.write(sample_data)
    
    noLZSS.write_factors_binary_file("plot_input.txt", "plot_factors.bin")
    noLZSS.plot_factor_lengths("plot_factors.bin", save_path="binary_plot.png", show_plot=False)
    print("Binary factor plot saved to binary_plot.png")
    
    # Clean up
    os.remove("plot_input.txt")
    os.remove("plot_factors.bin")
    if os.path.exists("factor_plot.png"):
        os.remove("factor_plot.png")
    if os.path.exists("binary_plot.png"):
        os.remove("binary_plot.png")
        
except ImportError:
    print("Plotting functionality requires: pip install noLZSS[plotting]")
except Exception as e:
    print(f"Plotting error: {e}")
```

### Performance Comparison

```python
import time

def benchmark_factorization_methods(text_size=10000):
    """Compare different factorization approaches."""
    # Generate test data with varying patterns
    test_text = "pattern" * (text_size // 7) + "unique_suffix"
    
    print(f"Benchmarking factorization (text size: {len(test_text):,} chars)")
    
    # Method 1: Direct string factorization
    start_time = time.time()
    factors1 = noLZSS.factorize(test_text)
    time1 = time.time() - start_time
    
    # Method 2: File-based factorization
    with open("benchmark.txt", "w") as f:
        f.write(test_text)
    
    start_time = time.time()
    factors2 = noLZSS.factorize_file("benchmark.txt")
    time2 = time.time() - start_time
    
    # Method 3: File-based with reserve hint
    start_time = time.time()
    factors3 = noLZSS.factorize_file("benchmark.txt", reserve_hint=len(factors1))
    time3 = time.time() - start_time
    
    # Method 4: Just counting (fastest)
    start_time = time.time()
    count = noLZSS.count_factors_file("benchmark.txt")
    time4 = time.time() - start_time
    
    print(f"Results:")
    print(f"  Direct string:     {time1:.4f}s, {len(factors1)} factors")
    print(f"  File-based:        {time2:.4f}s, {len(factors2)} factors")
    print(f"  With reserve hint: {time3:.4f}s, {len(factors3)} factors")
    print(f"  Count only:        {time4:.4f}s, {count} factors")
    
    # Verify all methods give same results
    assert len(factors1) == len(factors2) == len(factors3) == count
    assert factors1 == factors2 == factors3
    
    print(f"  All methods verified consistent!")
    
    os.remove("benchmark.txt")

# Run benchmark
benchmark_factorization_methods(50000)
```

### Advanced Genomics Example

```python
# Complete genomics workflow example
def genomics_analysis_workflow():
    """Demonstrate a complete genomics analysis workflow."""
    
    # Create realistic genomic data
    sequences = {
        "chr1_segment": "ATCGATCG" * 200,
        "chr2_segment": "GCTAGCTA" * 150,  
        "repetitive_region": "AAATTTCCCGGG" * 100,
        "variable_region": "ATCGATCGATCGATCGAAATTTCCCGGGATCGATCG" * 50
    }
    
    # Write multi-sequence FASTA
    fasta_content = ""
    for seq_id, sequence in sequences.items():
        fasta_content += f">{seq_id}\n{sequence}\n"
    
    with open("genome_analysis.fasta", "w") as f:
        f.write(fasta_content)
    
    print("Genomics Analysis Workflow")
    print("=" * 40)
    
    # Step 1: Process with Python wrapper
    print("1. Processing with Python wrapper...")
    results = noLZSS.genomics.read_nucleotide_fasta("genome_analysis.fasta")
    
    for seq_id, factors in results:
        original_len = len(sequences[seq_id])
        compression_ratio = len(factors) / original_len
        print(f"   {seq_id}: {original_len:,} bp → {len(factors)} factors "
              f"(ratio: {compression_ratio:.3f})")
    
    # Step 2: High-performance C++ processing
    print("\n2. High-performance C++ processing...")
    try:
        cpp_result = noLZSS.genomics.process_nucleotide_fasta("genome_analysis.fasta")
        total_factors = len(noLZSS.factorize(cpp_result["sequence"]))
        
        print(f"   Concatenated {cpp_result['num_sequences']} sequences")
        print(f"   Total length: {len(cpp_result['sequence']):,} characters")
        print(f"   Total factors: {total_factors:,}")
        
    except Exception as e:
        print(f"   C++ processing error: {e}")
    
    # Step 3: Analysis summary
    print("\n3. Analysis Summary:")
    total_bp = sum(len(seq) for seq in sequences.values())
    total_factors = sum(len(factors) for _, factors in results)
    overall_ratio = total_factors / total_bp
    
    print(f"   Total base pairs: {total_bp:,}")
    print(f"   Total factors: {total_factors:,}")
    print(f"   Overall compression ratio: {overall_ratio:.3f}")
    
    # Clean up
    os.remove("genome_analysis.fasta")

# Run the workflow
genomics_analysis_workflow()
```

This examples documentation provides comprehensive, working code samples that demonstrate all major features of the noLZSS library, from basic usage to advanced genomics applications and performance optimization techniques.
