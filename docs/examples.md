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

## Performance Optimization

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

## Benchmarking and Analysis

### Plotting and Visualization

### Performance Comparison

### Advanced Genomics Example

This examples documentation provides comprehensive, working code samples that demonstrate all major features of the noLZSS library, from basic usage to advanced genomics applications and performance optimization techniques.
