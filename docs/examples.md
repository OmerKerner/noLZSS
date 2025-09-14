# Examples and Usage

## Basic Usage

### String Factorization

```python
import noLZSS

# Simple string factorization
text = "abcabcabc"
factors = noLZSS.factorize(text)
print(factors)
# Output: [(0, 1, 0), (1, 1, 1), (2, 1, 2), (3, 3, 0), (6, 3, 0)]

# Get factorization with additional metadata
result = noLZSS.factorize_with_info(text)
factors = result['factors']
print(f"Number of factors: {result['num_factors']}")
print(f"Input size: {result['input_size']}")
print(f"Alphabet size: {result['alphabet_info']['size']}")
print(f"Alphabet characters: {result['alphabet_info']['characters']}")
print(f"Character distribution: {result['alphabet_info']['distribution']}")
print(f"Most common characters: {result['alphabet_info']['most_common']}")
print(f"Total length: {result['alphabet_info']['total_length']}")
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


## Performance Optimization

### Memory-Efficient Processing


### Batch Processing


## Advanced Features


### Binary Factor Storage


## Benchmarking and Analysis
