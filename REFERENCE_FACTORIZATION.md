# Reference Sequence Factorization

This document describes the new reference sequence factorization functionality added to noLZSS.

## Overview

The new `factorize_w_reference_seq` functions allow factorization of a target DNA sequence using a reference sequence. This is particularly useful for genomics applications where you want to compress or analyze a target sequence relative to a known reference (like a reference genome).

## Key Features

- **Reference-based factorization**: Concatenates reference and target sequences, then factorizes only the target part
- **Reverse complement awareness**: Finds matches in both forward and reverse complement orientations
- **Efficient starting position**: Factorization begins from where the target sequence starts, not from the beginning
- **Adjusted coordinates**: Factor positions are relative to the target sequence, not the combined reference+target string

## Implementation Details

### Core Algorithm Changes

1. **Modified `nolzss_multiple_dna_w_rc()`**: Added `start_pos` parameter to begin factorization from a specific position rather than always starting at position 0.

2. **New Functions**:
   - `factorize_w_reference_seq()`: In-memory factorization with reference
   - `factorize_w_reference_seq_file()`: File-based version that writes binary output

### Function Signatures

#### C++ API
```cpp
std::vector<Factor> factorize_w_reference_seq(
    const std::string& reference_seq, 
    const std::string& target_seq
);

size_t factorize_w_reference_seq_file(
    const std::string& reference_seq, 
    const std::string& target_seq, 
    const std::string& out_path
);
```

#### Python API
```python
def factorize_w_reference_seq(
    reference_seq: Union[str, bytes], 
    target_seq: Union[str, bytes], 
    validate: bool = True
) -> List[Tuple[int, int, int, bool]]

def factorize_w_reference_seq_file(
    reference_seq: Union[str, bytes], 
    target_seq: Union[str, bytes], 
    output_path: Union[str, Path], 
    validate: bool = True
) -> int
```

## Usage Examples

### Basic Usage
```python
import noLZSS

# Define sequences
reference = "ATCGATCGATCGATCGATCG"
target = "GATCGATCGATCGA"

# Factorize target using reference
factors = noLZSS.factorize_w_reference_seq(reference, target)

# Each factor is (start, length, ref_pos, is_reverse_complement)
for start, length, ref_pos, is_rc in factors:
    print(f"Position {start}-{start+length-1}: "
          f"{'RC' if is_rc else 'FWD'} match at ref pos {ref_pos}")
```

### File Output
```python
# Write factors to binary file
num_factors = noLZSS.factorize_w_reference_seq_file(
    reference, target, "output.bin"
)
print(f"Wrote {num_factors} factors to file")
```

## Technical Details

### String Format
The algorithm concatenates sequences as: `REFERENCE[sentinel]TARGET[sentinel]RC(TARGET)[sentinel]RC(REFERENCE)[sentinel]`

### Starting Position Calculation
For reference of length R and target of length T:
- Combined string length: ~2(R + T + 2) (including sentinels and reverse complements)
- Target starts at position: R + 1 (after reference and one sentinel)
- Factorization begins at this position

### Factor Coordinate Adjustment
All factor start positions are adjusted to be relative to the target sequence:
```cpp
adjusted_factor.start = original_factor.start - target_start_pos;
```

## Benefits

1. **Genomics Applications**: Perfect for read mapping, variant calling, and reference-based compression
2. **Efficiency**: Only factorizes the target sequence, not the entire reference
3. **Reverse Complement Support**: Handles both DNA strands automatically
4. **Memory Efficient**: File output version for large-scale processing

## Performance Considerations

- Reference sequence is included in suffix tree construction but not factorized
- Target sequence can reference patterns in both reference and its reverse complement
- Memory usage scales with combined reference + target size
- File output version recommended for large datasets

## Testing

Comprehensive tests are available in `tests/test_reference_seq.py` covering:
- Basic functionality
- File output
- Edge cases
- Reverse complement matching
- Input validation

Run tests with:
```bash
python tests/test_reference_seq.py
```

## Future Enhancements

Potential improvements could include:
- Multiple target sequences against single reference
- Streaming API for very large references
- Optional quality score integration for genomics applications
- Advanced filtering options for minimum match lengths