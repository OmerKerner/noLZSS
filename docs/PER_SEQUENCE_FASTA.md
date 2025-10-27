# Per-Sequence FASTA Factorization

## Overview

This document describes the new per-sequence FASTA factorization functionality added to noLZSS. Unlike the existing concatenation-based approach, these functions factorize each sequence in a FASTA file independently.

## Naming Schema

The new functions follow the pattern: `<operation>_fasta_dna_<rc_option>_per_sequence`

Examples:
- `factorize_fasta_dna_w_rc_per_sequence` - Factorize with reverse complement, per sequence
- `factorize_fasta_dna_no_rc_per_sequence` - Factorize without reverse complement, per sequence
- `write_factors_binary_file_fasta_dna_w_rc_per_sequence` - Write to binary file
- `count_factors_fasta_dna_w_rc_per_sequence` - Count factors only

## API Reference

### Python API

All functions are available via the `noLZSS._noLZSS` module:

```python
from noLZSS._noLZSS import (
    factorize_fasta_dna_w_rc_per_sequence,
    factorize_fasta_dna_no_rc_per_sequence,
    write_factors_binary_file_fasta_dna_w_rc_per_sequence,
    write_factors_binary_file_fasta_dna_no_rc_per_sequence,
    count_factors_fasta_dna_w_rc_per_sequence,
    count_factors_fasta_dna_no_rc_per_sequence,
    parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence,
    parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence,
)
```

### Function Signatures

#### Factorize and Return Results

```python
def factorize_fasta_dna_w_rc_per_sequence(fasta_path: str) -> tuple:
    """
    Factorize each sequence separately with reverse complement awareness.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        tuple: (per_sequence_factors, sequence_ids)
            - per_sequence_factors: List of factor lists, one per sequence
            - sequence_ids: List of sequence identifiers
            
    Each factor is a tuple: (start, length, ref, is_rc)
    """
```

#### Count Factors (Memory Efficient)

```python
def count_factors_fasta_dna_w_rc_per_sequence(fasta_path: str) -> int:
    """
    Count total factors across all sequences without storing them.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        int: Total number of factors across all sequences
    """
```

#### Write to Binary File

```python
def write_factors_binary_file_fasta_dna_w_rc_per_sequence(
    fasta_path: str, 
    out_path: str
) -> int:
    """
    Factorize and write to binary file.
    
    Args:
        fasta_path: Path to input FASTA file
        out_path: Path to output binary file
        
    Returns:
        int: Total number of factors written
    """
```

#### Parallel Processing

```python
def parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence(
    fasta_path: str, 
    out_path: str,
    num_threads: int = 0
) -> int:
    """
    Parallel version of write function.
    
    Args:
        fasta_path: Path to input FASTA file
        out_path: Path to output binary file
        num_threads: Number of threads (0 = auto-detect)
        
    Returns:
        int: Total number of factors written
    """
```

## Comparison with Existing Functions

### Existing Approach (`factorize_fasta_multiple_dna_*`)

```
Input: FASTA with sequences A, B, C

Step 1: Concatenate with sentinels
        "A!B@C$"

Step 2: Build single CST and factorize
        → Single factor list with sentinel markers
```

### New Approach (`factorize_fasta_dna_*_per_sequence`)

```
Input: FASTA with sequences A, B, C

Step 1: Process each independently
        Seq A → CST_A → Factors_A
        Seq B → CST_B → Factors_B  
        Seq C → CST_C → Factors_C

Step 2: Return per-sequence results
        → [[Factors_A], [Factors_B], [Factors_C]]
```

### Feature Comparison

| Feature | Existing | New (per-sequence) |
|---------|----------|-------------------|
| **Method** | Concatenate with sentinels | Independent processing |
| **CST building** | Single CST | One CST per sequence |
| **Cross-sequence matches** | Yes (via concatenation) | No |
| **Sequence limit** | 125 (w/ RC) or 250 (no RC) | Unlimited |
| **Result format** | Single list + sentinel indices | Per-sequence lists |
| **Parallelization** | Text-based splitting | Sequence-based distribution |
| **Memory usage** | All in memory at once | Can process incrementally |

## Use Cases

### When to Use Per-Sequence Functions

1. **Large numbers of sequences**: No sentinel limitations
2. **Independent analysis**: Each sequence analyzed separately
3. **Downstream per-sequence processing**: Results already organized by sequence
4. **Comparative analysis**: Easy to compare compression across sequences

### When to Use Concatenated Functions

1. **Cross-sequence pattern detection**: Need to find patterns spanning sequences
2. **Small number of sequences**: Under sentinel limits (125/250)
3. **Single compression ratio**: Want one metric for entire file

## Example Usage

```python
from noLZSS._noLZSS import factorize_fasta_dna_w_rc_per_sequence

# Factorize each sequence separately
per_seq_factors, sequence_ids = factorize_fasta_dna_w_rc_per_sequence("input.fasta")

# Process each sequence independently
for seq_id, factors in zip(sequence_ids, per_seq_factors):
    print(f"Sequence {seq_id}:")
    print(f"  Number of factors: {len(factors)}")
    
    # Analyze factors for this sequence
    total_length = sum(f[1] for f in factors)  # f[1] is length
    rc_count = sum(1 for f in factors if f[3])  # f[3] is is_rc
    
    print(f"  Total coverage: {total_length} bp")
    print(f"  RC matches: {rc_count}/{len(factors)}")
```

## Binary File Format

The binary output format includes per-sequence metadata:

```
[Factors for all sequences: 24 bytes each]
[Sequence IDs: null-terminated strings]
[Factors per sequence: uint64_t array]
[Footer: FactorFileFooter struct]
```

The footer includes:
- Total number of factors
- Number of sequences
- Number of sentinels (always 0 for per-sequence)
- Footer size
- Total length

## Testing

Comprehensive test suite in `tests/test_per_sequence_fasta.py`:

- Basic factorization (with/without RC)
- Count functions
- Binary file I/O
- Parallel processing
- Sequential vs parallel consistency
- Per-sequence independence verification
- Error handling

All tests passing (12/12).

## Performance Considerations

### Parallelization

- **Sequential** (`num_threads=1`): Processes sequences one by one
- **Parallel** (`num_threads>1`): Distributes sequences across threads
- **Auto-detect** (`num_threads=0`): Uses min(num_sequences, hardware_concurrency)

### Memory Usage

- Each sequence requires building a CST
- For very large sequences, memory usage can be significant
- Use count functions if only factor counts are needed

### When to Use Parallel

- Multiple sequences in FASTA file
- Sequences are reasonably sized (>10KB each)
- Available CPU cores
- I/O not bottleneck

## Implementation Notes

### Key Design Decisions

1. **Use existing preparation functions**: Reuses `prepare_multiple_dna_sequences_*` with single-element vectors for validation
2. **Use `factorize_multiple_dna_w_rc`**: Expects already-prepared strings (avoids double-preparation)
3. **Per-sequence metadata**: Binary format tracks factor counts per sequence for easy reconstruction

### C++ Implementation Details

- Located in `src/cpp/fasta_processor.{hpp,cpp}` and `src/cpp/parallel_fasta_processor.{hpp,cpp}`
- New struct: `FastaPerSequenceFactorizationResult`
- Parallel implementation uses work-stealing pattern with atomic counter
- Binary format compatible with existing footer structure
