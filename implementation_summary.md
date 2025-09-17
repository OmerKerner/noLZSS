# Enhanced plot_multiple_seq_self_lz_factor_plot_from_fasta() Implementation

## 🎯 Problem Statement (COMPLETED)
1. ✅ Add ability to read factors from binary factor file format
2. ✅ Add horizontal and vertical lines at sentinel positions  
3. ✅ Show sequence names on axes to indicate different sequences
4. ✅ Separate sequences with visual boundaries

## 🆕 New Function Signature
```python
def plot_multiple_seq_self_lz_factor_plot_from_fasta(
    fasta_filepath: Union[str, Path] = None,        # Original input
    factors_filepath: Union[str, Path] = None,      # NEW: Binary input  
    name: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
    return_panel: bool = False
) -> Optional["panel.viewable.Viewable"]:
```

## 🔧 Key Implementation Changes

### 1. Binary Factor File Reading
- Added `read_factors_binary_file_with_metadata()` function
- Reads enhanced binary format with:
  - Factor data (start, length, ref, is_rc)  
  - Sentinel factor indices
  - Sequence names from FASTA headers
  - Metadata (num_sequences, num_sentinels)

### 2. Dual Input Support
```python
# Validate input arguments
if (fasta_filepath is None) == (factors_filepath is None):
    raise ValueError("Exactly one of fasta_filepath or factors_filepath must be provided")

# Handle both input types
if input_type == "fasta":
    factors, sentinel_factor_indices = factorize_fasta_multiple_dna_w_rc(str(input_filepath))
    # Parse FASTA headers for sequence names
else:
    metadata = read_factors_binary_file_with_metadata(input_filepath)
    factors = metadata['factors']
    sentinel_factor_indices = metadata['sentinel_factor_indices']  
    sequence_names = metadata['sequence_names']
```

### 3. Sentinel Boundary Calculation
```python
# Calculate sentinel positions for lines
sentinel_positions = []
for idx in sentinel_factor_indices:
    if idx < len(factors):
        sentinel_start = factors[idx][0]
        sentinel_positions.append(sentinel_start)

# Calculate sequence boundaries  
sequence_boundaries = []
prev_pos = 0
for i, pos in enumerate(sentinel_positions):
    seq_name = sequence_names[i] if i < len(sequence_names) else f"seq_{i}"
    sequence_boundaries.append((prev_pos, pos, seq_name))
    prev_pos = pos + 1
```

### 4. Visual Enhancements
```python
# Add sentinel lines  
for pos in sentinel_positions:
    if min_val <= pos <= max_val:
        # Red vertical line at sequence boundary
        v_line = hv.VLine(pos).opts(line_color='red', line_width=2, alpha=0.7)
        # Red horizontal line at sequence boundary
        h_line = hv.HLine(pos).opts(line_color='red', line_width=2, alpha=0.7)
        
# Add sequence name labels on axes
for start_pos, end_pos, seq_name in sequence_boundaries:
    mid_pos = (start_pos + end_pos) / 2
    # X-axis label (bottom)
    x_label = hv.Text(mid_pos, min_val - offset, seq_name)
    # Y-axis label (left side, rotated)  
    y_label = hv.Text(min_val - offset, mid_pos, seq_name, angle=90)
```

## ✅ Validation Results

### Test Coverage
- ✅ FASTA factorization works correctly
- ✅ Binary file creation with metadata 
- ✅ Binary metadata reading
- ✅ Sentinel boundary calculation (N sequences → N-1 sentinels → N boundaries)
- ✅ Function imports and signature validation
- ✅ Backwards compatibility maintained
- ✅ Core tests still pass

### Example Usage
```python
# Original usage (backwards compatible)
plot_multiple_seq_self_lz_factor_plot_from_fasta(
    fasta_filepath="sequences.fasta"
)

# New usage with binary factors
plot_multiple_seq_self_lz_factor_plot_from_fasta(
    factors_filepath="factors.bin"  # Pre-computed factors with metadata
)
```

## 🎉 Features Delivered

1. **Binary Factor Input**: Reads pre-computed factors with metadata for faster plotting
2. **Sentinel Visualization**: Red lines clearly separate sequence boundaries  
3. **Sequence Labeling**: Names from FASTA headers appear on both axes
4. **Enhanced Titles**: Show sequence count and names in plot titles
5. **Backwards Compatibility**: Existing code continues to work unchanged
6. **Robust Error Handling**: Validates inputs and provides clear error messages

The enhanced plotting function now provides rich visualization capabilities for multi-sequence factor plots with clear sequence boundaries and comprehensive labeling.
