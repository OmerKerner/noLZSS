# Reference Sequence Factor Plotting

## Overview

This implementation adds new plotting functions to visualize DNA sequences that have been factorized with a reference sequence, compatible with the output of `factorize_w_reference_seq()`.

## New Functions

### 1. `plot_reference_seq_lz_factor_plot_simple()`

**Purpose**: Creates a simple matplotlib-based factor plot for reference sequence factorization.

**Key Features**:
- Simple matplotlib backend for basic plotting
- Distinct colors for reference vs target regions
- Clear sequence boundary visualization
- Legend and annotations
- PNG export capability

**Color Scheme**:
- **Red**: Forward factors in target region
- **Orange**: Reverse complement factors in target region  
- **Blue**: Forward factors in reference region (rare)
- **Dark Blue**: Reverse complement factors in reference region (rare)
- **Green**: Sequence boundary lines
- **Blue/Red backgrounds**: Reference/Target regions

### 2. `plot_reference_seq_lz_factor_plot()`

**Purpose**: Creates an advanced interactive Panel/Datashader plot with hover functionality and dynamic controls.

**Key Features**:
- Interactive Panel/Datashader backend
- Hover details for individual factors
- Dynamic length filtering
- Zoom/pan functionality
- High-performance rendering for large datasets
- Real-time plot updates

## Command-Line Interface

The plotting functions are available via command-line for easy access:

### Basic Usage

```bash
# Simple matplotlib plot
python -m noLZSS.genomics.plots reference-plot \
    "ATCGATCGATCGATCG" \
    "ATCGCCCCGATCGAAA" \
    --reference_name "Reference Genome" \
    --target_name "Query Sequence" \
    --save_path "reference_plot.png" \
    --show_plot
```

### Interactive Plot

```bash
# Interactive Panel/Datashader plot
python -m noLZSS.genomics.plots reference-plot \
    "ATCGATCGATCGATCG" \
    "ATCGCCCCGATCGAAA" \
    --interactive \
    --save_path "reference_plot.png"
```

### Using Pre-computed Factors

```bash
# Use binary factor file
python -m noLZSS.genomics.plots reference-plot \
    "ATCGATCGATCGATCG" \
    "ATCGCCCCGATCGAAA" \
    --factors_filepath "reference_factors.bin" \
    --save_path "reference_plot.png"
```

### Available Commands

```bash
# View all plot types
python -m noLZSS.genomics.plots --help

# View reference-plot specific help
python -m noLZSS.genomics.plots reference-plot --help
```

### Command-Line Arguments

**Required:**
- `reference_seq`: Reference DNA sequence
- `target_seq`: Target DNA sequence

**Optional:**
- `--factors_filepath`: Path to binary factors file (will compute if not provided)
- `--reference_name`: Name for reference sequence (default: "Reference")
- `--target_name`: Name for target sequence (default: "Target")
- `--save_path`: Path to save plot image
- `--show_plot`: Display the plot
- `--interactive`: Use interactive Panel plot instead of simple matplotlib
- `--return_panel`: Return Panel app object (interactive mode only)

## Python API Usage Examples

```python
from noLZSS.genomics.sequences import factorize_w_reference_seq
from noLZSS.genomics.plots import plot_reference_seq_lz_factor_plot_simple

# Define sequences
reference_seq = "ATCGATCGATCGATCG"
target_seq = "ATCGCCCCGATCGAAA"

# Create simple plot
plot_reference_seq_lz_factor_plot_simple(
    reference_seq=reference_seq,
    target_seq=target_seq,
    reference_name="Reference Genome",
    target_name="Query Sequence",
    save_path="reference_plot.png",
    show_plot=True
)
```

## Input Options

Both functions support multiple input methods:

1. **Direct factors**: Pass pre-computed factors from `factorize_w_reference_seq()`
2. **Auto-compute**: Functions will compute factors automatically if not provided
3. **Binary file**: Read factors from binary files created by `factorize_w_reference_seq_file()`

## Plot Interpretation

The plot shows:
- **X-axis**: Positions in concatenated sequence (reference + target)
- **Y-axis**: Reference positions that factors point to
- **Factors**: Lines from target positions to reference positions
- **Sequence regions**: Visually separated with colors and boundaries

This visualization helps understand:
- How much the target reuses patterns from the reference
- Which reference regions are most utilized
- The orientation of matches (forward vs reverse complement)
- The structure of the factorization

## Integration

The functions are:
- Available in `noLZSS.genomics.plots` module
- Exported via wildcard import in `__init__.py`
- Compatible with existing plotting infrastructure
- Follow established error handling patterns
- Include comprehensive documentation

## Dependencies

- **Simple version**: `matplotlib` (included in `[plotting]` extra)
- **Interactive version**: `panel`, `holoviews`, `datashader`, `bokeh`, etc. (included in `[panel]` extra)

## Files Added/Modified

1. **src/noLZSS/genomics/plots.py**: Added two new plotting functions
2. **examples/reference_sequence_plotting_example.py**: Comprehensive usage examples
3. Function exports automatically available via existing wildcard imports

## Compatibility

- Compatible with `factorize_w_reference_seq()` output format
- Follows existing noLZSS patterns for error handling and validation
- Supports both string and bytes input for sequences
- Works with binary factor files using established metadata format