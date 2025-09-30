# Command-Line Interface for Plotting

## Overview

The noLZSS package provides a convenient command-line interface for creating plots directly from the terminal. This is useful for quick visualization and batch processing.

## Available Plot Types

Use `python -m noLZSS.genomics.plots --help` to see all available plot types:

```
usage: plots.py [-h] {cumulative,self-factors-plot,reference-plot} ...

Run LZSS plots

positional arguments:
  {cumulative,self-factors-plot,reference-plot}
    cumulative          Plot cumulative factors
    self-factors-plot   Plot self-factors
    reference-plot      Plot target sequence factorized with reference sequence
```

## Reference Sequence Plotting

### Command Structure

```bash
python -m noLZSS.genomics.plots reference-plot <reference_seq> <target_seq> [options]
```

### Required Arguments

- `reference_seq`: The reference DNA sequence (A, C, T, G)
- `target_seq`: The target DNA sequence to be factorized

### Optional Arguments

- `--factors_filepath`: Path to binary factors file (will compute if not provided)
- `--reference_name`: Name for the reference sequence (default: "Reference")
- `--target_name`: Name for the target sequence (default: "Target")
- `--save_path`: Path to save the plot image
- `--show_plot`: Display the plot interactively
- `--interactive`: Use interactive Panel/Datashader plot instead of simple matplotlib
- `--return_panel`: Return Panel app object (interactive mode only)

### Examples

#### 1. Basic Usage - Simple Plot

```bash
python -m noLZSS.genomics.plots reference-plot \
    "ATCGATCGATCGATCG" \
    "ATCGCCCCGATCGAAA" \
    --save_path "my_plot.png"
```

This creates a simple matplotlib-based plot showing the factorization of the target sequence using the reference.

#### 2. With Custom Names

```bash
python -m noLZSS.genomics.plots reference-plot \
    "ATCGATCGATCGATCG" \
    "ATCGCCCCGATCGAAA" \
    --reference_name "Chromosome 1" \
    --target_name "Query Read" \
    --save_path "chr1_vs_query.png" \
    --show_plot
```

#### 3. Interactive Plot

```bash
python -m noLZSS.genomics.plots reference-plot \
    "ATCGATCGATCGATCGATCGATCGATCGATCG" \
    "ATCGCCCCGATCGAAAATTTTCCCCGGGG" \
    --interactive \
    --reference_name "Reference Genome" \
    --target_name "Query Sequence" \
    --save_path "interactive_plot.png"
```

The `--interactive` flag creates a Panel/Datashader plot with:
- Hover details for individual factors
- Dynamic length filtering controls
- Zoom/pan functionality
- High-performance rendering

#### 4. Using Pre-computed Factors

If you've already computed factors and saved them to a binary file:

```bash
python -m noLZSS.genomics.plots reference-plot \
    "ATCGATCGATCGATCG" \
    "ATCGCCCCGATCGAAA" \
    --factors_filepath "precomputed_factors.bin" \
    --save_path "plot_from_binary.png"
```

#### 5. Long Sequences from Files

For longer sequences stored in variables or files, you can use shell scripting:

```bash
REF_SEQ=$(cat reference.txt)
TARGET_SEQ=$(cat target.txt)

python -m noLZSS.genomics.plots reference-plot \
    "$REF_SEQ" \
    "$TARGET_SEQ" \
    --reference_name "Reference" \
    --target_name "Target" \
    --save_path "large_sequence_plot.png"
```

## Other Plot Types

### Cumulative Factors Plot

```bash
python -m noLZSS.genomics.plots cumulative \
    path/to/sequences.fasta \
    output_directory/ \
    --max_sequences 10 \
    --save_factors_binary
```

### Self-Factors Plot

```bash
python -m noLZSS.genomics.plots self-factors-plot \
    --fasta_filepath path/to/sequences.fasta \
    --name "My Analysis" \
    --save_path "self_factors.png" \
    --show_plot
```

Or using a binary factors file:

```bash
python -m noLZSS.genomics.plots self-factors-plot \
    --factors_filepath path/to/factors.bin \
    --name "My Analysis" \
    --save_path "self_factors.png"
```

## Tips and Best Practices

1. **For quick visualization**: Use the simple matplotlib plot (default)
2. **For large datasets**: Use `--interactive` for better performance and interactivity
3. **For batch processing**: Save plots without showing them (omit `--show_plot`)
4. **For publication**: Use simple plots with `--save_path` for clean, publication-ready figures
5. **For exploration**: Use interactive plots with hover details to examine specific factors

## Dependencies

- **Simple plots**: Requires `matplotlib` (install with `pip install 'noLZSS[plotting]'`)
- **Interactive plots**: Requires Panel stack (install with `pip install 'noLZSS[panel]'`)

## Troubleshooting

### Missing Dependencies

If you get an error about missing dependencies:

```bash
pip install 'noLZSS[plotting]'  # For matplotlib
pip install 'noLZSS[panel]'     # For interactive plots
```

### C++ Extension Not Built

If you get an error about the C++ extension:

```bash
pip install -e .  # Rebuild the package
```

## Getting Help

For detailed help on any command:

```bash
python -m noLZSS.genomics.plots --help
python -m noLZSS.genomics.plots reference-plot --help
python -m noLZSS.genomics.plots cumulative --help
python -m noLZSS.genomics.plots self-factors-plot --help
```