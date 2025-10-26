# Shuffle Analysis Feature

## Overview

The shuffle analysis feature in `batch_factorize.py` allows you to create shuffled versions of FASTA sequences as control datasets and compare their factorization patterns against the original sequences. This is useful for:

- **Statistical validation**: Demonstrating that observed compressibility is due to sequence structure, not random composition
- **Control experiments**: Establishing baseline compression for randomly ordered sequences
- **Biological insights**: Identifying whether sequence patterns are meaningful or random

## How It Works

### 1. Sequence Shuffling

When `--shuffle-analysis` is enabled, the tool:
- Reads each sequence from your FASTA file
- Randomizes the nucleotides within each sequence independently
- Preserves the original sequence headers
- Maintains the nucleotide composition (same A, C, T, G counts)
- Writes shuffled sequences to temporary files

### 2. Factorization

Both the original and shuffled sequences are factorized using the noLZSS algorithm:
- Original sequences → `output/with_reverse_complement/` or `output/without_reverse_complement/`
- Shuffled sequences → `output/shuffled/with_reverse_complement/` or `output/shuffled/without_reverse_complement/`

### 3. Comparison Plots

Automatically generated plots show:
- **Solid line (blue)**: Original sequence factorization (cumulative factor length vs. factor count)
- **Dotted line (black)**: Shuffled control factorization
- Saved to: `output/comparison_plots/`

## Usage

### Basic Usage

```bash
python -m noLZSS.genomics.batch_factorize \
    my_sequences.fasta \
    --output-dir results \
    --mode with_reverse_complement \
    --shuffle-analysis
```

### With Reproducible Shuffling

Use `--shuffle-seed` for reproducible results:

```bash
python -m noLZSS.genomics.batch_factorize \
    my_sequences.fasta \
    --output-dir results \
    --mode with_reverse_complement \
    --shuffle-analysis \
    --shuffle-seed 42
```

### Batch Processing

Process multiple files with shuffle analysis:

```bash
python -m noLZSS.genomics.batch_factorize \
    --file-list my_files.txt \
    --output-dir results \
    --mode both \
    --shuffle-analysis \
    --shuffle-seed 123
```

## Output Structure

```
results/
├── with_reverse_complement/
│   ├── file1.bin              # Original factorization
│   └── file2.bin
├── shuffled/
│   └── with_reverse_complement/
│       ├── file1_shuffled.bin  # Shuffled factorization
│       └── file2_shuffled.bin
└── comparison_plots/
    ├── file1_with_reverse_complement_comparison.png
    └── file2_with_reverse_complement_comparison.png
```

## Interpreting Results

### Factor Count Comparison

- **Fewer factors (original)**: Indicates structure/repetition in sequence
- **More factors (shuffled)**: Random sequences are less compressible
- **Similar factor counts**: Sequence may lack significant structure

### Example Interpretation

```
Original:  8 factors for 64 bp = 0.125 factors/bp (high compression)
Shuffled: 29 factors for 64 bp = 0.453 factors/bp (low compression)
```

This shows the original sequence has ~3.6× better compression, indicating strong internal structure (likely repetitive patterns).

## Command-Line Options

| Option | Description |
|--------|-------------|
| `--shuffle-analysis` | Enable shuffle analysis (creates shuffled versions and comparison plots) |
| `--shuffle-seed INT` | Random seed for reproducibility (optional) |
| `--mode MODE` | Factorization mode: `with_reverse_complement`, `without_reverse_complement`, or `both` |
| `--force` | Overwrite existing output files (by default, skips existing) |
| `--max-workers N` | Number of parallel workers for processing |

## Performance Considerations

- **Disk Space**: Shuffled files and additional factorizations require extra storage
- **Processing Time**: Roughly doubles processing time (original + shuffled)
- **Memory**: Shuffling happens per-file, so memory usage scales with individual file size

## Limitations

1. **Sequence Independence**: Each sequence in a multi-sequence FASTA is shuffled independently
2. **Composition Preservation**: Only nucleotide order changes, not composition
3. **Header Preservation**: Original headers are maintained in shuffled files
4. **Skip Existing**: By default, existing shuffled files/factors are reused (use `--force` to regenerate)

## Examples

### Example 1: Single FASTA File

```bash
# Create shuffled control and comparison plot
python -m noLZSS.genomics.batch_factorize \
    genome.fasta \
    --output-dir genome_analysis \
    --mode with_reverse_complement \
    --shuffle-analysis \
    --shuffle-seed 12345
```

### Example 2: Multiple Files with Both Modes

```bash
# Process multiple chromosomes with both factorization modes
python -m noLZSS.genomics.batch_factorize \
    chr1.fasta chr2.fasta chr3.fasta \
    --output-dir chromosome_analysis \
    --mode both \
    --shuffle-analysis \
    --shuffle-seed 99
```

### Example 3: Batch Processing from File List

```bash
# Create file list
cat > genomes.txt <<EOF
genome1.fasta
genome2.fasta
genome3.fasta
EOF

# Process all genomes with shuffle analysis
python -m noLZSS.genomics.batch_factorize \
    --file-list genomes.txt \
    --output-dir batch_results \
    --mode with_reverse_complement \
    --shuffle-analysis \
    --max-workers 4
```

## Troubleshooting

### No Plots Generated

**Issue**: Factorization succeeds but no plots are created.

**Solution**: Install matplotlib:
```bash
pip install matplotlib
```

### Shuffled Files Not Created

**Issue**: Using `--force` but shuffled files still not regenerated.

**Solution**: The `--force` flag applies to binary factorization outputs, not FASTA files. Delete the shuffled FASTA files manually if you want to regenerate them with a different seed.

### Memory Issues with Large Files

**Issue**: Out of memory errors with very large FASTA files.

**Solution**: Process files individually rather than in batch, or reduce `--max-workers`.

## Technical Details

### Shuffling Algorithm

- Uses Python's `random.shuffle()` with optional seed
- Operates on individual characters (nucleotides)
- Maintains original sequence length
- Each sequence shuffled independently

### Plot Format

- Format: PNG (300 DPI)
- Size: 10×6 inches
- Original: Blue solid line with circle markers
- Shuffled: Black dotted line with square markers
- X-axis: Cumulative factor length (bp)
- Y-axis: Factor index (count)

## References

For more information about noLZSS factorization, see:
- Main README: `../README.md`
- Batch factorization documentation: CLI help with `--help`
- Example scripts: `../examples/`
