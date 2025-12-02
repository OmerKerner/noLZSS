# Genomics Module

The genomics module provides specialized functions for biological sequence analysis using LZSS factorization.

## FASTA Processing

```{eval-rst}
.. automodule:: noLZSS.genomics.fasta
   :members:
   :undoc-members:
   :show-inheritance:
```

## Sequence Utilities

```{eval-rst}
.. automodule:: noLZSS.genomics.sequences
   :members:
   :undoc-members:
   :show-inheritance:
```

## Plotting and Visualization

```{eval-rst}
.. automodule:: noLZSS.genomics.plots
   :members:
   :undoc-members:
   :show-inheritance:
```

## Per-sequence Complexity Tables

`noLZSS.genomics.batch_factorize` now exposes a lightweight mode for computing the DNA LZSS complexity of each FASTA record with and without reverse complement awareness:

```bash
python -m noLZSS.genomics.batch_factorize my_sequences.fasta \
   --complexity-tsv results/complexity.tsv \
   --complexity-threads 8
```

The generated TSV contains three columns:

1. `sequence_id` – the exact FASTA header for the sequence
2. `complexity_w_rc` – factor count when reverse complements are allowed
3. `complexity_no_rc` – factor count without reverse complement matching

The command accepts local files or URLs (with optional `--download-dir`). No factor files are written when `--complexity-tsv` is supplied.

## Genomics Package

```{eval-rst}
.. automodule:: noLZSS.genomics
   :members:
   :undoc-members:
   :show-inheritance:
```