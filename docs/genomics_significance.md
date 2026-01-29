# Factor Length Significance Analysis

## Introduction

Factor length significance analysis is a statistical method to distinguish meaningful (signal) factors from noise in genomic LZSS factorizations. By comparing factorization results from a real genome with results from a shuffled version of the same genome, we can determine the minimal factor length threshold at which factors are confidently considered non-random.

### Why is this important?

In genomic analysis, not all LZSS factors represent biologically meaningful repetitive structures. Short factors may arise simply by chance, while longer factors are more likely to represent true biological repeats, conserved regions, or other significant genomic features. This analysis provides a principled, statistically-grounded way to separate signal from noise.

### Key Concepts

- **Real Genome**: The original, unmodified genome sequence
- **Shuffled Genome**: A randomized version that preserves nucleotide composition but destroys biological structure
- **Factor Length**: The length of each non-overlapping LZSS factor
- **Significance Threshold (L\*)**: The minimum factor length at which factors are confidently non-random
- **Rarity Score**: Per-factor probability of observing that length in the shuffled genome

## Statistical Methodology

### Clopper-Pearson Confidence Bounds

The analysis uses **Clopper-Pearson** confidence intervals to establish conservative upper bounds on the probability of observing factors of a given length in shuffled data. This is a robust method for binomial proportions that provides exact coverage probabilities.

For k factors of length ≥ L out of n total factors in the shuffled genome:
- We compute the empirical tail probability: S₀(L) = k/n
- We compute a conservative upper bound: S₀ᵁ(L) using Clopper-Pearson

The Clopper-Pearson upper bound is the (1-α) quantile of the Beta(k+1, n-k) distribution.

### Threshold Selection

Given N factors in the real genome and a desired false positive threshold τ:
1. For each unique length L in the shuffled data, compute expected false positives: E[FP] = N × S₀ᵁ(L)
2. Select L\* as the smallest length where E[FP] ≤ τ

This ensures that we expect at most τ false positives among factors with length ≥ L\*.

### Per-Factor Rarity Scores

Each factor in the real genome receives a rarity score s_i = S₀(L_i), representing the tail probability of observing a factor of that length or greater in the shuffled genome. Lower scores indicate rarer (more significant) factors.

### Genome-Wide Exceedance Probability

Using a Poisson approximation, we estimate the probability of observing at least one factor of length ≥ L:

P(any factor ≥ L) ≈ 1 - exp(-N × S₀(L))

## API Reference

### Core Statistical Functions

#### `clopper_pearson_upper(k: int, n: int, alpha: float = 0.05) -> float`

Compute one-sided (1-α) Clopper-Pearson upper confidence bound for Binomial(n, p).

**Parameters:**
- `k` (int): Number of successes (0 ≤ k ≤ n)
- `n` (int): Number of trials (must be positive)
- `alpha` (float): Confidence level (default: 0.05 for 95% confidence)

**Returns:**
- `float`: Upper confidence bound for the binomial probability

**Raises:**
- `ValueError`: If inputs are invalid

**Example:**
```python
from noLZSS.genomics import clopper_pearson_upper

# 95% upper confidence bound for 5 successes out of 100 trials
upper = clopper_pearson_upper(5, 100, 0.05)
print(f"Upper bound: {upper:.4f}")  # ~0.11
```

---

#### `extract_factor_lengths(factors) -> np.ndarray`

Extract factor lengths from either a list of factor tuples or a binary factor file.

**Parameters:**
- `factors` (Union[List[Tuple], str, Path]): Either a list of (pos, length, ref, ...) tuples or path to binary factor file

**Returns:**
- `np.ndarray`: Array of factor lengths (int64)

**Raises:**
- `ValueError`: If input type is invalid or list contains invalid tuples
- `NoLZSSError`: If binary file cannot be read

**Example:**
```python
from noLZSS.genomics import extract_factor_lengths

# From tuple list
factors = [(0, 5, 0), (5, 3, 2), (8, 10, 1)]
lengths = extract_factor_lengths(factors)
print(lengths)  # [5, 3, 10]

# From binary file
lengths = extract_factor_lengths("genome_factors.bin")
print(f"Extracted {len(lengths)} factor lengths")
```

---

#### `infer_length_significance(real_lengths, shuf_lengths, tau_expected_fp=1.0, alpha_cp=0.05) -> Dict`

Core statistical inference using ONE shuffled genome.

**Parameters:**
- `real_lengths` (array-like): Factor lengths from real genome
- `shuf_lengths` (array-like): Factor lengths from shuffled genome
- `tau_expected_fp` (float): Maximum expected false positives (default: 1.0)
- `alpha_cp` (float): Confidence level for Clopper-Pearson (default: 0.05)

**Returns:**
Dictionary with keys:
- `N_real` (int): Number of real factors
- `N_shuf` (int): Number of shuffled factors
- `L_star` (int or None): Threshold length (None if no L meets criterion)
- `tau_expected_fp` (float): Input parameter
- `alpha_cp` (float): Input parameter
- `rarity_scores_real` (np.ndarray): Per-factor tail probabilities
- `p_any_ge` (Callable): Function(L) -> genome-wide exceedance probability
- `uniq_L` (np.ndarray): Unique length values
- `S0` (np.ndarray): Empirical tail probabilities
- `S0_upper` (np.ndarray): Conservative upper bounds
- `expected_fp_upper` (np.ndarray): Expected FP counts

**Example:**
```python
from noLZSS.genomics import infer_length_significance
import numpy as np

real_lens = np.array([50, 60, 70, 80, 90])
shuf_lens = np.array([5, 10, 15, 20, 25, 30])

result = infer_length_significance(real_lens, shuf_lens, tau_expected_fp=1.0)
print(f"Threshold length: {result['L_star']}")
print(f"Number of real factors: {result['N_real']}")

# Per-factor rarity scores
for i, score in enumerate(result['rarity_scores_real']):
    print(f"Factor {i} (len={real_lens[i]}): rarity={score:.4f}")

# Genome-wide exceedance probability
p_any_ge = result['p_any_ge']
print(f"P(any factor >= 50): {p_any_ge(50):.4f}")
```

---

### High-Level API

#### `calculate_factor_length_threshold(real_factors_file, shuffled_factors_file, tau_expected_fp=1.0, alpha_cp=0.05, plot_output=None) -> Dict`

Main user-facing function for significance analysis.

**Parameters:**
- `real_factors_file` (str or Path): Path to binary factor file from real genome
- `shuffled_factors_file` (str or Path): Path to binary factor file from shuffled genome
- `tau_expected_fp` (float): Maximum expected false positives (default: 1.0)
- `alpha_cp` (float): Confidence level (default: 0.05)
- `plot_output` (str, Path, or None): Optional path to save visualization

**Returns:**
- `dict`: Result dictionary from `infer_length_significance()`

**Raises:**
- `FileNotFoundError`: If either input file doesn't exist
- `NoLZSSError`: If files have invalid format

**Example:**
```python
from noLZSS.genomics import calculate_factor_length_threshold

result = calculate_factor_length_threshold(
    "genome_real_factors.bin",
    "genome_shuffled_factors.bin",
    tau_expected_fp=1.0,
    plot_output="significance_analysis.png"
)

if result['L_star'] is not None:
    print(f"Minimum significant factor length: {result['L_star']}")
    
    # Count significant factors
    significant = result['rarity_scores_real'] < 0.05
    print(f"Number of significant factors (p < 0.05): {significant.sum()}")
else:
    print("No threshold meets the false positive criterion")
```

---

### Visualization

#### `plot_significance_analysis(result, save_path=None, show_plot=True) -> None`

Create visualization of significance analysis results.

**Parameters:**
- `result` (dict): Output dictionary from `infer_length_significance()`
- `save_path` (str, Path, or None): Optional path to save plot
- `show_plot` (bool): Whether to display plot interactively (default: True)

**Raises:**
- `ValueError`: If result dictionary is missing required keys

**Warnings:**
- `UserWarning`: If matplotlib is not available

**Example:**
```python
from noLZSS.genomics import infer_length_significance, plot_significance_analysis

result = infer_length_significance(real_lens, shuf_lens)

# Display plot
plot_significance_analysis(result)

# Save without displaying
plot_significance_analysis(result, save_path="analysis.png", show_plot=False)
```

## Workflow Guide

### Step-by-Step Analysis

#### 1. Prepare Your Data

Start with a FASTA file containing your genome sequence(s):

```python
import noLZSS

# Factorize the real genome
noLZSS.write_factors_binary_file_fasta_multiple_dna_w_rc(
    "genome.fasta",
    "genome_real_factors.bin"
)
```

#### 2. Generate Shuffled Genome

Create a shuffled version that preserves nucleotide composition but destroys biological structure. You can use external tools or write your own:

```python
from pathlib import Path
import random

def shuffle_fasta(input_path, output_path):
    """Simple shuffler - preserves per-sequence composition."""
    from noLZSS.genomics import read_fasta_auto
    
    sequences = read_fasta_auto(input_path)
    
    with open(output_path, 'w') as f:
        for seq_id, seq_data in sequences:
            # Convert to list and shuffle
            seq_list = list(seq_data)
            random.shuffle(seq_list)
            shuffled = ''.join(seq_list)
            
            # Write FASTA
            f.write(f">{seq_id}_shuffled\n")
            # Write in 80-char lines
            for i in range(0, len(shuffled), 80):
                f.write(shuffled[i:i+80] + '\n')

shuffle_fasta("genome.fasta", "genome_shuffled.fasta")
```

#### 3. Factorize Shuffled Genome

```python
# Factorize the shuffled genome
noLZSS.write_factors_binary_file_fasta_multiple_dna_w_rc(
    "genome_shuffled.fasta",
    "genome_shuffled_factors.bin"
)
```

#### 4. Calculate Significance Threshold

```python
from noLZSS.genomics import calculate_factor_length_threshold

result = calculate_factor_length_threshold(
    "genome_real_factors.bin",
    "genome_shuffled_factors.bin",
    tau_expected_fp=1.0,
    alpha_cp=0.05,
    plot_output="significance_analysis.png"
)

print(f"Analysis Results:")
print(f"  Real genome factors: {result['N_real']}")
print(f"  Shuffled genome factors: {result['N_shuf']}")
print(f"  Minimum significant length (L*): {result['L_star']}")
```

#### 5. Filter Significant Factors

```python
from noLZSS.utils import read_factors_binary_file

# Read real factors
factors = read_factors_binary_file("genome_real_factors.bin")

# Filter by threshold
if result['L_star'] is not None:
    significant_factors = [f for f in factors if f[1] >= result['L_star']]
    print(f"Significant factors: {len(significant_factors)}/{len(factors)}")
    
    # Or use rarity scores for finer control
    rarity_threshold = 0.05  # p < 0.05
    highly_significant = [
        f for f, score in zip(factors, result['rarity_scores_real'])
        if score < rarity_threshold
    ]
    print(f"Highly significant (p < 0.05): {len(highly_significant)}")
```

## Parameter Guide

### Choosing `tau_expected_fp`

The `tau_expected_fp` parameter controls the stringency of the threshold:

- **τ = 0.1**: Very stringent, expect ≤0.1 false positives
- **τ = 1.0**: Standard choice, expect ≤1 false positive (default)
- **τ = 5.0**: Relaxed, expect ≤5 false positives

**Recommendations:**
- Use τ = 0.1 to 1.0 for genomes with thousands of factors
- Use τ = 5.0 to 10.0 for smaller analyses with hundreds of factors
- Consider your tolerance for false discoveries in downstream analysis

### Choosing `alpha_cp`

The `alpha_cp` parameter controls the confidence level for the Clopper-Pearson bounds:

- **α = 0.05**: 95% confidence (default, standard)
- **α = 0.01**: 99% confidence (more conservative)
- **α = 0.10**: 90% confidence (less conservative)

**Recommendations:**
- α = 0.05 is appropriate for most analyses
- Use smaller α (more conservative) when false positives are very costly
- Larger α provides less conservative bounds but may increase false positives

## Interpretation

### Understanding the Results

#### L\* (Threshold Length)

- **L\* = 150**: Factors with length ≥150 are confidently non-random
- **L\* = None**: No length meets the criterion; factors are similar to shuffled
- Lower L\* indicates more distinctive real genome structure
- Higher L\* suggests less distinctive long-range structure

#### Rarity Scores

- **Score = 0.001**: Very rare, only 0.1% of shuffled factors this long
- **Score = 0.05**: Rare, only 5% of shuffled factors this long
- **Score = 0.50**: Common, 50% of shuffled factors this long or longer
- **Score = 1.0**: Length shorter than any shuffled factor

#### Expected False Positives

The plot shows E[FP] = N_real × S₀ᵁ(L):
- Where the curve crosses τ is approximately L\*
- Steep drop indicates clear separation between real and shuffled
- Gradual drop suggests overlapping length distributions

### Biological Interpretation

**High L\* (e.g., L\* > 200)**
- Strong biological signal in long factors
- Clear distinction between structured and random sequences
- Likely represents true biological repeats, gene families, etc.

**Moderate L\* (e.g., 50 < L\* < 200)**
- Moderate biological structure
- Mix of structured and random-like regions
- May represent moderate repeat content

**Low L\* or None**
- Weak distinction from shuffled genome
- May indicate:
  - Low repeat content
  - High sequence diversity
  - Short, less distinctive repeats
  - Need for different null model

## Visualization Guide

### Understanding the Plots

The `plot_significance_analysis()` function creates two panels:

#### Top Panel: Tail Probabilities

- **Blue line (S₀(L))**: Empirical probability from shuffled data
- **Red dashed (S₀ᵁ(L))**: Conservative upper bound
- **Green dotted (L\*)**: Significance threshold
- Y-axis is log scale to show full range

**Interpretation:**
- Steep drop = clear length separation
- Flat region = many factors in that length range
- Gap between S₀ and S₀ᵁ = statistical uncertainty

#### Bottom Panel: Expected False Positives

- **Purple line**: Expected FP = N_real × S₀ᵁ(L)
- **Orange dashed**: τ threshold
- **Green dotted**: L\* where curves intersect

**Interpretation:**
- Where purple crosses orange = L\*
- Curve below τ = acceptable FP region
- Steeper drop = easier to find low L\*

## Examples

### Example 1: Basic Analysis

```python
import noLZSS
from noLZSS.genomics import calculate_factor_length_threshold

# Factorize real and shuffled genomes
noLZSS.write_factors_binary_file_fasta_multiple_dna_w_rc(
    "ecoli.fasta",
    "ecoli_real.bin"
)

noLZSS.write_factors_binary_file_fasta_multiple_dna_w_rc(
    "ecoli_shuffled.fasta",
    "ecoli_shuffled.bin"
)

# Calculate threshold
result = calculate_factor_length_threshold(
    "ecoli_real.bin",
    "ecoli_shuffled.bin",
    tau_expected_fp=1.0,
    plot_output="ecoli_significance.png"
)

print(f"E. coli genome analysis:")
print(f"  Threshold length: {result['L_star']}")
print(f"  Real factors: {result['N_real']}")
print(f"  Shuffled factors: {result['N_shuf']}")
```

### Example 2: Filtering Significant Factors

```python
from noLZSS.utils import read_factors_binary_file
import numpy as np

# Read factors and perform analysis
factors = read_factors_binary_file("genome_real.bin")
result = calculate_factor_length_threshold(
    "genome_real.bin",
    "genome_shuffled.bin"
)

# Multiple filtering strategies
if result['L_star'] is not None:
    # Strategy 1: Length threshold
    length_filtered = [f for f in factors if f[1] >= result['L_star']]
    print(f"Length threshold: {len(length_filtered)} factors")
    
    # Strategy 2: Rarity score (more granular)
    rarity_scores = result['rarity_scores_real']
    rarity_filtered = [
        f for f, score in zip(factors, rarity_scores)
        if score < 0.05  # p < 0.05
    ]
    print(f"Rarity filtered: {len(rarity_filtered)} factors")
    
    # Strategy 3: Top N most significant
    top_n = 100
    top_indices = np.argsort(rarity_scores)[:top_n]
    top_factors = [factors[i] for i in top_indices]
    print(f"Top {top_n} factors by rarity")
```

### Example 3: Multiple Genomes Comparison

```python
from pathlib import Path
import pandas as pd

genomes = ["ecoli", "human_chr1", "yeast"]
results = []

for genome in genomes:
    real_file = f"{genome}_real.bin"
    shuf_file = f"{genome}_shuffled.bin"
    
    if Path(real_file).exists() and Path(shuf_file).exists():
        result = calculate_factor_length_threshold(
            real_file,
            shuf_file,
            tau_expected_fp=1.0
        )
        
        results.append({
            'genome': genome,
            'N_real': result['N_real'],
            'N_shuf': result['N_shuf'],
            'L_star': result['L_star'],
            'pct_significant': (
                100 * (result['rarity_scores_real'] < 0.05).sum() / result['N_real']
            )
        })

# Create summary table
df = pd.DataFrame(results)
print(df.to_string(index=False))
```

### Example 4: Custom Analysis with Direct Functions

```python
from noLZSS.genomics import (
    extract_factor_lengths,
    infer_length_significance,
    plot_significance_analysis
)

# Extract lengths
real_lengths = extract_factor_lengths("genome_real.bin")
shuf_lengths = extract_factor_lengths("genome_shuffled.bin")

# Run analysis with custom parameters
result = infer_length_significance(
    real_lengths,
    shuf_lengths,
    tau_expected_fp=0.5,  # Stringent
    alpha_cp=0.01  # Very conservative
)

# Generate custom visualization
plot_significance_analysis(
    result,
    save_path="custom_analysis.png",
    show_plot=True
)

# Access detailed statistics
print(f"Unique lengths in shuffled: {len(result['uniq_L'])}")
print(f"Min/Max shuffled length: {result['uniq_L'].min()}/{result['uniq_L'].max()}")

# Compute custom statistics
import numpy as np
median_rarity = np.median(result['rarity_scores_real'])
print(f"Median rarity score: {median_rarity:.4f}")
```

## Advanced Topics

### Limitations and Considerations

1. **Single Shuffled Genome**: This is a "Tier-0" analysis using one shuffled replicate. For more robust results, consider generating multiple shuffled genomes and averaging.

2. **Shuffling Strategy**: The shuffling method matters. Per-sequence shuffling preserves local composition. Whole-genome shuffling may be more appropriate for some applications.

3. **Length-Only Analysis**: This analysis considers only factor length, not position or sequence content. Factors at the same length are treated identically.

4. **Conservative Bounds**: Clopper-Pearson provides conservative (wide) confidence intervals. This is by design to avoid false significance claims.

### Performance Notes

- Memory usage is O(N_factors) for storing lengths
- Computation time is dominated by reading binary files
- Plot generation requires matplotlib (optional dependency)
- Exact Clopper-Pearson requires scipy (falls back to Wilson score if unavailable)

### Integration with Other Tools

This analysis integrates naturally with:
- Factor position analysis (see `noLZSS.genomics.plots`)
- Strand bias analysis (see strand-bias-heatmap)
- Per-sequence factorization workflows
- Downstream filtering pipelines

## References

1. Clopper, C.J. and Pearson, E.S. (1934). "The use of confidence or fiducial limits illustrated in the case of the binomial." Biometrika, 26(4), 404-413.

2. Wilson, E.B. (1927). "Probable inference, the law of succession, and statistical inference." Journal of the American Statistical Association, 22(158), 209-212.

3. Köppl, D. (2021). "Non-Overlapping LZ77 Factorization and LZ78 Substring Compression Queries with Suffix Trees." Algorithms, 14(2), 44.
