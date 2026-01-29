"""
Statistical significance analysis for factor lengths.

This module provides tools to determine the minimal factor length threshold from which
factors are confidently considered signal (not noise) by comparing real genome 
factorization results with shuffled genome factorization results.

The core methodology uses Clopper-Pearson confidence bounds for binomial tail 
probabilities to establish conservative thresholds while controlling for false positives.
"""

from typing import Union, Dict, Any, List, Tuple, Optional, Callable
import warnings
from pathlib import Path

import numpy as np

from ..utils import read_factors_binary_file, NoLZSSError


def clopper_pearson_upper(k: int, n: int, alpha: float = 0.05) -> float:
    """
    Compute one-sided (1-alpha) Clopper-Pearson upper confidence bound for Binomial(n, p).
    
    This provides a conservative upper bound on the true binomial success probability
    given k successes out of n trials.
    
    Args:
        k: Number of successes (must be 0 <= k <= n)
        n: Number of trials (must be positive)
        alpha: Confidence level (default: 0.05 for 95% confidence)
        
    Returns:
        Upper confidence bound for the binomial probability p
        
    Raises:
        ValueError: If inputs are invalid (k < 0, k > n, n <= 0, alpha not in (0,1))
        
    Examples:
        >>> clopper_pearson_upper(5, 100, 0.05)
        0.11...
        >>> clopper_pearson_upper(0, 100, 0.05)
        0.03...
    """
    # Input validation
    if n <= 0:
        raise ValueError(f"n must be positive, got {n}")
    if k < 0 or k > n:
        raise ValueError(f"k must be between 0 and n, got k={k}, n={n}")
    if alpha <= 0 or alpha >= 1:
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")
    
    # Edge cases
    if k == n:
        return 1.0
    if k == 0:
        # For k=0, CP upper bound is 1 - alpha^(1/n)
        return 1.0 - (alpha ** (1.0 / n))
    
    # Try to use scipy for exact Clopper-Pearson bound
    try:
        from scipy.stats import beta
        # CP upper bound is the (1-alpha) quantile of Beta(k+1, n-k)
        upper = beta.ppf(1.0 - alpha, k + 1, n - k)
        return float(upper)
    except ImportError:
        # Fallback to Wilson score interval (conservative approximation)
        warnings.warn(
            "scipy not available, using Wilson score approximation for Clopper-Pearson bound. "
            "Install scipy for exact bounds: pip install scipy",
            UserWarning
        )
        
        # Wilson score upper bound
        from math import sqrt
        z = 1.645 if alpha == 0.05 else 1.96 if alpha == 0.025 else 2.576  # Approximate z-scores
        
        p_hat = k / n
        denominator = 1 + z**2 / n
        center = (p_hat + z**2 / (2*n)) / denominator
        margin = z * sqrt((p_hat * (1 - p_hat) / n + z**2 / (4*n**2))) / denominator
        
        upper = min(center + margin, 1.0)
        return float(upper)


def extract_factor_lengths(
    factors: Union[List[Tuple[int, ...]], str, Path]
) -> np.ndarray:
    """
    Extract factor lengths from either a list of factor tuples or a binary factor file.
    
    Args:
        factors: Either a list of (pos, length, ref, ...) tuples or path to binary factor file
        
    Returns:
        numpy array of factor lengths (int64)
        
    Raises:
        ValueError: If input type is invalid or list contains invalid tuples
        NoLZSSError: If binary file cannot be read or has invalid format
        
    Examples:
        >>> factors = [(0, 5, 0), (5, 3, 2), (8, 10, 1)]
        >>> lengths = extract_factor_lengths(factors)
        >>> list(lengths)
        [5, 3, 10]
        
        >>> lengths = extract_factor_lengths("genome_factors.bin")
        >>> len(lengths)
        1000
    """
    if isinstance(factors, (str, Path)):
        # Read from binary file
        factors_list = read_factors_binary_file(factors)
        if not factors_list:
            return np.array([], dtype=np.int64)
        # Extract length (second element) from each tuple
        lengths = np.array([f[1] for f in factors_list], dtype=np.int64)
        return lengths
    
    elif isinstance(factors, list):
        if not factors:
            return np.array([], dtype=np.int64)
        
        # Validate that each element is a tuple with at least 2 elements
        for i, factor in enumerate(factors):
            if not isinstance(factor, tuple) or len(factor) < 2:
                raise ValueError(
                    f"Factor at index {i} must be a tuple with at least 2 elements "
                    f"(pos, length, ...), got {type(factor)}"
                )
        
        # Extract length (second element) from each tuple
        lengths = np.array([f[1] for f in factors], dtype=np.int64)
        return lengths
    
    else:
        raise ValueError(
            f"factors must be a list of tuples or a file path, got {type(factors)}"
        )


def infer_length_significance(
    real_lengths: Union[np.ndarray, List[int]],
    shuf_lengths: Union[np.ndarray, List[int]],
    tau_expected_fp: float = 1.0,
    alpha_cp: float = 0.05
) -> Dict[str, Any]:
    """
    Tier-0, length-only inference using ONE shuffled genome.
    
    Computes empirical tail CCDF from shuffled factor lengths: S0(L) = P0(len >= L)
    Computes conservative upper bound S0^U(L) via Clopper-Pearson.
    Chooses L* so that N_real * S0^U(L*) <= tau_expected_fp.
    
    Provides per-factor rarity score s_i = S0(L_i) and approximate genome-wide
    exceedance probability via Poisson approximation.
    
    Args:
        real_lengths: Array of factor lengths from real genome
        shuf_lengths: Array of factor lengths from shuffled genome  
        tau_expected_fp: Maximum expected false positives (default: 1.0)
        alpha_cp: Confidence level for Clopper-Pearson bounds (default: 0.05)
        
    Returns:
        Dictionary containing:
        - N_real: Number of real factors
        - N_shuf: Number of shuffled factors
        - L_star: Threshold length (None if no L meets criterion)
        - tau_expected_fp: Input parameter
        - alpha_cp: Input parameter
        - rarity_scores_real: Per-factor tail probabilities (length N_real)
        - p_any_ge: Function(L) -> genome-wide exceedance probability
        - uniq_L: Unique length values from shuffled data
        - S0: Empirical tail probabilities for each unique L
        - S0_upper: Conservative upper bounds for each unique L
        - expected_fp_upper: Expected false positive counts (N_real * S0_upper)
        
    Examples:
        >>> real_lens = np.array([5, 10, 15, 20, 25])
        >>> shuf_lens = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
        >>> result = infer_length_significance(real_lens, shuf_lens, tau_expected_fp=0.5)
        >>> result['L_star']  # Minimum significant length
        20
        >>> result['rarity_scores_real'][:3]  # Rarity of first 3 factors
        array([0.44..., 0.11..., 0.0])
    """
    # Convert to numpy arrays
    real_lengths = np.asarray(real_lengths, dtype=np.int64)
    shuf_lengths = np.asarray(shuf_lengths, dtype=np.int64)
    
    N_real = len(real_lengths)
    N_shuf = len(shuf_lengths)
    
    if N_shuf == 0:
        raise ValueError("Shuffled genome must have at least one factor")
    
    # Get unique lengths from shuffled data and sort
    uniq_L = np.unique(shuf_lengths)
    
    # Compute empirical tail probabilities S0(L) = P(len >= L) from shuffled data
    S0 = np.array([np.sum(shuf_lengths >= L) / N_shuf for L in uniq_L])
    
    # Compute Clopper-Pearson upper bounds S0^U(L)
    S0_upper = np.array([
        clopper_pearson_upper(int(np.sum(shuf_lengths >= L)), N_shuf, alpha_cp)
        for L in uniq_L
    ])
    
    # Expected false positive counts under the conservative bound
    expected_fp_upper = N_real * S0_upper
    
    # Find L* such that expected_fp_upper[L*] <= tau_expected_fp
    L_star = None
    valid_indices = np.where(expected_fp_upper <= tau_expected_fp)[0]
    if len(valid_indices) > 0:
        # Take the smallest L that meets the criterion
        L_star = int(uniq_L[valid_indices[0]])
    
    # Compute per-factor rarity scores: s_i = S0(L_i)
    # Use linear interpolation for lengths not in uniq_L
    # For L < min(uniq_L), use S0 = 1.0
    # For L > max(uniq_L), use S0 = 0.0
    rarity_scores_real = np.interp(
        real_lengths,
        uniq_L,
        S0,
        left=1.0,  # Lengths shorter than any shuffled factor
        right=0.0  # Lengths longer than any shuffled factor
    )
    
    # Genome-wide exceedance probability via Poisson approximation
    # P(at least one factor >= L) ≈ 1 - exp(-N_real * S0(L))
    def p_any_ge(L: float) -> float:
        """
        Approximate probability that at least one factor has length >= L.
        
        Uses Poisson approximation: P(X >= 1) ≈ 1 - exp(-lambda) where lambda = N_real * S0(L)
        """
        s0_L = np.interp(L, uniq_L, S0, left=1.0, right=0.0)
        lambda_val = N_real * s0_L
        return 1.0 - np.exp(-lambda_val)
    
    return {
        'N_real': N_real,
        'N_shuf': N_shuf,
        'L_star': L_star,
        'tau_expected_fp': tau_expected_fp,
        'alpha_cp': alpha_cp,
        'rarity_scores_real': rarity_scores_real,
        'p_any_ge': p_any_ge,
        'uniq_L': uniq_L,
        'S0': S0,
        'S0_upper': S0_upper,
        'expected_fp_upper': expected_fp_upper,
    }


def calculate_factor_length_threshold(
    real_factors_file: Union[str, Path],
    shuffled_factors_file: Union[str, Path],
    tau_expected_fp: float = 1.0,
    alpha_cp: float = 0.05,
    plot_output: Optional[Union[str, Path]] = None
) -> Dict[str, Any]:
    """
    Calculate the minimum significant factor length threshold.
    
    Main user-facing function that reads binary factor files, extracts lengths,
    and calls infer_length_significance() to determine the significance threshold.
    Optionally creates a visualization plot.
    
    Args:
        real_factors_file: Path to binary factor file from real genome
        shuffled_factors_file: Path to binary factor file from shuffled genome
        tau_expected_fp: Maximum expected false positives (default: 1.0)
        alpha_cp: Confidence level for Clopper-Pearson bounds (default: 0.05)
        plot_output: Optional path to save visualization plot
        
    Returns:
        Result dictionary from infer_length_significance() with all computed statistics
        
    Raises:
        FileNotFoundError: If either input file doesn't exist
        NoLZSSError: If files cannot be read or have invalid format
        
    Examples:
        >>> result = calculate_factor_length_threshold(
        ...     "genome_factors.bin",
        ...     "shuffled_factors.bin",
        ...     tau_expected_fp=1.0
        ... )
        >>> print(f"Minimum significant length: {result['L_star']}")
        Minimum significant length: 150
        
        >>> # With visualization
        >>> result = calculate_factor_length_threshold(
        ...     "genome_factors.bin",
        ...     "shuffled_factors.bin",
        ...     plot_output="significance.png"
        ... )
    """
    # Validate file paths
    real_path = Path(real_factors_file)
    shuf_path = Path(shuffled_factors_file)
    
    if not real_path.exists():
        raise FileNotFoundError(f"Real factors file not found: {real_path}")
    if not shuf_path.exists():
        raise FileNotFoundError(f"Shuffled factors file not found: {shuf_path}")
    
    # Extract factor lengths from both files
    real_lengths = extract_factor_lengths(real_path)
    shuf_lengths = extract_factor_lengths(shuf_path)
    
    # Perform significance analysis
    result = infer_length_significance(
        real_lengths,
        shuf_lengths,
        tau_expected_fp=tau_expected_fp,
        alpha_cp=alpha_cp
    )
    
    # Create plot if requested
    if plot_output is not None:
        plot_significance_analysis(result, save_path=plot_output, show_plot=False)
    
    return result


def plot_significance_analysis(
    result: Dict[str, Any],
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True
) -> None:
    """
    Create visualization of significance analysis results.
    
    Generates a multi-panel plot showing:
    - S0(L): Empirical tail probability from shuffled data
    - S0^U(L): Conservative upper bound via Clopper-Pearson
    - Expected FP curve: N_real * S0^U(L)
    - Vertical line at L* threshold (if found)
    - Horizontal line at tau_expected_fp
    
    Args:
        result: Output dictionary from infer_length_significance()
        save_path: Optional path to save the plot (e.g., 'significance.png')
        show_plot: Whether to display the plot interactively (default: True)
        
    Raises:
        ValueError: If result dictionary is missing required keys
        
    Warnings:
        UserWarning: If matplotlib is not available (function returns gracefully)
        
    Examples:
        >>> result = infer_length_significance(real_lens, shuf_lens)
        >>> plot_significance_analysis(result, save_path="analysis.png")
        
        >>> # Display without saving
        >>> plot_significance_analysis(result)
    """
    # Validate required keys
    required_keys = ['uniq_L', 'S0', 'S0_upper', 'expected_fp_upper', 
                     'L_star', 'tau_expected_fp', 'N_real', 'N_shuf']
    missing_keys = [k for k in required_keys if k not in result]
    if missing_keys:
        raise ValueError(f"Result dictionary missing required keys: {missing_keys}")
    
    # Try to import matplotlib
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        warnings.warn(
            "matplotlib is required for plotting. Install with: pip install matplotlib",
            UserWarning
        )
        return
    
    # Extract data
    uniq_L = result['uniq_L']
    S0 = result['S0']
    S0_upper = result['S0_upper']
    expected_fp_upper = result['expected_fp_upper']
    L_star = result['L_star']
    tau_expected_fp = result['tau_expected_fp']
    N_real = result['N_real']
    N_shuf = result['N_shuf']
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Top panel: Tail probabilities
    ax1.semilogy(uniq_L, S0, 'b-', linewidth=2, label='S0(L) - Empirical')
    ax1.semilogy(uniq_L, S0_upper, 'r--', linewidth=2, label='S0^U(L) - Upper bound')
    
    if L_star is not None:
        ax1.axvline(L_star, color='green', linestyle=':', linewidth=2, 
                   label=f'L* = {L_star}')
    
    ax1.set_ylabel('Tail Probability P(len ≥ L)', fontsize=12)
    ax1.set_title(
        f'Factor Length Significance Analysis\n'
        f'N_real = {N_real}, N_shuf = {N_shuf}, τ = {tau_expected_fp}',
        fontsize=13
    )
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Bottom panel: Expected false positives
    ax2.plot(uniq_L, expected_fp_upper, 'purple', linewidth=2, 
            label='Expected FP (upper)')
    ax2.axhline(tau_expected_fp, color='orange', linestyle='--', linewidth=2,
               label=f'τ = {tau_expected_fp}')
    
    if L_star is not None:
        ax2.axvline(L_star, color='green', linestyle=':', linewidth=2,
                   label=f'L* = {L_star}')
    
    ax2.set_xlabel('Factor Length L', fontsize=12)
    ax2.set_ylabel('Expected False Positives', fontsize=12)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save if requested
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    
    # Show if requested
    if show_plot:
        plt.show()
    else:
        plt.close()
