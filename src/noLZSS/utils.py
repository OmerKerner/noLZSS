"""
Utility functions for input validation, alphabet analysis, file I/O helpers, and visualization.

This module provides reusable utilities for the noLZSS package, including
input validation, sentinel handling, alphabet analysis, binary file I/O, and plotting functions.
"""

from typing import Union, Dict, Any, List, Tuple, Optional
import math
import struct
from pathlib import Path
from collections import Counter
import warnings


class NoLZSSError(Exception):
    """Base exception for noLZSS-related errors."""
    pass


class InvalidInputError(NoLZSSError):
    """Raised when input data is invalid for factorization."""
    pass


def validate_input(data: Union[str, bytes]) -> bytes:
    """
    Validate and normalize input data for factorization.
    
    Args:
        data: Input string or bytes to validate
        
    Returns:
        Normalized bytes data
        
    Raises:
        InvalidInputError: If input is invalid
        TypeError: If input type is not supported
    """
    if isinstance(data, str):
        # Convert string to bytes using ASCII encoding (1 byte per char)
        try:
            data = data.encode('ascii')
        except UnicodeEncodeError as e:
            raise InvalidInputError(f"Input string must contain only ASCII characters (1 byte each): {e}")
    elif isinstance(data, bytes):
        pass  # Already bytes
    else:
        raise TypeError(f"Input must be str or bytes, got {type(data)}")
    
    if len(data) == 0:
        raise InvalidInputError("Input data cannot be empty")
    
    # Check for null bytes in the middle (which might interfere with C++ processing)
    if b'\x00' in data[:-1]:  # Allow null byte only at the end (as potential sentinel)
        raise InvalidInputError("Input data contains null bytes")
    
    return data


def analyze_alphabet(data: Union[str, bytes]) -> Dict[str, Any]:
    """
    Analyze the alphabet of input data.
    
    Args:
        data: Input string or bytes to analyze
        
    Returns:
        Dictionary containing alphabet analysis:
        - 'size': Number of unique characters/bytes
        - 'characters': Set of unique characters/bytes
        - 'distribution': Counter of character/byte frequencies
        - 'entropy': Shannon entropy of the data
        - 'most_common': List of (char, count) tuples for most frequent characters
    """
    if isinstance(data, str):
        chars = data
        char_set = set(data)
    elif isinstance(data, bytes):
        chars = data
        char_set = set(data)
    else:
        raise TypeError(f"Input must be str or bytes, got {type(data)}")
    
    distribution = Counter(chars)
    total_chars = len(chars)
    
    # Calculate Shannon entropy
    entropy = 0.0
    if total_chars > 0:
        for count in distribution.values():
            if count > 0:
                p = count / total_chars
                entropy -= p * math.log2(p)
    
    return {
        'size': len(char_set),
        'characters': char_set,
        'distribution': distribution,
        'entropy': entropy,
        'most_common': distribution.most_common(10),  # Top 10 most frequent
        'total_length': total_chars
    }


def read_factors_binary_file(filepath: Union[str, Path]) -> List[Tuple[int, int, int]]:
    """
    Read factors from a binary file written by write_factors_binary_file.
    
    Args:
        filepath: Path to the binary factors file
        
    Returns:
        List of (position, length, ref) tuples
        
    Raises:
        NoLZSSError: If file cannot be read or has invalid format
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise NoLZSSError(f"File not found: {filepath}")
    
    try:
        with open(filepath, 'rb') as f:
            binary_data = f.read()
    except IOError as e:
        raise NoLZSSError(f"Error reading file {filepath}: {e}")
    
    if len(binary_data) % 24 != 0:
        raise NoLZSSError(f"Invalid binary file format: file size {len(binary_data)} is not a multiple of 24")
    
    factors = []
    for i in range(len(binary_data) // 24):
        start, length, ref = struct.unpack('<QQQ', binary_data[i*24:(i+1)*24])
        factors.append((start, length, ref))
    
    return factors


def plot_factor_lengths(
    factors_or_file: Union[List[Tuple[int, int, int]], str, Path],
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True
) -> None:
    """
    Plot the cumulative factor lengths vs factor index.
    
    Creates a scatter plot where:
    - X-axis: Cumulative sum of factor lengths
    - Y-axis: Factor index (number of factors)
    
    Args:
        factors_or_file: Either a list of (position, length, ref) tuples or path to binary factors file
        save_path: Optional path to save the plot image (e.g., 'plot.png')
        show_plot: Whether to display the plot (default: True)
        
    Raises:
        ImportError: If matplotlib is not installed
        NoLZSSError: If binary file cannot be read
        TypeError: If input type is invalid
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        warnings.warn("matplotlib is required for plotting. Install with: pip install matplotlib", UserWarning)
        return
    
    # Get factors
    if isinstance(factors_or_file, (str, Path)):
        factors = read_factors_binary_file(factors_or_file)
    elif isinstance(factors_or_file, list):
        factors = factors_or_file
    else:
        raise TypeError("factors_or_file must be a list of tuples or a path to a binary file")
    
    if not factors:
        raise ValueError("No factors to plot")
    
    # Compute cumulative lengths
    cumulative_lengths = []
    current_sum = 0
    for i, (_, length, _) in enumerate(factors):
        current_sum += length
        cumulative_lengths.append((i + 1, current_sum))  # y = factor index (1-based), x = cumulative
    
    # Extract x and y
    y_values, x_values = zip(*cumulative_lengths)
    
    # Create step (staircase) plot
    plt.figure(figsize=(10, 6))
    # 'where' controls alignment: 'post' holds the value until the next x,
    # 'pre' jumps before the x, 'mid' centers the step.
    plt.step(x_values, y_values, where='post', linewidth=1.5)
    # optional: show points at the step locations
    plt.plot(x_values, y_values, linestyle='', marker='o', markersize=4, alpha=0.6)
    plt.xlabel('Cumulative Factor Length')
    plt.ylabel('Factor Index')
    plt.title('Factor Length Accumulation (Step Plot)')
    plt.grid(True, alpha=0.3)
    
    # Save if requested
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    
    # Show plot
    if show_plot:
        plt.show()
