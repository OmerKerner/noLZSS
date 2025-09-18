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
        chars = data.decode('ascii')
        char_set = set(chars)
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


def read_factors_binary_file_with_metadata(filepath: Union[str, Path]) -> Dict[str, Any]:
    """
    Read factors from an enhanced binary file with metadata (sequence names and sentinel indices).
    
    This function reads binary files written by write_factors_binary_file_fasta_multiple_dna_*
    functions that contain metadata including sequence names and sentinel factor indices.
    
    Args:
        filepath: Path to the binary factors file with metadata
        
    Returns:
        Dictionary containing:
        - 'factors': List of (start, length, ref, is_rc) tuples 
        - 'sentinel_factor_indices': List of factor indices that are sentinels
        - 'sequence_names': List of sequence names from FASTA headers
        - 'num_sequences': Number of sequences
        - 'num_sentinels': Number of sentinel factors
        
    Raises:
        NoLZSSError: If file cannot be read or has invalid format
    """
    # RC_MASK is the MSB of uint64_t (defined in C++ as 1ULL << 63)
    RC_MASK = 1 << 63
    
    filepath = Path(filepath)
    if not filepath.exists():
        raise NoLZSSError(f"File not found: {filepath}")
    
    try:
        with open(filepath, 'rb') as f:
            # Read header
            header_data = f.read(48)  # FactorFileHeader is 8+8+8+8+8 = 40 bytes + padding
            if len(header_data) < 40:
                raise NoLZSSError("File too small to contain valid header")
            
            # Unpack header (magic is 8 chars, then 4 uint64_t)
            magic = header_data[:8]
            if magic != b'noLZSSv1':
                raise NoLZSSError("Invalid file format: missing noLZSS magic header")
            
            num_factors, num_sequences, num_sentinels, header_size = struct.unpack('<QQQQ', header_data[8:40])
            
            # Seek to beginning of header to read the full header
            f.seek(0)
            full_header = f.read(header_size)
            if len(full_header) != header_size:
                raise NoLZSSError(f"Could not read full header: expected {header_size}, got {len(full_header)}")
            
            # Skip the basic header structure
            offset = 40
            
            # Read sequence names
            sequence_names = []
            for i in range(num_sequences):
                # Find null terminator
                name_start = offset
                while offset < len(full_header) and full_header[offset] != 0:
                    offset += 1
                if offset >= len(full_header):
                    raise NoLZSSError("Invalid sequence name format")
                
                name = full_header[name_start:offset].decode('utf-8')
                sequence_names.append(name)
                offset += 1  # Skip null terminator
            
            # Read sentinel factor indices
            sentinel_indices = []
            for i in range(num_sentinels):
                if offset + 8 > len(full_header):
                    raise NoLZSSError("Insufficient data for sentinel indices")
                
                idx = struct.unpack('<Q', full_header[offset:offset+8])[0]
                sentinel_indices.append(idx)
                offset += 8
            
            # Read factors
            factors = []
            for i in range(num_factors):
                factor_data = f.read(24)  # Each factor is 3 * uint64_t = 24 bytes
                if len(factor_data) != 24:
                    raise NoLZSSError(f"Insufficient data for factor {i}")
                
                start, length, ref = struct.unpack('<QQQ', factor_data)
                
                # Extract is_rc flag and clean ref
                is_rc_flag = bool(ref & RC_MASK)
                clean_ref = ref & ~RC_MASK
                
                factors.append((start, length, clean_ref, is_rc_flag))
    
    except IOError as e:
        raise NoLZSSError(f"Error reading file {filepath}: {e}")
    except struct.error as e:
        raise NoLZSSError(f"Error unpacking binary data: {e}")
    
    return {
        'factors': factors,
        'sentinel_factor_indices': sentinel_indices,
        'sequence_names': sequence_names,
        'num_sequences': num_sequences,
        'num_sentinels': num_sentinels
    }


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
        NoLZSSError: If binary file cannot be read
        TypeError: If input type is invalid
        ValueError: If no factors to plot
        
    Warns:
        UserWarning: If matplotlib is not installed (function returns gracefully)
    """
    # Validate input and get factors BEFORE trying to import matplotlib
    if isinstance(factors_or_file, (str, Path)):
        factors = read_factors_binary_file(factors_or_file)
    elif isinstance(factors_or_file, list):
        factors = factors_or_file
    else:
        raise TypeError("factors_or_file must be a list of tuples or a path to a binary file")
    
    if not factors:
        raise ValueError("No factors to plot")
    
    # Now try to import matplotlib for plotting
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        warnings.warn("matplotlib is required for plotting. Install with: pip install matplotlib", UserWarning)
        return
    
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
