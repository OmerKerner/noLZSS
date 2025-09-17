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


class FactorizationData:
    """Complete factorization data from extended binary file."""
    
    def __init__(self, factors: List[Tuple[int, int, int]], 
                 sequence_ids: List[str], 
                 sentinel_factor_indices: List[int],
                 sequence_lengths: List[int],
                 sequence_positions: List[int]):
        self.factors = factors
        self.sequence_ids = sequence_ids
        self.sentinel_factor_indices = set(sentinel_factor_indices)
        self.sequence_lengths = sequence_lengths
        self.sequence_positions = sequence_positions
    
    def get_sequence_factors(self) -> List[Tuple[int, int, int]]:
        """Extract only sequence factors (excluding sentinels)."""
        return [factor for i, factor in enumerate(self.factors) 
                if i not in self.sentinel_factor_indices]
    
    def get_sentinel_factors(self) -> List[Tuple[int, int, int]]:
        """Extract only sentinel factors."""
        return [factor for i, factor in enumerate(self.factors) 
                if i in self.sentinel_factor_indices]


def read_factors_binary_file_extended(filepath: Union[str, Path]) -> FactorizationData:
    """
    Read factorization data from extended binary file with metadata.
    
    Args:
        filepath: Path to extended binary factor file
        
    Returns:
        FactorizationData with factors and metadata
        
    Raises:
        NoLZSSError: If file format is invalid or corrupted
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise NoLZSSError(f"Factor file not found: {filepath}")
    
    with open(filepath, 'rb') as f:
        # Read and validate header
        magic = f.read(8)
        if magic != b'noLZSSv1':
            raise NoLZSSError(f"Invalid file format. Expected noLZSS binary format.")
        
        # Read header fields
        header_data = f.read(32)  # 4 * uint64_t + 3 * uint64_t reserved
        if len(header_data) < 32:
            raise NoLZSSError("Incomplete header in extended binary file")
        
        num_factors, num_sequences, num_sentinels, header_size = struct.unpack('<QQQQ', header_data[:32])
        reserved = struct.unpack('<QQQ', f.read(24))  # Skip reserved fields
        
        # Read sequence names
        sequence_ids = []
        for _ in range(num_sequences):
            name_bytes = b''
            while True:
                byte = f.read(1)
                if not byte or byte == b'\0':
                    break
                name_bytes += byte
            sequence_ids.append(name_bytes.decode('utf-8'))
        
        # Read sentinel indices
        sentinel_indices = []
        for _ in range(num_sentinels):
            idx_data = f.read(8)
            if len(idx_data) < 8:
                raise NoLZSSError("Incomplete sentinel indices in file")
            idx, = struct.unpack('<Q', idx_data)
            sentinel_indices.append(int(idx))
        
        # Read sequence lengths
        sequence_lengths = []
        for _ in range(num_sequences):
            len_data = f.read(8)
            if len(len_data) < 8:
                raise NoLZSSError("Incomplete sequence lengths in file")
            length, = struct.unpack('<Q', len_data)
            sequence_lengths.append(int(length))
        
        # Read sequence positions
        sequence_positions = []
        for _ in range(num_sequences):
            pos_data = f.read(8)
            if len(pos_data) < 8:
                raise NoLZSSError("Incomplete sequence positions in file")
            position, = struct.unpack('<Q', pos_data)
            sequence_positions.append(int(position))
        
        # Read factors
        factors = []
        factor_size = 24  # 3 * uint64_t
        for _ in range(num_factors):
            factor_data = f.read(factor_size)
            if len(factor_data) != factor_size:
                raise NoLZSSError("Incomplete factor data in file")
            
            start_pos, length, reference = struct.unpack('<QQQ', factor_data)
            factors.append((int(start_pos), int(length), int(reference)))
    
    return FactorizationData(
        factors=factors,
        sequence_ids=sequence_ids,
        sentinel_factor_indices=sentinel_indices,
        sequence_lengths=sequence_lengths,
        sequence_positions=sequence_positions
    )


def detect_binary_file_format(filepath: Union[str, Path]) -> str:
    """
    Detect the format of a binary factor file.
    
    Args:
        filepath: Path to binary factor file
        
    Returns:
        'extended' for extended format with metadata, 'legacy' for standard format
        
    Raises:
        NoLZSSError: If file cannot be read or format is invalid
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise NoLZSSError(f"File not found: {filepath}")
    
    try:
        with open(filepath, 'rb') as f:
            # Check for extended format magic
            magic = f.read(8)
            if magic == b'noLZSSv1':
                return 'extended'
            else:
                # Reset and check if it's a valid legacy format (multiple of 24 bytes)
                f.seek(0)
                file_size = f.seek(0, 2)  # Seek to end to get size
                if file_size % 24 == 0:
                    return 'legacy'
                else:
                    raise NoLZSSError(f"Invalid binary file format: unrecognized magic and size not multiple of 24")
    except IOError as e:
        raise NoLZSSError(f"Error reading file {filepath}: {e}")


def read_factors_binary_file_auto(filepath: Union[str, Path]) -> Union[List[Tuple[int, int, int]], FactorizationData]:
    """
    Automatically detect and read binary factor file in any supported format.
    
    Args:
        filepath: Path to binary factor file
        
    Returns:
        For legacy format: List of (position, length, ref) tuples
        For extended format: FactorizationData object with metadata
        
    Raises:
        NoLZSSError: If file format is invalid or corrupted
    """
    format_type = detect_binary_file_format(filepath)
    
    if format_type == 'extended':
        return read_factors_binary_file_extended(filepath)
    else:  # legacy
        return read_factors_binary_file(filepath)


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
