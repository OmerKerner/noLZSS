"""
Utility functions for input validation, alphabet analysis, and file I/O helpers.

This module provides reusable utilities for the noLZSS package, including
input validation, sentinel handling, and alphabet analysis.
"""

from typing import Union, Dict, Any, Set
import re
import math
from pathlib import Path
from collections import Counter


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
        # Convert string to bytes using UTF-8 encoding
        try:
            data = data.encode('utf-8')
        except UnicodeEncodeError as e:
            raise InvalidInputError(f"Unable to encode string to UTF-8: {e}")
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


def ensure_sentinel(data: bytes, sentinel: bytes = b'$') -> bytes:
    """
    Ensure input data ends with a sentinel character.
    
    Args:
        data: Input bytes data
        sentinel: Sentinel character to append (default: '$')
        
    Returns:
        Data with sentinel appended if not already present
    """
    if not data.endswith(sentinel):
        data = data + sentinel
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


def is_dna_sequence(data: Union[str, bytes]) -> bool:
    """
    Check if data appears to be a DNA sequence (A, T, G, C, N).
    
    Args:
        data: Input data to check
        
    Returns:
        True if data contains only DNA nucleotides (case insensitive)
    """
    if isinstance(data, bytes):
        try:
            data = data.decode('ascii')
        except UnicodeDecodeError:
            return False
    
    # At this point data is guaranteed to be a string
    if not isinstance(data, str):
        return False
        
    # Allow standard DNA bases plus N (unknown), case insensitive
    # Remove sentinel character before checking
    clean_data = data.rstrip('$')
    dna_pattern = re.compile(r'^[ATGCN]+$', re.IGNORECASE)
    return bool(dna_pattern.match(clean_data))


def is_protein_sequence(data: Union[str, bytes]) -> bool:
    """
    Check if data appears to be a protein sequence (20 standard amino acids).
    
    Args:
        data: Input data to check
        
    Returns:
        True if data contains only standard amino acid codes
    """
    if isinstance(data, bytes):
        try:
            data = data.decode('ascii')
        except UnicodeDecodeError:
            return False
    
    # At this point data is guaranteed to be a string
    if not isinstance(data, str):
        return False
    
    # Standard 20 amino acids plus common extensions (B, J, O, U, X, Z)
    # Remove sentinel character before checking
    clean_data = data.rstrip('$')
    protein_pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWYBJOUX]+$', re.IGNORECASE)
    return bool(protein_pattern.match(clean_data))


def detect_sequence_type(data: Union[str, bytes]) -> str:
    """
    Detect the likely type of biological sequence.
    
    Args:
        data: Input data to analyze
        
    Returns:
        String indicating sequence type: 'dna', 'protein', 'text', or 'binary'
    """
    if isinstance(data, bytes):
        # Check if it's ASCII text
        try:
            text_data = data.decode('ascii')
        except UnicodeDecodeError:
            return 'binary'
        data = text_data
    
    if is_dna_sequence(data):
        return 'dna'
    elif is_protein_sequence(data):
        return 'protein'
    elif isinstance(data, str):
        return 'text'
    else:
        return 'binary'


def safe_file_reader(filepath: Union[str, Path], chunk_size: int = 8192):
    """
    Generator for safely reading large files in chunks.
    
    Args:
        filepath: Path to the file to read
        chunk_size: Size of each chunk in bytes
        
    Yields:
        Chunks of file data as bytes
    """
    try:
        with open(filepath, 'rb') as f:
            while True:
                chunk = f.read(chunk_size)
                if not chunk:
                    break
                yield chunk
    except IOError as e:
        raise NoLZSSError(f"Error reading file {filepath}: {e}")
