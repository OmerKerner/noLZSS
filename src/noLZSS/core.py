"""
Core Python wrappers for noLZSS C++ functionality.

This module provides enhanced Python wrappers around the C++ factorization functions,
adding input validation, error handling, and convenience features.
"""

from typing import List, Tuple, Union, Optional
import os
from pathlib import Path

from ._noLZSS import (
    factorize as _factorize,
    factorize_file as _factorize_file,
    count_factors as _count_factors,
    count_factors_file as _count_factors_file,
    write_factors_binary_file as _write_factors_binary_file,
    factorize_w_reference as _factorize_w_reference,
    factorize_w_reference_file as _factorize_w_reference_file,
)

from .utils import validate_input, analyze_alphabet


def factorize(data: Union[str, bytes], validate: bool = True) -> List[Tuple[int, int, int]]:
    """
    Factorize a string or bytes object into LZ factors.
    
    Args:
        data: Input string or bytes to factorize
        validate: Whether to perform input validation (default: True)
    
    Returns:
        List of (position, length, ref) tuples representing the factorization
        
    Raises:
        ValueError: If input is invalid (empty, etc.)
        TypeError: If input type is not supported
    """
    if validate:
        data = validate_input(data)
    
    return _factorize(data)


def factorize_file(filepath: Union[str, Path], reserve_hint: int = 0) -> List[Tuple[int, int, int]]:
    """
    Factorize the contents of a file into LZ factors.
    
    Args:
        filepath: Path to the input file
        reserve_hint: Optional hint for reserving space in output vector (0 = no hint)
        
    Returns:
        List of (position, length, ref) tuples representing the factorization
        
    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    return _factorize_file(str(filepath), reserve_hint)


def count_factors(data: Union[str, bytes], validate: bool = True) -> int:
    """
    Count the number of factors in a string without computing the full factorization.
    
    Args:
        data: Input string or bytes to analyze
        validate: Whether to perform input validation (default: True)
        
    Returns:
        Number of factors in the factorization
        
    Raises:
        ValueError: If input is invalid
        TypeError: If input type is not supported
    """
    if validate:
        data = validate_input(data)
    
    return _count_factors(data)


def count_factors_file(filepath: Union[str, Path], validate: bool = True) -> int:
    """
    Count the number of factors in a file without computing the full factorization.
    
    Args:
        filepath: Path to the input file
        validate: Whether to perform input validation (default: True)
        
    Returns:
        Number of factors in the factorization
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If file contents are invalid
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    return _count_factors_file(str(filepath))


def write_factors_binary_file(
    data: Union[str, bytes], 
    output_filepath: Union[str, Path]
) -> None:
    """
    Factorize input and write the factors to a binary file.
    
    Args:
        data: Input string or bytes to factorize
        output_filepath: Path where to write the binary factors
        
    Raises:
        ValueError: If input is invalid
        TypeError: If input type is not supported
        OSError: If unable to write to output file
    """
    data = validate_input(data)
    
    output_filepath = Path(output_filepath)
    # Ensure output directory exists
    output_filepath.parent.mkdir(parents=True, exist_ok=True)
    
    _write_factors_binary_file(data, str(output_filepath))


def factorize_with_info(data: Union[str, bytes], validate: bool = True) -> dict:
    """
    Factorize input and return both factors and additional information.
    
    Args:
        data: Input string or bytes to factorize
        validate: Whether to perform input validation (default: True)
        
    Returns:
        Dictionary containing:
        - 'factors': List of (position, length, ref) tuples
        - 'alphabet_info': Alphabet analysis results
        - 'input_size': Size of input data
        - 'num_factors': Number of factors
    """
    if validate:
        data = validate_input(data)
    
    factors = _factorize(data)
    alphabet_info = analyze_alphabet(data)
    
    return {
        'factors': factors,
        'alphabet_info': alphabet_info,
        'input_size': len(data),
        'num_factors': len(factors)
    }


def factorize_w_reference(reference_seq: Union[str, bytes], target_seq: Union[str, bytes], validate: bool = True) -> List[Tuple[int, int, int]]:
    """
    Factorize target sequence using a reference sequence without reverse complement.
    
    Concatenates a reference sequence and target sequence, then performs noLZSS factorization
    starting from where the target sequence begins. This allows the target sequence to reference
    patterns in the reference sequence without factorizing the reference itself. Suitable for
    general text or amino acid sequences.
    
    Args:
        reference_seq: Reference sequence (any text)
        target_seq: Target sequence to be factorized (any text)
        validate: Whether to perform input validation (default: True)
        
    Returns:
        List of (start, length, ref) tuples representing the factorization of target sequence
        
    Raises:
        ValueError: If sequences are empty
        TypeError: If input types are not supported
        RuntimeError: If processing errors occur
        
    Note:
        Factor start positions are absolute positions in the combined reference+target string.
        No reverse complement matching is performed - suitable for text or amino acid sequences.
        
    Warning:
        The sentinel character '\\x01' (ASCII 1) must not appear in either input sequence,
        as it is used internally to separate the reference and target sequences.
    """
    from .utils import validate_input
    
    if validate:
        reference_seq = validate_input(reference_seq)
        target_seq = validate_input(target_seq)
    
    # Convert to strings if they're bytes
    if isinstance(reference_seq, bytes):
        reference_seq = reference_seq.decode('ascii')
    if isinstance(target_seq, bytes):
        target_seq = target_seq.decode('ascii')
    
    return _factorize_w_reference(reference_seq, target_seq)


def factorize_w_reference_file(reference_seq: Union[str, bytes], target_seq: Union[str, bytes], 
                               output_path: Union[str, Path], validate: bool = True) -> int:
    """
    Factorize target sequence using a reference sequence and write factors to binary file.
    
    Concatenates a reference sequence and target sequence, then performs noLZSS factorization
    starting from where the target sequence begins, and writes the resulting factors to a binary file.
    Suitable for general text or amino acid sequences.
    
    Args:
        reference_seq: Reference sequence (any text)
        target_seq: Target sequence to be factorized (any text)
        output_path: Path to output file where binary factors will be written
        validate: Whether to perform input validation (default: True)
        
    Returns:
        Number of factors written to the output file
        
    Raises:
        ValueError: If sequences are empty
        TypeError: If input types are not supported
        RuntimeError: If unable to create output file or processing errors occur
        
    Note:
        Factor start positions are absolute positions in the combined reference+target string.
        No reverse complement matching is performed - suitable for text or amino acid sequences.
        This function overwrites the output file if it exists.
        
    Warning:
        The sentinel character '\\x01' (ASCII 1) must not appear in either input sequence,
        as it is used internally to separate the reference and target sequences.
    """
    from .utils import validate_input
    
    if validate:
        reference_seq = validate_input(reference_seq)
        target_seq = validate_input(target_seq)
    
    # Convert to strings if they're bytes
    if isinstance(reference_seq, bytes):
        reference_seq = reference_seq.decode('ascii')
    if isinstance(target_seq, bytes):
        target_seq = target_seq.decode('ascii')
    
    # Ensure output directory exists
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    return _factorize_w_reference_file(reference_seq, target_seq, str(output_path))
