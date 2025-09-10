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
