"""
Parallel factorization functionality for noLZSS.

This module provides functions for parallel factorization of strings using
multiple threads, which can significantly improve performance on large inputs.
"""

from pathlib import Path
from typing import Union, Optional, List
import os
import tempfile

from ._noLZSS import (
    parallel_factorize_to_file as _parallel_factorize_to_file,
    parallel_factorize_file_to_file as _parallel_factorize_file_to_file,
    parallel_factorize_dna_w_rc_to_file as _parallel_factorize_dna_w_rc_to_file,
    parallel_factorize_file_dna_w_rc_to_file as _parallel_factorize_file_dna_w_rc_to_file,
    Factor,
)
from .utils import validate_input, read_factors_binary_file_with_metadata


def parallel_factorize_to_file(text: Union[str, bytes], 
                              output_path: Union[str, Path], 
                              num_threads: int = 0,
                              validate: bool = True) -> int:
    """
    Factorize text in parallel and write factors to a binary file.
    
    Args:
        text: Input text to factorize
        output_path: Path to output binary file
        num_threads: Number of threads to use (0 = auto-detect)
        validate: Whether to validate input
        
    Returns:
        Number of factors produced
        
    Raises:
        ValueError: If input validation fails
        RuntimeError: If factorization fails
        
    Example:
        >>> from noLZSS.parallel import parallel_factorize_to_file
        >>> num_factors = parallel_factorize_to_file("CGACACGTA", "output.bin", num_threads=2)
        >>> print(f"Produced {num_factors} factors")
    """
    if validate:
        text = validate_input(text)
    
    output_path = Path(output_path)
    return _parallel_factorize_to_file(text, str(output_path), num_threads)


def parallel_factorize_file_to_file(input_path: Union[str, Path], 
                                   output_path: Union[str, Path], 
                                   num_threads: int = 0) -> int:
    """
    Factorize text from file in parallel and write factors to a binary file.
    
    Args:
        input_path: Path to input text file
        output_path: Path to output binary file
        num_threads: Number of threads to use (0 = auto-detect)
        
    Returns:
        Number of factors produced
        
    Raises:
        FileNotFoundError: If input file doesn't exist
        RuntimeError: If factorization fails
        
    Example:
        >>> from noLZSS.parallel import parallel_factorize_file_to_file
        >>> num_factors = parallel_factorize_file_to_file("input.txt", "output.bin", num_threads=4)
    """
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    output_path = Path(output_path)
    return _parallel_factorize_file_to_file(str(input_path), str(output_path), num_threads)


def parallel_factorize(text: Union[str, bytes], num_threads: int = 0, validate: bool = True) -> List[Factor]:
    """
    Factorize text in parallel and return the factors.
    
    Uses a temporary file internally and cleans it up after reading the factors.
    
    Args:
        text: Input text to factorize
        num_threads: Number of threads to use (0 = auto-detect)
        validate: Whether to validate input
        
    Returns:
        List of Factor objects
        
    Raises:
        ValueError: If input validation fails
        RuntimeError: If factorization fails
        
    Example:
        >>> from noLZSS.parallel import parallel_factorize
        >>> factors = parallel_factorize("CGACACGTA", num_threads=2)
        >>> print(f"First factor: start={factors[0].start}, length={factors[0].length}")
    """
    # Create temporary file
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.bin', delete=False) as tmp:
        temp_path = Path(tmp.name)
    
    try:
        # Factorize to temporary file
        parallel_factorize_to_file(text, temp_path, num_threads, validate)
        # Read factors back using the footer format
        import struct
        from collections import namedtuple
        
        factors = []
        with open(temp_path, 'rb') as f:
            # Read footer from the end (v2 format)
            # struct FactorFileFooter {
            #     char magic[8];           // "noLZSSv2"
            #     uint64_t num_factors;
            #     uint64_t num_sequences;
            #     uint64_t num_sentinels;
            #     uint64_t footer_size;
            # }
            f.seek(-40, 2)  # Seek to 40 bytes before end
            footer_data = f.read(40)  # 8 + 4 * 8 bytes
            magic = footer_data[:8]
            if magic != b'noLZSSv2':
                raise ValueError(f"Invalid file format: expected v2 footer, got {magic}")
            
            num_factors = struct.unpack('<Q', footer_data[8:16])[0]
            
            # Seek back to beginning to read factors
            f.seek(0)
            
            # Read factors
            FactorTuple = namedtuple('Factor', ['start', 'length', 'ref'])
            for _ in range(num_factors):
                factor_data = f.read(24)  # 3 * 8 bytes (start, length, ref)
                if len(factor_data) < 24:
                    break
                start, length, ref = struct.unpack('<QQQ', factor_data)
                factors.append(FactorTuple(start, length, ref))
        
        return factors
    finally:
        # Clean up temporary file
        if temp_path.exists():
            os.unlink(temp_path)


def parallel_factorize_dna_w_rc_to_file(text: Union[str, bytes], 
                                       output_path: Union[str, Path], 
                                       num_threads: int = 0,
                                       validate: bool = True) -> int:
    """
    Factorize DNA text in parallel with reverse complement and write factors to a binary file.
    
    Args:
        text: Input DNA text
        output_path: Path to output binary file
        num_threads: Number of threads to use (0 = auto-detect)
        validate: Whether to validate input as DNA sequence
        
    Returns:
        Number of factors produced
        
    Raises:
        ValueError: If input validation fails
        RuntimeError: If factorization fails
        
    Note:
        This function is currently under development and may not be fully implemented.
        
    Example:
        >>> from noLZSS.parallel import parallel_factorize_dna_w_rc_to_file
        >>> num_factors = parallel_factorize_dna_w_rc_to_file("ATCGATCG", "output.bin", num_threads=2)
    """
    if validate:
        text = validate_input(text)
    
    output_path = Path(output_path)
    return _parallel_factorize_dna_w_rc_to_file(text, str(output_path), num_threads)


def parallel_factorize_file_dna_w_rc_to_file(input_path: Union[str, Path], 
                                            output_path: Union[str, Path], 
                                            num_threads: int = 0) -> int:
    """
    Factorize DNA text from file in parallel with reverse complement and write factors to a binary file.
    
    Args:
        input_path: Path to input DNA text file
        output_path: Path to output binary file
        num_threads: Number of threads to use (0 = auto-detection)
        
    Returns:
        Number of factors produced
        
    Raises:
        FileNotFoundError: If input file doesn't exist
        RuntimeError: If factorization fails
        
    Note:
        This function is currently under development and may not be fully implemented.
        
    Example:
        >>> from noLZSS.parallel import parallel_factorize_file_dna_w_rc_to_file
        >>> num_factors = parallel_factorize_file_dna_w_rc_to_file("dna.txt", "output.bin", num_threads=4)
    """
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    output_path = Path(output_path)
    return _parallel_factorize_file_dna_w_rc_to_file(str(input_path), str(output_path), num_threads)
