"""
Sequence utilities for biological data.

This module provides functions for working with nucleotide and amino acid sequences,
including validation, transformation, and analysis functions.
"""

from typing import Union
import re


def is_dna_sequence(data: Union[str, bytes]) -> bool:
    """
    Check if data appears to be a DNA sequence (A, T, G, C).
    
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
        
    # Allow standard DNA bases only, case insensitive
    dna_pattern = re.compile(r'^[ATGC]+$', re.IGNORECASE)
    return bool(dna_pattern.match(data))


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
    protein_pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWYBJOUXZ]+$', re.IGNORECASE)
    return bool(protein_pattern.match(data))


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
    
    if not isinstance(data, str):
        return 'binary'
    
    # Convert to uppercase for analysis
    data_upper = data.upper()
    
    # Handle empty string
    if not data_upper:
        return 'text'
    
    # Check if all characters are alphabetic (no numbers, punctuation, etc.)
    if not all(c.isalpha() for c in data_upper):
        return 'text'
    
    # Check for characters that are amino-acid specific (not DNA nucleotides)
    amino_acid_only_chars = set('EFHIKLMNPQRSVWY')  # These are not DNA nucleotides
    has_amino_specific = any(c in amino_acid_only_chars for c in data_upper)
    
    # Check if all characters are valid DNA
    dna_chars = set('ACGT')
    all_dna_chars = all(c in dna_chars for c in data_upper)
    
    # Check if all characters are valid amino acids
    amino_chars = set('ACDEFGHIKLMNPQRSTVWY')
    all_amino_chars = all(c in amino_chars for c in data_upper)
    
    # Decision logic:
    # If it has amino-acid-specific characters, it's protein
    # If it only contains DNA characters and no amino-specific chars, it's DNA
    # If it contains other characters, check broader amino acid set
    if has_amino_specific and all_amino_chars:
        return 'protein'
    elif all_dna_chars and not has_amino_specific:
        return 'dna'
    elif all_amino_chars:
        return 'protein'
    else:
        return 'text'


def factorize_w_reference_seq(reference_seq: Union[str, bytes], target_seq: Union[str, bytes], validate: bool = True):
    """
    Factorize target DNA sequence using a reference sequence with reverse complement awareness.
    
    Concatenates a reference sequence and target sequence, then performs noLZSS factorization
    with reverse complement awareness starting from where the target sequence begins. This allows
    the target sequence to reference patterns in the reference sequence without factorizing the
    reference itself.
    
    Args:
        reference_seq: Reference DNA sequence (A, C, T, G - case insensitive)
        target_seq: Target DNA sequence to be factorized (A, C, T, G - case insensitive)
        validate: Whether to perform input validation (default: True)
        
    Returns:
        List of (start, length, ref, is_rc) tuples representing the factorization of target sequence
        
    Raises:
        ValueError: If sequences contain invalid nucleotides or are empty
        TypeError: If input types are not supported
        RuntimeError: If processing errors occur
        
    Note:
        Factor start positions are relative to the beginning of the target sequence.
        Both sequences are converted to uppercase before factorization.
        ref field has RC_MASK cleared. is_rc boolean indicates reverse complement matches.
    """
    from .._noLZSS import factorize_w_reference_seq as _factorize_w_reference_seq
    from ..utils import validate_input
    
    if validate:
        reference_seq = validate_input(reference_seq)
        target_seq = validate_input(target_seq)
        
        # Additional validation for DNA sequences
        if not is_dna_sequence(reference_seq):
            raise ValueError("Reference sequence must contain only DNA nucleotides (A, C, T, G)")
        if not is_dna_sequence(target_seq):
            raise ValueError("Target sequence must contain only DNA nucleotides (A, C, T, G)")
    
    # Convert to strings if they're bytes
    if isinstance(reference_seq, bytes):
        reference_seq = reference_seq.decode('ascii')
    if isinstance(target_seq, bytes):
        target_seq = target_seq.decode('ascii')
    
    return _factorize_w_reference_seq(reference_seq, target_seq)


def factorize_w_reference_seq_file(reference_seq: Union[str, bytes], target_seq: Union[str, bytes], 
                                   output_path: Union[str, 'Path'], validate: bool = True) -> int:
    """
    Factorize target DNA sequence using a reference sequence and write factors to binary file.
    
    Concatenates a reference sequence and target sequence, then performs noLZSS factorization
    with reverse complement awareness starting from where the target sequence begins, and writes
    the resulting factors to a binary file.
    
    Args:
        reference_seq: Reference DNA sequence (A, C, T, G - case insensitive)
        target_seq: Target DNA sequence to be factorized (A, C, T, G - case insensitive)
        output_path: Path to output file where binary factors will be written
        validate: Whether to perform input validation (default: True)
        
    Returns:
        Number of factors written to the output file
        
    Raises:
        ValueError: If sequences contain invalid nucleotides or are empty
        TypeError: If input types are not supported
        RuntimeError: If unable to create output file or processing errors occur
        
    Note:
        Factor start positions are relative to the beginning of the target sequence.
        Binary format follows the same structure as other DNA factorization binary outputs.
        This function overwrites the output file if it exists.
    """
    from .._noLZSS import factorize_w_reference_seq_file as _factorize_w_reference_seq_file
    from ..utils import validate_input
    from pathlib import Path
    
    if validate:
        reference_seq = validate_input(reference_seq)
        target_seq = validate_input(target_seq)
        
        # Additional validation for DNA sequences
        if not is_dna_sequence(reference_seq):
            raise ValueError("Reference sequence must contain only DNA nucleotides (A, C, T, G)")
        if not is_dna_sequence(target_seq):
            raise ValueError("Target sequence must contain only DNA nucleotides (A, C, T, G)")
    
    # Convert to strings if they're bytes
    if isinstance(reference_seq, bytes):
        reference_seq = reference_seq.decode('ascii')
    if isinstance(target_seq, bytes):
        target_seq = target_seq.decode('ascii')
    
    # Ensure output directory exists
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    return _factorize_w_reference_seq_file(reference_seq, target_seq, str(output_path))
