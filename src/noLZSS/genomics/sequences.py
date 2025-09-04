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
