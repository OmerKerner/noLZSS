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
    
    # Check DNA first (more restrictive)
    if is_dna_sequence(data):
        return 'dna'
    elif is_protein_sequence(data):
        return 'protein'
    elif isinstance(data, str):
        return 'text'
    else:
        return 'binary'
