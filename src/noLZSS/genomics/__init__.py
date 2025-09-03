"""
Genomics-specific functionality for noLZSS.

This subpackage provides specialized tools for working with biological sequences,
including FASTA file parsing, sequence validation, and genomics-aware compression.
"""

from .fasta import *
from .sequences import *

__all__ = [
    # From fasta module (will be populated when implemented)
    
    # From sequences module
    "is_dna_sequence",
    "is_protein_sequence", 
    "detect_sequence_type"
]
