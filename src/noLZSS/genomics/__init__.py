"""
Genomics-specific functionality for noLZSS.

This subpackage provides specialized tools for working with biological sequences,
including FASTA file parsing, sequence validation, and genomics-aware compression.
"""

from .fasta import *
from .sequences import *

__all__ = [
    # From fasta module
    "read_nucleotide_fasta",
    "read_protein_fasta", 
    "read_fasta_auto",
    "process_fasta_with_plots",
    "FASTAError",
    
    # From sequences module
    "is_dna_sequence",
    "is_protein_sequence", 
    "detect_sequence_type"
]
