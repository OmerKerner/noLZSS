"""
Genomics-specific functionality for noLZSS.

This subpackage provides specialized tools for working with biological sequences,
including FASTA file parsing, sequence validation, and genomics-aware compression.
"""

from .fasta import *
from .sequences import *
from .plots import *

__all__ = [
    # From fasta module
    "read_nucleotide_fasta",
    "read_protein_fasta", 
    "read_fasta_auto",
    "FASTAError",
    
    # From plots module
    "plot_single_seq_accum_factors_from_fasta",
    "plot_multiple_seq_self_weizmann_factor_plot_from_fasta",
    "PlotError",
    
    # From sequences module
    "is_dna_sequence",
    "is_protein_sequence", 
    "detect_sequence_type"
]
