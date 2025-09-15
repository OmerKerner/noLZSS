"""
Genomics-specific functionality for noLZSS.

This subpackage provides specialized tools for working with biological sequences,
including FASTA file parsing, sequence validation, and genomics-aware compression.
"""

from .._noLZSS import (
    # Import all genomics-related C++ bindings
    factorize_dna_w_rc,
    factorize_file_dna_w_rc,
    count_factors_dna_w_rc,
    count_factors_file_dna_w_rc,
    write_factors_binary_file_dna_w_rc,
    factorize_multiple_dna_w_rc,
    factorize_file_multiple_dna_w_rc,
    count_factors_multiple_dna_w_rc,
    count_factors_file_multiple_dna_w_rc,
    write_factors_binary_file_multiple_dna_w_rc,
    process_nucleotide_fasta,
    process_amino_acid_fasta,
    factorize_fasta_multiple_dna_w_rc,
    prepare_multiple_dna_sequences_w_rc,
)

from .fasta import *
from .sequences import *
from .plots import *
from . import batch_factorize
