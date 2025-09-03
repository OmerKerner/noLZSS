"""
noLZSS: Non-overlapping Lempel-Ziv-Storer-Szymanski factorization.

A high-performance Python package with C++ core for computing non-overlapping
LZ factorizations of strings and files.
"""

# Import C++ bindings
from ._noLZSS import (
    factorize as _factorize,
    factorize_file as _factorize_file,
    count_factors as _count_factors,
    count_factors_file as _count_factors_file,
    write_factors_binary_file as _write_factors_binary_file,
    __version__
)

# Import enhanced Python wrappers
from .core import (
    factorize,
    factorize_file, 
    count_factors,
    count_factors_file,
    write_factors_binary_file,
    factorize_with_info
)

# Import utilities
from .utils import (
    analyze_alphabet,
    read_factors_binary_file,
    plot_factor_lengths
)
from .genomics import (
    detect_sequence_type,
    is_dna_sequence,
    is_protein_sequence
)

__all__ = [
    # Enhanced wrapper functions (recommended for most users)
    "factorize",
    "factorize_file", 
    "count_factors",
    "count_factors_file",
    "write_factors_binary_file",
    "factorize_with_info",
    
    # Utility functions
    "analyze_alphabet",
    "read_factors_binary_file",
    "plot_factor_lengths",
    "detect_sequence_type", 
    "is_dna_sequence",
    "is_protein_sequence",
    
    # Version info
    "__version__"
]
