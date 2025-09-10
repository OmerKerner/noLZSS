"""
noLZSS: Non-overlapping Lempel-Ziv-Storer-Szymanski factorization.

A high-performance Python package with C++ core for computing non-overlapping
LZ factorizations of strings and files.
"""

# Import C++ bindings
from ._noLZSS import (
    __version__
)

from .core import *
from .utils import *
