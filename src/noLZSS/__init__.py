"""
noLZSS: Non-overlapping Lempel-Ziv-Storer-Szymanski factorization.

A high-performance Python package with C++ core for computing non-overlapping
LZ factorizations of strings and files.
"""

# Import version with robust fallback
try:
    from ._noLZSS import __version__
except ImportError:
    # Fallback to package metadata if C++ extension is not available
    from importlib.metadata import version as _pkg_version, PackageNotFoundError
    try:
        __version__ = _pkg_version("noLZSS")
    except PackageNotFoundError:
        __version__ = "0.0.0"

from .core import *
from .utils import *
