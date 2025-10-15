"""
Tests for parallel factorization functionality.
"""

import os
import tempfile
import unittest
from pathlib import Path

import pytest

try:
    from noLZSS.parallel import (
        parallel_factorize,
        parallel_factorize_to_file,
        parallel_factorize_file_to_file,
    )
    from noLZSS.utils import read_factors_binary_file_with_metadata
    from noLZSS import factorize
    HAS_CPP_EXTENSION = True
except ImportError:
    HAS_CPP_EXTENSION = False


@pytest.mark.skipif(not HAS_CPP_EXTENSION, reason="C++ extension not available")
class TestParallelFactorization(unittest.TestCase):
    """Test suite for parallel factorization functions."""

    def setUp(self):
        """Set up test cases."""
        self.test_string = "CGACACGTA"
        
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Create test input files
        self.input_path = self.temp_path / "input.txt"
        
        with open(self.input_path, 'w') as f:
            f.write(self.test_string)

    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()

    def test_parallel_factorize(self):
        """Test parallel factorization directly returning factors."""
        # Run with 2 threads
        factors = parallel_factorize(self.test_string, num_threads=2)
        
        # Compare with sequential factorization
        seq_factors = factorize(self.test_string)
        
        self.assertEqual(len(factors), len(seq_factors))
        
        # Verify factors cover the entire string
        positions = [f.start for f in factors]
        lengths = [f.length for f in factors]
        
        # Check if factors cover the entire string without gaps
        covered_positions = set()
        for pos, length in zip(positions, lengths):
            for i in range(pos, pos + length):
                covered_positions.add(i)
        
        self.assertEqual(len(covered_positions), len(self.test_string))

    def test_parallel_factorize_to_file(self):
        """Test parallel factorization of a string to a file."""
        output_path = self.temp_path / "output.bin"
        
        # Run with 2 threads
        num_factors = parallel_factorize_to_file(self.test_string, output_path, num_threads=2)
        
        # Verify the file exists
        self.assertTrue(output_path.exists())
        
        # Read the factors back and verify
        result = read_factors_binary_file_with_metadata(output_path)
        factors = result['factors']
        
        # Compare with sequential factorization
        seq_factors = factorize(self.test_string)
        
        self.assertEqual(len(factors), len(seq_factors))
        self.assertEqual(num_factors, len(seq_factors))

    def test_parallel_factorize_file_to_file(self):
        """Test parallel factorization from a file to another file."""
        output_path = self.temp_path / "output_from_file.bin"
        
        # Run with 2 threads
        num_factors = parallel_factorize_file_to_file(self.input_path, output_path, num_threads=2)
        
        # Verify the file exists
        self.assertTrue(output_path.exists())
        
        # Read the factors back and verify
        result = read_factors_binary_file_with_metadata(output_path)
        factors = result['factors']
        
        # Compare with sequential factorization
        seq_factors = factorize(self.test_string)
        
        self.assertEqual(len(factors), len(seq_factors))
        self.assertEqual(num_factors, len(seq_factors))

    def test_parallel_factorize_single_thread(self):
        """Test that parallel factorization works with a single thread."""
        factors_parallel = parallel_factorize(self.test_string, num_threads=1)
        factors_sequential = factorize(self.test_string)
        
        self.assertEqual(len(factors_parallel), len(factors_sequential))

    def test_parallel_factorize_auto_threads(self):
        """Test that parallel factorization works with auto-detected thread count."""
        factors = parallel_factorize(self.test_string, num_threads=0)
        seq_factors = factorize(self.test_string)
        
        self.assertEqual(len(factors), len(seq_factors))


if __name__ == '__main__':
    unittest.main()
