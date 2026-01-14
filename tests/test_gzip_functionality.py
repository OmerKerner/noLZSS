"""
Tests for gzip compression functionality in binary output files.

This module tests that binary files are written as gzipped and can be read
back correctly, supporting both gzipped and non-gzipped formats.
"""

import tempfile
import os
import gzip
import struct
from pathlib import Path

# Try importing noLZSS - gracefully handle if C++ extension not available
try:
    import noLZSS
    from noLZSS.utils import (
        read_factors_binary_file,
        read_binary_file_metadata,
        read_factors_binary_file_with_metadata,
        _is_gzipped_file
    )
    NOLZSS_AVAILABLE = True
except ImportError:
    NOLZSS_AVAILABLE = False
    import pytest
    pytestmark = pytest.mark.skip("noLZSS C++ extension not available")


class TestGzipFunctionality:
    """Test suite for gzip compression of binary output files."""
    
    def test_basic_write_creates_gzipped_file(self):
        """Test that write_factors_binary_file creates a gzipped file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = Path(tmpdir) / "input.txt"
            output_file = Path(tmpdir) / "output.bin.gz"
            
            # Create a simple input file
            input_file.write_text("ABRACADABRA")
            
            # Factorize and write to binary file
            num_factors = noLZSS.write_factors_binary_file(str(input_file), str(output_file))
            
            # Check that output file exists and is gzipped
            assert output_file.exists(), "Output file should exist"
            assert _is_gzipped_file(output_file), "Output file should be gzipped"
            assert num_factors > 0, "Should have at least one factor"
    
    def test_read_gzipped_file(self):
        """Test that we can read back factors from a gzipped file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = Path(tmpdir) / "input.txt"
            output_file = Path(tmpdir) / "output.bin.gz"
            
            # Create input and factorize
            test_text = "ABRACADABRA"
            input_file.write_text(test_text)
            
            num_factors = noLZSS.write_factors_binary_file(str(input_file), str(output_file))
            
            # Read back the factors
            factors = read_factors_binary_file(output_file)
            
            assert len(factors) == num_factors, "Should read same number of factors as written"
            assert all(len(f) == 3 for f in factors), "Each factor should be a 3-tuple"
            
            # Verify factors cover the input text
            total_length = sum(f[1] for f in factors)  # Sum of lengths
            assert total_length == len(test_text), "Factors should cover entire input text"
    
    def test_read_metadata_from_gzipped_file(self):
        """Test that we can read metadata from a gzipped file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = Path(tmpdir) / "input.txt"
            output_file = Path(tmpdir) / "output.bin.gz"
            
            test_text = "ABRACADABRA"
            input_file.write_text(test_text)
            
            num_factors = noLZSS.write_factors_binary_file(str(input_file), str(output_file))
            
            # Read metadata
            metadata = read_binary_file_metadata(output_file)
            
            assert metadata['num_factors'] == num_factors
            assert metadata['total_length'] == len(test_text)
    
    def test_gzipped_file_is_smaller(self):
        """Test that gzipped files are actually compressed (smaller than uncompressed)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = Path(tmpdir) / "input.txt"
            output_file = Path(tmpdir) / "output.bin.gz"
            
            # Create a longer input text with repetition (should compress well)
            test_text = "ABRACADABRA" * 100
            input_file.write_text(test_text)
            
            num_factors = noLZSS.write_factors_binary_file(str(input_file), str(output_file))
            
            # Get file sizes
            gzipped_size = output_file.stat().st_size
            
            # Decompress to check uncompressed size
            with gzip.open(output_file, 'rb') as f:
                uncompressed_data = f.read()
            
            uncompressed_size = len(uncompressed_data)
            
            # Gzipped should be smaller (at least 10% savings for this test)
            compression_ratio = gzipped_size / uncompressed_size
            print(f"Compression ratio: {compression_ratio:.2%}")
            assert compression_ratio < 0.9, "Gzipped file should be significantly smaller"
    
    def test_dna_w_rc_creates_gzipped_file(self):
        """Test that DNA with reverse complement write function creates gzipped files."""
        import noLZSS._noLZSS as cpp
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = Path(tmpdir) / "dna.txt"
            output_file = Path(tmpdir) / "dna_output.bin.gz"
            
            # Create a DNA sequence
            dna_seq = "ATCGATCGATCG"
            input_file.write_text(dna_seq)
            
            # Factorize with reverse complement awareness using C++ API directly
            num_factors = cpp.write_factors_binary_file_dna_w_rc(str(input_file), str(output_file))
            
            # Verify it's gzipped
            assert output_file.exists()
            assert _is_gzipped_file(output_file)
            assert num_factors > 0
            
            # Read back and verify
            factors = read_factors_binary_file(output_file)
            assert len(factors) == num_factors
    
    def test_file_without_gz_extension_still_gzipped(self):
        """Test that files are gzipped even without .gz extension."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = Path(tmpdir) / "input.txt"
            output_file = Path(tmpdir) / "output.bin"  # No .gz extension
            
            input_file.write_text("TESTDATA")
            
            num_factors = noLZSS.write_factors_binary_file(str(input_file), str(output_file))
            
            # Should still be gzipped
            assert output_file.exists()
            assert _is_gzipped_file(output_file), "File should be gzipped even without .gz extension"
            
            # Should still be readable
            factors = read_factors_binary_file(output_file)
            assert len(factors) == num_factors


if __name__ == "__main__":
    if NOLZSS_AVAILABLE:
        import sys
        test = TestGzipFunctionality()
        
        print("Running gzip functionality tests...")
        tests = [
            ("test_basic_write_creates_gzipped_file", test.test_basic_write_creates_gzipped_file),
            ("test_read_gzipped_file", test.test_read_gzipped_file),
            ("test_read_metadata_from_gzipped_file", test.test_read_metadata_from_gzipped_file),
            ("test_gzipped_file_is_smaller", test.test_gzipped_file_is_smaller),
            ("test_dna_w_rc_creates_gzipped_file", test.test_dna_w_rc_creates_gzipped_file),
            ("test_file_without_gz_extension_still_gzipped", test.test_file_without_gz_extension_still_gzipped),
        ]
        
        passed = 0
        failed = 0
        
        for test_name, test_func in tests:
            try:
                print(f"\n{test_name}...", end=" ")
                test_func()
                print("PASSED")
                passed += 1
            except Exception as e:
                print(f"FAILED: {e}")
                failed += 1
                import traceback
                traceback.print_exc()
        
        print(f"\n{'='*60}")
        print(f"Results: {passed} passed, {failed} failed")
        print(f"{'='*60}")
        
        sys.exit(0 if failed == 0 else 1)
    else:
        print("noLZSS not available - skipping tests")
        sys.exit(0)
