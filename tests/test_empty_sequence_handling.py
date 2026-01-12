"""
Test empty sequence handling in factorization to ensure SDSL warnings are eliminated.

This test specifically validates the fixes for empty sequence handling that prevent
SDSL int_vector warnings and ensure proper error propagation from worker threads.
"""

import pytest
import tempfile
import os
from pathlib import Path

try:
    import noLZSS
    import noLZSS._noLZSS as cpp
    factorize_fasta_multiple_dna_w_rc = cpp.factorize_fasta_multiple_dna_w_rc
    factorize_fasta_multiple_dna_no_rc = cpp.factorize_fasta_multiple_dna_no_rc
    cpp_available = True
except ImportError:
    cpp_available = False

pytestmark = pytest.mark.skipif(not cpp_available, reason="C++ extension not available")


class TestEmptySequenceHandling:
    """Test cases for empty sequence handling"""
    
    def test_empty_fasta_file_with_rc(self):
        """Test that completely empty FASTA files are handled gracefully with RC"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write("")  # Completely empty file
            fasta_path = f.name
        
        try:
            # Should raise RuntimeError with "No valid sequences found"
            with pytest.raises(RuntimeError, match="No valid sequences found"):
                factorize_fasta_multiple_dna_w_rc(fasta_path)
        finally:
            os.unlink(fasta_path)
    
    def test_empty_fasta_file_without_rc(self):
        """Test that completely empty FASTA files are handled gracefully without RC"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write("")  # Completely empty file
            fasta_path = f.name
        
        try:
            # Should raise RuntimeError with "No valid sequences found"
            with pytest.raises(RuntimeError, match="No valid sequences found"):
                factorize_fasta_multiple_dna_no_rc(fasta_path)
        finally:
            os.unlink(fasta_path)
    
    def test_fasta_with_only_headers(self):
        """Test FASTA file with headers but no sequence data"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")  # Header with no sequence
            f.write(">seq2\n")  # Another header with no sequence
            fasta_path = f.name
        
        try:
            # Should raise RuntimeError because no sequences have data
            with pytest.raises(RuntimeError, match="No valid sequences found"):
                factorize_fasta_multiple_dna_w_rc(fasta_path)
        finally:
            os.unlink(fasta_path)
    
    def test_fasta_with_mixed_empty_sequences(self):
        """Test FASTA file where some sequences are empty and some are not"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")  # Empty sequence
            f.write(">seq2\n")
            f.write("ATCGATCG\n")  # Valid sequence
            f.write(">seq3\n")  # Another empty sequence
            fasta_path = f.name
        
        try:
            # Should succeed - empty sequences are skipped by the parser
            result = factorize_fasta_multiple_dna_w_rc(fasta_path)
            # Result is a tuple: (factors, sentinel_indices, sequence_ids)
            factors, sentinel_indices, sequence_ids = result
            assert len(sequence_ids) == 1  # Only seq2 should be counted
            assert len(factors) > 0
        finally:
            os.unlink(fasta_path)
    
    def test_parallel_factorization_empty_string(self):
        """Test parallel factorization with empty string"""
        # Test with DNA w/ RC parallel function - use the file-based version
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bin', delete=False) as f:
            output_path = f.name
        
        try:
            # Should return 0 factors without throwing or warning
            result = cpp.parallel_factorize_dna_w_rc_to_file("", output_path, 2)
            assert result == 0
            
            # Output file should exist but be minimal/empty
            assert os.path.exists(output_path)
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)
    
    def test_prepare_empty_sequences_list(self):
        """Test that prepare_multiple_dna_sequences_w_rc handles empty list"""
        # Empty list should return empty prepared result
        # Result is a tuple: (prepared_string, original_length, sentinel_positions)
        result = cpp.prepare_multiple_dna_sequences_w_rc([])
        prepared_string, original_length, sentinel_positions = result
        assert prepared_string == ""
        assert original_length == 0
        assert sentinel_positions == []
    
    def test_factorize_with_very_short_sequence(self):
        """Test factorization with single character sequences"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")
            f.write("A\n")
            fasta_path = f.name
        
        try:
            # Should handle single character sequences without issues
            result = factorize_fasta_multiple_dna_w_rc(fasta_path)
            # Result is a tuple: (factors, sentinel_indices, sequence_ids)
            factors, sentinel_indices, sequence_ids = result
            assert len(sequence_ids) == 1
            assert len(factors) > 0  # Should have at least one factor
        finally:
            os.unlink(fasta_path)


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v"])
