"""
Additional tests for empty sequence handling fixes.

This test module validates the specific fixes implemented to address
worker process crashes in batch factorization.
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
    prepare_multiple_dna_sequences_w_rc = cpp.prepare_multiple_dna_sequences_w_rc
    prepare_multiple_dna_sequences_no_rc = cpp.prepare_multiple_dna_sequences_no_rc
    cpp_available = True
except ImportError:
    cpp_available = False

pytestmark = pytest.mark.skipif(not cpp_available, reason="C++ extension not available")


class TestEmptySequenceFixes:
    """Test cases for the empty sequence handling fixes"""
    
    def test_fasta_parser_filters_empty_with_warning(self):
        """Test that FASTA parser filters empty sequences and warns about them"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            # Create FASTA with some empty and some valid sequences
            f.write(">seq1\n")  # Empty sequence (will be filtered)
            f.write(">seq2\n")
            f.write("ATCGATCG\n")  # Valid sequence
            f.write(">seq3\n")  # Empty sequence (will be filtered)
            f.write(">seq4\n")
            f.write("GCTAGCTA\n")  # Valid sequence
            fasta_path = f.name
        
        try:
            # Should succeed with only the non-empty sequences
            result = factorize_fasta_multiple_dna_w_rc(fasta_path)
            factors, sentinel_indices, sequence_ids = result
            
            # Only seq2 and seq4 should be in the result
            assert len(sequence_ids) == 2
            assert sequence_ids[0] == "seq2"
            assert sequence_ids[1] == "seq4"
            assert len(factors) > 0
        finally:
            os.unlink(fasta_path)
    
    def test_prepare_sequences_filters_empty_strings_w_rc(self):
        """Test that prepare_multiple_dna_sequences_w_rc filters empty strings"""
        # Mix of empty and valid sequences
        sequences = ["", "ATCG", "", "GCTA", ""]
        
        # Should filter out empty strings and succeed
        result = prepare_multiple_dna_sequences_w_rc(sequences)
        prepared_string, original_length, sentinel_positions = result
        
        # Should only contain the 2 non-empty sequences and their reverse complements
        assert len(sentinel_positions) == 4  # 2 original + 2 RC sentinels
        assert "ATCG" in prepared_string
        assert "GCTA" in prepared_string
    
    def test_prepare_sequences_filters_empty_strings_no_rc(self):
        """Test that prepare_multiple_dna_sequences_no_rc filters empty strings"""
        # Mix of empty and valid sequences
        sequences = ["", "ATCG", "", "GCTA", ""]
        
        # Should filter out empty strings and succeed
        result = prepare_multiple_dna_sequences_no_rc(sequences)
        prepared_string, original_length, sentinel_positions = result
        
        # Should only contain the 2 non-empty sequences
        assert len(sentinel_positions) == 1  # 1 sentinel between sequences
        assert "ATCG" in prepared_string
        assert "GCTA" in prepared_string
    
    def test_prepare_sequences_all_empty_raises_error_w_rc(self):
        """Test that prepare_multiple_dna_sequences_w_rc raises error when all sequences are empty"""
        sequences = ["", "", ""]
        
        with pytest.raises(RuntimeError, match="All sequences are empty"):
            prepare_multiple_dna_sequences_w_rc(sequences)
    
    def test_prepare_sequences_all_empty_raises_error_no_rc(self):
        """Test that prepare_multiple_dna_sequences_no_rc raises error when all sequences are empty"""
        sequences = ["", "", ""]
        
        with pytest.raises(RuntimeError, match="All sequences are empty"):
            prepare_multiple_dna_sequences_no_rc(sequences)
    
    def test_very_short_sequence_single_nucleotide(self):
        """Test that single nucleotide sequences are handled correctly"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")
            f.write("A\n")
            f.write(">seq2\n")
            f.write("T\n")
            fasta_path = f.name
        
        try:
            # Should handle single nucleotide sequences
            result = factorize_fasta_multiple_dna_w_rc(fasta_path)
            factors, sentinel_indices, sequence_ids = result
            assert len(sequence_ids) == 2
            assert len(factors) > 0
        finally:
            os.unlink(fasta_path)
    
    def test_two_nucleotide_sequences(self):
        """Test that 2-nucleotide sequences are handled correctly"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")
            f.write("AT\n")
            f.write(">seq2\n")
            f.write("GC\n")
            fasta_path = f.name
        
        try:
            # Should handle 2-nucleotide sequences
            result = factorize_fasta_multiple_dna_w_rc(fasta_path)
            factors, sentinel_indices, sequence_ids = result
            assert len(sequence_ids) == 2
            assert len(factors) > 0
        finally:
            os.unlink(fasta_path)
    
    def test_fasta_only_empty_sequences_raises_error(self):
        """Test that FASTA files with only empty sequences raise clear error"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")  # Empty sequence
            f.write(">seq2\n")  # Empty sequence
            f.write(">seq3\n")  # Empty sequence
            fasta_path = f.name
        
        try:
            # Should raise RuntimeError with clear message
            with pytest.raises(RuntimeError, match="No valid sequences found"):
                factorize_fasta_multiple_dna_w_rc(fasta_path)
        finally:
            os.unlink(fasta_path)
    
    def test_mixed_empty_and_short_sequences(self):
        """Test handling of mix of empty and very short sequences"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">empty1\n")  # Empty
            f.write(">short1\n")
            f.write("A\n")  # Single nucleotide
            f.write(">empty2\n")  # Empty
            f.write(">short2\n")
            f.write("GC\n")  # Two nucleotides
            f.write(">normal\n")
            f.write("ATCGATCG\n")  # Normal length
            fasta_path = f.name
        
        try:
            # Should filter empty and handle short sequences
            result = factorize_fasta_multiple_dna_w_rc(fasta_path)
            factors, sentinel_indices, sequence_ids = result
            
            # Only non-empty sequences should be present
            assert len(sequence_ids) == 3
            assert "short1" in sequence_ids
            assert "short2" in sequence_ids
            assert "normal" in sequence_ids
            assert len(factors) > 0
        finally:
            os.unlink(fasta_path)


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v"])
