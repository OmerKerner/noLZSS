"""
Comprehensive tests for the utils module.
"""
import warnings
import pytest
import sys
import os
import tempfile
import struct
from pathlib import Path

from noLZSS.utils import (
    validate_input, analyze_alphabet, 
    read_factors_binary_file, plot_factor_lengths,
    InvalidInputError, NoLZSSError
)
from noLZSS.genomics import (
    is_dna_sequence, is_protein_sequence, detect_sequence_type
)


class TestValidateInput:
    """Test input validation functions."""
    
    def test_validate_string_input(self):
        """Test validation of string inputs."""
        result = validate_input("hello")
        assert isinstance(result, bytes)
        assert result == b"hello"
    
    def test_validate_bytes_input(self):
        """Test validation of bytes inputs."""
        result = validate_input(b"hello")
        assert isinstance(result, bytes)
        assert result == b"hello"
    
    def test_validate_unicode_string(self):
        """Test validation of unicode strings (should reject non-ASCII)."""
        with pytest.raises(InvalidInputError, match="Input string must contain only ASCII characters"):
            validate_input("héllo")
    
    def test_validate_empty_input_raises_error(self):
        """Test that empty input raises error."""
        with pytest.raises(InvalidInputError):
            validate_input("")
        
        with pytest.raises(InvalidInputError):
            validate_input(b"")
    
    def test_validate_null_bytes_raises_error(self):
        """Test that null bytes in middle raise error."""
        with pytest.raises(InvalidInputError):
            validate_input(b"hello\x00world")
    
    def test_validate_null_byte_at_end_allowed(self):
        """Test that null byte at end is allowed."""
        result = validate_input(b"hello\x00")
        assert result == b"hello\x00"
    
    def test_validate_invalid_type_raises_error(self):
        """Test that invalid input types raise TypeError."""
        with pytest.raises(TypeError):
            validate_input(123)
        
        with pytest.raises(TypeError):
            validate_input(None)
    
    def test_validate_invalid_unicode_raises_error(self):
        """Test that invalid unicode raises error."""
        # This is tricky to test directly, but we can create a mock scenario
        pass  # Most strings are valid UTF-8, so this is hard to trigger


class TestAnalyzeAlphabet:
    """Test alphabet analysis functions."""
    
    def test_analyze_string_alphabet(self):
        """Test alphabet analysis on strings."""
        result = analyze_alphabet("ATCGATCG")
        
        assert result['size'] == 4  # A, T, C, G
        assert result['total_length'] == 8
        assert result['characters'] == {'A', 'T', 'C', 'G'}
        assert 'distribution' in result
        assert 'entropy' in result
        assert 'most_common' in result
        
        # Check distribution
        assert result['distribution']['A'] == 2
        assert result['distribution']['T'] == 2
        assert result['distribution']['C'] == 2
        assert result['distribution']['G'] == 2
    
    def test_analyze_bytes_alphabet(self):
        """Test alphabet analysis on bytes."""
        result = analyze_alphabet(b"ATCGATCG")
        
        assert result['size'] == 4
        assert result['total_length'] == 8
        assert result['characters'] == {'A', 'T', 'C', 'G'}
    
    def test_analyze_empty_raises_error(self):
        """Test that analyzing empty data works."""
        result = analyze_alphabet("")
        assert result['size'] == 0
        assert result['total_length'] == 0
        assert result['entropy'] == 0.0
    
    def test_analyze_single_character(self):
        """Test analysis of single character."""
        result = analyze_alphabet("A")
        assert result['size'] == 1
        assert result['total_length'] == 1
        assert result['entropy'] == 0.0  # Single character has no entropy
    
    def test_entropy_calculation(self):
        """Test entropy calculation."""
        # Uniform distribution should have high entropy
        result_uniform = analyze_alphabet("ABCD")
        
        # Skewed distribution should have lower entropy
        result_skewed = analyze_alphabet("AAAB")
        
        assert result_uniform['entropy'] > result_skewed['entropy']


class TestSequenceDetection:
    """Test biological sequence detection functions."""
    
    def test_dna_sequence_detection(self):
        """Test DNA sequence detection."""
        assert is_dna_sequence("ATCGATCG")
        assert is_dna_sequence("atcgatcg")  # Case insensitive
        assert not is_dna_sequence("ATCGATCGN")  # With N (no longer allowed)
        assert is_dna_sequence(b"ATCG")  # Bytes input
        
        assert not is_dna_sequence("ATCGX")  # Invalid character
        assert not is_dna_sequence("PROTEIN")
        assert not is_dna_sequence("12345")
    
    def test_protein_sequence_detection(self):
        """Test protein sequence detection."""
        assert is_protein_sequence("ACDEFGHIKLMNPQRSTVWY")  # Standard amino acids
        assert is_protein_sequence("acdefg")  # Case insensitive
        assert is_protein_sequence("ACDEFGX")  # With X (unknown)
        assert is_protein_sequence(b"PROTEIN")  # Bytes input
        
        # Note: "ATCG" is actually valid for both DNA and protein sequences
        # since A, T, C, G are valid amino acid codes, so we test with
        # sequences that are clearly not DNA
        assert not is_protein_sequence("12345")
        assert not is_protein_sequence("PROTEINZ123")  # Invalid characters
    
    def test_sequence_type_detection(self):
        """Test sequence type detection."""
        assert detect_sequence_type("ATCGATCG") == 'dna'
        assert detect_sequence_type("ACDEFGHIKLMNPQRSTVWY") == 'protein'  # Valid amino acids
        assert detect_sequence_type("Hello World") == 'text'
        assert detect_sequence_type(b"\x80\x81\x82") == 'binary'
        
        # Edge cases
        assert detect_sequence_type("") == 'text'  # Empty string is text
    
    def test_sequence_detection_with_invalid_bytes(self):
        """Test sequence detection with invalid byte sequences."""
        # Non-ASCII bytes should return 'binary'
        assert detect_sequence_type(b"\xff\xfe\xfd") == 'binary'
        
        # Valid ASCII should be processed
        assert detect_sequence_type(b"ATCG") == 'dna'


class TestBinaryFileIO:
    """Test binary file reading functions."""
    
    def test_read_factors_binary_file_valid(self):
        """Test reading a valid binary factors file."""
        # Create mock binary data: 2 factors
        # Factor 1: pos=0, len=3, ref=1
        # Factor 2: pos=3, len=2, ref=2
        binary_data = struct.pack('<QQQ', 0, 3, 1) + struct.pack('<QQQ', 3, 2, 2)
        
        with tempfile.NamedTemporaryFile(delete=False) as f:
            f.write(binary_data)
            temp_path = f.name
        
        try:
            factors = read_factors_binary_file(temp_path)
            expected = [(0, 3, 1), (3, 2, 2)]
            assert factors == expected
        finally:
            os.unlink(temp_path)
    
    def test_read_factors_binary_file_invalid_size(self):
        """Test reading binary file with invalid size."""
        # Create invalid binary data (not multiple of 24)
        binary_data = b'12345678901234567890123'  # 23 bytes
        
        with tempfile.NamedTemporaryFile(delete=False) as f:
            f.write(binary_data)
            temp_path = f.name
        
        try:
            with pytest.raises(NoLZSSError):
                read_factors_binary_file(temp_path)
        finally:
            os.unlink(temp_path)
    
    def test_read_factors_binary_file_nonexistent(self):
        """Test reading non-existent binary file."""
        with pytest.raises(NoLZSSError):
            read_factors_binary_file("nonexistent.bin")


class TestPlotting:
    """Test plotting functions."""
    
    def test_plot_factor_lengths_with_list(self):
        """Test plotting with factor list."""
        factors = [(0, 3, 1), (3, 2, 2), (5, 4, 3)]
        
        # Test that function works when matplotlib is available
        # or gracefully handles when it's not available
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            plot_factor_lengths(factors, show_plot=False)
            # Function should complete without raising exceptions
    
    def test_plot_factor_lengths_empty_list(self):
        """Test plotting with empty factor list."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            with pytest.raises(ValueError):
                plot_factor_lengths([])
    
    def test_plot_factor_lengths_invalid_type(self):
        """Test plotting with invalid input type."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            with pytest.raises(TypeError):
                plot_factor_lengths(123)


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_very_long_sequences(self):
        """Test with very long sequences."""
        long_dna = "ATCG" * 1000
        assert is_dna_sequence(long_dna)
        
        result = analyze_alphabet(long_dna)
        assert result['size'] == 4
        assert result['total_length'] == 4000
    
    def test_mixed_case_sequences(self):
        """Test with mixed case sequences."""
        mixed_dna = "AtCgAtCg"
        assert is_dna_sequence(mixed_dna)
        
        mixed_protein = "PrOtEiN"
        assert is_protein_sequence(mixed_protein)
    
    def test_sequences_with_numbers(self):
        """Test sequences with numbers (should fail biological detection)."""
        assert not is_dna_sequence("ATCG123")
        assert not is_protein_sequence("PROTEIN123")
        assert detect_sequence_type("ATCG123") == 'text'
    
    def test_unicode_in_sequences(self):
        """Test unicode characters in sequences."""
        unicode_text = "héllo wørld"
        assert detect_sequence_type(unicode_text) == 'text'
        
        result = analyze_alphabet(unicode_text)
        assert result['size'] > len(set("hello world"))  # Should have more unique chars


if __name__ == "__main__":
    # Run tests without pytest
    import traceback
    
    test_classes = [
        TestValidateInput, TestAnalyzeAlphabet,
        TestSequenceDetection, TestBinaryFileIO, TestPlotting, TestEdgeCases
    ]
    
    total_tests = 0
    passed_tests = 0
    
    for test_class in test_classes:
        instance = test_class()
        methods = [method for method in dir(instance) if method.startswith('test_')]
        
        for method_name in methods:
            total_tests += 1
            try:
                method = getattr(instance, method_name)
                method()
                passed_tests += 1
                print(f"{test_class.__name__}.{method_name}")
            except Exception as e:
                print(f"Warning: {test_class.__name__}.{method_name}: {e}")
                traceback.print_exc()
    
    print(f"\nTests passed: {passed_tests}/{total_tests}")
    if passed_tests == total_tests:
        print("All tests passed!")
    else:
        print(f"Failed tests: {total_tests - passed_tests}")
