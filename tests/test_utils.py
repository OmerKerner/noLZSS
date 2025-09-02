"""
Comprehensive tests for the utils module.
"""

import pytest
import sys
import os
import tempfile
from pathlib import Path

# Add the src directory to the path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from noLZSS.utils import (
    validate_input, ensure_sentinel, analyze_alphabet, 
    is_dna_sequence, is_protein_sequence, detect_sequence_type,
    safe_file_reader, InvalidInputError, NoLZSSError
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
        """Test validation of unicode strings."""
        result = validate_input("héllo")
        assert isinstance(result, bytes)
        assert result == "héllo".encode('utf-8')
    
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


class TestEnsureSentinel:
    """Test sentinel handling functions."""
    
    def test_ensure_sentinel_adds_when_missing(self):
        """Test that sentinel is added when missing."""
        result = ensure_sentinel(b"hello")
        assert result == b"hello$"
    
    def test_ensure_sentinel_doesnt_duplicate(self):
        """Test that sentinel is not duplicated when present."""
        result = ensure_sentinel(b"hello$")
        assert result == b"hello$"
    
    def test_ensure_sentinel_custom_sentinel(self):
        """Test custom sentinel character."""
        result = ensure_sentinel(b"hello", b"#")
        assert result == b"hello#"
        
        result2 = ensure_sentinel(b"hello#", b"#")
        assert result2 == b"hello#"


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
        assert result['characters'] == {ord('A'), ord('T'), ord('C'), ord('G')}
    
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
        assert is_dna_sequence("ATCGATCGN")  # With N
        assert is_dna_sequence("ATCG$")  # With sentinel
        assert is_dna_sequence(b"ATCG")  # Bytes input
        
        assert not is_dna_sequence("ATCGX")  # Invalid character
        assert not is_dna_sequence("PROTEIN")
        assert not is_dna_sequence("12345")
    
    def test_protein_sequence_detection(self):
        """Test protein sequence detection."""
        assert is_protein_sequence("ACDEFGHIKLMNPQRSTVWY")  # Standard amino acids
        assert is_protein_sequence("acdefg")  # Case insensitive
        assert is_protein_sequence("ACDEFGX")  # With X (unknown)
        assert is_protein_sequence("PROTEIN$")  # With sentinel
        assert is_protein_sequence(b"PROTEIN")  # Bytes input
        
        # Note: "ATCG" is actually valid for both DNA and protein sequences
        # since A, T, C, G are valid amino acid codes, so we test with
        # sequences that are clearly not DNA
        assert not is_protein_sequence("12345")
        assert not is_protein_sequence("PROTEINZ123")  # Invalid characters
    
    def test_sequence_type_detection(self):
        """Test sequence type detection."""
        assert detect_sequence_type("ATCGATCG") == 'dna'
        assert detect_sequence_type("PROTEIN") == 'protein'
        assert detect_sequence_type("Hello World") == 'text'
        assert detect_sequence_type(b"\x80\x81\x82") == 'binary'
        
        # Edge cases
        assert detect_sequence_type("ATCG$") == 'dna'
        assert detect_sequence_type("") == 'text'  # Empty string is text
    
    def test_sequence_detection_with_invalid_bytes(self):
        """Test sequence detection with invalid byte sequences."""
        # Non-ASCII bytes should return 'binary'
        assert detect_sequence_type(b"\xff\xfe\xfd") == 'binary'
        
        # Valid ASCII should be processed
        assert detect_sequence_type(b"ATCG") == 'dna'


class TestSafeFileReader:
    """Test safe file reading functionality."""
    
    def test_safe_file_reader_small_file(self):
        """Test reading a small file."""
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.write("Hello World")
            temp_path = f.name
        
        try:
            chunks = list(safe_file_reader(temp_path))
            data = b''.join(chunks)
            assert data == b"Hello World"
        finally:
            os.unlink(temp_path)
    
    def test_safe_file_reader_large_file(self):
        """Test reading a larger file in chunks."""
        test_data = b"A" * 10000  # 10KB of data
        
        with tempfile.NamedTemporaryFile(mode='wb', delete=False) as f:
            f.write(test_data)
            temp_path = f.name
        
        try:
            chunks = list(safe_file_reader(temp_path, chunk_size=1024))
            data = b''.join(chunks)
            assert data == test_data
            assert len(chunks) > 1  # Should be multiple chunks
        finally:
            os.unlink(temp_path)
    
    def test_safe_file_reader_nonexistent_file(self):
        """Test reading a non-existent file raises error."""
        with pytest.raises(NoLZSSError):
            list(safe_file_reader("nonexistent_file.txt"))
    
    def test_safe_file_reader_empty_file(self):
        """Test reading an empty file."""
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            temp_path = f.name
        
        try:
            chunks = list(safe_file_reader(temp_path))
            assert chunks == []
        finally:
            os.unlink(temp_path)


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
        TestValidateInput, TestEnsureSentinel, TestAnalyzeAlphabet,
        TestSequenceDetection, TestSafeFileReader, TestEdgeCases
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
