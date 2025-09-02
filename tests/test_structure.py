"""
Test the new modular structure of noLZSS.
"""

import pytest
import sys
import os

# Add the src directory to the path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def test_utils_import():
    """Test that utils module can be imported and basic functions work."""
    from noLZSS.utils import validate_input, analyze_alphabet, detect_sequence_type
    
    # Test validate_input
    data = validate_input("ATCG")
    assert isinstance(data, bytes)
    assert data == b"ATCG"
    
    # Test analyze_alphabet
    result = analyze_alphabet("ATCGATCG")
    assert result['size'] == 4  # A, T, C, G
    assert result['total_length'] == 8
    assert 'distribution' in result
    
    # Test detect_sequence_type
    seq_type = detect_sequence_type("ATCGATCG")
    assert seq_type == 'dna'


def test_core_structure():
    """Test that core module has the expected structure."""
    from noLZSS import core
    
    # Check that the functions exist
    assert hasattr(core, 'factorize')
    assert hasattr(core, 'factorize_file')
    assert hasattr(core, 'count_factors')
    assert hasattr(core, 'factorize_with_info')


def test_genomics_import():
    """Test that genomics subpackage can be imported."""
    from noLZSS import genomics
    # Just test that it imports without error
    # The actual functionality will be implemented later


def test_main_package_imports():
    """Test that the main package exposes the expected API."""
    import noLZSS
    
    # Check that enhanced functions are available
    assert hasattr(noLZSS, 'analyze_alphabet')
    assert hasattr(noLZSS, 'detect_sequence_type')
    assert hasattr(noLZSS, 'is_dna_sequence')
    assert hasattr(noLZSS, 'is_protein_sequence')


if __name__ == "__main__":
    # Run basic tests
    test_utils_import()
    test_core_structure() 
    test_genomics_import()
    test_main_package_imports()
    print("All structure tests passed!")
