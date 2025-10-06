"""
Test module for reference sequence factorization functionality.

Tests the factorize_dna_w_reference_seq functions that allow factorization
of a target sequence using a reference sequence with reverse complement awareness,
as well as the general factorize_w_reference functions for text/amino acid sequences.
Factor positions are absolute positions in the combined reference+target string.
"""

import os
import tempfile
import sys
from pathlib import Path

# Add source to path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    import pytest
    PYTEST_AVAILABLE = True
except ImportError:
    PYTEST_AVAILABLE = False

def cpp_bindings_available():
    """Check that the C++ bindings are available."""
    try:
        # Try importing the built extension first
        build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
        if os.path.exists(build_path):
            sys.path.insert(0, build_path)
            import _noLZSS
            assert hasattr(_noLZSS, 'factorize_dna_w_reference_seq') and hasattr(_noLZSS, 'factorize_dna_w_reference_seq_file')
            assert hasattr(_noLZSS, 'factorize_w_reference') and hasattr(_noLZSS, 'factorize_w_reference_file')
        
        # Fallback to installed package
        from noLZSS._noLZSS import factorize_dna_w_reference_seq, factorize_dna_w_reference_seq_file
        from noLZSS._noLZSS import factorize_w_reference, factorize_w_reference_file
        # If we reach here, import was successful
        return True
    except ImportError:
        return False

def test_basic_dna_reference_factorization():
    """Test basic DNA reference sequence factorization with absolute factor positions."""
    if not cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_func = _noLZSS.factorize_dna_w_reference_seq
    else:
        from noLZSS._noLZSS import factorize_dna_w_reference_seq as factorize_func
    
    # Test case: reference contains patterns that target can reference
    reference = "ATCGATCGATCG"
    target = "GATCGATC"  # Should be able to reference patterns in reference
    
    factors = factorize_func(reference, target)
    
    # Verify we got factors
    assert len(factors) > 0, "Should produce at least one factor"
    
    # Verify factor structure and absolute positioning
    for factor in factors:
        assert len(factor) == 4, "Each factor should have 4 elements (start, length, ref, is_rc)"
        start, length, ref, is_rc = factor
        assert isinstance(start, int) and start >= 0, "Start should be non-negative integer (absolute position)"
        assert isinstance(length, int) and length > 0, "Length should be positive integer"
        assert isinstance(ref, int) and ref >= 0, "Ref should be non-negative integer"
        assert isinstance(is_rc, bool), "is_rc should be boolean"
        
        # Factor start positions should be in the target sequence range
        # Target starts at position len(reference) + 1 (after reference + sentinel)
        target_start_pos = len(reference) + 1
        target_end_pos = target_start_pos + len(target)
        assert target_start_pos <= start < target_end_pos, f"Factor start {start} should be within target range [{target_start_pos}, {target_end_pos})"
    
    print(f"✓ Basic reference factorization test passed: {len(factors)} factors")

def test_dna_file_output():
    """Test DNA file output functionality."""
    if not cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_file_func = _noLZSS.factorize_dna_w_reference_seq_file
    else:
        from noLZSS._noLZSS import factorize_dna_w_reference_seq_file as factorize_file_func
    
    reference = "ATCGATCGATCGATCG"
    target = "GATCGATCGATC"
    
    with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as f:
        output_path = f.name
    
    try:
        num_factors = factorize_file_func(reference, target, output_path)
        
        # Verify file was created and has content
        assert os.path.exists(output_path), "Output file should exist"
        file_size = os.path.getsize(output_path)
        assert file_size > 0, "Output file should have content"
        
        # Check that the file size makes sense (header + factors)
        # Each factor is 3 * 8 bytes = 24 bytes, plus header
        expected_min_size = 32  # Header size
        assert file_size >= expected_min_size, f"File size {file_size} seems too small"
        
        print(f"✓ File output test passed: {num_factors} factors, {file_size} bytes")
        
    finally:
        if os.path.exists(output_path):
            os.unlink(output_path)

def test_dna_edge_cases():
    """Test DNA edge cases and error conditions with absolute factor positions."""
    if not cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_func = _noLZSS.factorize_dna_w_reference_seq
    else:
        from noLZSS._noLZSS import factorize_dna_w_reference_seq as factorize_func
    
    # Test with minimal sequences
    reference = "A"
    target = "T"
    factors = factorize_func(reference, target)
    assert len(factors) > 0, "Should handle minimal sequences"
    
    # Verify factor positions are absolute
    for factor in factors:
        start, length, ref, is_rc = factor
        target_start_pos = len(reference) + 1  # After reference + sentinel
        assert start >= target_start_pos, f"Factor start {start} should be at or after target start position {target_start_pos}"
    
    # Test with identical sequences
    reference = "ATCG"
    target = "ATCG"
    factors = factorize_func(reference, target)
    assert len(factors) > 0, "Should handle identical sequences"
    
    # Verify factor positions for identical case
    for factor in factors:
        start, length, ref, is_rc = factor
        target_start_pos = len(reference) + 1
        assert start >= target_start_pos, f"Factor start {start} should be at or after target start position {target_start_pos}"
    
    print("✓ Edge cases test passed")

def test_dna_reverse_complement():
    """Test DNA reverse complement functionality with absolute factor positions."""
    if not cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_func = _noLZSS.factorize_dna_w_reference_seq
    else:
        from noLZSS._noLZSS import factorize_dna_w_reference_seq as factorize_func
    
    # Create a case where target should match reverse complement of reference
    reference = "ATCGATCG"
    target = "CGATCGAT"  # Reverse complement of reference
    
    factors = factorize_func(reference, target)
    
    # Check if any factors are reverse complement matches
    has_rc_match = any(factor[3] for factor in factors)  # factor[3] is is_rc
    
    # Verify factor positions are absolute
    for factor in factors:
        start, length, ref, is_rc = factor
        target_start_pos = len(reference) + 1
        assert start >= target_start_pos, f"Factor start {start} should be at or after target start position {target_start_pos}"
    
    print(f"✓ Reverse complement test: found RC matches = {has_rc_match}")


def test_general_reference_text():
    """Test general reference factorization with text sequences."""
    if not cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_func = _noLZSS.factorize_w_reference
        factorize_file_func = _noLZSS.factorize_w_reference_file
    else:
        from noLZSS._noLZSS import factorize_w_reference as factorize_func
        from noLZSS._noLZSS import factorize_w_reference_file as factorize_file_func
    
    # Test with text that has matching patterns
    reference = "hello world"
    target = "world hello"
    
    factors = factorize_func(reference, target)
    
    # Verify we got factors
    assert len(factors) > 0, "Should produce at least one factor"
    
    # Verify factor structure (no is_rc field for general factorization)
    for factor in factors:
        assert len(factor) == 3, "Each factor should have 3 elements (start, length, ref)"
        start, length, ref = factor
        assert isinstance(start, int) and start >= 0, "Start should be non-negative integer"
        assert isinstance(length, int) and length > 0, "Length should be positive integer"
        assert isinstance(ref, int) and ref >= 0, "Ref should be non-negative integer"
        
        # Factor start positions should be in the target sequence range
        target_start_pos = len(reference) + 1
        assert start >= target_start_pos, f"Factor start {start} should be at or after target start position {target_start_pos}"
    
    # Test file function
    with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as f:
        output_path = f.name
    
    try:
        num_factors = factorize_file_func(reference, target, output_path)
        assert num_factors == len(factors), "File function should produce same number of factors"
        assert os.path.exists(output_path), "Output file should exist"
        
        print(f"✓ General text test passed: {num_factors} factors")
        
    finally:
        if os.path.exists(output_path):
            os.unlink(output_path)


def test_general_reference_amino_acids():
    """Test general reference factorization with amino acid sequences."""
    if not cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_func = _noLZSS.factorize_w_reference
    else:
        from noLZSS._noLZSS import factorize_w_reference as factorize_func
    
    # Test with amino acid sequences
    reference = "ACDEFGHIKLMNPQRSTVWY"  # All 20 standard amino acids
    target = "ACDEFGHIKLMNPQR"  # Subset that should match reference
    
    factors = factorize_func(reference, target)
    
    # Verify we got factors and they make sense
    assert len(factors) > 0, "Should produce at least one factor"
    
    # Should have fewer factors than characters due to matching
    assert len(factors) <= len(target), "Should have compression"
    
    # Verify all factors are properly positioned
    for factor in factors:
        start, length, ref = factor
        target_start_pos = len(reference) + 1
        assert start >= target_start_pos, f"Factor start {start} should be in target region"
        assert start < target_start_pos + len(target), f"Factor should not exceed target region"
    
    print(f"✓ Amino acid test passed: {len(factors)} factors for {len(target)} characters")


def test_general_vs_dna_differences():
    """Test that general and DNA functions behave differently for DNA sequences."""
    if not cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        dna_func = _noLZSS.factorize_dna_w_reference_seq
        general_func = _noLZSS.factorize_w_reference
    else:
        from noLZSS._noLZSS import factorize_dna_w_reference_seq as dna_func
        from noLZSS._noLZSS import factorize_w_reference as general_func
    
    # Use DNA sequences where reverse complement could make a difference
    reference = "ATCGATCG"
    target = "CGATCGAT"  # Reverse complement of reference
    
    dna_factors = dna_func(reference, target)
    general_factors = general_func(reference, target)
    
    # DNA factors should have 4 elements (including is_rc)
    for factor in dna_factors:
        assert len(factor) == 4, "DNA factors should have 4 elements"
        assert isinstance(factor[3], bool), "Fourth element should be boolean is_rc"
    
    # General factors should have 3 elements (no is_rc)
    for factor in general_factors:
        assert len(factor) == 3, "General factors should have 3 elements"
    
    # Check if DNA function found reverse complement matches
    has_rc_matches = any(factor[3] for factor in dna_factors)
    
    print(f"✓ Function differences test: DNA has RC matches = {has_rc_matches}")
    print(f"  DNA factors: {len(dna_factors)}, General factors: {len(general_factors)}")


def test_pattern_matching_validation():
    """Test that factors correctly reference patterns in the reference sequence."""
    if not cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_func = _noLZSS.factorize_w_reference
    else:
        from noLZSS._noLZSS import factorize_w_reference as factorize_func
    
    reference = "abcdefghijk"
    target = "defghijkabc"  # Contains substrings from reference
    combined = reference + '\x01' + target  # How the C++ function combines them
    
    factors = factorize_func(reference, target)
    
    # Validate that each factor correctly references content
    target_start_pos = len(reference) + 1
    
    for factor in factors:
        start, length, ref = factor
        
        # Extract the factor content from the combined string
        factor_content = combined[start:start + length]
        
        # Extract the reference content
        ref_content = combined[ref:ref + length]
        
        # They should match
        assert factor_content == ref_content, f"Factor content '{factor_content}' should match reference content '{ref_content}'"
        
        # Factor should be in target region
        assert start >= target_start_pos, f"Factor start {start} should be in target region"
    
    print(f"✓ Pattern matching validation passed for {len(factors)} factors")


def test_fasta_reference_binary_output():
    """Test DNA reference sequence factorization from FASTA files with binary output."""
    if not cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Check if our new function is available
    try:
        from noLZSS.genomics.fasta import write_factors_dna_w_reference_fasta_files_to_binary
    except ImportError:
        print("Skipping test - write_factors_dna_w_reference_fasta_files_to_binary not available")
        return
    
    # Create test FASTA files
    with tempfile.TemporaryDirectory() as temp_dir:
        ref_fasta_path = os.path.join(temp_dir, "reference.fasta")
        target_fasta_path = os.path.join(temp_dir, "target.fasta")
        output_path = os.path.join(temp_dir, "factors.bin")
        
        # Create reference FASTA with multiple sequences
        with open(ref_fasta_path, 'w') as f:
            f.write(">ref1\n")
            f.write("ATCGATCGATCG\n")  # Reference sequence 1
            f.write(">ref2\n")
            f.write("GCTAGCTAGCTA\n")  # Reference sequence 2
        
        # Create target FASTA with sequences that should match parts of reference
        with open(target_fasta_path, 'w') as f:
            f.write(">target1\n")
            f.write("ATCGATCG\n")  # Should match beginning of ref1
            f.write(">target2\n")
            f.write("GCTAGCTA\n")  # Should match beginning of ref2
        
        # Test the function
        try:
            num_factors = write_factors_dna_w_reference_fasta_files_to_binary(
                ref_fasta_path, target_fasta_path, output_path
            )
            
            # Verify output file was created
            assert os.path.exists(output_path), "Output file was not created"
            
            # Verify file is not empty
            assert os.path.getsize(output_path) > 0, "Output file is empty"
            
            # Should have at least some factors
            assert num_factors > 0, f"Expected positive number of factors, got {num_factors}"
            
            print(f"FASTA reference binary output test: {num_factors} factors written successfully")
            
        except Exception as e:
            print(f"FASTA reference test failed: {e}")
            raise
    
    # Test with invalid files
    with tempfile.TemporaryDirectory() as temp_dir:
        invalid_ref = os.path.join(temp_dir, "nonexistent_ref.fasta")
        invalid_target = os.path.join(temp_dir, "nonexistent_target.fasta")
        output_path = os.path.join(temp_dir, "factors.bin")
        
        # Test with non-existent reference file
        try:
            write_factors_dna_w_reference_fasta_files_to_binary(
                invalid_ref, target_fasta_path, output_path
            )
            assert False, "Should have raised FileNotFoundError for missing reference file"
        except FileNotFoundError:
            pass  # Expected
        
        # Test with non-existent target file
        try:
            write_factors_dna_w_reference_fasta_files_to_binary(
                ref_fasta_path, invalid_target, output_path
            )
            assert False, "Should have raised FileNotFoundError for missing target file"
        except FileNotFoundError:
            pass  # Expected
    
    print("FASTA reference binary output test completed successfully")


def main():
    """Run all tests."""
    tests = [
        test_basic_dna_reference_factorization,
        test_dna_file_output,
        test_dna_edge_cases,
        test_dna_reverse_complement,
        test_general_reference_text,
        test_general_reference_amino_acids,
        test_general_vs_dna_differences,
        test_pattern_matching_validation,
        test_fasta_reference_binary_output,
    ]
    
    print("Running reference sequence factorization tests...")
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            test()  # Tests now use assertions and will raise exceptions on failure
            passed += 1
            print(f"{test.__name__} passed")
        except Exception as e:
            print(f"{test.__name__} failed with exception: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\nTest Results: {passed}/{total} passed")
    
    if passed == total:
        print("All tests passed!")
        return 0
    else:
        print("Some tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())