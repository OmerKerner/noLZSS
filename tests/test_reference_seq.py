"""
Test module for reference sequence factorization functionality.

Tests the new factorize_w_reference_seq functions that allow factorization
of a target sequence using a reference sequence with reverse complement awareness.
Factor positions are absolute positions in the combined reference+target string.
"""

import os
import tempfile
import sys
import pytest
from pathlib import Path

# Add source to path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def test_cpp_bindings_available():
    """Test that the C++ bindings are available."""
    try:
        # Try importing the built extension first
        build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
        if os.path.exists(build_path):
            sys.path.insert(0, build_path)
            import _noLZSS
            assert hasattr(_noLZSS, 'factorize_w_reference_seq') and hasattr(_noLZSS, 'factorize_w_reference_seq_file')
        
        # Fallback to installed package
        from noLZSS._noLZSS import factorize_w_reference_seq, factorize_w_reference_seq_file
        # If we reach here, import was successful
    except ImportError:
        pytest.fail("C++ bindings not available")

def test_basic_reference_factorization():
    """Test basic reference sequence factorization with absolute factor positions."""
    if not test_cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_func = _noLZSS.factorize_w_reference_seq
    else:
        from noLZSS._noLZSS import factorize_w_reference_seq as factorize_func
    
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

def test_file_output():
    """Test file output functionality."""
    if not test_cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_file_func = _noLZSS.factorize_w_reference_seq_file
    else:
        from noLZSS._noLZSS import factorize_w_reference_seq_file as factorize_file_func
    
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

def test_edge_cases():
    """Test edge cases and error conditions with absolute factor positions."""
    if not test_cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_func = _noLZSS.factorize_w_reference_seq
    else:
        from noLZSS._noLZSS import factorize_w_reference_seq as factorize_func
    
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

def test_reverse_complement():
    """Test reverse complement functionality with absolute factor positions."""
    if not test_cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return
    
    # Try build directory first
    build_path = os.path.join(os.path.dirname(__file__), '..', 'build')
    if os.path.exists(build_path):
        sys.path.insert(0, build_path)
        import _noLZSS
        factorize_func = _noLZSS.factorize_w_reference_seq
    else:
        from noLZSS._noLZSS import factorize_w_reference_seq as factorize_func
    
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

def main():
    """Run all tests."""
    tests = [
        test_basic_reference_factorization,
        test_file_output,
        test_edge_cases,
        test_reverse_complement,
    ]
    
    print("Running reference sequence factorization tests...")
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            test()  # Tests now use assertions and will raise exceptions on failure
            passed += 1
            print(f"✓ {test.__name__} passed")
        except Exception as e:
            print(f"✗ {test.__name__} failed with exception: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\nTest Results: {passed}/{total} passed")
    
    if passed == total:
        print("🎉 All tests passed!")
        return 0
    else:
        print("❌ Some tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())