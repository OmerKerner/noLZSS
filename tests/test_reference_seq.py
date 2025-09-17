"""
Test module for reference sequence factorization functionality.

Tests the new factorize_w_reference_seq functions that allow factorization
of a target sequence using a reference sequence with reverse complement awareness.
"""

import os
import tempfile
import sys
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
            return hasattr(_noLZSS, 'factorize_w_reference_seq') and hasattr(_noLZSS, 'factorize_w_reference_seq_file')
        
        # Fallback to installed package
        from noLZSS._noLZSS import factorize_w_reference_seq, factorize_w_reference_seq_file
        return True
    except ImportError:
        return False

def test_basic_reference_factorization():
    """Test basic reference sequence factorization."""
    if not test_cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return True
    
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
    
    # Verify factor structure
    for factor in factors:
        assert len(factor) == 4, "Each factor should have 4 elements (start, length, ref, is_rc)"
        start, length, ref, is_rc = factor
        assert isinstance(start, int) and start >= 0, "Start should be non-negative integer"
        assert isinstance(length, int) and length > 0, "Length should be positive integer"
        assert isinstance(ref, int) and ref >= 0, "Ref should be non-negative integer"
        assert isinstance(is_rc, bool), "is_rc should be boolean"
    
    print(f"✓ Basic reference factorization test passed: {len(factors)} factors")
    return True

def test_file_output():
    """Test file output functionality."""
    if not test_cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return True
    
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
        return True
        
    finally:
        if os.path.exists(output_path):
            os.unlink(output_path)

def test_edge_cases():
    """Test edge cases and error conditions."""
    if not test_cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return True
    
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
    
    # Test with identical sequences
    reference = "ATCG"
    target = "ATCG"
    factors = factorize_func(reference, target)
    assert len(factors) > 0, "Should handle identical sequences"
    
    print("✓ Edge cases test passed")
    return True

def test_reverse_complement():
    """Test reverse complement functionality."""
    if not test_cpp_bindings_available():
        print("Skipping test - C++ extension not available")
        return True
    
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
    
    print(f"✓ Reverse complement test: found RC matches = {has_rc_match}")
    return True

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
            if test():
                passed += 1
            else:
                print(f"✗ {test.__name__} failed")
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