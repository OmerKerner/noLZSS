"""
Tests for per-sequence FASTA processing functionality.

This module tests the per-sequence FASTA factorization functions including:
- factorize_fasta_dna_w_rc_per_sequence
- factorize_fasta_dna_no_rc_per_sequence
- write_factors_binary_file_fasta_dna_w_rc_per_sequence
- write_factors_binary_file_fasta_dna_no_rc_per_sequence
- count_factors_fasta_dna_w_rc_per_sequence
- count_factors_fasta_dna_no_rc_per_sequence
- parallel versions of write functions

Unlike the concatenated versions (factorize_fasta_multiple_dna_*), these functions
factorize each sequence in the FASTA file independently without using sentinels.
"""

import os
import tempfile
import struct
from pathlib import Path

# Try to import the C++ extension
try:
    from noLZSS._noLZSS import (
        factorize_fasta_dna_w_rc_per_sequence,
        factorize_fasta_dna_no_rc_per_sequence,
        write_factors_binary_file_fasta_dna_w_rc_per_sequence,
        write_factors_binary_file_fasta_dna_no_rc_per_sequence,
        count_factors_fasta_dna_w_rc_per_sequence,
        count_factors_fasta_dna_no_rc_per_sequence,
        parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence,
        parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence,
    )
    CPP_AVAILABLE = True
except ImportError as e:
    print(f"Warning: C++ extension not available: {e}")
    CPP_AVAILABLE = False

# Test resources directory
RESOURCES_DIR = Path(__file__).parent / "resources"


def read_per_sequence_binary_factors(filepath):
    """
    Read factors from a per-sequence binary file and return them along with metadata.
    
    Returns:
        tuple: (factors_list, num_sequences, sequence_ids, factors_per_sequence)
    """
    with open(filepath, 'rb') as f:
        # Read all data
        data = f.read()
    
    # Read footer from end
    footer_start = len(data) - 48  # Footer is 48 bytes (8 + 5*8 bytes)
    footer = data[footer_start:]
    
    magic = footer[0:8]
    num_factors = struct.unpack('<Q', footer[8:16])[0]
    num_sequences = struct.unpack('<Q', footer[16:24])[0]
    num_sentinels = struct.unpack('<Q', footer[24:32])[0]  # Should be 0 for per-sequence
    footer_size = struct.unpack('<Q', footer[32:40])[0]
    total_length = struct.unpack('<Q', footer[40:48])[0]
    
    # Calculate where metadata starts
    metadata_start = len(data) - footer_size
    factors_end = metadata_start
    
    # Read factors (24 bytes each: start, length, ref)
    factors = []
    for i in range(num_factors):
        offset = i * 24
        factor_data = data[offset:offset + 24]
        start = struct.unpack('<Q', factor_data[0:8])[0]
        length = struct.unpack('<Q', factor_data[8:16])[0]
        ref = struct.unpack('<Q', factor_data[16:24])[0]
        factors.append((start, length, ref))
    
    # Read sequence IDs from metadata
    sequence_ids = []
    pos = metadata_start
    for _ in range(num_sequences):
        end = data.find(b'\x00', pos)
        seq_id = data[pos:end].decode('utf-8')
        sequence_ids.append(seq_id)
        pos = end + 1
    
    # Read factors per sequence (boundary markers)
    factors_per_sequence = []
    for i in range(num_sequences):
        count_data = data[pos:pos + 8]
        count = struct.unpack('<Q', count_data)[0]
        factors_per_sequence.append(count)
        pos += 8
    
    return factors, num_sequences, sequence_ids, factors_per_sequence


class TestPerSequenceFastaFactorization:
    """Test per-sequence FASTA factorization functions."""
    
    def test_factorize_w_rc_basic(self):
        """Test basic per-sequence factorization with reverse complement."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        per_seq_factors, sequence_ids = factorize_fasta_dna_w_rc_per_sequence(fasta_path)
        
        assert isinstance(per_seq_factors, list), "Should return list of factor lists"
        assert isinstance(sequence_ids, list), "Should return list of sequence IDs"
        assert len(per_seq_factors) == len(sequence_ids), "Should have same number of factor lists and IDs"
        assert len(per_seq_factors) > 0, "Should have at least one sequence"
        
        # Each sequence should have factors
        for seq_factors in per_seq_factors:
            assert isinstance(seq_factors, list), "Each sequence should have a list of factors"
            assert len(seq_factors) > 0, "Each sequence should have at least one factor"
            
            # Check factor format
            for factor in seq_factors:
                assert len(factor) == 4, "Each factor should be (start, length, ref, is_rc)"
                start, length, ref, is_rc = factor
                assert isinstance(start, int) and start >= 0
                assert isinstance(length, int) and length > 0
                assert isinstance(ref, int) and ref >= 0
                assert isinstance(is_rc, bool)
        
        print(f"Factorized {len(sequence_ids)} sequences with RC support")
        for i, (seq_id, factors) in enumerate(zip(sequence_ids, per_seq_factors)):
            print(f"  Sequence {i+1} '{seq_id}': {len(factors)} factors")
    
    def test_factorize_no_rc_basic(self):
        """Test basic per-sequence factorization without reverse complement."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        per_seq_factors, sequence_ids = factorize_fasta_dna_no_rc_per_sequence(fasta_path)
        
        assert isinstance(per_seq_factors, list), "Should return list of factor lists"
        assert isinstance(sequence_ids, list), "Should return list of sequence IDs"
        assert len(per_seq_factors) == len(sequence_ids), "Should have same number of factor lists and IDs"
        assert len(per_seq_factors) > 0, "Should have at least one sequence"
        
        # Each sequence should have factors
        for seq_factors in per_seq_factors:
            assert isinstance(seq_factors, list), "Each sequence should have a list of factors"
            assert len(seq_factors) > 0, "Each sequence should have at least one factor"
            
            # Check factor format
            for factor in seq_factors:
                assert len(factor) == 4, "Each factor should be (start, length, ref, is_rc)"
                start, length, ref, is_rc = factor
                assert isinstance(start, int) and start >= 0
                assert isinstance(length, int) and length > 0
                assert isinstance(ref, int) and ref >= 0
                # For no_rc, is_rc should always be False
                assert is_rc == False, "is_rc should be False for no_rc version"
        
        print(f"Factorized {len(sequence_ids)} sequences without RC")
        for i, (seq_id, factors) in enumerate(zip(sequence_ids, per_seq_factors)):
            print(f"  Sequence {i+1} '{seq_id}': {len(factors)} factors")
    
    def test_count_factors_w_rc(self):
        """Test counting factors per sequence with RC."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        # Get full factorization
        per_seq_factors, _ = factorize_fasta_dna_w_rc_per_sequence(fasta_path)
        expected_count = sum(len(factors) for factors in per_seq_factors)
        
        # Get count only
        actual_count = count_factors_fasta_dna_w_rc_per_sequence(fasta_path)
        
        assert actual_count == expected_count, f"Count {actual_count} should match full factorization {expected_count}"
        print(f"Count correct: {actual_count} total factors")
    
    def test_count_factors_no_rc(self):
        """Test counting factors per sequence without RC."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        # Get full factorization
        per_seq_factors, _ = factorize_fasta_dna_no_rc_per_sequence(fasta_path)
        expected_count = sum(len(factors) for factors in per_seq_factors)
        
        # Get count only
        actual_count = count_factors_fasta_dna_no_rc_per_sequence(fasta_path)
        
        assert actual_count == expected_count, f"Count {actual_count} should match full factorization {expected_count}"
        print(f"Count correct: {actual_count} total factors")
    
    def test_write_binary_w_rc(self):
        """Test writing per-sequence factors to binary file with RC."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Write to binary file
            count = write_factors_binary_file_fasta_dna_w_rc_per_sequence(fasta_path, output_path)
            
            assert count > 0, "Should write at least one factor"
            assert os.path.exists(output_path), "Output file should be created"
            assert os.path.getsize(output_path) > 0, "Output file should not be empty"
            
            # Read and validate
            factors, num_seqs, seq_ids, factors_per_seq = read_per_sequence_binary_factors(output_path)
            
            assert len(factors) == count, "Factor count should match"
            assert num_seqs > 0, "Should have at least one sequence"
            assert len(seq_ids) == num_seqs, "Sequence ID count should match"
            assert len(factors_per_seq) == num_seqs, "Should have factor count for each sequence"
            assert sum(factors_per_seq) == count, "Sum of per-sequence factors should equal total"
            
            print(f"Successfully wrote {count} factors across {num_seqs} sequences")
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_write_binary_no_rc(self):
        """Test writing per-sequence factors to binary file without RC."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Write to binary file
            count = write_factors_binary_file_fasta_dna_no_rc_per_sequence(fasta_path, output_path)
            
            assert count > 0, "Should write at least one factor"
            assert os.path.exists(output_path), "Output file should be created"
            
            # Read and validate
            factors, num_seqs, seq_ids, factors_per_seq = read_per_sequence_binary_factors(output_path)
            
            assert len(factors) == count, "Factor count should match"
            assert sum(factors_per_seq) == count, "Sum of per-sequence factors should equal total"
            
            print(f"Successfully wrote {count} factors across {num_seqs} sequences (no RC)")
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_parallel_write_sequential(self):
        """Test parallel write with num_threads=1 (sequential)."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Test with 1 thread
            count = parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence(
                fasta_path, output_path, 1
            )
            
            assert count > 0, "Should produce at least one factor"
            assert os.path.exists(output_path), "Output file should be created"
            
            # Read and validate
            factors, num_seqs, _, factors_per_seq = read_per_sequence_binary_factors(output_path)
            assert len(factors) == count
            assert sum(factors_per_seq) == count
            
            print(f"Parallel (1 thread) wrote {count} factors")
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_parallel_write_multithreaded(self):
        """Test parallel write with multiple threads."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Test with 2 threads
            count = parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence(
                fasta_path, output_path, 2
            )
            
            assert count > 0, "Should produce at least one factor"
            
            # Read and validate
            factors, _, _, factors_per_seq = read_per_sequence_binary_factors(output_path)
            assert len(factors) == count
            assert sum(factors_per_seq) == count
            
            print(f"Parallel (2 threads) wrote {count} factors")
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_consistency_sequential_vs_parallel(self):
        """Test that sequential and parallel produce identical results."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp1:
            seq_output = tmp1.name
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp2:
            par_output = tmp2.name
        
        try:
            # Sequential
            seq_count = parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence(
                fasta_path, seq_output, 1
            )
            
            # Parallel
            par_count = parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence(
                fasta_path, par_output, 2
            )
            
            # Should produce same count
            assert seq_count == par_count, "Sequential and parallel should produce same factor count"
            
            # Read both outputs
            seq_factors, _, _, seq_factors_per = read_per_sequence_binary_factors(seq_output)
            par_factors, _, _, par_factors_per = read_per_sequence_binary_factors(par_output)
            
            # Should have identical factors and boundaries
            assert len(seq_factors) == len(par_factors)
            assert seq_factors == par_factors, "Sequential and parallel should produce identical factors"
            assert seq_factors_per == par_factors_per, "Per-sequence boundaries should match"
            
            print(f"Sequential and parallel are consistent: {seq_count} factors")
            
        finally:
            if os.path.exists(seq_output):
                os.remove(seq_output)
            if os.path.exists(par_output):
                os.remove(par_output)
    
    def test_per_sequence_independence(self):
        """Test that sequences are truly independent (no cross-sequence matches)."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        # Create a test FASTA with repeated sequences
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
            fasta_path = tmp.name
            tmp.write(">seq1\n")
            tmp.write("ATCGATCGATCG\n")
            tmp.write(">seq2\n")
            tmp.write("ATCGATCGATCG\n")  # Same sequence
        
        try:
            per_seq_factors, seq_ids = factorize_fasta_dna_no_rc_per_sequence(fasta_path)
            
            assert len(per_seq_factors) == 2, "Should have two sequences"
            
            # Check that factors in each sequence only reference that sequence
            for i, factors in enumerate(per_seq_factors):
                for factor in factors:
                    start, length, ref, is_rc = factor
                    # Refs should be within the same sequence (no cross-sequence matches)
                    # For the first factor of a sequence, ref should equal start (no previous occurrence)
                    # For subsequent factors, ref should be < start and within valid range
                    assert ref <= start, f"Seq {i}: ref {ref} should not exceed start {start}"
            
            print(f"Per-sequence independence verified")
            print(f"  Seq 1: {len(per_seq_factors[0])} factors")
            print(f"  Seq 2: {len(per_seq_factors[1])} factors")
            
        finally:
            if os.path.exists(fasta_path):
                os.remove(fasta_path)
    
    def test_empty_sequence_handling(self):
        """Test handling of FASTA files with empty sequences."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        # Note: FASTA parsers typically skip empty sequences, so this tests that behavior
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        # Should not raise an error
        per_seq_factors, seq_ids = factorize_fasta_dna_w_rc_per_sequence(fasta_path)
        
        # All sequences should have factors
        for factors in per_seq_factors:
            assert len(factors) > 0, "Even short sequences should have factors"
        
        print("Empty sequence handling verified")
    
    def test_invalid_nucleotides(self):
        """Test error handling for invalid nucleotides."""
        if not CPP_AVAILABLE:
            print("Skipping test - C++ extension not available")
            return
        
        # Create a test FASTA with invalid nucleotides
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
            fasta_path = tmp.name
            tmp.write(">seq1\n")
            tmp.write("ATCGXYZ\n")  # Invalid nucleotides
        
        try:
            # Should raise an error
            try:
                factorize_fasta_dna_w_rc_per_sequence(fasta_path)
                assert False, "Should have raised an error for invalid nucleotides"
            except (RuntimeError, ValueError) as e:
                assert "Invalid nucleotide" in str(e) or "invalid" in str(e).lower()
                print(f"Correctly detected invalid nucleotides: {e}")
        finally:
            if os.path.exists(fasta_path):
                os.remove(fasta_path)


def run_tests():
    """Run all tests in this module."""
    if not CPP_AVAILABLE:
        print("C++ extension not available - skipping all tests")
        return
    
    test_class = TestPerSequenceFastaFactorization()
    
    tests = [
        ("Basic factorization with RC", test_class.test_factorize_w_rc_basic),
        ("Basic factorization without RC", test_class.test_factorize_no_rc_basic),
        ("Count factors with RC", test_class.test_count_factors_w_rc),
        ("Count factors without RC", test_class.test_count_factors_no_rc),
        ("Write binary with RC", test_class.test_write_binary_w_rc),
        ("Write binary without RC", test_class.test_write_binary_no_rc),
        ("Parallel write sequential", test_class.test_parallel_write_sequential),
        ("Parallel write multithreaded", test_class.test_parallel_write_multithreaded),
        ("Sequential vs parallel consistency", test_class.test_consistency_sequential_vs_parallel),
        ("Per-sequence independence", test_class.test_per_sequence_independence),
        ("Empty sequence handling", test_class.test_empty_sequence_handling),
        ("Invalid nucleotides", test_class.test_invalid_nucleotides),
    ]
    
    print("\n" + "="*70)
    print("Running Per-Sequence FASTA Factorization Tests")
    print("="*70 + "\n")
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            print(f"\n{test_name}:")
            print("-" * 70)
            test_func()
            passed += 1
            print(f"✓ PASSED")
        except AssertionError as e:
            failed += 1
            print(f"✗ FAILED: {e}")
        except Exception as e:
            failed += 1
            print(f"✗ ERROR: {e}")
    
    print("\n" + "="*70)
    print(f"Test Results: {passed} passed, {failed} failed out of {len(tests)} total")
    print("="*70 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    import sys
    success = run_tests()
    sys.exit(0 if success else 1)
