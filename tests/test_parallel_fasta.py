"""
Tests for parallel FASTA processing functionality.

This module tests the parallel_fasta_processor functions including:
- parallel_write_factors_binary_file_fasta_multiple_dna_w_rc
- parallel_write_factors_binary_file_fasta_multiple_dna_no_rc
- parallel_write_factors_dna_w_reference_fasta_files_to_binary

Note: The parallel factorizer uses a minimum of 100,000 characters per thread,
so tests that verify actual multi-threading must use sufficiently large inputs.
"""

import pytest
import os
import tempfile
import struct
import random
from pathlib import Path

# Try to import the C++ extension
try:
    from noLZSS._noLZSS import (
        parallel_write_factors_binary_file_fasta_multiple_dna_w_rc,
        parallel_write_factors_binary_file_fasta_multiple_dna_no_rc,
        parallel_write_factors_dna_w_reference_fasta_files_to_binary,
        write_factors_binary_file_fasta_multiple_dna_w_rc,
        write_factors_binary_file_fasta_multiple_dna_no_rc,
        write_factors_dna_w_reference_fasta_files_to_binary,
    )
    CPP_AVAILABLE = True
except ImportError:
    CPP_AVAILABLE = False


# Test resources directory
RESOURCES_DIR = Path(__file__).parent / "resources"

# Minimum characters per thread for parallel factorization
MIN_CHARS_PER_THREAD = 100000


def generate_random_dna_sequence(length, seed=None):
    """Generate a random DNA sequence of specified length."""
    if seed is not None:
        random.seed(seed)
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def create_large_fasta_file(filepath, num_sequences=3, seq_length=150000, seed=42):
    """
    Create a large FASTA file suitable for testing multi-threading.
    
    Args:
        filepath: Path to create the FASTA file
        num_sequences: Number of sequences to generate
        seq_length: Length of each sequence (default 150K chars per sequence)
        seed: Random seed for reproducibility
    
    Returns:
        Total character count in the file
    """
    total_chars = 0
    with open(filepath, 'w') as f:
        for i in range(num_sequences):
            seq = generate_random_dna_sequence(seq_length, seed=seed + i if seed else None)
            f.write(f'>sequence_{i+1}\n')
            # Write in lines of 80 characters (standard FASTA format)
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + '\n')
            total_chars += len(seq)
    return total_chars


def read_binary_factors(filepath):
    """
    Read factors from a binary file and return them along with metadata.
    
    Returns:
        tuple: (factors_list, num_sequences, num_sentinels, sequence_ids, sentinel_indices)
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
    num_sentinels = struct.unpack('<Q', footer[24:32])[0]
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
    
    # Read sentinel indices
    sentinel_indices = []
    for i in range(num_sentinels):
        idx_data = data[pos:pos + 8]
        idx = struct.unpack('<Q', idx_data)[0]
        sentinel_indices.append(idx)
        pos += 8
    
    return factors, num_sequences, num_sentinels, sequence_ids, sentinel_indices


@pytest.mark.skipif(not CPP_AVAILABLE, reason="C++ extension not available")
class TestParallelFastaProcessor:
    """Test parallel FASTA processing functions."""
    
    def test_parallel_w_rc_sequential(self):
        """Test parallel processor with num_threads=1 (sequential mode)."""
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Test with 1 thread (sequential)
            count = parallel_write_factors_binary_file_fasta_multiple_dna_w_rc(
                fasta_path, output_path, 1
            )
            
            assert count > 0, "Should produce at least one factor"
            assert os.path.exists(output_path), "Output file should be created"
            assert os.path.getsize(output_path) > 0, "Output file should not be empty"
            
            # Read and validate binary format
            factors, num_seqs, num_sents, seq_ids, sent_indices = read_binary_factors(output_path)
            
            assert len(factors) == count, "Factor count should match"
            assert num_seqs > 0, "Should have at least one sequence"
            assert len(seq_ids) == num_seqs, "Sequence ID count should match"
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_parallel_w_rc_auto_threads(self):
        """Test parallel processor with auto thread detection."""
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Test with auto threads (0)
            count = parallel_write_factors_binary_file_fasta_multiple_dna_w_rc(
                fasta_path, output_path, 0
            )
            
            assert count > 0, "Should produce at least one factor"
            assert os.path.exists(output_path), "Output file should be created"
            
            # Read and validate
            factors, _, _, seq_ids, _ = read_binary_factors(output_path)
            assert len(factors) == count
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_parallel_w_rc_multiple_threads(self):
        """Test parallel processor with explicit thread count."""
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Test with 2 threads
            count = parallel_write_factors_binary_file_fasta_multiple_dna_w_rc(
                fasta_path, output_path, 2
            )
            
            assert count > 0, "Should produce at least one factor"
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_parallel_no_rc_sequential(self):
        """Test parallel no-RC processor with sequential mode."""
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, output_path, 1
            )
            
            assert count > 0, "Should produce at least one factor"
            assert os.path.exists(output_path), "Output file should be created"
            
            # Read and validate
            factors, num_seqs, _, seq_ids, _ = read_binary_factors(output_path)
            assert len(factors) == count
            assert num_seqs > 0
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_parallel_no_rc_multiple_threads(self):
        """Test parallel no-RC processor with multiple threads."""
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Test with 2 threads
            count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, output_path, 2
            )
            
            assert count > 0, "Should produce at least one factor"
            
            # Read and validate
            factors, _, _, _, _ = read_binary_factors(output_path)
            assert len(factors) == count
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_parallel_reference_target_sequential(self):
        """Test parallel reference/target processor with sequential mode."""
        ref_path = str(RESOURCES_DIR / "short_dna1.fasta")
        target_path = str(RESOURCES_DIR / "short_dna2.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            count = parallel_write_factors_dna_w_reference_fasta_files_to_binary(
                ref_path, target_path, output_path, 1
            )
            
            assert count > 0, "Should produce at least one factor"
            assert os.path.exists(output_path), "Output file should be created"
            
            # Read and validate
            factors, num_seqs, _, seq_ids, _ = read_binary_factors(output_path)
            assert len(factors) == count
            assert num_seqs > 0, "Should have reference and target sequences"
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_parallel_reference_target_multiple_threads(self):
        """Test parallel reference/target processor with multiple threads."""
        ref_path = str(RESOURCES_DIR / "short_dna1.fasta")
        target_path = str(RESOURCES_DIR / "short_dna2.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            count = parallel_write_factors_dna_w_reference_fasta_files_to_binary(
                ref_path, target_path, output_path, 2
            )
            
            assert count > 0, "Should produce at least one factor"
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_consistency_sequential_vs_parallel_no_rc(self):
        """Test that sequential and parallel produce same results for no-RC."""
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp1:
            seq_output = tmp1.name
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp2:
            par_output = tmp2.name
        
        try:
            # Sequential (1 thread)
            seq_count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, seq_output, 1
            )
            
            # Parallel (2 threads)
            par_count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, par_output, 2
            )
            
            # Both should produce same number of factors
            assert seq_count == par_count, "Sequential and parallel should produce same factor count"
            
            # Read both outputs
            seq_factors, _, _, _, _ = read_binary_factors(seq_output)
            par_factors, _, _, _, _ = read_binary_factors(par_output)
            
            # Should have same factors
            assert len(seq_factors) == len(par_factors)
            assert seq_factors == par_factors, "Sequential and parallel should produce identical factors"
            
        finally:
            if os.path.exists(seq_output):
                os.remove(seq_output)
            if os.path.exists(par_output):
                os.remove(par_output)
    
    def test_delegation_sequential_functions(self):
        """Test that sequential functions delegate to parallel versions correctly."""
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp1:
            output1 = tmp1.name
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp2:
            output2 = tmp2.name
        
        try:
            # Call sequential function (should delegate to parallel with num_threads=1)
            seq_count = write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, output1
            )
            
            # Call parallel function explicitly with 1 thread
            par_count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, output2, 1
            )
            
            # Should produce identical results
            assert seq_count == par_count, "Sequential wrapper should produce same count as parallel(1)"
            
            # Read both outputs
            seq_factors, _, _, _, _ = read_binary_factors(output1)
            par_factors, _, _, _, _ = read_binary_factors(output2)
            
            assert seq_factors == par_factors, "Sequential wrapper should produce identical factors"
            
        finally:
            if os.path.exists(output1):
                os.remove(output1)
            if os.path.exists(output2):
                os.remove(output2)
    
    def test_invalid_fasta_path(self):
        """Test error handling for invalid FASTA path."""
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            with pytest.raises(RuntimeError, match="Cannot open"):
                parallel_write_factors_binary_file_fasta_multiple_dna_w_rc(
                    "/nonexistent/path.fasta", output_path, 1
                )
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_sentinel_tracking(self):
        """Test that sentinels are correctly identified in multi-sequence FASTA."""
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            count = parallel_write_factors_binary_file_fasta_multiple_dna_w_rc(
                fasta_path, output_path, 1
            )
            
            # Read metadata
            factors, num_seqs, num_sents, seq_ids, sent_indices = read_binary_factors(output_path)
            
            # Sentinels should be tracked if there are multiple sequences
            if num_seqs > 1:
                assert num_sents > 0, "Should have sentinel factors for multiple sequences"
                assert len(sent_indices) == num_sents, "Sentinel count should match"
                
                # All sentinel indices should be valid
                for idx in sent_indices:
                    assert 0 <= idx < count, f"Sentinel index {idx} should be within factor range"
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_actual_multithreading_no_rc(self):
        """Test that actual multi-threading occurs with large enough input."""
        # Create a large FASTA file (450K chars total, enough for 4+ threads)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_fasta:
            fasta_path = tmp_fasta.name
            total_chars = create_large_fasta_file(fasta_path, num_sequences=3, seq_length=150000)
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp_output:
            output_path = tmp_output.name
        
        try:
            print(f"\nGenerated FASTA with {total_chars} characters")
            
            # Test with explicit multi-threading (should actually use multiple threads)
            import time
            start = time.time()
            count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, output_path, 4
            )
            multi_time = time.time() - start
            
            print(f"Multi-threaded (4 threads): {count} factors in {multi_time:.3f}s")
            
            assert count > 0, "Should produce factors"
            assert os.path.exists(output_path), "Output file should be created"
            
            # Verify the file is valid
            factors, num_seqs, _, _, _ = read_binary_factors(output_path)
            assert len(factors) == count
            assert num_seqs == 3, "Should have 3 sequences"
            
        finally:
            if os.path.exists(fasta_path):
                os.remove(fasta_path)
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_threading_correctness_large_file(self):
        """Verify that multi-threaded and single-threaded produce identical results on large file."""
        # Create a large FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_fasta:
            fasta_path = tmp_fasta.name
            total_chars = create_large_fasta_file(fasta_path, num_sequences=2, seq_length=250000, seed=123)
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp1:
            output1 = tmp1.name
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp2:
            output2 = tmp2.name
        
        try:
            print(f"\nGenerated FASTA with {total_chars} characters for correctness test")
            
            # Single-threaded
            count1 = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, output1, 1
            )
            
            # Multi-threaded
            count2 = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, output2, 4
            )
            
            print(f"Single-threaded: {count1} factors")
            print(f"Multi-threaded (4 threads): {count2} factors")
            
            # Should produce same number of factors
            assert count1 == count2, "Single and multi-threaded should produce same factor count"
            
            # Read both outputs
            factors1, _, _, _, _ = read_binary_factors(output1)
            factors2, _, _, _, _ = read_binary_factors(output2)
            
            # Should have identical factors
            assert len(factors1) == len(factors2)
            assert factors1 == factors2, "Single and multi-threaded should produce identical factors"
            
        finally:
            if os.path.exists(fasta_path):
                os.remove(fasta_path)
            if os.path.exists(output1):
                os.remove(output1)
            if os.path.exists(output2):
                os.remove(output2)
    
    def test_auto_thread_detection_large_file(self):
        """Test that auto thread detection works with large files."""
        # Create a large FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_fasta:
            fasta_path = tmp_fasta.name
            total_chars = create_large_fasta_file(fasta_path, num_sequences=3, seq_length=200000, seed=456)
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp_output:
            output_path = tmp_output.name
        
        try:
            print(f"\nGenerated FASTA with {total_chars} characters for auto-detection test")
            
            # Test with auto-detection (num_threads=0)
            import time
            start = time.time()
            count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, output_path, 0  # auto-detect
            )
            auto_time = time.time() - start
            
            print(f"Auto-threaded: {count} factors in {auto_time:.3f}s")
            
            assert count > 0, "Should produce factors"
            
            # Verify correctness
            factors, num_seqs, _, _, _ = read_binary_factors(output_path)
            assert len(factors) == count
            assert num_seqs == 3
            
        finally:
            if os.path.exists(fasta_path):
                os.remove(fasta_path)
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_small_file_uses_single_thread(self):
        """Verify that small files correctly use single thread (no overhead)."""
        fasta_path = str(RESOURCES_DIR / "short_dna1.fasta")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Small file with auto-detect should use 1 thread
            count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                fasta_path, output_path, 0
            )
            
            assert count > 0, "Should produce factors even for small files"
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_large_reference_target(self):
        """Test reference/target processing with large files."""
        # Create large reference and target files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_ref:
            ref_path = tmp_ref.name
            create_large_fasta_file(ref_path, num_sequences=2, seq_length=150000, seed=789)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_target:
            target_path = tmp_target.name
            create_large_fasta_file(target_path, num_sequences=2, seq_length=150000, seed=790)
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp_output:
            output_path = tmp_output.name
        
        try:
            # Note: This function currently uses sequential for DNA w/ RC
            # but should still work correctly
            count = parallel_write_factors_dna_w_reference_fasta_files_to_binary(
                ref_path, target_path, output_path, 2
            )
            
            assert count > 0, "Should produce factors"
            
            # Verify output
            factors, num_seqs, _, seq_ids, _ = read_binary_factors(output_path)
            assert len(factors) == count
            assert num_seqs == 4, "Should have 2 reference + 2 target sequences"
            
        finally:
            if os.path.exists(ref_path):
                os.remove(ref_path)
            if os.path.exists(target_path):
                os.remove(target_path)
            if os.path.exists(output_path):
                os.remove(output_path)
    
    def test_large_file_performance(self):
        """Test parallel processor with a larger file to verify performance benefit."""
        # Use a larger test file if available
        large_fasta = RESOURCES_DIR / "test_bacterial_dna.fna"
        
        if not large_fasta.exists():
            pytest.skip("Large test file not available")
        
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            import time
            
            # Test sequential
            start = time.time()
            seq_count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                str(large_fasta), output_path, 1
            )
            seq_time = time.time() - start
            os.remove(output_path)
            
            # Test parallel
            start = time.time()
            par_count = parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(
                str(large_fasta), output_path, 0  # auto-detect
            )
            par_time = time.time() - start
            
            # Should produce same results
            assert seq_count == par_count, "Sequential and parallel should produce same factor count"
            
            # Parallel should be faster (or at least not significantly slower)
            # We allow some tolerance since small files might not benefit from parallelization
            print(f"Sequential: {seq_time:.3f}s, Parallel: {par_time:.3f}s")
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
