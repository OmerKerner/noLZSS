"""
Tests for the lsf_batch_factorize module.

Note: These tests do not require actual LSF cluster access.
They test the helper functions and resource estimation logic.
"""

import os
import sys
import tempfile
import json
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Try to import the LSF batch factorization script
try:
    from noLZSS.genomics import lsf_batch_factorize
    LSF_BATCH_FACTORIZE_AVAILABLE = True
except ImportError:
    print("Warning: lsf_batch_factorize module not available")
    LSF_BATCH_FACTORIZE_AVAILABLE = False


class TestLSFBatchFactorize:
    """Test LSF batch factorization functionality."""
    
    def test_lsf_batch_factorize_import(self):
        """Test that lsf_batch_factorize can be imported."""
        if not LSF_BATCH_FACTORIZE_AVAILABLE:
            print("Skipping lsf_batch_factorize import test - module not available")
            return
        
        assert hasattr(lsf_batch_factorize, 'main')
        assert hasattr(lsf_batch_factorize, 'LSFBatchFactorizeError')
        assert hasattr(lsf_batch_factorize, 'estimate_fasta_nucleotides')
        assert hasattr(lsf_batch_factorize, 'decide_num_threads')
        print("lsf_batch_factorize import test passed")
    
    def test_estimate_fasta_nucleotides(self):
        """Test nucleotide estimation from FASTA file."""
        if not LSF_BATCH_FACTORIZE_AVAILABLE:
            print("Skipping estimate_fasta_nucleotides test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_fasta = temp_path / "test.fasta"
            
            # Create a test FASTA file with known content
            fasta_content = ">seq1\nATCGATCGATCGATCG\n>seq2\nGCTAGCTAGCTAGCTA\n"
            with open(test_fasta, 'w') as f:
                f.write(fasta_content)
            
            # Estimate nucleotides
            estimated = lsf_batch_factorize.estimate_fasta_nucleotides(test_fasta)
            
            # Should estimate around 32 nucleotides (16 + 16)
            # Allow some margin for header estimation
            assert estimated > 20, f"Estimated nucleotides too low: {estimated}"
            assert estimated < 100, f"Estimated nucleotides too high: {estimated}"
            print(f"Nucleotide estimation test passed: estimated={estimated}")
    
    def test_decide_num_threads(self):
        """Test thread count decision logic."""
        if not LSF_BATCH_FACTORIZE_AVAILABLE:
            print("Skipping decide_num_threads test - module not available")
            return
        
        # Test small file
        threads = lsf_batch_factorize.decide_num_threads(50_000, max_threads=16)
        assert threads == 1, f"Small file should use 1 thread, got {threads}"
        
        # Test medium file
        threads = lsf_batch_factorize.decide_num_threads(500_000, max_threads=16)
        assert threads <= 4, f"Medium file should use <=4 threads, got {threads}"
        
        # Test large file
        threads = lsf_batch_factorize.decide_num_threads(5_000_000, max_threads=16)
        assert threads <= 8, f"Large file should use <=8 threads, got {threads}"
        
        # Test very large file
        threads = lsf_batch_factorize.decide_num_threads(50_000_000, max_threads=16)
        assert threads == 16, f"Very large file should use max threads, got {threads}"
        
        # Test with low max_threads
        threads = lsf_batch_factorize.decide_num_threads(50_000_000, max_threads=4)
        assert threads == 4, f"Should respect max_threads, got {threads}"
        
        print("Thread count decision test passed")
    
    def test_estimate_resources_fallback(self):
        """Test fallback resource estimation."""
        if not LSF_BATCH_FACTORIZE_AVAILABLE:
            print("Skipping estimate_resources_fallback test - module not available")
            return
        
        # Test with 1 MB file
        file_size = 1_000_000
        estimate = lsf_batch_factorize.estimate_resources_fallback(
            file_size, num_threads=4, safety_factor=1.5
        )
        
        # Check required fields
        assert 'input_size_nucleotides' in estimate
        assert 'num_threads' in estimate
        assert 'estimated_time_seconds' in estimate
        assert 'estimated_memory_gb' in estimate
        assert 'cluster_memory_gb' in estimate
        assert 'safe_time_minutes' in estimate
        
        # Check reasonable values
        assert estimate['num_threads'] == 4
        assert estimate['safety_factor'] == 1.5
        assert estimate['estimated_time_seconds'] > 0
        assert estimate['estimated_memory_gb'] > 0
        assert estimate['cluster_memory_gb'] in [2, 4, 8, 16, 32, 64, 128]
        
        print("Fallback resource estimation test passed")
    
    def test_estimate_resources_from_trends(self):
        """Test resource estimation from benchmark trends."""
        if not LSF_BATCH_FACTORIZE_AVAILABLE:
            print("Skipping estimate_resources_from_trends test - module not available")
            return
        
        # Create mock trends
        trends = {
            'write_factors_binary_file_fasta_multiple_dna_w_rc_time': {
                'slope': 0.001,
                'intercept': 100,
                'log_scale': False
            },
            'write_factors_binary_file_fasta_multiple_dna_w_rc_memory': {
                'slope': 0.00003,  # 30 bytes per nucleotide
                'intercept': 10,
                'log_scale': False
            },
            'write_factors_binary_file_fasta_multiple_dna_w_rc_disk_space': {
                'slope': 0.00001,
                'intercept': 5,
                'log_scale': False
            }
        }
        
        from noLZSS.genomics.batch_factorize import FactorizationMode
        
        estimate = lsf_batch_factorize.estimate_resources_from_trends(
            file_size_nucleotides=1_000_000,
            trends=trends,
            mode=FactorizationMode.WITH_REVERSE_COMPLEMENT,
            num_threads=4,
            safety_factor=1.5
        )
        
        # Check required fields
        assert 'estimated_time_seconds' in estimate
        assert 'estimated_memory_gb' in estimate
        assert 'cluster_memory_gb' in estimate
        assert 'estimated_disk_gb' in estimate
        
        # Check that values are reasonable
        assert estimate['estimated_time_seconds'] > 0
        assert estimate['estimated_memory_gb'] > 0
        assert estimate['cluster_memory_gb'] in [2, 4, 8, 16, 32, 64, 128]
        
        print("Trend-based resource estimation test passed")
    
    def test_create_job_script(self):
        """Test job script creation."""
        if not LSF_BATCH_FACTORIZE_AVAILABLE:
            print("Skipping create_job_script test - module not available")
            return
        
        from noLZSS.genomics.batch_factorize import FactorizationMode
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            input_file = temp_path / "input.fasta"
            output_file = temp_path / "output.bin"
            script_path = temp_path / "job.sh"
            
            # Create dummy input file
            input_file.touch()
            
            # Create job script
            success = lsf_batch_factorize.create_job_script(
                input_file=input_file,
                output_file=output_file,
                mode=FactorizationMode.WITH_REVERSE_COMPLEMENT,
                num_threads=4,
                script_path=script_path
            )
            
            assert success, "Job script creation failed"
            assert script_path.exists(), "Job script not created"
            
            # Check script content
            with open(script_path, 'r') as f:
                script_content = f.read()
            
            assert '#!/bin/bash' in script_content
            assert 'parallel_write_factors_binary_file_fasta_multiple_dna_w_rc' in script_content
            assert str(input_file) in script_content
            assert str(output_file) in script_content
            assert 'num_threads = 4' in script_content
            
            # Check script is executable
            assert os.access(script_path, os.X_OK), "Script not executable"
            
            print("Job script creation test passed")
    
    def test_load_benchmark_trends_not_found(self):
        """Test loading benchmark trends when file not found."""
        if not LSF_BATCH_FACTORIZE_AVAILABLE:
            print("Skipping load_benchmark_trends test - module not available")
            return
        
        # Try to load non-existent file
        trends = lsf_batch_factorize.load_benchmark_trends(
            Path("/nonexistent/trends.pkl")
        )
        
        assert trends is None, "Should return None for non-existent file"
        print("Benchmark trends loading test passed")
    
    def test_check_job_output(self):
        """Test job output checking."""
        if not LSF_BATCH_FACTORIZE_AVAILABLE:
            print("Skipping check_job_output test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            output_file = temp_path / "output.bin"
            error_log = temp_path / "error.log"
            
            # Test: output file doesn't exist
            success, msg = lsf_batch_factorize.check_job_output(output_file, error_log)
            assert not success, "Should fail when output doesn't exist"
            assert "not created" in msg.lower()
            
            # Test: output file is empty
            output_file.touch()
            success, msg = lsf_batch_factorize.check_job_output(output_file, error_log)
            assert not success, "Should fail when output is empty"
            assert "empty" in msg.lower()
            
            # Test: output file has content, no error log
            with open(output_file, 'wb') as f:
                f.write(b'some binary data')
            success, msg = lsf_batch_factorize.check_job_output(output_file, error_log)
            assert success, "Should succeed when output has content"
            assert msg is None
            
            # Test: output file has content, error log with errors
            with open(error_log, 'w') as f:
                f.write("Error: something went wrong\n")
            success, msg = lsf_batch_factorize.check_job_output(output_file, error_log)
            assert not success, "Should fail when error log has errors"
            assert "error" in msg.lower()
            
            print("Job output checking test passed")
    
    def test_get_file_size(self):
        """Test file size retrieval."""
        if not LSF_BATCH_FACTORIZE_AVAILABLE:
            print("Skipping get_file_size test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_file = temp_path / "test.txt"
            
            # Create file with known size
            content = "A" * 1000
            with open(test_file, 'w') as f:
                f.write(content)
            
            size = lsf_batch_factorize.get_file_size(test_file)
            assert size == 1000, f"Expected size 1000, got {size}"
            
            print("File size test passed")


def run_all_tests():
    """Run all tests."""
    test_class = TestLSFBatchFactorize()
    
    tests = [
        test_class.test_lsf_batch_factorize_import,
        test_class.test_estimate_fasta_nucleotides,
        test_class.test_decide_num_threads,
        test_class.test_estimate_resources_fallback,
        test_class.test_estimate_resources_from_trends,
        test_class.test_create_job_script,
        test_class.test_load_benchmark_trends_not_found,
        test_class.test_check_job_output,
        test_class.test_get_file_size,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        test_name = test.__name__
        try:
            print(f"\nRunning {test_name}...")
            test()
            passed += 1
        except Exception as e:
            print(f"FAILED: {test_name}")
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*60)
    print(f"Test Summary: {passed} passed, {failed} failed")
    print("="*60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
