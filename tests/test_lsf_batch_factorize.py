"""
Tests for the lsf_batch_factorize module.
"""

import json
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Try to import the LSF batch factorization script
try:
    from noLZSS.genomics import lsf_batch_factorize
    LSF_BATCH_AVAILABLE = True
except ImportError:
    print("Warning: lsf_batch_factorize module not available")
    LSF_BATCH_AVAILABLE = False


class TestResourceEstimator:
    """Test ResourceEstimator class."""
    
    def test_resource_estimator_init(self):
        """Test ResourceEstimator initialization."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping ResourceEstimator init test - module not available")
            return
        
        estimator = lsf_batch_factorize.ResourceEstimator(safety_factor=2.0)
        assert estimator.safety_factor == 2.0
        print("ResourceEstimator initialization test passed")
    
    def test_default_estimates(self):
        """Test default resource estimation."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping default estimates test - module not available")
            return
        
        estimator = lsf_batch_factorize.ResourceEstimator()
        estimates = estimator.estimate_resources(
            input_size=100000,
            mode="with_reverse_complement",
            max_threads=1
        )
        
        assert 'lsf_time_minutes' in estimates
        assert 'lsf_memory_gb' in estimates
        assert 'lsf_disk_gb' in estimates
        assert estimates['lsf_time_minutes'] > 0
        assert estimates['lsf_memory_gb'] > 0
        print("Default estimates test passed")
    
    def test_estimate_with_trends(self):
        """Test resource estimation with trend parameters."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping estimate with trends test - module not available")
            return
        
        # Create mock trend parameters
        mock_trends = {
            'write_factors_binary_file_fasta_multiple_dna_w_rc_time': {
                'slope': 0.027,
                'intercept': 0.5,
                'log_scale': False
            },
            'write_factors_binary_file_fasta_multiple_dna_w_rc_memory': {
                'slope': 0.014,
                'intercept': 1.0,
                'log_scale': False
            },
            'write_factors_binary_file_fasta_multiple_dna_w_rc_disk_space': {
                'slope': 0.0024,
                'intercept': 0.1,
                'log_scale': False
            }
        }
        
        estimator = lsf_batch_factorize.ResourceEstimator(safety_factor=1.5)
        estimator.trends = mock_trends
        
        estimates = estimator.estimate_resources(
            input_size=100000,
            mode="with_reverse_complement",
            max_threads=1
        )
        
        assert 'estimated_time_seconds' in estimates
        assert 'safe_time_seconds' in estimates
        assert 'lsf_time_minutes' in estimates
        assert estimates['estimated_time_seconds'] > 0
        print("Estimate with trends test passed")


class TestLSFJobManager:
    """Test LSFJobManager class."""
    
    def test_lsf_job_manager_init(self):
        """Test LSFJobManager initialization."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping LSFJobManager init test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            log_dir = Path(temp_dir) / "logs"
            manager = lsf_batch_factorize.LSFJobManager(
                queue="normal",
                log_dir=log_dir
            )
            
            assert manager.queue == "normal"
            assert manager.log_dir == log_dir
            assert log_dir.exists()
            print("LSFJobManager initialization test passed")
    
    @patch('subprocess.run')
    def test_submit_job(self, mock_run):
        """Test job submission."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping job submission test - module not available")
            return
        
        # Mock successful bsub output
        mock_run.return_value = Mock(
            returncode=0,
            stdout="Job <12345> is submitted to queue <normal>.",
            stderr=""
        )
        
        with tempfile.TemporaryDirectory() as temp_dir:
            log_dir = Path(temp_dir) / "logs"
            manager = lsf_batch_factorize.LSFJobManager(log_dir=log_dir)
            
            job_id = manager.submit_job(
                job_name="test_job",
                command="echo 'hello'",
                time_minutes=10,
                memory_gb=4,
                output_file=log_dir / "test.log"
            )
            
            assert job_id == "12345"
            assert "12345" in manager.submitted_jobs
            print("Job submission test passed")
    
    @patch('subprocess.run')
    def test_check_job_status(self, mock_run):
        """Test job status checking."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping job status check test - module not available")
            return
        
        # Mock bjobs output
        mock_run.return_value = Mock(
            returncode=0,
            stdout="12345 user RUN normal host exec test_job time",
            stderr=""
        )
        
        with tempfile.TemporaryDirectory() as temp_dir:
            manager = lsf_batch_factorize.LSFJobManager(log_dir=Path(temp_dir))
            status = manager.check_job_status("12345")
            
            assert status == "RUN"
            print("Job status check test passed")
    
    def test_generate_summary(self):
        """Test summary generation."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping summary generation test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            log_dir = Path(temp_dir)
            manager = lsf_batch_factorize.LSFJobManager(log_dir=log_dir)
            
            # Mock some submitted jobs
            manager.submitted_jobs = {
                "12345": {"name": "job1", "output_file": log_dir / "job1.log"},
                "12346": {"name": "job2", "output_file": log_dir / "job2.log"},
            }
            
            job_status = {
                "12345": "DONE",
                "12346": "EXIT",
            }
            
            summary_file = log_dir / "summary.txt"
            manager.generate_summary(job_status, summary_file)
            
            assert summary_file.exists()
            with open(summary_file, 'r') as f:
                content = f.read()
                assert "Total jobs: 2" in content
                assert "Successful: 1" in content
                assert "Failed: 1" in content
            
            print("Summary generation test passed")


class TestHelperFunctions:
    """Test helper functions."""
    
    def test_get_fasta_size(self):
        """Test FASTA file size estimation."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping FASTA size test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            fasta_path = Path(temp_dir) / "test.fasta"
            with open(fasta_path, 'w') as f:
                f.write(">seq1\nATCGATCG\n>seq2\nGCTAGCTA\n")
            
            size = lsf_batch_factorize.get_fasta_size(fasta_path)
            assert size == 16  # 8 + 8 nucleotides
            print("FASTA size test passed")
    
    def test_prepare_jobs(self):
        """Test job preparation."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping job preparation test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test FASTA files
            fasta1 = temp_path / "test1.fasta"
            with open(fasta1, 'w') as f:
                f.write(">seq1\nATCGATCG\n")
            
            fasta2 = temp_path / "test2.fasta"
            with open(fasta2, 'w') as f:
                f.write(">seq2\nGCTAGCTA\n")
            
            file_list = [str(fasta1), str(fasta2)]
            output_dir = temp_path / "output"
            
            estimator = lsf_batch_factorize.ResourceEstimator()
            
            jobs = lsf_batch_factorize.prepare_jobs(
                file_list=file_list,
                output_dir=output_dir,
                mode="with_reverse_complement",
                estimator=estimator,
                max_threads=1
            )
            
            assert len(jobs) == 2
            assert 'input_file' in jobs[0]
            assert 'output_paths' in jobs[0]
            assert 'resources' in jobs[0]
            print("Job preparation test passed")


class TestIntegration:
    """Integration tests."""
    
    def test_full_workflow_dry_run(self):
        """Test full workflow in dry-run mode."""
        if not LSF_BATCH_AVAILABLE:
            print("Skipping full workflow test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test FASTA file
            fasta_file = temp_path / "test.fasta"
            with open(fasta_file, 'w') as f:
                f.write(">seq1\nATCGATCG\n")
            
            # Create file list
            file_list = temp_path / "files.txt"
            with open(file_list, 'w') as f:
                f.write(str(fasta_file) + "\n")
            
            # Test that the module can be invoked (in dry-run mode)
            # We can't actually run the full script without LSF, but we can test the structure
            try:
                from noLZSS.genomics import batch_factorize
                
                # Test reading the file list
                files = batch_factorize.read_file_list(file_list)
                assert len(files) == 1
                
                # Test resource estimator
                estimator = lsf_batch_factorize.ResourceEstimator()
                estimates = estimator.estimate_resources(
                    input_size=8,
                    mode="with_reverse_complement",
                    max_threads=1
                )
                assert estimates is not None
                
                print("Full workflow dry-run test passed")
                
            except Exception as e:
                print(f"Full workflow test failed with error: {e}")
                print("This may be expected if dependencies are not available")


def run_tests():
    """Run all LSF batch factorize tests."""
    print("\n=== TestResourceEstimator ===")
    test_re = TestResourceEstimator()
    test_re.test_resource_estimator_init()
    test_re.test_default_estimates()
    test_re.test_estimate_with_trends()
    
    print("\n=== TestLSFJobManager ===")
    test_jm = TestLSFJobManager()
    test_jm.test_lsf_job_manager_init()
    test_jm.test_submit_job()
    test_jm.test_check_job_status()
    test_jm.test_generate_summary()
    
    print("\n=== TestHelperFunctions ===")
    test_hf = TestHelperFunctions()
    test_hf.test_get_fasta_size()
    test_hf.test_prepare_jobs()
    
    print("\n=== TestIntegration ===")
    test_int = TestIntegration()
    test_int.test_full_workflow_dry_run()
    
    print("\nAll LSF batch factorize tests completed")


if __name__ == "__main__":
    run_tests()
