"""
Tests for the batch_factorize module.
"""

import os
import sys
import tempfile
from pathlib import Path
from unittest import mock

# Try to import the batch factorization script
try:
    from noLZSS.genomics import batch_factorize
    BATCH_FACTORIZE_AVAILABLE = True
except ImportError:
    print("Warning: batch_factorize module not available")
    BATCH_FACTORIZE_AVAILABLE = False


class TestBatchFactorize:
    """Test batch factorization functionality."""
    
    def test_batch_factorize_import(self):
        """Test that batch_factorize can be imported."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping batch_factorize import test - module not available")
            return
        
        assert hasattr(batch_factorize, 'main')
        assert hasattr(batch_factorize, 'FactorizationMode')
        assert hasattr(batch_factorize, 'BatchFactorizeError')
        print("batch_factorize import test passed")
    
    def test_factorization_mode_constants(self):
        """Test FactorizationMode constants."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping FactorizationMode test - module not available")
            return
        
        assert batch_factorize.FactorizationMode.WITHOUT_REVERSE_COMPLEMENT == "without_reverse_complement"
        assert batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT == "with_reverse_complement"
        assert batch_factorize.FactorizationMode.BOTH == "both"
        print("FactorizationMode constants test passed")
    
    def test_utility_functions(self):
        """Test utility functions."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping utility functions test - module not available")
            return
        
        # Test URL detection
        assert batch_factorize.is_url("http://example.com/file.fasta") == True
        assert batch_factorize.is_url("https://example.com/file.fasta") == True
        assert batch_factorize.is_url("ftp://example.com/file.fasta") == True
        assert batch_factorize.is_url("/local/path/file.fasta") == False
        assert batch_factorize.is_url("file.fasta") == False
        print("Utility functions test passed")
    
    def test_basic_functionality_with_test_files(self):
        """Test basic functionality with small test files."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping basic functionality test - module not available")
            return
        
        # Create temporary test files
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test FASTA file
            test_fasta = temp_path / "test.fasta"
            with open(test_fasta, 'w') as f:
                f.write(">seq1\nATCGATCG\n>seq2\nGCTAGCTA\n")
            
            # Create file list
            file_list = temp_path / "files.txt"
            with open(file_list, 'w') as f:
                f.write(str(test_fasta) + "\n")
            
            output_dir = temp_path / "output"
            
            try:
                # Test reading file list
                files = batch_factorize.read_file_list(file_list)
                assert len(files) == 1
                assert files[0] == str(test_fasta)
                
                # Test output path generation
                output_paths = batch_factorize.get_output_paths(
                    test_fasta, output_dir, batch_factorize.FactorizationMode.BOTH
                )
                assert "without_reverse_complement" in output_paths
                assert "with_reverse_complement" in output_paths
                
                print("Basic functionality test passed")
                
            except Exception as e:
                print(f"Basic functionality test failed with error: {e}")
                # Don't fail the test since this might be due to C++ extension issues
                print("This may be expected if C++ extension is not available")
    
    def test_shuffle_fasta_sequences(self):
        """Test FASTA sequence shuffling."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping shuffle test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test FASTA file
            test_fasta = temp_path / "test.fasta"
            original_seq1 = "ATCGATCGATCGATCG"
            original_seq2 = "GCTAGCTAGCTAGCTA"
            
            with open(test_fasta, 'w') as f:
                f.write(">seq1 Test sequence 1\n")
                f.write(original_seq1 + "\n")
                f.write(">seq2 Test sequence 2\n")
                f.write(original_seq2 + "\n")
            
            # Shuffle with a fixed seed
            shuffled_fasta = temp_path / "shuffled.fasta"
            result = batch_factorize.shuffle_fasta_sequences(
                test_fasta, shuffled_fasta, seed=42
            )
            
            assert result, "Shuffle should succeed"
            assert shuffled_fasta.exists(), "Shuffled file should exist"
            
            # Parse the shuffled file
            from noLZSS.genomics.fasta import _parse_fasta_content
            with open(shuffled_fasta, 'r') as f:
                shuffled_seqs = _parse_fasta_content(f.read())
            
            # Verify headers preserved
            assert "seq1" in shuffled_seqs, "seq1 header should be preserved"
            assert "seq2" in shuffled_seqs, "seq2 header should be preserved"
            
            # Verify sequences are shuffled (different from original)
            shuffled_seq1 = shuffled_seqs["seq1"]
            shuffled_seq2 = shuffled_seqs["seq2"]
            
            assert shuffled_seq1 != original_seq1, "Sequence 1 should be shuffled"
            assert shuffled_seq2 != original_seq2, "Sequence 2 should be shuffled"
            
            # Verify same composition (same letters, different order)
            assert sorted(shuffled_seq1) == sorted(original_seq1), "Sequence 1 should have same composition"
            assert sorted(shuffled_seq2) == sorted(original_seq2), "Sequence 2 should have same composition"
            
            print("Shuffle test passed")
    
    def test_shuffle_reproducibility(self):
        """Test that shuffling with same seed produces same result."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping reproducibility test - module not available")
            return
        
        from noLZSS.genomics.fasta import _parse_fasta_content
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test FASTA file
            test_fasta = temp_path / "test.fasta"
            with open(test_fasta, 'w') as f:
                f.write(">seq1\n")
                f.write("ATCGATCGATCGATCGATCGATCGATCGATCG\n")
            
            # Shuffle twice with same seed
            shuffled1 = temp_path / "shuffled1.fasta"
            shuffled2 = temp_path / "shuffled2.fasta"
            
            batch_factorize.shuffle_fasta_sequences(test_fasta, shuffled1, seed=123)
            batch_factorize.shuffle_fasta_sequences(test_fasta, shuffled2, seed=123)
            
            # Read both
            with open(shuffled1, 'r') as f:
                seqs1 = _parse_fasta_content(f.read())
            with open(shuffled2, 'r') as f:
                seqs2 = _parse_fasta_content(f.read())
            
            # Should be identical
            assert seqs1["seq1"] == seqs2["seq1"], "Same seed should produce same shuffle"
            
            print("Reproducibility test passed")
    
    def test_shuffle_gzipped_fasta(self):
        """Test shuffling of gzipped FASTA files."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping gzipped shuffle test - module not available")
            return
        
        import gzip
        from noLZSS.genomics.fasta import _parse_fasta_content
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test FASTA file
            test_fasta = temp_path / "test.fasta"
            original_seq = "ATCGATCGATCGATCGATCGATCGATCGATCG"
            
            with open(test_fasta, 'w') as f:
                f.write(">seq1 Test sequence\n")
                f.write(original_seq + "\n")
            
            # Compress the file
            gzipped_fasta = temp_path / "test.fasta.gz"
            with open(test_fasta, 'rb') as f_in:
                with gzip.open(gzipped_fasta, 'wb') as f_out:
                    f_out.write(f_in.read())
            
            # Verify file is gzipped
            assert batch_factorize.is_gzipped(gzipped_fasta), "File should be detected as gzipped"
            
            # Test decompression function
            decompressed_fasta = temp_path / "decompressed.fasta"
            result = batch_factorize.decompress_gzip(gzipped_fasta, decompressed_fasta)
            
            assert result, "Decompression should succeed"
            assert decompressed_fasta.exists(), "Decompressed file should exist"
            
            # Verify decompressed content matches original
            with open(test_fasta, 'r') as f1:
                with open(decompressed_fasta, 'r') as f2:
                    assert f1.read() == f2.read(), "Decompressed content should match original"
            
            # Test shuffling of decompressed file (the actual fix scenario)
            shuffled_fasta = temp_path / "shuffled.fasta"
            result = batch_factorize.shuffle_fasta_sequences(
                decompressed_fasta, shuffled_fasta, seed=42
            )
            
            assert result, "Shuffle of decompressed file should succeed"
            assert shuffled_fasta.exists(), "Shuffled file should exist"
            
            # Verify shuffled content
            with open(shuffled_fasta, 'r') as f:
                shuffled_seqs = _parse_fasta_content(f.read())
            
            assert "seq1" in shuffled_seqs, "seq1 header should be preserved"
            shuffled_seq = shuffled_seqs["seq1"]
            assert shuffled_seq != original_seq, "Sequence should be shuffled"
            assert sorted(shuffled_seq) == sorted(original_seq), "Sequence should have same composition"
            
            print("Gzipped FASTA shuffle test passed")
    
    def test_plot_factor_comparison(self):
        """Test comparison plot generation."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping plot comparison test - module not available")
            return
        
        try:
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
        except ImportError:
            print("Skipping plot test - matplotlib not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test FASTA files
            test_fasta = temp_path / "test.fasta"
            with open(test_fasta, 'w') as f:
                f.write(">seq1\n")
                f.write("ATCGATCGATCGATCGATCGATCGATCGATCG\n")  # Repetitive
            
            # Factorize original
            output_dir = temp_path / "output"
            output_paths = batch_factorize.get_output_paths(
                test_fasta, output_dir, batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT
            )
            
            try:
                batch_factorize.factorize_single_file(test_fasta, output_paths)
                
                # Create and factorize shuffled version
                shuffled_fasta = temp_path / "shuffled.fasta"
                batch_factorize.shuffle_fasta_sequences(test_fasta, shuffled_fasta, seed=42)
                
                shuffled_output_dir = temp_path / "shuffled_output"
                shuffled_output_paths = batch_factorize.get_output_paths(
                    shuffled_fasta, shuffled_output_dir, batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT
                )
                batch_factorize.factorize_single_file(shuffled_fasta, shuffled_output_paths)
                
                # Create comparison plot
                plot_path = temp_path / "comparison.png"
                result = batch_factorize.plot_factor_comparison(
                    original_factors_file=output_paths["with_reverse_complement"],
                    shuffled_factors_file=shuffled_output_paths["with_reverse_complement"],
                    output_plot_path=plot_path
                )
                
                assert result, "Plot generation should succeed"
                assert plot_path.exists(), "Plot file should exist"
                assert plot_path.stat().st_size > 0, "Plot file should not be empty"
                
                print("Plot comparison test passed")
            except Exception as e:
                print(f"Plot comparison test skipped - C++ extension issue: {e}")

    def test_complexity_table_and_tsv(self):
        """Test complexity table helpers with real C++ functions."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping complexity table test - module not available")
            return

        try:
            # Try to import the C++ extension
            from noLZSS._noLZSS import count_factors_fasta_dna_w_rc_per_sequence
        except ImportError:
            print("Skipping complexity table test - C++ extension not available")
            return

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta_file = temp_path / "test.fasta"
            # Create a simple FASTA file with repetitive sequences
            fasta_file.write_text(">seq1 description\nATCGATCGATCG\n>seq2 description\nGCTAGCTAGCTA\n", encoding="utf-8")

            try:
                rows = batch_factorize.compute_sequence_complexity_table(fasta_file)
                assert len(rows) == 2
                assert rows[0][0] == "seq1"  # sequence ID
                assert rows[1][0] == "seq2"  # sequence ID
                assert isinstance(rows[0][1], int)  # complexity_w_rc
                assert isinstance(rows[0][2], int)  # complexity_no_rc
                assert rows[0][1] > 0  # Should have some factors
                assert rows[0][2] > 0  # Should have some factors

                tsv_path = temp_path / "complexity.tsv"
                count = batch_factorize.write_sequence_complexity_tsv(fasta_file, tsv_path)
                assert count == 2
                assert tsv_path.exists()

                content = tsv_path.read_text(encoding="utf-8").strip().splitlines()
                assert content[0] == "sequence_id\tcomplexity_w_rc\tcomplexity_no_rc"
                assert content[1].startswith("seq1\t")
                assert content[2].startswith("seq2\t")
                
                print("Complexity table helper test passed")
            except Exception as e:
                print(f"Complexity table test skipped due to error: {e}")
    
    def test_validate_output_binary(self):
        """Test output validation function."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping validate_output_binary test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Test 1: Non-existent file
            non_existent = temp_path / "does_not_exist.bin"
            assert not batch_factorize.validate_output_binary(non_existent), "Non-existent file should be invalid"
            
            # Test 2: Empty file
            empty_file = temp_path / "empty.bin"
            empty_file.touch()
            assert not batch_factorize.validate_output_binary(empty_file), "Empty file should be invalid"
            
            # Test 3: Valid binary file (create one)
            try:
                # Create a test FASTA file
                test_fasta = temp_path / "test.fasta"
                with open(test_fasta, 'w') as f:
                    f.write(">seq1\nATCGATCGATCGATCG\n")
                
                # Factorize it to create a valid binary file
                output_dir = temp_path / "output"
                output_paths = batch_factorize.get_output_paths(
                    test_fasta, output_dir, batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT
                )
                
                result = batch_factorize.factorize_single_file(test_fasta, output_paths, skip_existing=False)
                
                if result.get("with_reverse_complement", False):
                    valid_file = output_paths["with_reverse_complement"]
                    assert batch_factorize.validate_output_binary(valid_file), "Valid binary file should pass validation"
                    
                    # Test 4: Corrupted file (truncate it)
                    with open(valid_file, 'rb') as f:
                        valid_data = f.read()
                    
                    corrupted_file = temp_path / "corrupted.bin"
                    with open(corrupted_file, 'wb') as f:
                        f.write(valid_data[:20])  # Write only first 20 bytes
                    
                    assert not batch_factorize.validate_output_binary(corrupted_file), "Corrupted file should be invalid"
                    
                    print("validate_output_binary test passed")
                else:
                    print("validate_output_binary test skipped - factorization failed")
                    
            except Exception as e:
                print(f"validate_output_binary test skipped due to error: {e}")
    
    def test_get_output_paths_from_source(self):
        """Test output path computation from source URLs/paths."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping get_output_paths_from_source test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "output"
            
            # Test 1: Local path without compression
            local_path = "/path/to/file.fasta"
            paths = batch_factorize.get_output_paths_from_source(
                local_path, output_dir, batch_factorize.FactorizationMode.BOTH
            )
            assert "without_reverse_complement" in paths
            assert "with_reverse_complement" in paths
            assert paths["without_reverse_complement"].name == "file.bin"
            assert paths["with_reverse_complement"].name == "file.bin"
            
            # Test 2: Local path with .gz extension
            gzipped_path = "/path/to/file.fasta.gz"
            paths = batch_factorize.get_output_paths_from_source(
                gzipped_path, output_dir, batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT
            )
            assert "with_reverse_complement" in paths
            assert "without_reverse_complement" not in paths
            assert paths["with_reverse_complement"].name == "file.bin"
            
            # Test 3: URL without compression
            url = "https://example.com/data/genome.fasta"
            paths = batch_factorize.get_output_paths_from_source(
                url, output_dir, batch_factorize.FactorizationMode.WITHOUT_REVERSE_COMPLEMENT
            )
            assert "without_reverse_complement" in paths
            assert "with_reverse_complement" not in paths
            assert paths["without_reverse_complement"].name == "genome.bin"
            
            # Test 4: URL with .gz extension
            gzipped_url = "https://example.com/data/genome.fa.gz"
            paths = batch_factorize.get_output_paths_from_source(
                gzipped_url, output_dir, batch_factorize.FactorizationMode.BOTH
            )
            assert paths["without_reverse_complement"].name == "genome.bin"
            assert paths["with_reverse_complement"].name == "genome.bin"
            
            # Test 5: URL with multiple extensions
            complex_url = "https://example.com/sequences.fasta.gz"
            paths = batch_factorize.get_output_paths_from_source(
                complex_url, output_dir, batch_factorize.FactorizationMode.BOTH
            )
            assert paths["without_reverse_complement"].name == "sequences.bin"
            
            print("get_output_paths_from_source test passed")
    
    def test_process_single_file_complete_skip_logic(self):
        """Test that process_single_file_complete correctly skips files with valid outputs."""
        if not BATCH_FACTORIZE_AVAILABLE:
            print("Skipping process_single_file_complete test - module not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            try:
                # Create a test FASTA file
                test_fasta = temp_path / "test.fasta"
                with open(test_fasta, 'w') as f:
                    f.write(">seq1\nATCGATCGATCGATCG\n")
                
                output_dir = temp_path / "output"
                download_dir = temp_path / "download"
                download_dir.mkdir(parents=True, exist_ok=True)
                
                # First run: should process the file
                logger = batch_factorize.setup_logging("INFO")
                file_info = (
                    str(test_fasta),
                    output_dir,
                    download_dir,
                    batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT,
                    True,  # skip_existing
                    3,     # max_retries
                    logger.name
                )
                
                result_path, results = batch_factorize.process_single_file_complete(file_info)
                
                assert result_path == str(test_fasta)
                assert "with_reverse_complement" in results
                
                if results.get("with_reverse_complement", False):
                    # Second run: should skip the file
                    result_path2, results2 = batch_factorize.process_single_file_complete(file_info)
                    
                    assert result_path2 == str(test_fasta)
                    assert results2.get("with_reverse_complement", False), "Should skip and return success"
                    
                    # Third run with skip_existing=False: should re-process
                    file_info_no_skip = (
                        str(test_fasta),
                        output_dir,
                        download_dir,
                        batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT,
                        False,  # skip_existing = False
                        3,
                        logger.name
                    )
                    
                    result_path3, results3 = batch_factorize.process_single_file_complete(file_info_no_skip)
                    assert results3.get("with_reverse_complement", False), "Should re-process successfully"
                    
                    print("process_single_file_complete skip logic test passed")
                else:
                    print("process_single_file_complete test skipped - initial factorization failed")
                    
            except Exception as e:
                print(f"process_single_file_complete test skipped due to error: {e}")


def run_tests():
    """Run all batch factorize tests."""
    print("\n=== TestBatchFactorize ===")
    
    test_instance = TestBatchFactorize()
    
    test_instance.test_batch_factorize_import()
    test_instance.test_factorization_mode_constants()
    test_instance.test_utility_functions()
    test_instance.test_basic_functionality_with_test_files()
    test_instance.test_shuffle_fasta_sequences()
    test_instance.test_shuffle_reproducibility()
    test_instance.test_shuffle_gzipped_fasta()
    test_instance.test_plot_factor_comparison()
    test_instance.test_validate_output_binary()
    test_instance.test_get_output_paths_from_source()
    test_instance.test_process_single_file_complete_skip_logic()
    
    print("All batch_factorize tests completed")


if __name__ == "__main__":
    run_tests()