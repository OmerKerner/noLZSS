"""
Tests for the genomics subpackage.
"""

import sys
import os


class TestGenomicsStructure:
    """Test the genomics subpackage structure."""
    
    def test_genomics_import(self):
        """Test that genomics subpackage can be imported."""
        try:
            from noLZSS import genomics
            print("Genomics subpackage imports successfully")
        except ImportError as e:
            print(f"Warning: Genomics import failed: {e}")
            assert False, f"Genomics import failed: {e}"
    
    def test_fasta_module_import(self):
        """Test that fasta module can be imported."""
        try:
            from noLZSS.genomics import fasta
            print("FASTA module imports successfully")
        except ImportError as e:
            print(f"Warning: FASTA module import failed: {e}")
            assert False, f"FASTA module import failed: {e}"
    
    def test_sequences_module_import(self):
        """Test that sequences module can be imported."""
        try:
            from noLZSS.genomics import sequences
            print("Sequences module imports successfully")
        except ImportError as e:
            print(f"Warning: Sequences module import failed: {e}")
            assert False, f"Sequences module import failed: {e}"
    
    def test_genomics_namespace(self):
        """Test the genomics namespace structure."""
        try:
            import noLZSS.genomics as genomics
            
            # Check that submodules are accessible
            assert hasattr(genomics, 'fasta')
            assert hasattr(genomics, 'sequences')
            
            print("Genomics namespace structure is correct")
        except Exception as e:
            print(f"Warning: Genomics namespace error: {e}")
            assert False, f"Genomics namespace error: {e}"


class TestFutureGenomicsFeatures:
    """Placeholder tests for future genomics functionality."""
    
    def test_fasta_placeholder(self):
        """Test that FASTA module exists and is ready for implementation."""
        try:
            from noLZSS.genomics import fasta
            
            # Check that it's a module
            assert hasattr(fasta, '__file__')
            
            print("FASTA module ready for implementation")
        except Exception as e:
            print(f"Warning: FASTA module not ready: {e}")
            assert False, f"FASTA module not ready: {e}"
    
    def test_sequences_placeholder(self):
        """Test that sequences module exists and is ready for implementation."""
        try:
            from noLZSS.genomics import sequences
            
            # Check that it's a module
            assert hasattr(sequences, '__file__')
            
            print("Sequences module ready for implementation")
        except Exception as e:
            print(f"Warning: Sequences module not ready: {e}")
            assert False, f"Sequences module not ready: {e}"


import tempfile
import os
from pathlib import Path
import pytest

from noLZSS.genomics.fasta import read_nucleotide_fasta, read_protein_fasta, read_fasta_auto, FASTAError
from noLZSS.genomics.plots import plot_single_seq_accum_factors_from_file, PlotError
from noLZSS._noLZSS import process_nucleotide_fasta, process_amino_acid_fasta  # C++ functions
from noLZSS._noLZSS import process_nucleotide_fasta  # C++ function


class TestFASTAFunctions:
    """Test FASTA file parsing and processing functions."""
    
    def test_read_nucleotide_fasta_valid(self):
        """Test reading valid nucleotide FASTA file."""
        fasta_content = """>seq1
ATCGATCG
>seq2
GCTAGCTA
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            results = read_nucleotide_fasta(temp_path)
            
            assert len(results) == 2
            assert results[0][0] == 'seq1'
            assert results[1][0] == 'seq2'
            
            # Check that factors are returned (we can't check exact values without the C++ extension)
            assert isinstance(results[0][1], list)
            assert isinstance(results[1][1], list)
            
            print("Nucleotide FASTA reading works correctly")
        except Exception as e:
            print(f"Warning: Nucleotide FASTA test failed: {e}")
            assert False, f"Nucleotide FASTA test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_read_nucleotide_fasta_invalid_nucleotides(self):
        """Test reading FASTA with invalid nucleotides."""
        fasta_content = """>seq1
ATCGNTCG
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            with pytest.raises(FASTAError):
                read_nucleotide_fasta(temp_path)
            print("Invalid nucleotide detection works correctly")
        except Exception as e:
            print(f"Warning: Invalid nucleotide test failed: {e}")
            assert False, f"Invalid nucleotide test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_read_protein_fasta_valid(self):
        """Test reading valid protein FASTA file."""
        fasta_content = """>protein1
MKVLWAALL
>protein2
ACDEFGHIK
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            results = read_protein_fasta(temp_path)
            
            assert len(results) == 2
            assert results[0][0] == 'protein1'
            assert results[0][1] == 'MKVLWAALL'
            assert results[1][0] == 'protein2'
            assert results[1][1] == 'ACDEFGHIK'
            
            print("Protein FASTA reading works correctly")
        except Exception as e:
            print(f"Warning: Protein FASTA test failed: {e}")
            assert False, f"Protein FASTA test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_read_protein_fasta_invalid_amino_acids(self):
        """Test reading FASTA with invalid amino acids."""
        fasta_content = """>protein1
MKVLWAALL1
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            with pytest.raises(FASTAError):
                read_protein_fasta(temp_path)
            print("Invalid amino acid detection works correctly")
        except Exception as e:
            print(f"Warning: Invalid amino acid test failed: {e}")
            assert False, f"Invalid amino acid test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_read_fasta_auto_nucleotide(self):
        """Test auto-detection of nucleotide FASTA."""
        fasta_content = """>seq1
ATCGATCG
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            results = read_fasta_auto(temp_path)
            
            # Should return factorized results for nucleotides
            assert len(results) == 1
            assert results[0][0] == 'seq1'
            assert isinstance(results[0][1], list)  # factors
            
            print("Auto-detection of nucleotide FASTA works correctly")
        except Exception as e:
            print(f"Warning: Auto-detection nucleotide test failed: {e}")
            assert False, f"Auto-detection nucleotide test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_read_fasta_auto_protein(self):
        """Test auto-detection of protein FASTA."""
        fasta_content = """>protein1
MKVLWAALL
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            results = read_fasta_auto(temp_path)
            
            # Should return sequence strings for proteins
            assert len(results) == 1
            assert results[0][0] == 'protein1'
            assert results[0][1] == 'MKVLWAALL'
            
            print("Auto-detection of protein FASTA works correctly")
        except Exception as e:
            print(f"Warning: Auto-detection protein test failed: {e}")
            assert False, f"Auto-detection protein test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_read_fasta_auto_invalid_type(self):
        """Test auto-detection with invalid sequence type."""
        fasta_content = """>seq1
123456789
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            with pytest.raises(FASTAError):
                read_fasta_auto(temp_path)
            print("Invalid sequence type detection works correctly")
        except Exception as e:
            print(f"Warning: Invalid type test failed: {e}")
            assert False, f"Invalid type test failed: {e}"
        finally:
            os.unlink(temp_path)


class TestCppFastaFunctions:
    """Test C++ FASTA processing functions for memory efficiency."""

    def test_process_nucleotide_fasta_valid(self):
        """Test C++ FASTA processing with valid nucleotide sequences."""
        fasta_content = """>seq1
ATCGATCG
>seq2
GCTAGCTA
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            # Try to import C++ function
            try:
                from noLZSS._noLZSS import process_nucleotide_fasta as cpp_process_fasta
                result = cpp_process_fasta(temp_path)
                
                # Verify result structure
                assert isinstance(result, dict)
                assert "sequence" in result
                assert "num_sequences" in result
                assert "sequence_ids" in result
                assert "sequence_lengths" in result
                assert "sequence_positions" in result
                
                # Verify content
                assert result["num_sequences"] == 2
                assert result["sequence_ids"] == ["seq1", "seq2"]
                assert result["sequence_lengths"] == [8, 8]
                assert len(result["sequence"]) == 16 + 1  # 8 + 8 + 1 sentinel
                
                print("C++ FASTA processing works correctly")
            except ImportError:
                print("C++ extension not available, skipping test")
                return  # Skip test if C++ extension not built
        except Exception as e:
            print(f"Warning: C++ FASTA test failed: {e}")
            assert False, f"C++ FASTA test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_process_nucleotide_fasta_invalid_nucleotides(self):
        """Test C++ FASTA processing with invalid nucleotides."""
        fasta_content = """>seq1
ATCGNTCG
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            try:
                from noLZSS._noLZSS import process_nucleotide_fasta as cpp_process_fasta
                with pytest.raises(RuntimeError):  # C++ throws runtime_error
                    cpp_process_fasta(temp_path)
                print("C++ FASTA invalid nucleotide detection works correctly")
            except ImportError:
                print("C++ extension not available, skipping test")
                return
        except Exception as e:
            print(f"Warning: C++ invalid nucleotide test failed: {e}")
            assert False, f"C++ invalid nucleotide test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_process_nucleotide_fasta_empty_file(self):
        """Test C++ FASTA processing with empty file."""
        fasta_content = ""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            try:
                from noLZSS._noLZSS import process_nucleotide_fasta as cpp_process_fasta
                with pytest.raises(RuntimeError):
                    cpp_process_fasta(temp_path)
                print("C++ FASTA empty file handling works correctly")
            except ImportError:
                print("C++ extension not available, skipping test")
                return
        except Exception as e:
            print(f"Warning: C++ empty file test failed: {e}")
            assert False, f"C++ empty file test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_process_amino_acid_fasta_valid(self):
        """Test C++ amino acid FASTA processing with valid sequences."""
        fasta_content = """>protein1
ACDEFGHIK
>protein2
LMNPQRSTV
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            # Try to import C++ function
            try:
                from noLZSS._noLZSS import process_amino_acid_fasta as cpp_process_aa_fasta
                result = cpp_process_aa_fasta(temp_path)
                
                # Verify result structure
                assert isinstance(result, dict)
                assert "sequence" in result
                assert "num_sequences" in result
                assert "sequence_ids" in result
                assert "sequence_lengths" in result
                assert "sequence_positions" in result
                
                # Verify content
                assert result["num_sequences"] == 2
                assert result["sequence_ids"] == ["protein1", "protein2"]
                assert result["sequence_lengths"] == [9, 9]
                assert len(result["sequence"]) == 18 + 1  # 9 + 9 + 1 sentinel
                
                print("C++ amino acid FASTA processing works correctly")
            except ImportError:
                print("C++ extension not available, skipping test")
                return  # Skip test if C++ extension not built
        except Exception as e:
            print(f"Warning: C++ amino acid FASTA test failed: {e}")
            assert False, f"C++ amino acid FASTA test failed: {e}"
        finally:
            os.unlink(temp_path)
    
    def test_process_amino_acid_fasta_invalid_amino_acids(self):
        """Test C++ amino acid FASTA processing with invalid amino acids."""
        fasta_content = """>protein1
ACDEFGHIK1
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            try:
                from noLZSS._noLZSS import process_amino_acid_fasta as cpp_process_aa_fasta
                with pytest.raises(RuntimeError):  # C++ throws runtime_error
                    cpp_process_aa_fasta(temp_path)
                print("C++ amino acid FASTA invalid character detection works correctly")
            except ImportError:
                print("C++ extension not available, skipping test")
                return
        except Exception as e:
            print(f"Warning: C++ invalid amino acid test failed: {e}")
            assert False, f"C++ invalid amino acid test failed: {e}"
        finally:
            os.unlink(temp_path)


class TestGenomicsIntegration:
    """Test integration with main package."""
    
    def test_genomics_from_main_package(self):
        """Test importing genomics from main noLZSS package."""
        try:
            import noLZSS
            
            # The genomics subpackage should be accessible
            # Note: It won't be in __all__ until we implement it
            import noLZSS.genomics
            
            print("Genomics accessible from main package")
        except Exception as e:
            print(f"Warning: Genomics not accessible from main package: {e}")
            assert False, f"Genomics not accessible from main package: {e}"
    
    def test_genomics_utils_integration(self):
        """Test that genomics can use utils functions."""
        try:
            from noLZSS.genomics import is_dna_sequence, is_protein_sequence
            
            # Test that genomics-related functions work
            assert is_dna_sequence("ATCG")
            assert is_protein_sequence("ACDEFG")
            
            print("Genomics utils integration works")
        except Exception as e:
            print(f"Warning: Genomics utils integration failed: {e}")
            assert False, f"Genomics utils integration failed: {e}"


if __name__ == "__main__":
    # Run tests without pytest
    test_classes = [TestGenomicsStructure, TestFutureGenomicsFeatures, TestFASTAFunctions, TestCppFastaFunctions, TestProcessFastaWithPlots, TestGenomicsIntegration]
    
    total_tests = 0
    passed_tests = 0
    
    for test_class in test_classes:
        print(f"\n=== {test_class.__name__} ===")
        instance = test_class()
        methods = [method for method in dir(instance) if method.startswith('test_')]
        
        for method_name in methods:
            total_tests += 1
            try:
                method = getattr(instance, method_name)
                success = method()
                if success is not False:  # None or True counts as success
                    passed_tests += 1
            except Exception as e:
                print(f"Warning: {method_name}: {e}")
    
    print(f"\n=== Summary ===")
    print(f"Tests passed: {passed_tests}/{total_tests}")
    if passed_tests == total_tests:
        print("All genomics tests passed!")
    else:
        print(f"Failed tests: {total_tests - passed_tests}")
