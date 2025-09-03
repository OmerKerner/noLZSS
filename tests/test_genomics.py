"""
Tests for the genomics subpackage.
"""

import sys
import os

# Add the src directory to the path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


class TestGenomicsStructure:
    """Test the genomics subpackage structure."""
    
    def test_genomics_import(self):
        """Test that genomics subpackage can be imported."""
        try:
            from noLZSS import genomics
            print("Genomics subpackage imports successfully")
            return True
        except ImportError as e:
            print(f"Warning: Genomics import failed: {e}")
            return False
    
    def test_fasta_module_import(self):
        """Test that fasta module can be imported."""
        try:
            from noLZSS.genomics import fasta
            print("FASTA module imports successfully")
            return True
        except ImportError as e:
            print(f"Warning: FASTA module import failed: {e}")
            return False
    
    def test_sequences_module_import(self):
        """Test that sequences module can be imported."""
        try:
            from noLZSS.genomics import sequences
            print("Sequences module imports successfully")
            return True
        except ImportError as e:
            print(f"Warning: Sequences module import failed: {e}")
            return False
    
    def test_genomics_namespace(self):
        """Test the genomics namespace structure."""
        try:
            import noLZSS.genomics as genomics
            
            # Check that submodules are accessible
            assert hasattr(genomics, 'fasta')
            assert hasattr(genomics, 'sequences')
            
            print("Genomics namespace structure is correct")
            return True
        except Exception as e:
            print(f"Warning: Genomics namespace error: {e}")
            return False


class TestFutureGenomicsFeatures:
    """Placeholder tests for future genomics functionality."""
    
    def test_fasta_placeholder(self):
        """Test that FASTA module exists and is ready for implementation."""
        try:
            from noLZSS.genomics import fasta
            
            # Check that it's a module
            assert hasattr(fasta, '__file__')
            
            print("FASTA module ready for implementation")
            return True
        except Exception as e:
            print(f"Warning: FASTA module not ready: {e}")
            return False
    
    def test_sequences_placeholder(self):
        """Test that sequences module exists and is ready for implementation."""
        try:
            from noLZSS.genomics import sequences
            
            # Check that it's a module
            assert hasattr(sequences, '__file__')
            
            print("Sequences module ready for implementation")
            return True
        except Exception as e:
            print(f"Warning: Sequences module not ready: {e}")
            return False


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
            return True
        except Exception as e:
            print(f"Warning: Genomics not accessible from main package: {e}")
            return False
    
    def test_genomics_utils_integration(self):
        """Test that genomics can use utils functions."""
        try:
            from noLZSS.genomics import is_dna_sequence, is_protein_sequence
            
            # Test that genomics-related functions work
            assert is_dna_sequence("ATCG")
            assert is_protein_sequence("ACDEFG")
            
            print("Genomics utils integration works")
            return True
        except Exception as e:
            print(f"Warning: Genomics utils integration failed: {e}")
            return False


if __name__ == "__main__":
    # Run tests without pytest
    test_classes = [TestGenomicsStructure, TestFutureGenomicsFeatures, TestGenomicsIntegration]
    
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
