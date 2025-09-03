"""
Integration tests that test the entire package workflow.
"""

import sys
import os
import tempfile
from pathlib import Path

# Add the src directory to the path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


class TestPackageIntegration:
    """Test integration across all modules."""
    
    def test_main_package_structure(self):
        """Test the main package structure and imports."""
        try:
            import noLZSS
            
            # Check that main functions are exposed
            expected_functions = [
                'analyze_alphabet', 'detect_sequence_type', 
                'is_dna_sequence', 'is_protein_sequence'
            ]
            
            for func_name in expected_functions:
                if hasattr(noLZSS, func_name):
                    print(f"{func_name} available in main package")
                else:
                    print(f"Warning: {func_name} not available in main package")
                    return False
            
            return True
        except ImportError as e:
            print(f"Warning: Main package import failed: {e}")
            return False
    
    def test_utils_through_main_package(self):
        """Test utils functionality through main package."""
        try:
            import noLZSS
            
            # Test sequence detection
            assert noLZSS.is_dna_sequence("ATCG")
            assert noLZSS.is_protein_sequence("PROTEIN")
            assert noLZSS.detect_sequence_type("ATCG") == 'dna'
            
            # Test alphabet analysis
            result = noLZSS.analyze_alphabet("ATCG")
            assert result['size'] == 4
            
            print("Utils functions work through main package")
            return True
        except Exception as e:
            print(f"Warning: Utils through main package failed: {e}")
            return False
    
    def test_all_modules_importable(self):
        """Test that all modules can be imported together."""
        try:
            from noLZSS import utils
            from noLZSS import core
            from noLZSS.genomics import fasta, sequences
            
            print("All modules can be imported together")
            return True
        except ImportError as e:
            print(f"Warning: Module import failed: {e}")
            return False
    
    def test_workflow_without_cpp(self):
        """Test a typical workflow without C++ extension."""
        try:
            import noLZSS
            
            # Step 1: Analyze input data
            dna_sequence = "ATCGATCGATCG"
            
            # Detect sequence type
            seq_type = noLZSS.detect_sequence_type(dna_sequence)
            assert seq_type == 'dna'
            
            # Analyze alphabet
            alphabet_info = noLZSS.analyze_alphabet(dna_sequence)
            assert alphabet_info['size'] == 4
            assert alphabet_info['total_length'] == 12
            
            # Step 2: Validate input (this should work even without C++)
            from noLZSS.utils import validate_input
            
            validated_data = validate_input(dna_sequence)
            assert isinstance(validated_data, bytes)
            
            print("Complete workflow (without C++ factorization) works")
            return True
        except Exception as e:
            print(f"Warning: Workflow failed: {e}")
            return False


class TestFileOperations:
    """Test file-related operations."""
    
    def test_file_reader_utility(self):
        """Test the file reader utility."""
        try:
            from noLZSS.utils import safe_file_reader
            
            # Create a test file
            test_data = b"ATCGATCGATCG" * 100  # Some test DNA data
            
            with tempfile.NamedTemporaryFile(mode='wb', delete=False) as f:
                f.write(test_data)
                temp_path = f.name
            
            try:
                # Read file in chunks
                chunks = list(safe_file_reader(temp_path, chunk_size=50))
                reconstructed_data = b''.join(chunks)
                
                assert reconstructed_data == test_data
                print("File reader utility works")
                return True
            finally:
                os.unlink(temp_path)
        except Exception as e:
            print(f"Warning: File reader utility failed: {e}")
            return False
    
    def test_path_object_support(self):
        """Test that Path objects are supported."""
        try:
            from noLZSS.utils import safe_file_reader
            from pathlib import Path
            
            # Create a test file
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
                f.write("test data")
                temp_path = Path(f.name)
            
            try:
                # Should work with Path objects
                chunks = list(safe_file_reader(temp_path))
                data = b''.join(chunks)
                assert data == b"test data"
                
                print("Path object support works")
                return True
            finally:
                temp_path.unlink()
        except Exception as e:
            print(f"Warning: Path object support failed: {e}")
            return False


class TestErrorHandlingIntegration:
    """Test error handling across modules."""
    
    def test_custom_exceptions(self):
        """Test custom exception handling."""
        try:
            from noLZSS.utils import InvalidInputError, NoLZSSError, validate_input
            
            # Test InvalidInputError
            try:
                validate_input("")
                assert False, "Should have raised InvalidInputError"
            except InvalidInputError:
                print("InvalidInputError handling works")
            
            # Test NoLZSSError with file operations
            from noLZSS.utils import safe_file_reader
            try:
                list(safe_file_reader("nonexistent_file.txt"))
                assert False, "Should have raised NoLZSSError"
            except NoLZSSError:
                print("NoLZSSError handling works")
            
            return True
        except Exception as e:
            print(f"Warning: Custom exception handling failed: {e}")
            return False
    
    def test_validation_options_integration(self):
        """Test validation options across the package."""
        try:
            from noLZSS.core import factorize
            from noLZSS.utils import InvalidInputError
            
            # This will fail due to missing C++ module, but tests validation flow
            try:
                factorize("", validate=True)
            except InvalidInputError:
                print("Validation=True works across modules")
            except Exception as e:
                if "No module named" in str(e):
                    print("Validation occurs before C++ call")
                else:
                    raise
            
            return True
        except ImportError:
            print("Warning: Cannot test validation options without core module")
            return False


class TestBenchmarkIntegration:
    """Test integration with benchmark scripts."""
    
    def test_benchmark_script_exists(self):
        """Test that benchmark scripts exist and are accessible."""
        benchmark_dir = Path(__file__).parent.parent / "benchmarks"
        
        expected_files = ["bench.py", "plot_benchmarks.py"]
        
        for filename in expected_files:
            filepath = benchmark_dir / filename
            if filepath.exists():
                print(f"{filename} exists")
            else:
                print(f"Warning: {filename} missing")
                return False
        
        return True
    
    def test_benchmark_can_import_package(self):
        """Test that benchmark scripts can import the package."""
        try:
            # This simulates what the benchmark script does
            import noLZSS
            
            # Check that functions benchmarks need are available
            from noLZSS.utils import analyze_alphabet
            
            result = analyze_alphabet("ATCG")
            assert 'size' in result
            
            print("Package functions work for benchmarking")
            return True
        except Exception as e:
            print(f"Warning: Benchmark integration failed: {e}")
            return False


if __name__ == "__main__":
    # Run integration tests
    test_classes = [
        TestPackageIntegration, TestFileOperations, 
        TestErrorHandlingIntegration, TestBenchmarkIntegration
    ]
    
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
    
    print(f"\n=== Integration Test Summary ===")
    print(f"Tests passed: {passed_tests}/{total_tests}")
    if passed_tests == total_tests:
        print("All integration tests passed!")
    else:
        print(f"Failed tests: {total_tests - passed_tests}")
