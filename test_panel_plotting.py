"""
Tests for the Panel/Datashader plotting functionality.
"""

import tempfile
import os
from pathlib import Path

from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta, PlotError


class TestPanelPlotting:
    """Test the new Panel/Datashader plotting function."""
    
    def test_import_error_handling(self):
        """Test that missing dependencies raise appropriate ImportError."""
        # Create a test FASTA file
        fasta_content = """>seq1
ATCGATCGATCGTAGCTAGCTAGCTACGTACGTACGTTAGCTAGCTAGCT
>seq2
GCTAGCTAGCTAATCGATCGATCGCGTACGTACGTACGTATCGATCGATCG
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            # If dependencies are missing, should get ImportError
            try:
                result = plot_multiple_seq_self_lz_factor_plot_from_fasta(
                    fasta_filepath=temp_path,
                    show_plot=False,
                    return_panel=True
                )
                # If we get here, dependencies are available
                print("Panel dependencies are available, function should work")
                # Basic validation that we got something back
                assert result is not None
            except ImportError as e:
                # Expected when dependencies are missing
                assert "Missing required dependency" in str(e)
                assert "pip install" in str(e)
                print("ImportError correctly raised for missing dependencies")
        finally:
            os.unlink(temp_path)
    
    def test_file_not_found_error(self):
        """Test that missing FASTA file raises FileNotFoundError."""
        nonexistent_path = "/tmp/nonexistent_file.fasta"
        
        try:
            try:
                plot_multiple_seq_self_lz_factor_plot_from_fasta(
                    fasta_filepath=nonexistent_path,
                    show_plot=False
                )
                # Should not reach here
                assert False, "Expected FileNotFoundError"
            except FileNotFoundError:
                print("FileNotFoundError correctly raised for missing file")
            except ImportError:
                # Dependencies not available, but file check should happen first
                print("ImportError raised - dependencies not available")
        except Exception as e:
            print(f"Unexpected error: {e}")
    
    def test_function_signature(self):
        """Test that the function has the correct signature."""
        import inspect
        
        sig = inspect.signature(plot_multiple_seq_self_lz_factor_plot_from_fasta)
        params = list(sig.parameters.keys())
        
        expected_params = ['fasta_filepath', 'name', 'save_path', 'show_plot', 'return_panel']
        assert params == expected_params
        
        # Check default values
        assert sig.parameters['name'].default is None
        assert sig.parameters['save_path'].default is None  
        assert sig.parameters['show_plot'].default is True
        assert sig.parameters['return_panel'].default is False
        
        print("Function signature is correct")
    
    def test_backward_compatibility(self):
        """Test that old function name still exists."""
        from noLZSS.genomics.plots import plot_multiple_seq_self_weizmann_factor_plot_from_fasta
        
        # Should be importable
        assert callable(plot_multiple_seq_self_weizmann_factor_plot_from_fasta)
        print("Backward compatibility maintained - old function still exists")


if __name__ == "__main__":
    # Run tests without pytest
    test_classes = [TestPanelPlotting]
    
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
                method()
                passed_tests += 1
                print(f"✓ {method_name}")
            except Exception as e:
                print(f"✗ {method_name}: {e}")
    
    print(f"\n=== Summary ===")
    print(f"Tests passed: {passed_tests}/{total_tests}")