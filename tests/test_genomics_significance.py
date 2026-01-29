"""
Tests for genomics significance analysis module.
"""

import os
import tempfile
from pathlib import Path

import pytest
import numpy as np

from noLZSS.genomics.significance import (
    clopper_pearson_upper,
    extract_factor_lengths,
    infer_length_significance,
    calculate_factor_length_threshold,
    plot_significance_analysis,
)
from noLZSS.utils import NoLZSSError


class TestClopperPearsonUpper:
    """Test Clopper-Pearson upper confidence bound calculation."""
    
    def test_clopper_pearson_upper_basic(self):
        """Test basic CP bound calculation."""
        # With k=5, n=100, alpha=0.05, upper bound should be around 0.11
        upper = clopper_pearson_upper(5, 100, 0.05)
        assert 0.08 < upper < 0.15, f"Expected upper bound ~0.11, got {upper}"
        
        # With k=50, n=100, alpha=0.05, upper bound should be around 0.60
        upper = clopper_pearson_upper(50, 100, 0.05)
        assert 0.55 < upper < 0.65, f"Expected upper bound ~0.60, got {upper}"
    
    def test_clopper_pearson_upper_edge_cases(self):
        """Test edge cases: k=0 and k=n."""
        # k=0: Should return 1 - alpha^(1/n)
        upper = clopper_pearson_upper(0, 100, 0.05)
        expected = 1.0 - (0.05 ** (1.0 / 100))
        assert abs(upper - expected) < 0.001, f"Expected {expected}, got {upper}"
        
        # k=n: Should return 1.0
        upper = clopper_pearson_upper(100, 100, 0.05)
        assert upper == 1.0, f"Expected 1.0 for k=n, got {upper}"
        
        # k=1, n=1: Should return 1.0
        upper = clopper_pearson_upper(1, 1, 0.05)
        assert upper == 1.0
    
    def test_clopper_pearson_upper_invalid_inputs(self):
        """Test error handling for invalid inputs."""
        # n must be positive
        with pytest.raises(ValueError, match="n must be positive"):
            clopper_pearson_upper(0, 0, 0.05)
        
        with pytest.raises(ValueError, match="n must be positive"):
            clopper_pearson_upper(0, -1, 0.05)
        
        # k must be between 0 and n
        with pytest.raises(ValueError, match="k must be between 0 and n"):
            clopper_pearson_upper(-1, 100, 0.05)
        
        with pytest.raises(ValueError, match="k must be between 0 and n"):
            clopper_pearson_upper(101, 100, 0.05)
        
        # alpha must be in (0, 1)
        with pytest.raises(ValueError, match="alpha must be in"):
            clopper_pearson_upper(5, 100, 0.0)
        
        with pytest.raises(ValueError, match="alpha must be in"):
            clopper_pearson_upper(5, 100, 1.0)
        
        with pytest.raises(ValueError, match="alpha must be in"):
            clopper_pearson_upper(5, 100, -0.1)


class TestExtractFactorLengths:
    """Test factor length extraction from lists and files."""
    
    def test_extract_factor_lengths_from_list(self):
        """Test extraction from list of tuples."""
        factors = [(0, 5, 0), (5, 3, 2), (8, 10, 1)]
        lengths = extract_factor_lengths(factors)
        
        assert isinstance(lengths, np.ndarray)
        assert lengths.dtype == np.int64
        np.testing.assert_array_equal(lengths, [5, 3, 10])
    
    def test_extract_factor_lengths_from_list_with_extra_fields(self):
        """Test extraction from tuples with extra fields (e.g., is_rc flag)."""
        factors = [(0, 5, 0, False), (5, 3, 2, True), (8, 10, 1, False)]
        lengths = extract_factor_lengths(factors)
        
        np.testing.assert_array_equal(lengths, [5, 3, 10])
    
    def test_extract_factor_lengths_empty_list(self):
        """Test extraction from empty list."""
        lengths = extract_factor_lengths([])
        
        assert isinstance(lengths, np.ndarray)
        assert len(lengths) == 0
        assert lengths.dtype == np.int64
    
    def test_extract_factor_lengths_from_file(self):
        """Test extraction from binary factor file."""
        # Create a temporary binary factor file
        with tempfile.TemporaryDirectory() as tmpdir:
            factor_file = Path(tmpdir) / "test_factors.bin"
            
            # Write binary factor file manually
            import struct
            factors = [(0, 5, 0), (5, 3, 2), (8, 10, 1)]
            
            with open(factor_file, 'wb') as f:
                # Write factors
                for pos, length, ref in factors:
                    f.write(struct.pack('<QQQ', pos, length, ref))
                
                # Write footer (magic + metadata)
                num_factors = len(factors)
                num_sequences = 0
                num_sentinels = 0
                footer_size = 48  # Just the basic footer
                total_length = sum(length for _, length, _ in factors)
                
                f.write(b'noLZSSv2')
                f.write(struct.pack('<QQQQQ', num_factors, num_sequences, 
                                   num_sentinels, footer_size, total_length))
            
            # Extract lengths
            lengths = extract_factor_lengths(factor_file)
            np.testing.assert_array_equal(lengths, [5, 3, 10])
    
    def test_extract_factor_lengths_invalid_input(self):
        """Test error handling for invalid input types."""
        # Invalid type
        with pytest.raises(ValueError, match="must be a list of tuples or a file path"):
            extract_factor_lengths(12345)
        
        with pytest.raises(ValueError, match="must be a list of tuples or a file path"):
            extract_factor_lengths({"factors": []})
    
    def test_extract_factor_lengths_invalid_tuples(self):
        """Test error handling for invalid tuple structure."""
        # Tuple with only one element
        with pytest.raises(ValueError, match="must be a tuple with at least 2 elements"):
            extract_factor_lengths([(0,)])
        
        # Not a tuple at all
        with pytest.raises(ValueError, match="must be a tuple with at least 2 elements"):
            extract_factor_lengths([[0, 5, 0]])
    
    def test_extract_factor_lengths_file_not_found(self):
        """Test error handling for non-existent file."""
        with pytest.raises(NoLZSSError, match="File not found"):
            extract_factor_lengths("/nonexistent/path/to/file.bin")


class TestInferLengthSignificance:
    """Test core significance inference function."""
    
    def test_infer_length_significance_basic(self):
        """Test basic significance inference."""
        # Real genome has longer factors
        real_lengths = np.array([10, 15, 20, 25, 30])
        # Shuffled genome has shorter factors (more samples for better statistics)
        shuf_lengths = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10] * 10)  # 90 samples
        
        result = infer_length_significance(real_lengths, shuf_lengths, tau_expected_fp=1.0)
        
        # Check structure
        assert 'N_real' in result
        assert 'N_shuf' in result
        assert 'L_star' in result
        assert 'rarity_scores_real' in result
        assert 'p_any_ge' in result
        
        # Check values
        assert result['N_real'] == 5
        assert result['N_shuf'] == 90
        # L_star may or may not be found depending on the data
        assert result['L_star'] is None or result['L_star'] > 0
        assert len(result['rarity_scores_real']) == 5
        
        # Longer factors should have lower rarity scores (more rare in shuffled)
        assert result['rarity_scores_real'][-1] < result['rarity_scores_real'][0]
    
    def test_infer_length_significance_no_threshold(self):
        """Test case where no L* meets the criterion."""
        # All real factors are short, within shuffled range
        real_lengths = np.array([5, 6, 7, 8, 9])
        shuf_lengths = np.array([5, 6, 7, 8, 9, 10, 11, 12])
        
        # Very stringent tau
        result = infer_length_significance(real_lengths, shuf_lengths, tau_expected_fp=0.01)
        
        # May not find a threshold
        # Just check that function completes without error
        assert result['L_star'] is None or result['L_star'] > 0
    
    def test_infer_length_significance_all_significant(self):
        """Test case where all real factors are significant."""
        # Real factors all very long
        real_lengths = np.array([50, 60, 70, 80, 90])
        # Shuffled factors all short (many samples for better statistics)
        shuf_lengths = np.array([1, 2, 3, 4, 5] * 20)  # 100 samples
        
        result = infer_length_significance(real_lengths, shuf_lengths, tau_expected_fp=2.0)
        
        # All rarity scores should be 0 (never seen in shuffled)
        np.testing.assert_array_equal(result['rarity_scores_real'], [0.0] * 5)
        
        # L_star should be found and small (since real factors are all far from shuffled)
        assert result['L_star'] is not None
        assert result['L_star'] <= 50
    
    def test_infer_length_significance_p_any_ge(self):
        """Test genome-wide exceedance probability function."""
        real_lengths = np.array([10, 15, 20, 25, 30])
        shuf_lengths = np.array([5, 10, 15, 20, 25])
        
        result = infer_length_significance(real_lengths, shuf_lengths)
        
        p_any_ge = result['p_any_ge']
        
        # Check that it's callable
        assert callable(p_any_ge)
        
        # Check that probabilities are valid
        p1 = p_any_ge(5)
        p2 = p_any_ge(25)
        p3 = p_any_ge(100)
        
        assert 0 <= p1 <= 1
        assert 0 <= p2 <= 1
        assert 0 <= p3 <= 1
        
        # Higher L should have lower probability
        assert p1 >= p2 >= p3
    
    def test_infer_length_significance_empty_shuffled(self):
        """Test error handling for empty shuffled genome."""
        real_lengths = np.array([10, 15, 20])
        shuf_lengths = np.array([])
        
        with pytest.raises(ValueError, match="Shuffled genome must have at least one factor"):
            infer_length_significance(real_lengths, shuf_lengths)


class TestCalculateFactorLengthThreshold:
    """Test main user-facing API function."""
    
    def test_calculate_factor_length_threshold_integration(self):
        """Test end-to-end integration with binary files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            real_file = Path(tmpdir) / "real_factors.bin"
            shuf_file = Path(tmpdir) / "shuf_factors.bin"
            
            # Create binary files
            import struct
            
            # Real genome factors (longer)
            real_factors = [(i*10, 20+i*5, i) for i in range(10)]
            with open(real_file, 'wb') as f:
                for pos, length, ref in real_factors:
                    f.write(struct.pack('<QQQ', pos, length, ref))
                f.write(b'noLZSSv2')
                f.write(struct.pack('<QQQQQ', len(real_factors), 0, 0, 48, 
                                   sum(l for _, l, _ in real_factors)))
            
            # Shuffled genome factors (shorter)
            shuf_factors = [(i*5, 5+i, i) for i in range(20)]
            with open(shuf_file, 'wb') as f:
                for pos, length, ref in shuf_factors:
                    f.write(struct.pack('<QQQ', pos, length, ref))
                f.write(b'noLZSSv2')
                f.write(struct.pack('<QQQQQ', len(shuf_factors), 0, 0, 48,
                                   sum(l for _, l, _ in shuf_factors)))
            
            # Calculate threshold
            result = calculate_factor_length_threshold(real_file, shuf_file)
            
            assert result['N_real'] == 10
            assert result['N_shuf'] == 20
            assert 'L_star' in result
            assert len(result['rarity_scores_real']) == 10
    
    def test_calculate_factor_length_threshold_with_plot(self):
        """Test with plot output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            real_file = Path(tmpdir) / "real_factors.bin"
            shuf_file = Path(tmpdir) / "shuf_factors.bin"
            plot_file = Path(tmpdir) / "plot.png"
            
            # Create minimal binary files
            import struct
            
            real_factors = [(0, 20, 0), (20, 25, 0)]
            with open(real_file, 'wb') as f:
                for pos, length, ref in real_factors:
                    f.write(struct.pack('<QQQ', pos, length, ref))
                f.write(b'noLZSSv2')
                f.write(struct.pack('<QQQQQ', len(real_factors), 0, 0, 48, 45))
            
            shuf_factors = [(i, 5+i, 0) for i in range(10)]
            with open(shuf_file, 'wb') as f:
                for pos, length, ref in shuf_factors:
                    f.write(struct.pack('<QQQ', pos, length, ref))
                f.write(b'noLZSSv2')
                f.write(struct.pack('<QQQQQ', len(shuf_factors), 0, 0, 48,
                                   sum(l for _, l, _ in shuf_factors)))
            
            # Calculate with plot
            result = calculate_factor_length_threshold(
                real_file, shuf_file, plot_output=plot_file
            )
            
            assert result['N_real'] == 2
            # Plot should be created (if matplotlib available)
            # We don't assert its existence since matplotlib might not be available
    
    def test_calculate_factor_length_threshold_file_not_found(self):
        """Test error handling for missing files."""
        with pytest.raises(FileNotFoundError, match="Real factors file not found"):
            calculate_factor_length_threshold(
                "/nonexistent/real.bin", 
                "/nonexistent/shuf.bin"
            )


class TestPlotSignificanceAnalysis:
    """Test visualization function."""
    
    def test_plot_significance_analysis(self):
        """Test plot generation (if matplotlib available)."""
        # Create a simple result dictionary
        real_lengths = np.array([10, 15, 20, 25, 30])
        shuf_lengths = np.array([5, 6, 7, 8, 9, 10, 11, 12])
        
        result = infer_length_significance(real_lengths, shuf_lengths)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            plot_path = Path(tmpdir) / "test_plot.png"
            
            # Should not raise an error (will warn if matplotlib not available)
            plot_significance_analysis(result, save_path=plot_path, show_plot=False)
            
            # If matplotlib is available, file should be created
            # We don't assert this because matplotlib might not be available
    
    def test_plot_significance_analysis_missing_keys(self):
        """Test error handling for invalid result dictionary."""
        incomplete_result = {
            'uniq_L': np.array([1, 2, 3]),
            'S0': np.array([0.5, 0.3, 0.1])
            # Missing other required keys
        }
        
        with pytest.raises(ValueError, match="missing required keys"):
            plot_significance_analysis(incomplete_result, show_plot=False)
    
    def test_plot_significance_analysis_no_save_no_show(self):
        """Test with no save and no show."""
        real_lengths = np.array([10, 15, 20])
        shuf_lengths = np.array([5, 6, 7, 8, 9])
        
        result = infer_length_significance(real_lengths, shuf_lengths)
        
        # Should complete without error
        plot_significance_analysis(result, save_path=None, show_plot=False)


if __name__ == "__main__":
    # Run tests when executed directly
    pytest.main([__file__, "-v"])
