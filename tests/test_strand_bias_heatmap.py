import math
import pytest

from noLZSS.genomics.plots import _compute_strand_bias_grid


def test_strand_bias_grid_splits_crossing_factors():
    # Two forward segments split across bins, two reverse-complement segments
    factors = [
        (0, 4, 0, False),  # forward diagonal from (0,0) -> (4,4)
        (0, 4, 0, True),   # reverse diagonal from (0,4) -> (4,0)
        (0, 2, 0, True),   # short reverse segment entirely in lower-left bin
    ]

    x_edges, y_edges, forward_grid, rc_grid, bias_grid = _compute_strand_bias_grid(factors, 2, total_length=4)

    assert forward_grid.shape == (2, 2)
    assert rc_grid.shape == (2, 2)

    # Forward coverage splits into two bins
    assert forward_grid[0, 0] == pytest.approx(2.0)
    assert forward_grid[1, 1] == pytest.approx(2.0)

    # Reverse coverage crosses boundaries and lands in three bins
    assert rc_grid[0, 0] == pytest.approx(2.0)
    assert rc_grid[0, 1] == pytest.approx(2.0)
    assert rc_grid[1, 0] == pytest.approx(2.0)

    # Totals reflect nucleotide coverage, not factor counts
    assert forward_grid.sum() == pytest.approx(4.0)
    assert rc_grid.sum() == pytest.approx(6.0)

    # Bias grid is finite everywhere because bins receive coverage
    assert bias_grid.mask.sum() == 0

    # Mixed bin has modest forward enrichment
    expected_bias = math.log2((0.5) / (1/3))  # (0.5 forward share) / (0.333 reverse share)
    assert bias_grid[0, 0] == pytest.approx(expected_bias, rel=1e-3)

    # Pure forward / pure reverse bins show strong signed bias
    assert bias_grid[1, 1] > 10
    assert bias_grid[0, 1] < -10
