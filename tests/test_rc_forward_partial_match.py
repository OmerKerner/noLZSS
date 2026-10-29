"""
Regression tests for forward "partial-edge" matches in the reverse-complement-aware factorization.

The RC-aware tree walk used to accept a forward candidate only from suffix-tree nodes whose
earliest occurrence fits entirely before the current position. When the next deeper node on
the path had an earliest occurrence that overlapped the current position, the shorter prefix
of that occurrence (which the base, non-RC algorithm handles on the edge between the two
nodes) was dropped, and a shorter factor was emitted.

These tests pin the fixed behaviour against a brute-force oracle that implements the
documented semantics directly: at each position take the longest non-overlapping previous
occurrence, forward preferred at equal length, reverse-complement matches must have
length >= 2, earliest start (forward) / earliest end (RC) among candidates of equal length.
"""

import random
import tempfile
from pathlib import Path

import pytest

try:
    import noLZSS._noLZSS as cpp
    from noLZSS.utils import read_factors_binary_file_with_metadata
    HAS_CPP_EXTENSION = True
except ImportError:  # pragma: no cover - exercised only when the extension is missing
    HAS_CPP_EXTENSION = False

pytestmark = pytest.mark.skipif(not HAS_CPP_EXTENSION, reason="C++ extension not available")

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def oracle_factorize_dna_w_rc(s: str):
    """Brute-force reference implementation; returns [(start, length, ref, is_rc), ...]."""
    n = len(s)
    i = 0
    out = []
    while i < n:
        candidates = []  # (length, kind, start) with kind 0 = forward, 1 = reverse complement
        for j in range(i):
            length = 0
            while i + length < n and j + length < i and s[j + length] == s[i + length]:
                length += 1
            if length > 0:
                candidates.append((length, 0, j))
        for end in range(i):  # RC occurrence T[start..end] with end < i
            length = 0
            while i + length < n and end - length >= 0 and COMPLEMENT[s[end - length]] == s[i + length]:
                length += 1
            if length > 1:  # an RC match must be strictly longer than a literal
                candidates.append((length, 1, end - length + 1))
        if not candidates:
            out.append((i, 1, i, False))
            i += 1
            continue
        best_len = max(c[0] for c in candidates)
        best = [c for c in candidates if c[0] == best_len]
        forward = [c for c in best if c[1] == 0]
        pick = min(forward, key=lambda c: c[2]) if forward else min(best, key=lambda c: c[2])
        length, kind, start = pick
        out.append((i, length, start, bool(kind)))
        i += length
    return out


def test_partial_edge_forward_match_regression():
    """'CT' at position 5 is a valid length-2 forward match for position 7; it used to be missed."""
    factors = cpp.factorize_dna_w_rc(b"GCGCACTCTCG")
    assert factors == oracle_factorize_dna_w_rc("GCGCACTCTCG")
    assert (7, 2, 5, False) in factors


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_matches_oracle_on_random_dna(seed):
    rng = random.Random(seed)
    for _ in range(150):
        n = rng.randint(1, 40)
        s = "".join(rng.choice("ACGT") for _ in range(n))
        assert cpp.factorize_dna_w_rc(s.encode()) == oracle_factorize_dna_w_rc(s), s


def test_rc_aware_never_shorter_than_plain_at_same_position():
    """Where both factorizations start a factor at the same position, the RC-aware one is at least as long."""
    rng = random.Random(11)
    for _ in range(200):
        n = rng.randint(1, 40)
        s = "".join(rng.choice("ACGT") for _ in range(n))
        plain = {start: length for start, length, _ in cpp.factorize(s.encode())}
        for start, length, _, _ in cpp.factorize_dna_w_rc(s.encode()):
            if start in plain:
                assert length >= plain[start], s


def test_parallel_rc_matches_sequential():
    rng = random.Random(5)
    for num_threads in (1, 2, 3):
        for _ in range(20):
            n = rng.randint(8, 60)
            s = "".join(rng.choice("ACGT") for _ in range(n))
            expected = cpp.factorize_dna_w_rc(s.encode())
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / "factors.bin"
                cpp.parallel_factorize_dna_w_rc_to_file(s.encode(), str(out), num_threads)
                got = read_factors_binary_file_with_metadata(out)["factors"]
            assert [tuple(f) for f in got] == expected, (s, num_threads)
