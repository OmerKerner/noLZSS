"""
Tests for the ``ref_policy`` option ("first" vs "last" reference selection).

Both policies emit the same factor boundaries: at each position the factor length is the maximum
over all earlier non-overlapping occurrences, so changing which occurrence is reported cannot change
the length or the next factor start. The policies differ only in the reported reference:

* "first": earliest occurrence (leftmost start for forward matches, leftmost end for reverse
  complement matches); at equal length a forward match beats a reverse-complement match.
* "last": closest occurrence (largest non-overlapping start/end); at equal length the strand
  whose occurrence ends closest to the factor wins, forward on an exact tie.

The brute-force oracles below implement these rules directly and are the reference semantics.
"""

import random
import tempfile
from pathlib import Path

import pytest

try:
    import noLZSS
    import noLZSS._noLZSS as cpp
    from noLZSS.utils import read_factors_binary_file_with_metadata, validate_ref_policy
    HAS_CPP_EXTENSION = True
except ImportError:  # pragma: no cover
    HAS_CPP_EXTENSION = False

pytestmark = pytest.mark.skipif(not HAS_CPP_EXTENSION, reason="C++ extension not available")

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
POLICIES = ("first", "last")


# ----------------------------------------------------------------------------------------------
# Oracles
# ----------------------------------------------------------------------------------------------

def oracle_plain(s: str, policy: str, start_pos: int = 0):
    """Brute-force plain factorization; returns [(start, length, ref), ...]."""
    n = len(s)
    i = start_pos
    out = []
    while i < n:
        cands = []
        for j in range(i):
            length = 0
            while i + length < n and j + length < i and s[j + length] == s[i + length]:
                length += 1
            if length > 0:
                cands.append((length, j))
        if not cands:
            out.append((i, 1, i))
            i += 1
            continue
        best_len = max(c[0] for c in cands)
        starts = [j for length, j in cands if length == best_len]
        ref = min(starts) if policy == "first" else max(starts)
        out.append((i, best_len, ref))
        i += best_len
    return out


def oracle_rc(s: str, policy: str, start_pos: int = 0):
    """
    Brute-force RC-aware factorization; returns [(start, length, ref, is_rc), ...].

    ``s`` may contain non-ACGT sentinel characters: they never take part in a reverse-complement
    match, which mirrors the prepared string used by the C++ code (sentinels in the RC half differ
    from those in the forward half).
    """
    n = len(s)
    i = start_pos
    out = []
    while i < n:
        cands = []  # (length, kind, start, end) kind 0 = forward, 1 = RC
        for j in range(i):
            length = 0
            while i + length < n and j + length < i and s[j + length] == s[i + length]:
                length += 1
            if length > 0:
                cands.append((length, 0, j, j + length - 1))
        for end in range(i):
            length = 0
            while (i + length < n and end - length >= 0
                   and COMPLEMENT.get(s[end - length]) is not None
                   and COMPLEMENT[s[end - length]] == s[i + length]):
                length += 1
            if length > 1:
                cands.append((length, 1, end - length + 1, end))
        if not cands:
            out.append((i, 1, i, False))
            i += 1
            continue
        best_len = max(c[0] for c in cands)
        best = [c for c in cands if c[0] == best_len]
        if policy == "first":
            fwd = [c for c in best if c[1] == 0]
            pick = min(fwd, key=lambda c: c[2]) if fwd else min(best, key=lambda c: c[3])
        else:
            # closest end wins; forward (kind 0) wins an exact tie
            pick = max(best, key=lambda c: (c[3], -c[1]))
        length, kind, start, _ = pick
        out.append((i, length, start, bool(kind)))
        i += length
    return out


def random_dna(rng, lo=1, hi=40):
    return "".join(rng.choice("ACGT") for _ in range(rng.randint(lo, hi)))


def random_text(rng, alphabet, lo=1, hi=40):
    return "".join(rng.choice(alphabet) for _ in range(rng.randint(lo, hi)))


def read_bin(path):
    return [tuple(f) for f in read_factors_binary_file_with_metadata(path)["factors"]]


def as4(factors):
    """Plain (start, length, ref) factors as the 4-tuples the no-RC FASTA/binary readers produce."""
    return [(a, b, c, False) for a, b, c in factors]


def prepared_text(prep):
    """(prepared_string, forward_part_without_its_final_sentinel) from a prepare_* result."""
    prepared = prep["prepared_string"] if isinstance(prep, dict) else prep[0]
    n = prep["original_length"] if isinstance(prep, dict) else prep[1]
    return prepared, prepared[: n - 1]


# ----------------------------------------------------------------------------------------------
# Fixed vectors
# ----------------------------------------------------------------------------------------------

def test_fixed_vectors_plain():
    first = cpp.factorize(b"abcabcabc")
    last = cpp.factorize(b"abcabcabc", ref_policy="last")
    assert first == cpp.factorize(b"abcabcabc", ref_policy="first")
    assert (6, 3, 0) in first
    assert (6, 3, 3) in last
    assert [(s, l) for s, l, _ in first] == [(s, l) for s, l, _ in last]

    first = cpp.factorize(b"abracadabra")
    last = cpp.factorize(b"abracadabra", ref_policy="last")
    assert (5, 1, 0) in first
    assert (5, 1, 3) in last


def test_fixed_vectors_rc():
    first = cpp.factorize_dna_w_rc(b"ACGTACGTACGT")
    last = cpp.factorize_dna_w_rc(b"ACGTACGTACGT", ref_policy="last")
    assert (8, 4, 0, False) in first
    assert (8, 4, 4, False) in last

    first = cpp.factorize_dna_w_rc(b"ACGTTTACGT")
    last = cpp.factorize_dna_w_rc(b"ACGTTTACGT", ref_policy="last")
    assert (5, 1, 3, False) in first
    assert (5, 1, 4, False) in last

    # Equal-length forward and RC candidates: "last" picks the closer RC occurrence
    first = cpp.factorize_dna_w_rc(b"TTTTCGAACTCGT")
    last = cpp.factorize_dna_w_rc(b"TTTTCGAACTCGT", ref_policy="last")
    assert (9, 3, 3, False) in first
    assert (9, 3, 4, True) in last

    # Exact palindromic tie (occurrence equals its own reverse complement) stays forward
    assert (2, 2, 0, False) in cpp.factorize_dna_w_rc(b"ATAT", ref_policy="last")


def test_default_policy_is_first():
    rng = random.Random(42)
    for _ in range(50):
        s = random_dna(rng).encode()
        assert cpp.factorize(s) == cpp.factorize(s, ref_policy="first")
        assert cpp.factorize_dna_w_rc(s) == cpp.factorize_dna_w_rc(s, ref_policy="first")


# ----------------------------------------------------------------------------------------------
# Oracle comparisons
# ----------------------------------------------------------------------------------------------

@pytest.mark.parametrize("policy", POLICIES)
@pytest.mark.parametrize("alphabet", ["ab", "abc", "ACGT"])
def test_plain_matches_oracle(policy, alphabet):
    rng = random.Random(hash((policy, alphabet)) & 0xFFFF)
    for _ in range(150):
        s = random_text(rng, alphabet)
        assert cpp.factorize(s.encode(), ref_policy=policy) == oracle_plain(s, policy), s


@pytest.mark.parametrize("policy", POLICIES)
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_rc_matches_oracle(policy, seed):
    rng = random.Random(seed)
    for _ in range(150):
        s = random_dna(rng)
        assert cpp.factorize_dna_w_rc(s.encode(), ref_policy=policy) == oracle_rc(s, policy), s


@pytest.mark.parametrize("policy", POLICIES)
def test_rc_low_entropy_inputs_match_oracle(policy):
    """Repetitive strings have many equal-length candidates on both strands."""
    rng = random.Random(7)
    for _ in range(100):
        unit = random_dna(rng, 1, 4)
        s = (unit * rng.randint(2, 8))[: rng.randint(4, 40)]
        assert cpp.factorize_dna_w_rc(s.encode(), ref_policy=policy) == oracle_rc(s, policy), s
        assert cpp.factorize(s.encode(), ref_policy=policy) == oracle_plain(s, policy), s


def test_boundaries_and_counts_identical_between_policies():
    rng = random.Random(9)
    for _ in range(100):
        s = random_dna(rng, 1, 60).encode()
        f, l = cpp.factorize(s), cpp.factorize(s, ref_policy="last")
        assert [(a, b) for a, b, _ in f] == [(a, b) for a, b, _ in l]
        assert cpp.count_factors(s) == len(l)
        f, l = cpp.factorize_dna_w_rc(s), cpp.factorize_dna_w_rc(s, ref_policy="last")
        assert [(a, b) for a, b, _, _ in f] == [(a, b) for a, b, _, _ in l]
        assert cpp.count_factors_dna_w_rc(s) == len(l)


def test_last_reference_is_never_further_than_first():
    rng = random.Random(13)
    for _ in range(100):
        s = random_dna(rng, 1, 60).encode()
        for f, l in zip(cpp.factorize(s), cpp.factorize(s, ref_policy="last")):
            assert l[2] >= f[2]
            if l[2] != l[0]:  # not a literal
                assert l[2] + l[1] <= l[0]  # still non-overlapping


# ----------------------------------------------------------------------------------------------
# Multiple sequences, reference sequences, start_pos
# ----------------------------------------------------------------------------------------------

@pytest.mark.parametrize("policy", POLICIES)
def test_multiple_dna_w_rc_matches_oracle(policy):
    rng = random.Random(21)
    for _ in range(60):
        seqs = [random_dna(rng, 1, 20) for _ in range(rng.randint(1, 3))]
        prepared, forward = prepared_text(cpp.prepare_multiple_dna_sequences_w_rc(seqs))
        got = cpp.factorize_multiple_dna_w_rc(prepared.encode("latin-1"), ref_policy=policy)
        assert got == oracle_rc(forward, policy), seqs


@pytest.mark.parametrize("policy", POLICIES)
def test_reference_functions_match_oracle(policy):
    rng = random.Random(33)
    for _ in range(60):
        ref, target = random_dna(rng, 1, 25), random_dna(rng, 1, 25)
        combined = ref + "\x01" + target
        got = cpp.factorize_w_reference(ref, target, ref_policy=policy)
        assert got == oracle_plain(combined, policy, start_pos=len(ref) + 1), (ref, target)

        got = cpp.factorize_dna_w_reference_seq(ref, target, ref_policy=policy)
        _, forward = prepared_text(cpp.prepare_multiple_dna_sequences_w_rc([ref, target]))
        assert got == oracle_rc(forward, policy, start_pos=len(ref) + 1), (ref, target)


# ----------------------------------------------------------------------------------------------
# Parallel implementations must agree with the sequential ones under both policies
# ----------------------------------------------------------------------------------------------

@pytest.mark.parametrize("policy", POLICIES)
def test_parallel_plain_matches_sequential(policy):
    rng = random.Random(5)
    for num_threads in (1, 2, 3):
        for _ in range(20):
            s = random_dna(rng, 8, 80)
            expected = as4(cpp.factorize(s.encode(), ref_policy=policy))
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / "f.bin"
                cpp.parallel_factorize_to_file(s, str(out), num_threads, 0, ref_policy=policy)
                assert read_bin(out) == expected, (s, num_threads)


@pytest.mark.parametrize("policy", POLICIES)
def test_parallel_rc_matches_sequential(policy):
    rng = random.Random(6)
    for num_threads in (1, 2, 3):
        for _ in range(20):
            s = random_dna(rng, 8, 80)
            expected = cpp.factorize_dna_w_rc(s.encode(), ref_policy=policy)
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / "f.bin"
                cpp.parallel_factorize_dna_w_rc_to_file(s, str(out), num_threads, ref_policy=policy)
                assert read_bin(out) == expected, (s, num_threads)


# ----------------------------------------------------------------------------------------------
# FASTA entry points
# ----------------------------------------------------------------------------------------------

@pytest.fixture
def fasta_file(tmp_path):
    rng = random.Random(77)
    seqs = {f"seq{i}": random_dna(rng, 10, 60) for i in range(3)}
    path = tmp_path / "in.fasta"
    path.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    return path, seqs


@pytest.mark.parametrize("policy", POLICIES)
def test_fasta_multiple_functions(fasta_file, tmp_path, policy):
    path, seqs = fasta_file
    prepared_rc, _ = prepared_text(cpp.prepare_multiple_dna_sequences_w_rc(list(seqs.values())))
    prepared_no, _ = prepared_text(cpp.prepare_multiple_dna_sequences_no_rc(list(seqs.values())))
    exp_rc = cpp.factorize_multiple_dna_w_rc(prepared_rc.encode("latin-1"), ref_policy=policy)
    exp_no = as4(cpp.factorize(prepared_no.encode("latin-1"), ref_policy=policy))

    assert cpp.factorize_fasta_multiple_dna_w_rc(str(path), ref_policy=policy)[0] == exp_rc
    assert cpp.factorize_fasta_multiple_dna_no_rc(str(path), ref_policy=policy)[0] == exp_no

    out = tmp_path / "rc.bin"
    cpp.write_factors_binary_file_fasta_multiple_dna_w_rc(str(path), str(out), ref_policy=policy)
    assert read_bin(out) == exp_rc
    out = tmp_path / "no.bin"
    cpp.write_factors_binary_file_fasta_multiple_dna_no_rc(str(path), str(out), ref_policy=policy)
    assert read_bin(out) == exp_no

    for threads in (1, 2):
        out = tmp_path / f"prc{threads}.bin"
        cpp.parallel_write_factors_binary_file_fasta_multiple_dna_w_rc(str(path), str(out), threads, ref_policy=policy)
        assert read_bin(out) == exp_rc
        out = tmp_path / f"pno{threads}.bin"
        cpp.parallel_write_factors_binary_file_fasta_multiple_dna_no_rc(str(path), str(out), threads, ref_policy=policy)
        assert read_bin(out) == exp_no


@pytest.mark.parametrize("policy", POLICIES)
def test_fasta_per_sequence_functions(fasta_file, tmp_path, policy):
    path, seqs = fasta_file
    exp_rc = [cpp.factorize_dna_w_rc(v.encode(), ref_policy=policy) for v in seqs.values()]
    exp_no = [as4(cpp.factorize(v.encode(), ref_policy=policy)) for v in seqs.values()]

    res = cpp.factorize_fasta_dna_w_rc_per_sequence(str(path), ref_policy=policy)
    assert [list(map(tuple, f)) for f in res[0]] == exp_rc
    res = cpp.factorize_fasta_dna_no_rc_per_sequence(str(path), ref_policy=policy)
    assert [list(map(tuple, f)) for f in res[0]] == exp_no

    for fn, exp, name in [
        (cpp.write_factors_binary_file_fasta_dna_w_rc_per_sequence, exp_rc, "w"),
        (cpp.write_factors_binary_file_fasta_dna_no_rc_per_sequence, exp_no, "n"),
    ]:
        out_dir = tmp_path / f"per_{name}"
        fn(str(path), str(out_dir), ref_policy=policy)
        for (seq_id, _), expected in zip(seqs.items(), exp):
            assert read_bin(out_dir / f"{seq_id}.bin") == expected
    for fn, exp, name in [
        (cpp.parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence, exp_rc, "pw"),
        (cpp.parallel_write_factors_binary_file_fasta_dna_no_rc_per_sequence, exp_no, "pn"),
    ]:
        out_dir = tmp_path / f"per_{name}"
        fn(str(path), str(out_dir), 2, ref_policy=policy)
        for (seq_id, _), expected in zip(seqs.items(), exp):
            assert read_bin(out_dir / f"{seq_id}.bin") == expected


@pytest.mark.parametrize("policy", POLICIES)
def test_fasta_reference_functions(tmp_path, policy):
    rng = random.Random(99)
    ref_seqs = [random_dna(rng, 10, 40) for _ in range(2)]
    tgt_seqs = [random_dna(rng, 10, 40) for _ in range(2)]
    ref_path, tgt_path = tmp_path / "ref.fa", tmp_path / "tgt.fa"
    ref_path.write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(ref_seqs)))
    tgt_path.write_text("".join(f">t{i}\n{s}\n" for i, s in enumerate(tgt_seqs)))

    _, forward = prepared_text(cpp.prepare_multiple_dna_sequences_w_rc(ref_seqs + tgt_seqs))
    target_start = sum(len(s) + 1 for s in ref_seqs)
    expected = oracle_rc(forward, policy, start_pos=target_start)

    got = cpp.factorize_dna_rc_w_ref_fasta_files(str(ref_path), str(tgt_path), ref_policy=policy)[0]
    assert got == expected
    assert all(start >= target_start for start, *_ in got)

    out = tmp_path / "ref.bin"
    cpp.write_factors_dna_w_reference_fasta_files_to_binary(str(ref_path), str(tgt_path), str(out), ref_policy=policy)
    assert read_bin(out) == got
    out = tmp_path / "pref.bin"
    cpp.parallel_write_factors_dna_w_reference_fasta_files_to_binary(str(ref_path), str(tgt_path), str(out), 2, ref_policy=policy)
    assert read_bin(out) == got


# ----------------------------------------------------------------------------------------------
# Python wrappers and validation
# ----------------------------------------------------------------------------------------------

def test_python_wrappers_pass_policy_through(tmp_path):
    s = "ACGTACGTACGT"
    assert noLZSS.factorize(s, ref_policy="last") == cpp.factorize(s.encode(), ref_policy="last")
    assert noLZSS.factorize(s, ref_policy="LAST") == cpp.factorize(s.encode(), ref_policy="last")
    assert noLZSS.factorize_with_info(s, ref_policy="last")["factors"] == cpp.factorize(s.encode(), ref_policy="last")

    p = tmp_path / "in.txt"
    p.write_text(s)
    assert noLZSS.factorize_file(p, ref_policy="last") == cpp.factorize(s.encode(), ref_policy="last")

    out = tmp_path / "out.bin"
    noLZSS.write_factors_binary_file(str(p), out, ref_policy="last")  # first argument is an input *path*
    assert read_bin(out) == as4(cpp.factorize(s.encode(), ref_policy="last"))

    from noLZSS.parallel import parallel_factorize, parallel_factorize_dna_w_rc_to_file
    got = parallel_factorize(s, num_threads=2, ref_policy="last")
    assert [tuple(f) for f in got] == cpp.factorize(s.encode(), ref_policy="last")
    parallel_factorize_dna_w_rc_to_file(s, out, num_threads=2, ref_policy="last")
    assert read_bin(out) == cpp.factorize_dna_w_rc(s.encode(), ref_policy="last")

    assert noLZSS.factorize_w_reference("ACGT", "ACGTACGT", ref_policy="last") == \
        cpp.factorize_w_reference("ACGT", "ACGTACGT", ref_policy="last")


def test_invalid_policy_rejected():
    with pytest.raises(ValueError):
        validate_ref_policy("middle")
    with pytest.raises(TypeError):
        validate_ref_policy(1)
    with pytest.raises(ValueError):
        noLZSS.factorize("ACGT", ref_policy="middle")
    with pytest.raises(ValueError):
        cpp.factorize(b"ACGT", ref_policy="middle")
    with pytest.raises(ValueError):
        cpp.factorize_dna_w_rc(b"ACGT", ref_policy="Last")  # the C++ layer is case-sensitive


def test_batch_factorize_single_file_uses_policy(fasta_file, tmp_path):
    from noLZSS.genomics.batch_factorize import factorize_single_file
    path, seqs = fasta_file
    outputs = {
        "without_reverse_complement": tmp_path / "no.bin",
        "with_reverse_complement": tmp_path / "rc.bin",
    }
    assert all(factorize_single_file(path, outputs, skip_existing=False, ref_policy="last").values())
    prepared_rc, _ = prepared_text(cpp.prepare_multiple_dna_sequences_w_rc(list(seqs.values())))
    assert read_bin(outputs["with_reverse_complement"]) == \
        cpp.factorize_multiple_dna_w_rc(prepared_rc.encode("latin-1"), ref_policy="last")


def test_lsf_job_script_carries_policy(tmp_path):
    from noLZSS.genomics.lsf_batch_factorize import create_job_script, FactorizationMode
    script = tmp_path / "job.sh"
    assert create_job_script(tmp_path / "in.fa", tmp_path / "out.bin", FactorizationMode.WITH_REVERSE_COMPLEMENT,
                             4, script, ref_policy="last")
    text = script.read_text()
    assert 'ref_policy = "last"' in text and "ref_policy=ref_policy" in text
    assert not create_job_script(tmp_path / "in.fa", tmp_path / "out.bin", FactorizationMode.WITH_REVERSE_COMPLEMENT,
                                 4, script, ref_policy="middle")
