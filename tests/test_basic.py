import noLZSS

def check_invariants(text: bytes):
    if not text.endswith(b"$"):
        text += b"$"
    factors = noLZSS.factorize(text)
    n = len(text) - 1
    covered = 0
    prev_end = 0
    for start, length in factors:
        assert 0 <= start < n
        assert length > 0
        assert start >= prev_end
        prev_end = start + length
        covered += length
    assert covered == n
    return factors

def test_repeated():
    check_invariants(b"aaaaa$")

def test_mixed():
    check_invariants(b"abracadabra$")

def test_short():
    check_invariants(b"a$")

def test_version():
    assert hasattr(noLZSS, "factorize")
