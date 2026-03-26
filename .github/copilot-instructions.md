# noLZSS Copilot Instructions

## Repository Purpose

`noLZSS` is a hybrid C++/Python library for non-overlapping LZSS factorization, with a strong genomics focus (FASTA processing and reverse-complement-aware factorization).

## Architecture at a Glance

- C++ core algorithms: `src/cpp/`
- Python package/wrappers: `src/noLZSS/`
- Genomics package: `src/noLZSS/genomics/`
- Tests: `tests/`
- Build + packaging config: `CMakeLists.txt`, `pyproject.toml`
- Docs: `docs/`
- Benchmarks: `benchmarks/`

Important coupling points:

- Python wrappers validate inputs, then call `_noLZSS` C++ bindings.
- Reverse-complement flag uses MSB in `ref` (`RC_MASK`).
- Binary factor file metadata is parsed in Python utilities and must stay consistent with C++ struct layout.

## Build & Setup (always run from repo root)

```bash
pip install -e .
```

If extension rebuild is stale:

```bash
pip install -e . --no-build-isolation --force-reinstall
```

## Validate Changes

Preferred full suite:

```bash
python -m pytest -q tests
```

Alternative script runner:

```bash
python tests/run_all_tests.py
```

Use focused tests when changing specific areas (for example `test_genomics_significance.py`, `test_reference_seq.py`, `test_parallel*.py`).

## Coding Conventions in This Repo

- Validate inputs in Python before crossing into C++.
- Use `pathlib.Path` in Python wrappers and convert to `str` only at binding boundaries.
- Preserve binary format compatibility when changing factor or footer structures.
- Keep forward and parallel C++ paths behaviorally aligned (often mirrored logic).
- For large data workflows, prefer file-based APIs and C++ implementations.

## Genomics-Specific Notes

- DNA validation is strict (`A/C/T/G` for DNA-oriented paths).
- Sentinels are used to separate concatenated sequences; do not introduce collisions.
- Reverse complement uses encoded flagging in `ref` and requires careful encode/decode parity.

## Optional Dependencies

- `matplotlib` for plotting (`[plotting]` extra)
- `scipy` improves genomics significance calculations (`[genomics]` extra; fallback exists)
- test dependencies available via `[test]`

## Key Workflows

Benchmarks:

```bash
python benchmarks/fasta_benchmark.py
python benchmarks/fasta_predictor.py benchmarks/fasta_results/trend_parameters.pkl --size 1000000
```

Docs (similar to CI):

```bash
pip install sphinx sphinx-rtd-theme myst-parser breathe exhale
pip install -e .
cd docs && sphinx-build -b html . _build/html
```

## Practical Change Checklist for Agents

When modifying APIs or algorithms, also check:

1. Python exports (`__init__.py`) and wrapper signatures
2. C++ bindings (`src/cpp/bindings.cpp`)
3. Tests under `tests/` for affected feature area
4. Relevant docs in `README.md` and `docs/*.md`
5. Packaging/install lists in `CMakeLists.txt` when adding Python modules
