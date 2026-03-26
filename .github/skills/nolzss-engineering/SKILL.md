---
name: nolzss-engineering
description: 'Work effectively in the noLZSS C++/Python hybrid genomics factorization codebase. Use when modifying LZSS algorithms, Python bindings/wrappers, genomics FASTA workflows, binary factor formats, tests, docs, or build/packaging logic.'
---

# noLZSS Engineering Skill

Use this skill for repository-aware implementation and review tasks in `noLZSS`.

## When to use this skill

- Changing factorization behavior in `src/cpp/`
- Updating Python wrappers in `src/noLZSS/`
- Touching genomics processing (`src/noLZSS/genomics/`)
- Modifying binary factor I/O or metadata parsing
- Updating packaging/build behavior (`CMakeLists.txt`, `pyproject.toml`)
- Validating with repository tests and docs

## Core repository map

- C++ core algorithms: `src/cpp/`
- pybind11 bindings: `src/cpp/bindings.cpp`
- Python wrappers/utilities: `src/noLZSS/core.py`, `src/noLZSS/utils.py`
- Genomics layer: `src/noLZSS/genomics/`
- Tests: `tests/`
- Docs: `README.md`, `docs/`
- Build/package: `CMakeLists.txt`, `pyproject.toml`

## API family map (use this instead of memorizing every function)

### 1) Core factorization (`noLZSS` top-level via `core.py`)

- In-memory text:
  - `factorize(data)`
  - `count_factors(data)`
  - `factorize_with_info(data)`
- File-based text:
  - `factorize_file(path, reserve_hint=0)`
  - `count_factors_file(path)`
- Binary output:
  - `write_factors_binary_file(data, output_path)`
- Reference-guided (no reverse complement):
  - `factorize_w_reference(reference_seq, target_seq)`
  - `factorize_w_reference_file(reference_seq, target_seq, output_path)`

### 2) Utilities (`noLZSS.utils`)

- Input and alphabet:
  - `validate_input(data)`
  - `analyze_alphabet(data)`
- Binary readers:
  - `read_factors_binary_file(path)`
  - `read_binary_file_metadata(path)`
  - `read_factors_binary_file_with_metadata(path)`
- Plot helper:
  - `plot_factor_lengths(...)`

### 3) Genomics APIs (`noLZSS.genomics`)

- DNA/RC factorization:
  - `factorize_dna_w_rc`, `factorize_file_dna_w_rc`
  - `count_factors_dna_w_rc`, `count_factors_file_dna_w_rc`
  - `write_factors_binary_file_dna_w_rc`
- Multi-sequence DNA:
  - `factorize_multiple_dna_w_rc`, `factorize_file_multiple_dna_w_rc`
  - `count_factors_multiple_dna_w_rc`, `count_factors_file_multiple_dna_w_rc`
  - `write_factors_binary_file_multiple_dna_w_rc`
  - `factorize_fasta_multiple_dna_w_rc`
  - `prepare_multiple_dna_sequences_w_rc`
- FASTA and sequence helpers:
  - `read_nucleotide_fasta`, `read_protein_fasta`, `read_fasta_auto`
  - sequence validators from `genomics.sequences`
- Significance/statistics:
  - `calculate_factor_length_threshold`
  - `infer_length_significance`
  - `extract_factor_lengths`
  - `plot_significance_analysis`

### 4) Parallel APIs (`noLZSS.parallel`)

- Generic parallel:
  - `parallel_factorize_to_file`, `parallel_factorize_file_to_file`, `parallel_factorize`
- DNA+RC parallel:
  - `parallel_factorize_dna_w_rc_to_file`
  - `parallel_factorize_file_dna_w_rc_to_file`

## Function selection guide (decision matrix)

- Need fastest path for very large input?
  - Prefer file-based or parallel APIs over loading everything in memory.
- Need only number of factors?
  - Use `count_*` functions instead of full `factorize_*`.
- Need durable output for later analysis?
  - Use `write_*_binary_file` variants.
- Working with DNA and reverse complements?
  - Use genomics `*_dna_w_rc` family.
- Processing FASTA directly?
  - Use `noLZSS.genomics` FASTA functions, especially `factorize_fasta_multiple_dna_w_rc`.
- Need to inspect binary outputs quickly?
  - Use metadata readers before loading all factors.
- Need generic reference-target matching without RC?
  - Use `factorize_w_reference*` family.

## Working rules

1. Preserve C++/Python behavioral parity.
2. Validate inputs in Python wrappers before C++ calls.
3. Keep reverse-complement encoding/decoding (`RC_MASK`) consistent across layers.
4. Preserve binary format compatibility (footer layout, metadata encoding).
5. When adding Python module files, ensure packaging/install logic includes them.
6. Update tests and docs for behavior changes.

## Build and test workflow

From repository root:

```bash
pip install -e .
python -m pytest -q tests
```

If extension appears stale:

```bash
pip install -e . --no-build-isolation --force-reinstall
```

Optional targeted runs:

```bash
python tests/test_genomics.py
python tests/test_reference_seq.py
python tests/test_parallel.py
```

## Genomics-specific cautions

- DNA-oriented functions expect strict nucleotide alphabets.
- Sentinel characters separate sequences; avoid introducing collisions.
- Reverse-complement factor references encode a flag in `ref` MSB.
- Significance analysis can use optional SciPy for exact bounds; fallback behavior should remain explicit and documented.

## Skill authoring guidance for library repositories

- Do not dump every function signature; that hurts discoverability and goes stale.
- Prefer API families + selection rules + caveats.
- Include representative entry points and where they live.
- Link to canonical docs/tests for full detail.
- Keep this file focused on choosing the right API and avoiding mistakes.

## Delivery checklist

- [ ] Change implemented in the correct layer(s)
- [ ] Bindings/wrapper exports updated if needed
- [ ] Tests updated and executed
- [ ] Docs updated for user-visible behavior
- [ ] Packaging/install behavior verified for new files
