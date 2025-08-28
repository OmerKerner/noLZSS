# noLZSS

[![Build Wheels](https://github.com/OmerKerner/noLZSS/actions/workflows/wheels.yml/badge.svg)](https://github.com/OmerKerner/noLZSS/actions/workflows/wheels.yml)

Non-overlapping Lempel–Ziv–Storer–Szymanski factorization (Python bindings over high-performance C++).

## Install (dev)

```bash
pip install -e .
```

## Usage

```python
import noLZSS
factors = noLZSS.factorize(b"abracadabra$")
print(factors)
```

## License

This project is licensed under the BSD 3‑Clause License (see `LICENSE`).

The repository vendors third-party components (notably SDSL v3). Third‑party
license texts and attribution are provided in `THIRD_PARTY_LICENSES.txt`.
