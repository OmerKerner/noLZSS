# noLZSS

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

## sdsl-lite (Succinct Data Structure Library 3.0) vendoring
By default tries system sdsl-lite (or the Succinct Data Structure Library 3.0 fork) then fetches pinned commit. Override:
 
```bash
pip install . --config-settings=cmake.args="-DNOLZSS_FORCE_VENDOR_SDSL=ON"
```

## License

GPL-3.0-or-later (includes sdsl-lite notice). See `LICENCE` and `LICENSE-GPL-3.0`.
