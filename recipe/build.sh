#!/usr/bin/env bash
set -euo pipefail

# Build wheel and install it into the build environment
# Avoid upgrading pip from PyPI inside the isolated conda-build env (network may be restricted).
python -m pip wheel . -w dist
python -m pip install dist/*.whl --no-deps -vv
