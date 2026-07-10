# Hydrogen S^3 Representation

[![CI](https://github.com/RQM-Technologies-dev/Hydrogens-Hypersphere/actions/workflows/ci.yml/badge.svg)](https://github.com/RQM-Technologies-dev/Hydrogens-Hypersphere/actions/workflows/ci.yml)

This repository provides an executable implementation of the standard S^3/SO(4)
representation-theoretic organization of nonrelativistic hydrogen bound-state
shells. It reproduces shell labels, pre-spin n^2 degeneracy, and branching into
ordinary angular-momentum sectors.

The project does not derive the Coulomb interaction, the Rydberg constant,
charge, Maxwell's equations, radial hydrogen wavefunctions, fine structure, or
new hydrogen physics. Its spectral energy function is a calibrated packaging of
the conventional E_n = -Ry/n^2 shell law.

## Mathematical Core

The retained construction is shell-level:

- scalar S^3 shell eigenvalues `K(K+2)`;
- shifted eigenvalues `(K+1)^2`;
- principal labels `n = K + 1`;
- shell dimensions `(K+1)^2 = n^2`;
- representation organization `H_K(S^3) ~= V_j tensor V_j*`, with `j = K/2`;
- branching into `direct sum_{ell=0}^{K} V_ell`.

The package constructs the branching transform `U_K` using SymPy
Clebsch-Gordan coefficients and verifies it against independently constructed
source and target angular-momentum generators.

## What Is Implemented

- `hydrogen_s3.spectrum`: S^3 shell identities and calibrated Rydberg shell energies.
- `hydrogen_s3.angular_momentum`: finite-dimensional spin matrices for integer and half-integer `j`.
- `hydrogen_s3.clebsch_gordan`: cached exact-wrapper access to SymPy Clebsch-Gordan coefficients.
- `hydrogen_s3.branching`: shell-wise unitary branching transform and diagnostics.
- `hydrogen_s3.angular_tensors`: angular Wigner-Eckart factors with explicit reduced matrix elements.
- `hydrogen_s3.reference_spectroscopy`: conventional Rydberg baseline with air/vacuum medium handling.

## What Is Not Claimed

This repository does not claim novelty for hydrogen's S^3/SO(4) organization.
It does not implement the full Pauli/Fock hydrogen construction, physical
measure, momentum-space stereographic transform, inverse transform, radial
Coulomb wavefunctions, oscillator strengths, fine structure, or experimental
predictions.

The branching transform is a unitary representation transform, not a physical
spatial projection.

See [docs/scope_and_claims.md](docs/scope_and_claims.md) for the full scope
matrix.

## Prior-Art Notice

Hydrogen's hidden SO(4) degeneracy and S^3/hyperspherical descriptions have
historical precedent, including Pauli's algebraic treatment and Fock's 1935
momentum-space formulation. This project is a reproducible software
organization of the shell-level representation structure.

See [docs/prior_art.md](docs/prior_art.md).

## Installation

```bash
python -m pip install -e ".[dev]"
```

Python 3.10 or newer is required.

## Quickstart

```python
from hydrogen_s3.branching import branching_diagnostics
from hydrogen_s3.spectrum import angular_labels, calibrated_spectral_energy_from_K

print(angular_labels(2))
print(calibrated_spectral_energy_from_K(2))
print(branching_diagnostics(2))
```

## Testing and Reports

```bash
ruff check .
ruff format --check .
mypy src
pytest -q
python -m build
python scripts/generate_report.py
```

The report is generated under `build/reports/` and is not committed.
