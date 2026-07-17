# Hydrogen Fock S^3 Bound States

[![CI](https://github.com/RQM-Technologies-dev/Hydrogens-Hypersphere/actions/workflows/ci.yml/badge.svg)](https://github.com/RQM-Technologies-dev/Hydrogens-Hypersphere/actions/workflows/ci.yml)

## RQM Technical Canon v2

This repository is a standard-compatible hydrogen spectral benchmark and
Fock-construction implementation. Its `S³` shell identities are established
harmonic mathematics, and its physical energy scale follows the conventional
Coulomb problem. It does not establish new hydrogen physics. See
[RQM_TECHNICAL_CANON_V2.md](RQM_TECHNICAL_CANON_V2.md).

This repository implements the standard nonrelativistic, spinless, bound-state Fock construction of the hydrogen atom. Starting from the Coulomb Hamiltonian, it derives the unitary-convention momentum integral equation, compactifies each negative-energy momentum space on an energy-dependent unit `S^3`, diagonalizes the Coulomb operator with normalized hyperspherical harmonics, derives the hydrogen spectrum, reconstructs normalized momentum states and ordinary associated-Laguerre radial orbitals, and verifies the fixed-shell `SO(4)` algebra.

The implementation supports atomic units and physical SI parameters with the electron-proton reduced mass. The legacy shell, branching, angular-tensor, and spectroscopy APIs remain available.

## Canon v2 boundary

Fock's hydrogen hypersphere is established standard momentum-space geometry.
It is **not AGQF**, does not use an anchor potential, does not require
`s^2=2n`, and is not evidence of an additional physical spatial dimension.
The principal quantum number and spectrum follow from the conventional Coulomb
eigenproblem. See [RQM_TECHNICAL_CANON_V2.md](RQM_TECHNICAL_CANON_V2.md).

## Quickstart

```python
import numpy as np

from hydrogen_s3.fock import CoulombSystem, momentum_wavefunction, quantization_from_fock

system = CoulombSystem.atomic_units()
state = quantization_from_fock(K=1, system=system)
phi_210 = momentum_wavefunction(2, 1, 0, np.array([0.1, 0.2, 0.3]), system)
print(state.kappa, state.energy, phi_210)
```

The focused implementation lives in `hydrogen_s3.fock`:

- conventions and Coulomb systems;
- stereographic geometry, measure, and chord identity;
- normalized complex `S^3` harmonics;
- Coulomb/Fock operators and quantization;
- normalized momentum and position radial states;
- inverse spherical-Bessel transforms;
- fixed-shell angular-momentum/Runge-Lenz generators;
- deterministic diagnostics and generated Markdown/JSON reports.

See [docs/fock_bound_state_construction.md](docs/fock_bound_state_construction.md) for the convention reconciliation and full derivation, [docs/scope_and_claims.md](docs/scope_and_claims.md) for the claim matrix, and [docs/prior_art.md](docs/prior_art.md) for historical attribution.

## Scope

Fock's `S^3` is an energy-dependent compactification of momentum space. It is not evidence for an extra physical spatial dimension or a quaternion-valued physical wavefunction, and this standard construction was not invented by RQM Technologies. Continuum/scattering states, relativity, spin, fine structure, the Lamb shift, QED, external fields, and many-electron atoms are outside scope.

## Development

```bash
python -m pip install -e ".[dev]"
ruff check .
ruff format --check .
mypy src
pytest -q
python scripts/generate_report.py
python -m build
```

Reports are written under `build/reports/` and are not committed.
