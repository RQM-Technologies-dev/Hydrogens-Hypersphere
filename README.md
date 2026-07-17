[![Hydrogen Bridge v1 CI](https://github.com/RQM-Technologies-dev/Hydrogen/actions/workflows/ci.yml/badge.svg)](https://github.com/RQM-Technologies-dev/Hydrogen/actions/workflows/ci.yml)

# Hydrogen Bridge v1

## RQM Technical Canon v2

This repository is a standard-compatible hydrogen spectral benchmark and
historical bridge implementation. Its `S³` shell identities are established
harmonic mathematics, and its physical energy scale is calibrated from the
conventional Rydberg law. It does not establish new hydrogen physics. See
[RQM_TECHNICAL_CANON_V2.md](RQM_TECHNICAL_CANON_V2.md).

## One-sentence claim
Hydrogen Bridge v1 implements a tested `S³` harmonic organization whose shell
labels reproduce the conventional principal-number and degeneracy identities
when paired with the calibrated Rydberg energy law.

## Central public claim
Scalar harmonics on `S³` provide a standard-compatible computational
organization of the conventional hydrogen shell labels.

## Core result
\[
-\Delta_{S^3}Y_K = K(K+2)Y_K
\]
\[
\hat N = \sqrt{-\Delta_{S^3}+1},\qquad \hat N Y_K = (K+1)Y_K
\]
\[
n = K+1,\qquad \dim\mathcal H_K(S^3) = (K+1)^2 = n^2
\]
\[
H_C = -\frac{\mathrm{Ry}}{-\Delta_{S^3}+1},\qquad H_C\Psi=E\Psi,\qquad E_n=-\frac{\mathrm{Ry}}{n^2}
\]

## Why this matters
Hydrogen Bridge v1 provides executable spectral bookkeeping and diagnostics
for `S³` harmonic shell identities.

- The ordinary hydrogen shell number appears as an S^3 spectral shell number.
- The same shell number controls both the energy denominator and the pre-spin shell degeneracy.
- The calibrated operator \(H_C=-\mathrm{Ry}/(-\Delta_{S^3}+1)\) is a compact
  reparameterization of the conventional `-Ry/n²` law, not a derivation of it
  from geometry alone.
- The angular bridge \(\Pi_K^{\mathrm{ang}}\) maps S^3 shell content to standard hydrogen angular labels.
- The repository includes executable tests and generated reports.

## What is implemented

### Core
- clean S^3 spectral equation (`notes/native_s3_spectral_hydrogen_equation.md`)
- S^3 scalar harmonic shell architecture (`notes/s3_scalar_harmonic_shell_architecture.md`)
- shell table and energy helper (`simulator/hydrogen_shell_simulator.py`)
- low-K numerical \(\Pi_K^{\mathrm{ang}}\) angular bridge and diagnostics (`notes/pi_k_angular_intertwiner.md`, `simulator/s3_s2_intertwiner.py`)
- claims matrix (`docs/claims_matrix.md`)

### Validation
- NIST ASD spectral comparison (`simulator/spectral_comparison.py`, `data/hydrogen_reference_lines.csv`)
- shell-locking numerical diagnostic (`simulator/shell_locking_test.py`)
- \(\Pi_K^{\mathrm{ang}}\) L² compatibility diagnostics and rank-1 angular transition-operator diagnostics (`simulator/angular_operators.py`, `reports/pi_k_l2_diagnostics.csv`, `reports/rank1_transition_diagnostics.csv`)
- report generation (`scripts/generate_reports.py`)
- tests/CI (`tests/`, `.github/workflows/ci.yml`)

### Appendices/support
- closure-action support (`notes/appendix_b_closure_geometry_inverse_square_action.md`)
- Hopf flux projection support (`notes/appendix_c_hopf_flux_projection_and_coulomb_field.md`)
- fine-structure benchmark (`simulator/fine_structure.py`, `reports/h_alpha_fine_structure.csv`)

## What this repository does not claim
- no first-principles Rydberg-constant derivation
- no derivation of \(\kappa\), electric charge, or Maxwell equations
- no full Schrödinger-Coulomb solution derivation
- no native fine-structure derivation yet
- no complete physical unitary \(\Pi_K\) operator yet
- no AGQF, anchor-well, or physical-slice derivation
- no requirement that physical space has an additional `S³` dimension
- no replacement for the conventional Coulomb eigenproblem

## Relationship to Fock hydrogen

Vladimir Fock's standard construction compactifies momentum space on an
energy-dependent `S³` and derives the Coulomb spectrum from the conventional
integral eigenproblem. That full construction is implemented in the sibling
`Hydrogens-Hypersphere` repository. Fock hydrogen is not AGQF and does not
require `s²=2n`.

## Start here
1. `docs/director_summary.md`
2. `notes/s3_scalar_harmonic_shell_architecture.md`
3. `notes/native_s3_spectral_hydrogen_equation.md`
4. `notes/pi_k_angular_intertwiner.md`
5. `docs/claims_matrix.md`
6. `reports/HYDROGEN_BRIDGE_V1_REPORT.md`

For appendices and support layers, see `notes/layer_separation.md`.

## Reviewer map

### Core reviewer path
- `docs/director_summary.md`
- `notes/s3_scalar_harmonic_shell_architecture.md`
- `notes/native_s3_spectral_hydrogen_equation.md`
- `notes/pi_k_angular_intertwiner.md`
- `docs/claims_matrix.md`

### Validation artifacts
- `reports/HYDROGEN_BRIDGE_V1_REPORT.md`
- `reports/shell_table.csv`
- `reports/angular_state_counts.csv`
- `reports/series_comparison.csv`
- `reports/shell_locking_validation.csv`

### Appendices/support
- `notes/layer_separation.md`
- `notes/appendix_coulomb_action_to_s3_operator.md`
- `notes/appendix_b_closure_geometry_inverse_square_action.md`
- `notes/appendix_c_hopf_flux_projection_and_coulomb_field.md`

## Quickstart
```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest
python scripts/generate_reports.py
python scripts/generate_plots.py
```

## Continuous Integration
The CI workflow runs tests and regenerates reports on pushes to `main` and pull requests targeting `main`.

## License
This repository is licensed under the MIT License. See [`LICENSE`](LICENSE).
