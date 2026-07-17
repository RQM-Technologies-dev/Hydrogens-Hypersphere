# Changelog

## 0.3.0

- Added the complete standard nonrelativistic spinless Coulomb bound-state Fock chain under `hydrogen_s3.fock`.
- Added convention-explicit stereographic geometry, normalized `S^3` harmonics, Coulomb operator diagonalization, spectrum derivation, normalized momentum states, inverse Hankel reconstruction, and fixed-shell SO(4) generators.
- Added deterministic focused validation for geometry, harmonics, kernel eigenvalues, quantization, low hydrogen states, inverse transforms, and branching-compatible SO(4) algebra.
- Extended report generation with Fock diagnostics and machine-readable JSON.
- Updated scientific scope and historical attribution; continuum, relativistic, spin, QED, and many-electron physics remain outside scope.

## 0.2.0

- Corrected the scientific scope to a narrow implementation of the standard S^3/SO(4) shell-level organization of nonrelativistic hydrogen bound states.
- Removed unsupported RQM, closure-action, shell-locking, Hopf-Coulomb, and fine-structure material from the active tree.
- Renamed the distribution to `hydrogen-s3` and the import package to `hydrogen_s3`.
- Added independent source and target angular-momentum generator construction for branching-transform diagnostics.
- Replaced the custom Clebsch-Gordan implementation with a cached wrapper around SymPy's standard convention.
- Corrected spectroscopy-medium handling so air-reference residuals are not computed against vacuum predictions.
- Replaced the documentation structure with scope, prior-art, mathematics, conventions, spectral-model, baseline, and limitations pages.
