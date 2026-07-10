# Changelog

## 0.2.0

- Corrected the scientific scope to a narrow implementation of the standard S^3/SO(4) shell-level organization of nonrelativistic hydrogen bound states.
- Removed unsupported RQM, closure-action, shell-locking, Hopf-Coulomb, and fine-structure material from the active tree.
- Renamed the distribution to `hydrogen-s3` and the import package to `hydrogen_s3`.
- Added independent source and target angular-momentum generator construction for branching-transform diagnostics.
- Replaced the custom Clebsch-Gordan implementation with a cached wrapper around SymPy's standard convention.
- Corrected spectroscopy-medium handling so air-reference residuals are not computed against vacuum predictions.
- Replaced the documentation structure with scope, prior-art, mathematics, conventions, spectral-model, baseline, and limitations pages.
