# Scope and Claims

This project implements a narrow shell-level representation structure for
nonrelativistic hydrogen bound states. It does not introduce a new physical
hydrogen theory.

| Statement | Status | Implementation | Test | External input | Limitation |
| --- | --- | --- | --- | --- | --- |
| Scalar harmonics on S^3 have Laplacian eigenvalue K(K+2). | established mathematical fact | `hydrogen_s3.spectrum.s3_laplacian_eigenvalue` | `tests/test_spectrum.py` | standard spectral geometry | Shell-level scalar harmonics only. |
| The shifted eigenvalue is (K+1)^2. | representation-theoretic derivation | `shifted_laplacian_eigenvalue` | `tests/test_spectrum.py` | S^3 eigenvalue formula | Does not derive a physical Hamiltonian. |
| Principal shell labels are represented by n = K + 1. | definition | `shell_number` | `tests/test_spectrum.py` | conventional hydrogen labels | Identification is shell-level. |
| Pre-spin shell dimension is n^2. | representation-theoretic derivation | `shell_dimension`, `angular_labels` | `tests/test_spectrum.py` | S^3 shell dimension and S^2 angular sectors | Spin, fine structure, and radial states are outside scope. |
| The angular content branches as V_j tensor V_j* to direct sum of V_ell for ell=0..K. | representation-theoretic derivation | `hydrogen_s3.branching.branching_transform` | `tests/test_branching.py` | SU(2)/SO(4) representation theory | It is a unitary basis transform, not a physical projection map. |
| H_spec(K) = -Ry/(K+1)^2 reproduces E_n = -Ry/n^2. | calibrated conventional model | `calibrated_spectral_energy_from_K`, `rydberg_shell_energy` | `tests/test_spectrum.py` | supplied Rydberg energy scale | Ry and the Coulomb interaction are not derived. |
| Generator intertwining holds for the constructed branch transform. | numerical consistency check | `source_generators`, `target_generators`, `branching_diagnostics` | `tests/test_branching.py` | finite-dimensional spin matrices | Low-dimensional numerical check through K=8. |
| Clebsch-Gordan coefficients follow the SymPy convention. | definition | `hydrogen_s3.clebsch_gordan` | `tests/test_clebsch_gordan.py` | SymPy Wigner implementation | Convention is delegated and documented; no independent tables are derived. |
| Rank-k angular tensor matrices can be assembled from Wigner-Eckart angular factors. | definition | `hydrogen_s3.angular_tensors` | `tests/test_angular_tensors.py` | angular momentum algebra | Radial matrix elements and oscillator strengths are not supplied. |
| The NIST line CSV can be used as a conventional Rydberg baseline. | conventional hydrogen reference calculation | `hydrogen_s3.reference_spectroscopy` | `tests/test_reference_spectroscopy.py` | raw NIST wavelengths and media in `data/` | Residuals are suppressed unless the reference medium matches the vacuum calculation. |
| Implementing Fock's complete momentum-space construction is future work. | future work | documentation only | not applicable | Fock/Pauli hydrogen symmetry literature | Current code does not implement the full stereographic map, physical measure, inverse transform, or radial Coulomb wavefunctions. |
