# Scope and Claims

The project implements the established Fock construction for conventional spinless Coulomb bound states. Its contribution is a convention-explicit, tested software implementation and exposition, not a new hydrogen theory.

| Statement | Status | Implementation and validation | Limitation |
| --- | --- | --- | --- |
| The unitary Fourier transform gives the Coulomb coefficient `g/(2 pi^2 hbar)`. | standard derivation | `fock.coulomb_kernel`; convention derivation document | Bound states only. |
| Negative-energy physical momentum compactifies on an energy-dependent `S^3` with scale `q=hbar*kappa`. | Fock construction | `fock.stereographic`; round-trip, Jacobian, and chord tests | The sphere is momentum-space geometry, not physical space. |
| Normalized scalar harmonics have eigenvalues `K(K+2)` and shell dimension `(K+1)^2`. | established harmonic analysis | `fock.harmonics`; quadrature and finite-difference tests | Scalar harmonics only. |
| The Coulomb operator eigenvalue is `2 pi^2/(K+1)`. | established Fock result | analytic spectral API plus independent endpoint-safe quadrature | Numerical quadrature is diagnostic, not the production operator. |
| Matching the transformed equation gives `n=K+1`, `kappa_n=1/(n a0)`, and `E_n=-mu g^2/(2 hbar^2 n^2)`. | derived conventional spectrum | `quantization_from_fock`; cross-check against legacy spectrum | No relativistic corrections. |
| Fock harmonics reconstruct normalized momentum states and associated-Laguerre position radial functions. | standard Coulomb eigenfunctions | `fock.bound_states`, `fock.transforms`; 1s, 2s, 2p, and n=3 tests | Oscillatory inverse quadrature has documented finite tolerance. |
| Fixed-shell generators realize `so(4)=su(2)+ direct sum su(2)-`. | Pauli/Fock symmetry | `fock.so4`; commutator, Casimir, and branching compatibility tests | Dimensionless scaled Runge-Lenz vector on bound shells. |
| Existing calibrated shell-energy helpers remain compatible. | compatibility behavior | `hydrogen_s3.spectrum`; reduced-mass-aware cross-check | The legacy `RYDBERG_EV` is the infinite-nuclear-mass scale; physical hydrogen uses reduced mass. |

Not implemented: continuum/scattering states, Dirac hydrogen, spin, fine structure, Lamb shift, QED, external fields, many-electron atoms, transition radial matrix elements, or novel RQM/QSG dynamics. A family of energy-scaled Fock spheres may motivate separate mathematical work, but is not part of the established implementation or evidence for extra physical dimensions.
