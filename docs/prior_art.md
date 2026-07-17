# Prior Art

Hydrogen's hidden `SO(4)` degeneracy and momentum-space `S^3` formulation are established physics and mathematics. Pauli derived the bound spectrum algebraically using angular momentum and the Runge-Lenz vector. Fock gave the momentum-space stereographic formulation in 1935. Bargmann further developed the dynamical-group interpretation. The implementation in this repository is not a new derivation invented by RQM Technologies.

The repository's contribution is software organization: explicit unitary Fourier constants, normalized harmonics and momentum states, independent numerical checks, inverse-transform comparisons, and compatibility with its pre-existing finite-dimensional branching machinery.

## Primary historical anchors

- W. Pauli Jr., “Über das Wasserstoffspektrum vom Standpunkt der neuen Quantenmechanik,” *Zeitschrift für Physik* **36**, 336–363 (1926), DOI: 10.1007/BF01450175.
- V. Fock, “Zur Theorie des Wasserstoffatoms,” *Zeitschrift für Physik* **98**, 145–154 (1935), DOI: 10.1007/BF01336904.
- V. Bargmann, “Zur Theorie des Wasserstoffatoms: Bemerkungen zur gleichnamigen Arbeit von V. Fock,” *Zeitschrift für Physik* **99**, 576–582 (1936), DOI: 10.1007/BF01798157.

Modern spectral-geometry, special-function, and representation-theory references remain useful for conventions, but the formulas implemented here belong to the established Pauli/Fock/Bargmann lineage.
