# Prior Art

Hydrogen's hidden SO(4) degeneracy structure is established physics and
mathematics. This repository does not claim novelty for that structure.

Pauli's algebraic treatment of the hydrogen spectrum used the enlarged symmetry
associated with angular momentum and the Runge-Lenz vector. Fock's 1935 work
placed the hydrogen problem into a hyperspherical momentum-space framework. In
modern language, these are part of the standard Pauli/Fock hydrogen symmetry
lineage.

This repository currently implements only the shell-level representation
organization:

- scalar S^3 shell labels K;
- the principal-shell identification n = K + 1;
- pre-spin shell dimension n^2;
- branching of V_j tensor V_j* into ordinary angular-momentum sectors V_ell;
- finite-dimensional numerical diagnostics for the branching transform.

It does not implement the complete Fock momentum-space stereographic map, the
physical integration measure, energy-dependent momentum scaling, an inverse
transform to position space, radial Coulomb wavefunctions, or transition
amplitudes.

Any novelty in the active project should be read as software organization,
explicit conventions, reproducible diagnostics, and exposition unless a
separate mathematical or physical result is demonstrated.

## Reference Anchors

The following entries are included only as conservative bibliographic anchors.
They are not page-level citation claims.

- W. Pauli Jr., "Uber das Wasserstoffspektrum vom Standpunkt der neuen Quantenmechanik," Zeitschrift fur Physik 36, 336-363 (1926), DOI 10.1007/BF01450175.
- V. Fock, "Zur Theorie des Wasserstoffatoms," Zeitschrift fur Physik 98, 145-154 (1935), DOI 10.1007/BF01336904.
- D. S. Chavel, Eigenvalues in Riemannian Geometry, Academic Press, 1984.
- S. Helgason, Groups and Geometric Analysis, Academic Press, 1984.
- B. C. Hall, Lie Groups, Lie Algebras, and Representations, 2nd ed., Springer, 2015.
- D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum Theory of Angular Momentum, World Scientific, 1988.
- L. D. Landau and E. M. Lifshitz, Quantum Mechanics: Non-Relativistic Theory, 3rd ed., Pergamon.
