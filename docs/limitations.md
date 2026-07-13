# Limitations

This project is intentionally limited to the standard nonrelativistic, spinless Coulomb bound-state problem. It derives the bound spectrum and radial eigenfunctions from that conventional Hamiltonian through Fock's momentum-space construction; it does not derive the Coulomb law itself, electric charge, or Maxwell's equations.

It does not include:

- continuum or scattering states;
- Dirac hydrogen, spin, or relativity;
- fine or hyperfine structure;
- the Lamb shift or QED;
- external fields;
- many-electron atoms;
- radial transition matrix elements or oscillator strengths;
- new measurable RQM/QSG hydrogen dynamics.

The adaptive inverse Hankel transform is a numerical validation path. Oscillatory-tail convergence becomes more expensive at small radii and tight tolerances; failures raise `TransformConvergenceError` rather than returning a silent low-accuracy result.

Fock's unit `S^3` is an energy-dependent compactification of momentum space. It is not a physical dimensional projection, evidence for a fourth spatial dimension, or evidence for a quaternion-valued physical wavefunction. The fixed-shell branching transform is a unitary representation transform, not a measurement map.
