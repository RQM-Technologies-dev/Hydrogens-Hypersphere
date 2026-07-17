# Spectral Energy Model

The package defines a calibrated spectral shell function

```text
H_spec(K) = -Ry / (-Delta_S3 + 1) = -Ry / (K+1)^2.
```

With `n = K + 1`, this reproduces the conventional nonrelativistic Rydberg
shell law

```text
E_n = -Ry / n^2.
```

In this legacy compatibility helper, the Rydberg energy scale is supplied
numerically. The function is useful for keeping shell labels, degeneracy, and
the conventional shell energy denominator in one executable representation. It
is not itself a derivation. The separate
`hydrogen_s3.fock.quantization_from_fock` path now derives the spinless bound
spectrum from the conventional Coulomb Hamiltonian and reconstructs radial
wavefunctions. Neither path derives the Coulomb law, electric charge, Maxwell
equations, fine structure, the Lamb shift, or new measurable hydrogen physics.

The code deliberately avoids naming this function as a derived Coulomb
Hamiltonian. Documentation and reports use `H_spec` or
`calibrated_spectral_energy_from_K`.

For physical hydrogen, `CoulombSystem.physical_hydrogen()` uses the
electron-proton reduced mass. The historical `RYDBERG_EV` compatibility value
is the infinite-nuclear-mass Rydberg scale, so their energies differ by the
standard reduced-mass ratio.
