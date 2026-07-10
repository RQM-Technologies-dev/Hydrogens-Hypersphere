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

The Rydberg energy scale is supplied numerically. The function is useful for
keeping shell labels, degeneracy, and the conventional shell energy denominator
in one executable representation. It is not a derivation of Ry, the Coulomb
interaction, electric charge, Maxwell equations, radial wavefunctions, fine
structure, the Lamb shift, or new measurable hydrogen physics.

The code deliberately avoids naming this function as a derived Coulomb
Hamiltonian. Documentation and reports use `H_spec` or
`calibrated_spectral_energy_from_K`.
