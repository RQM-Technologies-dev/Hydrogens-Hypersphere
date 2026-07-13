# Fock Bound-State Construction

## 1. Coulomb Schrödinger equation

For reduced mass `mu`, positive coupling `g`, and `hbar>0`,

```text
H = p^2/(2 mu) - g/r,       H psi = E psi.
```

This document concerns only `E<0` spinless bound states.

## 2. Fourier convention

The package uses physical momentum and the unitary pair

```text
phi(p) = (2 pi hbar)^(-3/2) integral exp(-i p.r/hbar) psi(r) d^3r,
psi(r) = (2 pi hbar)^(-3/2) integral exp(+i p.r/hbar) phi(p) d^3p.
```

Consequently Parseval normalization is identical in the two spaces. Since

```text
integral exp(-i q.r/hbar) / r d^3r = 4 pi hbar^2/q^2,
```

the transform of a product contributes a second factor `(2 pi hbar)^(-3/2)`. No constants are imported from a different Fourier convention.

## 3. Momentum-space integral equation

The resulting equation is

```text
(p^2/(2 mu) - E) phi(p)
  = g/(2 pi^2 hbar) integral phi(p')/|p-p'|^2 d^3p'.
```

In atomic units the coefficient is `1/(2 pi^2)`.

## 4. Bound-state parameter

Write

```text
E = -hbar^2 kappa^2/(2 mu),       q = hbar kappa.
```

`kappa` has inverse-length units, while `q` has the physical-momentum units of `p`. The stereographic module calls its same-dimension scale argument `kappa` for compatibility with the standard geometry formula; callers using SI physical momentum pass `q`. In atomic units `q=kappa`.

## 5. Stereographic compactification

For each negative energy,

```text
u = 2 q p/(p^2+q^2),       u4 = (p^2-q^2)/(p^2+q^2),
p = q u/(1-u4),            |u|^2+u4^2=1.
```

Thus each energy uses a different compactification of momentum space. Fock's `S^3` is not an additional dimension of ordinary physical space.

## 6. Measure and kernel transformation

Direct differentiation gives

```text
dOmega3 = (2q/(p^2+q^2))^3 d^3p,
|U-U'|^2 = 4q^2 |p-p'|^2 / ((p^2+q^2)(p'^2+q^2)).
```

Both identities are checked independently in finite differences and direct Euclidean coordinates.

## 7. Wavefunction rescaling

The convention-fixed Fock lift is

```text
Psi(U) = (p^2+q^2)^2 phi(p)/(4 q^(5/2)),
phi(p) = 4 q^(5/2) Psi(U)/(p^2+q^2)^2.
```

This is the equation-diagonalizing lift, not the pointwise unitary Jacobian lift. For a normalized single-shell harmonic, substitution of the measure gives

```text
integral |phi|^2 d^3p = integral (1-u4)|Psi|^2 dOmega3 = 1.
```

The last equality follows because `u4` maps degree `K` harmonics only into the adjacent `K-1` and `K+1` shells, so its diagonal shell expectation vanishes. The factor is also fixed by the normalized 1s transform; no amplitude fitting is used.

## 8. Harmonics and kernel eigenvalues

The package implements

```text
Y_Klm = N_Kl sin(chi)^l C_(K-l)^(l+1)(cos chi) Y_lm(theta,phi),
N_Kl = 2^l l! sqrt(2(K+1)(K-l)!/(pi(K+l+1)!)).
```

They are normalized on `S^3`, obey `-Delta Y=K(K+2)Y`, and diagonalize

```text
(Kcal f)(U) = integral f(U')/|U-U'|^2 dOmega3',
Kcal Y_Klm = 2 pi^2/(K+1) Y_Klm.
```

The production path applies this eigenvalue spectrally. A separate Gauss-Legendre diagnostic evaluates zonal harmonics at the south pole; rewriting the integrable endpoint factor as `(1-cos chi)/2` avoids sampling the singular point.

## 9. Quantization

Substitution of the lift, Jacobian, and chord identity gives

```text
Psi = mu g/(2 pi^2 hbar q) Kcal Psi.
```

Matching the harmonic eigenvalue yields

```text
q_n = mu g/(hbar n),       kappa_n = mu g/(hbar^2 n) = 1/(n a0),
n=K+1,                     a0=hbar^2/(mu g),
E_n = -q_n^2/(2mu) = -mu g^2/(2 hbar^2 n^2).
```

This route derives the spectrum from the Coulomb equation rather than calling the legacy calibrated energy helper.

## 10. Momentum-space states

With `K=n-1`, the normalized state is the inverse lift of `Y_Klm`, multiplied by the fixed forward-transform phase `(-i)^l`. The radial factor is

```text
F_nl(p) = (-i)^l 4 q_n^(5/2) Pi_(n-1,l)(chi)/(p^2+q_n^2)^2,
chi = 2 atan2(q_n,p).
```

The implementation supports scalar and batch evaluation and validates 1s, 2s, 2p, and `n=3` states without arbitrary rescaling.

## 11. Inverse transform to ordinary orbitals

Angular integration of the inverse Fourier transform gives

```text
R_nl(r) = 4 pi i^l/(2 pi hbar)^(3/2)
          integral_0^infinity p^2 F_nl(p) j_l(pr/hbar) dp.
```

Adaptive quadrature increases a finite momentum window until cutoff doubling meets the requested tolerance. It is compared with

```text
sqrt((2/(n a0))^3 (n-l-1)!/(2n(n+l)!))
exp(-rho/2) rho^l L_(n-l-1)^(2l+1)(rho),   rho=2r/(n a0).
```

The known `(-i)^l`/`i^l` phases cancel. The comparison API reports maximum absolute error, relative radial `L2` error, sampled normalization error, and phase alignment.

## 12. SO(4) and Runge-Lenz symmetry

On shell `K`, two commuting spin-`j=K/2` factors are built from the existing angular-momentum matrices. With

```text
L = J+ + J-,       M = J+ - J-,       J+/- = (L +/- M)/2,
```

the dimensionless scaled Runge-Lenz generators obey `[L,L]=i epsilon L`, `[L,M]=i epsilon M`, and `[M,M]=i epsilon L`. The checks include `L.M=0`, `L^2+M^2=K(K+2)`, dimension `(2j+1)^2=n^2`, and the existing Clebsch-Gordan branching to `l=0,...,K`.

## 13. Numerical validation

Exact finite-dimensional identities use tolerances near `1e-12`. Harmonic and kernel quadratures use method-appropriate tolerances near `1e-11`. Momentum and analytic radial normalizations are semi-infinite quadratures. Default inverse-transform tests use a `2e-7` integration request and require pointwise errors below `8e-7`; tighter tolerances require more oscillatory-tail work.

## 14. Limitations

The complete implementation here means the conventional nonrelativistic spinless Coulomb **bound-state** chain. It excludes continuum states, relativity, spin, fine structure, the Lamb shift, QED, external fields, and many-electron atoms. It is not proof of an extra physical spatial dimension, not proof of quaternion-valued physical wavefunctions, and not novel RQM/QSG dynamics. Energy-scaled families of Fock spheres may be studied separately but are not an established physical-space claim.
