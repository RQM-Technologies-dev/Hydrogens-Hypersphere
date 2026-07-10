# Conventions

## Basis Ordering

The source basis for shell K is ordered as `(a, b)` with `a` and `b` ascending
magnetic labels in `V_j` and `V_j*`, where `j = K/2`.

The target basis is ordered by increasing `ell`, then ascending `m`:

```text
(0, 0), (1, -1), (1, 0), (1, 1), ...
```

## Dual Representation

The algebraic dual generator is implemented as `-J_a.T` on `V_j*`. The source
generator is constructed independently as

```text
J_a_source = J_a tensor I + I tensor (-J_a.T).
```

The dual basis is identified with a second copy of `V_j` by

```text
e^b -> (-1)^(j-b) |j, -b>.
```

This phase is part of the basis convention.

## Clebsch-Gordan Convention

Clebsch-Gordan coefficients are obtained from
`sympy.physics.wigner.clebsch_gordan`. The wrapper converts integer and
half-integer inputs to exact rational values before evaluation. Numeric complex
values are created only at matrix assembly boundaries.

## Matrix Convention

All finite-dimensional matrices use `complex128`. Adjoint operations use
conjugate transpose, written in code as `.conj().T`.

## Angular Tensor Convention

The angular tensor module uses

```text
<ell_f m_f | T^(k)_q | ell_i m_i>
  = <ell_i m_i; k q | ell_f m_f>
    * <ell_f || T^(k) || ell_i> / sqrt(2 ell_f + 1).
```

The reduced matrix element is supplied by the caller. The package does not
provide radial matrix elements or oscillator strengths.
