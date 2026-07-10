# Mathematics

For scalar harmonics on S^3, the shell index K has Laplace-Beltrami eigenvalue

```text
-Delta_S3 Y_K = K(K+2) Y_K.
```

The shifted operator gives

```text
(-Delta_S3 + 1) Y_K = (K+1)^2 Y_K.
```

The project uses the shell label

```text
n = K + 1.
```

At fixed K, scalar harmonics on S^3 may be organized as

```text
H_K(S^3) ~= V_j tensor V_j*,  j = K/2.
```

Under the diagonal SU(2) action,

```text
V_j tensor V_j* ~= direct sum_{ell=0}^{K} V_ell.
```

The dimension identity is

```text
dim H_K(S^3) = (K+1)^2 = sum_{ell=0}^{K} (2 ell + 1).
```

The implementation verifies this finite-dimensional statement by constructing a
unitary branching transform U_K and independent source and target generators.
For each generator J_a with a in {x, y, z}, the numerical check is

```text
U_K J_a_source U_K^dagger ~= J_a_target.
```

This is a representation transform. It is not treated as an experimentally
derived spatial projection or measurement map.
