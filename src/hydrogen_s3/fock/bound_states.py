"""Normalized Fock momentum states and conventional position radial states."""

from __future__ import annotations

from math import exp, lgamma

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import eval_genlaguerre

from hydrogen_s3.fock.conventions import CoulombSystem
from hydrogen_s3.fock.harmonics import hyperspherical_radial
from hydrogen_s3.fock.stereographic import hyperspherical_angles, momentum_to_s3

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


def quantum_numbers(n: int, ell: int, m: int = 0) -> tuple[int, int, int]:
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("n must be a positive integer")
    if not isinstance(ell, int) or isinstance(ell, bool) or not 0 <= ell < n:
        raise ValueError("ell must satisfy 0 <= ell < n")
    if not isinstance(m, int) or isinstance(m, bool) or abs(m) > ell:
        raise ValueError("m must satisfy -ell <= m <= ell")
    return n, ell, m


def fock_lift_factor(p_magnitude: ArrayLike, momentum_scale: float) -> FloatArray:
    """Return (p^2+q^2)^2/(4 q^(5/2)) in Phi=lift*phi."""

    q = float(momentum_scale)
    if not np.isfinite(q) or q <= 0.0:
        raise ValueError("momentum_scale must be positive and finite")
    p = np.asarray(p_magnitude, dtype=np.float64)
    if np.any(~np.isfinite(p)) or np.any(p < 0.0):
        raise ValueError("p_magnitude must be non-negative and finite")
    return np.asarray((p * p + q * q) ** 2 / (4.0 * q**2.5), dtype=np.float64)


def momentum_radial(n: int, ell: int, p: ArrayLike, system: CoulombSystem | None = None) -> ComplexArray:
    """Return F_nell(p) where phi_nellm=F_nell(p)Y_ellm(p-hat).

    The phase (-i)^ell is fixed by the unitary forward Fourier transform.
    ``p`` is physical momentum (or dimensionless atomic-unit momentum).
    """

    quantum_numbers(n, ell)
    active = CoulombSystem.atomic_units() if system is None else system
    momentum = np.asarray(p, dtype=np.float64)
    if np.any(~np.isfinite(momentum)) or np.any(momentum < 0.0):
        raise ValueError("p must be non-negative and finite")
    q = active.momentum_scale_n(n)
    chi = 2.0 * np.arctan2(q, momentum)
    pi_factor = hyperspherical_radial(n - 1, ell, chi)
    amplitude = 4.0 * q**2.5 * pi_factor / (momentum * momentum + q * q) ** 2
    return np.asarray((-1.0j) ** ell * amplitude, dtype=np.complex128)


def momentum_wavefunction(n: int, ell: int, m: int, p: ArrayLike, system: CoulombSystem | None = None) -> ComplexArray:
    """Evaluate normalized phi_nellm for a final input shape (...,3)."""

    quantum_numbers(n, ell, m)
    vector = np.asarray(p, dtype=np.float64)
    if vector.ndim == 0 or vector.shape[-1] != 3 or np.any(~np.isfinite(vector)):
        raise ValueError("p must be finite with final dimension 3")
    active = CoulombSystem.atomic_units() if system is None else system
    magnitude = np.linalg.norm(vector, axis=-1)
    u = momentum_to_s3(vector, active.momentum_scale_n(n))
    _chi, theta, phi = hyperspherical_angles(u)
    from scipy.special import sph_harm_y

    angular = sph_harm_y(ell, m, theta, phi)
    return np.asarray(momentum_radial(n, ell, magnitude, active) * angular, dtype=np.complex128)


def analytic_position_radial(n: int, ell: int, r: ArrayLike, system: CoulombSystem | None = None) -> FloatArray:
    """Return the conventional normalized associated-Laguerre R_nell(r)."""

    quantum_numbers(n, ell)
    active = CoulombSystem.atomic_units() if system is None else system
    radius = np.asarray(r, dtype=np.float64)
    if np.any(~np.isfinite(radius)) or np.any(radius < 0.0):
        raise ValueError("r must be non-negative and finite")
    a0 = active.bohr_radius
    rho = 2.0 * radius / (float(n) * a0)
    log_normalization = 1.5 * np.log(2.0 / (float(n) * a0)) + 0.5 * (
        lgamma(n - ell) - np.log(2.0 * n) - lgamma(n + ell + 1)
    )
    values = exp(log_normalization) * np.exp(-rho / 2.0) * rho**ell * eval_genlaguerre(n - ell - 1, 2 * ell + 1, rho)
    return np.asarray(values, dtype=np.float64)
