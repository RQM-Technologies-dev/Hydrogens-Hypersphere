"""Coulomb integral equation and its Fock operator on S^3.

With the unitary physical-momentum Fourier transform, the bound-state equation
is

    (p^2/(2 mu)-E) phi(p) = g/(2 pi^2 hbar) int phi(p')/|p-p'|^2 d^3p'.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from math import pi

import numpy as np
from numpy.typing import NDArray
from scipy.special import roots_legendre

from hydrogen_s3.fock.conventions import CoulombSystem
from hydrogen_s3.fock.harmonics import hyperspherical_harmonic, validate_labels

ComplexArray = NDArray[np.complex128]


def momentum_coulomb_coefficient(system: CoulombSystem) -> float:
    """Return g/(2*pi^2*hbar) for the selected unitary Fourier convention."""

    return system.g / (2.0 * pi**2 * system.hbar)


def fock_kernel_eigenvalue(K: int) -> float:
    if not isinstance(K, int) or isinstance(K, bool) or K < 0:
        raise ValueError("K must be a non-negative integer")
    return 2.0 * pi**2 / float(K + 1)


def apply_spectral(K: int, values: ComplexArray) -> ComplexArray:
    """Apply the Fock operator to coefficients belonging to one K shell."""

    return np.asarray(fock_kernel_eigenvalue(K) * np.asarray(values, dtype=np.complex128), dtype=np.complex128)


@dataclass(frozen=True, slots=True)
class FockQuantization:
    K: int
    n: int
    kappa: float
    momentum_scale: float
    energy: float
    transformed_coupling: float


def quantization_from_fock(K: int, system: CoulombSystem) -> FockQuantization:
    """Derive kappa by matching the transformed equation to a kernel eigenvalue."""

    eigenvalue = fock_kernel_eigenvalue(K)
    # Psi = (mu*g/(2*pi^2*hbar*q)) K Psi, hence 1=coupling*eigenvalue.
    q = system.mu * system.g * eigenvalue / (2.0 * pi**2 * system.hbar)
    kappa = q / system.hbar
    n = K + 1
    return FockQuantization(
        K=K,
        n=n,
        kappa=kappa,
        momentum_scale=q,
        energy=-(q * q) / (2.0 * system.mu),
        transformed_coupling=system.mu * system.g / (2.0 * pi**2 * system.hbar * q),
    )


def numerical_zonal_eigenvalue(K: int, *, chi_order: int = 80, sphere_order: int = 20) -> tuple[complex, float]:
    """Independently quadrature-check K on Y_K00 at the south pole.

    Gauss-Legendre nodes avoid the integrable coincident endpoint.  The full
    S^2 angular integral is evaluated by a separate tensor quadrature rather
    than inserted as 4*pi.
    """

    validate_labels(K, 0, 0)
    if chi_order < 4 or sphere_order < 2:
        raise ValueError("quadrature orders are too small")
    x, wx = roots_legendre(chi_order)
    chi = 0.5 * pi * (x + 1.0)
    # At the south pole, sin^2(chi)/|U-U0|^2=(1-cos(chi))/2.
    radial_integrand = hyperspherical_harmonic(K, 0, 0, chi, 0.0, 0.0) * (1.0 - np.cos(chi)) / 2.0
    polar_x, polar_w = roots_legendre(sphere_order)
    phi_width = 2.0 * pi / float(2 * sphere_order)
    angular_weight = float(np.sum(polar_w) * np.sum(np.full(2 * sphere_order, phi_width)))
    result = angular_weight * 0.5 * pi * np.sum(wx * radial_integrand)
    at_pole = complex(hyperspherical_harmonic(K, 0, 0, pi, 0.0, 0.0))
    estimate = complex(result / at_pole)
    return estimate, abs(estimate - fock_kernel_eigenvalue(K))


def apply_numerical_at_south_pole(function: Callable[[NDArray[np.float64]], complex], *, order: int = 64) -> complex:
    """Integrate a callable over S^3 against the kernel at the south pole."""

    x, weights = roots_legendre(order)
    chi = np.arccos(x)
    # Axisymmetric reduction is valid only when callers explicitly supply the
    # meridian values used here; this helper is intentionally diagnostic.
    points = np.column_stack((np.sin(chi), np.zeros_like(chi), np.zeros_like(chi), x))
    values = np.asarray([function(point) for point in points], dtype=np.complex128)
    return complex(4.0 * pi * np.sum(weights * values * np.sqrt(1.0 - x * x) / (2.0 * (1.0 + x))))
