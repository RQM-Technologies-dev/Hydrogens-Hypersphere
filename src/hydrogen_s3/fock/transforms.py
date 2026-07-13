"""Numerical inverse spherical-Bessel transform to position space."""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import quad
from scipy.special import spherical_jn

from hydrogen_s3.fock.bound_states import analytic_position_radial, momentum_radial, quantum_numbers
from hydrogen_s3.fock.conventions import CoulombSystem

ComplexArray = NDArray[np.complex128]


class TransformConvergenceError(RuntimeError):
    """Raised when mapped semi-infinite quadrature misses the requested tolerance."""


def _integral_component(n: int, ell: int, r: float, system: CoulombSystem, tolerance: float) -> complex:
    scale = system.momentum_scale_n(n)

    def integrand(p: float) -> complex:
        return complex(p * p * momentum_radial(n, ell, p, system) * spherical_jn(ell, p * r / system.hbar))

    def integrate_to(cutoff: float) -> tuple[complex, float]:
        real, real_error = quad(
            lambda value: integrand(value).real, 0.0, cutoff, epsabs=tolerance / 4.0, epsrel=tolerance / 4.0, limit=500
        )
        imag, imag_error = quad(
            lambda value: integrand(value).imag, 0.0, cutoff, epsabs=tolerance / 4.0, epsrel=tolerance / 4.0, limit=500
        )
        return complex(real, imag), max(real_error, imag_error)

    if r == 0.0:
        real, real_error = quad(lambda value: integrand(value).real, 0.0, np.inf, epsabs=tolerance, epsrel=tolerance)
        imag, imag_error = quad(lambda value: integrand(value).imag, 0.0, np.inf, epsabs=tolerance, epsrel=tolerance)
        integral = complex(real, imag)
        estimated_error = max(real_error, imag_error)
    else:
        cutoff = 64.0 * scale + 16.0 * system.hbar / r
        previous, estimated_error = integrate_to(cutoff)
        integral = previous
        for _ in range(5):
            cutoff *= 2.0
            integral, quadrature_error = integrate_to(cutoff)
            estimated_error = max(quadrature_error, abs(integral - previous))
            if estimated_error <= 5.0 * tolerance:
                break
            previous = integral
        else:
            raise TransformConvergenceError(
                f"inverse transform cutoff-doubling error {estimated_error:.3e} exceeds requested tolerance"
            )
    coefficient = 4.0 * pi * (1.0j) ** ell / (2.0 * pi * system.hbar) ** 1.5
    if estimated_error > 20.0 * tolerance:
        raise TransformConvergenceError(
            f"inverse transform estimated error {estimated_error:.3e} exceeds requested tolerance"
        )
    return complex(coefficient * integral)


def inverse_hankel_radial(
    n: int,
    ell: int,
    r: ArrayLike,
    system: CoulombSystem | None = None,
    *,
    tolerance: float = 1e-9,
) -> ComplexArray:
    """Numerically reconstruct R_nell(r) from the Fock momentum radial state."""

    quantum_numbers(n, ell)
    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")
    active = CoulombSystem.atomic_units() if system is None else system
    radius = np.asarray(r, dtype=np.float64)
    if np.any(~np.isfinite(radius)) or np.any(radius < 0.0):
        raise ValueError("r must be non-negative and finite")
    flat = [_integral_component(n, ell, float(value), active, tolerance) for value in radius.ravel()]
    return np.asarray(flat, dtype=np.complex128).reshape(radius.shape)


@dataclass(frozen=True, slots=True)
class RadialComparison:
    max_absolute_error: float
    relative_l2_error: float
    normalization_error: float
    phase_alignment: complex
    reconstructed: ComplexArray
    analytic: NDArray[np.float64]


def compare_radial_reconstruction(
    n: int,
    ell: int,
    r: ArrayLike,
    system: CoulombSystem | None = None,
    *,
    tolerance: float = 1e-9,
) -> RadialComparison:
    """Compare numerical and associated-Laguerre radial functions on a grid."""

    active = CoulombSystem.atomic_units() if system is None else system
    radius = np.asarray(r, dtype=np.float64)
    if radius.ndim != 1 or radius.size < 2 or np.any(np.diff(radius) <= 0.0):
        raise ValueError("r must be a strictly increasing one-dimensional grid")
    reconstructed = inverse_hankel_radial(n, ell, radius, active, tolerance=tolerance)
    analytic = analytic_position_radial(n, ell, radius, active)
    overlap = np.trapz(np.conj(reconstructed) * analytic * radius**2, radius)
    phase = complex(overlap / abs(overlap)) if abs(overlap) > 0.0 else 1.0 + 0.0j
    aligned = reconstructed * phase
    difference = aligned - analytic
    reference_norm = float(np.trapz(np.abs(analytic) ** 2 * radius**2, radius))
    reconstructed_norm = float(np.trapz(np.abs(reconstructed) ** 2 * radius**2, radius))
    relative = sqrt(float(np.trapz(np.abs(difference) ** 2 * radius**2, radius)) / reference_norm)
    return RadialComparison(
        max_absolute_error=float(np.max(np.abs(difference))),
        relative_l2_error=relative,
        normalization_error=abs(reconstructed_norm - 1.0),
        phase_alignment=phase,
        reconstructed=reconstructed,
        analytic=analytic,
    )
