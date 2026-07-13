"""Compact deterministic diagnostics used by reports and CI tests."""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from hydrogen_s3.fock.bound_states import momentum_radial
from hydrogen_s3.fock.conventions import CoulombSystem
from hydrogen_s3.fock.coulomb_kernel import numerical_zonal_eigenvalue, quantization_from_fock
from hydrogen_s3.fock.harmonics import hyperspherical_harmonic
from hydrogen_s3.fock.so4 import so4_diagnostics
from hydrogen_s3.fock.stereographic import chord_distance_squared, measure_density, momentum_to_s3, s3_to_momentum
from hydrogen_s3.fock.transforms import compare_radial_reconstruction


def fock_diagnostics() -> dict[str, object]:
    system = CoulombSystem.atomic_units()
    points = np.asarray([[0.0, 0.0, 0.0], [0.2, -0.4, 1.1], [3.0, 2.0, -1.0]])
    mapped = momentum_to_s3(points, 0.7)
    roundtrip = s3_to_momentum(mapped, 0.7)
    chord = chord_distance_squared(points[1], points[2], 0.7)
    direct_chord = float(np.sum((mapped[1] - mapped[2]) ** 2))
    h = 1e-6
    jacobian = np.column_stack(
        [
            (momentum_to_s3(points[1] + h * np.eye(3)[i], 0.7) - momentum_to_s3(points[1] - h * np.eye(3)[i], 0.7))
            / (2.0 * h)
            for i in range(3)
        ]
    )
    numerical_density = float(np.sqrt(np.linalg.det(jacobian.T @ jacobian)))
    momentum_normalizations: dict[str, float] = {}
    for n, ell in ((1, 0), (2, 0), (2, 1), (3, 2)):
        integral, _ = quad(
            lambda p, n=n, ell=ell: float(abs(momentum_radial(n, ell, p, system)) ** 2 * p * p),
            0.0,
            np.inf,
        )
        momentum_normalizations[f"{n},{ell}"] = abs(integral - 1.0)
    kernel = {str(K): numerical_zonal_eigenvalue(K, chi_order=120)[1] for K in range(3)}
    chi_x, chi_w = roots_legendre(48)
    chi = 0.5 * np.pi * (chi_x + 1.0)
    harmonic_norm = float(
        4.0
        * np.pi
        * 0.5
        * np.pi
        * np.sum(chi_w * np.sin(chi) ** 2 * np.abs(hyperspherical_harmonic(2, 0, 0, chi, 0.0, 0.0)) ** 2)
    )
    quantization = {
        str(K): {
            "n": result.n,
            "kappa": result.kappa,
            "energy": result.energy,
        }
        for K in range(3)
        for result in (quantization_from_fock(K, system),)
    }
    radii = np.asarray([0.2, 1.0, 3.0, 7.0])
    radial = {
        f"{n},{ell}": {
            "max_absolute_error": comparison.max_absolute_error,
            "relative_l2_error": comparison.relative_l2_error,
            "sampled_normalization_error": comparison.normalization_error,
            "phase_alignment_real": comparison.phase_alignment.real,
            "phase_alignment_imag": comparison.phase_alignment.imag,
        }
        for n, ell in ((1, 0), (2, 0), (2, 1), (3, 2))
        for comparison in (compare_radial_reconstruction(n, ell, radii, system, tolerance=2e-7),)
    }
    return {
        "system": {"mu": system.mu, "g": system.g, "hbar": system.hbar, "a0": system.bohr_radius},
        "stereographic": {
            "roundtrip_max_error": float(np.max(np.abs(roundtrip - points))),
            "unit_s3_max_error": float(np.max(np.abs(np.sum(mapped * mapped, axis=-1) - 1.0))),
            "chord_error": abs(float(chord) - direct_chord),
            "jacobian_error": abs(float(measure_density(points[1], 0.7)) - numerical_density),
        },
        "harmonic_normalization_error": abs(harmonic_norm - 1.0),
        "kernel_eigenvalue_errors": kernel,
        "quantization": quantization,
        "momentum_normalization_errors": momentum_normalizations,
        "radial_reconstruction": radial,
        "so4": so4_diagnostics(2),
    }
