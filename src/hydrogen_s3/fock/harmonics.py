"""Normalized complex scalar hyperspherical harmonics on the unit S^3."""

from __future__ import annotations

from math import exp, lgamma, pi

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import eval_gegenbauer, sph_harm_y

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


def validate_labels(K: int, ell: int, m: int) -> None:
    if not isinstance(K, int) or isinstance(K, bool) or K < 0:
        raise ValueError("K must be a non-negative integer")
    if not isinstance(ell, int) or isinstance(ell, bool) or not 0 <= ell <= K:
        raise ValueError("ell must be an integer satisfying 0 <= ell <= K")
    if not isinstance(m, int) or isinstance(m, bool) or abs(m) > ell:
        raise ValueError("m must be an integer satisfying -ell <= m <= ell")


def normalization(K: int, ell: int) -> float:
    validate_labels(K, ell, 0)
    log_value = (
        ell * np.log(2.0)
        + lgamma(ell + 1.0)
        + 0.5 * (np.log(2.0 * (K + 1.0) / pi) + lgamma(K - ell + 1.0) - lgamma(K + ell + 2.0))
    )
    return float(exp(log_value))


def hyperspherical_radial(K: int, ell: int, chi: ArrayLike) -> FloatArray:
    """Evaluate Pi_Kell(chi), preserving the input array shape."""

    validate_labels(K, ell, 0)
    angle = np.asarray(chi, dtype=np.float64)
    if np.any(~np.isfinite(angle)) or np.any((angle < 0.0) | (angle > pi)):
        raise ValueError("chi must be finite and lie in [0, pi]")
    values = normalization(K, ell) * np.sin(angle) ** ell * eval_gegenbauer(K - ell, ell + 1, np.cos(angle))
    return np.asarray(values, dtype=np.float64)


def hyperspherical_harmonic(K: int, ell: int, m: int, chi: ArrayLike, theta: ArrayLike, phi: ArrayLike) -> ComplexArray:
    """Evaluate Y_Kellm=Pi_Kell Y_ellm with SciPy's Condon-Shortley phase."""

    validate_labels(K, ell, m)
    radial = hyperspherical_radial(K, ell, chi)
    polar = np.asarray(theta, dtype=np.float64)
    azimuth = np.asarray(phi, dtype=np.float64)
    if np.any(~np.isfinite(polar)) or np.any((polar < 0.0) | (polar > pi)):
        raise ValueError("theta must be finite and lie in [0, pi]")
    if np.any(~np.isfinite(azimuth)):
        raise ValueError("phi must be finite")
    return np.asarray(radial * sph_harm_y(ell, m, polar, azimuth), dtype=np.complex128)


def shell_labels(K: int) -> list[tuple[int, int, int]]:
    if not isinstance(K, int) or isinstance(K, bool) or K < 0:
        raise ValueError("K must be a non-negative integer")
    return [(K, ell, m) for ell in range(K + 1) for m in range(-ell, ell + 1)]


def shell_dimension(K: int) -> int:
    if not isinstance(K, int) or isinstance(K, bool) or K < 0:
        raise ValueError("K must be a non-negative integer")
    return (K + 1) ** 2


def laplacian_eigenvalue(K: int) -> int:
    if not isinstance(K, int) or isinstance(K, bool) or K < 0:
        raise ValueError("K must be a non-negative integer")
    return K * (K + 2)
