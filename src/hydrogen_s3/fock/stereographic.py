"""Energy-dependent stereographic compactification of momentum space."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def _positive_scale(kappa: float) -> float:
    value = float(kappa)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("kappa must be positive and finite")
    return value


def _last_dimension(value: ArrayLike, size: int, name: str) -> FloatArray:
    array = np.asarray(value, dtype=np.float64)
    if array.ndim == 0 or array.shape[-1] != size:
        raise ValueError(f"{name} must have final dimension {size}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def momentum_to_s3(p: ArrayLike, kappa: float) -> FloatArray:
    """Map (...,3) momentum coordinates to (...,4) unit vectors U."""

    scale = _positive_scale(kappa)
    momentum = _last_dimension(p, 3, "p")
    p2 = np.sum(momentum * momentum, axis=-1)
    denominator = p2 + scale * scale
    spatial = 2.0 * scale * momentum / denominator[..., None]
    u4 = (p2 - scale * scale) / denominator
    return np.concatenate((spatial, u4[..., None]), axis=-1)


def s3_to_momentum(u: ArrayLike, kappa: float, *, pole_tolerance: float = 1e-14) -> FloatArray:
    """Map (...,4) unit vectors U to momentum; reject the north pole."""

    scale = _positive_scale(kappa)
    point = _last_dimension(u, 4, "u")
    norm_error = np.abs(np.sum(point * point, axis=-1) - 1.0)
    if np.any(norm_error > 1e-10):
        raise ValueError("u must lie on the unit S^3")
    denominator = 1.0 - point[..., 3]
    if np.any(denominator <= pole_tolerance):
        raise ValueError("the compactification pole u4=1 represents infinite momentum")
    return scale * point[..., :3] / denominator[..., None]


def hyperspherical_angles(u: ArrayLike) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return chi, theta, phi for (...,4) S^3 points."""

    point = _last_dimension(u, 4, "u")
    if np.any(np.abs(np.sum(point * point, axis=-1) - 1.0) > 1e-10):
        raise ValueError("u must lie on the unit S^3")
    chi = np.arccos(np.clip(point[..., 3], -1.0, 1.0))
    radius = np.linalg.norm(point[..., :3], axis=-1)
    safe = np.where(radius > 0.0, radius, 1.0)
    theta = np.arccos(np.clip(point[..., 2] / safe, -1.0, 1.0))
    phi = np.mod(np.arctan2(point[..., 1], point[..., 0]), 2.0 * np.pi)
    return chi, theta, phi


def measure_density(p: ArrayLike, kappa: float) -> FloatArray:
    """Return dOmega_3/d^3p=(2*kappa/(p^2+kappa^2))^3."""

    scale = _positive_scale(kappa)
    momentum = _last_dimension(p, 3, "p")
    p2 = np.sum(momentum * momentum, axis=-1)
    return np.asarray((2.0 * scale / (p2 + scale * scale)) ** 3, dtype=np.float64)


def chord_distance_squared(p: ArrayLike, p_prime: ArrayLike, kappa: float) -> FloatArray:
    """Return the momentum-coordinate expression for |U-U'|^2."""

    scale = _positive_scale(kappa)
    left = _last_dimension(p, 3, "p")
    right = _last_dimension(p_prime, 3, "p_prime")
    difference2 = np.sum((left - right) ** 2, axis=-1)
    left2 = np.sum(left * left, axis=-1)
    right2 = np.sum(right * right, axis=-1)
    return np.asarray(
        4.0 * scale * scale * difference2 / ((left2 + scale * scale) * (right2 + scale * scale)),
        dtype=np.float64,
    )
