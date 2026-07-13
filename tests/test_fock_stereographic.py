import numpy as np
import pytest

from hydrogen_s3.fock.stereographic import (
    chord_distance_squared,
    measure_density,
    momentum_to_s3,
    s3_to_momentum,
)


def test_stereographic_roundtrip_norm_and_chord_identity() -> None:
    p = np.asarray([[0.0, 0.0, 0.0], [0.2, -1.0, 3.0], [12.0, 4.0, -2.0]])
    u = momentum_to_s3(p, 0.7)
    assert np.max(np.abs(np.sum(u * u, axis=1) - 1.0)) < 1e-12
    assert np.max(np.abs(s3_to_momentum(u, 0.7) - p)) < 1e-12
    assert chord_distance_squared(p[1], p[2], 0.7) == pytest.approx(np.sum((u[1] - u[2]) ** 2), abs=1e-12)


def test_measure_density_matches_radial_jacobian() -> None:
    p = np.asarray([0.4, -0.3, 0.8])
    h = 1e-6
    jacobian = np.column_stack(
        [
            (momentum_to_s3(p + h * np.eye(3)[i], 1.2) - momentum_to_s3(p - h * np.eye(3)[i], 1.2)) / (2 * h)
            for i in range(3)
        ]
    )
    gram_density = np.sqrt(np.linalg.det(jacobian.T @ jacobian))
    assert gram_density == pytest.approx(measure_density(p, 1.2), rel=2e-9)


def test_pole_and_scale_validation() -> None:
    with pytest.raises(ValueError, match="pole"):
        s3_to_momentum([0.0, 0.0, 0.0, 1.0], 1.0)
    with pytest.raises(ValueError):
        momentum_to_s3([0.0, 0.0, 0.0], 0.0)
