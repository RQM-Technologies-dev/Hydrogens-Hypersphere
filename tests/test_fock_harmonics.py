import numpy as np
from scipy.special import roots_legendre

from hydrogen_s3.fock.harmonics import (
    hyperspherical_harmonic,
    laplacian_eigenvalue,
    shell_dimension,
    shell_labels,
)


def test_shell_enumeration_and_eigenvalue() -> None:
    for K in range(6):
        assert len(shell_labels(K)) == shell_dimension(K) == (K + 1) ** 2
        assert laplacian_eigenvalue(K) == K * (K + 2)


def test_low_harmonic_orthonormality() -> None:
    x, wx = roots_legendre(48)
    chi = 0.5 * np.pi * (x + 1.0)
    wx = 0.5 * np.pi * wx
    z, wz = roots_legendre(24)
    theta = np.arccos(z)
    phi = np.linspace(0.0, 2.0 * np.pi, 49)[:-1]
    labels = shell_labels(2)
    values = []
    weights = []
    for i, angle in enumerate(chi):
        for j, polar in enumerate(theta):
            for azimuth in phi:
                values.append([hyperspherical_harmonic(*label, angle, polar, azimuth) for label in labels])
                weights.append(wx[i] * np.sin(angle) ** 2 * wz[j] * 2.0 * np.pi / len(phi))
    matrix = np.asarray(values)
    gram = matrix.conj().T @ (np.asarray(weights)[:, None] * matrix)
    assert np.linalg.norm(gram - np.eye(9), ord=np.inf) < 2e-12


def test_radial_laplacian_spot_check() -> None:
    K, ell, m = 3, 1, 0
    chi, theta, phi = 1.1, 0.9, 0.4
    h = 2e-5

    def value(c: float, t: float, f: float) -> complex:
        return complex(hyperspherical_harmonic(K, ell, m, c, t, f))

    f0 = value(chi, theta, phi)
    dc = (value(chi + h, theta, phi) - value(chi - h, theta, phi)) / (2 * h)
    d2c = (value(chi + h, theta, phi) - 2 * f0 + value(chi - h, theta, phi)) / h**2
    dt = (value(chi, theta + h, phi) - value(chi, theta - h, phi)) / (2 * h)
    d2t = (value(chi, theta + h, phi) - 2 * f0 + value(chi, theta - h, phi)) / h**2
    d2p = (value(chi, theta, phi + h) - 2 * f0 + value(chi, theta, phi - h)) / h**2
    laplacian = d2c + 2 / np.tan(chi) * dc + (d2t + dt / np.tan(theta) + d2p / np.sin(theta) ** 2) / np.sin(chi) ** 2
    assert abs(-laplacian - K * (K + 2) * f0) < 2e-6
