import math

import pytest

from hydrogen_s3.fock.conventions import CoulombSystem
from hydrogen_s3.fock.coulomb_kernel import (
    fock_kernel_eigenvalue,
    momentum_coulomb_coefficient,
    numerical_zonal_eigenvalue,
)


def test_physical_and_atomic_momentum_kernel_coefficient() -> None:
    assert momentum_coulomb_coefficient(CoulombSystem.atomic_units()) == pytest.approx(1.0 / (2.0 * math.pi**2))


@pytest.mark.parametrize("K", range(4))
def test_independent_numerical_kernel_eigenvalue(K: int) -> None:
    estimate, error = numerical_zonal_eigenvalue(K, chi_order=80, sphere_order=12)
    assert estimate.real == pytest.approx(fock_kernel_eigenvalue(K), abs=2e-11)
    assert abs(estimate.imag) < 1e-14
    assert error < 2e-11
