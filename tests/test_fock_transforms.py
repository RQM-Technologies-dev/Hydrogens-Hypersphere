import numpy as np
import pytest

from hydrogen_s3.fock.bound_states import analytic_position_radial
from hydrogen_s3.fock.transforms import inverse_hankel_radial


@pytest.mark.parametrize(("n", "ell"), [(1, 0), (2, 0), (2, 1), (3, 2)])
def test_inverse_transform_matches_laguerre_solution(n: int, ell: int) -> None:
    r = np.asarray([0.2, 1.0, 3.0, 7.0])
    reconstructed = inverse_hankel_radial(n, ell, r, tolerance=2e-7)
    expected = analytic_position_radial(n, ell, r)
    assert np.max(np.abs(reconstructed - expected)) < 8e-7
    assert np.max(np.abs(reconstructed.imag)) < 1e-12


def test_tighter_quadrature_improves_cutoff_result() -> None:
    expected = float(analytic_position_radial(1, 0, 0.2))
    coarse = complex(inverse_hankel_radial(1, 0, 0.2, tolerance=2e-5))
    fine = complex(inverse_hankel_radial(1, 0, 0.2, tolerance=2e-7))
    assert abs(fine - expected) < abs(coarse - expected)
