import numpy as np
import pytest
from scipy.integrate import quad

from hydrogen_s3.fock.bound_states import analytic_position_radial, momentum_radial, quantum_numbers


@pytest.mark.parametrize(("n", "ell"), [(1, 0), (2, 0), (2, 1), (3, 0), (3, 2)])
def test_momentum_and_position_radial_normalization(n: int, ell: int) -> None:
    p_norm, _ = quad(lambda p: abs(momentum_radial(n, ell, p)) ** 2 * p * p, 0.0, np.inf, epsabs=2e-11)
    r_norm, _ = quad(lambda r: analytic_position_radial(n, ell, r) ** 2 * r * r, 0.0, np.inf, epsabs=2e-11)
    assert p_norm == pytest.approx(1.0, abs=2e-11)
    assert r_norm == pytest.approx(1.0, abs=2e-11)


def test_known_1s_momentum_formula_and_phase() -> None:
    p = np.asarray([0.0, 0.3, 2.0])
    expected = 4.0 * np.sqrt(2.0 / np.pi) / (1.0 + p * p) ** 2
    assert momentum_radial(1, 0, p) == pytest.approx(expected)
    assert np.all(np.imag(momentum_radial(2, 1, p[1:])) != 0.0)


def test_quantum_number_validation() -> None:
    with pytest.raises(ValueError):
        quantum_numbers(2, 2, 0)
    with pytest.raises(ValueError):
        quantum_numbers(2, 1, 2)
