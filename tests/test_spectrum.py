import math

import pytest

from hydrogen_s3.spectrum import (
    RYDBERG_EV,
    angular_labels,
    calibrated_spectral_energy_from_K,
    rydberg_shell_energy,
    s3_laplacian_eigenvalue,
    shell_dimension,
    shell_number,
    shifted_laplacian_eigenvalue,
    transition_energy_ev,
    vacuum_transition_wavelength_nm,
)


def test_s3_shell_identities() -> None:
    for K in range(12):
        n = K + 1
        assert s3_laplacian_eigenvalue(K) == K * (K + 2)
        assert shifted_laplacian_eigenvalue(K) == n * n
        assert shell_number(K) == n
        assert shell_dimension(K) == n * n
        assert len(angular_labels(K)) == n * n


def test_angular_labels_order_and_counts() -> None:
    assert angular_labels(2) == [
        (0, 0),
        (1, -1),
        (1, 0),
        (1, 1),
        (2, -2),
        (2, -1),
        (2, 0),
        (2, 1),
        (2, 2),
    ]


def test_rydberg_energy_functions_are_consistent() -> None:
    for K in range(8):
        assert calibrated_spectral_energy_from_K(K) == pytest.approx(rydberg_shell_energy(K + 1), abs=1e-15)
    assert rydberg_shell_energy(1) == pytest.approx(-RYDBERG_EV)
    assert rydberg_shell_energy(2) == pytest.approx(-RYDBERG_EV / 4.0)


def test_transition_energy_and_vacuum_wavelength() -> None:
    energy = transition_energy_ev(3, 2)
    assert energy == pytest.approx(RYDBERG_EV * (1 / 4 - 1 / 9))
    assert vacuum_transition_wavelength_nm(3, 2) == pytest.approx(1239.8419843320026 / energy)


@pytest.mark.parametrize(
    ("function", "args"),
    [
        (s3_laplacian_eigenvalue, (-1,)),
        (shifted_laplacian_eigenvalue, (-1,)),
        (shell_number, (-1,)),
        (shell_dimension, (-1,)),
        (angular_labels, (-1,)),
        (rydberg_shell_energy, (0,)),
        (calibrated_spectral_energy_from_K, (-1,)),
        (transition_energy_ev, (2, 2)),
        (vacuum_transition_wavelength_nm, (2, 2)),
    ],
)
def test_input_validation(function: object, args: tuple[int, ...]) -> None:
    with pytest.raises(ValueError):
        function(*args)  # type: ignore[misc,operator]


def test_custom_rydberg_scale() -> None:
    assert math.isclose(rydberg_shell_energy(3, rydberg_ev=9.0), -1.0)
    assert math.isclose(calibrated_spectral_energy_from_K(2, rydberg_ev=9.0), -1.0)
