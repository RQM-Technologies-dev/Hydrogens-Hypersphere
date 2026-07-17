import numpy as np
import pytest

from hydrogen_s3.angular_momentum import casimir, dimension, jminus, jplus, jx, jy, jz, magnetic_labels


@pytest.mark.parametrize("j", [0, 0.5, 1, 1.5, 2, 4])
def test_dimensions_and_basis_ordering(j: float) -> None:
    labels = magnetic_labels(j)
    assert len(labels) == dimension(j)
    assert labels[0] == -j
    assert labels[-1] == j
    assert labels == tuple(sorted(labels))


@pytest.mark.parametrize("j", [0.5, 1, 1.5, 3])
def test_hermiticity(j: float) -> None:
    for matrix in [jx(j), jy(j), jz(j), casimir(j)]:
        assert np.linalg.norm(matrix - matrix.conj().T, ord=np.inf) < 1e-12


@pytest.mark.parametrize("j", [0.5, 1, 1.5, 3])
def test_commutation_relations(j: float) -> None:
    x = jx(j)
    y = jy(j)
    z = jz(j)
    assert np.linalg.norm(x @ y - y @ x - 1.0j * z, ord=np.inf) < 1e-12
    assert np.linalg.norm(y @ z - z @ y - 1.0j * x, ord=np.inf) < 1e-12
    assert np.linalg.norm(z @ x - x @ z - 1.0j * y, ord=np.inf) < 1e-12


@pytest.mark.parametrize("j", [0, 0.5, 1, 2.5, 4])
def test_casimir(j: float) -> None:
    expected = j * (j + 1.0) * np.eye(dimension(j), dtype=np.complex128)
    assert np.linalg.norm(casimir(j) - expected, ord=np.inf) < 1e-12


def test_ladder_matrix_elements() -> None:
    j = 2
    labels = magnetic_labels(j)
    plus = jplus(j)
    minus = jminus(j)
    for col, m in enumerate(labels):
        if m < j:
            row = labels.index(m + 1)
            expected = np.sqrt((j - m) * (j + m + 1))
            assert plus[row, col] == pytest.approx(expected)
        if m > -j:
            row = labels.index(m - 1)
            expected = np.sqrt((j + m) * (j - m + 1))
            assert minus[row, col] == pytest.approx(expected)
