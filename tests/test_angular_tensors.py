import numpy as np
import pytest

from hydrogen_s3.angular_tensors import (
    angular_tensor_matrix,
    rank1_selection_allowed,
    selection_rule_allowed,
    wigner_eckart_factor,
)


def test_selection_rules() -> None:
    assert selection_rule_allowed(1, 0, 2, 1, 1, 1)
    assert not selection_rule_allowed(1, 0, 2, 0, 1, 1)
    assert not selection_rule_allowed(0, 0, 2, 0, 1, 0)
    assert rank1_selection_allowed(1, -1, 0, 0, 1)


def test_wigner_eckart_reduced_matrix_element_scaling() -> None:
    base = wigner_eckart_factor(1, 0, 2, 1, 1, 1, reduced_matrix_element=1.0)
    scaled = wigner_eckart_factor(1, 0, 2, 1, 1, 1, reduced_matrix_element=3.5 - 2.0j)
    assert scaled == pytest.approx((3.5 - 2.0j) * base)


def test_angular_tensor_matrix_shape_and_entries() -> None:
    matrix = angular_tensor_matrix(1, 2, 1, 1, reduced_matrix_element=2.0)
    assert matrix.shape == (5, 3)
    for row, m_f in enumerate(range(-2, 3)):
        for col, m_i in enumerate(range(-1, 2)):
            expected = wigner_eckart_factor(1, m_i, 2, m_f, 1, 1, reduced_matrix_element=2.0)
            assert matrix[row, col] == pytest.approx(expected)
    assert np.count_nonzero(np.abs(matrix) > 1e-14) > 0


def test_invalid_rank_component() -> None:
    with pytest.raises(ValueError):
        angular_tensor_matrix(1, 1, 1, 2)
