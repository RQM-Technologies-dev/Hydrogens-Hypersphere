from fractions import Fraction

import pytest
from sympy import Rational, simplify, sqrt
from sympy.physics.wigner import clebsch_gordan as sympy_clebsch_gordan

from hydrogen_s3.clebsch_gordan import (
    clebsch_gordan_exact,
    clebsch_gordan_from_doubled,
    clebsch_gordan_value,
    m_values_doubled,
    quantum_number_to_doubled,
)


def test_known_highest_weight_coefficients_are_positive() -> None:
    assert clebsch_gordan_exact(Fraction(1, 2), Fraction(1, 2), Fraction(1, 2), Fraction(1, 2), 1, 1) == 1
    assert clebsch_gordan_exact(1, 1, 1, 1, 2, 2) == 1
    assert clebsch_gordan_value(1, 1, 1, 1, 2, 2) == pytest.approx(1.0 + 0.0j)


def test_known_standard_values() -> None:
    assert simplify(clebsch_gordan_exact(1, 1, 1, -1, 0, 0) - sqrt(3) / 3) == 0
    assert simplify(clebsch_gordan_exact(1, 0, 1, 0, 0, 0) + sqrt(3) / 3) == 0
    assert simplify(clebsch_gordan_exact(1, 1, 1, 0, 2, 1) - sqrt(2) / 2) == 0


def test_exhaustive_agreement_with_direct_sympy_through_j4() -> None:
    for j1_2 in range(0, 9):
        for j2_2 in range(0, 9):
            for m1_2 in m_values_doubled(j1_2):
                for m2_2 in m_values_doubled(j2_2):
                    min_j_2 = abs(j1_2 - j2_2)
                    max_j_2 = j1_2 + j2_2
                    for j_2 in range(min_j_2, max_j_2 + 1, 2):
                        for m_2 in m_values_doubled(j_2):
                            ours = clebsch_gordan_from_doubled(j1_2, m1_2, j2_2, m2_2, j_2, m_2)
                            direct = sympy_clebsch_gordan(
                                Rational(j1_2, 2),
                                Rational(j2_2, 2),
                                Rational(j_2, 2),
                                Rational(m1_2, 2),
                                Rational(m2_2, 2),
                                Rational(m_2, 2),
                            )
                            assert simplify(ours - direct) == 0


def test_invalid_coefficients_return_zero() -> None:
    assert clebsch_gordan_from_doubled(1, 1, 1, 1, 0, 0) == 0
    assert clebsch_gordan_from_doubled(1, 3, 1, 1, 2, 4) == 0


def test_quantum_number_validation() -> None:
    assert quantum_number_to_doubled(Fraction(3, 2), "j") == 3
    assert quantum_number_to_doubled(1.5, "j") == 3
    with pytest.raises(ValueError):
        quantum_number_to_doubled(0.25, "j")
