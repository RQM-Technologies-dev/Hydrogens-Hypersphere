"""Cached SymPy Clebsch-Gordan coefficients.

The wrapper uses SymPy's Condon-Shortley/Wigner convention through
``sympy.physics.wigner.clebsch_gordan``. Public arguments may be integers,
half-integer floats, or ``fractions.Fraction`` values; they are converted to
exact SymPy rationals before evaluation.
"""

from __future__ import annotations

from fractions import Fraction
from functools import cache
from numbers import Integral, Real

from sympy import Expr, Rational
from sympy.physics.wigner import clebsch_gordan as sympy_clebsch_gordan


def quantum_number_to_doubled(value: int | float | Fraction, name: str) -> int:
    """Convert an integer or half-integer quantum number to doubled form."""

    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer or half-integer")
    if isinstance(value, Integral):
        return int(value) * 2
    if isinstance(value, Fraction):
        doubled = value * 2
        if doubled.denominator != 1:
            raise ValueError(f"{name} must be an integer or half-integer")
        return int(doubled)
    if isinstance(value, Real):
        doubled_float = float(value) * 2.0
        doubled_int = round(doubled_float)
        if abs(doubled_float - doubled_int) > 1e-12:
            raise ValueError(f"{name} must be an integer or half-integer")
        return int(doubled_int)
    raise ValueError(f"{name} must be an integer or half-integer")


def rational_from_doubled(value_2: int) -> Rational:
    """Return an exact SymPy Rational from a doubled quantum number."""

    if not isinstance(value_2, int):
        raise ValueError("doubled quantum numbers must be integers")
    return Rational(value_2, 2)


def validate_jm_doubled(j_2: int, m_2: int) -> bool:
    """Return whether doubled (j, m) labels define an allowed state."""

    return isinstance(j_2, int) and isinstance(m_2, int) and j_2 >= 0 and abs(m_2) <= j_2 and (j_2 - m_2) % 2 == 0


def m_values_doubled(j_2: int) -> tuple[int, ...]:
    """Return doubled m labels in ascending order for doubled j."""

    if not isinstance(j_2, int) or j_2 < 0:
        raise ValueError("j_2 must be a non-negative integer")
    return tuple(range(-j_2, j_2 + 1, 2))


@cache
def clebsch_gordan_from_doubled(
    j1_2: int,
    m1_2: int,
    j2_2: int,
    m2_2: int,
    j_2: int,
    m_2: int,
) -> Expr:
    """Return exact <j1,m1; j2,m2 | j,m> using doubled quantum numbers."""

    labels = (j1_2, m1_2, j2_2, m2_2, j_2, m_2)
    if not all(isinstance(value, int) for value in labels):
        raise ValueError("all doubled quantum numbers must be integers")
    if not validate_jm_doubled(j1_2, m1_2):
        return Rational(0)
    if not validate_jm_doubled(j2_2, m2_2):
        return Rational(0)
    if not validate_jm_doubled(j_2, m_2):
        return Rational(0)
    if m_2 != m1_2 + m2_2:
        return Rational(0)
    if j_2 < abs(j1_2 - j2_2) or j_2 > j1_2 + j2_2:
        return Rational(0)

    return sympy_clebsch_gordan(
        rational_from_doubled(j1_2),
        rational_from_doubled(j2_2),
        rational_from_doubled(j_2),
        rational_from_doubled(m1_2),
        rational_from_doubled(m2_2),
        rational_from_doubled(m_2),
    )


def clebsch_gordan_exact(
    j1: int | float | Fraction,
    m1: int | float | Fraction,
    j2: int | float | Fraction,
    m2: int | float | Fraction,
    j: int | float | Fraction,
    m: int | float | Fraction,
) -> Expr:
    """Return exact <j1,m1; j2,m2 | j,m> in SymPy's convention."""

    return clebsch_gordan_from_doubled(
        quantum_number_to_doubled(j1, "j1"),
        quantum_number_to_doubled(m1, "m1"),
        quantum_number_to_doubled(j2, "j2"),
        quantum_number_to_doubled(m2, "m2"),
        quantum_number_to_doubled(j, "j"),
        quantum_number_to_doubled(m, "m"),
    )


def clebsch_gordan_value(
    j1: int | float | Fraction,
    m1: int | float | Fraction,
    j2: int | float | Fraction,
    m2: int | float | Fraction,
    j: int | float | Fraction,
    m: int | float | Fraction,
) -> complex:
    """Return a numerical complex Clebsch-Gordan coefficient."""

    return complex(clebsch_gordan_exact(j1, m1, j2, m2, j, m).evalf())
