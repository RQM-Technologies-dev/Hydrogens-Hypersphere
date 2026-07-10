"""Angular tensor matrices using an explicit Wigner-Eckart convention.

Convention:

    <ell_f m_f | T^(k)_q | ell_i m_i>
      = <ell_i m_i; k q | ell_f m_f>
        * <ell_f || T^(k) || ell_i> / sqrt(2 ell_f + 1).

The reduced matrix element is an explicit caller-supplied parameter. These
utilities provide angular factors only; they are not physical dipole operators
unless radial matrix elements and physical normalization are supplied.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

from hydrogen_s3.clebsch_gordan import clebsch_gordan_from_doubled

ComplexMatrix = NDArray[np.complex128]


def _validate_int(value: int, name: str, *, minimum: int = 0) -> None:
    if not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def magnetic_labels(ell: int) -> list[int]:
    """Return integer m labels for orbital angular momentum ell."""

    _validate_int(ell, "ell")
    return list(range(-ell, ell + 1))


def selection_rule_allowed(ell_i: int, m_i: int, ell_f: int, m_f: int, rank: int, component: int) -> bool:
    """Return whether angular-momentum tensor selection rules permit a matrix element."""

    _validate_int(ell_i, "ell_i")
    _validate_int(ell_f, "ell_f")
    _validate_int(rank, "rank")
    if not isinstance(component, int) or abs(component) > rank:
        raise ValueError("component must be an integer with |component| <= rank")
    if abs(m_i) > ell_i or abs(m_f) > ell_f:
        return False
    if m_f != m_i + component:
        return False
    return abs(ell_i - rank) <= ell_f <= ell_i + rank


def wigner_eckart_factor(
    ell_i: int,
    m_i: int,
    ell_f: int,
    m_f: int,
    rank: int,
    component: int,
    *,
    reduced_matrix_element: complex = 1.0,
) -> complex:
    """Return the angular tensor matrix element under the documented convention."""

    if not selection_rule_allowed(ell_i, m_i, ell_f, m_f, rank, component):
        return 0.0 + 0.0j
    cg = clebsch_gordan_from_doubled(2 * ell_i, 2 * m_i, 2 * rank, 2 * component, 2 * ell_f, 2 * m_f)
    return complex(reduced_matrix_element) * complex(cg.evalf()) / math.sqrt(2 * ell_f + 1)


def angular_tensor_matrix(
    ell_i: int,
    ell_f: int,
    rank: int,
    component: int,
    *,
    reduced_matrix_element: complex = 1.0,
) -> ComplexMatrix:
    """Return rows m_f, columns m_i for an angular tensor component."""

    labels_i = magnetic_labels(ell_i)
    labels_f = magnetic_labels(ell_f)
    out = np.zeros((len(labels_f), len(labels_i)), dtype=np.complex128)
    for row, m_f in enumerate(labels_f):
        for col, m_i in enumerate(labels_i):
            out[row, col] = wigner_eckart_factor(
                ell_i,
                m_i,
                ell_f,
                m_f,
                rank,
                component,
                reduced_matrix_element=reduced_matrix_element,
            )
    return out


def rank1_selection_allowed(ell_i: int, m_i: int, ell_f: int, m_f: int, component: int) -> bool:
    """Return rank-1 angular selection-rule status for q = component."""

    return selection_rule_allowed(ell_i, m_i, ell_f, m_f, 1, component)
