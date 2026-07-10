"""Finite-dimensional angular-momentum matrices."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from hydrogen_s3.clebsch_gordan import m_values_doubled, quantum_number_to_doubled

ComplexMatrix = NDArray[np.complex128]


def _j2(j: int | float) -> int:
    j_2 = quantum_number_to_doubled(j, "j")
    if j_2 < 0:
        raise ValueError("j must be non-negative")
    return j_2


def magnetic_labels(j: int | float) -> tuple[float, ...]:
    """Return m labels in ascending order: -j, -j+1, ..., j."""

    return tuple(value / 2.0 for value in m_values_doubled(_j2(j)))


def dimension(j: int | float) -> int:
    """Return dim V_j = 2j+1."""

    return _j2(j) + 1


def jz(j: int | float) -> ComplexMatrix:
    """Return the Jz matrix in ascending m-basis order."""

    labels = magnetic_labels(j)
    return np.diag(np.asarray(labels, dtype=np.complex128))


def jplus(j: int | float) -> ComplexMatrix:
    """Return the raising operator J+ with standard ladder coefficients."""

    j_2 = _j2(j)
    labels_2 = m_values_doubled(j_2)
    index_by_m = {m_2: index for index, m_2 in enumerate(labels_2)}
    out = np.zeros((len(labels_2), len(labels_2)), dtype=np.complex128)
    j_value = j_2 / 2.0
    for col, m_2 in enumerate(labels_2):
        next_m_2 = m_2 + 2
        if next_m_2 not in index_by_m:
            continue
        m_value = m_2 / 2.0
        out[index_by_m[next_m_2], col] = np.sqrt((j_value - m_value) * (j_value + m_value + 1.0))
    return out


def jminus(j: int | float) -> ComplexMatrix:
    """Return the lowering operator J- as the conjugate transpose of J+."""

    return jplus(j).conj().T


def jx(j: int | float) -> ComplexMatrix:
    """Return Jx = (J+ + J-) / 2."""

    return ((jplus(j) + jminus(j)) / 2.0).astype(np.complex128)


def jy(j: int | float) -> ComplexMatrix:
    """Return Jy = (J+ - J-) / (2i)."""

    return ((jplus(j) - jminus(j)) / (2.0j)).astype(np.complex128)


def casimir(j: int | float) -> ComplexMatrix:
    """Return J^2 = Jx^2 + Jy^2 + Jz^2."""

    x = jx(j)
    y = jy(j)
    z = jz(j)
    return (x @ x + y @ y + z @ z).astype(np.complex128)
