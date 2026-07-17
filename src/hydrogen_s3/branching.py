"""Shell-wise S^3/SO(4) branching transform.

For shell K, this module implements the unitary representation transform

    U_K: H_K(S^3) ~= V_j tensor V_j* -> direct sum_{ell=0}^K V_ell,
    j = K/2.

It is a branching transform, not a physical spatial projection or measurement
map. The source basis is ordered as (a, b), with a and b ascending magnetic
labels for V_j and V_j*. The target basis is ordered by increasing ell and then
ascending m. The dual basis is identified with V_j by
e^b -> (-1)^(j-b) |j,-b>, and Clebsch-Gordan coefficients use SymPy's
Condon-Shortley/Wigner convention.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from hydrogen_s3.angular_momentum import casimir, jx, jy, jz
from hydrogen_s3.clebsch_gordan import clebsch_gordan_from_doubled, m_values_doubled
from hydrogen_s3.spectrum import angular_labels, shell_dimension

ComplexMatrix = NDArray[np.complex128]


def _validate_K(K: int) -> None:
    if not isinstance(K, int) or K < 0:
        raise ValueError("K must be a non-negative integer")


def shell_j(K: int) -> float:
    """Return j = K/2 for the representation V_j tensor V_j*."""

    _validate_K(K)
    return K / 2.0


def source_basis(K: int) -> list[tuple[float, float]]:
    """Return source basis labels (a, b) in Kronecker-product order."""

    _validate_K(K)
    values = [m_2 / 2.0 for m_2 in m_values_doubled(K)]
    return [(a, b) for a in values for b in values]


def target_basis(K: int) -> list[tuple[int, int]]:
    """Return target labels (ell, m) ordered by ell and then m."""

    return angular_labels(K)


def branching_transform(K: int) -> ComplexMatrix:
    """Return U_K mapping source coefficients to target coefficients."""

    _validate_K(K)
    rows = target_basis(K)
    cols_2 = [(a_2, b_2) for a_2 in m_values_doubled(K) for b_2 in m_values_doubled(K)]
    out = np.zeros((len(rows), len(cols_2)), dtype=np.complex128)
    for row, (ell, m) in enumerate(rows):
        for col, (a_2, b_2) in enumerate(cols_2):
            phase_exponent = (K - b_2) // 2
            phase = -1.0 if phase_exponent % 2 else 1.0
            value = clebsch_gordan_from_doubled(K, a_2, K, -b_2, 2 * ell, 2 * m)
            out[row, col] = phase * complex(value.evalf())
    return out


def _dual_generator(matrix: ComplexMatrix) -> ComplexMatrix:
    return (-matrix.T).astype(np.complex128)


def source_generators(K: int) -> dict[str, ComplexMatrix]:
    """Construct source-side diagonal generators on V_j tensor V_j*."""

    _validate_K(K)
    j = shell_j(K)
    dim = K + 1
    identity = np.eye(dim, dtype=np.complex128)
    generators = {"x": jx(j), "y": jy(j), "z": jz(j)}
    return {
        axis: (np.kron(matrix, identity) + np.kron(identity, _dual_generator(matrix))).astype(np.complex128)
        for axis, matrix in generators.items()
    }


def _block_diagonal(blocks: list[ComplexMatrix]) -> ComplexMatrix:
    size = sum(block.shape[0] for block in blocks)
    out = np.zeros((size, size), dtype=np.complex128)
    offset = 0
    for block in blocks:
        next_offset = offset + block.shape[0]
        out[offset:next_offset, offset:next_offset] = block
        offset = next_offset
    return out


def target_generators(K: int) -> dict[str, ComplexMatrix]:
    """Construct target-side block generators on direct sum ell=0..K V_ell."""

    _validate_K(K)
    return {
        "x": _block_diagonal([jx(ell) for ell in range(K + 1)]),
        "y": _block_diagonal([jy(ell) for ell in range(K + 1)]),
        "z": _block_diagonal([jz(ell) for ell in range(K + 1)]),
    }


def source_casimir(K: int) -> ComplexMatrix:
    """Return source J^2 from independently constructed source generators."""

    generators = source_generators(K)
    return np.asarray(
        generators["x"] @ generators["x"] + generators["y"] @ generators["y"] + generators["z"] @ generators["z"],
        dtype=np.complex128,
    )


def target_casimir(K: int) -> ComplexMatrix:
    """Return target J^2 from block-diagonal spin-ell generators."""

    return _block_diagonal([casimir(ell) for ell in range(K + 1)])


def branching_diagnostics(K: int) -> dict[str, float | int]:
    """Return unitary, generator-intertwining, commutator, and Casimir errors."""

    _validate_K(K)
    U = branching_transform(K)
    dim = shell_dimension(K)
    identity = np.eye(dim, dtype=np.complex128)
    source = source_generators(K)
    target = target_generators(K)
    source_c = source_casimir(K)
    target_c = target_casimir(K)
    return {
        "K": K,
        "dimension": dim,
        "unitarity_U_Udagger_error": float(np.linalg.norm(U @ U.conj().T - identity, ord=np.inf)),
        "unitarity_Udagger_U_error": float(np.linalg.norm(U.conj().T @ U - identity, ord=np.inf)),
        "intertwining_Jx_error": float(np.linalg.norm(U @ source["x"] @ U.conj().T - target["x"], ord=np.inf)),
        "intertwining_Jy_error": float(np.linalg.norm(U @ source["y"] @ U.conj().T - target["y"], ord=np.inf)),
        "intertwining_Jz_error": float(np.linalg.norm(U @ source["z"] @ U.conj().T - target["z"], ord=np.inf)),
        "source_commutator_xy_error": float(
            np.linalg.norm(source["x"] @ source["y"] - source["y"] @ source["x"] - 1.0j * source["z"], ord=np.inf)
        ),
        "target_commutator_xy_error": float(
            np.linalg.norm(target["x"] @ target["y"] - target["y"] @ target["x"] - 1.0j * target["z"], ord=np.inf)
        ),
        "casimir_intertwining_error": float(np.linalg.norm(U @ source_c @ U.conj().T - target_c, ord=np.inf)),
    }
