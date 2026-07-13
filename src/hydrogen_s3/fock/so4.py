"""Fixed-shell SO(4) generators assembled from existing spin matrices."""

from __future__ import annotations

import numpy as np

from hydrogen_s3.angular_momentum import ComplexMatrix, jx, jy, jz
from hydrogen_s3.branching import branching_transform, shell_j


def _dual(matrix: ComplexMatrix) -> ComplexMatrix:
    return np.asarray(-matrix.T, dtype=np.complex128)


def so4_generators(K: int, *, branched: bool = False) -> dict[str, dict[str, ComplexMatrix]]:
    """Return dimensionless L, M, J_plus, and J_minus generator triples."""

    j = shell_j(K)
    dimension = K + 1
    identity = np.eye(dimension, dtype=np.complex128)
    spin = {"x": jx(j), "y": jy(j), "z": jz(j)}
    plus = {axis: np.kron(matrix, identity) for axis, matrix in spin.items()}
    minus = {axis: np.kron(identity, _dual(matrix)) for axis, matrix in spin.items()}
    angular = {axis: plus[axis] + minus[axis] for axis in spin}
    runge_lenz = {axis: plus[axis] - minus[axis] for axis in spin}
    result = {"L": angular, "M": runge_lenz, "J_plus": plus, "J_minus": minus}
    if not branched:
        return result
    transform = branching_transform(K)
    return {
        family: {axis: transform @ matrix @ transform.conj().T for axis, matrix in generators.items()}
        for family, generators in result.items()
    }


def so4_diagnostics(K: int) -> dict[str, float | int]:
    generators = so4_generators(K)
    angular = generators["L"]
    runge = generators["M"]
    plus = generators["J_plus"]
    minus = generators["J_minus"]
    identity = np.eye((K + 1) ** 2, dtype=np.complex128)

    def commutator(left: ComplexMatrix, right: ComplexMatrix) -> ComplexMatrix:
        return left @ right - right @ left

    l2 = sum((angular[axis] @ angular[axis] for axis in "xyz"), start=np.zeros_like(identity))
    m2 = sum((runge[axis] @ runge[axis] for axis in "xyz"), start=np.zeros_like(identity))
    dot = sum((angular[axis] @ runge[axis] for axis in "xyz"), start=np.zeros_like(identity))
    j = K / 2.0
    return {
        "K": K,
        "dimension": (K + 1) ** 2,
        "j_plus": j,
        "j_minus": j,
        "LL_error": float(np.linalg.norm(commutator(angular["x"], angular["y"]) - 1.0j * angular["z"], ord=np.inf)),
        "LM_error": float(np.linalg.norm(commutator(angular["x"], runge["y"]) - 1.0j * runge["z"], ord=np.inf)),
        "MM_error": float(np.linalg.norm(commutator(runge["x"], runge["y"]) - 1.0j * angular["z"], ord=np.inf)),
        "cross_casimir_error": float(np.linalg.norm(dot, ord=np.inf)),
        "so4_casimir_error": float(np.linalg.norm(l2 + m2 - K * (K + 2) * identity, ord=np.inf)),
        "commuting_su2_error": float(np.linalg.norm(commutator(plus["x"], minus["y"]), ord=np.inf)),
    }
