import numpy as np

from hydrogen_s3.branching import (
    branching_diagnostics,
    branching_transform,
    source_basis,
    source_casimir,
    source_generators,
    target_basis,
    target_casimir,
    target_generators,
)
from hydrogen_s3.spectrum import shell_dimension


def _expected_casimir_spectrum(K: int) -> np.ndarray:
    values = []
    for ell in range(K + 1):
        values.extend([ell * (ell + 1)] * (2 * ell + 1))
    return np.asarray(values, dtype=float)


def test_basis_dimensions_and_multiplicities() -> None:
    for K in range(9):
        assert len(source_basis(K)) == shell_dimension(K)
        assert len(target_basis(K)) == shell_dimension(K)
        multiplicities = {ell: 0 for ell in range(K + 1)}
        for ell, _m in target_basis(K):
            multiplicities[ell] += 1
        assert multiplicities == {ell: 2 * ell + 1 for ell in range(K + 1)}


def test_branching_unitarity_through_k8() -> None:
    for K in range(9):
        U = branching_transform(K)
        identity = np.eye(shell_dimension(K), dtype=np.complex128)
        assert np.linalg.norm(U @ U.conj().T - identity, ord=np.inf) < 1e-10
        assert np.linalg.norm(U.conj().T @ U - identity, ord=np.inf) < 1e-10


def test_independent_generator_intertwining_through_k8() -> None:
    for K in range(9):
        U = branching_transform(K)
        source = source_generators(K)
        target = target_generators(K)
        for axis in ["x", "y", "z"]:
            error = np.linalg.norm(U @ source[axis] @ U.conj().T - target[axis], ord=np.inf)
            assert error < 1e-10


def test_source_and_target_commutators_through_k8() -> None:
    for K in range(9):
        for generators in [source_generators(K), target_generators(K)]:
            x = generators["x"]
            y = generators["y"]
            z = generators["z"]
            assert np.linalg.norm(x @ y - y @ x - 1.0j * z, ord=np.inf) < 1e-10
            assert np.linalg.norm(y @ z - z @ y - 1.0j * x, ord=np.inf) < 1e-10
            assert np.linalg.norm(z @ x - x @ z - 1.0j * y, ord=np.inf) < 1e-10


def test_casimir_spectra_and_intertwining_through_k8() -> None:
    for K in range(9):
        U = branching_transform(K)
        source = source_casimir(K)
        target = target_casimir(K)
        expected = _expected_casimir_spectrum(K)
        assert np.max(np.abs(np.sort(np.linalg.eigvalsh(source).real) - expected)) < 1e-10
        assert np.max(np.abs(np.sort(np.linalg.eigvalsh(target).real) - expected)) < 1e-10
        assert np.linalg.norm(U @ source @ U.conj().T - target, ord=np.inf) < 1e-10


def test_diagnostics_report_errors_separately() -> None:
    diag = branching_diagnostics(4)
    for key in [
        "unitarity_U_Udagger_error",
        "unitarity_Udagger_U_error",
        "intertwining_Jx_error",
        "intertwining_Jy_error",
        "intertwining_Jz_error",
        "casimir_intertwining_error",
    ]:
        assert key in diag
        assert float(diag[key]) < 1e-10
