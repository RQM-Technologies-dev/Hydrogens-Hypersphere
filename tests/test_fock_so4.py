import numpy as np

from hydrogen_s3.branching import target_generators
from hydrogen_s3.fock.so4 import so4_diagnostics, so4_generators


def test_so4_commutators_casimirs_and_dimensions() -> None:
    for K in range(5):
        diagnostics = so4_diagnostics(K)
        assert diagnostics["dimension"] == (K + 1) ** 2
        assert diagnostics["j_plus"] == diagnostics["j_minus"] == K / 2
        for key in (
            "LL_error",
            "LM_error",
            "MM_error",
            "cross_casimir_error",
            "so4_casimir_error",
            "commuting_su2_error",
        ):
            assert float(diagnostics[key]) < 1e-12


def test_branched_angular_momentum_matches_existing_api() -> None:
    for K in range(4):
        branched = so4_generators(K, branched=True)["L"]
        existing = target_generators(K)
        for axis in "xyz":
            assert np.linalg.norm(branched[axis] - existing[axis], ord=np.inf) < 1e-12
