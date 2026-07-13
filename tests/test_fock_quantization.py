import pytest

from hydrogen_s3.fock.conventions import ELECTRON_MASS_KG, PROTON_MASS_KG, CoulombSystem
from hydrogen_s3.fock.coulomb_kernel import fock_kernel_eigenvalue, quantization_from_fock
from hydrogen_s3.spectrum import RYDBERG_EV, calibrated_spectral_energy_from_K


def test_quantization_is_derived_from_kernel_eigenvalue() -> None:
    system = CoulombSystem.atomic_units()
    for K in range(6):
        result = quantization_from_fock(K, system)
        assert result.n == K + 1
        assert result.kappa == pytest.approx(1.0 / result.n, abs=1e-15)
        assert result.energy == pytest.approx(-0.5 / result.n**2, abs=1e-15)
        assert result.transformed_coupling * fock_kernel_eigenvalue(K) == pytest.approx(1.0, abs=1e-15)


def test_physical_system_and_legacy_energy_cross_check() -> None:
    physical = CoulombSystem.physical_hydrogen()
    reduced_mass_ratio = PROTON_MASS_KG / (ELECTRON_MASS_KG + PROTON_MASS_KG)
    assert physical.energy_n(1) / 1.602176634e-19 == pytest.approx(-RYDBERG_EV * reduced_mass_ratio, rel=3e-8)
    for K in range(3):
        derived_ev = quantization_from_fock(K, physical).energy / 1.602176634e-19
        assert derived_ev == pytest.approx(calibrated_spectral_energy_from_K(K) * reduced_mass_ratio, rel=3e-8)
