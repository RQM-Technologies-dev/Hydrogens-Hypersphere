"""Physical constants and the unitary momentum Fourier convention."""

from __future__ import annotations

from dataclasses import dataclass
from math import pi

# CODATA 2022 values.  The package's legacy spectrum contains only an energy
# scale, so these are the first physical mass/coupling constants in the tree.
ELECTRON_MASS_KG = 9.109_383_713_9e-31
PROTON_MASS_KG = 1.672_621_925_95e-27
ELEMENTARY_CHARGE_C = 1.602_176_634e-19
VACUUM_PERMITTIVITY_F_M = 8.854_187_8188e-12
HBAR_J_S = 1.054_571_817e-34


def _positive_integer(value: int, name: str) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True, slots=True)
class CoulombSystem:
    """Parameters for H=p^2/(2 mu)-g/r.

    ``kappa_n`` is an inverse length. ``momentum_scale_n`` is the matching
    physical momentum hbar*kappa used by the stereographic map.
    """

    mu: float
    g: float
    hbar: float

    def __post_init__(self) -> None:
        for name, value in (("mu", self.mu), ("g", self.g), ("hbar", self.hbar)):
            if not isinstance(value, (int, float)) or value <= 0.0:
                raise ValueError(f"{name} must be a positive finite scalar")

    @classmethod
    def atomic_units(cls) -> CoulombSystem:
        """Return mu=g=hbar=1."""

        return cls(mu=1.0, g=1.0, hbar=1.0)

    @classmethod
    def physical_hydrogen(cls) -> CoulombSystem:
        """Return ordinary hydrogen with the electron-proton reduced mass."""

        reduced_mass = ELECTRON_MASS_KG * PROTON_MASS_KG / (ELECTRON_MASS_KG + PROTON_MASS_KG)
        coupling = ELEMENTARY_CHARGE_C**2 / (4.0 * pi * VACUUM_PERMITTIVITY_F_M)
        return cls(mu=reduced_mass, g=coupling, hbar=HBAR_J_S)

    @property
    def bohr_radius(self) -> float:
        return self.hbar**2 / (self.mu * self.g)

    @property
    def rydberg_energy(self) -> float:
        return self.mu * self.g**2 / (2.0 * self.hbar**2)

    def kappa_n(self, n: int) -> float:
        _positive_integer(n, "n")
        return 1.0 / (float(n) * self.bohr_radius)

    def momentum_scale_n(self, n: int) -> float:
        return self.hbar * self.kappa_n(n)

    def energy_n(self, n: int) -> float:
        kappa = self.kappa_n(n)
        return -((self.hbar * kappa) ** 2) / (2.0 * self.mu)


FOURIER_FORWARD_FACTOR = "(2*pi*hbar)^(-3/2)"
COULOMB_INTEGRAL_COEFFICIENT = "g/(2*pi^2*hbar)"
