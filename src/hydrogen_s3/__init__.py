"""Hydrogen S^3 Representation public API."""

from hydrogen_s3.spectrum import (
    HC_EV_NM,
    RYDBERG_EV,
    angular_labels,
    calibrated_spectral_energy_from_K,
    rydberg_shell_energy,
    s3_laplacian_eigenvalue,
    shell_dimension,
    shell_number,
    shifted_laplacian_eigenvalue,
)

__all__ = [
    "HC_EV_NM",
    "RYDBERG_EV",
    "angular_labels",
    "calibrated_spectral_energy_from_K",
    "rydberg_shell_energy",
    "s3_laplacian_eigenvalue",
    "shell_dimension",
    "shell_number",
    "shifted_laplacian_eigenvalue",
]

__version__ = "0.2.0"
