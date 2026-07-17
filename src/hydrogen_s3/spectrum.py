"""Shell-level spectral identities for the hydrogen S^3 representation.

The energy helpers implement the conventional nonrelativistic Rydberg shell
law. The numerical energy scale is supplied as an input constant; this module
does not derive the Rydberg constant, the Coulomb interaction, or radial
hydrogen wavefunctions.
"""

from __future__ import annotations

RYDBERG_EV: float = 13.605693122994
"""Rydberg energy in electronvolts used as a supplied conventional scale."""

HC_EV_NM: float = 1239.8419843320026
"""Planck constant times light speed in eV nm for wavelength baselines."""


def _validate_nonnegative_int(value: int, name: str) -> None:
    if not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer")


def _validate_positive_int(value: int, name: str) -> None:
    if not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def s3_laplacian_eigenvalue(K: int) -> int:
    """Return the scalar S^3 Laplacian eigenvalue K(K+2)."""

    _validate_nonnegative_int(K, "K")
    return K * (K + 2)


def shifted_laplacian_eigenvalue(K: int) -> int:
    """Return the shifted shell eigenvalue (-Delta_S3 + 1) = (K+1)^2."""

    _validate_nonnegative_int(K, "K")
    return (K + 1) ** 2


def shell_number(K: int) -> int:
    """Return the principal-shell label n = K + 1."""

    _validate_nonnegative_int(K, "K")
    return K + 1


def shell_dimension(K: int) -> int:
    """Return the pre-spin scalar shell dimension (K+1)^2."""

    return shifted_laplacian_eigenvalue(K)


def angular_labels(K: int) -> list[tuple[int, int]]:
    """Return all ordinary angular labels (ell, m) for a fixed shell K."""

    _validate_nonnegative_int(K, "K")
    return [(ell, m) for ell in range(K + 1) for m in range(-ell, ell + 1)]


def rydberg_shell_energy(n: int, rydberg_ev: float = RYDBERG_EV) -> float:
    """Return the conventional nonrelativistic shell energy -Ry/n^2.

    The Rydberg scale is supplied numerically. This is a calibrated reference
    energy law, not a first-principles derivation of the Coulomb Hamiltonian.
    """

    _validate_positive_int(n, "n")
    if rydberg_ev <= 0.0:
        raise ValueError("rydberg_ev must be positive")
    return -float(rydberg_ev) / float(n * n)


def calibrated_spectral_energy_from_K(K: int, rydberg_ev: float = RYDBERG_EV) -> float:
    """Package the conventional Rydberg shell law as H_spec(K) = -Ry/(K+1)^2.

    H_spec is a calibrated spectral function of the shifted S^3 Laplacian
    eigenvalue. The energy scale and physical Coulomb law are external inputs.
    """

    _validate_nonnegative_int(K, "K")
    return rydberg_shell_energy(shell_number(K), rydberg_ev=rydberg_ev)


def transition_energy_ev(n_i: int, n_f: int, rydberg_ev: float = RYDBERG_EV) -> float:
    """Return the conventional shell energy difference for n_i -> n_f."""

    _validate_positive_int(n_i, "n_i")
    _validate_positive_int(n_f, "n_f")
    if n_i <= n_f:
        raise ValueError("n_i must be greater than n_f for emission")
    if rydberg_ev <= 0.0:
        raise ValueError("rydberg_ev must be positive")
    return float(rydberg_ev) * ((1.0 / float(n_f * n_f)) - (1.0 / float(n_i * n_i)))


def vacuum_transition_wavelength_nm(
    n_i: int,
    n_f: int,
    *,
    rydberg_ev: float = RYDBERG_EV,
    hc_ev_nm: float = HC_EV_NM,
) -> float:
    """Return the conventional vacuum wavelength implied by the shell law."""

    if hc_ev_nm <= 0.0:
        raise ValueError("hc_ev_nm must be positive")
    return float(hc_ev_nm) / transition_energy_ev(n_i, n_f, rydberg_ev=rydberg_ev)
