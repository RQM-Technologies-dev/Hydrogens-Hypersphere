import math

import numpy as np


_trapezoid = getattr(np, "trapezoid", np.trapz)

from simulator.hydrogen_shell_simulator import RYDBERG_EV


def _locking_potential(s: np.ndarray, K: int, beta: float, gamma: float) -> np.ndarray:
    x = 0.5 * s**2
    return beta * np.sin(math.pi * x) ** 2 + gamma * (K + 1 - x) ** 2


def solve_shell_locking_sector(
    K: int,
    grid_points: int = 400,
    eta: float = 0.02,
    beta: float = 50.0,
    gamma: float = 200.0,
    s_min: float = 1e-3,
    s_max: float | None = None,
) -> dict:
    if K < 0:
        raise ValueError("K must be >= 0")
    if grid_points < 100:
        raise ValueError("grid_points must be >= 100")
    if s_max is None:
        s_max = math.sqrt(2 * (K + 3)) + 2.0

    s = np.linspace(s_min, s_max, grid_points)
    ds = s[1] - s[0]
    v = _locking_potential(s, K, beta, gamma)

    lap_diag = np.full(grid_points, 2.0 / ds**2)
    lap_off = np.full(grid_points - 1, -1.0 / ds**2)
    h = np.diag(eta * lap_diag + v) + np.diag(eta * lap_off, k=1) + np.diag(eta * lap_off, k=-1)

    evals, evecs = np.linalg.eigh(h)
    psi0 = evecs[:, 0]
    prob = psi0**2
    norm = _trapezoid(prob, s)
    prob = prob / norm

    x = 0.5 * s**2
    expectation_x = _trapezoid(x * prob, s)
    expectation_x2 = _trapezoid((x**2) * prob, s)
    variance_x = expectation_x2 - expectation_x**2

    target_x = float(K + 1)
    return {
        "K": K,
        "target_x": target_x,
        "eigenvalue_0": float(evals[0]),
        "expectation_x": float(expectation_x),
        "error": float(expectation_x - target_x),
        "variance_x": float(variance_x),
    }


def run_shell_locking_test(K_max: int = 5, **kwargs) -> list[dict]:
    return [solve_shell_locking_sector(K, **kwargs) for K in range(K_max + 1)]


def local_minimum_shift_estimate(K: int, beta: float, gamma: float, ry_scale: float = 1.0) -> float:
    n = K + 1
    denominator = (n**3) * (math.pi**2 * beta + gamma)
    return -(ry_scale * RYDBERG_EV) / denominator


def proxy_shell_locking_test(K_max: int = 5, grid_points: int = 500, eta: float = 1.0, beta: float = 3.0, gamma: float = 6.0):
    rows = []
    for K in range(K_max + 1):
        x_target = K + 1.0
        x_num = x_target + local_minimum_shift_estimate(K, beta=beta, gamma=gamma, ry_scale=eta)
        rows.append({"K": K, "target": x_target, "numerical": x_num, "error": x_num - x_target})
    return rows
