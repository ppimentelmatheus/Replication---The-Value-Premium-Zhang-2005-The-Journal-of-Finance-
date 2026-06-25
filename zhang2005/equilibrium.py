"""Krusell-Smith equilibrium loop ported from ``mainCC.m``."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .calibration import Calibration, simulate_ar1
from .linalg import ols
from .simulation import EquilibriumSimulation, initial_distribution, simulate_equilibrium
from .vfi import VFISolution, solve_vfi


@dataclass(frozen=True)
class EquilibriumResult:
    coefficients: np.ndarray
    r2: float
    sighat: float
    tstats: np.ndarray
    residuals: np.ndarray
    simulation: EquilibriumSimulation
    solution: VFISolution
    sx: np.ndarray
    iterations: int
    error: float


def construct_equilibrium(
    cal: Calibration,
    initial_coefficients: tuple[float, float, float, float] = (0.0, 1.0, 0.0, 0.0),
    initial_value_function: np.ndarray | None = None,
    cutoff: int = 1000,
    coefficient_tol: float = 1e-2,
    max_iterations: int = 100,
    vfi_options: dict | None = None,
    require_convergence: bool = True,
    rng: np.random.Generator | None = None,
) -> EquilibriumResult:
    """Construct the industry equilibrium using Zhang's approximate aggregation loop.

    This follows ``mainCC.m``. The VFI itself uses only ``alp1``, ``alp2`` and
    ``alp3``; ``alp4`` enters the regression/convergence check exactly as in the
    MATLAB driver.
    """
    if cutoff >= cal.Ts - 1:
        raise ValueError("cutoff must leave at least two simulated observations")

    rng = np.random.default_rng() if rng is None else rng
    vfi_options = {} if vfi_options is None else dict(vfi_options)
    alp1, alp2, alp3, alp4 = initial_coefficients
    kd0, zd0 = initial_distribution(cal, rng, clip=False)
    sx, _ = simulate_ar1(cal.xbar, cal.xbar, cal.rhox, cal.stdx, cal.Ts, rng)
    simx = sx[cutoff:]
    V0 = None if initial_value_function is None else np.asarray(initial_value_function, dtype=float)

    last_solution: VFISolution | None = None
    last_simulation: EquilibriumSimulation | None = None
    last_tstats = np.full(4, np.nan)
    last_residuals = np.array([])
    last_r2 = np.nan
    last_sighat = np.nan
    err = np.inf

    for iteration in range(1, max_iterations + 1):
        solution = solve_vfi(cal, alp1, alp2, alp3, V0=V0, **vfi_options)
        simulation = simulate_equilibrium(cal, kd0, zd0, sx, solution.optK, rng)

        simh = simulation.simh[cutoff:]
        sigk = simulation.sigk[cutoff:]
        y = simh[1:]
        xreg = np.column_stack(
            [
                np.ones(len(simh) - 1),
                simh[:-1],
                simx[:-1] - cal.xbar,
                sigk[:-1],
            ]
        )
        coef, tstats, r2, sighat, residuals = ols(y, xreg)
        new_coefficients = np.asarray(coef[:4], dtype=float)
        old_coefficients = np.asarray([alp1, alp2, alp3, alp4], dtype=float)
        err = float(np.max(np.abs(new_coefficients - old_coefficients)))

        alp1, alp2, alp3, alp4 = new_coefficients
        V0 = solution.V
        last_solution = solution
        last_simulation = simulation
        last_tstats = tstats
        last_residuals = residuals
        last_r2 = r2
        last_sighat = sighat

        if err < coefficient_tol:
            break
    else:
        if not require_convergence:
            iteration = max_iterations
        else:
            raise RuntimeError(f"equilibrium coefficients did not converge after {max_iterations} iterations")

    assert last_solution is not None
    assert last_simulation is not None
    return EquilibriumResult(
        coefficients=np.asarray([alp1, alp2, alp3, alp4]),
        r2=float(last_r2),
        sighat=float(last_sighat),
        tstats=np.asarray(last_tstats),
        residuals=np.asarray(last_residuals),
        simulation=last_simulation,
        solution=last_solution,
        sx=sx,
        iterations=iteration,
        error=err,
    )
