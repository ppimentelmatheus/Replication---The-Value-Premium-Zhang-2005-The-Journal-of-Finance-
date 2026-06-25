"""Simulation routines translated from the Fortran MEX files."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .calibration import Calibration, simulate_ar1
from .interpolation import interpolate_policy


@dataclass(frozen=True)
class EquilibriumSimulation:
    simh: np.ndarray
    sigk: np.ndarray
    sigz: np.ndarray
    kd: np.ndarray
    zd: np.ndarray


@dataclass(frozen=True)
class PanelSimulation:
    Pf: np.ndarray
    Bf: np.ndarray
    Df: np.ndarray
    Rf: np.ndarray
    In: np.ndarray
    iyr: float
    ikr: float
    theta: float
    dyr: float
    fyr: float
    Rm: np.ndarray
    GDPg: np.ndarray
    kd: np.ndarray
    zd: np.ndarray
    simh: np.ndarray


def antithetic_normals(n: int, rng: np.random.Generator) -> np.ndarray:
    half = n // 2
    draws = rng.standard_normal(half)
    shocks = np.concatenate([draws, -draws])
    if shocks.size < n:
        shocks = np.concatenate([shocks, rng.standard_normal(1)])
    return shocks


def update_z(zd: np.ndarray, rhoz: float, stdz: float, rng: np.random.Generator) -> np.ndarray:
    bound = 3.5 * stdz / np.sqrt(1.0 - rhoz**2)
    shocks = antithetic_normals(zd.size, rng)
    return np.clip(rhoz * zd + stdz * shocks, -bound, bound)


def initial_distribution(
    cal: Calibration,
    rng: np.random.Generator,
    n: int | None = None,
    clip: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    n = cal.N if n is None else n
    kd0 = np.ones(n)
    zd0 = cal.stdz / np.sqrt(1.0 - cal.rhoz**2) * antithetic_normals(n, rng)
    if not clip:
        return kd0, zd0
    bound = 3.5 * cal.stdz / np.sqrt(1.0 - cal.rhoz**2)
    return kd0, np.clip(zd0, -bound, bound)


def simulate_equilibrium(
    cal: Calibration,
    kd0: np.ndarray,
    zd0: np.ndarray,
    sx: np.ndarray,
    optK: np.ndarray,
    rng: np.random.Generator | None = None,
) -> EquilibriumSimulation:
    """Port of ``simIEfcn3.f90``."""
    rng = np.random.default_rng() if rng is None else rng
    kd = np.asarray(kd0, dtype=float).copy()
    zd = np.asarray(zd0, dtype=float).copy()
    periods = len(sx)
    simh = np.empty(periods)
    sigk = np.empty(periods)
    sigz = np.empty(periods)

    for t in range(periods):
        output = np.sum(np.exp(sx[t] + zd) * kd**cal.alpha)
        simh[t] = -cal.eta * np.log(output / kd.size)
        sigk[t] = np.std(kd, ddof=1)
        sigz[t] = np.std(zd, ddof=1)
        kd = interpolate_policy(optK, cal.k, cal.h, cal.x, cal.z, kd, zd, sx[t], simh[t])
        zd = update_z(zd, cal.rhoz, cal.stdz, rng)

    return EquilibriumSimulation(simh=simh, sigk=sigk, sigz=sigz, kd=kd, zd=zd)


def simulate_stationary_distribution(
    cal: Calibration,
    optK: np.ndarray,
    periods: int = 10_000,
    rng: np.random.Generator | None = None,
) -> EquilibriumSimulation:
    """Port of ``ssIE.m`` plus ``ssIEfcn.f90``."""
    rng = np.random.default_rng() if rng is None else rng
    kd0, zd0 = initial_distribution(cal, rng)
    sx, _ = simulate_ar1(cal.xbar, cal.xbar, cal.rhox, cal.stdx, periods, rng)
    return simulate_equilibrium(cal, kd0, zd0, sx, optK, rng)


def simulate_panel(
    cal: Calibration,
    kd0: np.ndarray,
    zd0: np.ndarray,
    sx: np.ndarray,
    optK: np.ndarray,
    V: np.ndarray,
    rng: np.random.Generator | None = None,
) -> PanelSimulation:
    """Port of ``panIEfcn.f90``.

    ``optK`` and ``V`` must have shape ``[nk, nh, nx, nz]``.
    """
    rng = np.random.default_rng() if rng is None else rng
    kd = np.asarray(kd0, dtype=float).copy()
    zd = np.asarray(zd0, dtype=float).copy()
    n = kd.size
    periods = len(sx)

    Bf = np.zeros((n, periods))
    In = np.zeros((n, periods))
    Rvf = np.zeros((n, periods))
    Af = np.zeros((n, periods))
    Df = np.zeros((n, periods))
    Vf = np.zeros((n, periods))
    simh = np.zeros(periods)

    for t in range(periods):
        Bf[:, t] = kd
        output = np.sum(np.exp(sx[t] + zd) * Bf[:, t] ** cal.alpha)
        simh[t] = -cal.eta * np.log(output / n)

        next_k = interpolate_policy(optK, cal.k, cal.h, cal.x, cal.z, Bf[:, t], zd, sx[t], simh[t])
        Vf[:, t] = interpolate_policy(V, cal.k, cal.h, cal.x, cal.z, Bf[:, t], zd, sx[t], simh[t])
        In[:, t] = next_k - (1.0 - cal.delta) * Bf[:, t]
        investment_rate = In[:, t] / Bf[:, t] - cal.istar
        Af[:, t] = np.where(investment_rate >= 0.0, cal.gP / 2.0, cal.gN / 2.0) * investment_rate**2 * Bf[:, t]
        Rvf[:, t] = np.exp(sx[t] + zd + simh[t]) * Bf[:, t] ** cal.alpha

        kd = next_k
        zd = update_z(zd, cal.rhoz, cal.stdz, rng)

    Df = Rvf - In - Af - cal.f
    Pf = Vf - Df
    Rf = (Pf[:, 1:] + Df[:, 1:]) / Pf[:, :-1]
    Rm = np.sum(Pf[:, :-1] * Rf, axis=0) / np.sum(Pf[:, :-1], axis=0)
    GDP = np.sum(Rvf, axis=0)
    GDPg = GDP[1:] / GDP[:-1]

    iyr = float(np.mean(np.sum(In, axis=0) / np.sum(Rvf, axis=0)))
    ikr = float(np.mean(np.sum(In, axis=0) / np.sum(Bf, axis=0)))
    theta = float(np.mean(np.sum(Af, axis=0) / np.sum(Rvf, axis=0)) / iyr)
    dyr = float(np.mean(np.sum(Df, axis=0) / np.sum(Rvf, axis=0)))
    fyr = float(np.mean(n * cal.f / np.sum(Rvf, axis=0)))

    return PanelSimulation(
        Pf=Pf,
        Bf=Bf,
        Df=Df,
        Rf=Rf,
        In=In,
        iyr=iyr,
        ikr=ikr,
        theta=theta,
        dyr=dyr,
        fyr=fyr,
        Rm=Rm,
        GDPg=GDPg,
        kd=kd,
        zd=zd,
        simh=simh,
    )
