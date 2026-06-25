"""Calibration routines ported from ``vpCode/CalibCC.m``."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Calibration:
    beta: float
    gamA: float
    gamB: float
    eta: float
    alpha: float
    delta: float
    gP: float
    gN: float
    istar: float
    f: float
    nz: int
    zbar: float
    rhoz: float
    stdz: float
    zdev: float
    Qz: np.ndarray
    z: np.ndarray
    nx: int
    rhox: float
    stdx: float
    xbar: float
    xdev: float
    Qx: np.ndarray
    x: np.ndarray
    rf: np.ndarray
    k: np.ndarray
    kp: np.ndarray
    h: np.ndarray
    N: int
    Ts: int


def rouwenhorst(rho: float, mean: float, step: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Discretize an AR(1) process using Zhang's Rouwenhorst convention.

    The returned transition matrix follows the MATLAB code: columns are current
    states and rows are next states, so each column sums to one.
    """
    if n % 2 == 0:
        raise ValueError("rouwenhorst requires an odd number of grid points")

    grid = np.linspace(mean - (n - 1) / 2 * step, mean + (n - 1) / 2 * step, n)
    p = (rho + 1.0) / 2.0
    q = p
    trans = np.array(
        [
            [p**2, p * (1 - q), (1 - q) ** 2],
            [2 * p * (1 - p), p * q + (1 - p) * (1 - q), 2 * q * (1 - q)],
            [(1 - p) ** 2, (1 - p) * q, q**2],
        ],
        dtype=float,
    ).T

    while trans.shape[0] <= n - 1:
        length = trans.shape[0]
        zcol = np.zeros((length, 1))
        zrow = np.zeros((1, length))
        trans = (
            p * np.block([[trans, zcol], [zrow, np.zeros((1, 1))]])
            + (1 - p) * np.block([[zcol, trans], [np.zeros((1, 1)), zrow]])
            + (1 - q) * np.block([[zrow, np.zeros((1, 1))], [trans, zcol]])
            + q * np.block([[np.zeros((1, 1)), zrow], [zcol, trans]])
        )
        trans[1:-1, :] /= 2.0

    trans = trans.T
    if np.max(np.abs(trans.sum(axis=0) - 1.0)) >= 1e-8:
        raise RuntimeError("Rouwenhorst transition matrix columns do not sum to one")
    return trans, grid


def simulate_ar1(
    initial: float,
    mean: float,
    rho: float,
    sigma: float,
    periods: int,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, float]:
    """Port of ``CspSimu.m`` for a bounded continuous AR(1)."""
    rng = np.random.default_rng() if rng is None else rng
    sx = np.zeros(periods + 1)
    sx[0] = initial
    innovations = rng.standard_normal(periods + 1)
    upper = mean + 3.5 * sigma / np.sqrt(1.0 - rho**2)
    lower = mean - 3.5 * sigma / np.sqrt(1.0 - rho**2)

    for t in range(1, periods + 1):
        sx[t] = np.clip(mean * (1.0 - rho) + rho * sx[t - 1] + sigma * innovations[t], lower, upper)

    return sx[:-1], float(sx[-1])


def get_rf_cc(x: np.ndarray, Qx: np.ndarray, beta: float, eta_1: float, eta_2: float) -> np.ndarray:
    """Closed-form real interest rate from ``getRfcc.m``."""
    x = np.asarray(x, dtype=float)
    Qx = np.asarray(Qx, dtype=float)
    xbar = float(np.mean(x))
    rf = np.empty_like(x)
    for ix, x_now in enumerate(x):
        det = beta * np.exp((eta_1 - eta_2 * xbar) * x_now + eta_2 * x_now**2)
        sto = np.exp(-(eta_1 - eta_2 * xbar + eta_2 * x_now) * x)
        rf[ix] = 1.0 / (det * (sto @ Qx[:, ix]))
    return rf


def annualized_average_sharpe(gamA: float, stdx: float) -> float:
    """Average annualized Sharpe ratio reported in ``CalibCC.m``."""
    exp_term = np.exp(gamA**2 * stdx**2)
    monthly = np.sqrt(exp_term * (exp_term - 1.0)) / np.exp(0.5 * gamA**2 * stdx**2)
    return float(monthly * np.sqrt(12.0))


def build_capital_grid(kmin: float = 0.01, kmax: float = 10.0, include_one: bool = True) -> np.ndarray:
    """Non-equally spaced capital grid used in ``CalibCC.m``."""
    values = [kmin]
    next_value = 1.0
    index = 1
    while next_value < kmax:
        index += 1
        next_value = values[-1] + 0.005 * np.exp(0.28165 * (index - 2))
        if next_value < kmax:
            values.append(float(next_value))

    k = np.asarray(values, dtype=float)
    if include_one and not np.any(np.isclose(k, 1.0)):
        loc = int(np.argmin(np.abs(k - 1.0)))
        if k[loc] < 1.0:
            k = np.insert(k, loc + 1, 1.0)
        else:
            k = np.insert(k, loc, 1.0)
    return k


def calibrate(N: int = 5000, Ts: int = 11000, nkp: int = 5000) -> Calibration:
    """Build the benchmark monthly calibration from Zhang's MATLAB code."""
    beta = 0.994
    gamA = 50.0
    gamB = -1000.0
    eta = 0.50
    alpha = 0.30
    delta = 0.01
    gP = 15.0
    gN = 150.0
    istar = 0.0

    nz = 15
    zbar = 0.0
    rhoz = 0.97
    stdz = 0.10
    zdev = 2.0 * stdz / np.sqrt((1.0 - rhoz**2) * (nz - 1))
    Qz, z = rouwenhorst(rhoz, zbar, zdev, nz)

    nx = 11
    rhox = 0.95 ** (1.0 / 3.0)
    stdx = 0.007 / 3.0
    xbar = (
        1.0
        / (1.0 - eta)
        * np.log(
            (
                1.0
                + gP * (delta - istar)
                - beta * ((gP / 2.0) * (delta - istar) * (delta + istar) + (1.0 - delta) * (1.0 + gP * (delta - istar)))
            )
            / (
                alpha
                * beta
                * np.exp(
                    0.5 * (1.0 - eta) ** 2 * stdx**2 / (1.0 - rhox**2)
                    + 0.5 * (1.0 - eta) ** 2 * stdz**2 / (1.0 - rhoz**2)
                )
            )
        )
    )
    xdev = 2.0 * stdx / np.sqrt((1.0 - rhox**2) * (nx - 1))
    Qx, x = rouwenhorst(rhox, xbar, xdev, nx)

    f = 0.0365
    rf = get_rf_cc(x, Qx, beta, gamA, gamB)
    k = build_capital_grid()
    kp = np.linspace(float(k.min()), float(k.max()), nkp)
    h = np.linspace(2.75, 3.25, 5)

    return Calibration(
        beta=beta,
        gamA=gamA,
        gamB=gamB,
        eta=eta,
        alpha=alpha,
        delta=delta,
        gP=gP,
        gN=gN,
        istar=istar,
        f=f,
        nz=nz,
        zbar=zbar,
        rhoz=rhoz,
        stdz=stdz,
        zdev=zdev,
        Qz=Qz,
        z=z,
        nx=nx,
        rhox=rhox,
        stdx=stdx,
        xbar=float(xbar),
        xdev=xdev,
        Qx=Qx,
        x=x,
        rf=rf,
        k=k,
        kp=kp,
        h=h,
        N=N,
        Ts=Ts,
    )
