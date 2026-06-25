"""Value function iteration translated from ``vfi3fcnIEccB.f90``."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .calibration import Calibration
from .interpolation import interp1_extrap


@dataclass(frozen=True)
class VFISolution:
    optK: np.ndarray
    V: np.ndarray
    I: np.ndarray
    div: np.ndarray
    iterations: int
    errV: float
    errK: float


def _h_transition(h: np.ndarray, x: np.ndarray, xbar: float, alp1: float, alp2: float, alp3: float) -> np.ndarray:
    Qh = np.zeros((h.size, h.size, x.size))
    for ix, x_now in enumerate(x):
        for ih, h_now in enumerate(h):
            hp = alp1 + alp2 * h_now + alp3 * (x_now - xbar)
            distance = np.abs(h - hp) + 1e-32
            Qh[:, ih, ix] = (1.0 / distance) / np.sum(1.0 / distance)
    return Qh


def solve_vfi(
    cal: Calibration,
    alp1: float,
    alp2: float,
    alp3: float,
    V0: np.ndarray | None = None,
    tol_v: float = 1e-3,
    tol_k: float = 1e-8,
    max_iter: int = 10_000,
    progress_every: int | None = 50,
    require_convergence: bool = True,
) -> VFISolution:
    """Solve the firm's dynamic program under the approximate law of motion.

    Arrays use shape ``[nk, nh, nx, nz]``. This is the same logical layout as
    the MATLAB code after ``reshape(optK, [nk nh nx nz])``.
    """
    k, kp, h, x, z = cal.k, cal.kp, cal.h, cal.x, cal.z
    nk, nh, nx, nz = k.size, h.size, x.size, z.size
    if V0 is None:
        V_old = np.broadcast_to(k[:, None, None, None], (nk, nh, nx, nz)).copy()
    else:
        V_old = np.asarray(V0, dtype=float).reshape((nk, nh, nx, nz))

    Qh = _h_transition(h, x, cal.xbar, alp1, alp2, alp3)
    detM = cal.beta * np.exp((cal.gamA - cal.gamB * cal.xbar) * x + cal.gamB * x**2)
    stoM = -(cal.gamA - cal.gamB * cal.xbar + cal.gamB * x)

    kk = k[:, None, None, None]
    hh = h[None, :, None, None]
    xx = x[None, None, :, None]
    zz = z[None, None, None, :]
    resources = np.exp(xx + zz + hh) * kk**cal.alpha + (1.0 - cal.delta) * kk - cal.f

    optK = np.zeros((nk, nh, nx, nz))
    optK_old = optK + 1.0
    V = V_old.copy()
    obj = np.zeros_like(V)

    for iteration in range(1, max_iter + 1):
        for ix in range(nx):
            discounted = np.exp(stoM[ix] * x)[None, None, :, None] * V_old
            expected = detM[ix] * np.einsum("kabc,ai,b,cd->kid", discounted, Qh[:, :, ix], cal.Qx[:, ix], cal.Qz)
            expected_minus_k = expected - k[:, None, None]
            expected_on_kp = interp1_extrap(k, expected_minus_k, kp)

            for ik, k_now in enumerate(k):
                investment_rate = (kp - (1.0 - cal.delta) * k_now) / k_now - cal.istar
                g = np.where(investment_rate >= 0.0, cal.gP, cal.gN)
                adjustment_cost = (g / 2.0) * investment_rate**2 * k_now
                candidate = expected_on_kp - adjustment_cost[:, None, None]
                argmax = np.argmax(candidate, axis=0)
                obj[ik, :, ix, :] = np.take_along_axis(candidate, argmax[None, :, :], axis=0)[0]
                optK[ik, :, ix, :] = kp[argmax]

        V = np.maximum(resources + obj, 1e-16)
        errK = float(np.max(np.abs(optK - optK_old)))
        errV = float(np.max(np.abs(V - V_old)))
        if progress_every and iteration % progress_every == 0:
            print(f"iteration={iteration:5d} errV={errV:.7f} errK={errK:.7f}")
        if errV <= tol_v and errK <= tol_k:
            break
        V_old = V.copy()
        optK_old = optK.copy()
    else:
        if require_convergence:
            raise RuntimeError(f"VFI did not converge after {max_iter} iterations")
        iteration = max_iter

    I = optK - (1.0 - cal.delta) * kk
    adj_rate = I / kk - cal.istar
    adj_cost = np.where(adj_rate >= 0.0, cal.gP / 2.0, cal.gN / 2.0) * adj_rate**2 * kk
    div = resources - optK - adj_cost
    return VFISolution(optK=optK, V=V, I=I, div=div, iterations=iteration, errV=errV, errK=errK)
