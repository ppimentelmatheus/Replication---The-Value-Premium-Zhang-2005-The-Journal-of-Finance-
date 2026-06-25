"""Data builders for the figures in Zhang (2005)."""

from __future__ import annotations

import numpy as np

from .calibration import Calibration
from .ff95 import FF95Result, ff95_profitability


def adjustment_cost_curve(cal: Calibration, investment_rate: np.ndarray, capital: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Figure 1 data: asymmetric adjustment cost as a function of ``i/k``."""
    investment_rate = np.asarray(investment_rate, dtype=float)
    theta = np.where(investment_rate >= 0.0, cal.gP, cal.gN)
    cost = (theta / 2.0) * investment_rate**2 * capital
    return investment_rate, cost


def pricing_kernel_moments(
    x: np.ndarray,
    beta: float,
    gamA: float,
    gamB: float,
    xbar: float,
    rhox: float,
    stdx: float,
) -> dict[str, np.ndarray]:
    """Figure 5 data, ported from ``vpCode/lambda_m.m``."""
    x = np.asarray(x, dtype=float)
    mu = (gamA + gamB * (x - xbar)) * (1.0 - rhox) * (x - xbar)
    sig = stdx * (gamA + gamB * (x - xbar))
    em = beta * np.exp(mu + 0.5 * sig**2)
    sigm = beta * np.exp(mu) * np.sqrt(np.exp(sig**2) * (np.exp(sig**2) - 1.0))
    sharpe = sigm / em
    lamb = sigm**2 / em
    return {"x": x, "em": em, "sigm": sigm, "sharpe": sharpe, "lambda": lamb}


def aggregation_quality(coefficients: np.ndarray, simh: np.ndarray, simx: np.ndarray, sigk: np.ndarray, xbar: float) -> dict[str, np.ndarray]:
    """Figure B.1 data: predicted versus actual output price and residuals."""
    coefficients = np.asarray(coefficients, dtype=float)
    simh = np.asarray(simh, dtype=float).reshape(-1)
    simx = np.asarray(simx, dtype=float).reshape(-1)
    sigk = np.asarray(sigk, dtype=float).reshape(-1)
    length = min(simh.size, simx.size, sigk.size)
    simh = simh[:length]
    simx = simx[:length]
    sigk = sigk[:length]
    predicted = coefficients[0] + coefficients[1] * simh[:-1] + coefficients[2] * (simx[:-1] - xbar) + coefficients[3] * sigk[:-1]
    actual = simh[1:]
    residual_pct = 100.0 * (actual - predicted) / actual
    return {
        "predicted_scaled": predicted / np.mean(predicted),
        "actual_scaled": actual / np.mean(actual),
        "residual_pct": residual_pct,
    }


def ff95_figure_data(Pf: np.ndarray, Bf: np.ndarray, Df: np.ndarray) -> FF95Result:
    """Figure 2 data."""
    return ff95_profitability(Pf, Bf, Df)


def investment_scatter_data(
    Pf: np.ndarray,
    Bf: np.ndarray,
    In: np.ndarray,
    sx: np.ndarray,
    gP: float,
    gN: float,
    max_points: int = 500,
) -> dict[str, np.ndarray]:
    """Figure 3 data: investment rate and adjustment cost for growth/value firms."""
    Pf = np.asarray(Pf, dtype=float)
    Bf = np.asarray(Bf, dtype=float)
    In = np.asarray(In, dtype=float)
    sx = np.asarray(sx, dtype=float).reshape(-1)
    periods = min(Pf.shape[1], Bf.shape[1], In.shape[1], sx.size)
    Pf, Bf, In, sx = Pf[:, :periods], Bf[:, :periods], In[:, :periods], sx[:periods]

    btm = Bf / Pf
    group_size = max(1, int(round(0.30 * Pf.shape[0])))
    bust_t = int(np.argmin(sx))
    boom_t = int(np.argmax(sx))

    def one_state(t: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        order = np.argsort(btm[:, t], kind="mergesort")
        growth = order[:group_size]
        value = order[-group_size:]
        growth = growth[:max_points]
        value = value[:max_points]

        ik_growth = In[growth, t] / Bf[growth, t]
        ik_value = In[value, t] / Bf[value, t]
        cost_growth = np.where(ik_growth >= 0.0, gP / 2.0, gN / 2.0) * ik_growth**2 * Bf[growth, t]
        cost_value = np.where(ik_value >= 0.0, gP / 2.0, gN / 2.0) * ik_value**2 * Bf[value, t]
        return ik_growth, cost_growth, ik_value, cost_value

    xg_bust, yg_bust, xv_bust, yv_bust = one_state(bust_t)
    xg_boom, yg_boom, xv_boom, yv_boom = one_state(boom_t)
    return {
        "bust_t": np.asarray(bust_t),
        "boom_t": np.asarray(boom_t),
        "xg_bust": xg_bust,
        "yg_bust": yg_bust,
        "xv_bust": xv_bust,
        "yv_bust": yv_bust,
        "xg_boom": xg_boom,
        "yg_boom": yg_boom,
        "xv_boom": xv_boom,
        "yv_boom": yv_boom,
    }
