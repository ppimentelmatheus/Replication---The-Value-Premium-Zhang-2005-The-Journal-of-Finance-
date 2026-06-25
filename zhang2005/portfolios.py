"""Portfolio construction and Fama-French moments."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ValuePremiumResult:
    table: np.ndarray
    SMB: np.ndarray
    HML: np.ndarray
    SH: np.ndarray
    BH: np.ndarray
    SL: np.ndarray
    BL: np.ndarray
    SM: np.ndarray
    BM: np.ndarray


def _weighted_return(values: np.ndarray, returns: np.ndarray, ids: np.ndarray) -> np.ndarray:
    if ids.size == 0:
        return np.full(values.shape[1], np.nan)
    weights = values[ids, :]
    return np.sum(weights * returns[ids, :], axis=0) / np.sum(weights, axis=0)


def value_premium(
    V_f: np.ndarray,
    R_f: np.ndarray,
    B_f: np.ndarray,
    R: np.ndarray,
    r: np.ndarray,
    lag: int = 60,
    matlab_compat: bool = True,
) -> ValuePremiumResult:
    """Port of ``ValPrem.m``.

    With ``matlab_compat=True`` the portfolio return windows follow Zhang's
    original file exactly. Set it to ``False`` to use the rebalance month as
    the first return month, matching the timing used in ``FF93.m``.
    """
    V_f = np.asarray(V_f, dtype=float)
    R_f = np.asarray(R_f, dtype=float)
    B_f = np.asarray(B_f, dtype=float)
    R = np.asarray(R, dtype=float).reshape(-1)
    r = np.asarray(r, dtype=float).reshape(-1)

    n, periods = R_f.shape
    btm = B_f / V_f
    shifted_btm = np.zeros_like(btm)
    shifted_btm[:, 6:] = btm[:, : periods - 6]

    small = int(np.round(0.50 * n))
    low = int(np.round(0.30 * n))
    lowmed = int(np.round(0.70 * n))
    out_periods = periods - lag

    pSL = np.zeros(out_periods)
    pSM = np.zeros(out_periods)
    pSH = np.zeros(out_periods)
    pBL = np.zeros(out_periods)
    pBM = np.zeros(out_periods)
    pBH = np.zeros(out_periods)

    for t in range(lag, periods):
        if t % 12 != 0:
            continue
        size_order = np.argsort(V_f[:, t], kind="mergesort")
        idS = size_order[:small]
        idB = size_order[small:]
        btm_order = np.argsort(shifted_btm[:, t], kind="mergesort")
        idL = btm_order[:low]
        idM = btm_order[low:lowmed]
        idH = btm_order[lowmed:]

        ids = {
            "SL": np.intersect1d(idS, idL, assume_unique=False),
            "SM": np.intersect1d(idS, idM, assume_unique=False),
            "SH": np.intersect1d(idS, idH, assume_unique=False),
            "BL": np.intersect1d(idB, idL, assume_unique=False),
            "BM": np.intersect1d(idB, idM, assume_unique=False),
            "BH": np.intersect1d(idB, idH, assume_unique=False),
        }

        out_start = t - lag
        out_stop = min(out_start + 12, out_periods)
        if matlab_compat:
            data_start = out_start
        else:
            data_start = t
        data_stop = data_start + (out_stop - out_start)
        target = slice(out_start, out_stop)
        source = slice(data_start, data_stop)

        pSL[target] = _weighted_return(V_f[:, source], R_f[:, source], ids["SL"])
        pSM[target] = _weighted_return(V_f[:, source], R_f[:, source], ids["SM"])
        pSH[target] = _weighted_return(V_f[:, source], R_f[:, source], ids["SH"])
        pBL[target] = _weighted_return(V_f[:, source], R_f[:, source], ids["BL"])
        pBM[target] = _weighted_return(V_f[:, source], R_f[:, source], ids["BM"])
        pBH[target] = _weighted_return(V_f[:, source], R_f[:, source], ids["BH"])

    SMB = (pSL + pSM + pSH) / 3.0 - (pBL + pBM + pBH) / 3.0
    HML = (pSH + pBH) / 2.0 - (pSL + pBL) / 2.0

    table = np.zeros((3, 9))
    table[0, :] = 100.0 * np.array(
        [
            np.mean(R - r),
            np.nanmean(SMB),
            np.nanmean(HML),
            np.nanmean(pSL - 1.0),
            np.nanmean(pSM - 1.0),
            np.nanmean(pSH - 1.0),
            np.nanmean(pBL - 1.0),
            np.nanmean(pBM - 1.0),
            np.nanmean(pBH - 1.0),
        ]
    )
    table[1, :] = 100.0 * np.array(
        [
            np.std(R - r, ddof=1),
            np.nanstd(SMB, ddof=1),
            np.nanstd(HML, ddof=1),
            np.nanstd(pSL, ddof=1),
            np.nanstd(pSM, ddof=1),
            np.nanstd(pSH, ddof=1),
            np.nanstd(pBL, ddof=1),
            np.nanstd(pBM, ddof=1),
            np.nanstd(pBH, ddof=1),
        ]
    )
    table[2, :] = np.sqrt(periods) * table[0, :] / table[1, :]
    return ValuePremiumResult(table=table, SMB=SMB, HML=HML, SH=pSH, BH=pBH, SL=pSL, BL=pBL, SM=pSM, BM=pBM)
