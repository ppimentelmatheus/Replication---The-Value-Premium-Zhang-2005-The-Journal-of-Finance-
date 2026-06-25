"""Fama-French (1995) profitability figures ported from ``vpCode/FF95.m``."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FF95Result:
    roevec: np.ndarray
    btmvec: np.ndarray
    ROElow: np.ndarray
    ROEhigh: np.ndarray


def _lag_book(book: np.ndarray) -> np.ndarray:
    lagged = np.full_like(book, np.nan, dtype=float)
    lagged[:, 1:] = book[:, :-1]
    return lagged


def ff95_profitability(V_f: np.ndarray, B_f: np.ndarray, D_f: np.ndarray) -> FF95Result:
    """Compute the Figure 2 profitability inputs from Zhang's simulated panel.

    ``V_f``, ``B_f`` and ``D_f`` are monthly firm panels with shape ``[N, T]``.
    The output matches the economic timing in ``FF95.m``: portfolio formation is
    annual, using low and high book-to-market firms.
    """
    V_f = np.asarray(V_f, dtype=float)
    B_f = np.asarray(B_f, dtype=float)
    D_f = np.asarray(D_f, dtype=float)
    if V_f.shape != B_f.shape or V_f.shape != D_f.shape:
        raise ValueError("V_f, B_f and D_f must have the same shape")

    n, periods = V_f.shape
    if periods % 12 != 0:
        periods = (periods // 12) * 12
        V_f = V_f[:, :periods]
        B_f = B_f[:, :periods]
        D_f = D_f[:, :periods]
    if periods < 132:
        raise ValueError("FF95 profitability plots require at least 132 monthly observations")

    lagged_book = _lag_book(B_f)
    roe = (B_f + D_f) / lagged_book - 1.0
    earnings = B_f + D_f - lagged_book
    btm = B_f / V_f

    group_size = int(round(0.30 * n))
    low_slice = slice(0, group_size)
    high_slice = slice(n - group_size, n)

    roe_windows: list[list[float]] = [[] for _ in range(11)]
    btm_windows: list[list[float]] = [[] for _ in range(11)]

    for t in range(61, periods - 71):
        if (t - 1) % 12 != 0:
            continue
        order = np.argsort(btm[:, t], kind="mergesort")
        low = order[low_slice]
        high = order[high_slice]
        for offset_index, offset in enumerate(range(-60, 72, 12)):
            start = t + offset
            stop = start + 12
            denom_start = start - 1
            denom_stop = stop - 1
            low_roe = 12.0 * np.nanmean(np.nansum(earnings[low, start:stop], axis=1) / np.nansum(B_f[low, denom_start:denom_stop], axis=1))
            high_roe = 12.0 * np.nanmean(
                np.nansum(earnings[high, start:stop], axis=1) / np.nansum(B_f[high, denom_start:denom_stop], axis=1)
            )
            roe_windows[offset_index].append([low_roe, high_roe])
            btm_windows[offset_index].append(
                [np.nanmean(btm[low, start:stop]), np.nanmean(btm[high, start:stop])]
            )

    roevec = np.asarray([np.nanmean(window, axis=0) for window in roe_windows])
    btmvec = np.asarray([np.nanmean(window, axis=0) for window in btm_windows])

    ROElow = np.zeros(periods)
    ROEhigh = np.zeros(periods)
    for t in range(12, periods - 11):
        if t % 12 != 0:
            continue
        order = np.argsort(btm[:, t], kind="mergesort")
        low = order[low_slice]
        high = order[high_slice]
        ROElow[t : t + 12] = np.nanmean(roe[low, t : t + 12], axis=0)
        ROEhigh[t : t + 12] = np.nanmean(roe[high, t : t + 12], axis=0)

    ROElow = 12.0 * np.nanmean(ROElow.reshape(periods // 12, 12), axis=1)
    ROEhigh = 12.0 * np.nanmean(ROEhigh.reshape(periods // 12, 12), axis=1)
    return FF95Result(roevec=roevec, btmvec=btmvec, ROElow=ROElow, ROEhigh=ROEhigh)
