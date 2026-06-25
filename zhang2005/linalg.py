"""Small regression helpers used by the replication code."""

from __future__ import annotations

import numpy as np


def ols(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float, np.ndarray]:
    """Port of ``vpCode/ols.m``.

    Returns coefficients, t-statistics, R-squared, residual standard deviation,
    and residuals.
    """
    y = np.asarray(y, dtype=float).reshape(-1)
    x = np.asarray(x, dtype=float)
    beta = np.linalg.solve(x.T @ x, x.T @ y)
    resid = y - x @ beta
    sighat = float(np.std(resid, ddof=1))
    r2 = 1.0 - float(np.var(resid, ddof=1) / np.var(y, ddof=1))
    cov = (resid @ resid) / (len(y) - x.shape[1]) * np.linalg.inv(x.T @ x)
    tstats = beta / np.sqrt(np.diag(cov))
    return beta, tstats, r2, sighat, resid
