"""Interpolation helpers matching the Fortran MEX routines."""

from __future__ import annotations

import numpy as np


def interp1_extrap(x: np.ndarray, y: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Linear interpolation with linear extrapolation at both tails.

    ``x`` is one-dimensional with length ``m``. ``y`` has ``m`` rows and any
    number of trailing dimensions. ``u`` is one-dimensional.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    u = np.asarray(u, dtype=float)
    if x.ndim != 1:
        raise ValueError("x must be one-dimensional")
    if y.shape[0] != x.size:
        raise ValueError("first dimension of y must match x")

    idx = np.searchsorted(x, u, side="right") - 1
    idx = np.clip(idx, 0, x.size - 2)
    weight = (u - x[idx]) / (x[idx + 1] - x[idx])
    out = (1.0 - weight)[(...,) + (None,) * (y.ndim - 1)] * y[idx] + weight[(...,) + (None,) * (y.ndim - 1)] * y[idx + 1]
    return out


def interpolate_x(policy: np.ndarray, x_grid: np.ndarray, sx: float) -> np.ndarray:
    """Interpolate policy/value array ``[k, h, x, z]`` along aggregate state."""
    x_grid = np.asarray(x_grid, dtype=float)
    if sx > x_grid[-1]:
        return policy[:, :, -1, :] + ((sx - x_grid[-1]) / (x_grid[-1] - x_grid[-2])) * (
            policy[:, :, -1, :] - policy[:, :, -2, :]
        )
    if sx < x_grid[0]:
        return policy[:, :, 0, :] - ((x_grid[0] - sx) / (x_grid[1] - x_grid[0])) * (
            policy[:, :, 1, :] - policy[:, :, 0, :]
        )
    ix = int(np.searchsorted(x_grid, sx, side="right") - 1)
    ix = min(ix, x_grid.size - 2)
    weight_left = (x_grid[ix + 1] - sx) / (x_grid[ix + 1] - x_grid[ix])
    return weight_left * policy[:, :, ix, :] + (1.0 - weight_left) * policy[:, :, ix + 1, :]


def interpolate_h(policy_x: np.ndarray, h_grid: np.ndarray, h_value: float) -> np.ndarray:
    """Interpolate ``[k, h, z]`` along log output price, clamped at h tails."""
    h_grid = np.asarray(h_grid, dtype=float)
    if h_value <= h_grid[0]:
        return policy_x[:, 0, :]
    if h_value >= h_grid[-1]:
        return policy_x[:, -1, :]
    ih = int(np.searchsorted(h_grid, h_value, side="right") - 1)
    ih = min(ih, h_grid.size - 2)
    weight_left = (h_grid[ih + 1] - h_value) / (h_grid[ih + 1] - h_grid[ih])
    return weight_left * policy_x[:, ih, :] + (1.0 - weight_left) * policy_x[:, ih + 1, :]


def interpolate_firms(policy_kz: np.ndarray, k_grid: np.ndarray, z_grid: np.ndarray, kd: np.ndarray, zd: np.ndarray) -> np.ndarray:
    """Interpolate ``[k, z]`` values for firm-level capital and productivity."""
    k_grid = np.asarray(k_grid, dtype=float)
    z_grid = np.asarray(z_grid, dtype=float)
    kd = np.asarray(kd, dtype=float)
    zd = np.asarray(zd, dtype=float)

    k_idx = np.searchsorted(k_grid, kd, side="left")
    k_idx = np.clip(k_idx, 1, k_grid.size - 1)
    lower = k_idx - 1
    upper = k_idx
    left_weight = (k_grid[upper] - kd) / (k_grid[upper] - k_grid[lower])
    by_z = left_weight[:, None] * policy_kz[lower, :] + (1.0 - left_weight)[:, None] * policy_kz[upper, :]

    z_idx = np.searchsorted(z_grid, zd, side="right") - 1
    z_idx = np.clip(z_idx, 0, z_grid.size - 2)
    z_weight = (zd - z_grid[z_idx]) / (z_grid[z_idx + 1] - z_grid[z_idx])
    return (1.0 - z_weight) * by_z[np.arange(kd.size), z_idx] + z_weight * by_z[np.arange(kd.size), z_idx + 1]


def interpolate_policy(
    policy: np.ndarray,
    k_grid: np.ndarray,
    h_grid: np.ndarray,
    x_grid: np.ndarray,
    z_grid: np.ndarray,
    kd: np.ndarray,
    zd: np.ndarray,
    sx: float,
    h_value: float,
) -> np.ndarray:
    """Full policy interpolation used by ``simIEfcn3``, ``ssIEfcn`` and ``panIEfcn``."""
    policy_x = interpolate_x(policy, x_grid, sx)
    policy_xh = interpolate_h(policy_x, h_grid, h_value)
    return interpolate_firms(policy_xh, k_grid, z_grid, kd, zd)
