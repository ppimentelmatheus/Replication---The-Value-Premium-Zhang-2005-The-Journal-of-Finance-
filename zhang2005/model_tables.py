"""Model-only table builders for Zhang (2005)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .calibration import Calibration, annualized_average_sharpe
from .linalg import ols

MIN_PRICE = 1e-6


@dataclass(frozen=True)
class PortfolioSortResult:
    returns: np.ndarray
    book_to_market: np.ndarray
    value_spread: np.ndarray
    labels: list[str]


def parameter_rows(cal: Calibration) -> list[dict[str, str | float]]:
    """Model parameter table, without empirical columns."""
    rows = [
        ("beta", "Subjective discount factor", cal.beta),
        ("gamma0", "Pricing kernel loading", cal.gamA),
        ("gamma1", "Countercyclical price-of-risk loading", cal.gamB),
        ("eta", "Inverse demand curvature", cal.eta),
        ("alpha", "Capital share", cal.alpha),
        ("delta", "Depreciation rate", cal.delta),
        ("theta_plus", "Adjustment cost for positive investment", cal.gP),
        ("theta_minus", "Adjustment cost for negative investment", cal.gN),
        ("f", "Fixed operating cost", cal.f),
        ("rho_x", "Aggregate productivity persistence", cal.rhox),
        ("sigma_x", "Aggregate productivity volatility", cal.stdx),
        ("rho_z", "Idiosyncratic productivity persistence", cal.rhoz),
        ("sigma_z", "Idiosyncratic productivity volatility", cal.stdz),
        ("N", "Number of simulated firms", cal.N),
        ("T", "Number of simulated monthly periods", cal.Ts),
    ]
    return [{"parameter": name, "description": desc, "value": value} for name, desc, value in rows]


def _safe_btm(Pf: np.ndarray, Bf: np.ndarray, min_price: float = MIN_PRICE) -> np.ndarray:
    return np.where((Pf > min_price) & np.isfinite(Pf) & np.isfinite(Bf), Bf / Pf, np.nan)


def _safe_industry_btm(Pf: np.ndarray, Bf: np.ndarray, min_price: float = MIN_PRICE) -> np.ndarray:
    valid = (Pf > min_price) & np.isfinite(Pf) & np.isfinite(Bf)
    book = np.where(valid, Bf, np.nan)
    value = np.where(valid, Pf, np.nan)
    return np.nansum(book, axis=0) / np.nansum(value, axis=0)


def aggregate_moment_rows(data: dict[str, np.ndarray]) -> list[dict[str, str | float]]:
    Pf, Bf, Rm, srf = data["Pf"], data["Bf"], data["Rm"], data["srf"]
    Rf = data["Rf"]
    HML = data.get("HML", np.array([]))
    SMB = data.get("SMB", np.array([]))

    industry_btm = _safe_industry_btm(Pf, Bf)
    firm_btm = _safe_btm(Pf, Bf)
    market_excess = Rm - srf
    rows = [
        ("Mean market excess return, annual pct", 1200.0 * np.mean(market_excess)),
        ("Market return volatility, annual pct", 100.0 * np.sqrt(12.0) * np.std(Rm, ddof=1)),
        ("Mean industry book-to-market", np.mean(industry_btm)),
        ("Volatility of industry book-to-market", np.std(industry_btm, ddof=1)),
        ("Mean cross-sectional B/M dispersion", np.nanmean(np.nanstd(firm_btm, axis=0, ddof=1))),
        ("Mean firm return, annual pct", 1200.0 * np.mean(Rf - 1.0)),
        ("Mean firm return volatility, annual pct", 100.0 * np.sqrt(12.0) * np.mean(np.std(Rf, axis=1, ddof=1))),
    ]
    if HML.size:
        rows.append(("HML mean return, annual pct", 1200.0 * np.nanmean(HML)))
        rows.append(("HML volatility, annual pct", 100.0 * np.sqrt(12.0) * np.nanstd(HML, ddof=1)))
    if SMB.size:
        rows.append(("SMB mean return, annual pct", 1200.0 * np.nanmean(SMB)))
        rows.append(("SMB volatility, annual pct", 100.0 * np.sqrt(12.0) * np.nanstd(SMB, ddof=1)))
    return [{"moment": name, "value": value} for name, value in rows]


def paper_table2_rows(cal: Calibration, data: dict[str, np.ndarray]) -> list[dict[str, str | float]]:
    """Moments matching Table II's model column definitions as closely as possible."""
    Rm = np.asarray(data["Rm"], dtype=float).reshape(-1)
    srf = np.asarray(data["srf"], dtype=float).reshape(-1)
    rows: list[tuple[str, float]] = [
        ("Average annual Sharpe ratio", annualized_average_sharpe(cal.gamA, cal.stdx)),
        ("Average annual real interest rate", float(np.nanmean(srf) ** 12.0 - 1.0)),
        ("Annual volatility of real interest rate", float(np.nanstd(srf, ddof=1) * np.sqrt(12.0))),
        ("Average annual value-weighted industry return", float(np.nanmean(Rm) ** 12.0 - 1.0)),
        ("Annual volatility of value-weighted industry return", float(np.nanstd(Rm, ddof=1) * np.sqrt(12.0))),
    ]

    if "Rf" in data:
        Rf = np.asarray(data["Rf"], dtype=float)
        rows.append(("Average volatility of individual stock return", float(np.nanmean(np.nanstd(Rf, axis=1, ddof=1)) * np.sqrt(12.0))))

    if {"Pf", "Bf"}.issubset(data):
        Pf = np.asarray(data["Pf"], dtype=float)
        Bf = np.asarray(data["Bf"], dtype=float)
        industry_btm = _safe_industry_btm(Pf, Bf)
        rows.append(("Average industry book-to-market ratio", float(np.nanmean(industry_btm))))
        rows.append(("Volatility of industry book-to-market ratio", float(np.nanstd(industry_btm, ddof=1))))

    if {"In", "Bf"}.issubset(data):
        In = np.asarray(data["In"], dtype=float)
        Bf = np.asarray(data["Bf"], dtype=float)
        investment_rate = In / Bf
        rows.append(("Annual average rate of investment", float(12.0 * np.nanmean(investment_rate[investment_rate > 0.0]))))
        rows.append(("Annual average rate of disinvestment", float(abs(12.0 * np.nanmean(investment_rate[investment_rate < 0.0])))))

    return [{"moment": name, "model": value} for name, value in rows]


def sort_by_book_to_market(
    Pf: np.ndarray,
    Rf: np.ndarray,
    Bf: np.ndarray,
    n_portfolios: int = 10,
    min_price: float = MIN_PRICE,
) -> PortfolioSortResult:
    """Value-weighted monthly returns for portfolios sorted on book-to-market."""
    Pf = np.asarray(Pf, dtype=float)
    Rf = np.asarray(Rf, dtype=float)
    Bf = np.asarray(Bf, dtype=float)
    periods = min(Pf.shape[1], Rf.shape[1], Bf.shape[1])
    Pf, Rf, Bf = Pf[:, :periods], Rf[:, :periods], Bf[:, :periods]
    n = Pf.shape[0]
    labels = [f"BM{j + 1}" for j in range(n_portfolios)]
    returns = np.full((n_portfolios, periods), np.nan)
    btm_series = np.full((n_portfolios, periods), np.nan)
    value_spread = np.full(periods, np.nan)

    btm = _safe_btm(Pf, Bf, min_price=min_price)
    for t in range(periods):
        valid = np.isfinite(btm[:, t]) & np.isfinite(Rf[:, t]) & np.isfinite(Pf[:, t])
        if np.sum(valid) < n_portfolios:
            continue
        order = np.argsort(btm[valid, t], kind="mergesort")
        valid_ids = np.flatnonzero(valid)
        ordered_ids = valid_ids[order]
        groups = np.array_split(ordered_ids, n_portfolios)
        for j, ids in enumerate(groups):
            weights = Pf[ids, t]
            returns[j, t] = np.sum(weights * Rf[ids, t]) / np.sum(weights)
            btm_series[j, t] = np.mean(btm[ids, t])
        value_spread[t] = np.log(btm_series[-1, t]) - np.log(btm_series[0, t])

    return PortfolioSortResult(returns=returns, book_to_market=btm_series, value_spread=value_spread, labels=labels)


def portfolio_summary_rows(sort: PortfolioSortResult, Rm: np.ndarray, srf: np.ndarray) -> list[dict[str, str | float]]:
    periods = min(sort.returns.shape[1], Rm.size, srf.size)
    portfolio_returns = sort.returns[:, :periods]
    market_excess = Rm[:periods] - srf[:periods]
    rows: list[dict[str, str | float]] = []
    for j, label in enumerate(sort.labels):
        gross = portfolio_returns[j]
        excess = gross - srf[:periods]
        xreg = np.column_stack([np.ones(periods), market_excess])
        valid = np.isfinite(excess) & np.all(np.isfinite(xreg), axis=1)
        beta, tstats, r2, _, _ = ols(excess[valid], xreg[valid])
        rows.append(
            {
                "portfolio": label,
                "mean_return_annual_pct": 100.0 * (np.nanmean(gross) ** 12.0 - 1.0),
                "volatility_annual_pct": 100.0 * np.sqrt(12.0) * np.nanstd(gross, ddof=1),
                "market_beta": beta[1],
                "market_beta_tstat": tstats[1],
                "capm_r2": r2,
                "mean_book_to_market": np.nanmean(sort.book_to_market[j, :periods]),
            }
        )

    hml = portfolio_returns[-1] - portfolio_returns[0]
    xreg = np.column_stack([np.ones(periods), market_excess])
    valid = np.isfinite(hml) & np.all(np.isfinite(xreg), axis=1)
    beta, tstats, r2, _, _ = ols(hml[valid], xreg[valid])
    rows.append(
        {
            "portfolio": "HML10",
            "mean_return_annual_pct": 1200.0 * np.nanmean(hml),
            "volatility_annual_pct": 100.0 * np.sqrt(12.0) * np.nanstd(hml, ddof=1),
            "market_beta": beta[1],
            "market_beta_tstat": tstats[1],
            "capm_r2": r2,
            "mean_book_to_market": np.nan,
        }
    )
    return rows


def predictive_regression_rows(data: dict[str, np.ndarray], sort: PortfolioSortResult, xbar: float) -> list[dict[str, str | float]]:
    Pf, Bf, Rm, sx = data["Pf"], data["Bf"], data["Rm"], data["sx"]
    periods = min(sort.returns.shape[1], Rm.size, sx.size, Pf.shape[1], Bf.shape[1])
    hml = sort.returns[-1, :periods] - sort.returns[0, :periods]
    value_spread = sort.value_spread[:periods]
    industry_btm = _safe_industry_btm(Pf[:, :periods], Bf[:, :periods])
    x_demeaned = sx[:periods] - xbar
    regressions = [
        ("HML on value spread", hml, np.column_stack([np.ones(periods), value_spread]), ["const", "value_spread"]),
        ("HML on aggregate productivity", hml, np.column_stack([np.ones(periods), x_demeaned]), ["const", "x_minus_xbar"]),
        (
            "HML on value spread and aggregate productivity",
            hml,
            np.column_stack([np.ones(periods), value_spread, x_demeaned]),
            ["const", "value_spread", "x_minus_xbar"],
        ),
        ("Market return on industry B/M", Rm[:periods], np.column_stack([np.ones(periods), np.log(industry_btm)]), ["const", "log_industry_btm"]),
        ("Market return on value spread", Rm[:periods], np.column_stack([np.ones(periods), value_spread]), ["const", "value_spread"]),
    ]
    rows: list[dict[str, str | float]] = []
    for name, y, x, labels in regressions:
        valid = np.isfinite(y) & np.all(np.isfinite(x), axis=1)
        y = y[valid]
        x = x[valid]
        if y.size <= x.shape[1]:
            continue
        coef, tstats, r2, _, _ = ols(y, x)
        row: dict[str, str | float] = {"regression": name, "r2": r2}
        for label, c, t in zip(labels, coef, tstats, strict=True):
            row[f"{label}_coef"] = c
            row[f"{label}_tstat"] = t
        rows.append(row)
    return rows
