#!/usr/bin/env python3
"""Build data files for Zhang (2005) article figures."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zhang2005.calibration import calibrate
from zhang2005.figure_data import (
    adjustment_cost_curve,
    aggregation_quality,
    ff95_figure_data,
    investment_scatter_data,
    pricing_kernel_moments,
)
from zhang2005.io import save_npz
from zhang2005.simulation import initial_distribution, simulate_equilibrium
from zhang2005.calibration import simulate_ar1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-output-dir",
        type=Path,
        default=None,
        help="directory produced by scripts/run_full_replication.py; preferred for full-spec figures",
    )
    parser.add_argument(
        "--replication-npz",
        type=Path,
        default=Path("outputs/mini_replication_ks5/mini_replication_results.npz"),
        help="NPZ produced by scripts/run_mini_replication.py",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/figures/data"))
    parser.add_argument("--x-points", type=int, default=200)
    parser.add_argument(
        "--recompute-b1",
        action="store_true",
        help="if full output lacks equilibrium_diagnostics.npz, simulate aggregation diagnostics from saved optK",
    )
    parser.add_argument("--b1-seed", type=int, default=20250625)
    return parser.parse_args()


def _load_first_full_panel_arrays(full_output_dir: Path) -> dict[str, np.ndarray] | None:
    panel_dir = full_output_dir / "panels"
    candidates = sorted(panel_dir.glob("panel_*_arrays.npz"))
    if not candidates:
        print(f"No full panel arrays found in {panel_dir}; skipping figures 2, 3 and panel-based value spread")
        return None
    print(f"using panel arrays: {candidates[0]}")
    return dict(np.load(candidates[0]))


def _load_first_full_factor_file(full_output_dir: Path) -> dict[str, np.ndarray] | None:
    panel_dir = full_output_dir / "panels"
    candidates = sorted(panel_dir.glob("panel_*_factors.npz"))
    if not candidates:
        print(f"No panel factor files found in {panel_dir}; skipping factor-based figures")
        return None
    print(f"using panel factors: {candidates[0]}")
    return dict(np.load(candidates[0]))


def _recompute_b1_diagnostics(full_output_dir: Path, seed: int) -> dict[str, np.ndarray] | None:
    summary_path = full_output_dir / "full_replication_summary.json"
    solution_path = full_output_dir / "equilibrium_solution.npz"
    if not summary_path.exists() or not solution_path.exists():
        return None
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    config = summary.get("config", {})
    cal = calibrate(
        N=int(config.get("n_firms", 5000)),
        Ts=int(config.get("equilibrium_periods", 11000)),
        nkp=int(config.get("nkp", 5000)),
    )
    cutoff = int(config.get("cutoff", 1000))
    solution = np.load(solution_path)
    coefficients = np.asarray(solution["coefficients"], dtype=float)
    optK = np.asarray(solution["optK"], dtype=float)

    print("Recomputing Figure B.1 diagnostics from saved equilibrium solution...")
    rng = np.random.default_rng(seed)
    kd0, zd0 = initial_distribution(cal, rng, clip=False)
    sx, _ = simulate_ar1(cal.xbar, cal.xbar, cal.rhox, cal.stdx, cal.Ts, rng)
    sim = simulate_equilibrium(cal, kd0, zd0, sx, optK, rng)
    diagnostics = {
        "coefficients": coefficients,
        "eq_sx": sx,
        "eq_simh": sim.simh,
        "eq_sigk": sim.sigk,
        "eq_sigz": sim.sigz,
        "cutoff": np.asarray(cutoff),
    }
    save_npz(full_output_dir / "equilibrium_diagnostics.npz", **diagnostics)
    print(f"wrote {full_output_dir / 'equilibrium_diagnostics.npz'}")
    return diagnostics


def _load_full_data(
    full_output_dir: Path,
    recompute_b1: bool = False,
    b1_seed: int = 20250625,
) -> tuple[dict[str, np.ndarray] | None, dict[str, np.ndarray] | None, dict[str, np.ndarray] | None]:
    if not full_output_dir.exists():
        raise FileNotFoundError(f"full output directory not found: {full_output_dir}")
    summary_path = full_output_dir / "full_replication_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        config = summary.get("config", {})
        print(
            "full output config: "
            f"N={config.get('n_firms')}, panel_periods={config.get('panel_periods')}, "
            f"panel_simulations={config.get('panel_simulations')}"
        )
    diagnostics_path = full_output_dir / "equilibrium_diagnostics.npz"
    diagnostics = dict(np.load(diagnostics_path)) if diagnostics_path.exists() else None
    if diagnostics is None:
        if recompute_b1:
            diagnostics = _recompute_b1_diagnostics(full_output_dir, b1_seed)
        else:
            print("No equilibrium_diagnostics.npz found; Figure B.1 will be skipped for this full run")
    return _load_first_full_panel_arrays(full_output_dir), _load_first_full_factor_file(full_output_dir), diagnostics


def _load_panel(path: Path) -> dict[str, np.ndarray] | None:
    if not path.exists():
        print(f"Panel file not found, skipping panel figures: {path}")
        return None
    return dict(np.load(path))


def _save_figure4_proxy(data: dict[str, np.ndarray], out: Path) -> None:
    pf = data["Pf"]
    bf = data["Bf"]
    sx = data["sx"][: pf.shape[1]]
    hml = data["HML"]
    lag = pf.shape[1] - hml.size

    btm = bf / pf
    n = pf.shape[0]
    group_size = max(1, int(round(0.30 * n)))
    value_spread = np.full(pf.shape[1], np.nan)
    for t in range(pf.shape[1]):
        order = np.argsort(btm[:, t], kind="mergesort")
        growth = order[:group_size]
        value = order[-group_size:]
        value_spread[t] = np.nanmean(np.log(btm[value, t])) - np.nanmean(np.log(btm[growth, t]))

    save_npz(
        out / "figure4_simulated_spreads.npz",
        x_for_hml=sx[lag : lag + hml.size],
        hml=hml,
        x_for_value_spread=sx,
        value_spread=value_spread,
    )
    print("wrote figure4_simulated_spreads.npz")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cal = calibrate()

    inv_rate = np.linspace(-0.30, 0.30, 400)
    x_fig1, cost_fig1 = adjustment_cost_curve(cal, inv_rate)
    save_npz(args.output_dir / "figure1_adjustment_cost.npz", investment_rate=x_fig1, adjustment_cost=cost_fig1)
    print("wrote figure1_adjustment_cost.npz")

    x_grid = np.linspace(cal.xbar - 3.5 * cal.stdx / np.sqrt(1.0 - cal.rhox**2), cal.xbar + 3.5 * cal.stdx / np.sqrt(1.0 - cal.rhox**2), args.x_points)
    benchmark = pricing_kernel_moments(x_grid, cal.beta, cal.gamA, cal.gamB, cal.xbar, cal.rhox, cal.stdx)
    constant = pricing_kernel_moments(x_grid, cal.beta, cal.gamA, 0.0, cal.xbar, cal.rhox, cal.stdx)
    save_npz(
        args.output_dir / "figure5_pricing_kernel.npz",
        x=x_grid,
        benchmark_em=benchmark["em"],
        benchmark_sigm=benchmark["sigm"],
        benchmark_sharpe_annual=np.sqrt(12.0) * benchmark["sharpe"],
        constant_em=constant["em"],
        constant_sigm=constant["sigm"],
        constant_sharpe_annual=np.sqrt(12.0) * constant["sharpe"],
    )
    print("wrote figure5_pricing_kernel.npz")

    if args.full_output_dir is not None:
        data, factor_data, diagnostics = _load_full_data(args.full_output_dir, args.recompute_b1, args.b1_seed)
        if data is not None and factor_data is not None:
            data = {**factor_data, **data}
        elif data is None:
            data = factor_data
    else:
        data = _load_panel(args.replication_npz)
        factor_data = data
        diagnostics = data
    if data is None:
        return

    if {"Pf", "Bf", "Df"}.issubset(data):
        try:
            ff95 = ff95_figure_data(data["Pf"], data["Bf"], data["Df"])
        except ValueError as exc:
            print(f"Skipping Figure 2: {exc}")
        else:
            save_npz(
                args.output_dir / "figure2_ff95_profitability.npz",
                roevec=ff95.roevec,
                btmvec=ff95.btmvec,
                ROElow=ff95.ROElow,
                ROEhigh=ff95.ROEhigh,
            )
            print("wrote figure2_ff95_profitability.npz")

    if {"Pf", "Bf", "In", "sx"}.issubset(data):
        fig3 = investment_scatter_data(data["Pf"], data["Bf"], data["In"], data["sx"], cal.gP, cal.gN)
        save_npz(args.output_dir / "figure3_investment_scatter.npz", **fig3)
        print("wrote figure3_investment_scatter.npz")

    if {"HML", "Pf", "Bf", "sx"}.issubset(data):
        _save_figure4_proxy(data, args.output_dir)

    b1_data = diagnostics if diagnostics is not None else data
    if b1_data is not None and {"coefficients", "eq_simh", "eq_sx", "eq_sigk"}.issubset(b1_data):
        cutoff = int(np.asarray(b1_data.get("cutoff", 0)))
        figb1 = aggregation_quality(
            b1_data["coefficients"],
            b1_data["eq_simh"][cutoff:],
            b1_data["eq_sx"][cutoff:],
            b1_data["eq_sigk"][cutoff:],
            cal.xbar,
        )
        save_npz(args.output_dir / "figureB1_aggregation_quality.npz", **figb1)
        print("wrote figureB1_aggregation_quality.npz")
    else:
        print("Skipping Figure B.1: rerun full replication with the updated script to save equilibrium_diagnostics.npz")


if __name__ == "__main__":
    main()
