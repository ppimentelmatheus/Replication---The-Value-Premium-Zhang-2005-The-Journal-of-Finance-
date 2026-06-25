#!/usr/bin/env python3
"""Run a small end-to-end Zhang (2005) replication pipeline.

This script is intentionally small. It is a plumbing check before running the
full benchmark: calibration, approximate equilibrium, panel simulation, and
value-premium portfolio construction.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zhang2005.calibration import Calibration, calibrate, get_rf_cc, rouwenhorst, simulate_ar1
from zhang2005.equilibrium import construct_equilibrium
from zhang2005.interpolation import interp1_extrap
from zhang2005.io import load_mat, save_npz
from zhang2005.portfolios import value_premium
from zhang2005.simulation import simulate_panel


def _odd_grid_size(value: int, name: str) -> int:
    if value < 3 or value % 2 == 0:
        raise argparse.ArgumentTypeError(f"{name} must be an odd integer >= 3")
    return value


def make_mini_calibration(
    n_firms: int,
    periods: int,
    nkp: int,
    nx: int,
    nz: int,
    capital_grid: str,
) -> Calibration:
    """Build a reduced calibration that preserves Zhang's primitive parameters."""
    base = calibrate(N=n_firms, Ts=periods, nkp=nkp)

    qx, x = rouwenhorst(base.rhox, base.xbar, base.xdev, nx)
    qz, z = rouwenhorst(base.rhoz, base.zbar, base.zdev, nz)
    if capital_grid == "benchmark":
        k_grid = base.k
        kp_grid = np.linspace(float(base.k.min()), float(base.k.max()), nkp)
    elif capital_grid == "compact":
        k_grid = np.array([0.50, 0.75, 1.00, 1.30, 1.70], dtype=float)
        kp_grid = np.linspace(float(k_grid.min()), float(k_grid.max()), nkp)
    else:
        raise ValueError(f"unknown capital grid preset: {capital_grid}")

    return replace(
        base,
        N=n_firms,
        Ts=periods,
        nx=nx,
        nz=nz,
        Qx=qx,
        Qz=qz,
        x=x,
        z=z,
        rf=get_rf_cc(x, qx, base.beta, base.gamA, base.gamB),
        k=k_grid,
        kp=kp_grid,
        h=np.array([2.75, 3.00, 3.25], dtype=float),
    )


def summarize_table(table: np.ndarray) -> dict[str, dict[str, float]]:
    columns = ["Mkt-Rf", "SMB", "HML", "SL", "SM", "SH", "BL", "BM", "BH"]
    rows = ["mean_pct", "std_pct", "t_stat"]
    return {
        row: {col: float(table[i, j]) for j, col in enumerate(columns)}
        for i, row in enumerate(rows)
    }


def load_initial_coefficients(path: Path | None) -> tuple[float, float, float, float]:
    if path is None:
        return (0.0, 1.0, 0.0, 0.0)
    data = load_mat(path)
    return tuple(float(data[name]) for name in ("alp1", "alp2", "alp3", "alp4"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-firms", type=int, default=120, help="number of simulated firms")
    parser.add_argument("--periods", type=int, default=96, help="number of monthly panel periods")
    parser.add_argument("--nkp", type=int, default=201, help="candidate capital grid points")
    parser.add_argument(
        "--capital-grid",
        choices=("benchmark", "compact"),
        default="benchmark",
        help="capital grid preset: full benchmark k grid or a tiny compact grid",
    )
    parser.add_argument("--nx", type=lambda x: _odd_grid_size(int(x), "nx"), default=3, help="aggregate shock grid size")
    parser.add_argument("--nz", type=lambda x: _odd_grid_size(int(x), "nz"), default=3, help="idiosyncratic shock grid size")
    parser.add_argument("--cutoff", type=int, default=12, help="burn-in observations for the law-of-motion regression")
    parser.add_argument("--ks-iterations", type=int, default=1, help="Krusell-Smith coefficient iterations")
    parser.add_argument("--vfi-iterations", type=int, default=120, help="value-function iterations per KS iteration")
    parser.add_argument(
        "--initial-coefficients",
        type=Path,
        default=Path("vpCode/coefIS.mat"),
        help="MAT file with alp1-alp4; use 'none' to start from (0,1,0,0)",
    )
    parser.add_argument("--lag", type=int, default=60, help="portfolio formation lag in months")
    parser.add_argument("--seed", type=int, default=20250622, help="random seed")
    parser.add_argument("--matlab-compat", action="store_true", help="use the original ValPrem.m return timing")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"), help="directory for saved outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.periods <= args.lag + 12:
        raise ValueError("--periods must exceed --lag by at least 13 months")

    rng = np.random.default_rng(args.seed)
    coefficient_path = None if str(args.initial_coefficients).lower() == "none" else args.initial_coefficients
    initial_coefficients = load_initial_coefficients(coefficient_path)
    cal = make_mini_calibration(args.n_firms, args.periods, args.nkp, args.nx, args.nz, args.capital_grid)

    print("Mini replication setup")
    print(f"  firms: {cal.N}")
    print(f"  periods: {cal.Ts}")
    print(f"  grids: nk={cal.k.size}, nkp={cal.kp.size}, nh={cal.h.size}, nx={cal.nx}, nz={cal.nz}")
    print(f"  initial coefficients: {np.array2string(np.asarray(initial_coefficients), precision=6)}")

    print("\nSolving approximate equilibrium...")
    eq = construct_equilibrium(
        cal,
        initial_coefficients=initial_coefficients,
        cutoff=args.cutoff,
        coefficient_tol=0.0,
        max_iterations=args.ks_iterations,
        vfi_options={
            "max_iter": args.vfi_iterations,
            "tol_v": 1e-3,
            "tol_k": 1e-8,
            "progress_every": None,
            "require_convergence": False,
        },
        require_convergence=False,
        rng=rng,
    )
    print(f"  KS iterations: {eq.iterations}")
    print(f"  coefficients: {np.array2string(eq.coefficients, precision=6)}")
    print(f"  R2: {eq.r2:.6f}")
    print(f"  final coefficient error: {eq.error:.6g}")
    print(f"  VFI iterations in final step: {eq.solution.iterations}")
    print(f"  VFI final errV/errK: {eq.solution.errV:.6g} / {eq.solution.errK:.6g}")

    print("\nSimulating firm panel...")
    sx, _ = simulate_ar1(cal.xbar, cal.xbar, cal.rhox, cal.stdx, cal.Ts, rng)
    panel = simulate_panel(cal, eq.simulation.kd, eq.simulation.zd, sx, eq.solution.optK, eq.solution.V, rng)

    pf = np.maximum(panel.Pf[:, :-1], 1e-8)
    bf = panel.Bf[:, :-1]
    rf = panel.Rf
    rm = panel.Rm
    srf = interp1_extrap(cal.x, cal.rf, sx).reshape(-1)[:-1]

    print("  panel arrays:")
    print(f"    Pf: {pf.shape}")
    print(f"    Rf: {rf.shape}")
    print(f"    Bf: {bf.shape}")
    print(f"    Rm: {rm.shape}")
    print(f"    Bf min/mean/max: {bf.min():.6f} / {bf.mean():.6f} / {bf.max():.6f}")
    print(f"    Bf cross-sectional std, avg over time: {np.mean(np.std(bf, axis=0, ddof=1)):.6g}")

    print("\nConstructing value/growth portfolios...")
    vp = value_premium(pf, rf, bf, rm, srf, lag=args.lag, matlab_compat=args.matlab_compat)
    print("  Table rows: percent mean, percent std, t-stat")
    print(np.array2string(vp.table, precision=4, suppress_small=False))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = args.output_dir / "mini_replication_results.npz"
    json_path = args.output_dir / "mini_replication_summary.json"
    save_npz(
        npz_path,
        coefficients=eq.coefficients,
        tstats=eq.tstats,
        residuals=eq.residuals,
        optK=eq.solution.optK,
        V=eq.solution.V,
        Pf=pf,
        Bf=bf,
        Df=panel.Df[:, :-1],
        In=panel.In[:, :-1],
        Rf=rf,
        Rm=rm,
        srf=srf,
        sx=sx,
        eq_sx=eq.sx,
        eq_simh=eq.simulation.simh,
        eq_sigk=eq.simulation.sigk,
        eq_sigz=eq.simulation.sigz,
        cutoff=np.asarray(args.cutoff),
        panel_simh=panel.simh,
        SMB=vp.SMB,
        HML=vp.HML,
        table=vp.table,
    )

    summary = {
        "seed": args.seed,
        "settings": {
            "n_firms": cal.N,
            "periods": cal.Ts,
            "nkp": cal.kp.size,
            "capital_grid": args.capital_grid,
            "nx": cal.nx,
            "nz": cal.nz,
            "cutoff": args.cutoff,
            "ks_iterations": args.ks_iterations,
            "vfi_iterations": args.vfi_iterations,
            "lag": args.lag,
            "matlab_compat": args.matlab_compat,
            "initial_coefficients": None if coefficient_path is None else str(coefficient_path),
        },
        "equilibrium": {
            "coefficients": eq.coefficients.tolist(),
            "r2": float(eq.r2),
            "sighat": float(eq.sighat),
            "error": float(eq.error),
            "iterations": int(eq.iterations),
            "vfi_iterations_final": int(eq.solution.iterations),
            "vfi_errV_final": float(eq.solution.errV),
            "vfi_errK_final": float(eq.solution.errK),
        },
        "panel_diagnostics": {
            "capital_min": float(bf.min()),
            "capital_mean": float(bf.mean()),
            "capital_max": float(bf.max()),
            "average_cross_sectional_capital_std": float(np.mean(np.std(bf, axis=0, ddof=1))),
            "firm_value_min": float(pf.min()),
            "firm_value_max": float(pf.max()),
            "firm_return_min": float(rf.min()),
            "firm_return_max": float(rf.max()),
        },
        "aggregate_ratios": {
            "iyr": panel.iyr,
            "ikr": panel.ikr,
            "theta": panel.theta,
            "dyr": panel.dyr,
            "fyr": panel.fyr,
        },
        "value_premium_table": summarize_table(vp.table),
    }
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("\nSaved outputs")
    print(f"  {npz_path}")
    print(f"  {json_path}")


if __name__ == "__main__":
    main()
