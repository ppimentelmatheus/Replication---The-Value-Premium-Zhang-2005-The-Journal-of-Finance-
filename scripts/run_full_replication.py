#!/usr/bin/env python3
"""Run the full Zhang (2005) model-simulation pipeline.

The defaults mirror the benchmark model in the original MATLAB code:

* equilibrium: N=5000, Ts=11000, cutoff=1000, nkp=5000
* stationary distribution: 10000 months
* panel simulation: 20 panels with 421 monthly observations

This script is intentionally explicit about phases because the full benchmark is
large. Use ``--dry-run`` first.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zhang2005.calibration import Calibration, calibrate, simulate_ar1
from zhang2005.equilibrium import construct_equilibrium
from zhang2005.interpolation import interp1_extrap
from zhang2005.io import load_mat, save_npz
from zhang2005.portfolios import value_premium
from zhang2005.simulation import simulate_panel, simulate_stationary_distribution


@dataclass(frozen=True)
class FullRunConfig:
    n_firms: int
    equilibrium_periods: int
    stationary_periods: int
    panel_periods: int
    panel_simulations: int
    nkp: int
    cutoff: int
    ks_max_iterations: int
    coefficient_tol: float
    vfi_max_iterations: int
    vfi_tol_v: float
    vfi_tol_k: float
    seed: int
    use_matlab_equilibrium: bool
    matlab_solution: str
    matlab_coefficients: str
    stationary_mat: str | None
    save_panel_arrays: str
    output_dir: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-firms", type=int, default=5000)
    parser.add_argument("--equilibrium-periods", type=int, default=11000)
    parser.add_argument("--stationary-periods", type=int, default=10000)
    parser.add_argument("--panel-periods", type=int, default=421, help="421 for panIE.m; 901 for panIE900.m")
    parser.add_argument("--panel-simulations", type=int, default=20, help="20 for panIE.m; 10 for panIE900.m")
    parser.add_argument("--nkp", type=int, default=5000)
    parser.add_argument("--cutoff", type=int, default=1000)
    parser.add_argument("--ks-max-iterations", type=int, default=100)
    parser.add_argument("--coefficient-tol", type=float, default=1e-2)
    parser.add_argument("--vfi-max-iterations", type=int, default=10_000)
    parser.add_argument("--vfi-tol-v", type=float, default=1e-3)
    parser.add_argument("--vfi-tol-k", type=float, default=1e-8)
    parser.add_argument("--vfi-progress-every", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20250625)
    parser.add_argument(
        "--use-matlab-equilibrium",
        action="store_true",
        help="skip Python equilibrium solve and load vpCode/vfi3Mat.mat plus vpCode/coefIS.mat",
    )
    parser.add_argument("--matlab-solution", type=Path, default=Path("vpCode/vfi3Mat.mat"))
    parser.add_argument("--matlab-coefficients", type=Path, default=Path("vpCode/coefIS.mat"))
    parser.add_argument(
        "--stationary-mat",
        type=Path,
        default=None,
        help="optional MAT file with kd0/zd0, e.g. vpCode/distrSS.mat; skips stationary simulation",
    )
    parser.add_argument(
        "--save-panel-arrays",
        choices=("none", "first", "all"),
        default="first",
        help="full panel arrays are large; compact factor files are always saved",
    )
    parser.add_argument("--price-floor", type=float, default=1e-8)
    parser.add_argument("--ff93-timing", action="store_true", help="use rebalance-month timing instead of original ValPrem.m timing")
    parser.add_argument("--allow-nonconvergence", action="store_true", help="return last equilibrium iteration if coefficients do not converge")
    parser.add_argument("--dry-run", action="store_true", help="print configuration and memory estimates without running")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/full_replication"))
    return parser.parse_args()


def config_from_args(args: argparse.Namespace) -> FullRunConfig:
    return FullRunConfig(
        n_firms=args.n_firms,
        equilibrium_periods=args.equilibrium_periods,
        stationary_periods=args.stationary_periods,
        panel_periods=args.panel_periods,
        panel_simulations=args.panel_simulations,
        nkp=args.nkp,
        cutoff=args.cutoff,
        ks_max_iterations=args.ks_max_iterations,
        coefficient_tol=args.coefficient_tol,
        vfi_max_iterations=args.vfi_max_iterations,
        vfi_tol_v=args.vfi_tol_v,
        vfi_tol_k=args.vfi_tol_k,
        seed=args.seed,
        use_matlab_equilibrium=args.use_matlab_equilibrium,
        matlab_solution=str(args.matlab_solution),
        matlab_coefficients=str(args.matlab_coefficients),
        stationary_mat=None if args.stationary_mat is None else str(args.stationary_mat),
        save_panel_arrays=args.save_panel_arrays,
        output_dir=str(args.output_dir),
    )


def format_gb(bytes_count: float) -> str:
    return f"{bytes_count / 1024**3:.2f} GB"


def print_dry_run(cal: Calibration, args: argparse.Namespace) -> None:
    shape = (cal.k.size, cal.h.size, cal.x.size, cal.z.size)
    panel_cols = args.panel_periods - 1
    one_panel_array = args.n_firms * panel_cols * 8
    saved_panel_arrays = 5 * one_panel_array
    all_saved_panel_arrays = args.panel_simulations * saved_panel_arrays
    print("Full replication dry run")
    print(f"  calibration: N={cal.N}, equilibrium Ts={cal.Ts}, nkp={cal.kp.size}")
    print(f"  grids: nk={cal.k.size}, nh={cal.h.size}, nx={cal.x.size}, nz={cal.z.size}")
    print(f"  policy/value shape: {shape}")
    print(f"  stationary periods: {args.stationary_periods}")
    print(f"  panel periods x simulations: {args.panel_periods} x {args.panel_simulations}")
    print(f"  one N x (T-1) panel array: {format_gb(one_panel_array)}")
    print(f"  Pf/Bf/Df/In/Rf for one panel: {format_gb(saved_panel_arrays)}")
    print(f"  Pf/Bf/Df/In/Rf for all panels: {format_gb(all_saved_panel_arrays)}")
    print(f"  save_panel_arrays: {args.save_panel_arrays}")
    print(f"  use_matlab_equilibrium: {args.use_matlab_equilibrium}")
    print()
    print("Recommendation:")
    print("  1. First run with --dry-run.")
    print("  2. Then try --use-matlab-equilibrium to test the full panel simulator.")
    print("  3. Run from-scratch equilibrium only after accelerating VFI/simulation.")


def reshape_matlab_state(array: np.ndarray, cal: Calibration) -> np.ndarray:
    shape = (cal.k.size, cal.h.size, cal.x.size, cal.z.size)
    array = np.asarray(array, dtype=float)
    if array.shape == shape:
        return array
    if array.size != np.prod(shape):
        raise ValueError(f"MATLAB array has size {array.size}, expected {np.prod(shape)} for shape {shape}")
    return array.reshape(shape, order="F")


def load_initial_coefficients(path: Path) -> tuple[float, float, float, float]:
    data = load_mat(path)
    return tuple(float(data[name]) for name in ("alp1", "alp2", "alp3", "alp4"))


def load_matlab_equilibrium(cal: Calibration, solution_path: Path, coefficient_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    solution = load_mat(solution_path)
    coefficients = np.asarray(load_initial_coefficients(coefficient_path), dtype=float)
    optK = reshape_matlab_state(solution["optK"], cal)
    value_name = "V0" if "V0" in solution else "V"
    V = reshape_matlab_state(solution[value_name], cal)
    return coefficients, optK, V


def load_stationary_distribution(path: Path) -> tuple[np.ndarray, np.ndarray, float | None]:
    data = load_mat(path)
    kd_name = "kd0" if "kd0" in data else "kd"
    zd_name = "zd0" if "zd0" in data else "zd"
    xinl = float(data["xinl"]) if "xinl" in data else None
    return np.asarray(data[kd_name], dtype=float).reshape(-1), np.asarray(data[zd_name], dtype=float).reshape(-1), xinl


def save_equilibrium_outputs(output_dir: Path, coefficients: np.ndarray, optK: np.ndarray, V: np.ndarray, metadata: dict) -> None:
    save_npz(output_dir / "equilibrium_solution.npz", coefficients=coefficients, optK=optK, V=V)
    (output_dir / "equilibrium_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def panel_should_save(mode: str, index: int) -> bool:
    return mode == "all" or (mode == "first" and index == 1)


def run_panels(
    cal: Calibration,
    optK: np.ndarray,
    V: np.ndarray,
    kd0: np.ndarray,
    zd0: np.ndarray,
    args: argparse.Namespace,
    rng: np.random.Generator,
) -> dict:
    panels_dir = args.output_dir / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)

    tables = []
    ratios = []
    factor_files = []
    x_initial = cal.xbar
    matlab_compat = not args.ff93_timing

    kd = np.asarray(kd0, dtype=float).copy()
    zd = np.asarray(zd0, dtype=float).copy()
    for sim in range(1, args.panel_simulations + 1):
        print(f"\nPanel simulation {sim}/{args.panel_simulations}")
        sx, x_initial = simulate_ar1(x_initial, cal.xbar, cal.rhox, cal.stdx, args.panel_periods, rng)
        srf = interp1_extrap(cal.x, cal.rf, sx).reshape(-1)[:-1]
        panel = simulate_panel(cal, kd, zd, sx, optK, V, rng)
        kd, zd = panel.kd, panel.zd

        Pf = np.maximum(panel.Pf[:, :-1], args.price_floor)
        Bf = panel.Bf[:, :-1]
        Df = panel.Df[:, :-1]
        In = panel.In[:, :-1]
        Rf = panel.Rf
        Rm = panel.Rm

        vp = value_premium(Pf, Rf, Bf, Rm, srf, matlab_compat=matlab_compat)
        tables.append(vp.table)
        ratios.append([panel.iyr, panel.ikr, panel.theta, panel.dyr, panel.fyr])

        compact_path = panels_dir / f"panel_{sim:03d}_factors.npz"
        save_npz(
            compact_path,
            table=vp.table,
            SMB=vp.SMB,
            HML=vp.HML,
            Rm=Rm,
            srf=srf,
            sx=sx,
            GDPg=panel.GDPg,
            simh=panel.simh,
            kd_final=kd,
            zd_final=zd,
            ratios=np.asarray(ratios[-1]),
        )
        factor_files.append(str(compact_path))

        if panel_should_save(args.save_panel_arrays, sim):
            panel_path = panels_dir / f"panel_{sim:03d}_arrays.npz"
            save_npz(panel_path, Pf=Pf, Bf=Bf, Df=Df, In=In, Rf=Rf, Rm=Rm, srf=srf, sx=sx)
            print(f"  saved full panel arrays: {panel_path}")

        print(f"  value premium table mean row: {np.array2string(vp.table[0], precision=4)}")

    tables_array = np.asarray(tables)
    ratios_array = np.asarray(ratios)
    save_npz(
        args.output_dir / "panel_summary_arrays.npz",
        value_premium_tables=tables_array,
        aggregate_ratios=ratios_array,
        mean_value_premium_table=np.nanmean(tables_array, axis=0),
        mean_aggregate_ratios=np.nanmean(ratios_array, axis=0),
        kd_final=kd,
        zd_final=zd,
    )
    return {
        "factor_files": factor_files,
        "mean_value_premium_table": np.nanmean(tables_array, axis=0).tolist(),
        "mean_aggregate_ratios": {
            "iyr": float(np.nanmean(ratios_array[:, 0])),
            "ikr": float(np.nanmean(ratios_array[:, 1])),
            "theta": float(np.nanmean(ratios_array[:, 2])),
            "dyr": float(np.nanmean(ratios_array[:, 3])),
            "fyr": float(np.nanmean(ratios_array[:, 4])),
        },
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    cal = calibrate(N=args.n_firms, Ts=args.equilibrium_periods, nkp=args.nkp)

    if args.dry_run:
        print_dry_run(cal, args)
        (args.output_dir / "full_replication_config.json").write_text(
            json.dumps(asdict(config_from_args(args)), indent=2) + "\n",
            encoding="utf-8",
        )
        return

    print("Full replication setup")
    print(json.dumps(asdict(config_from_args(args)), indent=2))

    if args.use_matlab_equilibrium:
        print("\nLoading MATLAB equilibrium solution...")
        coefficients, optK, V = load_matlab_equilibrium(cal, args.matlab_solution, args.matlab_coefficients)
        equilibrium_metadata = {
            "source": "matlab",
            "coefficients": coefficients.tolist(),
            "solution_path": str(args.matlab_solution),
            "coefficient_path": str(args.matlab_coefficients),
        }
    else:
        print("\nSolving equilibrium in Python...")
        initial_coefficients = load_initial_coefficients(args.matlab_coefficients) if args.matlab_coefficients.exists() else (0.0, 1.0, 0.0, 0.0)
        initial_value = None
        if args.matlab_solution.exists():
            initial_value = reshape_matlab_state(load_mat(args.matlab_solution)["V0"], cal)
            print(f"  using initial V0 from {args.matlab_solution}")
        eq = construct_equilibrium(
            cal,
            initial_coefficients=initial_coefficients,
            initial_value_function=initial_value,
            cutoff=args.cutoff,
            coefficient_tol=args.coefficient_tol,
            max_iterations=args.ks_max_iterations,
            vfi_options={
                "max_iter": args.vfi_max_iterations,
                "tol_v": args.vfi_tol_v,
                "tol_k": args.vfi_tol_k,
                "progress_every": args.vfi_progress_every,
                "require_convergence": not args.allow_nonconvergence,
            },
            require_convergence=not args.allow_nonconvergence,
            rng=rng,
        )
        coefficients, optK, V = eq.coefficients, eq.solution.optK, eq.solution.V
        equilibrium_metadata = {
            "source": "python",
            "coefficients": coefficients.tolist(),
            "r2": eq.r2,
            "sighat": eq.sighat,
            "tstats": eq.tstats.tolist(),
            "iterations": eq.iterations,
            "error": eq.error,
            "vfi_iterations": eq.solution.iterations,
            "vfi_errV": eq.solution.errV,
            "vfi_errK": eq.solution.errK,
        }
        save_npz(
            args.output_dir / "equilibrium_diagnostics.npz",
            coefficients=coefficients,
            eq_sx=eq.sx,
            eq_simh=eq.simulation.simh,
            eq_sigk=eq.simulation.sigk,
            eq_sigz=eq.simulation.sigz,
            cutoff=np.asarray(args.cutoff),
            residuals=eq.residuals,
        )

    save_equilibrium_outputs(args.output_dir, coefficients, optK, V, equilibrium_metadata)

    if args.stationary_mat is not None:
        print("\nLoading stationary distribution...")
        kd0, zd0, _ = load_stationary_distribution(args.stationary_mat)
    else:
        print("\nSimulating stationary distribution...")
        stationary = simulate_stationary_distribution(cal, optK, periods=args.stationary_periods, rng=rng)
        kd0, zd0 = stationary.kd, stationary.zd
        save_npz(
            args.output_dir / "stationary_distribution.npz",
            kd0=kd0,
            zd0=zd0,
            simh=stationary.simh,
            sigk=stationary.sigk,
            sigz=stationary.sigz,
        )

    print("\nSimulating panels...")
    panel_summary = run_panels(cal, optK, V, kd0, zd0, args, rng)

    summary = {
        "config": asdict(config_from_args(args)),
        "equilibrium": equilibrium_metadata,
        "panels": panel_summary,
    }
    summary_path = args.output_dir / "full_replication_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\nSaved summary: {summary_path}")


if __name__ == "__main__":
    main()
