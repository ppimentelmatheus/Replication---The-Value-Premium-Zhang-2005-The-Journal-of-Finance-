#!/usr/bin/env python3
"""Generate model-only tables for Zhang (2005)."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zhang2005.calibration import calibrate
from zhang2005.model_tables import (
    aggregate_moment_rows,
    paper_table2_rows,
    parameter_rows,
    portfolio_summary_rows,
    predictive_regression_rows,
    sort_by_book_to_market,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-output-dir",
        type=Path,
        default=None,
        help="directory produced by scripts/run_full_replication.py; preferred for full-spec tables",
    )
    parser.add_argument(
        "--replication-npz",
        type=Path,
        default=Path("outputs/replication_for_figures/mini_replication_results.npz"),
        help="NPZ produced by scripts/run_mini_replication.py",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=Path("outputs/replication_for_figures/mini_replication_summary.json"),
        help="JSON summary produced by scripts/run_mini_replication.py",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/tables/model"))
    return parser.parse_args()


def _format_value(value: Any) -> str:
    if isinstance(value, (float, np.floating)):
        if np.isnan(value):
            return ""
        return f"{float(value):.6g}"
    return str(value)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, title: str, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    lines = [f"# {title}", ""]
    lines.append("| " + " | ".join(columns) + " |")
    lines.append("| " + " | ".join(["---"] * len(columns)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(_format_value(row.get(col, "")) for col in columns) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_table(base: Path, stem: str, title: str, rows: list[dict[str, Any]]) -> None:
    write_csv(base / "csv" / f"{stem}.csv", rows)
    write_markdown(base / "markdown" / f"{stem}.md", title, rows)
    print(f"wrote {stem}.csv/md")


def load_first_panel_arrays(full_output_dir: Path) -> tuple[dict[str, np.ndarray] | None, Path | None]:
    candidates = sorted((full_output_dir / "panels").glob("panel_*_arrays.npz"))
    if not candidates:
        return None, None
    return dict(np.load(candidates[0])), candidates[0]


def load_factor_series(full_output_dir: Path) -> dict[str, np.ndarray]:
    files = sorted((full_output_dir / "panels").glob("panel_*_factors.npz"))
    if not files:
        return {}
    series_names = ["Rm", "srf", "HML", "SMB", "GDPg"]
    out: dict[str, list[np.ndarray]] = {name: [] for name in series_names}
    tables = []
    ratios = []
    for path in files:
        data = np.load(path)
        for name in series_names:
            if name in data:
                out[name].append(np.asarray(data[name]).reshape(-1))
        if "table" in data:
            tables.append(np.asarray(data["table"]))
        if "ratios" in data:
            ratios.append(np.asarray(data["ratios"]).reshape(-1))
    result = {name: np.concatenate(values) for name, values in out.items() if values}
    if tables:
        result["value_premium_tables"] = np.asarray(tables)
        result["mean_value_premium_table"] = np.nanmean(np.asarray(tables), axis=0)
    if ratios:
        result["aggregate_ratios"] = np.asarray(ratios)
        result["mean_aggregate_ratios"] = np.nanmean(np.asarray(ratios), axis=0)
    return result


def value_premium_table_rows(table: np.ndarray) -> list[dict[str, Any]]:
    columns = ["Mkt-Rf", "SMB", "HML", "SL", "SM", "SH", "BL", "BM", "BH"]
    row_names = ["percent_mean", "percent_std", "t_stat"]
    return [
        {"statistic": row_name, **{col: float(table[i, j]) for j, col in enumerate(columns)}}
        for i, row_name in enumerate(row_names)
    ]


def aggregate_ratio_rows(ratios: np.ndarray) -> list[dict[str, Any]]:
    names = ["investment_output", "investment_capital", "adjustment_cost_share", "dividend_output", "fixed_cost_output"]
    ratios = np.asarray(ratios, dtype=float).reshape(-1)
    return [{"ratio": name, "value": float(value)} for name, value in zip(names, ratios, strict=True)]


def load_full_output(full_output_dir: Path) -> tuple[dict[str, Any], dict[str, np.ndarray], dict[str, np.ndarray] | None, Path | None]:
    if not full_output_dir.exists():
        raise FileNotFoundError(f"full output directory not found: {full_output_dir}")
    summary_path = full_output_dir / "full_replication_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else {}
    factors = load_factor_series(full_output_dir)
    panel_arrays, panel_path = load_first_panel_arrays(full_output_dir)
    if panel_path is not None:
        print(f"using full panel arrays for portfolio sorts/regressions: {panel_path}")
    else:
        print("no full panel arrays found; portfolio sorts/regressions that need Pf/Bf/Rf will be skipped")
    return summary, factors, panel_arrays, panel_path


def main() -> None:
    args = parse_args()
    if args.full_output_dir is not None:
        summary, factors, panel_arrays, _ = load_full_output(args.full_output_dir)
        config = summary.get("config", {})
        n_firms = int(config.get("n_firms", panel_arrays["Pf"].shape[0] if panel_arrays else 5000))
        periods = int(config.get("equilibrium_periods", 11000))
        nkp = int(config.get("nkp", 5000))
        cal = calibrate(N=n_firms, Ts=periods, nkp=nkp)
        if panel_arrays is None:
            data = factors
        else:
            data = {**panel_arrays, **factors}
    else:
        if not args.replication_npz.exists():
            raise FileNotFoundError(f"Run scripts/run_mini_replication.py first: {args.replication_npz}")
        data = dict(np.load(args.replication_npz))
        summary = json.loads(args.summary_json.read_text(encoding="utf-8")) if args.summary_json.exists() else {}
        settings = summary.get("settings", {})
        cal = calibrate(
            N=int(settings.get("n_firms", data["Pf"].shape[0])),
            Ts=int(settings.get("periods", data["Pf"].shape[1] + 1)),
            nkp=int(settings.get("nkp", 5000)),
        )

    save_table(args.output_dir, "table1_model_parameters", "Table I - Model Parameters", parameter_rows(cal))
    if {"Pf", "Bf", "Rf", "Rm", "srf"}.issubset(data):
        sort = sort_by_book_to_market(data["Pf"], data["Rf"], data["Bf"])
        save_table(args.output_dir, "table2_model_moments", "Table II - Key Moments, Model Only", paper_table2_rows(cal, data))
        save_table(args.output_dir, "table2_additional_model_moments", "Additional Model Moments", aggregate_moment_rows(data))
        save_table(
            args.output_dir,
            "table3_bm_portfolios_model",
            "Table III - Book-to-Market Portfolios, Model Only",
            portfolio_summary_rows(sort, data["Rm"], data["srf"]),
        )
        save_table(
            args.output_dir,
            "table5_6_predictive_regressions_model",
            "Tables V-VI - Predictive Regressions, Model Only",
            predictive_regression_rows(data, sort, cal.xbar),
        )
    else:
        print("skipping moments/portfolio/predictive tables: Pf/Bf/Rf/Rm/srf not all available")

    if "mean_value_premium_table" in data:
        save_table(
            args.output_dir,
            "table_value_premium_full_model",
            "Full Simulation - Mean Value Premium Table",
            value_premium_table_rows(data["mean_value_premium_table"]),
        )
    if "mean_aggregate_ratios" in data:
        save_table(
            args.output_dir,
            "table_aggregate_ratios_full_model",
            "Full Simulation - Mean Aggregate Ratios",
            aggregate_ratio_rows(data["mean_aggregate_ratios"]),
        )


if __name__ == "__main__":
    main()
