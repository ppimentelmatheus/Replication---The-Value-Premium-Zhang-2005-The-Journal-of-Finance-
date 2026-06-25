#!/usr/bin/env python3
"""Plot Zhang (2005) article figures from precomputed figure data."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-zhang2005")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("outputs/figures/data"))
    parser.add_argument("--png-dir", type=Path, default=Path("outputs/figures/png"))
    parser.add_argument("--pdf-dir", type=Path, default=Path("outputs/figures/pdf"))
    return parser.parse_args()


def save(fig: plt.Figure, name: str, png_dir: Path, pdf_dir: Path) -> None:
    png_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(png_dir / f"{name}.png", dpi=180)
    fig.savefig(pdf_dir / f"{name}.pdf")
    plt.close(fig)
    print(f"wrote {name}.png/pdf")


def maybe_load(path: Path) -> dict[str, np.ndarray] | None:
    if not path.exists():
        print(f"missing {path.name}, skipping")
        return None
    return dict(np.load(path))


def plot_figure1(data: dict[str, np.ndarray], png_dir: Path, pdf_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(data["investment_rate"], data["adjustment_cost"], color="black", linewidth=2)
    ax.axvline(0.0, color="0.65", linewidth=1)
    ax.set_xlabel("Investment rate, i/k")
    ax.set_ylabel("Adjustment cost")
    ax.set_title("Figure 1. Asymmetric Adjustment Cost")
    save(fig, "figure1_adjustment_cost", png_dir, pdf_dir)


def plot_figure5(data: dict[str, np.ndarray], png_dir: Path, pdf_dir: Path) -> None:
    x = data["x"]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    panels = [
        ("benchmark_sharpe_annual", "constant_sharpe_annual", "Annualized Sharpe ratio"),
        ("benchmark_sigm", "constant_sigm", r"$\sigma_t[M_{t+1}]$"),
        ("benchmark_em", "constant_em", r"$E_t[M_{t+1}]$"),
    ]
    for ax, (bench_key, const_key, ylabel) in zip(axes, panels, strict=True):
        ax.plot(x, data[bench_key], color="black", linewidth=2, label=r"$\gamma_1=-1000$")
        ax.plot(x, data[const_key], color="0.35", linewidth=2, linestyle="--", label=r"$\gamma_1=0$")
        ax.set_xlabel("Aggregate productivity, x")
        ax.set_ylabel(ylabel)
        ax.grid(True, color="0.9")
    axes[0].legend(frameon=False)
    fig.suptitle("Figure 5. Properties of the Pricing Kernel")
    save(fig, "figure5_pricing_kernel", png_dir, pdf_dir)


def plot_figure2(data: dict[str, np.ndarray], png_dir: Path, pdf_dir: Path) -> None:
    years = np.arange(-5, 6)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    axes[0].plot(years, data["roevec"][:, 0], color="black", linewidth=2, label="Growth")
    axes[0].plot(years, data["roevec"][:, 1], color="black", linewidth=2, linestyle="--", label="Value")
    axes[0].set_xlabel("Formation year")
    axes[0].set_ylabel("Profitability")
    axes[0].set_title("Panel A: Return on Equity")
    axes[0].legend(frameon=False)
    t = np.arange(data["ROElow"].size)
    axes[1].plot(t, data["ROElow"], color="black", linewidth=2, label="Growth")
    axes[1].plot(t, data["ROEhigh"], color="black", linewidth=2, linestyle="--", label="Value")
    axes[1].set_xlabel("Year")
    axes[1].set_ylabel("Profitability")
    axes[1].set_title("Panel B: Time-Series of ROE")
    axes[1].legend(frameon=False)
    fig.suptitle("Figure 2. The Value Factor in Profitability")
    save(fig, "figure2_ff95_profitability", png_dir, pdf_dir)


def plot_figure3(data: dict[str, np.ndarray], png_dir: Path, pdf_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    panels = [
        ("Panel A: Bad Times", "xg_bust", "yg_bust", "xv_bust", "yv_bust"),
        ("Panel B: Good Times", "xg_boom", "yg_boom", "xv_boom", "yv_boom"),
    ]
    for ax, (title, xg, yg, xv, yv) in zip(axes, panels, strict=True):
        ax.scatter(data[xg], data[yg], s=18, marker="o", facecolors="none", edgecolors="black", label="Growth")
        ax.scatter(data[xv], data[yv], s=24, marker="+", color="black", label="Value")
        ax.set_xlabel("i/k")
        ax.set_title(title)
        ax.grid(True, color="0.92")
    axes[0].set_ylabel("Adjustment cost")
    axes[0].legend(frameon=False)
    fig.suptitle("Figure 3. The Value Factor in Corporate Investment")
    save(fig, "figure3_investment_scatter", png_dir, pdf_dir)


def plot_figure4(data: dict[str, np.ndarray], png_dir: Path, pdf_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    axes[0].scatter(data["x_for_hml"], data["hml"], s=16, color="black", alpha=0.75)
    axes[0].set_xlabel("Aggregate productivity, x")
    axes[0].set_ylabel("HML realized return")
    axes[0].set_title("Panel A: Simulated HML")
    axes[1].scatter(data["x_for_value_spread"], data["value_spread"], s=16, color="black", alpha=0.75)
    axes[1].set_xlabel("Aggregate productivity, x")
    axes[1].set_ylabel("Log B/M spread")
    axes[1].set_title("Panel B: Value Spread")
    fig.suptitle("Figure 4. Simulated Spreads")
    save(fig, "figure4_simulated_spreads", png_dir, pdf_dir)


def plot_figure_b1(data: dict[str, np.ndarray], png_dir: Path, pdf_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    axes[0].scatter(data["predicted_scaled"], data["actual_scaled"], s=14, color="black", alpha=0.7)
    lo = min(float(data["predicted_scaled"].min()), float(data["actual_scaled"].min()))
    hi = max(float(data["predicted_scaled"].max()), float(data["actual_scaled"].max()))
    axes[0].plot([lo, hi], [lo, hi], color="0.45", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Predicted output price")
    axes[0].set_ylabel("Actual output price")
    axes[0].set_title("Panel A: Predicted versus Actual")
    axes[1].hist(data["residual_pct"], bins=40, color="0.25")
    axes[1].set_xlabel("100 x excess demand / output")
    axes[1].set_ylabel("Frequency")
    axes[1].set_title("Panel B: Excess Demand")
    fig.suptitle("Figure B.1. Quality of Aggregation")
    save(fig, "figureB1_aggregation_quality", png_dir, pdf_dir)


def main() -> None:
    args = parse_args()
    plotters = [
        ("figure1_adjustment_cost.npz", plot_figure1),
        ("figure5_pricing_kernel.npz", plot_figure5),
        ("figure2_ff95_profitability.npz", plot_figure2),
        ("figure3_investment_scatter.npz", plot_figure3),
        ("figure4_simulated_spreads.npz", plot_figure4),
        ("figureB1_aggregation_quality.npz", plot_figure_b1),
    ]
    for filename, plotter in plotters:
        data = maybe_load(args.data_dir / filename)
        if data is not None:
            plotter(data, args.png_dir, args.pdf_dir)


if __name__ == "__main__":
    main()
