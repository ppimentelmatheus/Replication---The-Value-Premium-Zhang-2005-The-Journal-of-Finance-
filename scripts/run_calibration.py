#!/usr/bin/env python3
"""Print Zhang (2005) benchmark calibration objects."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zhang2005.calibration import annualized_average_sharpe, calibrate


def main() -> None:
    cal = calibrate()
    print(f"xbar: {cal.xbar:.10f}")
    print(f"capital grid points: {cal.k.size}")
    print(f"k min/max: {cal.k.min():.4f} / {cal.k.max():.4f}")
    print(f"kp grid points: {cal.kp.size}")
    print(f"h grid: {cal.h}")
    print(f"annualized average Sharpe ratio: {annualized_average_sharpe(cal.gamA, cal.stdx):.6f}")


if __name__ == "__main__":
    main()
