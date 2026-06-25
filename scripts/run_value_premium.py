#!/usr/bin/env python3
"""Compute SMB/HML from simulated panel arrays saved in a MATLAB file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zhang2005.io import load_mat
from zhang2005.portfolios import value_premium


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mat_file", help="MAT file containing Pf, Rf, Bf, Rm and srf/rf arrays")
    parser.add_argument("--matlab-compat", action="store_true", help="use Zhang's original ValPrem.m timing")
    args = parser.parse_args()

    data = load_mat(args.mat_file)
    rate_name = "srf" if "srf" in data else "rf"
    result = value_premium(data["Pf"], data["Rf"], data["Bf"], data["Rm"], data[rate_name], matlab_compat=args.matlab_compat)
    print(result.table)


if __name__ == "__main__":
    main()
