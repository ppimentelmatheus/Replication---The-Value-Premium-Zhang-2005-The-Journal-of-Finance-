# Zhang (2005) Replication Porting Notes

This repository now keeps the original MATLAB/Fortran files in `vpCode/` and
starts a Python port in `zhang2005/`.

## Original Pipeline

1. `CalibCC.m`: calibrates the monthly model, builds Rouwenhorst grids for
   aggregate and idiosyncratic shocks, creates the capital and price grids, and
   saves `Params.mat`.
2. `mainCC.m`: solves the firm dynamic program and the industry equilibrium
   using Krusell and Smith approximate aggregation. The original speed-critical
   routines are `vfi3fcnIEccB.f90` and `simIEfcn3.f90`.
3. `ssIE.m`: simulates firms until the stationary cross-sectional distribution,
   using `ssIEfcn.f90`.
4. `panIE.m`: simulates panel data, computes Fama-French tests, and constructs
   value and size premia using `panIEfcn.f90`, `FF92.m`, `FF93.m`, `FF95.m`,
   and `ValPrem.m`.

## Python Modules

- `zhang2005.calibration`: port of `CalibCC.m`, `rouwTrans.m`, `CspSimu.m`,
  and `getRfcc.m`.
- `zhang2005.vfi`: direct NumPy translation of `vfi3fcnIEccB.f90`.
- `zhang2005.simulation`: translations of `simIEfcn3.f90`, `ssIEfcn.f90`, and
  `panIEfcn.f90`.
- `zhang2005.equilibrium`: port of the `mainCC.m` Krusell-Smith loop.
- `zhang2005.portfolios`: port of `ValPrem.m` for SMB/HML and summary moments.
- `zhang2005.linalg`: port of `ols.m`.
- `zhang2005.io`: MATLAB `.mat` interoperability.

## Notes

The benchmark model is large: 5,000 firms, 11,000 aggregate periods, 5,000
candidate next-period capital choices, and a four-dimensional policy grid. The
Python code is written to match the original logic first. For production-size
runs, the next step is to add Numba acceleration to `vfi.py` and
`simulation.py`, or to wrap the original Fortran with `f2py`.

`ValPrem.m` appears to use portfolio return windows indexed as `t-lag:t-lag+11`
after sorting at month `t`; this has been preserved with
`matlab_compat=True`. Passing `matlab_compat=False` uses the rebalance month
itself as the first return month, which is closer to the timing in `FF93.m`.
