from __future__ import annotations

import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from zhang2005.calibration import calibrate, get_rf_cc, rouwenhorst
from zhang2005.equilibrium import construct_equilibrium
from zhang2005.portfolios import value_premium
from zhang2005.simulation import simulate_panel
from zhang2005.vfi import solve_vfi


class CalibrationSmokeTest(unittest.TestCase):
    def test_benchmark_calibration_shapes(self) -> None:
        cal = calibrate(N=10, Ts=20, nkp=25)
        self.assertEqual(cal.Qx.shape, (cal.nx, cal.nx))
        self.assertEqual(cal.Qz.shape, (cal.nz, cal.nz))
        self.assertTrue(np.allclose(cal.Qx.sum(axis=0), 1.0))
        self.assertTrue(np.allclose(cal.Qz.sum(axis=0), 1.0))
        self.assertTrue(np.any(np.isclose(cal.k, 1.0)))
        self.assertEqual(cal.kp.size, 25)

    def test_calibration_matches_original_params_file(self) -> None:
        params_path = Path("vpCode/Params.mat")
        self.assertTrue(params_path.exists())
        matlab = loadmat(params_path, squeeze_me=True)
        cal = calibrate()

        self.assertAlmostEqual(float(matlab["xbar"]), cal.xbar, places=12)
        self.assertTrue(np.allclose(matlab["Qx"], cal.Qx))
        self.assertTrue(np.allclose(matlab["Qz"], cal.Qz))
        self.assertTrue(np.allclose(matlab["x"], cal.x))
        self.assertTrue(np.allclose(matlab["z"], cal.z))
        self.assertTrue(np.allclose(matlab["k"], cal.k))
        self.assertTrue(np.allclose(np.ravel(matlab["rf"]), cal.rf))


class PortfolioSmokeTest(unittest.TestCase):
    def test_value_premium_synthetic_panel(self) -> None:
        rng = np.random.default_rng(123)
        n, periods = 240, 84
        values = np.exp(rng.normal(size=(n, periods)))
        books = values * np.exp(rng.normal(scale=0.4, size=(n, periods)))
        firm_returns = 1.0 + rng.normal(loc=0.01, scale=0.04, size=(n, periods))
        market = 1.0 + rng.normal(loc=0.008, scale=0.03, size=periods)
        rf = np.full(periods, 1.002)

        result = value_premium(values, firm_returns, books, market, rf, matlab_compat=False)
        self.assertEqual(result.table.shape, (3, 9))
        self.assertEqual(result.SMB.shape, (periods - 60,))
        self.assertTrue(np.isfinite(result.table).all())


class VFISmokeTest(unittest.TestCase):
    def test_tiny_vfi_and_panel_shapes(self) -> None:
        base = calibrate(N=6, Ts=5, nkp=4)
        Qx, x = rouwenhorst(0.2, base.xbar, 0.01, 3)
        Qz, z = rouwenhorst(0.3, 0.0, 0.10, 3)
        cal = replace(
            base,
            N=6,
            Ts=5,
            nx=3,
            nz=3,
            Qx=Qx,
            Qz=Qz,
            x=x,
            z=z,
            rf=get_rf_cc(x, Qx, base.beta, base.gamA, base.gamB),
            k=np.array([0.5, 1.0, 1.5]),
            kp=np.array([0.5, 0.8, 1.1, 1.5]),
            h=np.array([2.8, 3.0, 3.2]),
        )
        sol = solve_vfi(cal, alp1=0.0, alp2=1.0, alp3=0.0, max_iter=1, tol_v=1e12, tol_k=1e12, progress_every=None)
        self.assertEqual(sol.optK.shape, (3, 3, 3, 3))
        self.assertEqual(sol.V.shape, (3, 3, 3, 3))

        rng = np.random.default_rng(456)
        panel = simulate_panel(cal, np.ones(cal.N), np.zeros(cal.N), np.full(cal.Ts, cal.xbar), sol.optK, sol.V, rng)
        self.assertEqual(panel.Pf.shape, (cal.N, cal.Ts))
        self.assertEqual(panel.Rf.shape, (cal.N, cal.Ts - 1))
        self.assertEqual(panel.Rm.shape, (cal.Ts - 1,))

    def test_tiny_equilibrium_loop_shapes(self) -> None:
        base = calibrate(N=6, Ts=8, nkp=4)
        Qx, x = rouwenhorst(0.2, base.xbar, 0.01, 3)
        Qz, z = rouwenhorst(0.3, 0.0, 0.10, 3)
        cal = replace(
            base,
            N=6,
            Ts=8,
            nx=3,
            nz=3,
            Qx=Qx,
            Qz=Qz,
            x=x,
            z=z,
            rf=get_rf_cc(x, Qx, base.beta, base.gamA, base.gamB),
            k=np.array([0.5, 1.0, 1.5]),
            kp=np.array([0.5, 0.8, 1.1, 1.5]),
            h=np.array([2.8, 3.0, 3.2]),
        )
        result = construct_equilibrium(
            cal,
            cutoff=1,
            max_iterations=1,
            coefficient_tol=1e12,
            vfi_options={"max_iter": 1, "tol_v": 1e12, "tol_k": 1e12, "progress_every": None},
            require_convergence=False,
            rng=np.random.default_rng(789),
        )
        self.assertEqual(result.coefficients.shape, (4,))
        self.assertEqual(result.solution.optK.shape, (3, 3, 3, 3))


if __name__ == "__main__":
    unittest.main()
