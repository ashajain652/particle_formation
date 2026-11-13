"""The shared similarity boundary-layer solver (Lees local similarity, total-enthalpy energy equation, real frozen-air
properties), ported from analysis/gas_loads/reference/similarity_heat.py and checked against that reference's committed
outputs (fixtures in tests/fixtures/gas_loads/) and against the Falkner-Skan and cold-wall stagnation values."""
import csv
import math
import os

import numpy as np
import pytest

from reentry_model.shared import gas, similarity

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "gas_loads")
R_AIR, T_INF, V, RADIUS, T_W = 287.05, 200.0, 7108.0, 0.1524, 1500.0
CASES = [("0.002", 1.007e-4, 151.0), ("0.01", 2.014e-5, 151.0), ("0.05", 3.987e-6, 151.0), ("0.25", 8.057e-7, 158.0)]


@pytest.fixture(scope="module")
def air():
    return gas.UnreactedAir()


@pytest.fixture(scope="module")
def props(air):
    return similarity.FrozenAirProperties(air)


def _rows(name):
    with open(os.path.join(FIX, name)) as fh:
        return list(csv.DictReader(fh))


def test_falkner_skan():
    solver = similarity.SimilaritySolver(similarity.ConstantProperties())
    for beta, expect in ((0.5, 0.9277), (0.0, 0.4696), (1.0, 1.2326)):
        solver.reset()
        pt = solver.solve(beta, 1.0e6, 0.0, 1.0e6 / 1005.0)
        assert pt.fpp0 == pytest.approx(expect, abs=5e-4), beta
        assert pt.g_w == pytest.approx(1.0)


def test_cold_wall_stagnation(props):
    solver = similarity.SimilaritySolver(props)
    pt = solver.solve(0.5, 25162715.7, 0.0, 1500.0)
    assert pt.fpp0 == pytest.approx(0.5024, abs=2e-3)
    assert pt.gp0 == pytest.approx(0.3313, abs=2e-3)
    solver.reset()
    assert solver.solve(0.5, 25162715.7, 0.0, 16395.0).fpp0 == pytest.approx(0.9273, abs=1e-3)


def test_march_reproduces_reference(air, props):
    heat, lee = _rows("similarity_heat.csv"), _rows("similarity_lee.csv")
    th = np.arange(0.0, 181.0)
    mu_w = float(props.mu_of_T(T_W))
    for kn, rho, sep in CASES:
        p_inf = rho * R_AIR * T_INF
        p_w = np.array([float(r["p_over_pinf"]) for r in lee if r["case_kn"] == kn]) * p_inf
        stag = air.stagnation(rho, T_INF, V)
        u_e = np.array([math.sqrt(max(0.0, 2.0 * (stag.h - air.expand(stag, float(p)).h))) for p in p_w])
        u_e[0] = 0.0
        sol = similarity.march(th, RADIUS, p_w, u_e, stag.h, T_W, p_w / (R_AIR * T_W) * mu_w,
                               similarity.SimilaritySolver(props))
        assert sol.separation_deg == sep
        rows = [r for r in heat if r["case_kn"] == kn]
        for i, r in enumerate(rows):
            if r["Cq_continuum"] == "":
                assert sol.q_c[i] == 0.0 and sol.tau_c[i] == 0.0 or i == 180
                continue
            assert sol.q_c[i] / (0.5 * rho * V**3) == pytest.approx(float(r["Cq_continuum"]), rel=1e-6)
            assert sol.tau_c[i] / (0.5 * rho * V**2) == pytest.approx(float(r["Ctau_continuum"]), rel=1e-6, abs=1e-12)


def test_failed_solve_is_separation(props, monkeypatch):
    th = np.arange(0.0, 181.0)
    p_w = 1000.0 * (1.0 + 0.1 * np.cos(np.radians(th)))
    u_e = 3000.0 * np.sin(np.radians(th)) + 1.0
    u_e[0] = 0.0
    solver = similarity.SimilaritySolver(props)
    real = similarity.SimilaritySolver.solve
    args = (th, RADIUS, p_w, u_e, 2.0e7, T_W, p_w / (R_AIR * T_W) * 5e-5)
    ref = similarity.march(*args, solver)
    assert solver._sol is not None
    # now force failure from the bin at 100 degrees on
    n = {"i": 0}

    def failing(self, beta, H_e, u, T_w):
        n["i"] += 1
        if n["i"] > 1 + 99:          # call 1 is the stagnation point, then one call per bin; bin 100 is call 101
            raise similarity.SimilarityError("forced")
        return real(self, beta, H_e, u, T_w)
    monkeypatch.setattr(similarity.SimilaritySolver, "solve", failing)
    solver.reset()
    sol = similarity.march(*args, solver)
    assert sol.separation_deg == 100.0
    assert np.all(sol.tau_c[100:] == 0.0) and np.all(sol.q_c[100:] == 0.0)
    assert np.all(sol.q_c[:100] > 0.0)
    # reset() sends the next solve back through the homotopy
    monkeypatch.setattr(similarity.SimilaritySolver, "solve", real)
    solver.reset()
    assert solver._sol is None
    again = similarity.march(*args, solver)
    assert again.separation_deg == ref.separation_deg
    np.testing.assert_allclose(again.q_c, ref.q_c, rtol=1e-9)
