"""The shared gas-loads provider (reentry_model.shared.gas_loads): the joined pressure and the meridian loads, checked
against the reference calculation's committed outputs (fixtures in tests/fixtures/gas_loads/)."""
import csv
import math
import os
import pickle
import time

import numpy as np
import pytest

from reentry_model.shared import aero, atmosphere, gas_loads, surface_flow, trajectory
from reentry_model.shared.constants import M_BAR_AIR
from reentry_model.shared.gas import UnreactedAir

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "gas_loads")
R_AIR, T_INF, V, RADIUS, T_W = 287.05, 200.0, 7108.0, 0.1524, 1500.0
CASES = [("0.002", 1.007e-4, 150.0), ("0.01", 2.014e-5, 150.0), ("0.05", 3.987e-6, 150.0), ("0.25", 8.057e-7, 157.0)]


def _rows(name, kn):
    with open(os.path.join(FIX, name)) as fh:
        return [r for r in csv.DictReader(fh) if r["case_kn"] == kn]


def aero_state(rho, v=V, T=T_INF):
    p = rho * R_AIR * T
    fs = atmosphere.Freestream(rho, T, p, M_BAR_AIR, np.zeros(3))
    kn = aero.knudsen(rho, M_BAR_AIR, 2 * RADIUS)
    return trajectory.AeroState(float("nan"), 0.0, 0.0, fs, np.array([v, 0.0, 0.0]), v, kn,
                                aero.mach(v, T, M_BAR_AIR), 1.0, np.zeros(3), 0.5 * rho * v**2)


@pytest.fixture(scope="module")
def fine():
    return gas_loads.GasLoads(bin_deg=1.0)


@pytest.fixture(scope="module")
def fine_loads(fine):
    return fine.meridian(aero_state(1.007e-4), RADIUS, T_W)


def test_joined_pressure_reproduces_reference():
    air = UnreactedAir()
    th = np.arange(0.0, 181.0)
    for kn, rho, meet_expect in CASES:
        p_inf, q_inf = rho * R_AIR * T_INF, 0.5 * rho * V**2
        ds = _rows("ds_model.csv", kn)
        p_ds = p_inf + np.array([float(r["Cp"]) for r in ds]) * q_inf
        stag = air.stagnation(rho, T_INF, V)
        p, meet = gas_loads.joined_pressure(th, stag.p, p_inf, q_inf, surface_flow.GAMMA_PM, p_ds)
        ref = np.array([float(r["p_over_pinf"]) for r in _rows("similarity_lee.csv", kn)])
        # ahead of the meet angle the pressure is the model's own; beyond it it is built from the 6-figure Cp
        assert p[:int(meet)] / p_inf == pytest.approx(ref[:int(meet)], rel=1e-9)
        assert p / p_inf == pytest.approx(ref, rel=1e-5)
        assert meet == meet_expect


def test_meridian_matches_reference_on_fine_bins(fine_loads):
    assert fine_loads.shear_source == "similarity"
    rho = 1.007e-4
    rows = _rows("similarity_heat.csv", "0.002")
    q3, q2 = 0.5 * rho * V**3, 0.5 * rho * V**2
    for i, r in enumerate(rows):
        if r["Cq_continuum"] == "":
            # blank fixture rows are the separated region (and the last bin, which is not marched)
            assert i >= fine_loads.separation_deg, i
            assert fine_loads.q_c[i] == 0.0 and fine_loads.tau[i] == 0.0, i
            continue
        tol = 1e-6 if i < fine_loads.meet_deg else 1e-5     # R5: 1e-5 only downstream of the D&S lee pressure
        assert fine_loads.q_c[i] / q3 == pytest.approx(float(r["Cq_continuum"]), rel=tol), i
        assert fine_loads.tau[i] / q2 == pytest.approx(float(r["Ctau_continuum"]), rel=tol, abs=1e-12), i
    assert fine_loads.q_c0 == fine_loads.q_c[0]
    assert fine_loads.separation_deg == 151.0
    assert fine_loads.meet_deg == 150.0


def test_default_bins_close_to_fine(fine_loads):
    rho = 1.007e-4
    coarse = gas_loads.GasLoads().meridian(aero_state(rho), RADIUS, T_W)
    th = np.radians(np.arange(0.0, 140.5, 0.5))
    for field in ("tau", "q_c"):
        a, b = coarse.at(th, field), fine_loads.at(th, field)
        assert a[1:] == pytest.approx(b[1:], rel=0.03), field
    assert abs(coarse.separation_deg - fine_loads.separation_deg) <= 2.0
    assert coarse.theta_deg[1] - coarse.theta_deg[0] == 2.0


def test_switch_to_daniels_schwartzentruber():
    rho = 3.987e-6
    state = aero_state(rho)
    loads = gas_loads.GasLoads(bin_deg=1.0).meridian(state, RADIUS, T_W)
    assert 0.05 <= state.kn and loads.kn_body == state.kn
    assert loads.shear_source == "daniels-schwartzentruber"
    ds = _rows("ds_model.csv", "0.05")
    expect = np.array([float(r["Ctau"]) for r in ds])
    assert loads.tau / loads.q_inf == pytest.approx(expect, rel=1e-5, abs=1e-9)
    assert loads.q_c0 > 0.0           # the similarity flux is computed at every Kn_body


def test_cache_shared():
    g = gas_loads.GasLoads()
    state = aero_state(1.007e-4)
    a = g.meridian(state, RADIUS, T_W)
    assert g.meridian(state, RADIUS, T_W) is a
    b = g.meridian(state, RADIUS, 1400.0)
    assert b is not a and b.T_wall == 1400.0


def test_low_speed_state():
    state = aero_state(0.0184, v=600.0, T=226.5)     # about 30 km, Mach about 2
    loads = gas_loads.GasLoads().meridian(state, 0.05, 400.0)
    for f in (loads.p_w, loads.tau, loads.q_c, loads.dp_ds):
        assert np.all(np.isfinite(f))
    assert loads.q_c0 > 0.0


def test_total_enthalpy_above_the_property_table_fails_the_meridian_with_the_freestream_and_the_limit():
    """At 8200 m/s the unreacted total enthalpy, 33.5 MJ/kg, is above the frozen-air table's top (32.2 MJ/kg at
    20 000 K): the meridian raises SimilarityError giving the freestream state (ruling R6) and both enthalpies, so a flight
    fails that step cleanly instead of flying on clipped properties."""
    g = gas_loads.GasLoads()
    state = aero_state(1e-5, v=8200.0)
    with pytest.raises(gas_loads.SimilarityError) as exc:
        g.meridian(state, RADIUS, T_W)
    msg = str(exc.value)
    assert "V=8200" in msg and "Kn_body={:.4g}".format(state.kn) in msg
    assert "exceeds the frozen-air property table" in msg and "{:.4g} J/kg".format(g._props.h_max) in msg
    assert g._cache == {}


def test_a_non_positive_joined_pressure_is_floored_before_cantera(monkeypatch):
    """The joined pressure is clamped to `P_FLOOR` times the freestream pressure before Cantera's isentropic expansion (and
    so in the loads returned): a non-positive value, which Daniels and Schwartzentruber's lee pressure could give where
    their c_p is below -p_inf/q_inf, no longer reaches Cantera, and every other bin is unchanged. The floor is far below
    any joined pressure met so far, so it changes no current result. Cantera itself cannot expand to the floor (the
    isentropic temperature goes below zero), so the recording stand-in below expands such a bin to the freestream
    pressure instead: the test asks only what pressure reaches Cantera."""
    state = aero_state(0.0184, v=600.0, T=226.5)     # about 30 km, Mach about 2: separated by 106 degrees
    plain = gas_loads.GasLoads().meridian(state, 0.05, 400.0)
    real_join, real_expand = gas_loads.joined_pressure, UnreactedAir.expand
    p_inf, seen = state.freestream.p, []

    def join(*a, **k):
        p, meet = real_join(*a, **k)
        p[-2:] = (0.0, -5.0)                          # 178 and 180 degrees
        return p, meet

    def expand(self, stag, p):
        seen.append(p)
        return real_expand(self, stag, p if p > 1e-3 * p_inf else p_inf)
    monkeypatch.setattr(gas_loads, "joined_pressure", join)
    monkeypatch.setattr(UnreactedAir, "expand", expand)
    loads = gas_loads.GasLoads().meridian(state, 0.05, 400.0)
    assert 0.0 < gas_loads.P_FLOOR <= 1e-20
    assert min(seen) == gas_loads.P_FLOOR * p_inf
    assert list(loads.p_w[-2:]) == [gas_loads.P_FLOOR * p_inf] * 2
    np.testing.assert_array_equal(loads.p_w[:-2], plain.p_w[:-2])
    np.testing.assert_array_equal(loads.tau, plain.tau)


def test_zero_speed_returns_zero_loads():
    loads = gas_loads.GasLoads().meridian(aero_state(1e-4, v=0.0), RADIUS, T_W)
    assert np.all(loads.tau == 0.0) and np.all(loads.q_c == 0.0) and np.all(loads.p_w == loads.p_inf)


def test_pickle_round_trip():
    g = gas_loads.GasLoads()
    state = aero_state(1.007e-4)
    a = g.meridian(state, RADIUS, T_W)
    g2 = pickle.loads(pickle.dumps(g))
    b = g2.meridian(state, RADIUS, T_W)
    assert b is not a
    for f in ("p_w", "tau", "q_c"):
        assert getattr(b, f) == pytest.approx(getattr(a, f), rel=1e-9)


def test_at_interpolates_and_clips():
    g = gas_loads.GasLoads()
    loads = g.meridian(aero_state(1.007e-4), RADIUS, T_W)
    th = loads.theta_deg
    mid = np.radians(0.5 * (th[10] + th[11]))
    assert loads.at(mid, "p_w") == pytest.approx(0.5 * (loads.p_w[10] + loads.p_w[11]))
    assert loads.at(-1.0, "p_w") == pytest.approx(loads.p_w[0])
    assert loads.at(10.0, "p_w") == pytest.approx(loads.p_w[-1])


def test_meridian_cost():
    g = gas_loads.GasLoads()
    g.meridian(aero_state(1.007e-4), RADIUS, T_W)
    t0 = time.time()
    g.meridian(aero_state(1.1e-4), RADIUS, T_W)
    dt = time.time() - t0
    print(f"meridian cost {dt:.2f} s")
    assert dt < 5.0


def test_meridian_cost_near_separation_failure():
    """The Kn 0.01 state (Kn_body 0.013) has a separation bin whose solve fails; it must fail fast."""
    g = gas_loads.GasLoads()
    g.meridian(aero_state(1.007e-4), RADIUS, T_W)
    t0 = time.time()
    g.meridian(aero_state(2.014e-5), RADIUS, T_W)
    dt = time.time() - t0
    print(f"Kn 0.01 meridian cost {dt:.2f} s")
    assert dt < 5.0
