"""Checks of the droplet model (from the prototype's verify_droplet_demise.py), the release= entry point, the parent
snapshot, and the reproduction of the 2026-10-22 demise sweep."""
import csv
import math
import os

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from reentry_model.droplet_reentry import model as dd
from reentry_model.droplet_reentry import parent_trajectory

PARENT_SHA256 = "0ba1cb3238654de5d1474600a32f4bc9f1626d705eac53decd0f322f73a27b7f"   # DRAMA run of 2026-09-13


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def old_sweep(fixtures_dir, column):
    """The 2026-10-14 free-molecule sweep (eps 0.1, accommodation 0.8) for 10 um, interpolated to 80 km."""
    rows = sorted(read_rows(os.path.join(fixtures_dir, "droplet_heating_sweep_2026-10-14_d10um.csv")),
                  key=lambda r: float(r["h0_km"]))
    return float(np.interp(80.0, [float(r["h0_km"]) for r in rows], [float(r[column]) for r in rows]))


def test_parent_snapshot_is_the_committed_drama_run(parent):
    assert parent.sha256 == PARENT_SHA256
    assert parent.object_name == "Main Body"


@pytest.mark.parametrize("el", dd.VOLATILE, ids=lambda e: e.name)
def test_vapour_pressure_is_one_atmosphere_at_the_boiling_point(el):
    assert dd.p_sat(el, el.T_boil) / 101325.0 == pytest.approx(1.0, abs=0.01)


@pytest.fixture(scope="module")
def sim80(parent):
    return dd.DropletSimulator(dd.parent_release(parent, 80.0), 10.0)


def initial_volatiles(sim):
    return np.array([sim.m0 * e.mass_fraction for e in dd.VOLATILE])


def test_rarefied_aluminium_evaporation_is_hertz_knudsen(sim80):
    fs = sim80.atmosphere.state(0.0, 80e3, 0.0, 0.0)
    m_vol = initial_volatiles(sim80)
    mdot, _ = sim80.evaporation(fs, 1.0, 1500.0, m_vol, 10e-6)
    x_al = sim80.mole_fractions(m_vol)[0]
    ref = math.pi * 1e-10 * x_al * dd.p_sat(dd.VOLATILE[0], 1500.0) * math.sqrt(
        0.026982 / (2 * math.pi * dd.R_GAS * 1500.0))
    assert mdot[0] / ref == pytest.approx(1.0, abs=0.01)


def test_zinc_and_magnesium_boil_at_release_at_80_km(sim80):
    fs = sim80.atmosphere.state(0.0, 80e3, 0.0, 0.0)
    assert sim80.evaporation(fs, 1.0, dd.T_RELEASE, initial_volatiles(sim80), 10e-6)[1] > 1.0


def test_diffusion_throttles_evaporation_below_boiling_at_sea_level(sim80):
    fs = sim80.atmosphere.state(0.0, 0.0, 0.0, 0.0)
    m_vol = initial_volatiles(sim80)
    mdot, y = sim80.evaporation(fs, 1.0, 1100.0, m_vol, 100e-6)
    hk = math.pi * 1e-8 * sim80.mole_fractions(m_vol)[1] * dd.p_sat(dd.VOLATILE[1], 1100.0) * math.sqrt(
        0.06538 / (2 * math.pi * dd.R_GAS * 1100.0))
    assert y < 1.0 and mdot[1] / hk < 0.5


def test_weber_precheck_skips_an_unstable_droplet(parent, monkeypatch):
    monkeypatch.setattr(dd, "WE_CRIT", 1e-3)
    s, rows = dd.run_case(10.0, 80.0, parent)
    assert s["skipped"] and rows is None


def test_without_evaporation_peak_temperature_matches_the_free_molecule_sweep(parent, fixtures_dir):
    s, _ = dd.run_case(10.0, 80.0, parent, alpha_evap=0.0, drag_law="blend")   # same drag as that sweep
    assert s["T_peak_K"] / old_sweep(fixtures_dir, "T_peak_K") == pytest.approx(1.0, abs=0.05)
    assert abs(s["mass_lost_frac"]) < 1e-9


def test_heat_models_agree_where_blended(parent):
    """Where a flight blends PhysicsHeating with the exchange model, the gap is under 30 % of the incident energy
    flux alpha 1/2 rho V^3 A (the relative gap is meaningless near recovery)."""
    gaps = []
    for d, h in ((10.0, 80.0), (100.0, 60.0), (200.0, 52.0)):
        _, rows = dd.run_case(d, h, parent)
        sim = dd.DropletSimulator(dd.parent_release(parent, h), d)
        for r in rows:
            if 0.0 < r["w_physics_heating"] < 1.0:
                fs = sim.atmosphere.state(0.0, r["h_km"] * 1e3, 0.0, 0.0)
                dm = r["d_um"] * 1e-6
                scale = dd.ALPHA_HEAT * 0.5 * fs.rho * r["V_ms"] ** 3 * math.pi * dm * dm / 4.0
                gaps.append(abs(sim.q_physics(fs, r["V_ms"], r["T_K"], dm)
                                - sim.q_exchange(fs, r["V_ms"], r["T_K"], dm)) / scale)
    assert gaps and max(gaps) < 0.3


def test_final_state_is_steady_and_element_losses_add_up(parent):
    s, _ = dd.run_case(10.0, 80.0, parent)
    assert s["end_reason"] == "steady"
    assert s["V_final_ms"] / s["v_terminal_final_ms"] == pytest.approx(1.0, abs=0.01)
    assert abs(s["T_final_K"] - s["T_eq_final_K"]) < 1.0
    lost = sum(s["lost_{}_kg".format(e.name)] for e in dd.VOLATILE)
    assert lost / s["mass_lost_kg"] == pytest.approx(1.0, abs=1e-9)


def test_release_state_given_directly_matches_the_parent_lookup(parent):
    """The entry point for the fragment model: an InitialState instead of a release altitude."""
    by_altitude, _ = dd.run_case(9.69, 70.0, parent)
    by_state, _ = dd.run_case(9.69, release=dd.parent_release(parent, 70.0))
    assert by_state["h0_km"] == pytest.approx(70.0)
    for key in ("t_final_s", "V_final_ms", "mass_lost_frac", "T_peak_K"):
        assert by_state[key] == by_altitude[key]


def test_settling_quadrature_matches_an_ode(parent):
    s, rows = dd.run_case(20.6, 70.0, parent, w_meso=0.05)
    sim = dd.DropletSimulator(dd.parent_release(parent, 70.0), 20.6,
                              dd.US76TableAtmosphere(wind=dd.MeanVerticalWind(0.05)))
    last = rows[-1]
    lat, lon = math.radians(last["lat_deg"]), math.radians(last["lon_deg"])

    def sink(t, y):
        h = max(y[0], 0.0)
        fs = sim.atmosphere.state(0.0, h, lat, lon)
        return [-(sim.terminal_velocity(fs, dd.earth.geodetic_to_ecef(h, lat, lon), fs.T, last["d_um"] * 1e-6,
                                        last["mass_kg"]) - sim.atmosphere.wind.w(h))]

    def ground(t, y):
        return y[0]
    ground.terminal, ground.direction = True, -1
    ode = solve_ivp(sink, (0.0, 1e9), [last["h_km"] * 1e3], method="LSODA", rtol=1e-8, atol=1e-3, events=ground)
    assert s["t_fall_s"] / ode.t_events[0][0] == pytest.approx(1.0, abs=1e-3)


def test_time_to_ground_matches_the_free_molecule_sweep(parent, fixtures_dir):
    s, _ = dd.run_case(10.0, 80.0, parent, w_meso=0.0, drag_law="blend", alpha_evap=0.0)
    assert s["t_ground_s"] / old_sweep(fixtures_dir, "t_ground_s") == pytest.approx(1.0, abs=0.01)


def test_droplet_is_held_where_terminal_velocity_meets_the_updraft(parent):
    s, _ = dd.run_case(1.0, 80.0, parent, w_meso=0.05, w_strat=3e-4)
    assert s["end_settle"] == "suspended"
    h = s["h_suspended_km"] * 1e3
    fs = dd.US76TableAtmosphere().state(0.0, h, 0.0, 0.0)
    sim = dd.DropletSimulator(dd.parent_release(parent, 80.0), 1.0)
    vt = sim.terminal_velocity(fs, dd.earth.geodetic_to_ecef(h, 0.0, 0.0), fs.T, s["d_final_um"] * 1e-6,
                               s["mass_final_kg"])
    assert vt / dd.MeanVerticalWind(0.05, 3e-4).w(h) == pytest.approx(1.0, abs=0.02)


def reference_cases():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "demise_sweep_2026-10-22_cases.csv")
    return read_rows(path)


@pytest.mark.parametrize("ref", reference_cases(), ids=lambda r: "{}um-{}km".format(r["d0_um"], r["h0_km"]))
def test_reproduces_the_2026_10_22_demise_sweep(parent, ref):
    """Cases of the prototype's sweep (runs/demise_sweep/final_states.csv), re-run from the package."""
    s, _ = dd.run_case(float(ref["d0_um"]), float(ref["h0_km"]), parent)
    assert s["end_reason"] == ref["end_reason"]
    for key in ref:
        if key not in ("d0_um", "h0_km", "end_reason"):
            assert s[key] == pytest.approx(float(ref[key]), rel=1e-9), key


def test_default_parent_is_the_data_snapshot():
    assert parent_trajectory.DEFAULT_PARENT.endswith(os.path.join("data", "generic_satellite_reentry",
                                                                  "dmf_output.json"))
    assert os.path.isfile(parent_trajectory.DEFAULT_PARENT)
