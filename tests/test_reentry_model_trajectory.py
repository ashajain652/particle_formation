"""body.py and trajectory.py: the rotating-Earth 3-DOF integrator and its outputs."""
import csv
import json
import math
from datetime import datetime

import numpy as np
import pytest

from reentry_model import aero, atmosphere, body, earth
from reentry_model import trajectory as tj
from reentry_model.constants import MU_EARTH, OMEGA_EARTH

EPOCH = datetime(2024, 8, 1, 12, 53, 7)
R100 = tj.InitialState(7500.0, 77500.133, math.radians(-0.959331), math.radians(347.168296),
                       math.radians(29.546067), math.radians(-82.134333), EPOCH)
ORBIT = tj.InitialState(7800.0, 300e3, 0.0, math.radians(90.0), 0.0, 0.0, EPOCH)
MASS_100MM = 2813.0 * 4.0 / 3.0 * math.pi * 0.05 ** 3


def make(initial, atm, mass=MASS_100MM, **settings):
    s = tj.Settings(diameter=0.1, **settings)
    return tj.Simulator(initial, body.ConstantBody(mass, 300.0), atm, aero.SphereDragTables.from_json(), aero.SesamTable(), s)


def test_constant_body():
    assert body.sphere_mass(0.1, 2813.0) == pytest.approx(MASS_100MM)
    b = body.ConstantBody(1.5, 310.0)
    assert b.mass(0.0) == 1.5 and b.mass(100.0) == 1.5 and b.mean_temperature() == 310.0
    assert b.advance(0.0, 0.5, None) is None and b.field() is None and b.energy() == 0.0
    assert b.surface_temperature().tolist() == [310.0]


def test_initial_state_vector_round_trips_the_inputs():
    y0 = make(R100, atmosphere.VacuumAtmosphere()).initial_state_vector()
    h, lat, lon = earth.ecef_to_geodetic(y0[:3])
    assert h == pytest.approx(77500.133, abs=1e-6) and lat == pytest.approx(R100.lat) and lon == pytest.approx(R100.lon)
    V, gamma, heading = earth.flight_angles(y0[3:], lat, lon)
    assert V == pytest.approx(7500.0) and gamma == pytest.approx(R100.flight_path) and heading == pytest.approx(R100.heading)


def test_vacuum_point_mass_non_rotating_conserves_energy_and_angular_momentum():
    sim = make(ORBIT, atmosphere.VacuumAtmosphere(), gravity="point", rotating_frame=False,
               cadence=60.0, t_max=600.0, escape_altitude=1e9)
    hist = sim.run()
    assert hist.end_reason == "t_max" and hist.states.shape == (11, 6)
    r, v = hist.states[:, :3], hist.states[:, 3:]
    energy = 0.5 * np.sum(v * v, axis=1) - MU_EARTH / np.linalg.norm(r, axis=1)
    momentum = np.linalg.norm(np.cross(r, v), axis=1)
    # samples come from DOP853's dense output, whose interpolation error at rtol 1e-9 is ~2e-9
    assert np.max(np.abs(energy - energy[0]) / abs(energy[0])) < 1e-8
    assert np.max(np.abs(momentum - momentum[0]) / momentum[0]) < 1e-8


def test_vacuum_rotating_frame_conserves_the_jacobi_integral():
    sim = make(ORBIT, atmosphere.VacuumAtmosphere(), gravity="point", rotating_frame=True,
               cadence=60.0, t_max=600.0, escape_altitude=1e9)
    hist = sim.run()
    r, v = hist.states[:, :3], hist.states[:, 3:]
    jacobi = 0.5 * np.sum(v * v, axis=1) - MU_EARTH / np.linalg.norm(r, axis=1) \
        - 0.5 * OMEGA_EARTH ** 2 * (r[:, 0] ** 2 + r[:, 1] ** 2)
    # samples come from DOP853's dense output, whose interpolation error at rtol 1e-9 is ~2e-9
    assert np.max(np.abs(jacobi - jacobi[0]) / abs(jacobi[0])) < 1e-8


def test_flight_path_angle_rate_matches_sesam_and_vinh():
    # SESAM: -0.00790 deg/s at the R100 start; Vinh's rotating-Earth equations: -0.00792 deg/s (facts note s.7).
    sim = make(R100, atmosphere.VacuumAtmosphere(), cadence=2.5, t_max=2.5, escape_altitude=1e9)
    hist = sim.run()
    gamma = hist.columns["flight_path_deg"]
    rate = (gamma[-1] - gamma[0]) / hist.columns["time_s"][-1]
    assert rate == pytest.approx(-0.0079, abs=0.0002)


def test_us76_flight_ends_on_the_ground_and_samples_as_requested():
    sim = make(R100, atmosphere.US76TableAtmosphere(), cadence=10.0)
    hist = sim.run(extra_times=[3.3, 7.7])
    t = hist.columns["time_s"]
    assert hist.end_reason == "ground"
    assert abs(hist.columns["altitude_km"][-1]) < 1e-6                 # within 1 mm of h = 0
    assert 300.0 < t[-1] < 450.0                                        # SESAM: 366 s on US76
    assert t[0] == 0.0 and np.all(np.diff(t) > 0)
    for wanted in (3.3, 7.7, 10.0, 20.0):
        assert np.any(np.isclose(t, wanted))
    assert hist.results["impact_time_s"] == t[-1] and hist.results["end_reason"] == "ground"
    assert 60.0 < hist.results["knudsen_crossings"]["0.01"] < 76.0     # km; SESAM: 70.1 km
    assert np.all(hist.columns["mass_kg"] == MASS_100MM) and np.all(hist.columns["temperature_K"] == 300.0)
    assert hist.columns["drag"][0] == pytest.approx(aero.drag_coefficient(hist.columns["knudsen"][0], hist.columns["mach"][0],
                                                                            aero.SphereDragTables.from_json(), aero.SesamTable()))
    assert 6.0 < hist.results["max_deceleration_g"] < 10.0 and 35.0 < hist.results["altitude_of_max_deceleration_km"] < 50.0   # SESAM: 8.2 g at 41 km


def test_escape_event_ends_the_flight_at_the_escape_altitude():
    climbing = tj.InitialState(7800.0, 140e3, math.radians(2.0), math.radians(90.0), 0.0, 0.0, EPOCH)
    sim = make(climbing, atmosphere.VacuumAtmosphere(), cadence=5.0, t_max=600.0)
    hist = sim.run()
    assert hist.end_reason == "escape"
    assert abs(hist.columns["altitude_km"][-1] - 150.0) < 1e-6          # within 1 mm of h = 150 km
    assert 20.0 < hist.columns["time_s"][-1] < 60.0                     # climb rate ~= 272 m/s
    assert hist.results["impact_time_s"] is None
    assert hist.results["final_altitude_km"] == pytest.approx(150.0, abs=1e-6)


def test_start_at_or_above_escape_altitude_is_rejected():
    for h in (150e3, 200e3):
        above = tj.InitialState(7800.0, h, math.radians(2.0), math.radians(90.0), 0.0, 0.0, EPOCH)
        with pytest.raises(ValueError):
            make(above, atmosphere.VacuumAtmosphere())


def test_aero_state_tolerates_trial_stage_altitudes_outside_the_table():
    # DOP853's adaptive RK stages can trial-evaluate the RHS a little beyond the ground or the top
    # of a hard-bounded table (e.g. US76's [0, 150000] m) before backing off; aero_state must not
    # raise on those excursions, and must still report the TRUE (unclamped) altitude.
    sim = make(R100, atmosphere.US76TableAtmosphere())
    v = earth.velocity_from_flight_angles(R100.velocity, R100.flight_path, R100.heading, R100.lat, R100.lon)
    for h in (-5.0, 150_500.0):
        r = earth.geodetic_to_ecef(h, R100.lat, R100.lon)
        a = sim.aero_state(0.0, r, v)
        assert a.h == pytest.approx(h, abs=1e-3)
        assert np.all(np.isfinite(a.a_drag))


def test_csv_and_json_round_trip(tmp_path):
    sim = make(R100, atmosphere.US76TableAtmosphere(), cadence=50.0, t_max=100.0)
    hist = sim.run()
    csv_path = tmp_path / "h.csv"
    tj.write_history_csv(hist, str(csv_path))
    with open(csv_path) as fh:
        rows = list(csv.DictReader(fh))
    assert list(rows[0].keys()) == tj.CSV_COLUMNS and len(rows) == len(hist.columns["time_s"])
    back = tj.read_history_csv(str(csv_path))
    assert np.allclose(back.columns["velocity_kms"], hist.columns["velocity_kms"], rtol=1e-6)
    tj.write_run_json(str(tmp_path / "h.json"), {"results": hist.results, "epoch": EPOCH})
    doc = json.load(open(tmp_path / "h.json"))
    assert doc["results"]["end_reason"] == "t_max" and doc["epoch"].startswith("2024-08-01")


def test_advance_reproduces_run_and_truncates_at_the_ground():
    """The stepper (a fresh DOP853 solve per macro step) vs the single dense-output solve of run(): < 0.01 m/s, < 0.1 m."""
    atm = atmosphere.US76TableAtmosphere()
    hist = make(R100, atm, cadence=0.5, t_max=60.0).run()
    sim = make(R100, atm, t_max=60.0)
    assert sim.t == 0.0 and sim.end_reason is None
    rows = [sim.sample_row(sim.t, sim.y, sim.aero_state(sim.t, sim.y[:3], sim.y[3:]))]
    while sim.end_reason is None:
        assert sim.advance(0.5) == pytest.approx(0.5)
        rows.append(sim.sample_row(sim.t, sim.y, sim.aero_state(sim.t, sim.y[:3], sim.y[3:])))
    assert sim.end_reason == "t_max" and sim.t == 60.0 and len(rows) == len(hist) and sim.nfev > 0
    V = np.array([r["velocity_kms"] for r in rows]) * 1e3
    h = np.array([r["altitude_km"] for r in rows]) * 1e3
    assert np.abs(V - hist.columns["velocity_kms"] * 1e3).max() < 0.01
    assert np.abs(h - hist.columns["altitude_km"] * 1e3).max() < 0.1
    with pytest.raises(RuntimeError):
        sim.advance(0.5)
    # ground event truncates the final step
    low = tj.InitialState(300.0, 200.0, math.radians(-60.0), R100.heading, R100.lat, R100.lon, EPOCH)
    sim = make(low, atm)
    advanced = sim.advance(10.0)
    assert sim.end_reason == "ground" and 0.5 < advanced < 1.5 and abs(earth.ecef_to_geodetic(sim.y[:3])[0]) < 1e-3


def test_history_csv_round_trips_extra_columns(tmp_path):
    hist = make(R100, atmosphere.US76TableAtmosphere(), cadence=5.0, t_max=10.0).run()
    hist.columns["convective_heat_W"] = np.array([1.0, 2.0, 3.0])
    tj.write_history_csv(hist, tmp_path / "h.csv")
    with open(tmp_path / "h.csv") as fh:
        header = fh.readline().strip().split(",")
    assert header[:len(tj.CSV_COLUMNS)] == tj.CSV_COLUMNS and header[-1] == "convective_heat_W"
    back = tj.read_history_csv(tmp_path / "h.csv")
    assert back.columns["convective_heat_W"].tolist() == [1.0, 2.0, 3.0] and len(back) == 3
