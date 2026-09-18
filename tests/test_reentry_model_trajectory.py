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
    assert b.mass(0.0) == 1.5 and b.mass(100.0) == 1.5 and b.temperature(5.0) == 310.0
    assert b.on_step(0.0, None, None, None) is None


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
    assert np.max(np.abs(energy - energy[0]) / abs(energy[0])) < 1e-9
    assert np.max(np.abs(momentum - momentum[0]) / momentum[0]) < 1e-9


def test_vacuum_rotating_frame_conserves_the_jacobi_integral():
    sim = make(ORBIT, atmosphere.VacuumAtmosphere(), gravity="point", rotating_frame=True,
               cadence=60.0, t_max=600.0, escape_altitude=1e9)
    hist = sim.run()
    r, v = hist.states[:, :3], hist.states[:, 3:]
    jacobi = 0.5 * np.sum(v * v, axis=1) - MU_EARTH / np.linalg.norm(r, axis=1) \
        - 0.5 * OMEGA_EARTH ** 2 * (r[:, 0] ** 2 + r[:, 1] ** 2)
    assert np.max(np.abs(jacobi - jacobi[0]) / abs(jacobi[0])) < 1e-9


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
