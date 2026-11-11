"""earth.py: WGS-84 conversions, local frame, gravity."""
import math

import numpy as np
import pytest

from reentry_model import constants as c
from reentry_model import earth

R100 = dict(h=77500.133, lat=math.radians(29.546067), lon=math.radians(-82.134333))
R50 = dict(h=115000.0, lat=math.radians(29.546067), lon=math.radians(-82.134333))


@pytest.mark.parametrize("state", [R100, R50, dict(h=0.0, lat=0.0, lon=0.0), dict(h=1e5, lat=math.radians(-60.0), lon=math.radians(170.0))])
def test_geodetic_ecef_round_trip(state):
    r = earth.geodetic_to_ecef(state["h"], state["lat"], state["lon"])
    h, lat, lon = earth.ecef_to_geodetic(r)
    assert abs(h - state["h"]) < 1e-6
    assert abs(lat - state["lat"]) < 1e-12 and abs(lon - state["lon"]) < 1e-12


def test_equator_and_pole_radii():
    assert np.linalg.norm(earth.geodetic_to_ecef(0.0, 0.0, 0.0)) == pytest.approx(c.WGS84_A)
    b = c.WGS84_A * (1.0 - c.WGS84_F)
    assert earth.geodetic_to_ecef(0.0, math.pi / 2, 0.0)[2] == pytest.approx(b)


def test_enu_basis_is_orthonormal_and_up_is_radial_at_equator():
    e, n, u = earth.enu_basis(R100["lat"], R100["lon"])
    for a in (e, n, u):
        assert np.linalg.norm(a) == pytest.approx(1.0)
    assert abs(e @ n) < 1e-12 and abs(e @ u) < 1e-12 and abs(n @ u) < 1e-12
    assert np.allclose(np.cross(e, n), u)
    e0, n0, u0 = earth.enu_basis(0.0, 0.0)
    assert np.allclose(u0, [1, 0, 0]) and np.allclose(n0, [0, 0, 1]) and np.allclose(e0, [0, 1, 0])


def test_velocity_from_flight_angles_round_trip():
    V, gamma, heading = 7500.0, math.radians(-0.959331), math.radians(347.168296)
    v = earth.velocity_from_flight_angles(V, gamma, heading, R100["lat"], R100["lon"])
    assert np.linalg.norm(v) == pytest.approx(V)
    V2, gamma2, heading2 = earth.flight_angles(v, R100["lat"], R100["lon"])
    assert V2 == pytest.approx(V) and gamma2 == pytest.approx(gamma, abs=1e-12)
    assert heading2 == pytest.approx(heading, abs=1e-12)


def test_heading_convention_north_and_east():
    e, n, u = earth.enu_basis(R100["lat"], R100["lon"])
    north = earth.velocity_from_flight_angles(100.0, 0.0, 0.0, R100["lat"], R100["lon"])
    east = earth.velocity_from_flight_angles(100.0, 0.0, math.pi / 2, R100["lat"], R100["lon"])
    assert north @ n == pytest.approx(100.0) and abs(north @ e) < 1e-9
    assert east @ e == pytest.approx(100.0) and abs(east @ n) < 1e-9


def test_gravity_point_mass_and_j2_at_equator_and_pole():
    r_eq = earth.geodetic_to_ecef(0.0, 0.0, 0.0)
    g_point = earth.gravity(r_eq, "point")
    assert np.linalg.norm(g_point) == pytest.approx(c.MU_EARTH / c.WGS84_A ** 2)
    g_j2 = earth.gravity(r_eq, "j2")
    assert np.linalg.norm(g_j2) == pytest.approx(c.MU_EARTH / c.WGS84_A ** 2 * (1.0 + 1.5 * c.J2), rel=1e-9)
    assert g_j2 @ r_eq < 0                                  # points inward
    r_pole = earth.geodetic_to_ecef(0.0, math.pi / 2, 0.0)
    g_pole = earth.gravity(r_pole, "j2")
    b = np.linalg.norm(r_pole)
    expected = c.MU_EARTH / b ** 2 * (1.0 - 3.0 * c.J2 * (c.WGS84_A / b) ** 2)
    assert np.linalg.norm(g_pole) == pytest.approx(expected, rel=1e-9)
    assert abs(np.linalg.norm(g_pole) - 9.832) < 0.002       # textbook polar gravity


def test_gravity_j4_is_a_small_correction():
    r = earth.geodetic_to_ecef(77500.0, R100["lat"], R100["lon"])
    g2, g4 = earth.gravity(r, "j2"), earth.gravity(r, "j2j4")
    assert np.linalg.norm(g4 - g2) < 3e-5 * np.linalg.norm(g2)


def test_gravity_rejects_unknown_model():
    with pytest.raises(ValueError):
        earth.gravity(np.array([c.WGS84_A, 0.0, 0.0]), "j6")


def test_great_circle_distance():
    quarter = earth.great_circle_distance(0.0, 0.0, 0.0, math.pi / 2)
    assert quarter == pytest.approx(math.pi / 2 * c.WGS84_A)
    assert earth.great_circle_distance(R100["lat"], R100["lon"], R100["lat"], R100["lon"]) == 0.0
