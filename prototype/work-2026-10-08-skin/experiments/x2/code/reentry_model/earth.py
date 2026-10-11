"""WGS-84 geodesy, the local East-North-Up frame, gravity with J2/J4, great-circle distance.

Positions are ECEF metres; angles radians; heading is clockwise from north.
"""
import math

import numpy as np

from .constants import J2, J4, MU_EARTH, WGS84_A, WGS84_E2

GRAVITY_MODELS = ("point", "j2", "j2j4")


def geodetic_to_ecef(h, lat, lon):
    """Geodetic altitude [m], latitude and longitude [rad] -> ECEF position [m]."""
    s, c = math.sin(lat), math.cos(lat)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * s * s)
    return np.array([(n + h) * c * math.cos(lon), (n + h) * c * math.sin(lon), (n * (1.0 - WGS84_E2) + h) * s])


def ecef_to_geodetic(r):
    """ECEF position [m] -> (geodetic altitude [m], latitude [rad], longitude [rad]).

    Fixed-point iteration on the latitude (converges to < 1e-9 m in a few iterations for
    |lat| < 89 deg, which covers every trajectory we run)."""
    x, y, z = float(r[0]), float(r[1]), float(r[2])
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1.0 - WGS84_E2))
    h = 0.0
    for _ in range(10):
        s = math.sin(lat)
        n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * s * s)
        h = p / math.cos(lat) - n
        lat_new = math.atan2(z, p * (1.0 - WGS84_E2 * n / (n + h)))
        if abs(lat_new - lat) < 1e-15:
            lat = lat_new
            break
        lat = lat_new
    return h, lat, lon


def enu_basis(lat, lon):
    """Unit vectors (east, north, up) of the local frame at geodetic (lat, lon), in ECEF."""
    sl, cl = math.sin(lat), math.cos(lat)
    so, co = math.sin(lon), math.cos(lon)
    e = np.array([-so, co, 0.0])
    n = np.array([-sl * co, -sl * so, cl])
    u = np.array([cl * co, cl * so, sl])
    return e, n, u


def velocity_from_flight_angles(V, gamma, heading, lat, lon):
    """Speed [m/s], flight-path angle [rad, + up] and heading [rad, clockwise from north] -> ECEF velocity."""
    e, n, u = enu_basis(lat, lon)
    return V * (math.cos(gamma) * math.sin(heading) * e + math.cos(gamma) * math.cos(heading) * n + math.sin(gamma) * u)


def flight_angles(v, lat, lon):
    """ECEF velocity -> (speed, flight-path angle [rad], heading [rad] in [0, 2 pi))."""
    e, n, u = enu_basis(lat, lon)
    V = float(np.linalg.norm(v))
    if V == 0.0:
        return 0.0, 0.0, 0.0
    ve, vn, vu = float(v @ e), float(v @ n), float(v @ u)
    gamma = math.asin(max(-1.0, min(1.0, vu / V)))
    heading = math.atan2(ve, vn) % (2.0 * math.pi)
    return V, gamma, heading


def gravity(r, model="j2"):
    """Gravitational acceleration [m/s^2] at ECEF position r for a point mass, J2, or J2+J4 field."""
    if model not in GRAVITY_MODELS:
        raise ValueError("gravity model must be one of {}, got {!r}".format(GRAVITY_MODELS, model))
    x, y, z = float(r[0]), float(r[1]), float(r[2])
    r2 = x * x + y * y + z * z
    rr = math.sqrt(r2)
    base = -MU_EARTH / (rr * r2)
    if model == "point":
        return base * np.array([x, y, z])
    z2 = z * z / r2
    k2 = 1.5 * J2 * WGS84_A * WGS84_A / r2
    g = base * np.array([x * (1.0 - k2 * (5.0 * z2 - 1.0)),
                         y * (1.0 - k2 * (5.0 * z2 - 1.0)),
                         z * (1.0 - k2 * (5.0 * z2 - 3.0))])
    if model == "j2j4":
        c4 = 15.0 * MU_EARTH * J4 * WGS84_A ** 4 / (8.0 * rr ** 7)
        g = g + c4 * np.array([x * (1.0 - 14.0 * z2 + 21.0 * z2 * z2),
                               y * (1.0 - 14.0 * z2 + 21.0 * z2 * z2),
                               z * (5.0 - 70.0 * z2 / 3.0 + 21.0 * z2 * z2)])
    return g


def great_circle_distance(lat1, lon1, lat2, lon2, radius=WGS84_A):
    """Haversine distance [m] between two (lat, lon) [rad] points on a sphere of `radius`."""
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2.0 * radius * math.asin(min(1.0, math.sqrt(a)))
