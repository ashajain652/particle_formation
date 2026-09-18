"""Rotating-Earth 3-DOF trajectory of a sphere, integrated in ECEF Cartesian coordinates.

d r/dt = v ;  d v/dt = g(r) - 2 w x v - w x (w x r) + a_drag,  a_drag = -1/2 rho |v_rel| v_rel C_D A / m
(spec section 6). Positions/velocities are relative to the rotating Earth; v_rel subtracts the wind.
"""
import csv
import json
import math
import time
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
from scipy.integrate import solve_ivp

from . import aero, earth
from .constants import G0, OMEGA_EARTH

CSV_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "temperature_K", "mass_kg", "thick_mm", "lat_deg",
               "lon_deg", "downrange_km", "flight_path_deg", "heading_deg", "drag", "lift", "side", "knudsen",
               "mach", "density_kgm3", "dynamic_pressure_Pa", "load_factor_g"]
KNUDSEN_THRESHOLDS = (10.0, 1.0, 0.1, 0.01)


@dataclass(frozen=True)
class InitialState:
    velocity: float      # m/s relative to the rotating atmosphere
    altitude: float      # m, geodetic
    flight_path: float   # rad, negative = descending
    heading: float       # rad, clockwise from north
    lat: float           # rad, geodetic
    lon: float           # rad
    epoch: datetime


@dataclass
class Settings:
    diameter: float                  # m
    gravity: str = "j2"              # point | j2 | j2j4
    rotating_frame: bool = True
    rtol: float = 1e-9
    atol_position: float = 1e-6      # m
    atol_velocity: float = 1e-9      # m/s
    cadence: float = 1.0             # s between history samples
    t_max: float = 3600.0            # s
    ground_altitude: float = 0.0     # m
    escape_altitude: float = 150e3   # m


@dataclass
class AeroState:
    h: float
    lat: float
    lon: float
    freestream: object
    v_rel: np.ndarray
    V: float
    kn: float
    ma: float
    cd: float
    a_drag: np.ndarray
    q_dyn: float


@dataclass
class History:
    columns: dict
    states: np.ndarray
    end_reason: str
    results: dict = field(default_factory=dict)

    def __len__(self):
        return len(self.columns["time_s"])


class Simulator:
    def __init__(self, initial, body, atmosphere, tables, bridging, settings):
        if initial.altitude >= settings.escape_altitude:
            raise ValueError("initial altitude {:.0f} m must be below the escape altitude {:.0f} m".format(
                initial.altitude, settings.escape_altitude))
        self.initial, self.body, self.atmosphere = initial, body, atmosphere
        self.tables, self.bridging, self.settings = tables, bridging, settings
        self.area = math.pi * settings.diameter ** 2 / 4.0
        self.omega = np.array([0.0, 0.0, OMEGA_EARTH]) if settings.rotating_frame else np.zeros(3)
        self.r0 = earth.geodetic_to_ecef(initial.altitude, initial.lat, initial.lon)
        self.v0 = earth.velocity_from_flight_angles(initial.velocity, initial.flight_path, initial.heading,
                                                    initial.lat, initial.lon)

    def initial_state_vector(self):
        return np.concatenate([self.r0, self.v0])

    def aero_state(self, t, r, v):
        h, lat, lon = earth.ecef_to_geodetic(r)
        # DOP853's adaptive stages can trial-evaluate the RHS a little beyond the ground or the
        # top of a hard-bounded table (e.g. US76's [0, 150000] m) before backing off. Clamp ONLY
        # the altitude used for the atmosphere.state() lookup to [ground_altitude, escape_altitude]
        # so that transient excursion doesn't raise; `h` itself (and thus AeroState.h, used for
        # the reported/history altitude) stays the TRUE, unclamped geodetic altitude.
        # The upper clamp only ever sees this trial-stage overshoot: __init__ rejects any initial
        # altitude at or above escape_altitude, so a real (non-trial) state can't reach here already
        # past the top of the table.
        h_atm = min(max(h, self.settings.ground_altitude), self.settings.escape_altitude)
        fs = self.atmosphere.state(t, h_atm, lat, lon)
        e, n, u = earth.enu_basis(lat, lon)
        wind = fs.wind_enu[0] * e + fs.wind_enu[1] * n + fs.wind_enu[2] * u
        v_rel = v - wind
        V = float(np.linalg.norm(v_rel))
        ma = aero.mach(V, fs.T, fs.m_bar)
        if fs.rho > 0.0:
            kn = aero.knudsen(fs.rho, fs.m_bar, self.settings.diameter)
            cd = aero.drag_coefficient(kn, ma, self.tables, self.bridging)
            a_drag = -0.5 * fs.rho * V * v_rel * cd * self.area / self.body.mass(t)
        else:
            kn, cd, a_drag = math.inf, self.tables.cd_free_molecular(ma), np.zeros(3)
        return AeroState(h, lat, lon, fs, v_rel, V, kn, ma, cd, a_drag, 0.5 * fs.rho * V * V)

    def rhs(self, t, y):
        r, v = y[:3], y[3:]
        a = earth.gravity(r, self.settings.gravity) - 2.0 * np.cross(self.omega, v) \
            - np.cross(self.omega, np.cross(self.omega, r)) + self.aero_state(t, r, v).a_drag
        return np.concatenate([v, a])

    def _events(self):
        s = self.settings

        def ground(t, y):
            return earth.ecef_to_geodetic(y[:3])[0] - s.ground_altitude

        def escape(t, y):
            return s.escape_altitude - earth.ecef_to_geodetic(y[:3])[0]

        ground.terminal, ground.direction = True, -1
        escape.terminal, escape.direction = True, -1
        return [ground, escape]

    def run(self, extra_times=None):
        s = self.settings
        started = time.perf_counter()
        sol = solve_ivp(self.rhs, (0.0, s.t_max), self.initial_state_vector(), method="DOP853",
                        rtol=s.rtol, atol=[s.atol_position] * 3 + [s.atol_velocity] * 3,
                        dense_output=True, events=self._events())
        if not sol.success:
            raise RuntimeError("integration failed: {}".format(sol.message))
        t_end = float(sol.t[-1])
        if sol.t_events[0].size:
            end_reason = "ground"
        elif sol.t_events[1].size:
            end_reason = "escape"
        else:
            end_reason = "t_max"
        times = np.arange(0.0, t_end, s.cadence)
        if extra_times is not None:
            extra = np.asarray(extra_times, dtype=float)
            times = np.union1d(times, extra[extra <= t_end])
        times = np.union1d(times, [t_end])
        cols = {k: [] for k in CSV_COLUMNS}
        states = []
        for t in times:
            y = sol.sol(t)
            r, v = y[:3], y[3:]
            a = self.aero_state(t, r, v)
            self.body.on_step(t, y, a.freestream, a)
            # velocity_kms/flight_path_deg/heading_deg are kinematic and use the ground-relative
            # (rotating-frame) velocity v -- the quantity SESAM reports relative to the rotating
            # atmosphere -- while mach/knudsen/drag/dynamic_pressure_Pa below use aero_state's
            # wind-relative v_rel; the two coincide unless a wind model is active.
            V, gamma, heading = earth.flight_angles(v, a.lat, a.lon)
            cols["time_s"].append(t)
            cols["altitude_km"].append(a.h / 1e3)
            cols["velocity_kms"].append(V / 1e3)
            cols["temperature_K"].append(self.body.temperature(t))
            cols["mass_kg"].append(self.body.mass(t))
            cols["thick_mm"].append(s.diameter * 500.0)
            cols["lat_deg"].append(math.degrees(a.lat))
            cols["lon_deg"].append(math.degrees(a.lon))
            cols["downrange_km"].append(earth.great_circle_distance(self.initial.lat, self.initial.lon, a.lat, a.lon) / 1e3)
            cols["flight_path_deg"].append(math.degrees(gamma))
            cols["heading_deg"].append(math.degrees(heading))
            cols["drag"].append(a.cd)
            cols["lift"].append(0.0)
            cols["side"].append(0.0)
            cols["knudsen"].append(a.kn)
            cols["mach"].append(a.ma)
            cols["density_kgm3"].append(a.freestream.rho)
            cols["dynamic_pressure_Pa"].append(a.q_dyn)
            cols["load_factor_g"].append(float(np.linalg.norm(a.a_drag)) / G0)     # = SESAM's column (drag / m g0, verified)
            states.append(y)
        columns = {k: np.array(v, dtype=float) for k, v in cols.items()}
        history = History(columns, np.array(states), end_reason)
        history.results = self._results(history, sol, time.perf_counter() - started)
        return history

    def _results(self, history, sol, runtime):
        c = history.columns
        i_dec = int(np.argmax(c["load_factor_g"]))
        crossings = {}
        for thr in KNUDSEN_THRESHOLDS:
            below = np.nonzero(c["knudsen"] < thr)[0]
            crossings["{:g}".format(thr)] = float(c["altitude_km"][below[0]]) if below.size else None
        results = {
            "end_reason": history.end_reason,
            "final_time_s": float(c["time_s"][-1]),
            "impact_time_s": float(c["time_s"][-1]) if history.end_reason == "ground" else None,
            "final_velocity_kms": float(c["velocity_kms"][-1]),
            "final_altitude_km": float(c["altitude_km"][-1]),
            "max_deceleration_g": float(c["load_factor_g"][i_dec]),
            "altitude_of_max_deceleration_km": float(c["altitude_km"][i_dec]),
            "time_of_max_deceleration_s": float(c["time_s"][i_dec]),
            "max_dynamic_pressure_Pa": float(c["dynamic_pressure_Pa"].max()),
            "knudsen_start": float(c["knudsen"][0]),
            "knudsen_crossings": crossings,
            "n_samples": len(history),
            "rhs_evaluations": int(sol.nfev),
            "runtime_s": runtime,
        }
        return results


def write_history_csv(history, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLUMNS)
        for i in range(len(history)):
            w.writerow(["{:.9g}".format(history.columns[k][i]) for k in CSV_COLUMNS])


def read_history_csv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    columns = {k: np.array([float(r[k]) for r in rows]) for k in CSV_COLUMNS if k in rows[0]}
    return History(columns, np.zeros((len(rows), 6)), "unknown")


def write_run_json(path, doc):
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2, default=str)
