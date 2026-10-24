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
        # stepper state (advance): time, ECEF state vector, end reason once an event or t_max is reached
        self.t, self.y, self.end_reason, self.nfev = 0.0, self.initial_state_vector(), None, 0

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

    def advance(self, dt):
        """One macro step from the stored state: DOP853 over [t, t + dt] at the Step 1 tolerances, a fresh solve_ivp
        per step (spec section 7). A ground/escape event or t_max truncates the step and sets `end_reason`.
        Returns the time actually advanced."""
        if self.end_reason is not None:
            raise RuntimeError("the flight already ended ({})".format(self.end_reason))
        s = self.settings
        t_target = min(self.t + dt, s.t_max)
        sol = solve_ivp(self.rhs, (self.t, t_target), self.y, method="DOP853", rtol=s.rtol,
                        atol=[s.atol_position] * 3 + [s.atol_velocity] * 3, events=self._events())
        if not sol.success:
            raise RuntimeError("integration failed: {}".format(sol.message))
        self.nfev += int(sol.nfev)
        t_new, self.y = float(sol.t[-1]), sol.y[:, -1].copy()
        if sol.t_events[0].size:
            self.end_reason = "ground"
        elif sol.t_events[1].size:
            self.end_reason = "escape"
        elif t_new >= s.t_max:
            self.end_reason = "t_max"
        advanced, self.t = t_new - self.t, t_new
        return advanced

    def sample_row(self, t, y, a):
        """The Step 1 history columns for state y at time t with its AeroState a (the body supplies temperature/mass)."""
        r, v = y[:3], y[3:]
        # velocity_kms/flight_path_deg/heading_deg are kinematic and use the ground-relative
        # (rotating-frame) velocity v -- the quantity SESAM reports relative to the rotating
        # atmosphere -- while mach/knudsen/drag/dynamic_pressure_Pa use aero_state's
        # wind-relative v_rel; the two coincide unless a wind model is active.
        V, gamma, heading = earth.flight_angles(v, a.lat, a.lon)
        return {
            "time_s": t, "altitude_km": a.h / 1e3, "velocity_kms": V / 1e3,
            "temperature_K": self.body.mean_temperature(), "mass_kg": self.body.mass(t),
            "thick_mm": self.settings.diameter * 500.0,
            "lat_deg": math.degrees(a.lat), "lon_deg": math.degrees(a.lon),
            "downrange_km": earth.great_circle_distance(self.initial.lat, self.initial.lon, a.lat, a.lon) / 1e3,
            "flight_path_deg": math.degrees(gamma), "heading_deg": math.degrees(heading),
            "drag": a.cd, "lift": 0.0, "side": 0.0, "knudsen": a.kn, "mach": a.ma,
            "density_kgm3": a.freestream.rho, "dynamic_pressure_Pa": a.q_dyn,
            "load_factor_g": float(np.linalg.norm(a.a_drag)) / G0,     # = SESAM's column (drag / m g0, verified)
        }

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
            row = self.sample_row(t, y, self.aero_state(t, y[:3], y[3:]))
            for k in CSV_COLUMNS:
                cols[k].append(row[k])
            states.append(y)
        columns = {k: np.array(v, dtype=float) for k, v in cols.items()}
        history = History(columns, np.array(states), end_reason)
        history.results = self.results(history, int(sol.nfev), time.perf_counter() - started)
        return history

    def results(self, history, nfev, runtime):
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
            "rhs_evaluations": nfev,
            "runtime_s": runtime,
        }
        return results


def write_history_csv(history, path):
    """All of the history's columns, the Step 1 columns first, then any thermal columns (coupled runs)."""
    names = CSV_COLUMNS + [k for k in history.columns if k not in CSV_COLUMNS]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(names)
        for i in range(len(history)):
            w.writerow(["{:.9g}".format(history.columns[k][i]) for k in names])


def read_history_csv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    columns = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return History(columns, np.zeros((len(rows), 6)), "unknown")


def write_run_json(path, doc):
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2, default=str)
