"""Adapters between the repository's integrators and the overlap engine: the twin of the DRAMA parent (spec 4.3).

The twin parent
---------------
The DRAMA file (`data/generic_satellite_reentry/dmf_output.json`, object `Main Body`) is too coarse to place a sphere
relative to the plume: its rows are 5 s apart and its latitude and longitude are rounded to 0.001 degree (spec 4.1).
The parent track is therefore recomputed by this repository's own rotating-Earth integrator (`trajectory.Simulator`:
ECEF state, J2 gravity, Coriolis and centrifugal terms, DOP853), the same code that integrates the follower spheres,
so that parent and follower are integrated alike. The file calibrates the twin; it does not supply its geometry.

What is taken from the file, and what is assumed
  * Start state: `model.parent_release` at the top of the descending branch (about 111.2 km, 7.918 km/s, flight-path
    angle -0.226 degrees), i.e. the file's own row there, with t = 0 of the twin at that row's time.
  * Drag: a_drag = -(1/2) rho V v_rel / beta(h), with beta(h) the parent's effective ballistic coefficient
    m / (C_D A) tabulated against altitude as dynamicPressure / (loadFactor * G0) from the file's rows on the
    descending branch (`loadFactor` is SESAM's drag deceleration over g0; checked in trajectory.sample_row). This
    folds the parent's ablation (its mass falls from 740 kg at 81 km to about 216 kg at 66 km) and its changing
    drag coefficient into one number, so the twin needs neither a mass model nor a drag table.
  * Smoothing of beta: rows with loadFactor below 1e-4 are dropped (the ratio of two rounded small numbers), and the
    rest are smoothed with a centred running median of 5 rows (the window shrinks at the ends of the branch). The
    median removes single-row jumps from the rounding of `loadFactor` (4 significant figures) and `dynamicPressure`
    without shifting the steep fall of beta between 80 and 60 km the way a running mean would. beta, mass and
    C_D A = mass / beta are then linear in altitude between rows and clamped at the ends of the table.
  * Air: the US76 table with no wind (spec 4.3), which is also what DRAMA used for this run (its `density` column
    agrees with the US76 table to within 0.5 percent from 111 to 42 km, measured 2026-10-29). The plume frame is the
    Earth-fixed frame.
  * Flow arrays carried by the track (for the plume laws): the air's rho, T and p; the viscosity from the
    repository's Sutherland law `particle_trajectory.viscosity` (the form the repository's drag code uses: 1.716e-5
    Pa s at 273.15 K, S = 110.4 K; Feldman's form in `laws.sutherland_mu` uses 1.72e-5, 273 K and 113 K and differs
    by under 1 percent between 180 and 300 K); the hard-sphere mean free path `aero.mean_free_path`; the parent's
    equal-area diameter d and radius r0 from the constant frontal area 6.215 m2 (spec 4.2); and C_D A = m / beta
    from the table.
  * The twin is integrated once, with dense output, from the branch top down to `h_stop_km` (default 45 km, below the
    lowest release altitude of 50 km), and sampled every `cadence` seconds into a `Track`. A release altitude outside
    that range raises ValueError (Review Focus 5); nothing is extrapolated.

`twin_vs_drama` is the consistency check of spec 4.3: speed against altitude (target within 1 percent above 60 km)
and altitude against time.
"""
import json
import math
import time
from dataclasses import dataclass, field, replace
from datetime import timedelta

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from reentry_model import aero, earth
from reentry_model.atmosphere import US76TableAtmosphere
from reentry_model.body import ConstantBody
from reentry_model.constants import G0
from reentry_model.droplet_reentry import parent_trajectory
from reentry_model.droplet_reentry.model import parent_release
from reentry_model.droplet_reentry.particle_trajectory import viscosity
from reentry_model.trajectory import InitialState, Settings, Simulator
from reentry_model.wake_residence.track import FLOW_KEYS, Track

PARENT_FRONTAL_AREA = 6.215     # m2, DRAMA's impacting cross-section of the parent (spec 4.2)
LOAD_FACTOR_MIN = 1e-4          # rows with a smaller loadFactor (g) are dropped from the beta table
TWIN_T_MAX = 3600.0             # s, integration limit for the twin (the flight to 45 km takes about 950 s)


# ---------------------------------------------------------------- the ballistic coefficient from the DRAMA file

@dataclass
class BetaTable:
    """Parent ballistic coefficient beta = m / (C_D A) [kg/m2] and mass [kg] against altitude h_m [m] (increasing),
    linear in altitude and clamped at the ends."""
    h_m: np.ndarray
    beta: np.ndarray
    mass: np.ndarray

    def beta_at(self, h):
        return np.interp(h, self.h_m, self.beta)

    def mass_at(self, h):
        return np.interp(h, self.h_m, self.mass)

    def cda_at(self, h):
        """Effective drag area C_D A = mass / beta [m2]."""
        return self.mass_at(h) / self.beta_at(h)


def _running_median(x, window):
    """Centred running median of `window` (odd) points; the window shrinks symmetrically at the ends."""
    if window < 1 or window % 2 == 0:
        raise ValueError(f"smooth_window must be a positive odd number of rows, got {window!r}")
    k = window // 2
    n = len(x)
    out = np.empty(n)
    for i in range(n):
        w = min(k, i, n - 1 - i)
        out[i] = np.median(x[i - w:i + w + 1])
    return out


def _aerothermal_history(parent):
    """The `_AeroThermalHistory.txt` rows of the parent's object, read from the parent's own file."""
    with open(parent.path) as fh:
        doc = json.load(fh)
    results = doc["singleModuleOutputs"]["satellites"][0]["missionPhases"][0]["epochs"][0]["analysisModules"][0][
        "results"]
    wanted = parent.object_name.replace(" ", "_")
    keys = [k for k in results if k.endswith("_AeroThermalHistory.txt")
            and len(k.split(".")) >= 2 and k.split(".")[1] == wanted]
    if len(keys) != 1:
        raise ValueError(f"expected one aerothermal history for object {parent.object_name!r} in {parent.path}, "
                         f"found {keys}")
    return results[keys[0]]


def drama_beta_table(parent, smooth_window=5):
    """BetaTable from the descending branch of the DRAMA parent: beta = dynamicPressure / (loadFactor * G0), rows
    with loadFactor below LOAD_FACTOR_MIN dropped, a centred running median of `smooth_window` rows (module
    docstring), and the mass from the same object's aerothermal history interpolated to each row's time."""
    rows = [r for r in parent.branch if r["loadFactor"] >= LOAD_FACTOR_MIN]
    if len(rows) < 2:
        raise ValueError("fewer than 2 rows of the descending branch have a usable loadFactor")
    t = np.array([r["time"] for r in rows], dtype=float)
    h = np.array([r["altitude"] for r in rows], dtype=float) * 1e3
    beta = _running_median(np.array([r["dynamicPressure"] / (r["loadFactor"] * G0) for r in rows]), smooth_window)
    hist = _aerothermal_history(parent)
    t_hist = np.array([r["time"] for r in hist], dtype=float)
    m_hist = np.array([r["mass"] for r in hist], dtype=float)
    mass = np.interp(t, t_hist, m_hist)
    order = np.argsort(h, kind="stable")
    h, beta, mass = h[order], beta[order], mass[order]
    keep = np.concatenate([[True], np.diff(h) > 0.0])          # strictly increasing altitude nodes
    return BetaTable(h[keep], beta[keep], mass[keep])


# ---------------------------------------------------------------- the twin's integrator

class TwinSimulator(Simulator):
    """`trajectory.Simulator` with the drag replaced by -(1/2) rho V v_rel / (beta_scale * beta(h)).

    Everything else (gravity, rotating frame, atmosphere, events) is the parent class's. The body is a
    `ConstantBody(1.0)` and the diameter is `settings.diameter`; the parent class uses them only for its own drag,
    which is replaced here. In the returned AeroState, `kn` is the Knudsen number on `settings.diameter` and `cd` the
    effective drag coefficient m / (beta A) on that diameter's area.
    """

    def __init__(self, initial, beta_table, atmosphere, settings, beta_scale=1.0):
        if not (math.isfinite(beta_scale) and beta_scale > 0.0):
            raise ValueError(f"beta_scale must be a finite positive number, got {beta_scale!r}")
        super().__init__(initial, ConstantBody(1.0), atmosphere, aero.SphereDragTables.from_json(),
                         aero.bridging_by_name("sesam-table"), settings)
        self.beta_table, self.beta_scale = beta_table, float(beta_scale)

    def beta(self, h):
        return self.beta_scale * float(self.beta_table.beta_at(h))

    def aero_state(self, t, r, v):
        a = super().aero_state(t, r, v)
        beta = self.beta(a.h)
        a_drag = -0.5 * a.freestream.rho * a.V * a.v_rel / beta
        cd = float(self.beta_table.mass_at(a.h)) / (beta * self.area)
        return replace(a, a_drag=a_drag, cd=cd)


@dataclass
class DenseRun:
    sol: object            # scipy OdeResult with dense output (sol.sol(t) -> (6,) or (6, N))
    t_end: float           # s
    end_reason: str        # name of the terminal event that ended the run, or "t_max"
    runtime_s: float = math.nan


def dense_solve(sim, t_max, events=None, method="DOP853", rtol=None, atol=None):
    """One solve_ivp of the simulator's right-hand side from t = 0 with dense output, at the simulator's tolerances by
    default. `events` defaults to the simulator's ground and escape events."""
    s = sim.settings
    rtol = s.rtol if rtol is None else rtol
    atol = [s.atol_position] * 3 + [s.atol_velocity] * 3 if atol is None else atol
    events = events or sim._events()
    started = time.perf_counter()
    sol = solve_ivp(sim.rhs, (0.0, float(t_max)), sim.initial_state_vector(), method=method, rtol=rtol, atol=atol,
                    dense_output=True, events=events)
    runtime = time.perf_counter() - started
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    end_reason = "t_max"
    for ev, te in zip(events, sol.t_events):
        if getattr(ev, "terminal", False) and te.size:
            end_reason = getattr(ev, "__name__", "event")
            break
    return DenseRun(sol, float(sol.t[-1]), end_reason, runtime)


def flow_arrays(t, states, atmosphere, d, cda, r0):
    """Per-sample flow state for a Track: rho, T, p, mu (Sutherland, `particle_trajectory.viscosity`), mfp
    (`aero.mean_free_path`), d, cda (scalar or array), r0, plus `speed`, the air-relative speed. `states` is (N, 6)
    ECEF position and velocity."""
    t = np.asarray(t, dtype=float)
    states = np.asarray(states, dtype=float)
    n = len(t)
    out = {k: np.empty(n) for k in ("rho", "T", "p", "mu", "mfp", "speed")}
    for i in range(n):
        r, v = states[i, :3], states[i, 3:]
        h, lat, lon = earth.ecef_to_geodetic(r)
        fs = atmosphere.state(float(t[i]), h, lat, lon)
        e, nn, u = earth.enu_basis(lat, lon)
        wind = fs.wind_enu[0] * e + fs.wind_enu[1] * nn + fs.wind_enu[2] * u
        out["rho"][i], out["T"][i], out["p"][i] = fs.rho, fs.T, fs.p
        out["mu"][i] = viscosity(fs.T)
        out["mfp"][i] = aero.mean_free_path(fs.rho, fs.m_bar)
        out["speed"][i] = float(np.linalg.norm(v - wind))
    for k, val in (("d", d), ("cda", cda), ("r0", r0)):
        out[k] = np.broadcast_to(np.asarray(val, dtype=float), (n,)).copy()
    return out


# ---------------------------------------------------------------- the twin parent

@dataclass
class TwinParent:
    """The twin's Track (t = 0 at the top of DRAMA's descending branch, at `epoch`), its dense run, the beta table it
    was integrated with, and the SHA-256 of the DRAMA file."""
    track: Track
    run: DenseRun
    beta: BetaTable
    epoch: object
    sha256: str
    _grid: tuple = field(default=None, init=False, repr=False, compare=False)

    def altitude(self, t):
        """Geodetic altitude [m] of the track at time t [s]."""
        return earth.ecef_to_geodetic(self.track.position(float(t)))[0]

    def t_at_altitude(self, h_km):
        """Time [s] at which the descending track first reaches altitude h_km. ValueError naming the track's altitude
        range if h_km is outside it (nothing is extrapolated)."""
        h = float(h_km) * 1e3
        if self._grid is None:
            tg = np.linspace(self.track.t0, self.track.t1, max(int(self.track.t1 - self.track.t0) + 2, 2))
            self._grid = (tg, np.array([self.altitude(x) for x in tg]))
        tg, hg = self._grid
        h_top, h_bottom = hg[0], hg[-1]
        tol = 1e-6                                        # m: round-off of the geodetic conversion at the ends
        if not (math.isfinite(h) and h_bottom - tol <= h <= h_top + tol):
            raise ValueError(f"altitude {h_km!r} km is outside the twin's integrated range {h_bottom / 1e3:.3f} to "
                             f"{h_top / 1e3:.3f} km")
        if h >= h_top:
            return float(tg[0])
        if h <= h_bottom:
            return float(tg[-1])
        i =int(np.argmax(hg <= h))                       # first grid point at or below h (descending track)
        return float(brentq(lambda x: self.altitude(x) - h, tg[i - 1], tg[i], xtol=1e-10, rtol=1e-14))


def _stop_event(h_stop_m):
    def stop_altitude(t, y):
        return earth.ecef_to_geodetic(y[:3])[0] - h_stop_m
    stop_altitude.terminal, stop_altitude.direction = True, -1
    return stop_altitude


def twin_parent(parent=None, h_stop_km=45.0, atmosphere=None, frontal_area=PARENT_FRONTAL_AREA, beta_scale=1.0,
                cadence=0.02):
    """Integrate the twin from the DRAMA state at the top of the descending branch down to h_stop_km, and sample its
    dense output every `cadence` seconds into a Track carrying the flow arrays (module docstring)."""
    parent = parent if parent is not None else parent_trajectory.load_parent()
    atmosphere = atmosphere if atmosphere is not None else US76TableAtmosphere()
    if not (cadence > 0.0 and math.isfinite(cadence)):
        raise ValueError(f"cadence must be a finite positive number of seconds, got {cadence!r}")
    beta = drama_beta_table(parent)
    initial = parent_release(parent, parent.branch[0]["altitude"])
    d = math.sqrt(4.0 * frontal_area / math.pi)
    r0 = math.sqrt(frontal_area / math.pi)
    settings = Settings(diameter=d, t_max=TWIN_T_MAX)
    sim = TwinSimulator(initial, beta, atmosphere, settings, beta_scale=beta_scale)
    run = dense_solve(sim, TWIN_T_MAX, events=sim._events() + [_stop_event(h_stop_km * 1e3)])
    if run.end_reason != "stop_altitude":
        raise RuntimeError(f"the twin ended by {run.end_reason!r} at t = {run.t_end:.3f} s before reaching "
                           f"{h_stop_km} km")
    t = np.arange(0.0, run.t_end, cadence)
    if run.t_end - t[-1] < 1e-3 * cadence:
        t = t[:-1]
    t = np.append(t, run.t_end)
    states = run.sol.sol(t).T
    h = np.array([earth.ecef_to_geodetic(s[:3])[0] for s in states])
    cda = np.asarray(beta.mass_at(h)) / (beta_scale * np.asarray(beta.beta_at(h)))
    flow = flow_arrays(t, states, atmosphere, d, cda, r0)
    speed = flow.pop("speed")
    track = Track.from_arrays(t, states[:, :3], states[:, 3:], speed=speed, flow={k: flow[k] for k in FLOW_KEYS})
    return TwinParent(track, run, beta, initial.epoch, parent.sha256)


def release_state(track, t, epoch):
    """InitialState of a follower released at time t exactly on the track (position and velocity of the track),
    with epoch + t as its epoch."""
    r, v = track.position(float(t)), track.velocity(float(t))
    h, lat, lon = earth.ecef_to_geodetic(r)
    V, gamma, heading = earth.flight_angles(v, lat, lon)
    return InitialState(velocity=V, altitude=h, flight_path=gamma, heading=heading, lat=lat, lon=lon,
                        epoch=epoch + timedelta(seconds=float(t)))


def twin_vs_drama(twin, parent):
    """One dict per DRAMA row on the descending branch: the DRAMA time, altitude, speed and downrange, the twin's speed
    at the DRAMA altitude, the twin's altitude at the DRAMA time, and speed_rel_diff = twin / DRAMA speed - 1. Twin
    values outside the twin's range are NaN."""
    offset = (twin.epoch - parent.epoch).total_seconds()      # DRAMA time of the twin's t = 0
    out = []
    for r in parent.branch:
        try:
            v_twin = float(twin.track.speed(twin.t_at_altitude(r["altitude"])))
        except ValueError:
            v_twin = math.nan
        tt = r["time"] - offset
        h_twin = twin.altitude(tt) / 1e3 if twin.track.t0 <= tt <= twin.track.t1 else math.nan
        v_drama = r["velocity"] * 1e3
        out.append({"time_s": r["time"], "drama_altitude_km": r["altitude"], "drama_speed_ms": v_drama,
                    "drama_downrange_km": r["downrange"], "twin_speed_ms": v_twin, "twin_altitude_km": h_twin,
                    "speed_rel_diff": v_twin / v_drama - 1.0})
    return out
