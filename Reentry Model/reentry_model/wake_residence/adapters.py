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

Followers in the wake (spec 4.3 and 4.3a)
-----------------------------------------
  * `WakeAtmosphere` wraps any atmosphere. Its `state(t, h, lat, lon)` (t on the follower's clock) asks the plume for
    the wake gas at the follower (`Plume.wake_state` at parent time t + t_offset, emission window anchored at the
    release, t_release = t_offset) and returns the wrapped atmosphere's own `Freestream` object when the mode is "off"
    (the plume is then not queried at all) or the gas there is plain ambient air (beyond the profiles' cutoff, or no
    law at the best emission point); otherwise the wake's density, temperature and pressure, the wrapped atmosphere's
    mean molecular mass, and its wind plus the wake gas velocity converted from Earth-fixed to east-north-up.
    `trajectory.Simulator` and `DropletSimulator` read the air only through `state()`, so neither is modified. A
    follower time beyond the end of the parent's track raises ValueError naming both ranges (the plume is undefined
    there); `follower_t_max` is the longest follower integration the parent track supports.
  * `sphere_simulator` builds a sphere follower exactly as `cli.py` builds a `--thermal none` run.
  * `solution_track` samples a dense run every `cadence` seconds (default 0.002 s, plus the end time) into a Track on
    the parent's clock; `sphere_parent_track` adds the flow arrays a sphere parent's plume needs, with
    C_D A = C_D(Kn, Ma) times the sphere's area from the simulator's own drag tables and bridging.

Rarefied followers in the parent's shock layer (spec 4.3d)
---------------------------------------------------------
  * `WakeAtmosphere(..., shock_layer=field)` also applies a `shock_layer.ShockLayerField`, in any mode (mode "off"
    applies the field and no wake). The field lives in the parent frame (`parent_frame`: x along the parent's
    velocity at the parent time, y the distance from the axis) and moves and turns with the twin track; its gas is
    blended with the wake-or-ambient gas by the field's weight w: rho, T, p and m_bar linearly, and the Earth-fixed
    gas velocity as w (parent velocity + field velocity) + (1 - w) (wind + wake velocity). `pressure_gradient(t, x)` is
    the gradient of the field's excess pressure in ECEF (zeros without a field); `last_field_weight` is w at the last
    `state()` call.
  * `RarefiedSphereSimulator` (`rarefied_sphere_simulator`) is a sphere follower whose drag is the composed law of
    `droplet_reentry.particle_trajectory.drag_coefficient` (Schaaf and Chambre free-molecular, a continuum
    correlation, Cunningham slip, blended on the speed ratio) relative to the local gas, with the wall-to-gas
    temperature ratio from a named wall temperature (300 K, assumed), zero at zero relative speed, plus the
    pressure-gradient acceleration -grad(p) / rho_f. Its AeroState's `a_drag` is the sum of the two (the whole
    aerodynamic acceleration the integrator uses), so `Simulator.sample_row`'s drag and load-factor columns include the
    pressure-gradient force for this simulator; `max_a_drag` and `max_a_pressure` record the largest of each seen
    inside the field (weight above 0.5), counting trial evaluations of the integrator.
  * `parent_collision_event` stops a run when the follower's distance from the parent's centroid falls below R (the
    parent as an obstacle); `shock_layer_exit` classifies the run: the first of `through_shock` (x >= x_s(y)) and
    `past_base_plane` (x below -1 micrometre inside the envelope) on the dense output, else `hit_parent` if the
    collision event ended the run, else `inside_at_end`.
"""
import json
import math
import time
from dataclasses import dataclass, field, replace
from datetime import timedelta

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from reentry_model import aero, body, earth
from reentry_model.atmosphere import Freestream, US76TableAtmosphere
from reentry_model.body import ConstantBody
from reentry_model.constants import G0
from reentry_model.droplet_reentry import parent_trajectory
from reentry_model.droplet_reentry.model import parent_release
from reentry_model.droplet_reentry.particle_trajectory import drag_coefficient as composed_drag_coefficient
from reentry_model.droplet_reentry.particle_trajectory import viscosity
from reentry_model.trajectory import AeroState, InitialState, Settings, Simulator
from reentry_model.wake_residence.track import FLOW_KEYS, Track
from reentry_model.wake_residence.wake_gas import Ambient, check_mode

PARENT_FRONTAL_AREA = 6.215     # m2, DRAMA's impacting cross-section of the parent (spec 4.2)
AA7075_DENSITY = 2813.0         # kg/m3, AA7075 solid (DRAMA Sphere Experiments/data/user_materials/AA7075_nomelt.json)
SOLUTION_CADENCE = 0.002        # s, default sampling of a follower's dense run into a Track
LOAD_FACTOR_MIN = 1e-4          # rows with a smaller loadFactor (g) are dropped from the beta table
TWIN_T_MAX = 3600.0             # s, integration limit for the twin (the flight to 45 km takes about 950 s)
HIT_DEPTH = 1e-7                # the collision event fires below a distance R (1 - HIT_DEPTH) from the centroid: well
                                # above ECEF round-off (about 1e-9 m on 6.4e6 m), so a follower released on the surface
                                # never starts with the event already negative
BASE_PLANE_TOL = 1e-6           # m: `past_base_plane` needs x below -BASE_PLANE_TOL, so a release on the base plane
                                # (theta = 90 degrees) is not labelled by round-off at its first sample


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


# ---------------------------------------------------------------- followers in the wake (spec 4.3a)

class WakeAtmosphere:
    """An atmosphere whose `state()` is the wrapped atmosphere's, with the parent's wake gas inside the plume and,
    optionally, the parent's shock-layer gas around it.

    base         the atmosphere wrapped (anything with `state(t, h, lat, lon) -> Freestream`)
    plume        the parent's `overlap.Plume`
    t_offset     the follower's release time on the parent's clock [s]: follower time t is parent time t + t_offset,
                 and the plume's emission window is anchored there (spec 4.4)
    mode         `wake_flow`: "off", "velocity", "density" or "both" (spec 4.3a)
    k_E          scale on the line-blast radius of the energy profile (`wake_gas.gas_state`)
    shock_layer  None, or the `shock_layer.ShockLayerField` of this release (spec 4.3d), applied in every mode

    `n_calls` counts `state()` calls; `last_margin` is the `MarginResult` of the last plume query (None in mode "off");
    `last_field_weight` is the shock-layer field's weight at the last `state()` call (0 without a field).
    """

    def __init__(self, base, plume, t_offset, mode, k_E=1.0, shock_layer=None):
        check_mode(mode)
        t_offset, k_E = float(t_offset), float(k_E)
        if not math.isfinite(t_offset):
            raise ValueError(f"t_offset must be a finite number of seconds, got {t_offset!r}")
        if not (math.isfinite(k_E) and k_E > 0.0):
            raise ValueError(f"k_E must be a finite positive number, got {k_E!r}")
        self.base, self.plume, self.t_offset, self.mode, self.k_E = base, plume, t_offset, mode, k_E
        self.shock_layer = shock_layer
        self.n_calls = 0
        self.last_margin = None
        self.last_field_weight = 0.0

    @property
    def follower_t_max(self):
        """Longest follower time [s] for which the parent's track defines the plume (parent t1 - t_offset)."""
        return self.plume.parent.t1 - self.t_offset

    def _parent_time(self, t):
        """Parent time of follower time t; ValueError naming both ranges outside the parent's track."""
        tp = float(t) + self.t_offset
        par = self.plume.parent
        slack = 1e-9 * (par.t1 - par.t0)
        if not (par.t0 - slack <= tp <= par.t1 + slack):
            raise ValueError(f"follower time {float(t)!r} s is parent time {tp!r} s, outside the parent track's range "
                             f"[{par.t0!r}, {par.t1!r}] s (follower range [{par.t0 - self.t_offset!r}, "
                             f"{self.follower_t_max!r}] s)")
        return tp

    def parent_frame(self, t, x_ecef):
        """(x, y, e_x, e_y) of the ECEF point x_ecef at follower time t: x [m] along the parent's velocity from its
        centroid, y [m] the distance from that axis, e_x the parent's direction of flight and e_y the unit lateral
        vector towards the point (any perpendicular to e_x on the axis)."""
        tp = self._parent_time(t)
        par = self.plume.parent
        v = par.velocity(tp)
        e_x = v / np.linalg.norm(v)
        d = np.asarray(x_ecef, dtype=float) - par.position(tp)
        x = float(d @ e_x)
        perp = d - x * e_x
        y = float(np.linalg.norm(perp))
        if y > 0.0:
            e_y = perp / y
        else:
            a = np.cross(e_x, [0.0, 0.0, 1.0])
            if not np.linalg.norm(a) > 1e-6:
                a = np.cross(e_x, [1.0, 0.0, 0.0])
            e_y = a / np.linalg.norm(a)
        return x, y, e_x, e_y

    def pressure_gradient(self, t, x_ecef):
        """Gradient [Pa/m] in ECEF of the shock-layer field's excess pressure, weight * (p - p_inf), at the ECEF point
        x_ecef and follower time t (zeros without a field)."""
        if self.shock_layer is None:
            return np.zeros(3)
        x, y, e_x, e_y = self.parent_frame(t, x_ecef)
        gx, gy = self.shock_layer.excess_pressure_gradient(x, y)
        return gx * e_x + gy * e_y

    def state(self, t, h, lat, lon):
        self.n_calls += 1
        fs = self.base.state(t, h, lat, lon)
        gas = fs if self.mode == "off" else self._wake_state(t, h, lat, lon, fs)
        if self.shock_layer is None:
            return gas
        return self._with_shock_layer(t, h, lat, lon, gas)

    def _with_shock_layer(self, t, h, lat, lon, gas):
        """`gas` (wake or ambient) blended with the shock-layer field by the field's weight w: rho, T, p and m_bar
        linearly, the Earth-fixed velocity as w (parent velocity + field velocity) + (1 - w) (gas velocity)."""
        x, y, e_x, e_y = self.parent_frame(t, earth.geodetic_to_ecef(h, lat, lon))
        st = self.shock_layer.state(x, y)
        w = st.weight
        self.last_field_weight = w
        if w == 0.0:
            return gas
        e, n, up = earth.enu_basis(lat, lon)
        v_field = self.plume.parent.velocity(float(t) + self.t_offset) + st.w_x * e_x + st.w_y * e_y
        v_gas = gas.wind_enu[0] * e + gas.wind_enu[1] * n + gas.wind_enu[2] * up
        v = w * v_field + (1.0 - w) * v_gas
        mix = lambda a, b: w * a + (1.0 - w) * b
        return Freestream(mix(st.rho, gas.rho), mix(st.T, gas.T), mix(st.p, gas.p), mix(st.m_bar, gas.m_bar),
                          np.array([v @ e, v @ n, v @ up]))

    def _wake_state(self, t, h, lat, lon, fs):
        """The wake gas at the point (Task 6), or `fs` itself where the gas there is plain ambient air."""
        tp = self._parent_time(t)
        x = earth.geodetic_to_ecef(h, lat, lon)
        gas, m = self.plume.wake_state(x, tp, Ambient(fs.p, fs.T, fs.rho), self.mode, self.k_E,
                                       t_release=self.t_offset)
        self.last_margin = m
        u = np.asarray(gas.u, dtype=float)
        if gas.rho == fs.rho and gas.T == fs.T and gas.p == fs.p and not np.any(u):
            return fs                                    # plain ambient air: beyond the cutoff, or no law there
        e, n, up = earth.enu_basis(lat, lon)
        u_enu = np.array([u @ e, u @ n, u @ up])
        return Freestream(gas.rho, gas.T, gas.p, fs.m_bar, np.asarray(fs.wind_enu, dtype=float) + u_enu)


def sphere_simulator(initial, d_m, atmosphere, density=AA7075_DENSITY, bridging="sesam-table", t_max=600.0):
    """A solid sphere of diameter d_m [m] and density [kg/m3] from `initial`, built as `cli.py` builds a
    `--thermal none` run: ConstantBody of the sphere's mass, the SESAM sphere drag tables, the named bridging, J2
    gravity, rotating frame, and the repository's default tolerances."""
    d_m = float(d_m)
    if not (math.isfinite(d_m) and d_m > 0.0):
        raise ValueError(f"d_m must be a finite positive diameter in metres, got {d_m!r}")
    settings = Settings(diameter=d_m, gravity="j2", rotating_frame=True, t_max=float(t_max))
    return Simulator(initial, ConstantBody(body.sphere_mass(d_m, density)), atmosphere,
                     aero.SphereDragTables.from_json(), aero.bridging_by_name(bridging), settings)


class RarefiedSphereSimulator(Simulator):
    """A solid sphere follower for the shock layer (spec 4.3d): `trajectory.Simulator` with the drag replaced by
    -(1/2) rho V v_rel C_D A / m, C_D from the composed law `particle_trajectory.drag_coefficient` with the local gas
    and the wall-to-gas temperature ratio wall_T_K / T (zero force at zero relative speed), plus the pressure-gradient
    acceleration -grad(p) / density read from `atmosphere.pressure_gradient(t, r)` when `pressure_gradient` is on and
    the atmosphere has that method. The AeroState's `a_drag` is their sum (the whole aerodynamic acceleration, so
    history columns labelled drag or load factor include a_p), `cd` the composed C_D (infinite at zero relative
    speed) and `kn` the hard-sphere Knudsen number. `max_a_drag` and
    `max_a_pressure` [m/s2] are the largest drag and pressure-gradient accelerations seen where the field's weight
    exceeds 0.5 (trial evaluations of the integrator included). The SESAM drag tables are kept on the simulator
    (`tables`, `bridging`) for consistency checks only.
    """

    def __init__(self, initial, d_m, atmosphere, density=AA7075_DENSITY, wall_T_K=300.0, pressure_gradient=True,
                 t_max=600.0):
        d_m, density, wall_T_K = float(d_m), float(density), float(wall_T_K)
        for name, val in (("d_m", d_m), ("density", density), ("wall_T_K", wall_T_K)):
            if not (math.isfinite(val) and val > 0.0):
                raise ValueError(f"{name} must be a finite positive number, got {val!r}")
        settings = Settings(diameter=d_m, gravity="j2", rotating_frame=True, t_max=float(t_max))
        super().__init__(initial, ConstantBody(body.sphere_mass(d_m, density)), atmosphere,
                         aero.SphereDragTables.from_json(), aero.bridging_by_name("sesam-table"), settings)
        self.density, self.wall_T_K, self.pressure_gradient = density, wall_T_K, bool(pressure_gradient)
        self.max_a_drag = 0.0
        self.max_a_pressure = 0.0

    def aero_state(self, t, r, v):
        s = self.settings
        h, lat, lon = earth.ecef_to_geodetic(r)
        fs = self.atmosphere.state(t, min(max(h, s.ground_altitude), s.escape_altitude), lat, lon)
        in_field = getattr(self.atmosphere, "last_field_weight", 0.0) > 0.5
        e, n, u = earth.enu_basis(lat, lon)
        v_rel = v - (fs.wind_enu[0] * e + fs.wind_enu[1] * n + fs.wind_enu[2] * u)
        V = float(np.linalg.norm(v_rel))
        if V > 0.0 and fs.rho > 0.0:
            cd = composed_drag_coefficient(V, fs.rho, fs.T, fs.p, s.diameter,
                                           temperature_ratio=self.wall_T_K / fs.T)[0]
            a_drag = -0.5 * fs.rho * V * v_rel * cd * self.area / self.body.mass(t)
        else:
            cd, a_drag = math.inf, np.zeros(3)
        a_aero = a_drag
        if self.pressure_gradient and hasattr(self.atmosphere, "pressure_gradient"):
            a_p = -np.asarray(self.atmosphere.pressure_gradient(t, r), dtype=float) / self.density
            a_aero = a_drag + a_p
            if in_field:
                self.max_a_pressure = max(self.max_a_pressure, float(np.linalg.norm(a_p)))
        if in_field:
            self.max_a_drag = max(self.max_a_drag, float(np.linalg.norm(a_drag)))
        kn = aero.knudsen(fs.rho, fs.m_bar, s.diameter) if fs.rho > 0.0 else math.inf
        return AeroState(h, lat, lon, fs, v_rel, V, kn, aero.mach(V, fs.T, fs.m_bar), cd, a_aero,
                         0.5 * fs.rho * V * V)


def rarefied_sphere_simulator(initial, d_m, atmosphere, density=AA7075_DENSITY, wall_T_K=300.0,
                              pressure_gradient=True, t_max=600.0):
    """A `RarefiedSphereSimulator` of diameter d_m [m] and density [kg/m3] from `initial`."""
    return RarefiedSphereSimulator(initial, d_m, atmosphere, density=density, wall_T_K=wall_T_K,
                                   pressure_gradient=pressure_gradient, t_max=t_max)


def parent_collision_event(wake_atmosphere, R):
    """solve_ivp event `hit_parent`: the follower's distance from the parent's centroid minus R (1 - HIT_DEPTH),
    terminal, on inward crossings only. A run it ends has `end_reason` "hit_parent"."""
    track, t_offset, R = wake_atmosphere.plume.parent, wake_atmosphere.t_offset, float(R)

    def hit_parent(t, y):
        return float(np.linalg.norm(y[:3] - track.position(t + t_offset))) - R * (1.0 - HIT_DEPTH)
    hit_parent.terminal, hit_parent.direction = True, -1
    return hit_parent


@dataclass
class ShockLayerExit:
    kind: str       # "through_shock", "past_base_plane", "hit_parent" or "inside_at_end"
    t: float        # follower time [s]
    x: float        # parent frame: along the flight direction from the centroid [m]
    y: float        # parent frame: distance from the axis [m]
    v_rel: float    # follower speed relative to the parent [m/s]


def shock_layer_exit(run, wake_atmosphere, field, hit=False, tol=1e-3):
    """How the follower of a dense run left the parent's shock layer (spec 4.3d): the first of `through_shock`
    (x >= x_s(y)) or `past_base_plane` (x < -BASE_PLANE_TOL inside the envelope), found on the dense output sampled at
    most `tol` seconds apart and refined by root-finding; else `hit_parent` if `hit` (the collision event ended the
    run); else `inside_at_end`. An exit followed by a hit is reported as the exit (the integration continues after an
    exit). The tolerance keeps a release on the base plane (theta = 90 degrees) from being labelled at t0 by
    round-off."""
    wa, shape = wake_atmosphere, field.shape
    track = wa.plume.parent
    t0 = float(run.sol.t[0])
    ts = np.linspace(t0, run.t_end, max(int(math.ceil((run.t_end - t0) / tol)), 1) + 1)
    states = np.atleast_2d(run.sol.sol(ts).T)[:, :6]
    v_par = track.velocity(ts + wa.t_offset)
    e_x = v_par / np.linalg.norm(v_par, axis=1)[:, None]
    d = states[:, :3] - track.position(ts + wa.t_offset)
    xs = np.einsum("ij,ij->i", d, e_x)
    ys = np.linalg.norm(d - xs[:, None] * e_x, axis=1)

    def frame(t):
        st = run.sol.sol(t)[:6]
        x, y, _, _ = wa.parent_frame(t, st[:3])
        return x, y, float(np.linalg.norm(st[3:] - track.velocity(t + wa.t_offset)))

    def found(kind, t):
        return ShockLayerExit(kind, float(t), *frame(t))

    for i in range(len(ts)):
        if xs[i] >= shape.x_s(ys[i]):
            kind, f = "through_shock", lambda t: (lambda x, y, _: x - shape.x_s(y))(*frame(t))
        elif xs[i] < -BASE_PLANE_TOL:
            kind, f = "past_base_plane", lambda t: -frame(t)[0] - BASE_PLANE_TOL
        else:
            continue
        return found(kind, ts[0] if i == 0 else brentq(f, ts[i - 1], ts[i], xtol=1e-9))
    return found("hit_parent" if hit else "inside_at_end", run.t_end)


def _solution_samples(run, cadence):
    """(t, states): the dense run sampled every `cadence` s from its start, plus its end time; states (N, 6) are the
    first six components (ECEF position and velocity), so a droplet's longer state vector works too."""
    cadence = SOLUTION_CADENCE if cadence is None else float(cadence)
    if not (math.isfinite(cadence) and cadence > 0.0):
        raise ValueError(f"cadence must be a finite positive number of seconds, got {cadence!r}")
    t0 = float(run.sol.t[0])
    t = t0 + np.arange(0.0, run.t_end - t0, cadence)
    if t.size and run.t_end - t[-1] < 1e-3 * cadence:
        t = t[:-1]
    t = np.append(t, run.t_end)
    states = np.atleast_2d(run.sol.sol(t).T)[:, :6]
    return t, states


def solution_track(run, t_offset, flow=None, cadence=None):
    """Track of a dense run on the parent's clock (times shifted by t_offset), sampled every `cadence` s (default
    0.002 s, plus the end time), with the run's positions and velocities as a Hermite spline.

    `flow` is None, a dict of per-sample arrays for exactly the Track's flow keys (optionally with `speed`, the
    air-relative speed), or a callable `flow(t, states)` returning such a dict, called with the sample times on the
    run's own clock and the (N, 6) states. Without `speed` the Track's speed is the Earth-fixed speed.
    """
    t, states = _solution_samples(run, cadence)
    speed = None
    if flow is not None:
        flow = dict(flow(t, states) if callable(flow) else flow)
        speed = flow.pop("speed", None)
    return Track.from_arrays(t + float(t_offset), states[:, :3], states[:, 3:6], speed=speed, flow=flow)


def sphere_parent_track(run, sim, t_offset, atmosphere, cadence=None):
    """Track of a sphere integrated by `sim` (a `sphere_simulator`), as the parent of a plume: `solution_track` with
    flow arrays from `atmosphere` (rho, T, p, viscosity, mean free path, air-relative speed), d the sphere diameter,
    r0 = d / 2 and C_D A = C_D(Kn, Ma) times the sphere's area from the simulator's own tables and bridging. The
    altitude used for the air is clamped to the simulator's ground and escape altitudes, as `Simulator.aero_state`
    does, so the sample at a ground event (a fraction of a millimetre below 0 m) is valid."""
    s = sim.settings
    d, area = float(s.diameter), float(sim.area)

    def flow(t, states):
        n = len(t)
        out = {k: np.empty(n) for k in ("rho", "T", "p", "mu", "mfp", "speed", "cda")}
        for i in range(n):
            r, v = states[i, :3], states[i, 3:6]
            h, lat, lon = earth.ecef_to_geodetic(r)
            fs = atmosphere.state(float(t[i]), min(max(h, s.ground_altitude), s.escape_altitude), lat, lon)
            e, nn, up = earth.enu_basis(lat, lon)
            wind = fs.wind_enu[0] * e + fs.wind_enu[1] * nn + fs.wind_enu[2] * up
            V = float(np.linalg.norm(v - wind))
            ma = aero.mach(V, fs.T, fs.m_bar)
            if fs.rho > 0.0:
                cd = aero.drag_coefficient(aero.knudsen(fs.rho, fs.m_bar, d), ma, sim.tables, sim.bridging)
            else:
                cd = sim.tables.cd_free_molecular(ma)
            out["rho"][i], out["T"][i], out["p"][i] = fs.rho, fs.T, fs.p
            out["mu"][i] = viscosity(fs.T)
            out["mfp"][i] = aero.mean_free_path(fs.rho, fs.m_bar)
            out["speed"][i], out["cda"][i] = V, cd * area
        out["d"], out["r0"] = np.full(n, d), np.full(n, 0.5 * d)
        return out

    return solution_track(run, t_offset, flow=flow, cadence=cadence)
