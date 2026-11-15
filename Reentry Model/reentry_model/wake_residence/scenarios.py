"""Scenarios: one follower released from one parent, integrated, and measured against the parent's plume (spec 4.3b-4.4).

`sphere_in_spacecraft_plume(d_mm, h_km, theta_deg)` releases a solid AA7075 sphere from the twin of the Generic
Spacecraft (`adapters.twin_parent`) at the time the twin first reaches h_km, at release angle theta (90 to 180 degrees,
the windward side only, Asha). The order of decisions is that of spec 4.3c:

  1. The parent's regime at release (`laws.classify`). A rarefied or too-slow parent has no plume law, no shock layer
     and no AFFIx outcome: nothing is integrated and the reason is `rarefied_no_law` or `too_slow_no_law`.
  2. The parent's shock layer is built from its free stream at release (`shock_layer.ShockLayerField`, R = D1/2 with
     D1 the parent's equal-area diameter). The follower's Knudsen number Kn_f is the mean free path of the shock-layer
     gas at the release point (`field.surface_state(theta)`) over the follower's diameter, and sets the follower's
     regime (`follower_regime`). AFFIx (`affix.affix_regime`) classifies every release; its status is `governing`,
     `indicative` or `not_applicable` by that regime (`AFFIX_STATUS`).
  3. A rarefied (transitional or free-molecular) follower starts at its windward surface point with the parent's
     velocity and is integrated through the shock layer (`adapters.rarefied_sphere_simulator`, the parent as an
     obstacle), whatever AFFIx says. A transitional follower is flagged (`transitional_follower`).
  4. A continuum follower follows AFFIx: far wake on the axis two parent diameters behind the centroid, escape at the
     windward surface point (both with the parent's velocity, integrated by `adapters.sphere_simulator` with no shock
     layer), shock band and trapped not integrated (reasons `shock_band` and `near_wake_trapped`).

Integrated followers run to the smaller of `ScenarioSettings.t_max` and the end of the parent's track (a release that
leaves less than MIN_FLIGHT_S of that track is refused), and `overlap.residence` measures them against
`Plume(parent track, law, bound)` from the release time. The far-wake start does not depend on theta, so the far-wake
integration and its residence are memoised in this process by (twin, diameter, altitude, settings) and shared by every
far-wake theta; each theta still has its own run name and writes its own files, since its AFFIx record differs.

`droplet_in_sphere_plume(sphere_d_mm, h_km, droplet_d_um)` is the worked example of spec 4.0: a solid sphere released
from the twin at the far-wake start (integrated with no wake gas) becomes the parent, and a molten droplet released at
the sphere's own release state (`release="centroid"`: the droplet-from-sphere release geometry is not decided) is the
follower, with the droplet model's own heating and evaporation (`DropletSimulator`, not modified).

Clocks. Times on the parent's clock (the twin's t = 0 at the top of DRAMA's descending branch): `t_release`,
`follower_t_end`, and the history's `t`. Times measured from release: `residence_time`, `t_first_entry`,
`first_interval`, `total_inside`, `t_shock_layer_exit` and the history's `t_since_release`.

Output: `CaseResult` holds the run name (`run_name`, which encodes every setting), the configuration, the
`ResidenceResult` (None when nothing was integrated), the summary dict and the history rows; `write_case` writes
`<out_dir>/runs/<name>.json` (the summary, NaN as null) and `<name>.csv` (the history, when there is one).
"""
from collections import OrderedDict
from dataclasses import asdict, dataclass, fields
import csv
import functools
import json
import math
import os
import re
import time
from datetime import timedelta

import numpy as np

from reentry_model.shared import aero, body, earth
from reentry_model.shared.atmosphere import US76TableAtmosphere
from reentry_model.droplet_reentry.model import T_RELEASE, VOLATILE, WE_CRIT, DropletSimulator, weber
from reentry_model.droplet_reentry.particle_trajectory import drag_coefficient as composed_drag_coefficient
from reentry_model.shared.gas import EquilibriumAir
from reentry_model.wake_column.run import _clean, _git_commit, _git_dirty
from reentry_model.wake_residence.adapters import (
    WakeAtmosphere, dense_solve, parent_collision_event, rarefied_sphere_simulator, shock_layer_exit,
    solution_track, sphere_parent_track, sphere_simulator, state_from_ecef, twin_parent)
from reentry_model.wake_residence.affix import (
    AffixSettings, affix_regime, far_wake_offset, release_offset)
from reentry_model.wake_residence.laws import (
    LawSettings, classify, knudsen, plume_radius, pressure_parameter, reynolds)
from reentry_model.wake_residence.overlap import Plume, residence
from reentry_model.wake_residence.shock_layer import (
    FreeStream, ShockLayerField, ShockLayerSettings)
from reentry_model.wake_residence.wake_gas import check_mode

AFFIX_STATUS = {"continuum": "governing", "transitional": "indicative", "free_molecular": "not_applicable"}
RELEASES = ("far_wake", "escape", "surface", "shock_band", "trapped", "no_law", "centroid")
NO_THETA_RELEASES = ("centroid",)             # releases whose run name carries no theta (theta is undefined there)
NO_LAW = {"rarefied": "rarefied_no_law", "too_slow": "too_slow_no_law"}
IN_LAYER_WEIGHT = 0.5         # the follower is "in the shock layer" where the field's weight exceeds this
LAYER_SAMPLE_DT = 1e-4        # s, sampling of the dense run up to the shock-layer exit for the smallest Kn_f ...
LAYER_MAX_SAMPLES = 20_000    # ... widened so that this grid never has more samples
MIN_FLIGHT_S = 1.0            # s, least flight time the parent's track must leave after a release
FAR_WAKE_MEMO_SIZE = 8        # far-wake integrations kept in memory (each serves every far-wake theta of its case)
HISTORY_FINE_DT = 0.02        # s, history spacing over the first HISTORY_FINE_S after release ...
HISTORY_FINE_S = 10.0
HISTORY_LOG_POINTS = 200      # ... plus this many points spaced geometrically from HISTORY_LOG_START to the end
HISTORY_LOG_START = 1e-4      # s
PARENT_AT_RELEASE_KEYS = ("regime_at_release", "laws_at_release", "kn_at_release", "re_at_release", "p_at_release",
                          "parent_speed_release")
DROPLET_LIMITATION = ("the plume laws are metre-scale wake results applied to a centimetre-scale sphere parent, and "
                      "the droplet's release geometry from the sphere is not decided (released at its centroid)")


@dataclass(frozen=True)
class ScenarioSettings:
    """Every choice a case depends on; each non-default value is encoded in the run name."""
    law: LawSettings = LawSettings()
    bound: str = "low"                 # "low" or "high": the smaller or larger radius where several laws apply
    wake_flow: str = "off"             # the wake gas the follower feels: "off", "velocity", "density" or "both"
    density: float = 2813.0            # kg/m3, the follower sphere (AA7075 solid)
    frontal_area: float = 6.215        # m2, the parent's (sets D1 and r0; the twin must be built with it)
    beta_scale: float = 1.0            # scale on the parent's ballistic-coefficient table (the twin's)
    k_E: float = 1.0                   # scale on the line-blast radius of the wake's energy profile
    t_max: float = 600.0               # s, longest follower integration (also capped by the parent's track)
    scan_dt: object = None             # s, finest residence scan step (None: overlap.residence's default)
    affix: AffixSettings = AffixSettings()
    kn_free_molecular: float = 10.0    # follower Kn_f at or above which it is free-molecular (spec 4.3c)
    shock: ShockLayerSettings = ShockLayerSettings()

    def __post_init__(self):
        if self.bound not in ("low", "high"):
            raise ValueError(f"bound must be 'low' or 'high', got {self.bound!r}")
        check_mode(self.wake_flow)
        for name in ("density", "frontal_area", "beta_scale", "k_E", "t_max", "kn_free_molecular"):
            v = getattr(self, name)
            if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0.0):
                raise ValueError(f"{name} must be a finite positive number, got {v!r}")
        if self.scan_dt is not None and not (math.isfinite(self.scan_dt) and self.scan_dt > 0.0):
            raise ValueError(f"scan_dt must be None or a finite positive number of seconds, got {self.scan_dt!r}")


@dataclass
class CaseResult:
    name: str
    config: dict
    residence: object                  # overlap.ResidenceResult, or None when nothing was integrated
    summary: dict
    history: list


def follower_regime(kn_f, kn_continuum, kn_free_molecular):
    """`continuum` below kn_continuum, `free_molecular` at or above kn_free_molecular, `transitional` between."""
    if not math.isfinite(kn_f) or kn_f < 0.0:
        raise ValueError(f"the follower's Knudsen number must be finite and non-negative, got {kn_f!r}")
    if kn_f < kn_continuum:
        return "continuum"
    return "free_molecular" if kn_f >= kn_free_molecular else "transitional"


# ---------------------------------------------------------------------------------------------------- run names

def _tag(v):
    """Text of a setting's value for a run name: floats by repr, sequences joined by '-', unsafe characters as '-'."""
    if isinstance(v, (tuple, list)):
        s = "-".join(_tag(x) for x in v)
    elif isinstance(v, float):
        s = repr(v)
    else:
        s = str(v)
    return re.sub(r"[^A-Za-z0-9.+-]", "-", s)


def _non_default(obj, default, skip=()):
    return {f.name: getattr(obj, f.name) for f in fields(obj)
            if f.name not in skip and getattr(obj, f.name) != getattr(default, f.name)}


def run_name(kind, d_mm, h_km, settings, droplet_d_um=None, release="far_wake", theta_deg=None):
    """The case's name: kind and diameters, release altitude, release (with theta, except for centroid releases; a
    far-wake name carries theta too, because its summary holds theta's AFFIx record, although its integration is
    shared by every far-wake theta), bound, wake_flow and turbulent law, then `_<field><value>` for every other non-default setting
    (law, scenario, `affix_`-prefixed AFFIx and shock-layer fields), sorted by field name."""
    if kind == "sphere":
        head = f"sphere_d{d_mm:06.1f}mm"
    elif kind == "droplet":
        if droplet_d_um is None:
            raise ValueError("a droplet run name needs droplet_d_um")
        head = f"droplet_s{d_mm:06.1f}mm_d{droplet_d_um:06.1f}um"
    else:
        raise ValueError(f"kind must be 'sphere' or 'droplet', got {kind!r}")
    if release not in RELEASES:
        raise ValueError(f"release must be one of {RELEASES}, got {release!r}")
    rel = release.replace("_", "")
    if release not in NO_THETA_RELEASES:
        if theta_deg is None:
            raise ValueError(f"a {release!r} run name needs theta_deg")
        rel += f"_th{theta_deg:05.1f}"
    extra = _non_default(settings.law, LawSettings(), skip=("turbulent_law",))
    extra.update(_non_default(settings, ScenarioSettings(), skip=("law", "bound", "wake_flow", "affix", "shock")))
    extra.update({f"affix_{k}": v for k, v in _non_default(settings.affix, AffixSettings()).items()})
    extra.update(settings.shock.non_default())
    suffix = "".join(f"_{k}{_tag(extra[k])}" for k in sorted(extra))
    return (f"{head}_h{h_km:05.1f}km_{rel}_{settings.bound}_{settings.wake_flow}_{settings.law.turbulent_law}"
            f"{suffix}")


# ------------------------------------------------------------------------------------- shared, expensive objects

@functools.lru_cache(maxsize=None)
def _base_atmosphere():
    return US76TableAtmosphere()


@functools.lru_cache(maxsize=None)
def _air():
    """One EquilibriumAir per process: building it takes a while."""
    return EquilibriumAir()


@functools.lru_cache(maxsize=None)
def _drag_tables():
    return aero.SphereDragTables.from_json(), aero.bridging_by_name("sesam-table")


@functools.lru_cache(maxsize=4)
def _built_twin(frontal_area, beta_scale):
    return twin_parent(frontal_area=frontal_area, beta_scale=beta_scale)


def _twin_for(twin, settings):
    """The twin to use: built for the settings when None, else checked against settings.frontal_area and
    settings.beta_scale (ValueError if it was built with others, so a run name never claims a parent it did not use)."""
    if twin is None:
        return _built_twin(settings.frontal_area, settings.beta_scale)
    tr = twin.track
    e = tr.emission_state(tr.t0)
    d = math.sqrt(4.0 * settings.frontal_area / math.pi)
    if not math.isclose(e.d, d, rel_tol=1e-9):
        raise ValueError(f"the twin's parent diameter is {e.d!r} m, but settings.frontal_area = "
                         f"{settings.frontal_area!r} m2 gives {d!r} m; pass a twin built with that frontal_area")
    cda = float(twin.beta.cda_at(twin.altitude(tr.t0))) / settings.beta_scale
    if not math.isclose(e.cda, cda, rel_tol=1e-6):
        raise ValueError(f"the twin's C_D A at its start is {e.cda!r} m2, but settings.beta_scale = "
                         f"{settings.beta_scale!r} gives {cda!r} m2; pass a twin built with that beta_scale")
    return twin


@dataclass
class _Release:
    t_r: float              # release time on the parent's clock, s
    P: np.ndarray           # parent's centroid, ECEF m
    V: np.ndarray           # parent's velocity, ECEF m/s
    v_hat: np.ndarray
    up: np.ndarray          # local geodetic up at P
    h: float                # parent's geodetic altitude, m
    epoch: object           # epoch of the release
    e: object               # parent's EmissionState at release
    m_bar: float            # free-stream mean molecular mass, kg


def _release(twin, h_km):
    t_r = twin.t_at_altitude(h_km)
    P, V = twin.track.position(t_r), twin.track.velocity(t_r)
    h, lat, lon = earth.ecef_to_geodetic(P)
    up = earth.enu_basis(lat, lon)[2]
    m_bar = _base_atmosphere().state(0.0, h, lat, lon).m_bar
    return _Release(t_r, P, V, V / np.linalg.norm(V), up, h, twin.epoch + timedelta(seconds=t_r),
                    twin.track.emission_state(t_r), m_bar)


def _check_flight_time(twin, rel, h_km):
    """ValueError when a release at h_km leaves the follower less than MIN_FLIGHT_S of the parent's track."""
    left = twin.track.t1 - rel.t_r
    if left < MIN_FLIGHT_S:
        raise ValueError(f"a release at {h_km!r} km (t = {rel.t_r:.3f} s) leaves {left:.3f} s of the twin's track, which "
                         f"ends at t = {twin.track.t1:.3f} s ({twin.altitude(twin.track.t1) / 1e3:.3f} km); at least "
                         f"{MIN_FLIGHT_S} s is needed")


def classify_release(d_mm, h_km, theta_deg, twin, settings):
    """AFFIx outcome (`affix.AffixResult`) of a d_mm sphere released at theta from the twin when it reaches h_km: the
    parent's EmissionState at the release time, the sphere's density, the parent's mass there (`twin.beta.mass_at`),
    and D1 the parent's equal-area diameter."""
    return _affix(d_mm * 1e-3, theta_deg, _release(twin, h_km), twin, settings)


def _affix(d, theta_deg, rel, twin, settings):
    return affix_regime(d, rel.e.d, theta_deg, rel.e, settings.density,
                        parent_mass=float(twin.beta.mass_at(rel.h)), settings=settings.affix)


# --------------------------------------------------------------------------------------------- summary pieces

def _provenance(sha256):
    return {"parent_sha256": sha256, "git_commit": _git_commit(), "git_dirty": _git_dirty()}


def _parent_at_release(e, settings):
    """The parent's plume regime, laws and flow numbers at the release emission point."""
    rr = plume_radius(0.0, e, settings.law, settings.bound)
    return {"regime_at_release": classify(e, settings.law), "laws_at_release": list(rr.laws),
            "kn_at_release": knudsen(e), "re_at_release": reynolds(e), "p_at_release": pressure_parameter(e),
            "parent_speed_release": e.U}


def _not_integrated(reason, censored=False):
    return {"reason": reason, "residence_reason": None, "residence_time": math.nan, "censored": censored,
            "total_inside": math.nan, "n_intervals": 0, "t_first_entry": math.nan, "first_interval": math.nan,
            "follower_t_end": math.nan, "follower_end_reason": None, "interp_error_m": math.nan}


def _air_speed_and_kn(wa, t, state, d):
    """(speed relative to the local gas [m/s], local Kn_f) of a follower at follower time t, state (6,) ECEF, from
    `wa.state` (wake, shock-layer or ambient gas)."""
    h, lat, lon = earth.ecef_to_geodetic(state[:3])
    fs = wa.state(float(t), max(h, 0.0), lat, lon)
    e, n, u = earth.enu_basis(lat, lon)
    wind = fs.wind_enu[0] * e + fs.wind_enu[1] * n + fs.wind_enu[2] * u
    return float(np.linalg.norm(state[3:6] - wind)), aero.mean_free_path(fs.rho, fs.m_bar) / d


def _residence_summary(res, run, wa, plume, t_r, d):
    """Residence numbers, the first exit's state and the local Kn_f at that exit (or at the end) and at the end."""
    out = {"residence_reason": res.reason, "residence_time": res.residence_time, "censored": res.censored,
           "total_inside": res.total_inside, "n_intervals": len(res.intervals),
           "follower_t_end": t_r + run.t_end, "follower_end_reason": run.end_reason,
           "interp_error_m": res.interp_error_m}
    if res.intervals:
        start, end = res.intervals[0]
        out["t_first_entry"] = start - t_r
        out["first_interval"] = (t_r + run.t_end if end is None else end) - start
    else:
        out["t_first_entry"] = out["first_interval"] = math.nan
    flags = {}
    for v in res.validity:
        for k, raised in v.flags.items():
            flags[k] = flags.get(k, False) or bool(raised)
    out.update({f"plume_{k}": v for k, v in flags.items()})
    keys = ("exit_altitude_m", "exit_speed", "exit_lag_m", "exit_distance_m", "exit_r_b", "exit_s", "exit_through",
            "regime_at_exit", "laws_at_exit")
    out.update(dict.fromkeys(keys))
    t_exit = run.t_end
    if res.exits:
        ex = res.exits[0]
        t_exit = ex.t - t_r
        m = plume.margin(ex.position, ex.t, t_release=t_r)
        out.update({"exit_altitude_m": earth.ecef_to_geodetic(ex.position)[0], "exit_lag_m": ex.lag_m,
                    "exit_distance_m": ex.distance_m, "exit_r_b": ex.r_b, "exit_s": ex.s_behind,
                    "exit_through": ex.through, "regime_at_exit": m.regime, "laws_at_exit": list(m.laws)})
    speed, out["kn_f_exit"] = _air_speed_and_kn(wa, t_exit, run.sol.sol(t_exit), d)
    if res.exits:
        out["exit_speed"] = speed
    out["kn_f_end"] = _air_speed_and_kn(wa, run.t_end, run.sol.sol(run.t_end), d)[1]
    return out


def _history(plume, follower, t_r, extra_times=()):
    """Rows of the follower and the parent at the same times: every HISTORY_FINE_DT over the first HISTORY_FINE_S,
    HISTORY_LOG_POINTS more spaced geometrically to the end, and `extra_times` (the exits), all from release."""
    span = follower.t1 - t_r
    ts = [np.arange(0.0, min(span, HISTORY_FINE_S), HISTORY_FINE_DT), [span], np.asarray(extra_times) - t_r]
    if span > HISTORY_LOG_START:
        ts.append(np.geomspace(HISTORY_LOG_START, span, HISTORY_LOG_POINTS))
    ts = np.unique(np.clip(np.concatenate(ts), 0.0, span))
    rows = []
    for dt in ts:
        t = t_r + float(dt)
        x, P, v = follower.position(t), plume.parent.position(t), plume.parent.velocity(t)
        m = plume.margin(x, t, t_release=t_r)
        rows.append({"t": t, "t_since_release": float(dt),
                     "follower_x": x[0], "follower_y": x[1], "follower_z": x[2],
                     "parent_x": P[0], "parent_y": P[1], "parent_z": P[2],
                     "lag_m": float(np.dot(P - x, v / np.linalg.norm(v))), "distance_m": m.distance,
                     "r_b": m.r_b, "margin": m.margin, "regime": m.regime})
    return rows


def _layer_sample_times(step_times, t_until):
    """The integrator's step times and a grid from 0 to t_until every LAYER_SAMPLE_DT, widened so the grid has at most
    LAYER_MAX_SAMPLES points (a follower held in the layer for a long time would otherwise give millions)."""
    n = min(int(math.ceil(t_until / LAYER_SAMPLE_DT)), LAYER_MAX_SAMPLES - 1) + 1 if t_until > 0.0 else 1
    return np.unique(np.concatenate([np.asarray(step_times, dtype=float), np.linspace(0.0, max(t_until, 0.0), n)]))


def _residence_exceeds_t_fz(res, t_fz):
    """True when the follower spends more than 3 t_fz inside the plume in all (spec 4.3b): `total_inside`, because an
    escape release starts outside the plume, so its first-interval `residence_time` is 0."""
    return bool(res.total_inside > 3.0 * t_fz)


def _shock_layer_summary(run, sim, wa, field, ex, twin, rel, d, settings):
    """The shock-layer record of a surface release (spec 4.3c-d)."""
    hit = ex.kind == "hit_parent"
    hit_theta = math.nan
    if hit:
        x, y, _, _ = wa.parent_frame(run.t_end, run.sol.sol(run.t_end)[:3])
        hit_theta = math.degrees(math.atan2(y, -x))
    # smallest local Kn_f where the field's weight exceeds IN_LAYER_WEIGHT, on the integrator's steps and a fine grid
    # up to the shock-layer exit
    ts = _layer_sample_times(run.sol.t, min(ex.t, run.t_end))
    states = np.atleast_2d(run.sol.sol(ts).T)
    kn_min = math.inf
    for t, st in zip(ts, states):
        xx, yy, _, _ = wa.parent_frame(t, st[:3])
        if field.state(xx, yy).weight > IN_LAYER_WEIGHT:
            kn_min = min(kn_min, _air_speed_and_kn(wa, t, st, d)[1])
    kn_min = kn_min if math.isfinite(kn_min) else math.nan
    q0 = 0.5 * rel.e.rho * rel.e.U ** 2
    e1 = twin.track.emission_state(rel.t_r + ex.t)
    q_change = 0.5 * e1.rho * e1.U ** 2 / q0 - 1.0
    fs = field.fs
    tables, bridging = _drag_tables()
    out = {"shock_layer_exit": ex.kind, "t_shock_layer_exit": ex.t, "x_exit_parent_frame": ex.x,
           "y_exit_parent_frame": ex.y, "v_rel_shock_layer_exit": ex.v_rel, "hit_theta_deg": hit_theta,
           "kn_f_min_in_layer": kn_min, "flux_mismatch": field.flux_mismatch, "flux_iterations": field.flux_iterations,
           "epsilon": field.eps, "standoff_m": field.standoff, "standoff_billig_m": field.standoff_billig,
           "mach": field.mach, "a_p_max": sim.max_a_pressure, "a_drag_max_shock_layer": sim.max_a_drag,
           "q_change": q_change,
           "cd_composed_release": composed_drag_coefficient(fs.U, fs.rho, fs.T, fs.p, d,
                                                            temperature_ratio=settings.shock.follower_wall_T_K
                                                            / fs.T)[0],
           "cd_sesam_release": aero.drag_coefficient(aero.knudsen(fs.rho, fs.m_bar, d),
                                                     aero.mach(fs.U, fs.T, fs.m_bar), tables, bridging),
           "field_build_time_s": field.build_time_s,
           "follower_continuum_in_layer": bool(kn_min < settings.affix.kn_continuum),
           "shock_layer_nonequilibrium": field.flags["kn_parent_above_flag"],
           "q_change_large": abs(q_change) > settings.shock.q_change_flag}
    out.update(field.flags)
    return out


# -------------------------------------------------------------------------------------------------- scenarios

# Far-wake integrations by (id of the twin, d_mm, h_km, settings): (twin, dense run, follower Track, ResidenceResult).
# The far-wake start does not depend on theta (spec 4.3b), so every far-wake theta of a case reuses one integration and
# writes its own files. The twin is kept in the entry so a recycled id never matches another twin.
_FAR_WAKE_RUNS = OrderedDict()


def _case(name, config, res, summary, history, started):
    summary["name"] = name
    summary["wall_time_s"] = time.perf_counter() - started
    return CaseResult(name, config, res, summary, history)


def sphere_in_spacecraft_plume(d_mm, h_km, theta_deg=90.0, settings=ScenarioSettings(), twin=None):
    """A d_mm solid sphere released at theta_deg (90 to 180 degrees) from the Generic Spacecraft's twin when it reaches
    h_km (module docstring). `twin` defaults to one built for settings.frontal_area and settings.beta_scale."""
    if not 90.0 <= theta_deg <= 180.0:
        raise ValueError(f"spacecraft releases are windward only: theta_deg must lie in [90, 180], got {theta_deg!r}")
    started = time.perf_counter()
    twin = _twin_for(twin, settings)
    d = d_mm * 1e-3
    rel = _release(twin, h_km)
    _check_flight_time(twin, rel, h_km)
    e, D1 = rel.e, rel.e.d
    base = _base_atmosphere()
    plume = Plume(twin.track, settings.law, settings.bound)
    tables, bridging = _drag_tables()
    cd_sesam = aero.drag_coefficient(aero.knudsen(e.rho, rel.m_bar, d), aero.mach(e.U, e.T, rel.m_bar), tables,
                                     bridging)
    summary = {"kind": "sphere", "diameter_mm": d_mm, "altitude_km": h_km, "theta_deg": theta_deg,
               "t_release": rel.t_r, "bound": settings.bound, "wake_flow": settings.wake_flow,
               "beta_sphere_release": body.sphere_mass(d, settings.density) / (cd_sesam * math.pi * d * d / 4.0),
               "beta_parent_release": float(twin.beta.mass_at(rel.h)) / e.cda,
               "parent_interp_error_m": twin.track.interpolation_error(),
               "kn_f_free_stream": aero.mean_free_path(e.rho, rel.m_bar) / d,
               "transitional_follower": False, **_parent_at_release(e, settings), **_provenance(twin.sha256)}
    config = {"kind": "sphere", "d_mm": d_mm, "h_km": h_km, "theta_deg": theta_deg, "settings": asdict(settings)}

    def done(release, res, history):
        name = run_name("sphere", d_mm, h_km, settings, release=release, theta_deg=theta_deg)
        return _case(name, config, res, summary, history, started)

    # 1. the parent's regime
    if summary["regime_at_release"] in NO_LAW:
        summary.update(_not_integrated(NO_LAW[summary["regime_at_release"]]), release_kind=None, affix_zone=None,
                       follower_regime=None, affix_status=None, kn_f_release=math.nan)
        return done("no_law", None, [])

    # 2. the shock layer, the follower's regime and AFFIx
    fs = FreeStream(rho=e.rho, T=e.T, p=e.p, U=e.U, mfp=e.mfp, m_bar=rel.m_bar)
    field = ShockLayerField.build(fs, 0.5 * D1, e.r0, settings.shock, air=_air())
    surf = field.surface_state(theta_deg)
    kn_f = aero.mean_free_path(surf.rho, surf.m_bar) / d
    regime = follower_regime(kn_f, settings.affix.kn_continuum, settings.kn_free_molecular)
    af = _affix(d, theta_deg, rel, twin, settings)
    summary.update({"kn_f_release": kn_f, "follower_regime": regime, "affix_status": AFFIX_STATUS[regime],
                    "affix_zone": af.zone, "affix_r": af.r, "theta_inn": af.theta_inn, "theta_bif": af.theta_bif,
                    "t_hat": af.t_hat, "t_fz": af.t_fz, "transitional_follower": regime == "transitional",
                    "residence_exceeds_t_fz": False, **{f"affix_{k}": v for k, v in af.flags.items()},
                    "surface_match_phi_deg": field.surface_match_phi_deg,
                    "shoulder_pressure_ratio": field.shoulder_pressure_ratio,
                    "surface_npm_no_match": field.flags["npm_no_match"],
                    "surface_below_p_inf": field.flags["surface_below_p_inf"]})

    cached = None                  # the far-wake memo entry this case reuses, if any
    if regime != "continuum":
        # 3. rarefied follower: the surface point, through the shock layer
        release = "surface"
        wa = WakeAtmosphere(base, plume, rel.t_r, settings.wake_flow, k_E=settings.k_E, shock_layer=field)
        t_max = min(settings.t_max, wa.follower_t_max)
        initial = state_from_ecef(rel.P + release_offset(theta_deg, D1, rel.v_hat, rel.up), rel.V, rel.epoch)
        sim = rarefied_sphere_simulator(initial, d, wa, density=settings.density,
                                        wall_T_K=settings.shock.follower_wall_T_K,
                                        pressure_gradient=settings.shock.pressure_gradient, t_max=t_max)
        run = dense_solve(sim, t_max, events=sim._events() + [parent_collision_event(wa, 0.5 * D1)])
        ex = shock_layer_exit(run, wa, field, hit=run.end_reason == "hit_parent")
        summary.update(_shock_layer_summary(run, sim, wa, field, ex, twin, rel, d, settings))
    elif af.zone in ("far_wake", "escape"):
        # 4. continuum follower, AFFIx governs
        release = af.zone
        wa = WakeAtmosphere(base, plume, rel.t_r, settings.wake_flow, k_E=settings.k_E)
        key = (id(twin), d_mm, h_km, settings)
        cached = _FAR_WAKE_RUNS.get(key) if release == "far_wake" else None
        if cached is not None and cached[0] is not twin:
            cached = None
        if cached is not None:
            _FAR_WAKE_RUNS.move_to_end(key)
            run, follower, res = cached[1:]
        else:
            offset = (far_wake_offset(D1, rel.v_hat) if release == "far_wake"
                      else release_offset(theta_deg, D1, rel.v_hat, rel.up))
            t_max = min(settings.t_max, wa.follower_t_max)
            sim = sphere_simulator(state_from_ecef(rel.P + offset, rel.V, rel.epoch), d, wa,
                                   density=settings.density, t_max=t_max)
            run = dense_solve(sim, t_max)
    else:
        summary["release_kind"] = None
        if af.zone == "shock_band":
            summary.update(_not_integrated("shock_band"))
        else:
            summary.update(_not_integrated("near_wake_trapped", censored=True))
        return done(af.zone, None, [])

    if cached is None:
        follower = solution_track(run, rel.t_r)
        res = residence(plume, follower, t_start=rel.t_r, t_end=follower.t1, scan_dt=settings.scan_dt)
        if release == "far_wake":
            _FAR_WAKE_RUNS[key] = (twin, run, follower, res)
            while len(_FAR_WAKE_RUNS) > FAR_WAKE_MEMO_SIZE:
                _FAR_WAKE_RUNS.popitem(last=False)
    summary["release_kind"] = release
    summary.update(_residence_summary(res, run, wa, plume, rel.t_r, d))
    summary["reason"] = "hit_parent" if run.end_reason == "hit_parent" else res.reason
    if release == "escape":
        summary["residence_exceeds_t_fz"] = _residence_exceeds_t_fz(res, af.t_fz)
    return done(release, res, _history(plume, follower, rel.t_r, [x.t for x in res.exits]))


def _droplet_events(sim):
    """The two terminal events of `DropletSimulator.run`: the ground, and the droplet's mass falling to VAPORISED of
    its initial mass (named `vaporised`, the end reason `run` reports for it)."""
    s = sim.settings

    def ground(tt, yy):
        return earth.ecef_to_geodetic(yy[:3])[0] - s.ground_altitude
    ground.terminal, ground.direction = True, -1

    def vaporised(tt, yy):
        return float(np.maximum(yy[7:], 0.0).sum()) + sim.m_inert - sim.VAPORISED * sim.m0
    vaporised.terminal, vaporised.direction = True, -1
    return [ground, vaporised]


def droplet_in_sphere_plume(sphere_d_mm, h_km, droplet_d_um, settings=ScenarioSettings(), twin=None):
    """A droplet_d_um molten droplet in the plume of a sphere_d_mm solid sphere released from the twin at h_km (module
    docstring). The sphere is released at the far-wake start with the parent's velocity and integrated with no wake
    gas; the droplet starts at the sphere's release state. A droplet whose release Weber number is at least WE_CRIT is
    not run (reason `skipped_weber`, as `run_case` does); nor is one whose sphere parent has no plume law at release."""
    started = time.perf_counter()
    twin = _twin_for(twin, settings)
    rel = _release(twin, h_km)
    _check_flight_time(twin, rel, h_km)
    base = _base_atmosphere()
    d_s, d0 = sphere_d_mm * 1e-3, droplet_d_um * 1e-6
    initial = state_from_ecef(rel.P + far_wake_offset(rel.e.d, rel.v_hat), rel.V, rel.epoch)
    fs = base.state(0.0, initial.altitude, initial.lat, initial.lon)
    we0 = weber(fs.rho, initial.velocity, d0, T_RELEASE)        # known before the sphere is integrated
    summary = {"kind": "droplet", "sphere_diameter_mm": sphere_d_mm, "droplet_diameter_um": droplet_d_um,
               "altitude_km": h_km, "release": "centroid", "t_release": rel.t_r, "bound": settings.bound,
               "wake_flow": settings.wake_flow, "weber_release": we0,
               "beta_parent_release": float(twin.beta.mass_at(rel.h)) / rel.e.cda,
               "kn_f_free_stream": aero.mean_free_path(fs.rho, fs.m_bar) / d0, "limitation": DROPLET_LIMITATION,
               **_provenance(twin.sha256)}
    config = {"kind": "droplet", "sphere_d_mm": sphere_d_mm, "h_km": h_km, "droplet_d_um": droplet_d_um,
              "release": "centroid", "settings": asdict(settings)}
    name = run_name("droplet", sphere_d_mm, h_km, settings, droplet_d_um=droplet_d_um, release="centroid")
    if we0 >= WE_CRIT:
        # the sphere is not integrated either: its keys are kept, empty, so every droplet summary has the same keys
        summary.update({"sphere_t_end": math.nan, "sphere_end_reason": None, "parent_interp_error_m": math.nan,
                        "beta_sphere_release": math.nan, **dict.fromkeys(PARENT_AT_RELEASE_KEYS)},
                       **_not_integrated("skipped_weber"))
        return _case(name, config, None, summary, [], started)

    wa_s = WakeAtmosphere(base, Plume(twin.track, settings.law, settings.bound), rel.t_r, "off")
    t_max_s = min(settings.t_max, wa_s.follower_t_max)
    sim_s = sphere_simulator(initial, d_s, wa_s, density=settings.density, t_max=t_max_s)
    run_s = dense_solve(sim_s, t_max_s)
    sphere = sphere_parent_track(run_s, sim_s, rel.t_r, base)
    plume = Plume(sphere, settings.law, settings.bound)
    e = sphere.emission_state(rel.t_r)
    summary.update({"sphere_t_end": rel.t_r + run_s.t_end, "sphere_end_reason": run_s.end_reason,
                    "parent_interp_error_m": sphere.interpolation_error(),
                    "beta_sphere_release": body.sphere_mass(d_s, settings.density) / e.cda,
                    **_parent_at_release(e, settings)})
    if summary["regime_at_release"] in NO_LAW:
        summary.update(_not_integrated(NO_LAW[summary["regime_at_release"]]))
        return _case(name, config, None, summary, [], started)

    wa = WakeAtmosphere(base, plume, rel.t_r, settings.wake_flow, k_E=settings.k_E)
    sim = DropletSimulator(initial, droplet_d_um, atmosphere=wa)
    s = sim.settings
    atol = [s.atol_position] * 3 + [s.atol_velocity] * 3 + [1e-6] + [1e-9 * sim.m0] * len(VOLATILE)
    t_max = min(settings.t_max, wa.follower_t_max)
    run = dense_solve(sim, t_max, events=_droplet_events(sim), method="LSODA", rtol=1e-8, atol=atol)
    follower = solution_track(run, rel.t_r)
    res = residence(plume, follower, t_start=rel.t_r, t_end=follower.t1, scan_dt=settings.scan_dt)
    summary.update(_residence_summary(res, run, wa, plume, rel.t_r, d0))
    summary["reason"] = res.reason
    return _case(name, config, res, summary, _history(plume, follower, rel.t_r, [x.t for x in res.exits]), started)


# ------------------------------------------------------------------------------------------------------- files

def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (set, frozenset)):
        return sorted(_jsonable(v) for v in x)
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    return _clean(x)


def write_case(case, out_dir):
    """Write `<out_dir>/runs/<name>.json` (the summary, with the configuration under "config"; NaN as null) and, when
    the case was integrated, `<name>.csv` (the history). Returns (csv_path or None, json_path)."""
    runs = os.path.join(os.fspath(out_dir), "runs")
    os.makedirs(runs, exist_ok=True)
    csv_path = None
    if case.residence is not None and case.history:
        csv_path = os.path.join(runs, case.name + ".csv")
        with open(csv_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(case.history[0]))
            w.writeheader()
            for row in case.history:
                w.writerow({k: _clean(v) for k, v in row.items()})
    json_path = os.path.join(runs, case.name + ".json")
    tmp = json_path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(_jsonable({**case.summary, "config": case.config}), fh, indent=2, allow_nan=False)
        fh.write("\n")
    os.replace(tmp, json_path)
    return csv_path, json_path
