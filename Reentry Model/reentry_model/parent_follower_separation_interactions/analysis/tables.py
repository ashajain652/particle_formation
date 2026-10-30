"""The supporting tables (spec 5, Task 9), all from the twin and the laws, with no follower integrated:

    "$PY" -m reentry_model.parent_follower_separation_interactions.analysis.tables [--out DIR]

  plume_size.csv           per release altitude: the parent's plume regime, Knudsen number, Reynolds number on its
                           diameter, Taylor's pressure parameter, and r_b (low and high bound) at s = 10 ... 30,000 m
                           behind it, with the time since passage s / U and a flag for s beyond Park's distance limit.
  ballistic_coefficient.csv  the twin's ballistic coefficient and C_D A at each altitude, and each default sphere's, so
                           the direction of the lag is visible (a sphere with the larger coefficient falls behind less).
  twin_vs_drama.csv        `adapters.twin_vs_drama` plus the twin's downrange at the DRAMA time and its difference from
                           DRAMA's downrange since the top of the descending branch.
  affix_regimes.csv        the AFFIx classification of every (diameter, altitude, theta): zone, angles, times, flags.
  drag_ratio_check.csv     per altitude and bound: rho_w (V - u)^2 / (rho_inf V^2) at the far-wake start point from
                           the wake gas model (`wake_flow` "both"), next to Register et al.'s 0.14.
  wake_gas_sanity.csv      per altitude and bound: the model's core gas speed 30 km behind the parent and core
                           temperature 300 m behind, next to the values spec 4.3a assumed (about 10 m/s, 1,000-2,000 K).

Rows with no plume law (a rarefied parent: 100 km in the default set) have blank radii and gas values.
"""
import argparse
import math
import os
import sys

import numpy as np

from reentry_model import aero, body, earth
from reentry_model.droplet_reentry import parent_trajectory
from reentry_model.parent_follower_separation_interactions import OUTPUT
from reentry_model.parent_follower_separation_interactions.adapters import twin_vs_drama
from reentry_model.parent_follower_separation_interactions.affix import far_wake_offset
from reentry_model.parent_follower_separation_interactions.analysis.sweep import (
    DEFAULT_ALTITUDES, DEFAULT_BOUNDS, DEFAULT_DIAMETERS, DEFAULT_THETAS, get_twin, kept_altitudes, write_csv)
from reentry_model.parent_follower_separation_interactions.laws import (
    LawSettings, classify, knudsen, plume_radius, pressure_parameter, reynolds)
from reentry_model.parent_follower_separation_interactions.overlap import Plume
from reentry_model.parent_follower_separation_interactions.scenarios import (
    ScenarioSettings, _base_atmosphere, _drag_tables, _release, classify_release)
from reentry_model.parent_follower_separation_interactions.wake_gas import Ambient

PLUME_DISTANCES_M = (10.0, 100.0, 1000.0, 3000.0, 10000.0, 30000.0)
REGISTER_DRAG_RATIO = 0.14           # Register et al. 2020: the follower's drag in the wake relative to free stream
CORE_SPEED_DISTANCE_M = 30e3         # m behind the parent, where the core gas speed is read
CORE_TEMPERATURE_DISTANCE_M = 300.0  # m behind the parent, where the core temperature is read
ASSUMED_CORE_SPEED_MS = 10.0         # m/s, spec 4.3a
ASSUMED_CORE_T_K = "1000-2000"       # K, spec 4.3a


def plume_size_rows(twin, altitudes, law=LawSettings()):
    rows = []
    for h in altitudes:
        e = _release(twin, h).e
        regime = classify(e, law)
        for s in PLUME_DISTANCES_M:
            lo, hi = (plume_radius(s, e, law, b) for b in DEFAULT_BOUNDS)
            rows.append({"altitude_km": h, "regime": regime, "kn": knudsen(e), "re_d": reynolds(e),
                         "p_cmhg_cm": pressure_parameter(e), "parent_speed_ms": e.U, "s_m": s,
                         "t_since_passage_s": s / e.U, "r_b_low_m": lo.r, "r_b_high_m": hi.r,
                         "laws_low": list(lo.laws), "laws_high": list(hi.laws),
                         "beyond_s_park_limit": s > law.s_park_limit,
                         **{f"flag_{k}": v for k, v in lo.flags.items()}})
    return rows


def ballistic_rows(twin, altitudes, diameters, density=ScenarioSettings().density):
    tables, bridging = _drag_tables()
    rows = []
    for h in altitudes:
        rel = _release(twin, h)
        e = rel.e
        mass = float(twin.beta.mass_at(rel.h))
        beta_parent = mass / e.cda
        for d_mm in diameters:
            d = d_mm * 1e-3
            cd = aero.drag_coefficient(aero.knudsen(e.rho, rel.m_bar, d), aero.mach(e.U, e.T, rel.m_bar), tables,
                                       bridging)
            beta = body.sphere_mass(d, density) / (cd * math.pi * d * d / 4.0)
            rows.append({"altitude_km": h, "parent_mass_kg": mass, "parent_cda_m2": e.cda,
                         "parent_beta_kg_m2": beta_parent, "diameter_mm": d_mm, "sphere_cd": cd,
                         "sphere_beta_kg_m2": beta, "sphere_over_parent_beta": beta / beta_parent,
                         "sphere_falls_behind_less": beta > beta_parent})
    return rows


def drama_rows(twin, parent):
    """`twin_vs_drama` rows with the twin's downrange (great-circle distance from the twin's start point) at the DRAMA
    time, and its difference from DRAMA's downrange measured from the top of the descending branch (DRAMA's own downrange
    counts from the start of the whole flight, 18,185 km at the top)."""
    offset = (twin.epoch - parent.epoch).total_seconds()
    _, lat0, lon0 = earth.ecef_to_geodetic(twin.track.position(twin.track.t0))
    top = parent.branch[0]["downrange"]
    rows = []
    for row, r in zip(twin_vs_drama(twin, parent), parent.branch):
        tt = r["time"] - offset
        if twin.track.t0 <= tt <= twin.track.t1:
            _, lat, lon = earth.ecef_to_geodetic(twin.track.position(tt))
            down = earth.great_circle_distance(lat0, lon0, lat, lon) / 1e3
        else:
            down = math.nan
        rows.append({**row, "drama_downrange_from_top_km": r["downrange"] - top, "twin_downrange_km": down,
                     "downrange_diff_km": down - (r["downrange"] - top)})
    return rows


def affix_rows(twin, diameters, altitudes, thetas, settings=ScenarioSettings()):
    rows = []
    for h in altitudes:
        for d_mm in diameters:
            for th in thetas:
                a = classify_release(d_mm, h, th, twin, settings)
                rows.append({"diameter_mm": d_mm, "altitude_km": h, "theta_deg": th, "zone": a.zone, "r": a.r,
                             "theta_inn": a.theta_inn, "theta_bif": a.theta_bif, "t_hat": a.t_hat, "t_fz": a.t_fz,
                             **{f"flag_{k}": v for k, v in a.flags.items()}})
    return rows


def _wake_gas(twin, rel, bound, x, mode="both"):
    """(GasState or None, MarginResult) of the wake gas at ECEF x at the release time on the parent's clock."""
    h, lat, lon = earth.ecef_to_geodetic(x)
    fs = _base_atmosphere().state(0.0, max(h, 0.0), lat, lon)
    plume = Plume(twin.track, LawSettings(), bound)
    gas, m = plume.wake_state(x, rel.t_r, Ambient(fs.p, fs.T, fs.rho), mode)
    return (None if not math.isfinite(m.margin) else gas), m, fs


def drag_ratio_rows(twin, altitudes):
    rows = []
    for h in altitudes:
        rel = _release(twin, h)
        x = rel.P + far_wake_offset(rel.e.d, rel.v_hat)
        for b in DEFAULT_BOUNDS:
            gas, m, fs = _wake_gas(twin, rel, b, x)
            row = {"altitude_km": h, "bound": b, "register_drag_ratio": REGISTER_DRAG_RATIO,
                   "parent_speed_ms": rel.e.U, "regime": m.regime}
            if gas is not None:
                u = float(np.linalg.norm(gas.u))
                ratio = gas.rho * (rel.e.U - u) ** 2 / (fs.rho * rel.e.U ** 2)
                row.update({"drag_ratio": ratio, "wake_density_over_ambient": gas.rho / fs.rho,
                            "wake_gas_speed_ms": u, "r_axis_m": m.distance, "r_b_m": m.r_b})
            rows.append(row)
    return rows


def wake_gas_rows(twin, altitudes):
    rows = []
    for h in altitudes:
        rel = _release(twin, h)
        for b in DEFAULT_BOUNDS:
            row = {"altitude_km": h, "bound": b, "assumed_core_speed_ms": ASSUMED_CORE_SPEED_MS,
                   "assumed_core_temperature_K": ASSUMED_CORE_T_K}
            for key, axis_key, s in (("core_speed_30km_ms", "core_speed_r_axis_m", CORE_SPEED_DISTANCE_M),
                                     ("core_T_300m_K", "core_T_r_axis_m", CORE_TEMPERATURE_DISTANCE_M)):
                # a point on the parent's own track s metres behind its position at release
                x = twin.track.position(rel.t_r - s / rel.e.U)
                gas, m, _ = _wake_gas(twin, rel, b, x)
                row["regime"] = m.regime
                if gas is not None:
                    row[key] = float(np.linalg.norm(gas.u)) if key.startswith("core_speed") else gas.T
                    row[axis_key] = m.distance
            rows.append(row)
    return rows


def _parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], prog="tables")
    p.add_argument("--diameters", nargs="+", type=float, default=list(DEFAULT_DIAMETERS))
    p.add_argument("--altitudes", nargs="+", type=float, default=list(DEFAULT_ALTITUDES))
    p.add_argument("--thetas", nargs="+", type=float, default=list(DEFAULT_THETAS))
    p.add_argument("--out", default=OUTPUT, help="output folder")
    return p


def _columns(rows, first):
    return list(first) + sorted({k for r in rows for k in r} - set(first))


def main(argv=None):
    args = _parser().parse_args(argv)
    out = os.fspath(args.out)
    os.makedirs(out, exist_ok=True)
    twin = get_twin()
    altitudes = kept_altitudes(args.altitudes, twin)

    rows = plume_size_rows(twin, altitudes)
    write_csv(os.path.join(out, "plume_size.csv"), rows, _columns(rows, (
        "altitude_km", "regime", "kn", "re_d", "p_cmhg_cm", "parent_speed_ms", "s_m", "t_since_passage_s", "r_b_low_m",
        "r_b_high_m", "laws_low", "laws_high", "beyond_s_park_limit")))
    rows = ballistic_rows(twin, altitudes, args.diameters)
    write_csv(os.path.join(out, "ballistic_coefficient.csv"), rows, _columns(rows, ("altitude_km", "diameter_mm")))
    rows = drama_rows(twin, parent_trajectory.load_parent())
    write_csv(os.path.join(out, "twin_vs_drama.csv"), rows, _columns(rows, ("time_s",)), comment=[
        "The twin runs about 6.5 s behind DRAMA (its speed at a given altitude is within a fraction of a percent, but it",
        "reaches each altitude later); the engine uses the twin's own clock and track, not DRAMA's.",
        "twin_downrange_km is the great-circle distance from the twin's start point at the DRAMA time (time_s minus the",
        "twin's start time); downrange_diff_km is that minus DRAMA's downrange since the top of the descending branch."])
    rows = affix_rows(twin, args.diameters, altitudes, args.thetas)
    write_csv(os.path.join(out, "affix_regimes.csv"), rows, _columns(rows, (
        "diameter_mm", "altitude_km", "theta_deg", "zone", "theta_inn", "theta_bif", "t_hat", "t_fz")))
    rows = drag_ratio_rows(twin, altitudes)
    write_csv(os.path.join(out, "drag_ratio_check.csv"), rows, _columns(rows, ("altitude_km", "bound", "drag_ratio",
                                                                                 "register_drag_ratio")))
    rows = wake_gas_rows(twin, altitudes)
    write_csv(os.path.join(out, "wake_gas_sanity.csv"), rows, _columns(rows, (
        "altitude_km", "bound", "regime", "core_speed_30km_ms", "assumed_core_speed_ms", "core_T_300m_K",
        "assumed_core_temperature_K")))
    print(f"tables written to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
