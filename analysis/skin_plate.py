"""Sub-plan 20 Task 9: level 2 of the layered-skin design's verification ladder (section 11.1, amended 2026-10-12), flow
off: the plate.

The q6p box of sub-plan 18 cannot be followed by corner columns: its 90-degree edges fold past the fixed surface
triangulation's limit (1.4 s on 2 mm elements, 0.9 s on 1 mm, found while planning), so the level-2 test is a plate
whose heated region tapers to nothing well inside its edges (decided 2026-10-12, Asha). The plate is 8 mm deep (the
box's length) and 110 x 110 mm across, of 2 mm elements at 700 K, Scheil material, no radiation. Its x = 0 face is
heated at 6 MW/m2 within 32 mm of its centre, the flux falling to zero at 45 mm along a half cosine, every column
receiving what the projection of its top dual cell on the face receives (`projected_flux`); liquid is removed as it
forms. The centre is one-dimensional while sideways conduction from the taper cannot reach it (the diffusion length
2 sqrt(alpha t) is about 19 mm at 2 s), so its recession -- read on the front corners within 3 mm of
the centre, as the top of each column with all its metal solid (`solid_top`) -- is compared at 1.0, 1.5 and 2.0 s with
the one-dimensional 1 um reference (`analysis/ablation_1d.py`) and the stage-1 fine box
(`data/reference_values/box_q6p_fine.json`). Criterion: within 5 % of both; mass closed to 1e-12 and energy to 1e-8.
The crater's wall is where the surface folds under heating-driven recession: the fold, the triangle and element
quality, rebuilds, faults, repeats and the wall time are recorded every step.

The hot plate (`--case hot`, added 2026-10-14, Asha): the same plate started at 880 K and heated at 2 MW/m2, where the
metal between the melt front and T_rigid -- the non-rigid layer the skins must hold -- is 1.5-1.7 mm thick in the exact
solution, so the skins' absorption is exercised as it is not at 700 K. Its reference is the one-dimensional solution at
1 um and 0.5 ms, recession and non-rigid depth (`data/reference_values/plate_hot_1d.json`, `analysis/ablation_1d.py
--t0 880 --q 2e6`). Criterion: recession within 5 % at 1.5 and 2.0 s (the onset at 1.0 s is reported, not judged: the
melting onset is open), the skins' flowing depth at the centre within 10 % of the exact non-rigid depth at 1.0, 1.5
and 2.0 s, and the centre's base no more than 1 K above T_rigid -- the skins hold the non-rigid metal, none is left in
the 3D model; mass closed to 1e-12 and energy to 1e-8.

    python analysis/skin_plate.py [--case cold|hot] [--dt 0.05] [--preset dev|production] [--h-mm 2] [--t-end 2]
                                  [--out PATH.json]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reentry_model import heating, material, mesh, skin, skinbody  # noqa: E402
from reentry_model.thermal import thermal_solver  # noqa: E402

Q0, T0 = 6.0e6, 700.0
LENGTH, WIDTH = 8e-3, 110e-3                 # m: the plate's depth along the heating and its extent across
R_FLAT, R_ZERO, R_READ = 32e-3, 45e-3, 3e-3  # m: full flux within R_FLAT, none beyond R_ZERO; recession read within R_READ
# R_FLAT - R_READ is half again the diffusion length 2 sqrt(alpha t), about 19 mm at 2 s: with full flux only to 12 mm
# the centre lost heat sideways and fell 9 % behind the exact solution by 2 s (measured while planning)
REFERENCE_1D = {1.0: 0.7890, 1.5: 2.7120, 2.0: 5.2190}       # mm, analysis/ablation_1d.py at 1 um and 0.5 ms (2026-10-12)
CASES = {"cold": (T0, Q0), "hot": (880.0, 2.0e6)}            # K, W/m2: the level-2 plate and the hot plate (2026-10-14)
HOT_FLOWING_TOL, HOT_BASE_TOL_K = 0.10, 1.0                    # the hot plate's flowing depth (relative) and base excess


def hot_reference():
    """The hot plate's one-dimensional reference at 1.0, 1.5 and 2.0 s: {t: (recession mm, non-rigid depth mm)}."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    doc = json.load(open(os.path.join(here, "data", "reference_values", "plate_hot_1d.json")))
    return {t: (float(np.interp(t, doc["t_s"], doc["recession_mm"])), float(np.interp(t, doc["t_s"], doc["nonrigid_mm"])))
            for t in (1.0, 1.5, 2.0)}


def fine_box():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    doc = json.load(open(os.path.join(here, "data", "reference_values", "box_q6p_fine.json")))
    return {t: float(np.interp(t, doc["t_s"], doc["r_total_mm"])) for t in (1.0, 1.5, 2.0)}


def taper(r, r_flat=R_FLAT, r_zero=R_ZERO):
    """1 within r_flat, a half cosine down to 0 at r_zero, 0 beyond."""
    x = np.clip((r - r_flat) / (r_zero - r_flat), 0.0, 1.0)
    return 0.5 * (1.0 + np.cos(np.pi * x))


def projected_flux(body, radius, r_flat=R_FLAT, r_zero=R_ZERO, q0=Q0):
    """Each column's flux [W/m2 over its top area]: q0 times the taper at its corner's starting radius `radius`, times
    the projection of its top dual cell on the face -- a third of each touching top triangle times that triangle's
    -n_x, floored at zero -- over its area. A column on a tilted crater wall receives the heat its projection catches."""
    s = body.surface
    a, b, c = (s.points[s.faces[:, i]] for i in range(3))
    n2 = np.cross(b - a, c - a)                                   # twice the area times the outward normal
    proj = np.zeros(len(s.areas))
    for i in range(3):
        np.add.at(proj, s.faces[:, i], np.clip(-n2[:, 0], 0.0, None) / 6.0)
    return q0 * taper(radius, r_flat, r_zero) * proj / s.areas


def solid_top(body):
    """Each column's top as if all its metal were solid: the base plus its mass over the solid density and its area,
    along its direction. The references hold their metal at one density, so this is the recession they measure; the
    skins' own top also swells as the mushy layer's density falls (by up to 0.15 mm before melting, measured)."""
    c = body.columns
    return c.points + (body.skins.column_mass() / (body.material.rho * c.area))[:, None] * c.normal


def run(dt=0.05, preset="dev", deform=4e-3, t_end=2.0, solver="skfem", verbose=True, h=2e-3, width=WIDTH, r_flat=R_FLAT,
        r_zero=R_ZERO, repeat_rule=skinbody.DEFAULT_REPEAT_RULE, T_start=T0, q0=Q0):
    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
    m = mesh.box_mesh(LENGTH, width, width, h)
    m.points[:, 1:] -= 0.5 * width                                 # the face's centre on the x axis
    body = skinbody.SkinBody(m, mat, thermal_solver(solver), mat.rho * m.volume(),
                             skin.SkinSettings(cell=skin.SKIN_PRESETS[preset]), T0=T_start, emissivity=0.0,
                             v_hat=(-1.0, 0.0, 0.0), deform_depth=deform, repeat_rule=repeat_rule)
    p = body._top0                                                 # the top surface at time zero: the plate's own
    front = np.abs(p[:, 0]) < 0.05e-3
    radius = np.where(front, np.hypot(p[:, 1], p[:, 2]), np.inf)
    centre = radius <= R_READ
    solid0, area0 = solid_top(body), body.surface.areas.copy()
    M0, rows, t, started, failure = body.mass(), [], 0.0, time.perf_counter(), None
    for k in range(int(round(t_end / dt))):
        try:
            q = projected_flux(body, radius, r_flat, r_zero, q0)
            body.advance(t + dt, dt, heating.HeatingResult(q, q0, 0.0, 0.0, 0.0))
        except RuntimeError as exc:                        # a move no rebuild can absorb: record where, and stop
            failure = {"t": t + dt, "error": str(exc), "fold_deg": body.moving.surface_fold()}
            print(json.dumps(failure), flush=True)
            break
        t += dt
        st = body.skin_stats()
        rows.append({"t": t, "recession_mm": float(((solid_top(body) - solid0) @ np.array([1.0, 0.0, 0.0]))[centre].mean() * 1e3),
                     "top_recession_mm": float(body._recession()[centre].mean() * 1e3),
                     "removed_recession_mm": float(body.column_removed[centre].sum() / (body.material.rho * area0[centre].sum()) * 1e3),
                     "centre_area_ratio": float(body.surface.areas[centre].sum() / area0[centre].sum()),
                     "T_top_centre": float(body.surface_temperature()[centre].mean()),
                     "flowing_mm": float(body.skins.depth_above(mat.T_rigid)[centre].mean() * 1e3),
                     "base_excess_K": float((body.solver.temperature()[body.moving.corners][centre] - mat.T_rigid).max()),
                     "mass_error": (body.mass() + body.removed_mass) / M0 - 1.0,
                     "energy_residual": body.energy_balance_residual(), "fold_deg": st["surface_fold_deg"],
                     "quality": st["surface_quality_min"], "zone_ratio": st["zone_ratio_min"], "rebuilds": st["rebuilds"],
                     "faults": st["jump_faults"], "repeats": st["interface_repeats"], "cells": st["skin_cells"],
                     "columns_over": st["repeat_columns_over"], "worst_ratio": st["mismatch_worst_ratio"],
                     "total_ratio": st["mismatch_total_ratio"],
                     "absorbed_kg": st["absorbed_mass_kg"], "handed_kg": st["handed_mass_kg"],
                     "wall_s": time.perf_counter() - started})
        if verbose and (k + 1) % max(1, int(round(0.25 / dt))) == 0:
            print(json.dumps({key: (round(v, 6) if isinstance(v, float) else v) for key, v in rows[-1].items()}), flush=True)
    return body, rows, int(centre.sum()), failure


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--case", default="cold", choices=tuple(CASES))
    a.add_argument("--dt", type=float, default=0.05)
    a.add_argument("--preset", default="dev", choices=tuple(skin.SKIN_PRESETS))
    a.add_argument("--deform-mm", type=float, default=4.0)
    a.add_argument("--t-end", type=float, default=2.0)
    a.add_argument("--solver", default="skfem")
    a.add_argument("--h-mm", type=float, default=2.0, help="3D element size [mm]")
    a.add_argument("--width-mm", type=float, default=WIDTH * 1e3)
    a.add_argument("--repeat-rule", default=skinbody.DEFAULT_REPEAT_RULE, choices=skinbody.REPEAT_RULES)
    a.add_argument("--r-flat-mm", type=float, default=R_FLAT * 1e3)
    a.add_argument("--r-zero-mm", type=float, default=R_ZERO * 1e3)
    a.add_argument("--out", default="skin_plate.json")
    args = a.parse_args()
    T_start, q0 = CASES[args.case]
    body, rows, n_centre, failure = run(args.dt, args.preset, args.deform_mm * 1e-3, args.t_end, args.solver,
                                        h=args.h_mm * 1e-3, width=args.width_mm * 1e-3, r_flat=args.r_flat_mm * 1e-3,
                                        r_zero=args.r_zero_mm * 1e-3, repeat_rule=args.repeat_rule, T_start=T_start, q0=q0)
    at = lambda t_, key: float(np.interp(t_, [x["t"] for x in rows], [x[key] for x in rows]))
    compare, verdict = {}, []
    fine = fine_box() if args.case == "cold" else None
    hot = hot_reference() if args.case == "hot" else None
    for t_ in (1.0, 1.5, 2.0):
        if not rows or t_ > rows[-1]["t"] + 1e-9:
            continue
        r, g = at(t_, "recession_mm"), at(t_, "removed_recession_mm")
        if fine is not None:
            compare[str(t_)] = {"centre_mm": r, "removed_mm": g, "fine_mm": fine[t_], "ref1d_mm": REFERENCE_1D[t_],
                                "rel_fine": r / fine[t_] - 1.0, "rel_1d": r / REFERENCE_1D[t_] - 1.0,
                                "removed_rel_1d": g / REFERENCE_1D[t_] - 1.0}
            verdict.append(abs(r / fine[t_] - 1.0) <= 0.05 and abs(r / REFERENCE_1D[t_] - 1.0) <= 0.05)
        else:
            ref_r, ref_nr = hot[t_]
            flow = at(t_, "flowing_mm")
            compare[str(t_)] = {"centre_mm": r, "ref1d_mm": ref_r, "rel_1d": r / ref_r - 1.0, "flowing_mm": flow,
                                "nonrigid_ref_mm": ref_nr, "flowing_rel": flow / ref_nr - 1.0,
                                "base_excess_K": at(t_, "base_excess_K")}
            verdict.append((t_ < 1.25 or abs(r / ref_r - 1.0) <= 0.05) and abs(flow / ref_nr - 1.0) <= HOT_FLOWING_TOL
                           and at(t_, "base_excess_K") <= HOT_BASE_TOL_K)
    summary = {"case": args.case, "T0_K": T_start, "q_Wm2": q0, "passes": bool(verdict) and all(verdict) and failure is None,
               "failure": failure, "completed_s": rows[-1]["t"] if rows else 0.0, "h_mm": args.h_mm, "dt": args.dt,
               "width_mm": args.width_mm, "r_flat_mm": args.r_flat_mm, "r_zero_mm": args.r_zero_mm,
               "repeat_rule": args.repeat_rule,
               "preset": args.preset, "deform_mm": args.deform_mm, "n_columns": body.columns.n, "n_centre": n_centre,
               "compare": compare, "worst_mass_error": max((abs(x["mass_error"]) for x in rows), default=0.0),
               "worst_energy_residual": max((abs(x["energy_residual"]) for x in rows), default=0.0),
               "fold_max_deg": max((x["fold_deg"] for x in rows), default=0.0),
               "rebuilds": rows[-1]["rebuilds"] if rows else 0, "faults": rows[-1]["faults"] if rows else 0,
               "repeats": rows[-1]["repeats"] if rows else 0, "wall_s": rows[-1]["wall_s"] if rows else 0.0,
               "wall_split": body.wall}
    print(json.dumps(summary, indent=1))
    json.dump({"summary": summary, "rows": rows}, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
