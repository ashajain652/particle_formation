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

Both plates have the skin runs' fine band (`--fine-mm 4 --fine-size-mm 0.7`, the defaults; Asha, 2026-10-14, option 1):
the 3D model's elements are finer to 4 mm below the surface, which keeps its 2 mm triangulation, and every rebuild makes
the band again; on 2 mm elements throughout the plates melted 2-6 % late.

The skins have their defaults: a 0.8 mm minimum with the finest cells in the top 0.4 mm (option B', Asha,
2026-10-16), so that on the cold plate, whose non-rigid layer is 0.39-0.49 mm thick, the skin-3D interface lies in rigid
metal below the front; with 0.4 mm for both the cold plate melted 3.6 % slow from 1.0 to 2.0 s, with B' 1.7 %.

Both cases run with the skins held at the solid's density (`--swell off`, the default here; Asha, 2026-10-14): the
references and the fine box hold their metal at one density, while the skins' own law (spec 5.2) swells the liquid by
15 %, which alone melts the exact plates sooner -- by 0.033 mm (cold) and 0.094 mm (hot) at 1.0 s -- and had hidden a
numerical delay of the same size. `--swell on` runs the skins' own law; the swelling itself is verified at level 1
(`tests/test_reentry_model_skin.py`, against the one-dimensional reference with the skins' density law).

With the flow on (`--flow on`, sub-plan 21 Task 8; spec 11.1's last sentence: "with flow on and a push along the
face: mass, momentum and energy close exactly and the result is physically sensible") the liquid stays in the skins
and flows; nothing is removed as it forms. Two settings drive it:

- `--push flight`: the face as a patch of the 100 mm sphere's shoulder PUSH_ALPHA (45 degrees) from the stagnation
  point at 93 s of its flight, where the surface flow's pressure gradient is largest. The flight direction is tilted by
  that angle about the z axis (`inclined`), so that along the face the meridian runs to -y, and every column of the
  face (`face_columns`: its corners at x = 0 whose normal is the face's, the edges' columns not) takes that point's
  forcing (`plate_forcing`): dp/ds -3.42e5 Pa/m and 80.3 Pa of gas shear down the meridian, and the deceleration of
  29.8 m/s2 along the flight direction -- 21.1 m/s2 out of the face, the front face's Rayleigh-Taylor pull, and 21.1
  m/s2 along it toward the nose, against the gas. On the plate's own flight direction, normal to the face, neither the
  pressure gradient nor the shear has a direction along it and the deceleration pulls straight out. The flight
  direction feeds nothing else here: the heating is the plate's own projected flux, and the stagnation and rear
  patches it picks are not read. `--push none` (the default) drives nothing but surface tension.
- `--release all`: Task 7's release from the top cells at RELEASE_ALL, far more than any column holds in the shortest
  flow sub-step, so that the liquid leaves as it forms -- at the liquidus (Task 7's ruling 3) rather than at T_feed,
  2 K higher, where instant removal takes it. With no push this is the flow-off plate's physics with the flow on, the
  regression against the references.

The flow sub-steps are chosen by the skins (Task 6d) or fixed (`--flow-cycles`). With the flow on the history gains the
skin body's flow columns (the momentum residual, the fastest cell, the roughness, the flow sub-steps' worst Courant
number, growth and linearisations, the released mass), and `--columns PATH.npz` keeps every column's mass, thickness,
top, top temperature, depths, top speed and volume flux at every macro step, for the sensibility checks;
`--save-state T` keeps the skins, the columns, the heating and the base at the first macro step ending at or after T
(`<out>.state.npz`), for the flow sub-step study and the cost split on the skins alone.

    python analysis/skin_plate.py [--case cold|hot] [--dt 0.05] [--preset dev|production] [--h-mm 2] [--t-end 2]
                                  [--swell off|on] [--fine-mm 4] [--fine-size-mm 0.7] [--out PATH.json]
                                  [--flow off|on] [--push none|flight] [--release none|all] [--flow-cycles N]
                                  [--columns PATH.npz] [--save-state T]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reentry_model import heating, material, mesh, movemesh, skin, skinbody  # noqa: E402
from reentry_model.thermal import thermal_solver  # noqa: E402

Q0, T0 = 6.0e6, 700.0
LENGTH, WIDTH = 8e-3, 110e-3                 # m: the plate's depth along the heating and its extent across
R_FLAT, R_ZERO, R_READ = 32e-3, 45e-3, 3e-3  # m: full flux within R_FLAT, none beyond R_ZERO; recession read within R_READ
# R_FLAT - R_READ is half again the diffusion length 2 sqrt(alpha t), about 19 mm at 2 s: with full flux only to 12 mm
# the centre lost heat sideways and fell 9 % behind the exact solution by 2 s (measured while planning)
REFERENCE_1D = {1.0: 0.7890, 1.5: 2.7120, 2.0: 5.2190}       # mm, analysis/ablation_1d.py at 1 um and 0.5 ms (2026-10-12)
# The references with the liquid leaving at the liquidus, as the release from the top cells takes it (sub-plan 21 Task 7's
# ruling 3), rather than at T_feed (stage 2's removal): analysis/ablation_1d.py --remove-at liquidus at 1 um and 0.5 ms
# (2026-10-24). The cold plate's recession [mm], and the hot plate's recession and non-rigid depth [mm]; the fine box
# is moved by the one-dimensional difference (`references`)
REFERENCE_1D_LIQUIDUS = {1.0: 0.8370, 1.5: 2.7720, 2.0: 5.2980}
HOT_LIQUIDUS = {1.0: (0.9280, 1.3530), 1.5: (1.9660, 1.4330), 2.0: (3.0750, 1.5370)}
CASES = {"cold": (T0, Q0), "hot": (880.0, 2.0e6)}            # K, W/m2: the level-2 plate and the hot plate (2026-10-14)
HOT_FLOWING_TOL, HOT_BASE_TOL_K = 0.10, 1.0                    # the hot plate's flowing depth (relative) and base excess
FINE = (4.0e-3, 0.7e-3)                                         # m: the fine band's depth and size, as a skin run's
# The flow (sub-plan 21 Task 8): the push along the face, the face as a patch of the 100 mm sphere's shoulder at
# PUSH_ALPHA from the stagnation point at 93 s of the flight -- the surface flow there (dp/ds [Pa/m], gas shear [Pa]) and
# the deceleration [m/s2], evaluated 2026-10-24 (measure/task8/forcing45.py): the largest pressure gradient of 6c's
# forcing, 3.4e5 Pa/m between 44 and 46.5 degrees, the shear 87 Pa at its largest at 49 degrees
PUSH_ALPHA = 45.0                                               # deg
PUSHES = {"none": None, "flight": (PUSH_ALPHA, -3.424e5, 80.29, 29.77)}
RELEASES = {"none": None, "all": 1.0e9}  # kg/(m2 s): a top area of 4 mm2 then releases 4 kg in a 1 ms flow sub-step,
                                         # some 1e5 times what a column holds: the liquid leaves as it forms
FACE_TILT = 5.0                          # deg: a column of the face has its corner at x = 0 and its normal within this
                                         # of the face's; the edges' columns, tilted 45 degrees, take no gas forcing


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


def references(case, at_liquidus=False):
    """The references a plate is judged against at 1.0, 1.5 and 2.0 s: (one-dimensional, fine box) recession [mm] for the
    cold plate, {t: (recession, non-rigid depth)} [mm] and None for the hot one; `at_liquidus`: with the liquid leaving
    at the liquidus (REFERENCE_1D_LIQUIDUS, HOT_LIQUIDUS), the fine box moved by the one-dimensional difference -- the
    liquid layer that stage 2's T_feed keeps on the front."""
    if case == "hot":
        return (HOT_LIQUIDUS if at_liquidus else hot_reference()), None
    fine = fine_box()
    if not at_liquidus:
        return REFERENCE_1D, fine
    return REFERENCE_1D_LIQUIDUS, {t: fine[t] + REFERENCE_1D_LIQUIDUS[t] - REFERENCE_1D[t] for t in fine}


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


def inclined(alpha_deg):
    """The flight direction tilted by alpha [deg] from the plate's own (-1, 0, 0) about the z axis: the heated face then
    lies alpha from the stagnation point, which is off the plate toward +y, and its meridian runs to -y."""
    a = np.radians(alpha_deg)
    return (-np.cos(a), np.sin(a), 0.0)


def face_columns(body, tilt=FACE_TILT):
    """The heated face's columns, (n,) bool: corners at x = 0 at time zero whose normal is within `tilt` [deg] of the
    face's (-1, 0, 0); the columns on the plate's edges, whose normals lean 45 degrees, are not the face's."""
    return (np.abs(body._top0[:, 0]) < 0.05e-3) & (-body._normal0[:, 0] >= np.cos(np.radians(tilt)))


def plate_forcing(body, face, push):
    """The gas side's push on the plate (skinbody.FlowForcing): with `push` (alpha [deg], dp/ds [Pa/m], gas shear [Pa],
    deceleration [m/s2]; PUSHES) the surface flow of that point of the shoulder on every column of the face and the
    deceleration along the body's flight direction, which `inclined(alpha)` gave it; None: nothing (surface tension
    alone drives the flow)."""
    n = body.columns.n
    if push is None:
        return skinbody.FlowForcing(np.zeros(n), np.zeros(n), np.zeros(3))
    _, dp_ds, tau, decel = push
    return skinbody.FlowForcing(np.where(face, dp_ds, 0.0), np.where(face, tau, 0.0), decel * body.v_hat)


def column_record(body):
    """Every column's state for the sensibility checks: its mass [kg], thickness [m], top [m, (n, 3)], top temperature
    [K], flowing and liquid depths [m], top speed [m/s] and volume flux [m2/s, (n, 3)] -- the sum of its cells'
    thicknesses times their velocities, what crosses a unit width of the surface."""
    sk, mat = body.skins, body.material
    dz = sk.dz()
    return {"mass": sk.column_mass(), "thickness": dz.sum(axis=1), "top": body.surface.points.copy(),
            "T_top": sk.T[:, 0].copy(), "flowing": sk.depth_above(mat.T_rigid), "liquid": sk.depth_above(mat.T_liquidus),
            "top_speed": np.linalg.norm(sk.U[:, 0], axis=1), "flux": (dz[..., None] * sk.U).sum(axis=1)}


def save_state(path, body, q):
    """The skins, their columns (base corners, triangles, top areas), the heating q [W/m2], the base temperatures and the
    3D model's gradient below them: what `skin.advance` needs to run a part on the skins alone."""
    np.savez(path, **dict(zip(("area", "m", "H", "T", "count", "U", "flow_memory"), body.skins.state())),
             points=body.columns.points, faces=body.columns.faces, A_top=body.surface.areas, q=q,
             T_b=body.solver.temperature()[body.moving.corners], g_base=body._normal_gradient(), v_hat=body.v_hat,
             rho_l=body.skins.rho_l, face=face_columns(body))


def run(dt=0.05, preset="dev", deform=4e-3, t_end=2.0, solver="skfem", verbose=True, h=2e-3, width=WIDTH, r_flat=R_FLAT,
        r_zero=R_ZERO, repeat_rule=skinbody.DEFAULT_REPEAT_RULE, T_start=T0, q0=Q0, swell=False, fine=FINE, flow=False,
        push=None, release=None, flow_cycles=None, record=False, save_at=None, save_path=None):
    """The plate (module docstring) to t_end at macro steps dt; returns (body, rows, centre columns, failure) and, with
    `record`, every macro step's `column_record` as a fifth item. `flow` turns the flow on (the liquid kept), driven by
    `push` (PUSHES) with the release `release` [kg/(m2 s), RELEASES] and `flow_cycles` (None: chosen); `save_at` [s]
    keeps the skins' state at the first macro step ending at or after it in `save_path`."""
    if not flow and (push is not None or release is not None or flow_cycles is not None):
        raise ValueError("a push, a release and the flow sub-steps need the flow on")
    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
    m = mesh.box_mesh(LENGTH, width, width, h)
    m.points[:, 1:] -= 0.5 * width                                 # the face's centre on the x axis
    if fine is not None and fine[0] > 0.0:                         # the fine band (depth, size) m, as a skin run's
        m = movemesh.refined_mesh(m, fine[0], min(fine[1], h))
    settings = skin.SkinSettings(cell=skin.SKIN_PRESETS[preset], swell=swell)
    if flow:                                                       # the liquid stays and flows (sub-plan 21 Task 8)
        body = skinbody.SkinBody(m, mat, thermal_solver(solver), mat.rho * m.volume(), settings, T0=T_start,
                                 emissivity=0.0, v_hat=(-1.0, 0.0, 0.0) if push is None else inclined(push[0]),
                                 deform_depth=deform, repeat_rule=repeat_rule, remove=False, flow=True,
                                 flow_cycles=flow_cycles, release=None if release is None else skin.prescribed_release(release))
    else:
        body = skinbody.SkinBody(m, mat, thermal_solver(solver), mat.rho * m.volume(), settings, T0=T_start,
                                 emissivity=0.0, v_hat=(-1.0, 0.0, 0.0), deform_depth=deform, repeat_rule=repeat_rule)
    p = body._top0                                                 # the top surface at time zero: the plate's own
    front = np.abs(p[:, 0]) < 0.05e-3
    radius = np.where(front, np.hypot(p[:, 1], p[:, 2]), np.inf)
    centre = radius <= R_READ
    face = face_columns(body)
    solid0, area0 = solid_top(body), body.surface.areas.copy()
    M0, rows, t, started, failure = body.mass(), [], 0.0, time.perf_counter(), None
    records = [column_record(body)] if record else None
    for k in range(int(round(t_end / dt))):
        try:
            q = projected_flux(body, radius, r_flat, r_zero, q0)
            if flow:
                body.advance(t + dt, dt, heating.HeatingResult(q, q0, 0.0, 0.0, 0.0),
                             forcing=plate_forcing(body, face, push))
            else:
                body.advance(t + dt, dt, heating.HeatingResult(q, q0, 0.0, 0.0, 0.0))
        except RuntimeError as exc:                        # a move no rebuild can absorb: record where, and stop
            failure = {"t": t + dt, "error": str(exc), "fold_deg": body.moving.surface_fold()}
            print(json.dumps(failure), flush=True)
            break
        t += dt
        if save_at is not None and t >= save_at - 1e-9 and save_path is not None:
            save_state(save_path, body, q)
            save_at = None
        if record:
            records.append(column_record(body))
        st = body.skin_stats()
        rows.append({"t": t, "recession_mm": float(((solid_top(body) - solid0) @ np.array([1.0, 0.0, 0.0]))[centre].mean() * 1e3),
                     "top_recession_mm": float(body._recession()[centre].mean() * 1e3),
                     "removed_recession_mm": float(body.column_removed[centre].sum() / (body.material.rho * area0[centre].sum()) * 1e3),
                     "centre_area_ratio": float(body.surface.areas[centre].sum() / area0[centre].sum()),
                     "T_top_centre": float(body.surface_temperature()[centre].mean()),
                     "flowing_mm": float(body.skins.depth_above(mat.T_rigid)[centre].mean() * 1e3),
                     "base_excess_K": float((body.solver.temperature()[body.moving.corners][centre] - mat.T_rigid).max()),
                     "mass_error": (body.mass() + body.removed_mass + body.released_mass) / M0 - 1.0,
                     "energy_residual": body.energy_balance_residual(), "fold_deg": st["surface_fold_deg"],
                     "quality": st["surface_quality_min"], "zone_ratio": st["zone_ratio_min"], "rebuilds": st["rebuilds"],
                     "faults": st["jump_faults"], "repeats": st["interface_repeats"], "cells": st["skin_cells"],
                     "columns_over": st["repeat_columns_over"], "worst_ratio": st["mismatch_worst_ratio"],
                     "total_ratio": st["mismatch_total_ratio"],
                     "absorbed_kg": st["absorbed_mass_kg"], "handed_kg": st["handed_mass_kg"],
                     "wall_s": time.perf_counter() - started})
        if flow:                                           # the skin body's flow columns (sub-plan 21 Tasks 6c-7)
            sk = body.skins
            rows[-1].update({key: st[key] for key in ("momentum_residual", "speed_max_ms", "roughness_max",
                                                      "flow_courant_max", "flow_growth_max", "flow_passes_max",
                                                      "flow_cycles_max", "released_mass_kg", "wall_flow_s",
                                                      "liquid_depth_max_mm", "flowing_depth_max_mm")})
            rows[-1].update({"liquid_mm": float(sk.depth_above(mat.T_liquidus)[centre].mean() * 1e3),
                             "kinetic_J": 0.5 * float((sk.m * (sk.U ** 2).sum(axis=2)).sum()),
                             "wall_3d_s": st["wall_3d_s"], "wall_skins_s": st["wall_skins_s"],
                             "wall_move_s": st["wall_move_s"], "flow_seconds": dict(skin.FLOW_SECONDS)})
        if verbose and (k + 1) % max(1, int(round(0.25 / dt))) == 0:
            print(json.dumps({key: (round(v, 6) if isinstance(v, float) else v) for key, v in rows[-1].items()}), flush=True)
    if record:
        return body, rows, int(centre.sum()), failure, records
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
    a.add_argument("--fine-mm", type=float, default=FINE[0] * 1e3, help="the fine band's depth [mm]; 0: none")
    a.add_argument("--fine-size-mm", type=float, default=FINE[1] * 1e3, help="the fine band's element size asked [mm]")
    a.add_argument("--swell", choices=("off", "on"), default="off",
                   help="off: the skins at the solid's density, as the references hold their metal (default); on: the "
                        "skins' own law, each cell at its phase's density")
    a.add_argument("--out", default="skin_plate.json")
    a.add_argument("--flow", choices=("off", "on"), default="off",
                   help="on: the skins' liquid stays and flows (sub-plan 21 Task 8); off: removed as it forms (default)")
    a.add_argument("--push", choices=tuple(PUSHES), default="none",
                   help="with --flow on, flight: the face as the 100 mm sphere's shoulder at 45 degrees at 93 s")
    a.add_argument("--release", choices=tuple(RELEASES), default="none",
                   help="with --flow on, all: the liquid released from the top cells as it forms")
    a.add_argument("--flow-cycles", type=int, default=None, help="with --flow on: flow sub-steps fixed (default chosen)")
    a.add_argument("--columns", default=None, help="keep every column's state at every macro step in this npz")
    a.add_argument("--save-state", type=float, default=None, help="keep the skins' state at this time [s]")
    args = a.parse_args()
    if args.flow == "off" and (args.push != "none" or args.release != "none" or args.flow_cycles is not None):
        a.error("--push, --release and --flow-cycles need --flow on")
    T_start, q0 = CASES[args.case]
    out = run(args.dt, args.preset, args.deform_mm * 1e-3, args.t_end, args.solver,
              h=args.h_mm * 1e-3, width=args.width_mm * 1e-3, r_flat=args.r_flat_mm * 1e-3,
              r_zero=args.r_zero_mm * 1e-3, repeat_rule=args.repeat_rule, T_start=T_start, q0=q0,
              swell=args.swell == "on", fine=(args.fine_mm * 1e-3, args.fine_size_mm * 1e-3), flow=args.flow == "on",
              push=PUSHES[args.push], release=RELEASES[args.release], flow_cycles=args.flow_cycles,
              record=args.columns is not None, save_at=args.save_state,
              save_path=os.path.splitext(args.out)[0] + ".state.npz")
    body, rows, n_centre, failure = out[:4]
    if args.columns is not None:
        np.savez(args.columns, t=np.r_[0.0, [x["t"] for x in rows]], face=face_columns(body), top0=body._top0,
                 area=body.columns.area, **{key: np.array([rec[key] for rec in out[4]]) for key in out[4][0]})
    at = lambda t_, key: float(np.interp(t_, [x["t"] for x in rows], [x[key] for x in rows]))
    compare, verdict = {}, []
    ref, fine = references(args.case, at_liquidus=args.release == "all")     # the release takes liquid at the liquidus
    hot = ref if args.case == "hot" else None
    for t_ in (1.0, 1.5, 2.0):
        if not rows or t_ > rows[-1]["t"] + 1e-9:
            continue
        r, g = at(t_, "recession_mm"), at(t_, "removed_recession_mm")
        if fine is not None:
            compare[str(t_)] = {"centre_mm": r, "removed_mm": g, "fine_mm": fine[t_], "ref1d_mm": ref[t_],
                                "rel_fine": r / fine[t_] - 1.0, "rel_1d": r / ref[t_] - 1.0,
                                "removed_rel_1d": g / ref[t_] - 1.0}
            verdict.append(abs(r / fine[t_] - 1.0) <= 0.05 and abs(r / ref[t_] - 1.0) <= 0.05)
        else:
            ref_r, ref_nr = hot[t_]
            flow = at(t_, "flowing_mm")
            compare[str(t_)] = {"centre_mm": r, "ref1d_mm": ref_r, "rel_1d": r / ref_r - 1.0, "flowing_mm": flow,
                                "nonrigid_ref_mm": ref_nr, "flowing_rel": flow / ref_nr - 1.0,
                                "base_excess_K": at(t_, "base_excess_K")}
            verdict.append((t_ < 1.25 or abs(r / ref_r - 1.0) <= 0.05) and abs(flow / ref_nr - 1.0) <= HOT_FLOWING_TOL
                           and at(t_, "base_excess_K") <= HOT_BASE_TOL_K)
    judged = args.push == "none"                     # a push along the face has no exact reference (spec 11.1)
    summary = {"case": args.case, "T0_K": T_start, "q_Wm2": q0,
               "passes": judged and bool(verdict) and all(verdict) and failure is None,
               "failure": failure, "completed_s": rows[-1]["t"] if rows else 0.0, "h_mm": args.h_mm, "dt": args.dt,
               "width_mm": args.width_mm, "r_flat_mm": args.r_flat_mm, "r_zero_mm": args.r_zero_mm,
               "repeat_rule": args.repeat_rule, "swell": args.swell, "fine_mm": args.fine_mm,
               "fine_size_mm": args.fine_size_mm, "n_elements": body.mesh.n_elements,
               "preset": args.preset, "deform_mm": args.deform_mm, "n_columns": body.columns.n, "n_centre": n_centre,
               "compare": compare, "worst_mass_error": max((abs(x["mass_error"]) for x in rows), default=0.0),
               "worst_energy_residual": max((abs(x["energy_residual"]) for x in rows), default=0.0),
               "fold_max_deg": max((x["fold_deg"] for x in rows), default=0.0),
               "rebuilds": rows[-1]["rebuilds"] if rows else 0, "faults": rows[-1]["faults"] if rows else 0,
               "repeats": rows[-1]["repeats"] if rows else 0, "wall_s": rows[-1]["wall_s"] if rows else 0.0,
               "wall_split": body.wall}
    if args.flow == "on":
        summary.update({"flow": "on", "push": args.push, "release": args.release, "flow_cycles": args.flow_cycles,
                        "judged": judged, "released_kg": body.released_mass,
                        "worst_momentum_residual": max(abs(x["momentum_residual"]) for x in rows) if rows else 0.0,
                        "roughness_max": max(x["roughness_max"] for x in rows) if rows else 0.0,
                        "speed_max_ms": max(x["speed_max_ms"] for x in rows) if rows else 0.0,
                        "flow_seconds": dict(skin.FLOW_SECONDS)})
    print(json.dumps(summary, indent=1))
    json.dump({"summary": summary, "rows": rows}, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
