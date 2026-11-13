#!/usr/bin/env python3
"""The probe of a saved skin flight under each near-liquidus viscosity law (level 4 of the near-liquidus viscosity
design, docs/superpowers/specs/2026-11-12-near-liquidus-viscosity-design.md, section 8): the skins' flow alone, run from
the saved state with the gas loads held at the saved step's, and the spraying modes evaluated along the way without
releasing any metal.

    "$PY" analysis/near_liquidus_probe.py --state <pkl> --case {bridge,composite,composite_fi1,composite_fi025}
        [--flow-time 0.1] [--out output/near_liquidus_viscosity/probe] [--single-column]

The cases: `bridge`, the present law; `composite`, the new law with its internal solid fraction 0.5; `composite_fi1` and
`composite_fi025`, the same law with the internal solid fraction set to 1 or 0.25 in this process only
(`rheology.INTERNAL_SOLID_FRACTION`, read at call time). Nothing in the package is changed.

What a run does:

- loads the state (`cli.load_state`) and sets the skins' law (`SkinSettings.rheology`, the body's and the skins' one
  object);
- recomputes the gas loads with the current code, `SkinBody.flow_forcing` at the saved step's aero state
  (`Simulator.aero_state` at the saved time and trajectory state, as `CoupledRun.loads_at` calls it), and records their
  difference from the forcing the body saved for that step's end (`_forcing[1]`);
- builds the flow (`skin.Flow`) from those loads as the body's parts do, with no release and the adaptive count of flow
  sub-steps, and calls `skin.flow_substeps` once per heat sub-step (the skins' 10 ms), each followed by the heat sub-step's
  regrid, until `--flow-time`. No heat, no 3D conduction, no trajectory step and no release run;
- at the start, every 20 ms and at the end, on the windward columns (the angle of the top surface's normal to the flight
  direction, `SkinBody.load_angles`, below 60 degrees), records: the liquid's and the exposed slurry's spraying modes
  (`skinbody.girin_plan(..., slurry=True)`, as the body calls it at a macro step's start) under two sets of release rules,
  the saved flight's (gates "before", slip band off) and the recent d050 work's (gates "2026-11-10", slip band on), set on
  the body's spray model and restored after -- each mode's Weber number against the critical value, the gate that decided
  it and its release rate; the closure check, the slurry mode recomputed through `skinbody.slurry_modes` with its flowing
  depth ended at the 10 % solid depth (`Skins.depth_above` at the temperature where the material's Scheil curve stands at
  10 % solid, found by inverting `liquid_fraction`), beside the momentum step's top-cell speeds; and the top three
  cells' solid fraction, viscosity over the liquid's and the stress each carries (the mean over its two faces of the face
  stresses, as `skin._momentum_columns`' `rates()` computes them, at the cells' present velocities); and where the skins'
  metal is, in bands of that angle: all of it, the flowing cells' and the cells under 10 % solid. Each heat sub-step's
  flow sub-steps and wall time are recorded too;
- with `--single-column`, instead of the flow: the windward column with the median flowing depth among those holding
  slurry cells, run by `skin.momentum_step` alone on `Skins.rows([i])` with that column's push and shear, 1 ms sub-steps,
  until the largest velocity change is below 1e-9 of the largest velocity or 2000 sub-steps have run; its profiles of
  velocity, viscosity and stress are written.

The spraying modes' inputs and results are read by recording the arguments and results of every `SprayModel.evaluate` call
that `girin_plan` and `slurry_modes` make (the class's method is wrapped for the duration of the call and restored), so
the gate that decided each mode is read from the model's own numbers; the closure variant's arguments are those
`girin_plan` passes, and at every sample the slurry mode with the real flowing depth is checked to reproduce the plan's
rate exactly before the variant is recomputed.

Outputs (git-ignored): `<out>/<case>.npz` (per sample, per windward column), `<out>/<case>.json` (the medians and 10th to
90th percentiles over the windward columns, the gates' counts, the costs and the loads' difference); with
`--single-column`, `<out>/<case>_single.npz` and `<case>_single.json`.

The composite law's suspension part below 10 % solid and its internal solid fraction are assumed (design 3.7: no aluminium
measurement exists between the liquid and 10 % solid, and the internal fraction is unmeasured for 7075); the jump at 10 %
solid is unresolved. See the module docstring of reentry_model/shared/rheology.py."""
import argparse
import collections
import contextlib
import json
import math
import os
import subprocess
import sys
import time

import numpy as np
from scipy.optimize import brentq

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)
from reentry_model import OUTPUT_DIR  # noqa: E402
from reentry_model import cli  # noqa: E402
from reentry_model.shared import rheology, spray as spray_mod  # noqa: E402
from reentry_model.skyscrapers import skin as skin_mod, skinbody  # noqa: E402

CASES = {"bridge": ("bridge", None), "composite": ("composite", None), "composite_fi1": ("composite", 1.0),
         "composite_fi025": ("composite", 0.25)}           # (law, internal solid fraction set in this process or None)
RULES = {"saved": ("before", False), "recent": ("2026-11-10", True)}     # (spray_gates, slip_band)
WINDWARD_DEG = 60.0               # columns whose top normal is within this of the flight direction are sampled
SAMPLE_EVERY = 0.02               # s of flow between samples
SOLID_LIMIT = 0.1                 # the closure variant ends the slurry layer where the Scheil curve stands at this solid
SINGLE_DT, SINGLE_MAX, SINGLE_TOL = 1.0e-3, 2000, 1.0e-9   # the single column: sub-step [s], most sub-steps, steady test
TOP = 3                           # top cells recorded
BANDS_DEG = (0.0, 30.0, 60.0, 90.0, 180.01)   # bands of the angle to the flight direction for the metal's whereabouts
BRANCH_NAMES = {spray_mod.BRANCH_THICK: "thick", spray_mod.BRANCH_THIN: "thin", spray_mod.BRANCH_RAREFIED: "rarefied",
                spray_mod.BRANCH_RT: "Rayleigh-Taylor", spray_mod.BRANCH_REGIME2: "regime 2"}


# -- reading the spraying modes ----------------------------------------------------------------------------------------
@contextlib.contextmanager
def recorded():
    """Record every `SprayModel.evaluate` call (its model, arguments and result) made inside the block."""
    calls, original = [], spray_mod.SprayModel.evaluate

    def evaluate(self, flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius=0.05, **kw):
        res = original(self, flow, state, b, delta_m, v_s, windward, dt, areas, m_f, radius, **kw)
        calls.append({"model": self, "flow": flow, "state": state, "b": np.array(b, dtype=float),
                      "delta_m": np.array(delta_m, dtype=float), "v_s": np.array(v_s, dtype=float),
                      "windward": np.array(windward, dtype=bool), "kw": kw, "res": res})
        return res

    spray_mod.SprayModel.evaluate = evaluate
    try:
        yield calls
    finally:
        spray_mod.SprayModel.evaluate = original


def gates(call):
    """Per patch of one recorded `evaluate` call, the gate that decided its mode, in the order the model applies them: a
    label "sprays: <branch>" where it releases, else the first test that stopped it."""
    m, res, kw, flow, state = call["model"], call["res"], call["kw"], call["flow"], call["state"]
    b, dm, v_s, w = call["b"], call["delta_m"], call["v_s"], call["windward"]
    extent = kw.get("extent")
    on_slurry, cap = kw.get("on_slurry"), kw.get("cap_depth")
    out = []
    for i in range(b.size):
        br = int(res.branch[i])
        name = BRANCH_NAMES.get(br, "none")
        no_cap = cap is not None and not np.asarray(cap, dtype=float)[i] > 0.0
        L = None if extent is None else np.asarray(extent, dtype=float)[i:i + 1]
        if not w[i]:
            label = "leeward"
        elif not b[i] >= spray_mod.B_MIN:
            label = "no layer to spray"
        elif res.mdot[i] > 0.0:
            label = "sprays: " + name
        elif br == -1:
            if on_slurry is not None and bool(np.asarray(on_slurry)[i]) and not np.isfinite(dm[i]):
                label = "held: liquid on slurry with no conjugate depth"
            elif m.spray_gates != "before" and not np.isfinite(dm[i]):
                label = "continuum gate: no conjugate depth"
            else:
                label = "none: unexplained"
        elif br == spray_mod.BRANCH_THICK:
            if not res.we_s[i] > m.we_critical:
                label = "thick: We_s at or below We_cr"
            else:
                delta_f, im_f, _ = m.table(np.array([res.we_s[i]]))
                lam = 2.0 * np.pi * dm[i] / delta_f[0] if delta_f[0] > 0.0 else np.inf
                if not (im_f[0] > 0.0 and v_s[i] > 0.0):
                    label = "thick: no growth"
                elif not spray_mod.wave_fits(np.array([lam]), L)[0]:
                    label = "thick: wave does not fit"
                elif no_cap:
                    label = "no depth to release"
                else:
                    label = "thick: unexplained"
        elif br in (spray_mod.BRANCH_THIN, spray_mod.BRANCH_RAREFIED, spray_mod.BRANCH_REGIME2):
            if br == spray_mod.BRANCH_RAREFIED:
                mach, flux = state.ma, state.freestream.rho * state.V ** 2
            else:
                mach, flux = flow.mach_e[i], flow.rho_e[i] * flow.u_eff[i] ** 2
            lam, tau = spray_mod.thin_film_mode(np.array([mach], dtype=float), np.array([flux], dtype=float), m.liquid)
            if not (np.isfinite(lam[0]) and lam[0] > 0.0 and tau[0] > 0.0):
                label = name + ": no wavelength"
            elif res.gamma_tau[i] >= 1.0:
                label = name + ": viscous damping (gamma tau* >= 1)"
            elif res.outside is not None and res.outside[i]:
                label = name + ": film inside the floor's Stokes layer"
            elif not spray_mod.wave_fits(lam, L)[0]:
                label = name + ": wave does not fit"
            elif no_cap:
                label = "no depth to release"
            else:
                label = name + ": unexplained"
        else:
            label = "no depth to release" if no_cap else name + ": unexplained"
        out.append(label)
    return out


def flow_counts(sk, mat):
    """(flowing cells, liquid cells, flowing run per column, liquid run per column), as `girin_plan` counts them."""
    flowing = sk.valid & np.isfinite(skin_mod.cell_viscosity(sk, 0.0))
    liquid = flowing & (sk.T >= mat.T_liquidus)
    return flowing, liquid, np.cumprod(flowing, axis=1).sum(axis=1), np.cumprod(liquid, axis=1).sum(axis=1)


def slurry_rows(calls, rows, n):
    """Per column (n,), from `slurry_modes`' recorded one-row calls (one per row of `rows`, in order): rate, We_s, branch,
    surface speed, gate (nan, -1 and "" off the rows)."""
    if len(calls) != len(rows):
        raise RuntimeError("expected {} recorded slurry calls, got {}".format(len(rows), len(calls)))
    rate, we, v = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    branch, gate = np.full(n, -1), np.full(n, "", dtype=object)
    for i, c in zip(rows, calls):
        rate[i], we[i], v[i], branch[i] = c["res"].mdot[0], c["res"].we_s[0], c["v_s"][0], c["res"].branch[0]
        gate[i] = gates(c)[0]
    return rate, we, branch, v, gate


def modes(body, forcing, state, dt_macro, T10, idx):
    """The spraying modes on the columns idx under both rule sets, without releasing anything: a dict of arrays."""
    sk, mat, model = body.skins, body.material, body.spray
    theta = body.load_angles()
    surface, edges, n = body.surface, body.columns.edges, body.skins.n
    radius, nose = body.transverse_radius(), body.nose_radius()
    flow = forcing.result
    windward = theta <= 0.5 * np.pi
    w_n = flow.deceleration * np.maximum(np.cos(theta), 0.0)
    _, _, n_flow, n_liq = flow_counts(sk, mat)
    rows = np.flatnonzero((n_flow > n_liq) & windward)
    F, F10 = sk.depth_above(mat.T_rigid), sk.depth_above(T10)
    own = {key: key in vars(model) for key in ("spray_gates", "slip_band")}
    kept = (model.spray_gates, model.slip_band)
    out = {}
    try:
        for rule, (spray_gates, band) in RULES.items():
            model.spray_gates, model.slip_band = spray_gates, band
            with recorded() as calls:
                plan = skinbody.girin_plan(sk, flow, state, model, theta, surface.areas, edges, radius, dt_macro, True,
                                           points=surface.points.copy(),
                                           centre=spray_mod.enclosed_centroid(surface.points, surface.faces),
                                           v_hat=body.v_hat, nose_radius=nose)
            liq = calls[0]
            if liq["b"].size != n:
                raise RuntimeError("the first recorded call is not the liquid's mode over every column")
            s_rate, s_we, s_branch, s_v, s_gate = slurry_rows(calls[1:], rows, n)
            if not np.array_equal(np.nan_to_num(s_rate), plan.rate_slurry):
                raise RuntimeError("the recorded slurry rates differ from the plan's")
            args = (sk, flow, state, model, windward, w_n, surface.areas, edges, radius, dt_macro)
            with recorded() as again:
                check, _ = skinbody.slurry_modes(*args, F, n_liq, n_flow, theta=theta, nose_radius=nose)
            if not np.array_equal(check, plan.rate_slurry):
                raise RuntimeError("slurry_modes with girin_plan's arguments does not reproduce the plan's slurry rate")
            with recorded() as closure:
                rate10, out10 = skinbody.slurry_modes(*args, F10, n_liq, n_flow, theta=theta, nose_radius=nose)
            c_rate, c_we, c_branch, c_v, c_gate = slurry_rows(closure, rows, n)
            lg = np.array(gates(liq), dtype=object)
            mu_l = mat.liquid.mu
            got = {"liquid_rate": plan.rate, "liquid_we": liq["res"].we_s, "liquid_branch": liq["res"].branch,
                   "liquid_v": liq["v_s"], "liquid_delta_m": liq["delta_m"], "liquid_gate": lg,
                   "slurry_rate": s_rate, "slurry_we": s_we, "slurry_branch": s_branch, "slurry_v": s_v, "slurry_gate": s_gate,
                   "slurry_mu_ratio": plan.slurry["mu"] / mu_l, "slurry_gamma": plan.slurry["gamma"],
                   "slurry_delta_m": plan.slurry["delta_m"], "slurry_depth": plan.slurry["depth"],
                   "slurry_r": plan.slurry["r"],
                   "closure_rate": c_rate, "closure_we": c_we, "closure_branch": c_branch, "closure_v": c_v,
                   "closure_gate": c_gate, "closure_depth": out10["depth"], "closure_delta_m": out10["delta_m"]}
            for key, value in got.items():
                out["{}_{}".format(rule, key)] = np.asarray(value)[idx]
    finally:
        model.spray_gates, model.slip_band = kept
        for key, mine in own.items():
            if not mine:
                vars(model).pop(key, None)                 # back to the class default, as the saved model had it
    out["F"], out["F10"], out["n_flow"], out["n_liq"] = F[idx], F10[idx], n_flow[idx], n_liq[idx]
    return out


# -- the cells' stresses -----------------------------------------------------------------------------------------------
def stresses(sk, shear):
    """Every cell's viscosity [Pa s], face stresses [Pa, (n, K+1), the top face first], the stress it carries (the mean
    of its two faces' [Pa]) and its shear rate [1/s], at the cells' present velocities, as `_momentum_columns` builds them:
    a floating crust at CRUST_VISCOSITY, the drags in series, the gas shear on a moving top, the anchored metal's drag on
    a flowing bottom cell; the viscosities at their own rates by the same fixed point, the velocities held."""
    n, K = sk.m.shape
    rows = np.arange(n)
    last = np.maximum(sk.count - 1, 0)
    valid, dz = sk.valid, sk.dz()
    flowing = valid & np.isfinite(skin_mod.cell_viscosity(sk, 0.0))
    anchored = skin_mod.anchored_cells(sk, flowing)
    moving = valid & ~anchored
    crust = moving & ~flowing
    top_moves = moving[:, 0]
    bottom_flows = flowing[rows, last] & (sk.count > 0)
    U = np.where(moving[..., None], sk.U, 0.0)
    shear = np.broadcast_to(np.asarray(shear, dtype=float), (n, 3))

    def viscosities(gamma):
        return np.where(crust, skin_mod.CRUST_VISCOSITY, skin_mod.cell_viscosity(sk, gamma))

    def faces(eta):
        D = skin_mod.series_drag(sk, eta)
        Db = np.where(bottom_flows, eta[rows, last] / (0.5 * np.maximum(dz[rows, last], skin_mod.DZ_FLOOR)), 0.0)
        Db = np.where(np.isfinite(Db), Db, 0.0)
        s = np.zeros((n, K + 1))
        s[:, 1:K] = D * np.linalg.norm(U[:, :-1] - U[:, 1:], axis=2)
        s[:, 0] = np.where(top_moves, np.linalg.norm(shear, axis=1), 0.0)
        s[rows, last + 1] = Db * np.linalg.norm(U[rows, last], axis=1)
        return s

    def rates(eta):
        s = faces(eta)
        face = 0.5 * (s[:, :-1] + s[:, 1:])
        return np.where(flowing, np.divide(face, eta, out=np.zeros_like(face), where=np.isfinite(eta)), 0.0)

    eta = viscosities(0.0)
    eta = viscosities(rates(eta))
    for _ in range(skin_mod.SHEAR_PASSES):
        new = viscosities(rates(eta))
        change = np.abs(new[flowing] - eta[flowing]) / eta[flowing]
        eta = new
        if change.size == 0 or change.max() <= skin_mod.SHEAR_TOL:
            break
    s = faces(eta)
    return eta, s, 0.5 * (s[:, :-1] + s[:, 1:]), rates(eta)


def cells(sk, mat, flow, idx):
    """The top TOP cells of the columns idx: temperature, solid fraction, viscosity over the liquid's, stress carried,
    shear rate, speed; and the top cell's speed. nan beyond a column's cells."""
    sub = sk.rows(idx)
    eta, _, carried, gamma = stresses(sub, flow.shear[idx])
    k = min(TOP, sub.m.shape[1])
    pad = lambda x: np.concatenate([x[:, :k], np.full((len(idx), TOP - k), np.nan)], axis=1)
    on = pad(sub.valid.astype(float)) > 0.0
    f_s = 1.0 - sub.props.liquid_fraction(sub.T)
    out = {"top_T": sub.T, "top_fs": f_s, "top_mu_ratio": eta / mat.liquid.mu, "top_stress": carried, "top_gamma": gamma,
           "top_speed": np.linalg.norm(sub.U, axis=2)}
    return {key: np.where(on, pad(value), np.nan) for key, value in out.items()}


def masses(sk, mat, theta, T10):
    """Where the skins' metal is [kg]: in bands of the angle to the flight direction (BANDS_DEG), all of it, the flowing
    cells' and the cells under 10 % solid (warmer than T10); and the total."""
    flowing = sk.valid & np.isfinite(skin_mod.cell_viscosity(sk, 0.0))
    warm = sk.valid & (sk.T > T10)
    deg = np.degrees(theta)
    out = {"mass_total_kg": float(sk.m.sum())}
    for lo, hi in zip(BANDS_DEG[:-1], BANDS_DEG[1:]):
        rows = (deg >= lo) & (deg < hi)
        for name, cells_in in (("all", sk.valid), ("flowing", flowing), ("under10", warm)):
            out["mass_{}_{:g}_{:g}_kg".format(name, lo, hi)] = float((sk.m * cells_in)[rows].sum())
    return out


# -- statistics --------------------------------------------------------------------------------------------------------
def stats(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if not x.size:
        return {"n": 0, "median": None, "p10": None, "p90": None, "max": None}
    p10, med, p90 = np.percentile(x, [10, 50, 90])
    return {"n": int(x.size), "median": float(med), "p10": float(p10), "p90": float(p90), "max": float(x.max())}


def summary(sample, A_top):
    """The JSON summary of one sample: stats of every numeric array (per top cell for the (n, TOP) ones), the gates'
    counts, and per mode the spraying columns and their total release over the windward columns [kg/s]."""
    out = {}
    for key, value in sample.items():
        value = np.asarray(value)
        if value.dtype == object or value.dtype.kind == "U":
            out[key] = dict(collections.Counter(str(v) for v in value))
        elif value.ndim == 2:
            out[key] = [stats(value[:, j]) for j in range(value.shape[1])]
        elif value.ndim == 1:
            out[key] = stats(value)
            if key.endswith("_rate"):
                r = np.nan_to_num(value)
                out[key]["spraying"] = int((r > 0.0).sum())
                out[key]["total_kg_s"] = float((r * A_top).sum())
        else:
            out[key] = float(value)
    return out


def git_commit():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:                                               # noqa: BLE001 -- a record, not a requirement
        return None


def loads_change(new, saved, windward):
    """The recomputed forcing's difference from the saved one, per field (tau, dp_ds): the median and largest absolute
    and relative change, over every column and over the windward ones."""
    if saved is None:
        return None
    old = saved[1]
    out = {}
    for key in ("tau", "dp_ds"):
        a, b = np.asarray(getattr(new, key), dtype=float), np.asarray(getattr(old, key), dtype=float)
        if a.shape != b.shape:
            out[key] = {"shapes": [list(a.shape), list(b.shape)]}
            continue
        d = np.abs(a - b)
        rel = np.divide(d, np.abs(b), out=np.full_like(d, np.nan), where=np.abs(b) > 0.0)
        for where, mask in (("all", np.ones(a.shape, dtype=bool)), ("windward", windward)):
            out["{}_{}".format(key, where)] = {"abs": stats(d[mask]), "rel": stats(rel[mask]),
                                               "new": stats(a[mask]), "saved": stats(b[mask])}
    return out


def solid_temperature(mat, f_s):
    """The temperature [K] at which the material's Scheil curve stands at the solid fraction f_s."""
    g = lambda T: 1.0 - float(np.asarray(mat.liquid_fraction(np.array([T])))[0]) - f_s
    return brentq(g, mat.T_rigid, mat.T_liquidus, xtol=1e-10)


def write(out_dir, name, arrays, record):
    np.savez_compressed(os.path.join(out_dir, name + ".npz"),
                        **{k: (np.asarray(v, dtype=str) if np.asarray(v).dtype == object else np.asarray(v))
                           for k, v in arrays.items()})
    with open(os.path.join(out_dir, name + ".json"), "w") as fh:
        json.dump(record, fh, indent=1)


# -- the single column -------------------------------------------------------------------------------------------------
def single_column(body, flow, idx, out_dir, case, record):
    sk, mat = body.skins, body.material
    _, _, n_flow, n_liq = flow_counts(sk, mat)
    F = sk.depth_above(mat.T_rigid)
    cands = idx[n_flow[idx] > n_liq[idx]]
    if not cands.size:
        raise RuntimeError("no windward column holds slurry cells")
    i = int(cands[np.argsort(F[cands], kind="stable")[len(cands) // 2]])
    G, _ = skin_mod.push(flow.columns, sk.thickness(), flow.dp_ds, flow.g, flow.v_hat, flow.rho_l, flow.sigma,
                         wet=skin_mod.flowing_tops(sk))
    sub = sk.rows([i])
    push, shear, A_top = G[[i]], flow.shear[[i]], flow.A_top[[i]]
    U0 = sub.U.copy()
    started = time.perf_counter()
    change, steps = np.inf, 0
    for steps in range(1, SINGLE_MAX + 1):
        before = sub.U.copy()
        skin_mod.momentum_step(sub, SINGLE_DT, push, shear, A_top)
        top = np.abs(sub.U).max()
        change = float(np.abs(sub.U - before).max() / top) if top > 0.0 else 0.0
        if change < SINGLE_TOL:
            break
    wall = time.perf_counter() - started
    eta, face, carried, gamma = stresses(sub, shear)
    on = sub.valid[0]
    dz = sub.dz()[0]
    arrays = {"depth_mid": (np.cumsum(dz) - 0.5 * dz)[on], "dz": dz[on], "T": sub.T[0][on],
              "fs": (1.0 - sub.props.liquid_fraction(sub.T))[0][on], "speed": np.linalg.norm(sub.U[0], axis=1)[on],
              "speed_start": np.linalg.norm(U0[0], axis=1)[on], "mu_ratio": (eta / mat.liquid.mu)[0][on],
              "stress": carried[0][on], "gamma": gamma[0][on], "face_stress": face[0][:int(sub.count[0]) + 1]}
    record.update({"column": i, "theta_deg": float(np.degrees(body.load_angles()[i])), "flowing_depth_m": float(F[i]),
                   "flowing_cells": int(n_flow[i]), "liquid_cells": int(n_liq[i]), "cells": int(sub.count[0]),
                   "candidates": int(cands.size), "push_N_m3": push[0].tolist(), "shear_Pa": shear[0].tolist(),
                   "steps": steps, "converged": bool(change < SINGLE_TOL), "last_relative_change": change,
                   "wall_s": wall, "top_speed_m_s": float(arrays["speed"][0]),
                   "flow_per_width_m2_s": float((arrays["speed"] * arrays["dz"]).sum())})
    write(out_dir, case + "_single", arrays, record)
    print("{} single column {}: {} sub-steps, converged {}, last change {:.3g}, top speed {:.4g} m/s, {:.1f} s".format(
        case, i, steps, record["converged"], change, record["top_speed_m_s"], wall), flush=True)


# -- the run -----------------------------------------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--state", required=True, help="a saved skin flight (--save-state-at), read only")
    parser.add_argument("--case", required=True, choices=sorted(CASES))
    parser.add_argument("--flow-time", type=float, default=0.1, help="s of flow alone (a whole number of heat sub-steps)")
    parser.add_argument("--out", default=os.path.join(OUTPUT_DIR, "near_liquidus_viscosity", "probe"))
    parser.add_argument("--single-column", action="store_true", help="run the median windward slurry column to steady "
                                                                    "state instead of the flow")
    args = parser.parse_args(argv)
    law, internal = CASES[args.case]
    state = cli.load_state(args.state, parser)
    run = state["run"]
    body = run.body
    if not isinstance(body, skinbody.SkinBody) or not body.flow_on or body.spray is None:
        parser.error("{} is not a flight of flowing skins with Girin's spraying".format(args.state))
    substep = body.settings.substep
    n_heat = int(round(args.flow_time / substep))
    if not (math.isfinite(args.flow_time) and n_heat >= 1 and abs(n_heat * substep - args.flow_time) < 1e-9):
        parser.error("--flow-time must be a positive whole number of the skins' {} s heat sub-steps, got {!r}".format(
            substep, args.flow_time))
    if internal is not None:
        rheology.INTERNAL_SOLID_FRACTION = internal                 # this process only (read at call time)
    body.settings.rheology = law
    body.skins.settings.rheology = law                              # the same object; set both should they ever differ
    os.makedirs(args.out, exist_ok=True)
    mat, sk = body.material, body.skins
    T10 = solid_temperature(mat, SOLID_LIMIT)
    theta = body.load_angles()
    idx = np.flatnonzero(theta < np.radians(WINDWARD_DEG))
    a = run.sim.aero_state(run.sim.t, run.sim.y[:3], run.sim.y[3:])
    started = time.perf_counter()
    forcing = body.flow_forcing(a)
    loads_s = time.perf_counter() - started
    sigma = mat.liquid.sigma if mat.liquid is not None else 0.0
    A_top = body.surface.areas
    flow = skin_mod.Flow(body.columns, forcing.dp_ds, skin_mod.along_meridian(body.columns.normal, body.v_hat, forcing.tau),
                         forcing.g, body.v_hat, sk.rho_l, sigma, A_top, None, release=None, slurry=body.slurry_spray)
    res = forcing.result
    record = {"case": args.case, "law": law, "internal_solid_fraction": rheology.INTERNAL_SOLID_FRACTION,
              "state": os.path.abspath(args.state), "time_s": state["time_s"], "step": state["step"], "commit": git_commit(),
              "T10_K": T10, "T_rigid_K": mat.T_rigid, "T_liquidus_K": mat.T_liquidus, "mu_liquid_Pa_s": mat.liquid.mu,
              "we_critical": body.spray.we_critical, "saved_rules": [body.spray.spray_gates, body.spray.slip_band],
              "windward_deg": WINDWARD_DEG, "windward_columns": int(idx.size), "columns": int(sk.n),
              "altitude_km": float(a.h / 1e3), "speed_m_s": float(a.V), "kn_body": float(res.kn_body),
              "surface_flow_branch": res.branch_name, "in_slip_band": bool(res.in_slip_band),
              "closures_windward": collections.Counter(int(c) for c in np.asarray(res.closure)[idx]),
              "loads_change": loads_change(forcing, body._forcing, theta < np.radians(WINDWARD_DEG)),
              "loads_wall_s": loads_s, "macro_step_s": run.settings.dt}
    record["closures_windward"] = {str(k): v for k, v in record["closures_windward"].items()}
    if args.single_column:
        single_column(body, flow, idx, args.out, args.case, record)
        return 0

    record.update({"flow_time_s": args.flow_time, "heat_substeps": n_heat, "samples": [], "costs": []})
    samples, arrays_t = [], []
    every = max(1, int(round(SAMPLE_EVERY / substep)))

    def sample(t):
        began = time.perf_counter()
        s = modes(body, forcing, a, run.settings.dt, T10, idx)
        s.update(cells(sk, mat, flow, idx))
        s.update(masses(sk, mat, theta, T10))
        s["theta_deg"] = np.degrees(theta[idx])
        samples.append(s)
        arrays_t.append(t)
        summ = summary(s, A_top[idx])
        summ["t_s"], summ["wall_s"] = t, time.perf_counter() - began
        record["samples"].append(summ)
        keys = sorted(s)
        arrays = {"t_s": np.array(arrays_t), "columns": idx, "A_top": A_top[idx]}
        for key in keys:
            arrays[key] = np.stack([np.asarray(x[key]) for x in samples])
        write(args.out, args.case, arrays, record)
        print("{} t = {:.3f} s sampled in {:.1f} s: slurry rate median {} (saved) / {} (recent)".format(
            args.case, t, summ["wall_s"], summ["saved_slurry_rate"]["median"], summ["recent_slurry_rate"]["median"]),
            flush=True)

    sample(0.0)
    seconds0 = dict(skin_mod.FLOW_SECONDS)
    for k in range(1, n_heat + 1):
        began = time.perf_counter()
        fb = skin_mod.flow_substeps(sk, flow, substep)
        sk.regrid()                                                 # the heat sub-step's own regrid follows its flow
        wall = time.perf_counter() - began
        cost = {"heat_substep": k, "flow_substeps": int(fb.sub_steps), "cycles": int(fb.cycles), "wall_s": wall,
                "wall_per_flow_substep_s": wall / fb.sub_steps if fb.sub_steps else None, "courant": float(fb.courant),
                "growth": float(fb.growth), "thickness_passes": int(fb.passes)}
        record["costs"].append(cost)
        print("{} heat sub-step {}/{}: {} flow sub-steps in {:.1f} s (Courant {:.2f})".format(
            args.case, k, n_heat, fb.sub_steps, wall, fb.courant), flush=True)
        if k % every == 0 or k == n_heat:
            sample(k * substep)
    record["flow_seconds"] = {key: skin_mod.FLOW_SECONDS[key] - seconds0.get(key, 0.0) for key in skin_mod.FLOW_SECONDS}
    flows = sum(c["flow_substeps"] for c in record["costs"])
    walls = sum(c["wall_s"] for c in record["costs"])
    record["totals"] = {"flow_substeps": flows, "flow_wall_s": walls, "wall_per_flow_substep_s": walls / flows if flows else None}
    write(args.out, args.case, {**{"t_s": np.array(arrays_t), "columns": idx, "A_top": A_top[idx]},
                                **{key: np.stack([np.asarray(x[key]) for x in samples]) for key in sorted(samples[0])}},
          record)
    print("{} done: {} flow sub-steps, {:.1f} s of flow wall time".format(args.case, flows, walls), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
