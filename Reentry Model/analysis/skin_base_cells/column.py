#!/usr/bin/env python3
"""The single-column test of a thick skin's base (open item 7, 2026-11-14; Asha: "Run the single-column test with coarse
and fine bottom cells", then "Go with option 1, test it on the column first").

Why: the jump-fault probe of the 50 mm flight (hand-down off, 178.5-186.5 s) found the bottom cell of thick (about 28 mm)
skins about 10 K above the rigid-front metal crossing their base, the skin's bottom two cells at the same temperature. A
skin holds its base at the 3D model's temperature through its bottom half-cell (`skin._conduct`), so the step is the drop
across that half-cell, and graded from the top alone a 28 mm skin's bottom cell is 2-4 mm thick.

The problem, in one dimension: a pool of metal `pool_mm` thick at `T_pool` on rigid metal `sub_mm` deep at `T_sub`,
AA7075 with Scheil's curve at one density (the skins' `swell=False`), the top heated at `q` (adiabatic for 0, no
radiation), the bottom adiabatic. Case A: 906 K on 894 K, no heating; case B: 906 K on 885 K, heated at 2 MW/m2.

- `skins`: the pool is one real skin column (`reentry_model.skyscrapers.skin`: its conduction, split and merge, trial and
  real passes), graded by `growth` and the base cell `base_cell_mm` (None: the settings' default, 0: none); the rigid
  metal is a uniform grid of `dx_um` cells solved here (backward Euler in the enthalpy, Newton, banded). Every coupling
  step `dt_ms` (the flight's 12.5 ms) follows `skinbody`'s damped rule: the skin's trial pass at a base temperature gives
  the base law E/dt = a + b T_b, the rigid metal takes it as its top flux, the real pass holds the base at the rigid
  metal's new top temperature, and the pass repeats re-linearised part of the way there (Aitken's factor) until the two
  base energies agree within `skin.INTERFACE_TOLERANCE` or the base moves at most `DAMPED_TOL_K`, at most
  `DAMPED_MAX_PASSES` times; what is left is a load on the rigid metal's top in the next step. The rigid metal's top cells
  at or above T_rigid then join the skin's bottom cell with their mass and enthalpy (the metal crossing the rigid front).
  Mass and energy close exactly.
- `reference`: the whole column on one uniform grid of `dx_um` cells, sharing nothing with the skins but the material.

Reported every `every` seconds: the rigid front's advance below the original interface (skins: the rigid metal taken
in; reference: the deepest metal at or above T_rigid) and, for the skins, the bottom cell's temperature and thickness
and its step above the base.

Results (2026-11-14, `compare`): the reference is converged (5 and 10 um, 1 and 0.5 ms agree to 5 um). Graded from the
top alone (growth 1.10, bottom cell 2.2-4.3 mm) the front is 17-64 % too far in case A and, in case B, stuck for 0.1 s,
then 28 % too far at 2 s, with steps of 7-11 K at first; the 12.5 and 2.5 ms coupling steps agree, so the error is the
cell's. With 10 um cells throughout the skins reproduce the reference to its grid. With a base cell of 0.1, 0.2, 0.3 and
0.5 mm the front is within 0.3, 1.2, 2.3 and 6.2 % from 0.5 s on (case A; case B within 0.8, 0.8, 1.5 and 0.5 %), at
130-137, 120-126, 114-123 and 110-115 cells (case A) against 100-103 without; 0.2 mm is the default.

    python -I -B analysis/skin_base_cells/column.py skins --out output/skin_base_cells/A_g110.json
    python -I -B analysis/skin_base_cells/column.py skins --base-cell-mm 0 --out output/skin_base_cells/A_g110_off.json
    python -I -B analysis/skin_base_cells/column.py reference --dx-um 5 --dt-ms 1 --out output/skin_base_cells/A_ref_5um_1ms.json
    (case B: add --T-sub 885 --q 2e6 and name the files B_...; then python -I -B analysis/skin_base_cells/compare.py)
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import solve_banded

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from reentry_model.shared import material  # noqa: E402
from reentry_model.skyscrapers import skin  # noqa: E402
from reentry_model.skyscrapers import skinbody as sb  # noqa: E402

MAT = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
T_RIGID = MAT.T_rigid
RHO = MAT.rho
PROFILE_STEP = 20e-6              # m, the spacing of the reference's saved temperature profile


def newton_column(T, H_old, m, dt, dx, props, a=0.0, b=0.0, q_top=0.0, tol=1e-7, max_it=60):
    """One backward-Euler step of a uniform column of cells (mass m per m2 each, thickness dx), top flux a + b T[0] + q_top
    [W/m2], bottom adiabatic; `props(T)` -> (h, cp, k). Newton in T with a backtracking line search."""
    def resid(T):
        h, cp, k = props(T)
        G = 1.0 / (0.5 * dx / k[:-1] + 0.5 * dx / k[1:])
        flux = G * (T[:-1] - T[1:])
        net = np.zeros_like(T)
        net[:-1] -= flux
        net[1:] += flux
        net[0] += a + b * T[0] + q_top
        return (m * h - H_old) / dt - net, cp, G
    R, cp, G = resid(T)
    for _ in range(max_it):
        diag = m * cp / dt
        diag[:-1] += G
        diag[1:] += G
        diag[0] -= b
        ab = np.zeros((3, len(T)))
        ab[0, 1:] = -G
        ab[1] = diag
        ab[2, :-1] = -G
        dT = solve_banded((1, 1), ab, -R)
        r0 = np.linalg.norm(R)
        lam = 1.0
        for _ in range(20):
            T_try = T + lam * dT
            R_try, cp_try, G_try = resid(T_try)
            if np.linalg.norm(R_try) <= (1.0 - 1e-4 * lam) * r0 or lam * np.abs(dT).max() <= tol:
                break
            lam *= 0.5
        T, R, cp, G = T_try, R_try, cp_try, G_try
        if lam * np.abs(dT).max() <= tol:
            return T
    raise RuntimeError("column step did not converge (last change {:.2e} K)".format(lam * np.abs(dT).max()))


def material_props(T):
    return MAT.enthalpy(T), MAT.cp_eff(T), MAT.k(T)


def run_reference(a, log=print):
    """The whole column on one grid; rows of {t, front_mm, profile}."""
    dx = a.dx_um * 1e-6
    n_pool, n_sub = int(round(a.pool_mm * 1e-3 / dx)), int(round(a.sub_mm * 1e-3 / dx))
    T = np.concatenate([np.full(n_pool, a.T_pool), np.full(n_sub, a.T_sub)])
    m = RHO * dx
    H = m * MAT.enthalpy(T)
    dt, t, rows = a.dt_ms * 1e-3, 0.0, []
    next_report, every = a.every, max(1, int(round(PROFILE_STEP / dx)))
    t0 = time.time()
    while t < a.t_end - 1e-12:
        T = newton_column(T, H, m, dt, dx, material_props, q_top=a.q)
        H = m * MAT.enthalpy(T)
        t += dt
        if t >= next_report - 1e-9:
            nonrigid = np.flatnonzero(T >= T_RIGID)
            deepest = (nonrigid.max() + 1) * dx if nonrigid.size else 0.0
            rows.append({"t": round(t, 6), "front_mm": (deepest - a.pool_mm * 1e-3) * 1e3,
                         "profile": {"z_mm": ((np.arange(len(T)) + 0.5) * dx * 1e3)[::every].tolist(),
                                     "T": T[::every].tolist()}})
            next_report += a.every
            log("t {:.3f} s  front {:+.4f} mm  ({:.0f} s)".format(t, rows[-1]["front_mm"], time.time() - t0))
    return rows


def build_skin(a):
    """One skin column holding the pool at T_pool, its cells as the settings' targets make them."""
    extra = {} if a.base_cell_mm is None else {"base_cell": a.base_cell_mm * 1e-3 if a.base_cell_mm > 0.0 else None}
    s = skin.SkinSettings(growth=a.growth, swell=False, **extra)
    total = a.pool_mm * 1e-3
    edges, z = [0.0], 0.0
    while z < total - 1e-12:                         # graded from the top; the split and merge then apply the base's rule
        t = s.cell + (s.growth - 1.0) * max(z - s.fine_zone, 0.0)
        z = min(z + t, total)
        if total - z < 0.5 * t:
            z = total
        edges.append(z)
    dz = np.diff(edges)
    sk = skin.Skins(MAT, [1.0], s)
    m = RHO * dz[None, :]
    sk._install(m, m * MAT.enthalpy(np.full(dz.shape, a.T_pool))[None, :], np.array([len(dz)]))
    sk.regrid()
    return sk


def run_skins(a, log=print):
    """The pool as a skin on a grid of rigid metal; (rows, closure)."""
    sk = build_skin(a)
    props = sk.props
    dx = a.dx_um * 1e-6
    n_sub = int(round(a.sub_mm * 1e-3 / dx))
    m = RHO * dx
    T_s = np.full(n_sub, a.T_sub)
    H_s = m * props.enthalpy(T_s)
    top = 0                                            # the first cell of rigid metal not yet taken in
    E0, M0 = sk.energy() + H_s.sum(), sk.mass() + m * n_sub
    q = np.array([a.q])
    dt, t, rows = a.dt_ms * 1e-3, 0.0, []
    next_report = a.every
    mismatch_max = E_top_total = deferred = 0.0      # deferred: J/m2 of last step's mismatch, a load on the top cell
    passes_total = steps = unconverged = 0
    t0 = time.time()
    while t < a.t_end - 1e-12:
        saved = sk.state()
        x, omega, r_prev = float(T_s[top]), sb.DAMPED_OMEGA0, None
        for p in range(sb.DAMPED_MAX_PASSES):
            sk.set_state(saved)
            law_a, law_b = skin.advance(sk, dt, q, np.array([1.0]), 0.0, 0.0, np.array([x]), 0.0, remove=False)
            T_new = newton_column(T_s[top:].copy(), H_s[top:], m, dt, dx, props.conduction,
                                  a=float(law_a[0]) + deferred / dt, b=float(law_b[0]))
            H_new = m * props.enthalpy(T_new)
            E_in = H_new.sum() - H_s[top:].sum() - deferred
            book = skin.advance(sk, dt, q, np.array([1.0]), 0.0, 0.0, np.array([x]), 0.0,
                                T_b1=np.array([T_new[0]]), remove=False)
            E_base = float(book.E_base[0])
            mis, r = E_base - E_in, float(T_new[0]) - x
            ok = abs(mis) <= skin.INTERFACE_TOLERANCE * (abs(float(book.E_top[0])) + abs(E_base))
            if ok or abs(r) <= sb.DAMPED_TOL_K or p == sb.DAMPED_MAX_PASSES - 1:
                break
            if r_prev is not None and r != r_prev:
                omega = float(np.clip(-omega * r_prev / (r - r_prev), *sb.DAMPED_OMEGA_RANGE))
            x, r_prev = x + omega * r, r
        passes_total, steps = passes_total + p + 1, steps + 1
        unconverged += int(not ok and abs(r) > sb.DAMPED_TOL_K)
        T_s[top:], H_s[top:] = T_new, H_new
        mismatch_max = max(mismatch_max, abs(mis) / max(abs(E_base), 1e-30))
        deferred = mis
        E_top_total += float(book.E_top[0])
        k = 0
        while top + k < n_sub and T_s[top + k] >= T_RIGID:
            k += 1
        if k:
            dm, dH = m * k, H_s[top:top + k].sum()
            sk.add_bottom(np.array([dm]), np.array([dH / dm]))
            sk.regrid()
            top += k
        t += dt
        if t >= next_report - 1e-9:
            b = sk.count[0] - 1
            dz = sk.dz()[0]
            rows.append({"t": round(t, 6), "front_mm": top * dx * 1e3, "T_bottom": float(sk.T[0, b]),
                         "T_above": float(sk.T[0, b - 1]), "T_base": float(T_s[top]),
                         "step_K": float(sk.T[0, b] - T_s[top]), "dz_bottom_mm": float(dz[b] * 1e3),
                         "dz_above_mm": float(dz[b - 1] * 1e3), "cells": int(sk.count[0]),
                         "pool_mm": float(dz[:b + 1].sum() * 1e3)})
            next_report += a.every
            log("t {:.3f} s  front {:.4f} mm  bottom {:.2f} K ({:.3f} mm)  base {:.2f} K  step {:.2f} K  cells {}  "
                "({:.0f} s)".format(t, rows[-1]["front_mm"], rows[-1]["T_bottom"], rows[-1]["dz_bottom_mm"],
                                    rows[-1]["T_base"], rows[-1]["step_K"], rows[-1]["cells"], time.time() - t0))
    E1 = sk.energy() + H_s[top:].sum() + deferred     # the last mismatch is still owed to the rigid metal
    M1 = sk.mass() + m * (n_sub - top)
    closure = {"mass_rel": (M1 - M0) / M0, "energy_rel": float((E1 - E0 - E_top_total) / abs(E0)),
               "interface_mismatch_max_rel": float(mismatch_max), "passes_per_step": passes_total / max(steps, 1),
               "steps_unconverged": unconverged, "steps": steps}
    return rows, closure


def build_parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("model", choices=("skins", "reference"))
    p.add_argument("--growth", type=float, default=1.10)
    p.add_argument("--base-cell-mm", type=float, default=None, help="the skin's base cell (default: the settings'; 0: none)")
    p.add_argument("--dx-um", type=float, default=5.0)
    p.add_argument("--dt-ms", type=float, default=12.5)
    p.add_argument("--pool-mm", type=float, default=28.0)
    p.add_argument("--sub-mm", type=float, default=30.0)
    p.add_argument("--T-pool", type=float, default=906.0)
    p.add_argument("--T-sub", type=float, default=894.0)
    p.add_argument("--q", type=float, default=0.0)
    p.add_argument("--t-end", type=float, default=3.0)
    p.add_argument("--every", type=float, default=0.1)
    p.add_argument("--out", default=None)
    return p


def options(model, **overrides):
    """The command line's defaults for `model`, with `overrides` (keyword names as the attributes: base_cell_mm=0.0)."""
    a = build_parser().parse_args([model])
    for key, value in overrides.items():
        if not hasattr(a, key):
            raise AttributeError(key)
        setattr(a, key, value)
    return a


def main(argv=None):
    a = build_parser().parse_args(argv)
    if a.out is None:
        raise SystemExit("--out is needed")
    t0 = time.time()
    print("T_rigid {:.3f} K, liquidus {:.1f} K, args {}".format(T_RIGID, MAT.T_liquidus, vars(a)), flush=True)
    if a.model == "reference":
        rows, closure = run_reference(a), None
    else:
        rows, closure = run_skins(a)
        print("closure", closure, flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump({"args": vars(a), "T_rigid": T_RIGID, "rows": rows, "closure": closure, "wall_s": time.time() - t0}, fh)
    print("done in {:.0f} s -> {}".format(time.time() - t0, a.out), flush=True)


if __name__ == "__main__":
    main()
