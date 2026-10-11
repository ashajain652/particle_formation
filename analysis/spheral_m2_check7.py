"""Spheral M2 check 7: tension behaviour and the choice between SolidSPH and SolidCRKSPH (M2 plan, Task 9; spec §12
check 7, §14.1 conditions 1 and 2).

    "$PY" analysis/spheral_m2_check7.py [--only run|damage|analyse] [--max-load 2.0] [--force]

**The runs** (`CASES`). Checks 4, 5 and 6 rerun with `--scheme crksph` (run.py "Schemes": SolidCRKSPH, integrated
density, mass-over-density volumes; in RZ SolidCRKSPHRZ without compatible energy) with exactly their drivers'
settings, at the spacings the SPH studies already ran:

- check 4 (the heated sphere, `spheral_m2_check4.py`): 10 crossings, damping 0.05, frames 0-3, at RZ 2.2 and 1.1
  mm and 3D 3.0 mm (3D 2.2 mm would cost about 50 min for three intervals);
- check 5 (the bar held at both ends, `spheral_m2_check5.py`): frames 0-10 at RZ 2.2 mm (SolidSPHRZ's collapsed
  rod: does CRKSPH hold it?), RZ 1.1 mm and 3D 2.2 mm;
- check 6 (frame 182's loads, `spheral_m2_check6.py`): `--area-shares patch` (Task 8: the volume shares misplace
  10-13 % of the load on a lattice; passed explicitly, as the SPH patch rows; the default since decision 9),
  RZ 2.2 and 1.1 mm and 3D 2.2 mm, over frames 0-3 instead of 0-6: the SPH runs had settled by k 2 (the stress change
  between frames 2e-4 of the peak), and the CRKSPH 3D step costs about five times SolidSPH's.

- the free uniform expansion ("fx", added after check 4's RZ CRKSPH runs came out stressed): the 100 mm sphere heated
  uniformly by 50 K over the first interval and held for a second (body `--synthetic uniform --uniform-delta-k 50
  --uniform-ramp-frames 1`), elastic, loads off, at RZ 2.2 mm and 3D 3.0 mm with *both* schemes (run here): the patch
  test, exact answer zero stress and an outer radial strain alpha_eff dT; `err_*` the von Mises over E alpha_eff dT,
  `fx_p_mean` the mean pressure over it, `fx_dilatation_ratio` the shell's mean u_r / r over alpha_eff dT. (Its KE
  ratios and monitor estimates are relative to a stress that should be zero: not meaningful.)

The SPH side of every comparison is the existing SPH run of the same case (checks 4-6's run directories, read, not
rerun), measured here by the same functions, so both schemes' rows come from one code path.

**Measures** per run and frame (`measure_*`), each the check's own (see its driver's docstring):
check 4 `err_interior` / `err_shell` / `err_interior_rms` over sigma_0 (binned interior and surface-shell error
against Timoshenko and Goodier, the particle RMS); check 5 `sxx_err_over_Y` (the middle third's axial stress
against fe_ratio (-E alpha_eff dT) capped at -Y), `lat_max_over_Y`, `ps_mean` and `ps_kinematic`; check 6
`err_interior` / `err_interior_rms` / `err_core` / `err_shell` / `err_shell_rms` over VM0 (vM(sigma_p -
sigma_FE)), `drag_rel` (the last step's sum F against the table's integral over the exact sphere), `settle_vm`.
For all: the internal KE ratio, the net momentum over the interval's drag impulse, the re-centring, the run's
crossings and the monitor's estimate, the cost.

**Clumping** (Step 2, `clumping`), at every compared frame (the start is the first frame after K0, the end the
last; at K0 the body is stress free): the particles under tension are those whose largest principal stress
(of the full 3 x 3 sigma, RZ with the hoop element) is positive; over pairs of them closer than the kernel's support
2 h = 2 x 2.01 dx (current positions; RZ in the meridional plane), `tension_min_pair_dx` the smallest pair
distance over the spacing and `tension_share_lt_half` the share of those pairs closer than dx / 2;
`all_min_pair_dx` the smallest over every pair, for context (a lattice starts at 1 dx; RZ rows likewise).

**Damage** (Step 3): `analysis/spheral_m2_check7_damage.py` in the container, both schemes and forms, written to
`spheral_output/m2/check7/damage/<scheme>_<form>.json` and summarised in `data/spheral/m2_check7_damage.csv`.

**Measured** (2026-10-11; `data/spheral/m2_check7.csv`, thresholds in tests/test_spheral_frag_runner_check7.py,
the choice in decision 8 of the M2 plan's "Decisions at review").

- *SolidCRKSPHRZ is not usable at 116c71f.* It cannot step with compatible energy (run.py "Schemes"), and without it
  it fails the patch test: the freely, uniformly heated sphere ends stressed by 1.07 E alpha dT (von Mises, interior)
  and over-expanded (outer radial strain 1.23-1.24 alpha_eff dT) where SolidSPHRZ gives 5.7e-4 and 0.9996. Hence
  check 4's RZ interior errors of 0.38 (2.2 mm) and 0.62 (1.1 mm) sigma_0 against SPH's 0.045 and 0.021, check 5's
  lateral stresses of 0.97-0.99 Y (SPH 1.1 mm 0.018 Y) -- although SolidCRKSPHRZ holds the 2.2 mm rod that collapses
  under SolidSPHRZ -- check 6's RZ 2.2 mm interior error 0.60 VM0 (SPH 0.28), and its RZ 1.1 mm run diverging (1e5
  VM0 by frame 2). Not diagnosed further: the free sphere's over-expansion points at SolidCRKSPHRZ's velocity
  gradient or hoop term (its own S_zz equals -(S_xx + S_rr), so the stress read-out is not the cause).
- *SolidCRKSPH in 3D passes the patch test* (3.0e-4 E alpha dT held, the same as SolidSPH) and runs at the same cost
  per step (51-117 ms on 18 processes). Check 4 (3.0 mm): interior 1.2-1.6 % of sigma_0 against SPH's 2.6-3.6 %, but
  the surface shell 7.1-7.4 % against 5.9-6.2 %: **CRKSPH does not fix the free-surface error.** Check 5 (2.2 mm):
  the elastic middle -0.93 % of Y against -1.46 %, the plateau within 0.13 % of Y to k 6, then drifting off (0.27,
  0.56 and 1.1 % at k 7-9, 8.1 % at k 10 after the monitor estimated 130 % on its pressure and switched to 20 crossings), lateral stresses
  up to 0.35 Y (SPH 0.11 Y, the corners). Check 6 (2.2 mm, patch shares): interior max 0.95 VM0 against 0.65, RMS
  0.058 against 0.040, core 0.042 against 0.040, shell 2.58 against 2.99 VM0 (both: the loads act on one particle
  layer); drag identical (+7.8e-4); its monitor estimated 19.7 % on the first interval (SPH 1.2 %) and switched to 20
  crossings.
- *No clumping in either scheme*: in every run and frame no pair of particles under tension is closer than half a
  spacing (share 0); the smallest tension pair distance is 0.91 dx (CRKSPH RZ 1.1 mm, checks 4 and 6), 0.93 dx
  (CRKSPH 3D bar at k 10), 0.94 dx (SPH RZ 1.1 mm bar at k 10) and 1.0 dx elsewhere.
- *Damage* runs with both schemes in both forms (400 steps each, damage growing: largest 1.0 in SPH RZ and 3D and
  CRKSPH 3D, 0.0028 in CRKSPH RZ), once SpheralController's initial H iteration is off.

**Running it.** Under `caffeinate -i -s`. Before each run on more than 4 processes the driver waits for the VM's
load to fall to half of --max-load (as checks 4-6). Resumable: runs whose summary.json exists are skipped; runs go to
`spheral_output/m2/check7/runs/` in naming.run_dir's layout, logs to `spheral_output/m2/check7/logs/`.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "analysis"))

import spheral_m2_check4 as c4  # noqa: E402
import spheral_m2_check5 as c5  # noqa: E402
import spheral_m2_check6 as c6  # noqa: E402
from spheral_frag import material, naming  # noqa: E402
from spheral_frag.runner import checkpoint, nodes as rn, run  # noqa: E402

PREPARED = c4.PREPARED
OUT = os.path.join("spheral_output", "m2", "check7")
RUNS_OUT = os.path.join(OUT, "runs")
DAMAGE_OUT = os.path.join(OUT, "damage")
CSV = os.path.join("data", "spheral", "m2_check7.csv")
DAMAGE_CSV = os.path.join("data", "spheral", "m2_check7_damage.csv")
LAUNCHER = c4.LAUNCHER
SCHEMES = ("sph", "crksph")
C6_FRAMES = (0, 3)
FX_FRAMES = (0, 2)
FX_DT = 50                                   # K: the free expansion's uniform rise, ramped over one interval
FX_SAMPLER = "uniformdt{}r1".format(FX_DT)
CASES = [("fx", "rz", 2.2, 1), ("fx", "3d", 3.0, 18),
         ("c4", "rz", 2.2, 1), ("c4", "rz", 1.1, 4), ("c4", "3d", 3.0, 18),
         ("c5", "rz", 2.2, 1), ("c5", "rz", 1.1, 4), ("c5", "3d", 2.2, 18),
         ("c6", "rz", 2.2, 1), ("c6", "rz", 1.1, 4), ("c6", "3d", 2.2, 18)]
SUPPORT_DX = 2.0*rn.NPERH                   # the kernel's support 2 h in spacings

COLUMNS = ["check", "scheme", "form", "dx_mm", "k", "processes", "n_particles", "err_interior", "err_interior_rms",
           "err_shell", "err_shell_rms", "err_core", "sxx_err_over_Y", "lat_max_over_Y", "ps_mean", "ps_kinematic",
           "drag_rel", "settle_vm", "fx_p_mean", "fx_dilatation_ratio", "ke_internal_ratio", "momentum_over_impulse", "recentre_shift_dx",
           "max_disp_ref_dx", "run_crossings", "crossings_err_limit", "n_tension", "tension_min_pair_dx",
           "tension_share_lt_half", "n_tension_pairs", "all_min_pair_dx", "steps", "wall_per_step_ms", "run"]
DAMAGE_COLUMNS = ["scheme", "form", "n_particles", "constructed", "stepped", "steps_done", "error", "time_s",
                  "volumetric_strain", "damage_max", "share_damage_gt_0p01", "share_damage_gt_0p99",
                  "wall_per_step_s", "packages"]


# ------------------------------------------------------------------------------------------------- run directories
def frames(check):
    return {"c4": c4.FRAMES, "c5": c5.FRAMES, "c6": C6_FRAMES, "fx": FX_FRAMES}[check]


def brackets(check, scheme):
    if check == "c4":
        b = {"crossings": 10, "damping": 0.05, "strength": "elastic", "sampler": "radial", "loads": "off",
             "mass_scaling": "on"}
    elif check == "c5":
        b = {"strength": c5.STRENGTH, "sampler": "uniform", "loads": "off", "geometry": "bar"}
    elif check == "fx":
        b = {"strength": "elastic", "sampler": FX_SAMPLER, "loads": "off"}
    else:
        b = {"strength": c6.STRENGTH, "sampler": c6.SAMPLER, "loads_frame": c6.LOADS_FRAME, "area_shares": "patch"}
    b["scheme"] = scheme
    return naming.canonical_brackets(b)


def run_name(check, form, dx, scheme):
    return naming.run_name(os.path.basename(PREPARED), "replay", frames(check), form, dx, brackets(check, scheme), 1)


def run_dir(check, form, dx, scheme):
    """CRKSPH (and both schemes of the free expansion): this study's run; SPH: the check's own existing run (check
    6's over frames 0-6)."""
    if scheme == "crksph" or check == "fx":
        return naming.run_dir(os.path.join(REPO, RUNS_OUT), run_name(check, form, dx, scheme))
    if check == "c4":
        return c4.run_dir(form, dx, 10, 0.05, True)
    if check == "c5":
        return c5.run_dir(form, dx)
    return c6.run_dir(form, dx, "patch")


def fx_body_dir(form, dx):
    return naming.body_dir(os.path.join("spheral_output", "body"),
                           naming.body_name(os.path.basename(PREPARED), *FX_FRAMES, form, dx, synthetic=FX_SAMPLER))


def build_fx_bodies():
    for check, form, dx, _ in CASES:
        if check == "fx" and not os.path.isfile(os.path.join(REPO, fx_body_dir(form, dx), "body.json")):
            subprocess.run([sys.executable, "-m", "spheral_frag", "body", "--prepared", PREPARED, "--frames",
                            "{}-{}".format(*FX_FRAMES), "--form", form, "--dx-mm", str(dx), "--synthetic", "uniform",
                            "--uniform-delta-k", str(FX_DT), "--uniform-ramp-frames", "1", "--quiet"],
                           cwd=REPO, check=True)


def body_dir(check, form, dx):
    if check == "fx":
        return fx_body_dir(form, dx)
    return {"c4": c4.body_dir, "c5": c5.body_dir, "c6": c6.body_dir}[check](form, dx)


def command(check, form, dx, out=RUNS_OUT, fr=None, max_load=2.0, scheme="crksph"):
    """run.py's arguments for the CRKSPH run of a case (the check's own settings + --scheme crksph)."""
    fr = fr or frames(check)
    base = ["spheral_frag/runner/run.py", "replay", "--body", body_dir(check, form, dx), "--prepared", PREPARED,
            "--frames", "{}-{}".format(*fr), "--ke-fraction", "1e300", "--out", out, "--max-load", repr(max_load),
            "--scheme", scheme]
    if check == "fx":
        return base + ["--strength", "elastic", "--loads", "off"]
    if check == "c4":
        return base + ["--crossings", "10", "--max-crossings", "10", "--damping", "0.05", "--strength", "elastic",
                       "--loads", "off", "--mass-scaling", "on"]
    if check == "c5":
        return base + ["--strength", c5.STRENGTH, "--loads", "off"]
    return base + ["--strength", c6.STRENGTH, "--loads-frame", str(c6.LOADS_FRAME), "--area-shares", "patch"]


def launch(check, form, dx, procs, max_load, force, scheme="crksph"):
    d = run_dir(check, form, dx, scheme)
    if os.path.isfile(os.path.join(d, "summary.json")):
        return "done"
    if procs > 4 and not force:
        c4.wait_for_load(0.5*max_load)
    cmd = [LAUNCHER] + (["-n", str(procs)] if procs > 1 else []) + command(check, form, dx, max_load=max_load,
                                                                            scheme=scheme)
    if os.path.isdir(d):
        cmd.append("--resume")
    if force:
        cmd.append("--force")
    t = time.perf_counter()
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    os.makedirs(os.path.join(REPO, OUT, "logs"), exist_ok=True)
    with open(os.path.join(REPO, OUT, "logs", "{}_{}_{}_dx{}.log".format(check, scheme, form, dx)), "w") as fh:
        fh.write(" ".join(cmd) + "\n" + r.stdout + "\n--- stderr\n" + r.stderr)
    if r.returncode != 0:
        print("  FAILED ({}): {}".format(r.returncode, r.stderr.strip().splitlines()[-1:] if r.stderr else ""),
              flush=True)
        return "failed"
    return "ran in {:.0f} s".format(time.perf_counter() - t)


def run_damage():
    for scheme in SCHEMES:
        for form in ("rz", "3d"):
            if os.path.isfile(os.path.join(REPO, DAMAGE_OUT, "{}_{}.json".format(scheme, form))):
                continue
            subprocess.run([LAUNCHER, "analysis/spheral_m2_check7_damage.py", "--scheme", scheme, "--form", form,
                            "--out", DAMAGE_OUT], cwd=REPO, capture_output=True, text=True)


# ------------------------------------------------------------------------------------------------------ measuring
def principal_max(sigma):
    return np.linalg.eigvalsh(0.5*(sigma + np.transpose(sigma, (0, 2, 1))))[:, -1]


def clumping(x, sigma, dx):
    """(n under tension, smallest tension pair distance / dx, share of tension pairs < dx/2, n tension pairs,
    smallest pair distance of all / dx): module docstring, "Clumping"."""
    from scipy.spatial import cKDTree
    x = np.asarray(x, dtype=np.float64)
    tree = cKDTree(x)
    dd, _ = tree.query(x, k=2)
    all_min = float(dd[:, 1].min())/dx
    ten = principal_max(sigma) > 0.0
    nt = int(ten.sum())
    if nt < 2:
        return nt, math.nan, math.nan, 0, all_min
    xt = x[ten]
    pairs = cKDTree(xt).query_pairs(SUPPORT_DX*dx, output_type="ndarray")
    if not len(pairs):
        return nt, math.nan, math.nan, 0, all_min
    d = np.linalg.norm(xt[pairs[:, 0]] - xt[pairs[:, 1]], axis=1)/dx
    return nt, float(d.min()), float(np.mean(d < 0.5)), int(len(d)), all_min


def interval_values(iv):
    def num(c):
        v = iv.get(c, "")
        return float(v) if v not in ("", None) else math.nan
    return dict(ke_internal_ratio=num("ke_internal_ratio"), momentum_over_impulse=num("momentum_over_impulse"),
                recentre_shift_dx=num("recentre_shift_dx"), max_disp_ref_dx=num("max_disp_ref_dx"),
                run_crossings=int(num("crossings") if math.isnan(num("run_crossings")) else num("run_crossings")),
                crossings_err_limit=num("crossings_err_limit"),
                steps=int(num("steps")), wall_per_step_ms=1e3*num("wall_per_step_s"))


class Context:
    def __init__(self):
        self.table = material.MaterialTable.load(os.path.join(REPO, PREPARED))
        self.c4 = c4.constants(self.table)                       # E, nu, alpha_eff
        self.c5 = c5.constants(self.table)
        self.c5_fe = c5.fe_reference(self.c5[0], self.c5[1])
        fes = c6.fe_summary()
        self.drag_ref = float(fes["tables"]["lee"]["continuum_drag_N"])
        self.fine = c6.reference(c6.FE_FINE)


def measure_run(ctx, check, form, dx_mm, scheme):
    d = run_dir(check, form, dx_mm, scheme)
    if not os.path.isfile(os.path.join(d, "summary.json")):
        print("  missing: {}".format(naming.run_name_from_dir(d)))
        return []
    with open(os.path.join(d, "summary.json")) as fh:
        summ = json.load(fh)
    ivs = {int(r["k1"]): r for r in run.read_intervals(os.path.join(d, "intervals.csv"))}
    dx = dx_mm*1e-3
    k0, k1 = frames(check)
    rows, prev = [], None
    if check == "c6":
        _, prev = c6.state(d, k0, ctx.table, form)
    for k in range(k0 + 1, k1 + 1):
        row = {c: math.nan for c in COLUMNS}
        row.update(check=check, scheme=scheme, form=form, dx_mm=dx_mm, k=k, processes=summ["processes"],
                   run=naming.run_name_from_dir(d), **interval_values(ivs[k]))
        if check == "fx":
            E, nu, alpha = ctx.c4
            ck = checkpoint.read(os.path.join(d, "checkpoints", checkpoint.filename(k)))
            st = run.stress_state(ck, ctx.table)
            sref = E*alpha*FX_DT
            p3 = c6.points3(form, st["x_ref"])
            r = np.linalg.norm(p3, axis=1)
            rc = np.linalg.norm(c6.points3(form, st["x"]), axis=1)
            inner, outer = r <= c4.R - c4.SHELL_DX*dx, r > c4.R - c4.SHELL_DX*dx
            vm = st["von_mises"]/sref
            row.update(err_interior=float(vm[inner].max()), err_shell=float(vm[outer].max()),
                       err_interior_rms=c6.wrms(vm[inner], st["V"][inner]),
                       fx_p_mean=float(np.sum(st["V"]*st["P"])/np.sum(st["V"]))/sref,
                       fx_dilatation_ratio=float(np.mean((rc[outer] - r[outer])/r[outer]))/(alpha*FX_DT*min(k, 1)))
        elif check == "c4":
            E, nu, alpha = ctx.c4
            mm = c4.measure(d, k, ctx.table, E, nu, alpha, dx)
            st = mm["st"]
            row.update(err_interior=mm["interior"][0], err_shell=mm["shell"][0], err_interior_rms=mm["rms"])
        elif check == "c5":
            E, nu, alpha = ctx.c5
            m = c5.measure(d, k, ctx.table, E, alpha, dx, form)
            ratio = ctx.c5_fe["{}_{}".format(form, dx_mm)]["ratio"]
            dT = c5.DT_PER_FRAME*k
            sxx_exp = -min(ratio*E*alpha*dT, c5.Y)
            row.update(sxx_err_over_Y=(m["sxx_mean"] - sxx_exp)/c5.Y, lat_max_over_Y=m["lat_max"]/c5.Y,
                       ps_mean=m["ps_mean"], ps_kinematic=m["ps_kin"])
            st = run.stress_state(checkpoint.read(os.path.join(d, "checkpoints", checkpoint.filename(k))), ctx.table)
        else:
            m = c6.measure(d, k, ctx.table, form, ctx.fine, None, prev)
            st = prev = m.pop("_state")
            m.pop("_fe")
            h = np.load(os.path.join(d, "history", "interval_{:05d}.npz".format(k - 1)))["history"]
            row.update(err_interior=m["int_vm_max"], err_interior_rms=m["int_vm_rms"], err_core=m["core_vm_max"],
                       err_shell=m["shell_vm_max"], err_shell_rms=m["shell_vm_rms"], settle_vm=m["settle_vm"],
                       drag_rel=-float(h[-1, 1])/ctx.drag_ref - 1.0)
        nt, tmin, tshare, npairs, amin = clumping(st["x"], st["sigma"], dx)
        row.update(n_particles=len(st["pid"]), n_tension=nt, tension_min_pair_dx=tmin, tension_share_lt_half=tshare,
                   n_tension_pairs=npairs, all_min_pair_dx=amin)
        rows.append(row)
    return rows


def write_csv(path, rows, columns):
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(columns)
        for row in rows:
            wr.writerow([repr(float(row[c])) if isinstance(row[c], (float, np.floating)) else row[c]
                         for c in columns])


def analyse():
    ctx = Context()
    rows = []
    for check, form, dx, _ in CASES:
        for scheme in SCHEMES:
            rows += measure_run(ctx, check, form, dx, scheme)
    write_csv(os.path.join(REPO, CSV), rows, COLUMNS)
    print("wrote {} ({} rows)".format(CSV, len(rows)))
    drows = []
    for scheme in SCHEMES:
        for form in ("rz", "3d"):
            p = os.path.join(REPO, DAMAGE_OUT, "{}_{}.json".format(scheme, form))
            if not os.path.isfile(p):
                print("  missing damage probe {} {}".format(scheme, form))
                continue
            with open(p) as fh:
                j = json.load(fh)
            r = {c: j.get(c, math.nan) for c in DAMAGE_COLUMNS}
            r["packages"] = "+".join(j.get("packages", []))
            r["error"] = (j.get("error") or "").replace("\n", " ")[:300]
            for c in ("constructed", "stepped"):
                r[c] = int(bool(j.get(c)))
            drows.append(r)
    write_csv(os.path.join(REPO, DAMAGE_CSV), drows, DAMAGE_COLUMNS)
    print("wrote {} ({} rows)".format(DAMAGE_CSV, len(drows)))
    report(rows)
    return rows


def report(rows):
    print("\n{:3s} {:6s} {:3s} {:4s} {:>2s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>6s} {:>8s} {:>8s} {:>8s} "
          "{:>7s}".format("chk", "scheme", "frm", "dx", "k", "interior", "int_rms", "shell", "sxx/Y", "drag_rel",
                          "KEint/E", "nten", "tmin_dx", "t<half", "allmin", "ms/stp"))
    for r in rows:
        print("{:3s} {:6s} {:3s} {:4.1f} {:2d} {:9.3e} {:9.3e} {:9.3e} {:9.2e} {:9.2e} {:9.2e} {:6d} {:8.3f} {:8.2e} "
              "{:8.3f} {:7.1f}".format(r["check"], r["scheme"], r["form"], r["dx_mm"], r["k"], r["err_interior"],
                                       r["err_interior_rms"], r["err_shell"], r["sxx_err_over_Y"], r["drag_rel"],
                                       r["ke_internal_ratio"], r["n_tension"], r["tension_min_pair_dx"],
                                       r["tension_share_lt_half"], r["all_min_pair_dx"], r["wall_per_step_ms"]))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--only", choices=("run", "damage", "analyse"))
    p.add_argument("--max-load", type=float, default=2.0)
    p.add_argument("--force", action="store_true", help="start runs even above --max-load")
    a = p.parse_args(argv)
    os.chdir(REPO)
    if a.only in (None, "run"):
        build_fx_bodies()
        for check, form, dx, procs in CASES:
            for scheme in (SCHEMES if check == "fx" else ("crksph",)):
                print("{} {} {} mm {}: ".format(check, form, dx, scheme), end="", flush=True)
                print(launch(check, form, dx, procs, a.max_load, a.force, scheme), flush=True)
    if a.only in (None, "damage"):
        run_damage()
    if a.only in (None, "analyse"):
        analyse()
    return 0


if __name__ == "__main__":
    sys.exit(main())
