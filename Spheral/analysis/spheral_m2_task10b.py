"""Spheral M2 Task 10b: a remedy for SolidSPH's tensile pairing on the MVP flight's replay (check 13's finding 1).

    "$PY" analysis/spheral_m2_task10b.py [--only probe|replay|checks|analyse] [--max-load 2.0] [--force]

**The reproducer** (`PROBES`). Task 10 found the core of the MVP flight's replay in 40-50 MPa axial tension and
SolidSPHRZ's particles there pairing up axially at RZ 1.1 mm: whole lattice planes x = const shift alternately
along the axis (the same sign at every r out to about 33 mm; `pairing`, check 13's measure: the largest second
difference of the axial displacement along a lattice line, in spacings), growing x3.2 per interval at 10 crossings
from k 3, and the core's axial stress lost (`core_sxx`, the volume-weighted mean over |x_ref| <= 25 mm). The
reproducer reruns check 13's RZ 1.1 mm interval 4 -> 5 from its checkpoint (pairing 0.0029 dx at its start) with each
candidate at 10, 20 and 40 crossings on 4 processes (about 20 s per 10 crossings; analysis/spheral_m2_task10b_probe.py):
the pairing's growth factor over the interval and the core's axial stress at k 5 against RZ 2.2 mm's and 3D's
(38.09 and 36.52 MPa, check 13's runs, where the pairing is weak). Confirmed on RZ 2.2 mm's interval 39 -> 40 and 3D
2.2 mm's 39 -> 40 (after check 13's switches, where the pairing grew at 2.2 mm), and by replays from k 0 (`REPLAYS`).

**Candidates** (`PROBES`, measured 2026-10-12; the CSV):
- the runner as built: WendlandC4 (Spheral's Wendland kernels have extent 1: the support is h) at nPerh 2.01, h
  ~2.2 dx after IdealH, so every neighbour lies beyond the kernel's inflection point (q 0.235 against the nearest
  neighbour's 0.45-0.50): the regime of Swegle's tensile instability. x3.20, 7.85, 17.2 at 10, 20, 40 crossings.
- (a) kernel and neighbours: WendlandC6 at 2.01 x8.9 (worse, and rings: internal KE/E 0.37); nPerh 3.01 x0.69, 0.46,
  0.65 (stable; 4.01 x0.56) but +25-35 % per step in RZ and about (3/2)^3 the neighbours in 3D, and it moves the core
  stress (36.9 MPa against 38.1); WendlandC6 at 3.01 x0.93; WendlandC2 at 3.01 likewise stable.
- (b) Spheral's Monaghan artificial stress (SPH's epsTensile, nTensile 4): at 0.3 Spheral's step falls 21 times
  (7.0e-7 s against 1.5e-5 s over one crossing) and each step costs 1.5 times more, the body thrown into motion
  (internal KE/E 1.9 after one crossing) -- Task 10's "beyond 20 minutes" is a ~30 times slower interval; 0.1 still
  4.5 times the steps; at 0.02 and 0.05 the step holds but the pairing grows faster (x3.9 and x9.7) and the core's
  axial stress falls (35.2 and 7.9 MPa): no affordable working setting.
- (c) damping: the replay's damping at zeta 0.2 x3.14 (it acts at the body's lowest frequency); XSPH x1.66;
  the unlimited Monaghan-Gingold viscosity x3.45. A growing exponential mode is only slowed by velocity damping.
- (d) the hourglass stiffness (packages/hourglass.py, run.py --hourglass KAPPA): a stiffness against the non-affine
  part of each particle's displacement from its reference position, with a viscous part on the same mode; zero on
  every affine field, momentum-conserving. kappa 1e-4 x3.11, 1e-3 x2.4/4.6/11, 3e-3 x1.3/1.45/1.6 (slowed, still
  growing), 1e-2 x0.15/0.08/0.08 (removed) with the core stress 38.07-38.16 MPa; stiffness alone at 1e-2 lets an
  axial oscillation of the two rows next to the axis grow at 40 crossings (internal KE/E 9.9e-3), which the viscous
  part (beta 1) removes (9.2e-6), x0.18/0.08/0.08.

**Checks 4-6 with the hourglass** (`CHECKS`): the cheapest representative runs of each rerun with --hourglass 0.01 and
measured with the checks' own functions against their committed rows (`data/m2_task10b_checks.csv`).

Outputs: `data/m2_task10b_probe.csv` (one row per probe and crossings number, and the epsTensile timings),
`m2_task10b_replay.csv` (the hourglass replays against check 13's runs per frame), `m2_task10b_checks.csv`; runs and
probe states in `spheral_output/m2_probe/task10b/` (git-ignored)."""
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

from spheral_frag import material, naming  # noqa: E402
from spheral_frag.runner import checkpoint, run  # noqa: E402

import spheral_m2_check13 as c13  # noqa: E402

OUT = os.path.join("spheral_output", "m2_probe", "task10b")
PROBE_SCRIPT = os.path.join("analysis", "spheral_m2_task10b_probe.py")
LAUNCHER = os.path.join("spheral_frag", "container", "spheral")
PROBE_CSV = os.path.join("data", "m2_task10b_probe.csv")
REPLAY_CSV = os.path.join("data", "m2_task10b_replay.csv")
CHECKS_CSV = os.path.join("data", "m2_task10b_checks.csv")
KAPPA = 0.01                                   # the hourglass stiffness chosen (module docstring)

# name, form, dx, processes, interval start k, crossings, config (analysis/spheral_m2_task10b_probe.py)
PROBES = [
    ("base", "rz", 1.1, 4, 4, (10, 20, 40), {}),
    ("nperh3", "rz", 1.1, 4, 4, (10, 20, 40), {"nperh": 3.01}),
    ("nperh4", "rz", 1.1, 4, 4, (10,), {"nperh": 4.01}),
    ("c6", "rz", 1.1, 4, 4, (10,), {"kernel": "C6"}),
    ("c6_nperh3", "rz", 1.1, 4, 4, (10,), {"kernel": "C6", "nperh": 3.01}),
    ("c2_nperh3", "rz", 1.1, 4, 4, (10,), {"kernel": "C2", "nperh": 3.01}),
    ("xsph", "rz", 1.1, 4, 4, (10,), {"xsph": True}),
    ("q_mg", "rz", 1.1, 4, 4, (10,), {"Q": "mg"}),
    ("zeta0.2", "rz", 1.1, 4, 4, (10,), {"zeta": 0.2}),
    ("eps0.02", "rz", 1.1, 4, 4, (10, 40), {"eps": 0.02}),
    ("eps0.05", "rz", 1.1, 4, 4, (10,), {"eps": 0.05}),
    ("hg1e-4", "rz", 1.1, 4, 4, (10,), {"hg": 1e-4, "beta": 0.0}),
    ("hg1e-3", "rz", 1.1, 4, 4, (10, 20, 40), {"hg": 1e-3, "beta": 0.0}),
    ("hg3e-3", "rz", 1.1, 4, 4, (10, 20, 40), {"hg": 3e-3, "beta": 0.0}),
    ("hg1e-2_stiff", "rz", 1.1, 4, 4, (10, 20, 40), {"hg": 1e-2, "beta": 0.0}),
    ("hg1e-2", "rz", 1.1, 4, 4, (10, 20, 40), {"hg": 1e-2}),
    ("rz22_k39_base", "rz", 2.2, 1, 39, (10, 40), {}),
    ("rz22_k39_hg1e-2", "rz", 2.2, 1, 39, (10, 40), {"hg": 1e-2}),
    ("3d22_k39_base", "3d", 2.2, 18, 39, (20,), {}),
    ("3d22_k39_hg1e-2", "3d", 2.2, 18, 39, (20,), {"hg": 1e-2}),
]
TIMINGS = [("eps0", {}), ("eps0.02", {"eps": 0.02}), ("eps0.1", {"eps": 0.1}), ("eps0.3", {"eps": 0.3})]
# check 13's probes whose states are the baseline here (same interval, same configuration): reused, not rerun
REUSE = {"rz22_k39_base": os.path.join("spheral_output", "m2_probe", "task10", "rz22"),
         "3d22_k39_base": os.path.join("spheral_output", "m2_probe", "task10", "3d22")}
REPLAYS = [("rz", 1.1, 4, 12), ("rz", 2.2, 1, 12), ("3d", 2.2, 18, 3)]       # form, dx, processes, last frame

PROBE_COLUMNS = ["probe", "form", "dx_mm", "processes", "config", "k0", "k1", "crossings", "steps", "wall_s",
                 "ms_per_step", "dt_mean_s", "ke_internal_ratio", "pairing_start_dx", "pairing_end_dx", "growth",
                 "core_sxx_Pa", "core_sxx_rz22_Pa", "core_sxx_3d22_Pa", "vm_peak_Pa"]
REPLAY_COLUMNS = ["form", "dx_mm", "k", "kappa", "pairing_dx", "pairing_base_dx", "core_sxx_Pa", "core_sxx_base_Pa",
                  "core_sxx_rz22_base_Pa", "core_sxx_3d22_base_Pa", "vm_peak_base_Pa", "diff_vm_max", "diff_vm_rms",
                  "diff_p_max", "diff_core_vm_max", "diff_core_p_max", "hourglass_energy_ratio", "ke_internal_ratio",
                  "steps", "ms_per_step", "ms_per_step_base"]


def probe_dir(name):
    return REUSE.get(name, os.path.join(OUT, "probes", name))


def src_run(form, dx):
    return c13.run_dir(form, dx)


def launch_probe(name, form, dx, procs, k, Ns, cfg):
    d = probe_dir(name)
    if all(os.path.isfile(os.path.join(REPO, d, "k{}_N{}.npz".format(k, n))) for n in Ns):
        return
    cmd = [LAUNCHER] + (["-n", str(procs)] if procs > 1 else []) + [
        PROBE_SCRIPT, "interval", json.dumps(cfg), os.path.relpath(src_run(form, dx), REPO), str(k), d,
        ",".join(str(n) for n in Ns)]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=REPO, check=True, stdout=subprocess.DEVNULL)


def launch_timing(name, cfg):
    d = os.path.join(OUT, "timing", name)
    if os.path.isfile(os.path.join(REPO, d, "timing.json")):
        return
    cmd = [LAUNCHER, PROBE_SCRIPT, "timing", json.dumps(cfg), os.path.relpath(src_run("rz", 1.1), REPO), "4", d]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=REPO, check=True, stdout=subprocess.DEVNULL)


def replay_name(form, dx, kappa=KAPPA):
    b = naming.canonical_brackets({"monitor_every": 0, "hourglass": kappa})
    return naming.run_name(os.path.basename(c13.PREPARED), "replay", c13.FRAMES, form, dx, b, 1)


def replay_dir(form, dx, kappa=KAPPA):
    return naming.run_dir(os.path.join(REPO, OUT, "replays"), replay_name(form, dx, kappa))


def launch_replay(form, dx, procs, k1, max_load, force):
    d = replay_dir(form, dx)
    if os.path.isfile(os.path.join(d, "checkpoints", checkpoint.filename(k1))):
        return "done"
    if procs > 4 and not force:
        import spheral_m2_check4 as c4
        c4.wait_for_load(0.5*max_load)
    cmd = [LAUNCHER] + (["-n", str(procs)] if procs > 1 else []) + [
        "spheral_frag/runner/run.py", "replay", "--body", c13.body_dir(form, dx), "--prepared", c13.PREPARED,
        "--frames", "{}-{}".format(*c13.FRAMES), "--hourglass", repr(KAPPA), "--monitor-every", "0",
        "--stop-at", str(k1), "--out", os.path.join(OUT, "replays"), "--max-load", repr(max_load), "--resume"]
    if force:
        cmd.append("--force")
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    return "ok" if r.returncode == 0 else "FAILED: " + r.stderr.strip()[-300:]


# ----------------------------------------------------------------------------------------------- the checks 4-6
def checks_runs():
    """(check, form, dx, processes, run command args without the launcher, run dir) for the reruns of checks 4-6
    with the hourglass: the cheapest representative of each form (module docstring)."""
    import spheral_m2_check4 as c4
    import spheral_m2_check5 as c5
    import spheral_m2_check6 as c6
    out = os.path.join(OUT, "checks")
    prep = os.path.basename(c4.PREPARED)
    runs = []
    for form, dx, procs in (("rz", 2.2, 1), ("rz", 1.1, 4), ("3d", 3.0, 18)):
        b = naming.canonical_brackets({"crossings": 10, "damping": 0.05, "strength": "elastic", "sampler": "radial",
                                       "loads": "off", "hourglass": KAPPA})
        d = naming.run_dir(os.path.join(REPO, out), naming.run_name(prep, "replay", c4.FRAMES, form, dx, b, 1))
        cmd = ["spheral_frag/runner/run.py", "replay", "--body", c4.body_dir(form, dx), "--prepared", c4.PREPARED,
               "--frames", "{}-{}".format(*c4.FRAMES), "--crossings", "10", "--max-crossings", "10", "--damping",
               "0.05", "--strength", "elastic", "--loads", "off", "--ke-fraction", "1e300", "--hourglass",
               repr(KAPPA)]
        runs.append(("check4", form, dx, procs, cmd, d))
    # the RZ 2.2 mm bar is check 5's collapsed case: under Task 10b's rule its monitor runs to the cap at interval 0
    # (17.5 %, 195 %, 230 % at 10, 20 and 40 crossings), so it runs with the monitor off, at 10 crossings
    for form, dx, procs, mon in (("rz", 2.2, 1, 0), ("rz", 1.1, 4, 10)):
        b = naming.canonical_brackets({"strength": c5.STRENGTH, "sampler": "uniform", "loads": "off",
                                       "geometry": "bar", "hourglass": KAPPA, "monitor_every": mon})
        d = naming.run_dir(os.path.join(REPO, out), naming.run_name(prep, "replay", c5.FRAMES, form, dx, b, 1))
        cmd = ["spheral_frag/runner/run.py", "replay", "--body", c5.body_dir(form, dx), "--prepared", c5.PREPARED,
               "--frames", "{}-{}".format(*c5.FRAMES), "--strength", c5.STRENGTH, "--loads", "off",
               "--ke-fraction", "1e300", "--hourglass", repr(KAPPA), "--monitor-every", str(mon)]
        runs.append(("check5", form, dx, procs, cmd, d))
    for form, dx, procs in (("rz", 2.2, 1), ("rz", 1.1, 4), ("3d", 3.0, 18)):
        b = dict(c6.brackets("patch"), hourglass=KAPPA)
        d = naming.run_dir(os.path.join(REPO, out), naming.run_name(prep, "replay", c6.FRAMES, form, dx,
                                                                   naming.canonical_brackets(b), 1))
        cmd = c6.run_command(form, dx, "patch", out)[:-2] + ["--hourglass", repr(KAPPA)]
        runs.append(("check6", form, dx, procs, cmd, d))
    return runs


def launch_checks(max_load, force):
    out = os.path.join(OUT, "checks")
    for check, form, dx, procs, cmd, d in checks_runs():
        if os.path.isfile(os.path.join(d, "summary.json")):
            continue
        full = [LAUNCHER] + (["-n", str(procs)] if procs > 1 else []) + cmd
        if "--out" not in full:
            full += ["--out", out]
        full += ["--max-load", repr(max_load)] + (["--resume"] if os.path.isdir(d) else []) + \
            (["--force"] if force else [])
        if procs > 4 and not force:
            import spheral_m2_check4 as c4
            c4.wait_for_load(0.5*max_load)       # as the check drivers: 18 ranks starting raise the load
        print(check, form, dx, flush=True)
        t = time.perf_counter()
        r = subprocess.run(full, cwd=REPO, capture_output=True, text=True)
        os.makedirs(os.path.join(REPO, OUT, "logs"), exist_ok=True)
        with open(os.path.join(REPO, OUT, "logs", "{}_{}_{}.log".format(check, form, dx)), "w") as fh:
            fh.write(" ".join(full) + "\n" + r.stdout + "\n--- stderr\n" + r.stderr)
        print("  {} in {:.0f} s".format("ok" if r.returncode == 0 else "FAILED", time.perf_counter() - t),
              flush=True)


# --------------------------------------------------------------------------------------------------- analysis
def load_state(path, table):
    z = dict(np.load(path))
    st = run.stress_state(dict(z, scale=float(z["scale"])), table)
    return z, st


def probe_rows(table):
    rows = []
    for name, form, dx, procs, k, Ns, cfg in PROBES:
        ck0 = checkpoint.read(os.path.join(src_run(form, dx), "checkpoints", checkpoint.filename(k)))
        p0 = c13.pairing(form, ck0["x"], ck0["x_ref"], dx*1e-3)
        refs = {}
        for f, dd in (("rz", 2.2), ("3d", 2.2)):
            ck = checkpoint.read(os.path.join(src_run(f, dd), "checkpoints", checkpoint.filename(k + 1)))
            refs[f] = c13.core_sxx(run.stress_state(ck, table))
        for n in Ns:
            path = os.path.join(REPO, probe_dir(name), "k{}_N{}.npz".format(k, n))
            if not os.path.isfile(path):
                print("  missing", path)
                continue
            z, st = load_state(path, table)
            p1 = c13.pairing(form, z["x"], z["x_ref"], dx*1e-3)
            steps = int(z["steps"])
            rows.append(dict(probe=name, form=form, dx_mm=dx, processes=procs, config=json.dumps(cfg, sort_keys=True),
                             k0=k, k1=k + 1, crossings=n, steps=steps, wall_s=float(z["wall"]),
                             ms_per_step=1e3*float(z["wall"])/steps,
                             dt_mean_s=float(z["dt_mean"]) if "dt_mean" in z else math.nan,
                             ke_internal_ratio=float(z["ke_internal_ratio"]), pairing_start_dx=p0,
                             pairing_end_dx=p1, growth=p1/p0, core_sxx_Pa=c13.core_sxx(st),
                             core_sxx_rz22_Pa=refs["rz"], core_sxx_3d22_Pa=refs["3d"],
                             vm_peak_Pa=float(st["von_mises"].max())))
    for name, cfg in TIMINGS:
        path = os.path.join(REPO, OUT, "timing", name, "timing.json")
        if not os.path.isfile(path):
            print("  missing", path)
            continue
        with open(path) as fh:
            t = json.load(fh)
        rows.append(dict(probe="timing_" + name, form="rz", dx_mm=1.1, processes=t["processes"],
                         config=json.dumps(cfg, sort_keys=True), k0=4, k1=5, crossings=1, steps=t["steps"],
                         wall_s=t["wall_s"], ms_per_step=t["ms_per_step"], dt_mean_s=t["dt_mean_s"],
                         ke_internal_ratio=t["ke_internal_ratio"]))
    return rows


def replay_rows(table):
    rows = []
    for form, dx, procs, k1 in REPLAYS:
        d = replay_dir(form, dx)
        base = src_run(form, dx)
        iv = {int(r["k1"]): r for r in run.read_intervals(os.path.join(d, "intervals.csv"))}
        ivb = {int(r["k1"]): r for r in run.read_intervals(os.path.join(base, "intervals.csv"))}
        for k in range(1, k1 + 1):
            p = os.path.join(d, "checkpoints", checkpoint.filename(k))
            if not os.path.isfile(p):
                print("  missing", p)
                continue
            ck = checkpoint.read(p)
            ckb = checkpoint.read(os.path.join(base, "checkpoints", checkpoint.filename(k)))
            st, sb = run.stress_state(ck, table), run.stress_state(ckb, table)
            assert np.array_equal(st["pid"], sb["pid"])
            top = float(sb["von_mises"].max())
            dvm = run.von_mises(st["sigma"] - sb["sigma"])
            dp = np.abs(st["P"] - sb["P"])
            core = np.sqrt(np.einsum("ij,ij->i", st["x_ref"], st["x_ref"])) <= c13.CORE
            refs = {}
            for f, dd in (("rz", 2.2), ("3d", 2.2)):
                c = checkpoint.read(os.path.join(src_run(f, dd), "checkpoints", checkpoint.filename(k)))
                refs[f] = c13.core_sxx(run.stress_state(c, table))
            r, rb = iv[k], ivb[k]
            rows.append(dict(form=form, dx_mm=dx, k=k, kappa=KAPPA,
                             pairing_dx=c13.pairing(form, ck["x"], ck["x_ref"], dx*1e-3),
                             pairing_base_dx=c13.pairing(form, ckb["x"], ckb["x_ref"], dx*1e-3),
                             core_sxx_Pa=c13.core_sxx(st), core_sxx_base_Pa=c13.core_sxx(sb),
                             core_sxx_rz22_base_Pa=refs["rz"], core_sxx_3d22_base_Pa=refs["3d"], vm_peak_base_Pa=top,
                             diff_vm_max=float(dvm.max())/top, diff_vm_rms=math.sqrt(float(np.mean(dvm*dvm)))/top,
                             diff_p_max=float(dp.max())/top, diff_core_vm_max=float(dvm[core].max())/top,
                             diff_core_p_max=float(dp[core].max())/top,
                             hourglass_energy_ratio=float(r.get("hourglass_energy_ratio", "nan") or "nan"),
                             ke_internal_ratio=float(r["ke_internal_ratio"]), steps=int(r["steps"]),
                             ms_per_step=1e3*float(r["wall_per_step_s"]),
                             ms_per_step_base=1e3*float(rb["wall_per_step_s"])))
    return rows


CHECK_METRICS = {"check4": ("interior_err", "shell_err", "interior_rms_particle", "ke_internal_ratio",
                            "wall_per_step_ms"),
                 "check5": ("sxx_err_over_Y", "sxx_sd_Pa", "lat_mean_over_Y", "lat_max_over_Y", "ps_mean",
                            "ps_kinematic", "ke_internal_ratio", "wall_per_step_ms"),
                 "check6": ("int_vm_max", "int_vm_rms", "int_p_max", "core_vm_max", "core_vm_rms", "shell_vm_max",
                            "shell_vm_rms", "drag_rel", "ke_internal_ratio", "crossings_err_limit",
                            "wall_per_step_ms")}
CHECK_COLUMNS = ["check", "form", "dx_mm", "k", "metric", "committed", "hourglass", "difference"]


def check_rows(table):
    import spheral_m2_check4 as c4
    import spheral_m2_check5 as c5
    import spheral_m2_check6 as c6
    committed = {}
    for check, path in (("check4", c4.CSV), ("check5", c5.CSV), ("check6", c6.CSV)):
        with open(os.path.join(REPO, path), newline="") as fh:
            committed[check] = list(csv.DictReader(fh))
    rows = []
    E, nu, alpha = c4.constants(table)
    fes = c6.fe_summary()
    drag_ref = float(fes["tables"]["lee"]["continuum_drag_N"])
    fe_decel = -float(fes["meshes"]["1.6667"]["cases"]["lee"]["a_ms2"][0])
    fine, coarse = c6.reference(c6.FE_FINE), c6.reference(c6.FE_COARSE)
    for check, form, dx, procs, cmd, d in checks_runs():
        if not os.path.isfile(os.path.join(d, "summary.json")):
            print("  missing", d)
            continue
        iv = {int(r["k1"]): r for r in run.read_intervals(os.path.join(d, "intervals.csv"))}
        got = {}
        if check == "check4":
            for k in c4.COMPARE_K:
                m = c4.measure(d, k, table, E, nu, alpha, dx*1e-3)
                got[k] = dict(interior_err=m["interior"][0], shell_err=m["shell"][0],
                              interior_rms_particle=m["rms"], ke_internal_ratio=float(iv[k]["ke_internal_ratio"]),
                              wall_per_step_ms=1e3*float(iv[k]["wall_per_step_s"]))
            ref = {int(r["k"]): r for r in committed[check] if r["form"] == form and float(r["dx_mm"]) == dx
                   and r["crossings"] == "10" and float(r["damping"]) == 0.05 and r["scaled"] == "1"}
        elif check == "check5":
            E5, nu5, alpha5 = c5.constants(table)
            for k in range(c5.FRAMES[0] + 1, c5.FRAMES[1] + 1):
                m = c5.measure(d, k, table, E5, alpha5, dx*1e-3, form)
                ratio = c5.fe_reference(E5, nu5)["{}_{}".format(form, dx)]["ratio"]
                sxx_exp = -min(ratio*E5*alpha5*c5.DT_PER_FRAME*k, c5.Y)
                got[k] = dict(sxx_err_over_Y=(m["sxx_mean"] - sxx_exp)/c5.Y, sxx_sd_Pa=m["sxx_sd"],
                              lat_mean_over_Y=m["lat_mean"]/c5.Y, lat_max_over_Y=m["lat_max"]/c5.Y,
                              ps_mean=m["ps_mean"], ps_kinematic=m["ps_kin"],
                              ke_internal_ratio=float(iv[k]["ke_internal_ratio"]),
                              wall_per_step_ms=1e3*float(iv[k]["wall_per_step_s"]))
            ref = {int(r["k"]): r for r in committed[check] if r["form"] == form and float(r["dx_mm"]) == dx}
        else:
            with open(os.path.join(d, "summary.json")) as fh:
                summ = json.load(fh)
            for r in c6.measure_run(d, form, dx, "patch", table, fine, coarse, drag_ref, fe_decel,
                                    summ["processes"]):
                got[int(r["k"])] = r
            ref = {int(r["k"]): r for r in committed[check] if r["form"] == form and float(r["dx_mm"]) == dx
                   and r["area_shares"] == "patch"}
        for k in sorted(got):
            for metric in CHECK_METRICS[check]:
                a = float(ref[k][metric]) if k in ref and ref[k].get(metric, "") not in ("", None) else math.nan
                b = float(got[k].get(metric, math.nan))
                rows.append(dict(check=check, form=form, dx_mm=dx, k=k, metric=metric, committed=a, hourglass=b,
                                 difference=b - a))
    return rows


def write_rows(path, columns, rows):
    with open(os.path.join(REPO, path), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(columns)
        for row in rows:
            w.writerow([repr(float(row.get(c, math.nan))) if isinstance(row.get(c, math.nan), (float, np.floating))
                        else row.get(c, "") for c in columns])
    print("wrote {} ({} rows)".format(path, len(rows)))


def analyse(parts=("probe", "replay", "checks")):
    table = material.MaterialTable.load(os.path.join(REPO, c13.PREPARED))
    if "probe" in parts:
        rows = probe_rows(table)
        write_rows(PROBE_CSV, PROBE_COLUMNS, rows)
        for r in rows:
            if not r["probe"].startswith("timing"):
                print("{:18s} {:3s} {:3.1f} N {:3d}  pairing {:.4f} -> {:.4f} x{:5.2f}  core sxx {:6.2f} MPa (RZ 2.2 "
                      "{:6.2f}, 3D {:6.2f})  KE/E {:.2e}  {:.2f} ms/step".format(
                          r["probe"], r["form"], r["dx_mm"], r["crossings"], r["pairing_start_dx"],
                          r["pairing_end_dx"], r["growth"], r["core_sxx_Pa"]/1e6, r["core_sxx_rz22_Pa"]/1e6,
                          r["core_sxx_3d22_Pa"]/1e6, r["ke_internal_ratio"], r["ms_per_step"]))
            else:
                print("{:18s} steps {:6d}  dt {:.3e} s  {:.2f} ms/step  KE/E {:.2e}".format(
                    r["probe"], r["steps"], r["dt_mean_s"], r["ms_per_step"], r["ke_internal_ratio"]))
    if "replay" in parts:
        rows = replay_rows(table)
        write_rows(REPLAY_CSV, REPLAY_COLUMNS, rows)
        for r in rows:
            print("{:3s} {:3.1f} k {:2d}  pairing {:.4f} (base {:.4f})  core sxx {:6.2f} (base {:6.2f}; RZ 2.2 {:6.2f}, "
                  "3D {:6.2f})  vs base vM max {:.3f} rms {:.4f} P {:.3f} core {:.3f}/{:.3f}  E_hg/E {:.2e}  "
                  "{:.2f} ({:.2f}) ms/step".format(
                      r["form"], r["dx_mm"], r["k"], r["pairing_dx"], r["pairing_base_dx"], r["core_sxx_Pa"]/1e6,
                      r["core_sxx_base_Pa"]/1e6, r["core_sxx_rz22_base_Pa"]/1e6, r["core_sxx_3d22_base_Pa"]/1e6,
                      r["diff_vm_max"], r["diff_vm_rms"], r["diff_p_max"], r["diff_core_vm_max"],
                      r["diff_core_p_max"], r["hourglass_energy_ratio"], r["ms_per_step"], r["ms_per_step_base"]))
    if "checks" in parts:
        rows = check_rows(table)
        write_rows(CHECKS_CSV, CHECK_COLUMNS, rows)
        for r in rows:
            print("{} {:3s} {:3.1f} k {:2d} {:22s} {:+.4e} -> {:+.4e} ({:+.2e})".format(
                r["check"], r["form"], r["dx_mm"], r["k"], r["metric"], r["committed"], r["hourglass"],
                r["difference"]))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--only", choices=("probe", "replay", "checks", "analyse"))
    p.add_argument("--max-load", type=float, default=2.0)
    p.add_argument("--force", action="store_true")
    a = p.parse_args(argv)
    os.chdir(REPO)
    if a.only in (None, "probe"):
        for pr in PROBES:
            launch_probe(*pr)
        for name, cfg in TIMINGS:
            launch_timing(name, cfg)
    if a.only in (None, "replay"):
        for form, dx, procs, k1 in REPLAYS:
            print("replay", form, dx, launch_replay(form, dx, procs, k1, a.max_load, a.force), flush=True)
    if a.only in (None, "checks"):
        launch_checks(a.max_load, a.force)
    if a.only in (None, "analyse"):
        analyse()
    return 0


if __name__ == "__main__":
    sys.exit(main())
