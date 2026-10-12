"""Spheral M2 check 4: the heated sphere against Timoshenko and Goodier, with and without mass scaling, and the
crossings per interval and damping studies (M2 plan, Task 6; spec §12 check 4, §9.1, §14).

    "$PY" analysis/spheral_m2_check4.py [--only run|analyse] [--max-load 2.0] [--force]

**The case.** The 100 mm sphere of the MVP flight's frames 0-3 (the intact sphere; no particle is removed), its
temperatures replaced by the body stage's `--synthetic radial` profile T = 300 K + 50 K min(k, 1) (r/R)^2: ramped
in linearly over the first 0.5 s interval and held for two more. Elastic (`--strength elastic`), loads off
(`--loads off`: pure thermal stress), mass-scaled replay as spec §9.1, compared at the end of interval 1 (k = 1,
just after the ramp) and of interval 3 (k = 3, held for 1 s). Exact: `analytic.sphere_quadratic_stress` with the
material table's E, nu (constant below the solidus) and the free density's expansion coefficient referred to 300 K
(`analytic.expansion_strain_alpha`); peak sigma_0 = 2 alpha E dT / (5 (1 - nu)) = 50.07 MPa.

**The runs** (`DESIGN`). Per form and spacing (RZ 2.2 and 1.1 mm, 3D 3.0 and 2.2 mm): the crossings sweep 5, 10, 20
at the default damping 0.05, and the damping sweep 0, 0.2 at the default 10 crossings (the 0.05 run is shared): five
runs, not the 3 x 3 product, because the two devices act on different things -- the crossings set how far the
quasi-static ramp lags (an error during the ramp), the damping how fast the ringing it leaves decays (the held
state) -- and the plan's answers ask about each at the other's default. Without mass scaling (`--mass-scaling off`,
s = 1, each interval in N_c D / c_l seconds): crossings 5, 10, 20 at damping 0.05 for RZ 2.2 mm and 3D 3.0 mm, each
the twin of a scaled run. Every run uses `--ke-fraction 1e300` (no automatic repeat at more crossings: the study
measures each crossing number as given and records the KE ratios), so the runs live in their own directory
(`spheral_output/m2/check4/runs/`), not with replays at the default slowness test.

**Measures** per run and compared frame (`measure`): the particles' stresses from the checkpoint
(`run.stress_state`), resolved into the radial and mean tangential components about the sphere's centre at each
particle's reference position (`analytic.spherical_components`); the error e = sigma - sigma_exact(r_ref) per
particle; its volume-weighted mean in radial bins one spacing wide counted inward from R (`binned_profile`). The
**interior error** is the largest |binned mean e| of either component over the bins inside R - 2 dx, the **surface
shell error** the same over the two outer bins, both over sigma_0; `interior_rms_particle` the particle scatter
(RMS of e over interior particles, both components). The **scaled-unscaled difference** is the same binned measure
of sigma_scaled - sigma_unscaled per particle (the same pids), interior bins; `unscaled_diff_particle_max` the
largest per-particle difference of any component of the full tensor. The **damping difference** is the binned
interior measure of sigma(zeta) - sigma(0.05) at 10 crossings. KE ratios are intervals.csv's at the compared
frame; cost the run's steps, wall seconds and processes.

**The floor of agreement** (what the exact small-strain solution cannot resolve, a priori): the table's free
density gives the free strain alpha_eff (T - 300 K) exactly (checked by `constants` to 1e-9), but the body takes
it as a volume ratio, (1 + alpha dT)^3 - 1 against 3 alpha dT, relative 1.2e-3; and P = K (rho/rho_free - 1) is
-K e (1 - e + ...) in the mechanical volumetric strain e ~ sigma_0 / K = 7.1e-4. Together about 2e-3 of sigma_0, an
order below every measured interior error.

**Measured** (2026-10-10, 34 runs, the CSV). *Option A is sound*: every scaled run and its unscaled twin agree to
at most 3.1e-5 sigma_0 in the binned interior (3.5e-8 to 3.1e-5; per particle at most 3.1e-5), against interior
errors of 1.5 % or more: three orders below. Reruns of three runs that had been slowed by the host's sleep gave
bitwise-identical stresses. *Errors at the working values* (10-20 crossings, damping 0.05-0.2): interior 1.3-3.7 %
of sigma_0, the surface shell 5.5-9.8 % (the outermost bin's hoop stress, short of the exact -sigma_0 by SPH's
free-surface deficiency, and not falling with spacing: RZ 9.8 -> 8.5 %, 3D 6.2 -> 7.9 %). The settled interior
error (20 crossings, k 3) is a near-uniform isotropic offset of the core, +1.4 % in RZ and -1.6 to -2.4 % in 3D,
falling slowly with spacing (RZ 2.33 -> 1.69 %, 3D 2.64 -> 2.07 %); the same in scaled and unscaled runs, so a
property of the particle discretisation, not of mass scaling (a candidate for Task 9's SolidCRKSPH comparison).
*Crossings*: 5 never qualifies (3-11 %); the plan's rule gives 10 at every form, spacing and frame except RZ 1.1 mm
right after the ramp (k 1: 2.09 % against 1.2 x 1.48 %, the quasi-static lag of the ramp, 0.6 % of sigma_0), where
it gives 20; at the held state 10 everywhere. Internal KE/E at 10 crossings: at most 3.7e-3 at the end of the ramp,
1.5e-6 to 7.7e-6 held. *Damping*: 0.2 against 0.05 moves the held state by 0.03-0.49 % (below the interior
error); 0 leaves the body ringing (internal KE/E 2e-3 to 1.1e-2 after 1 s of hold, the held state off by 2-17 %).
*Cost* per 0.5 s interval at 10 crossings: RZ 2.2 mm 6 s (1 process, 3.7 ms/step), RZ 1.1 mm 19 s (4 processes),
3D 3.0 mm 62 s and 3D 2.2 mm 189 s (18 processes, 52 and 118 ms/step); the 34 runs 1.5 h of compute.

**Running it.** The host must not sleep: Docker's VM is suspended with it, which stretches wall times (not results);
run the driver under `caffeinate -i -s`. Before each run on more than 4 processes the driver waits for the VM's load
to fall to half of --max-load, since starting 18 ranks raises it by about 0.5 before run.py tests it.

Outputs: `data/spheral/m2_check4.csv` (one row per run and compared frame; committed), the plots and profiles in
`spheral_output/m2/check4/`. Resumable: bodies and runs whose outputs exist are skipped (run.py's summary.json).
Runs on more than 4 processes check the Linux VM's load as bench.py does (run.py refuses above --max-load).
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

from spheral_frag import analytic, body as body_stage, material, naming  # noqa: E402
from spheral_frag.runner import checkpoint, run  # noqa: E402

PREPARED = os.path.join("spheral_output", "prepare", "prep_model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_"
                        "none_fem-physics_melt-girin_material-AA7075_scheil_deeprunoff-off_k00000-01254")
FRAMES = (0, 3)
COMPARE_K = (1, 3)
OUT = os.path.join("spheral_output", "m2", "check4")
RUNS_OUT = os.path.join(OUT, "runs")
CSV = os.path.join("data", "spheral", "m2_check4.csv")
LAUNCHER = os.path.join("spheral_frag", "container", "spheral")
R = body_stage.RADIAL_R_M
DELTA_K = body_stage.RADIAL_DELTA_K
T_BASE = body_stage.RADIAL_T_BASE_K
# the hourglass stiffness (Task 10b): these runs were made before it, without it (_hourglass-0 in their names)
HOURGLASS = 0.0
SHELL_DX = 2.0                      # the surface shell: r > R - 2 dx (plan Task 6 Step 3)

FORMS = [("rz", 2.2, 1), ("rz", 1.1, 4), ("3d", 3.0, 18), ("3d", 2.2, 18)]     # form, dx mm, processes
SWEEP = [(5, 0.05), (10, 0.05), (20, 0.05), (10, 0.0), (10, 0.2)]               # (crossings, damping)
UNSCALED = {("rz", 2.2): [5, 10, 20], ("3d", 3.0): [5, 10, 20]}                # at damping 0.05

COLUMNS = ["form", "dx_mm", "crossings", "damping", "scaled", "k", "t_s", "n_particles", "processes", "scale",
           "time_scale", "peak_exact_Pa", "interior_err", "interior_err_r", "interior_err_t", "shell_err",
           "shell_err_r", "shell_err_t", "interior_rms_particle", "ke_ratio", "ke_internal_ratio",
           "max_disp_ref_dx", "steps_run", "wall_run_s", "wall_per_step_ms", "unscaled_diff_interior",
           "unscaled_diff_particle_max", "damping_diff_interior", "run"]


def design():
    """[(form, dx_mm, procs, crossings, damping, scaled)] in the order run."""
    out = []
    for form, dx, procs in FORMS:
        for n, z in SWEEP:
            out.append((form, dx, procs, n, z, True))
        for n in UNSCALED.get((form, dx), []):
            out.append((form, dx, procs, n, 0.05, False))
    return out


def body_dir(form, dx):
    prep = os.path.basename(PREPARED)
    return naming.body_dir(os.path.join("spheral_output", "body"), naming.body_name(prep, *FRAMES, form, dx,
                                                                                     synthetic="radial"))


def run_name(form, dx, n, z, scaled):
    b = naming.canonical_brackets({"crossings": n, "damping": z, "strength": "elastic", "sampler": "radial",
                                   "loads": "off", "mass_scaling": "on" if scaled else "off",
                                   "hourglass": HOURGLASS})
    return naming.run_name(os.path.basename(PREPARED), "replay", FRAMES, form, dx, b, 1)


def run_dir(form, dx, n, z, scaled):
    return naming.run_dir(os.path.join(REPO, RUNS_OUT), run_name(form, dx, n, z, scaled))


# ------------------------------------------------------------------------------------------------------- running
def build_bodies():
    for form, dx, _ in FORMS:
        if os.path.isfile(os.path.join(REPO, body_dir(form, dx), "body.json")):
            continue
        cmd = [sys.executable, "-m", "spheral_frag", "body", "--prepared", PREPARED, "--frames", "{}-{}".format(*FRAMES),
               "--form", form, "--dx-mm", str(dx), "--synthetic", "radial", "--quiet"]
        subprocess.run(cmd, cwd=REPO, check=True)


def vm_load():
    """The Linux VM's 1-minute load average (what run.py tests), read through the launcher."""
    r = subprocess.run([LAUNCHER, "-c", "import os; print(os.getloadavg()[0])"], cwd=REPO, capture_output=True,
                       text=True)
    try:
        return float(r.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return math.inf


def wait_for_load(max_load, timeout_s=1800.0, poll_s=20.0):
    """Wait until the VM's load is at most max_load (the previous run's own load decays over a minute or two)."""
    t0 = time.perf_counter()
    while True:
        load = vm_load()
        if load <= max_load or time.perf_counter() - t0 > timeout_s:
            return load
        time.sleep(poll_s)


def launch(form, dx, procs, n, z, scaled, max_load, force):
    d = run_dir(form, dx, n, z, scaled)
    if os.path.isfile(os.path.join(d, "summary.json")):
        return "done"
    if procs > 4 and not force:
        # half the threshold: starting 18 ranks that import Spheral raises the load by ~0.5 before run.py tests it
        wait_for_load(0.5*max_load)
    cmd = [LAUNCHER] + (["-n", str(procs)] if procs > 1 else []) + [
        "spheral_frag/runner/run.py", "replay", "--body", body_dir(form, dx), "--prepared", PREPARED,
        "--frames", "{}-{}".format(*FRAMES), "--crossings", str(n), "--max-crossings", str(n), "--damping", repr(z),
        "--strength", "elastic", "--loads", "off", "--mass-scaling", "on" if scaled else "off",
        "--ke-fraction", "1e300", "--out", RUNS_OUT, "--max-load", repr(max_load), "--hourglass", repr(HOURGLASS)]
    if os.path.isdir(d):
        cmd.append("--resume")
    if force:
        cmd.append("--force")
    t = time.perf_counter()
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    os.makedirs(os.path.join(REPO, OUT, "logs"), exist_ok=True)
    with open(os.path.join(REPO, OUT, "logs", naming.run_name_from_dir(d) + ".log"), "w") as fh:
        fh.write(" ".join(cmd) + "\n" + r.stdout + "\n--- stderr\n" + r.stderr)
    if r.returncode != 0:
        print("  FAILED ({}): {}".format(r.returncode, r.stderr.strip().splitlines()[-1:] if r.stderr else ""),
              flush=True)
        return "failed"
    return "ran in {:.0f} s".format(time.perf_counter() - t)


# ----------------------------------------------------------------------------------------------------- measuring
def constants(table):
    """(E, nu, alpha_eff): the table's E and nu at 300 K, checked constant to 350 K, and the free density's linear
    expansion coefficient referred to 300 K (module docstring)."""
    T = np.linspace(300.0, 350.0, 51)
    E, nu = table.young(T), table.poisson(T)
    if np.ptp(E) != 0.0 or np.ptp(nu) != 0.0:
        raise RuntimeError("E or nu is not constant between 300 and 350 K in this table")
    alpha = float(table.header["mechanical"]["entries"]["alpha"]["value"])
    a_eff = analytic.expansion_strain_alpha(alpha, T_BASE, material.T_REF_EXPANSION)
    # the table's rho_free is rho_s / (1 + alpha (T - 293))^3 on a 1 K grid: its free strain must be a_eff (T - 300)
    lin = (table.rho_free(T_BASE)/table.rho_free(T))**(1.0/3.0) - 1.0
    dev = float(np.max(np.abs(lin - a_eff*(T - T_BASE))))
    if dev > 1e-9*a_eff*50.0:
        raise RuntimeError("the table's free density is not the expected expansion law ({:.2e})".format(dev))
    return float(E[0]), float(nu[0]), a_eff


def edges_for(dx):
    n = int(math.floor(R/dx))
    e = R - dx*np.arange(n + 1)[::-1]
    e[0] = 0.0                                             # the innermost bin absorbs the remainder
    e[-1] = R*(1.0 + 1e-9)
    return e


def components(ck, table):
    st = run.stress_state(ck, table)
    r, srr, stt = analytic.spherical_components(st["x_ref"], st["sigma"])
    return st, r, srr, stt


def binned_errors(r, er, et, V, dx):
    """(interior max over bins and components, (r, t) parts, shell likewise, bin centres, binned (er, et))."""
    edges = edges_for(dx)
    c, m, cnt = analytic.binned_profile(r, np.stack([er, et], axis=1), V, edges)
    hi = edges[1:]
    interior = (hi <= R - SHELL_DX*dx + 1e-12) & (cnt > 0)
    shell = (edges[:-1] >= R - SHELL_DX*dx - 1e-12) & (cnt > 0)
    ir = float(np.nanmax(np.abs(m[interior, 0])))
    it = float(np.nanmax(np.abs(m[interior, 1])))
    sr = float(np.nanmax(np.abs(m[shell, 0])))
    stt = float(np.nanmax(np.abs(m[shell, 1])))
    return (max(ir, it), ir, it), (max(sr, stt), sr, stt), c, m


def measure(d, k, table, E, nu, alpha, dx):
    ck = checkpoint.read(os.path.join(d, "checkpoints", checkpoint.filename(k)))
    st, r, srr, stt = components(ck, table)
    ex_r, ex_t = analytic.sphere_quadratic_stress(r, R, E, nu, alpha, DELTA_K)
    s0 = analytic.quadratic_peak(E, nu, alpha, DELTA_K)
    er, et = srr - ex_r, stt - ex_t
    (ie, ir, it), (se, sr_, st_), c, m = binned_errors(r, er, et, st["V"], dx)
    inner = r <= R - SHELL_DX*dx
    w = st["V"][inner]
    rms = math.sqrt(float(np.sum(w*0.5*(er[inner]**2 + et[inner]**2))/np.sum(w)))/s0
    return dict(ck=ck, st=st, r=r, srr=srr, stt=stt, s0=s0, interior=(ie/s0, ir/s0, it/s0),
                shell=(se/s0, sr_/s0, st_/s0), rms=rms, centres=c, binned=m)


def pair_diff(a, b, dx):
    """Binned interior difference and the largest per-particle tensor difference of two measured runs, over sigma_0."""
    if not np.array_equal(a["st"]["pid"], b["st"]["pid"]):
        raise RuntimeError("the two runs hold different particles")
    (ie, _, _), _, _, _ = binned_errors(a["r"], a["srr"] - b["srr"], a["stt"] - b["stt"], a["st"]["V"], dx)
    pmax = float(np.abs(a["st"]["sigma"] - b["st"]["sigma"]).max())
    return ie/a["s0"], pmax/a["s0"]


def analyse():
    table = material.MaterialTable.load(os.path.join(REPO, PREPARED))
    E, nu, alpha = constants(table)
    s0 = analytic.quadratic_peak(E, nu, alpha, DELTA_K)
    print("exact: E {:.4g} Pa, nu {}, alpha_eff {:.6e} 1/K, sigma_0 {:.6g} MPa".format(E, nu, alpha, s0/1e6))
    rows, measured = [], {}
    for form, dx_mm, procs, n, z, scaled in design():
        d = run_dir(form, dx_mm, n, z, scaled)
        if not os.path.isfile(os.path.join(d, "summary.json")):
            print("  missing: {}".format(naming.run_name_from_dir(d)))
            continue
        with open(os.path.join(d, "summary.json")) as fh:
            summ = json.load(fh)
        iv = {int(r["k1"]): r for r in run.read_intervals(os.path.join(d, "intervals.csv"))}
        for k in COMPARE_K:
            mm = measure(d, k, table, E, nu, alpha, dx_mm*1e-3)
            measured[(form, dx_mm, n, z, scaled, k)] = mm
            r = iv[k]
            rows.append(dict(form=form, dx_mm=dx_mm, crossings=n, damping=z, scaled=int(scaled), k=k,
                             t_s=float(r["t1_s"]), n_particles=len(mm["r"]), processes=summ["processes"],
                             scale=float(r["scale"]), time_scale=float(r["time_scale"]), peak_exact_Pa=s0,
                             interior_err=mm["interior"][0], interior_err_r=mm["interior"][1],
                             interior_err_t=mm["interior"][2], shell_err=mm["shell"][0], shell_err_r=mm["shell"][1],
                             shell_err_t=mm["shell"][2], interior_rms_particle=mm["rms"],
                             ke_ratio=float(r["ke_ratio"]), ke_internal_ratio=float(r["ke_internal_ratio"]),
                             max_disp_ref_dx=float(r["max_disp_ref_dx"]), steps_run=summ["steps"],
                             wall_run_s=summ["wall_s_intervals"],
                             wall_per_step_ms=1e3*summ["wall_s_intervals"]/summ["steps"],
                             unscaled_diff_interior=math.nan, unscaled_diff_particle_max=math.nan,
                             damping_diff_interior=math.nan, run=naming.run_name_from_dir(d)))
    for row in rows:
        key = (row["form"], row["dx_mm"], row["crossings"], row["damping"])
        twin = measured.get(key + (not bool(row["scaled"]), row["k"]))
        if twin is not None:
            a = measured[key + (bool(row["scaled"]), row["k"])]
            row["unscaled_diff_interior"], row["unscaled_diff_particle_max"] = pair_diff(a, twin, row["dx_mm"]*1e-3)
        if row["scaled"] and row["crossings"] == 10 and row["damping"] != 0.05:
            ref = measured.get((row["form"], row["dx_mm"], 10, 0.05, True, row["k"]))
            if ref is not None:
                a = measured[key + (True, row["k"])]
                row["damping_diff_interior"] = pair_diff(a, ref, row["dx_mm"]*1e-3)[0]
    path = os.path.join(REPO, CSV)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COLUMNS)
        for row in rows:
            w.writerow([repr(float(row[c])) if isinstance(row[c], (float, np.floating)) else row[c] for c in COLUMNS])
    print("wrote {} ({} rows)".format(CSV, len(rows)))
    report(rows)
    plots(measured, s0, E, nu, alpha)
    return rows


def report(rows):
    print("\n{:3s} {:>4s} {:>3s} {:>5s} {:>3s} {:>2s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>7s}".format(
        "frm", "dx", "N", "zeta", "sc", "k", "interior", "shell", "rms_part", "KEint/E", "unsc.diff", "damp.diff",
        "steps", "wall_s"))
    for r in rows:
        print("{:3s} {:4.1f} {:3d} {:5.2f} {:3d} {:2d} {:9.3e} {:9.3e} {:9.3e} {:9.3e} {:9.3e} {:9.3e} {:9d} {:7.0f}"
              .format(r["form"], r["dx_mm"], r["crossings"], r["damping"], r["scaled"], r["k"], r["interior_err"],
                      r["shell_err"], r["interior_rms_particle"], r["ke_internal_ratio"],
                      r["unscaled_diff_interior"], r["damping_diff_interior"], r["steps_run"], r["wall_run_s"]))


# -------------------------------------------------------------------------------------------------------- plots
INK, MUTED = "#1f2328", "#6e7781"
RAMP = {5: "#9ecae1", 10: "#4292c6", 20: "#08519c"}                 # crossings: one hue, light -> dark
DAMP = {0.0: "#c6dbef", 0.05: "#4292c6", 0.2: "#08306b"}


def plots(measured, s0, E, nu, alpha):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = os.path.join(REPO, OUT)
    os.makedirs(out, exist_ok=True)
    rr = np.linspace(0.0, R, 200)
    exr, ext = analytic.sphere_quadratic_stress(rr, R, E, nu, alpha, DELTA_K)
    for form, dx_mm, _ in FORMS:
        for k in COMPARE_K:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
            for ax, sweep, colours, label in ((axes[0], [(n, 0.05) for n in (5, 10, 20)], RAMP, "crossings"),
                                              (axes[1], [(10, z) for z in (0.0, 0.05, 0.2)], DAMP, "damping")):
                ax.plot(rr*1e3, exr/1e6, color=INK, lw=1.5, label="exact sigma_r")
                ax.plot(rr*1e3, ext/1e6, color=INK, lw=1.5, ls="--", label="exact sigma_theta")
                for n, z in sweep:
                    mm = measured.get((form, dx_mm, n, z, True, k))
                    if mm is None:
                        continue
                    key = n if label == "crossings" else z
                    c = mm["centres"]
                    ex_c = analytic.sphere_quadratic_stress(c, R, E, nu, alpha, DELTA_K)
                    ax.plot(c*1e3, (ex_c[0] + mm["binned"][:, 0])/1e6, color=colours[key], lw=2, marker="o", ms=3,
                            label="{} {} sigma_r".format(label, key))
                    ax.plot(c*1e3, (ex_c[1] + mm["binned"][:, 1])/1e6, color=colours[key], lw=2, ls="--",
                            marker="s", ms=3, label="{} {} sigma_theta".format(label, key))
                ax.axvline((R - SHELL_DX*dx_mm*1e-3)*1e3, color=MUTED, lw=0.8, ls=":")
                ax.set_xlabel("reference radius (mm)")
                ax.grid(color="#d0d7de", lw=0.5)
                ax.set_title("{} sweep".format(label), fontsize=10)
                ax.legend(fontsize=7, frameon=False)
            axes[0].set_ylabel("stress (MPa, binned)")
            fig.suptitle("Check 4: {} {} mm, k = {} (t = {} s)".format(form.upper(), dx_mm, k, 0.5*k), fontsize=11)
            fig.tight_layout()
            fig.savefig(os.path.join(out, "profiles_{}_dx{}_k{}.png".format(form, dx_mm, k)), dpi=130)
            plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    marks = {("rz", 2.2): "o", ("rz", 1.1): "s", ("3d", 3.0): "^", ("3d", 2.2): "D"}
    shades = {("rz", 2.2): "#9ecae1", ("rz", 1.1): "#2171b5", ("3d", 3.0): "#a1d99b", ("3d", 2.2): "#238b45"}
    for ax, k in zip(axes, COMPARE_K):
        for form, dx_mm, _ in FORMS:
            xs, ys = [], []
            for n in (5, 10, 20):
                mm = measured.get((form, dx_mm, n, 0.05, True, k))
                if mm is not None:
                    xs.append(n)
                    ys.append(mm["interior"][0])
            ax.plot(xs, ys, marker=marks[(form, dx_mm)], color=shades[(form, dx_mm)], lw=2, ms=8,
                    label="{} {} mm".format(form.upper(), dx_mm))
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xticks([5, 10, 20], ["5", "10", "20"])
        ax.set_xlabel("sound crossings per interval")
        ax.set_title("k = {}".format(k), fontsize=10)
        ax.grid(color="#d0d7de", lw=0.5)
        ax.legend(fontsize=8, frameon=False)
    axes[0].set_ylabel("interior error / sigma_0 (binned)")
    fig.suptitle("Check 4: interior error against crossings (damping 0.05)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "interior_error_vs_crossings.png"), dpi=130)
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--only", choices=("run", "analyse"))
    p.add_argument("--max-load", type=float, default=2.0)
    p.add_argument("--force", action="store_true", help="start runs even above --max-load")
    a = p.parse_args(argv)
    os.chdir(REPO)
    if a.only != "analyse":
        build_bodies()
        for form, dx, procs, n, z, scaled in design():
            print("{} {} mm N_c {} zeta {} {}: ".format(form, dx, n, z, "scaled" if scaled else "UNSCALED"), end="",
                  flush=True)
            print(launch(form, dx, procs, n, z, scaled, a.max_load, a.force), flush=True)
    if a.only != "run":
        analyse()
    return 0


if __name__ == "__main__":
    sys.exit(main())
