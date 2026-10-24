"""Spheral M2 check 6: the loads on an elastic sphere, the deceleration and the finite-element reference (M2 plan,
Task 8 Steps 2-3; spec §7, §12 check 6).

    "$PY" analysis/spheral_m2_check6.py [--only run|analyse] [--max-load 2.0] [--force]

**The case.** The 100 mm sphere of the MVP flight's frame 0 (the intact sphere; frames 0-FRAMES[1] are before the
first film, nothing is removed), held at a uniform 300 K at every frame (the body stage's `--synthetic uniform
--uniform-delta-k 0`, sampler "uniformdt0": no thermal stress), `--strength elastic`, the surface loads and the
deceleration on with **frame 182's load table held** (`--loads-frame 182`: peak dynamic pressure, plan fact 8; the
table ramps in linearly from zero over the first 0.5 s interval and is held after, run.py "One frame's loads
held"), the mass-scaled replay at the default crossings (10) and damping (0.05), the crossings monitor at its
default (it estimates interval 0 -> 1). RZ at 2.2 and 1.1 mm, 3D at 3.0 and 2.2 mm. `--ke-fraction 1e300` as in
checks 4 and 5: the study measures the default crossings as given and records the KE ratios.

**Why the ramp** (measured on RZ 2.2 mm, 2026-10-10, `spheral_output/m2_probe/check6_sudden/`): applied at full
strength at the first step the loads set the body ringing; at 10 crossings the damping leaves it ringing at the end
of the interval, the monitor's 20-crossing rerun ends at another phase of it, the estimate is 32 % and the run
switches to 20 crossings. With the ramp the monitor estimates 2.3 % and the run stays at 10; the held state is the
same (linear elasticity: the static solution does not depend on the path).

**How long** (`FRAMES`): until the internal KE is below 1 % of the elastic energy and the stress field no longer
changes. Measured per frame: `settle_vm` = max over particles of vM(sigma_k - sigma_(k-1)) / VM0 and `settle_p` the
same for the pressure. Over six intervals (frames 0-6: the ramp and five held) both fall by a factor 0.1-0.4 per
interval (the damping's e^(-zeta pi N_c) = 0.21 at 10 crossings), and the comparison is made at the last frame
(`COMPARE_K`).

**The reference.** `analysis/spheral_m2_elastic_fe.py` (Step 1): P2 tetrahedra on the ball of R = 50 mm, the same
tractions -p(theta) n + tau(theta) t from frame 182's table (with the lee extension, which sets no bin at frame 182),
balanced by the inertial body force -rho (a + alpha x (x - x_c)); the fine (1.667 mm) solution, sampled at each
particle's reference position (`fe_stress`: the containing element by `geometry.TetLocator`, the exact gradient of
the P2 field there; a point outside the faceted ball -- a lattice point within a few micrometres of the true sphere
-- takes the nearest element's polynomial, extrapolated by that distance), and the coarse (2.5 mm) solution at the
same points for the reference's own uncertainty (`fe_conv_*`). VM0 = 3,510.5 Pa, the reference's largest sampled
von Mises stress (Step 1). RZ: a ring particle's reference (x, r) is the point (x, r, 0) of the x-y plane, its
stress rows/columns (axial, radial, hoop) the Cartesian (xx, yy, zz) there.

**Measures** per run and frame k (`measure`):

- *stresses*: the particle's full stress sigma = S - P I (`run.stress_state`) against the reference's; the von Mises
  stress of the difference (deviatoric part) and the pressure difference |tr(sigma_p - sigma_FE)| / 3 (the part von
  Mises does not see), each over VM0, their maximum and volume-weighted RMS over the interior (|x_ref| <= 45 mm) and
  the surface shell (|x_ref| > 45 mm) separately; the signed volume-weighted mean of tr(sigma_p - sigma_FE) / 3 / VM0
  over the interior (`int_p_mean`; negative: the particles more compressed); the same von Mises max, RMS and mean
  pressure offset over the core |x_ref| <= 25 mm (`core_*`: the interior at 45 mm is only 1.7-4.5 spacings deep, so
  the core separates the surface layer's reach from the body's); the particles' largest von Mises over VM0;
- *settling*: `settle_vm`, `settle_p` (above);
- *forces and motion* (intervals.csv): the drag of the particle loads sum F (`drag_N`, at the held state), the
  deceleration sum F / M (`decel_ms2`, the run's physical -s a_frame) against the table's integral over the exact
  sphere divided by the body's mass (`decel_ref_ms2`; 30.717 N / M), the finite element's a_x (`fe_decel_ms2`,
  20.869 m/s^2) and the history's (`flight_decel_ms2`, 86.69 m/s^2 at k 182 -- a diagnostic only: the history's
  body at k 182 has lost 61 % of its mass and its drag is the trajectory's drag table, not the frame's loads, so
  it is not like for like); the net momentum left after the interval over its drag impulse
  (`momentum_over_impulse`), the RZ axial drift the boundary removes (`v_com_removed_ms`, `recentre_shift_dx`);
- *surface detection* (`surface`, recomputed in drama_env with the package's own functions from the checkpoint the
  held interval started from, and checked against the run's counts and drag): n surface, n exposed, the area
  shares' sum over 4 pi R^2 (`area_cover`: calibrated to frame 0's faceted surface area), the shares' closure
  |sum A n| over 4 pi R^2, the windward projected area sum_{n_x > 0} A n_x over pi R^2, the normals' angle from
  radial (max and RMS), the drag of the same loads with exact radial normals at the same particles
  (`drag_radial_normals_N`: what the normals' error costs), the loads' moment about the mass centre over drag R
  (`moment_rel`, 3D; the Deceleration package removes the net force, not a moment, so a moment turns the body
  slowly against the damping);
- *section equilibrium* (`section`): the particles' axial force through the mid-plane x = 0 (sigma_xx over the
  lattice plane x_ref = 0 times each particle's section) against the load the nose half transmits there (its
  surface forces less its inertia); `section_ratio` is their ratio, 1 in the continuum;
- the crossings monitor's estimates on interval 0 -> 1; KE ratios; cost.

**The area shares** (`SHARES`). Every run is made twice: with Task 4's shares ("volume", V^(2/3) with one
calibration factor; the default then) and with the frame's patch areas ("patch", the default since decision 9:
each finite-element surface patch's area to the nearest surface particle by reference position; packages/loads.py,
"Area shares"). The first RZ runs showed why:
on a lattice the volume shares give every surface particle the same area, but the staircase carries fewer surface
particles per unit area at 45 degrees than near the axis and the equator, so against the true (spherical Voronoi)
area each particle stands for, the shares are 10-11 % too large near 0, 90 and 180 degrees and 11-13 % too small
near 45 and 135 degrees, at 2.2 and 1.1 mm alike (a property of the lattice, not of the resolution).

**Measured** (2026-10-10/11, 8 runs, the CSV; thresholds in tests/test_spheral_frag_runner_check6.py). Stress
errors over VM0 = 3,510.5 Pa, at k 6 unless stated.

- *Settling.* Internal KE/E 1.8-6.5e-4 at the end of the ramp, at most 1.5e-5 held (k >= 2); the stress change
  between frames falls from 2-8 % (k 2) to 0.2-2e-4 of VM0 (k 6), so six frames settle it (by k 4 to 0.2-1.8e-3).
- *At rest.* 3D: net momentum after an interval at most 1.3e-16 of its drag impulse; RZ (SolidSPHRZ does not
  conserve axial momentum) at most 1.9e-3, the axial drift the boundary removes at most 3.5e-7 spacings per
  interval. With the patch shares the 3D loads keep a small moment (4-6e-5 of drag x R; the volume shares' is
  1e-16), which the Deceleration (a force, not a moment) leaves: the body turns slowly against the damping, 0.5-1.4e-8
  spacings per interval (3.0, 2.2 mm), internal KE/E staying at 0.9-4.4e-7.
- *Deceleration and drag.* sum F / M is the deceleration to round-off; against the table over the exact sphere
  (30.717 N / 1.47202 kg = 20.867 m/s^2; the finite element's 20.869): volume shares RZ -2.9 % and -3.0 % (2.2,
  1.1 mm), 3D -0.61 % and -0.17 % (3.0, 2.2 mm); patch shares RZ +0.35 % and +0.20 %, 3D -0.54 % and +0.08 %.
  The history's 86.69 m/s^2 at k 182 is 4.1-4.3 times larger: its body has lost 61 % of its mass (0.572 kg) and its
  drag is the trajectory's drag table (49.6 N, against frame 182's own table integral 32.78 N over its receded
  body) -- not like for like, a diagnostic.
- *Surface detection.* RZ 63 and 127 surface particles, 3D 3,350 and 6,194, every one exposed (convex body); the
  shares sum to frame 0's faceted area, 0.99967 of 4 pi R^2; normals within 4.8 and 5.7 degrees of radial in RZ
  (RMS 2.1-2.2), 1.2 and 2.4 in 3D (RMS 0.7-0.9); with exact radial normals at the same particles the drag changes
  by at most 0.6 %: the inclination is not where the error is. The patch shares leave 583 and 963 of the 3D surface
  particles (the staircase's inner corners) without area.
- *Stresses against the reference.* Interior (|x_ref| <= 45 mm) max and RMS: volume shares RZ 0.62 / 0.10 and
  0.35 / 0.087, 3D 1.19 / 0.13 and 0.58 / 0.068; patch shares RZ 0.28 / 0.051 and 0.19 / 0.041, 3D 0.46 / 0.056
  and 0.65 / 0.039 (the 45 mm interior starts 1.7-4.5 spacings in, inside the surface layer's reach). Core
  (|x_ref| <= 25 mm) max: volume RZ 0.21 and 0.21, 3D 0.13 and 0.066; patch RZ 0.058 and 0.065, 3D 0.083 and
  0.039. Shell (|x_ref| > 45 mm) max 1.4-1.7 (volume) and 1.8-3.0 (patch), RMS 0.18-0.31: the loads act as a body
  force on the outermost particle layer, whose own stress (and the next layer's) SPH does not resolve -- the
  free-surface limitation of check 4, here at the size of the peak stress itself, and not falling with spacing.
  The reference's own mesh difference at the particles is 0.3-0.8 % (interior) and 2-4 % (shell): it resolves far
  finer than the comparison.
- *RZ's section.* Through the mid-plane x = 0 the RZ particles carry 5.0-7.7 % more axial force than the nose half's
  load (net of its inertia) requires, at both spacings and with both shares; 3D carries 1.6-2.6 % less (its outer
  layer's deficit). This is the RZ core's near-uniform offset (sigma_xx 145-270 Pa more compressive, mean pressure
  offset 1.5-4.9 % of VM0), absent in 3D (core mean offset at most 0.34 %): a property of SolidSPHRZ under axial
  load, for Task 10's RZ-against-3D comparison.
- *Monitor* (interval 0 -> 1): patch shares 1.0-1.5 %, volume RZ 1.2-2.3 %; volume 3D 9.0 % and 6.5 % (von Mises),
  above the 5 % limit, so those two runs switched to 20 crossings.
- *Cost* per 0.5 s interval: RZ 2.2 mm 6.2 s (1 process, 3.8 ms/step), RZ 1.1 mm 19-20 s (4 processes), 3D 3.0 mm
  63 s at 10 crossings and 131 s at 20, 3D 2.2 mm 189 s at 10 and 387 s at 20 (18 processes, 54-56 and 118-120
  ms/step); the eight runs 1.6 h.

**The answer** (plan Task 8: do inclination tables work on particles?). The inclination does: the colour-gradient
normals are within 1-6 degrees, and with correct areas the particles' drag and deceleration match the table's
integral to 0.1-0.5 %. The area shares are what decide it: Task 4's volume shares misplace up to 13 % of the
load locally on a lattice (RZ drag -3 %, every core error 1.5-3.6 times larger than with the patch shares); the patch
shares give the load and the body's stresses to 4-8 % of the peak in the core and 4-6 % RMS over the interior. The
outer one or two spacings are not resolved by either (errors up to 1.4-3 times the peak there), and RZ adds its
5-8 % section excess.

**Running it.** Under `caffeinate -i -s` (the host's sleep suspends Docker's VM). Before each run on more than 4
processes the driver waits for the VM's load to fall to half of --max-load (as checks 4 and 5).

Outputs: `data/m2_check6.csv` (one row per run and frame; committed), the runs in
`spheral_output/m2/check6/runs/`, plots in `spheral_output/m2/check6/`. Resumable: bodies and runs whose outputs exist
are skipped.
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

from spheral_frag import material, naming  # noqa: E402
from spheral_frag.runner import checkpoint, run  # noqa: E402

PREPARED = os.path.join("spheral_output", "prepare", "prep_model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_"
                        "none_fem-physics_melt-girin_material-AA7075_scheil_deeprunoff-off_k00000-01254")
LOADS_FRAME = 182
FRAMES = (0, 6)
COMPARE_K = FRAMES[1]
SAMPLER = "uniformdt0"
# the hourglass stiffness (Task 10b): these runs were made before it, without it (_hourglass-0 in their names)
HOURGLASS = 0.0
STRENGTH = "elastic"
OUT = os.path.join("spheral_output", "m2", "check6")
RUNS_OUT = os.path.join(OUT, "runs")
FE_DIR = os.path.join(OUT, "fe")
FE_FINE = os.path.join(FE_DIR, "solution_lee_h1.6667mm.npz")
FE_COARSE = os.path.join(FE_DIR, "solution_lee_h2.5mm.npz")
CSV = os.path.join("data", "m2_check6.csv")
LAUNCHER = os.path.join("spheral_frag", "container", "spheral")
R = 0.05                                    # m
INTERIOR = 0.045                            # m: |x_ref| <= this is the interior
CORE = 0.025                                # m: |x_ref| <= this is the core (half the radius; a depth measure)
VM0 = 3510.5                                # Pa: the reference's largest sampled von Mises (Step 1)

FORMS = [("rz", 2.2, 1), ("rz", 1.1, 4), ("3d", 3.0, 18), ("3d", 2.2, 18)]     # form, dx mm, processes
SHARES = ("volume", "patch")                # the area shares: Task 4's, and the frame's patches (decision 9's default)

COLUMNS = ["form", "dx_mm", "area_shares", "k", "n_particles", "n_interior", "n_shell", "processes", "scale",
           "int_vm_max", "int_vm_rms", "int_p_max", "int_p_rms", "int_p_mean", "core_vm_max", "core_vm_rms",
           "core_p_mean", "shell_vm_max", "shell_vm_rms", "shell_p_max",
           "shell_p_rms", "particle_vm_max", "fe_conv_int_vm_max", "fe_conv_shell_vm_max", "settle_vm", "settle_p",
           "ke_ratio", "ke_internal_ratio", "drag_N", "drag_ref_N", "drag_rel", "frame_table_drag_N", "decel_ms2",
           "decel_ref_ms2", "decel_rel", "fe_decel_ms2", "flight_decel_ms2", "mass_kg", "momentum_axial_Ns",
           "momentum_over_impulse", "v_com_removed_ms", "recentre_shift_dx", "max_disp_ref_dx", "n_surface",
           "n_exposed", "area_cover", "closure", "projected_windward", "normal_err_max_deg", "normal_err_rms_deg",
           "drag_surface_N", "drag_radial_normals_N", "n_zero_area", "moment_rel", "section_force_N", "section_applied_N",
           "section_ratio", "monitored", "crossings_err_est", "crossings_err_pressure",
           "crossings_err_limit", "run_crossings", "steps", "wall_s", "wall_per_step_ms", "run"]


def body_dir(form, dx):
    prep = os.path.basename(PREPARED)
    return naming.body_dir(os.path.join("spheral_output", "body"),
                           naming.body_name(prep, *FRAMES, form, dx, synthetic=SAMPLER))


def brackets(shares="volume"):
    return naming.canonical_brackets({"strength": STRENGTH, "sampler": SAMPLER, "loads_frame": LOADS_FRAME,
                                      "area_shares": shares, "hourglass": HOURGLASS})


def run_name(form, dx, shares="volume", frames=FRAMES):
    return naming.run_name(os.path.basename(PREPARED), "replay", frames, form, dx, brackets(shares), 1)


def run_dir(form, dx, shares="volume", out=None, frames=FRAMES):
    return naming.run_dir(out or os.path.join(REPO, RUNS_OUT), run_name(form, dx, shares, frames))


# ------------------------------------------------------------------------------------------------------- running
def build_bodies():
    for form, dx, _ in FORMS:
        if os.path.isfile(os.path.join(REPO, body_dir(form, dx), "body.json")):
            continue
        cmd = [sys.executable, "-m", "spheral_frag", "body", "--prepared", PREPARED, "--frames",
               "{}-{}".format(*FRAMES), "--form", form, "--dx-mm", str(dx), "--synthetic", "uniform",
               "--uniform-delta-k", "0", "--quiet"]
        subprocess.run(cmd, cwd=REPO, check=True)


def run_command(form, dx, shares, out, frames=FRAMES, max_load=2.0):
    return ["spheral_frag/runner/run.py", "replay", "--body", body_dir(form, dx), "--prepared", PREPARED,
            "--frames", "{}-{}".format(*frames), "--strength", STRENGTH, "--loads-frame", str(LOADS_FRAME),
            "--area-shares", shares, "--ke-fraction", "1e300", "--hourglass", repr(HOURGLASS), "--out", out,
            "--max-load", repr(max_load)]


def launch(form, dx, shares, procs, max_load, force):
    import spheral_m2_check4 as c4
    d = run_dir(form, dx, shares)
    if os.path.isfile(os.path.join(d, "summary.json")):
        return "done"
    if procs > 4 and not force:
        c4.wait_for_load(0.5*max_load)
    cmd = [LAUNCHER] + (["-n", str(procs)] if procs > 1 else []) + run_command(form, dx, shares, RUNS_OUT,
                                                                                max_load=max_load)
    if os.path.isdir(d):
        cmd.append("--resume")
    if force:
        cmd.append("--force")
    t = time.perf_counter()
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    os.makedirs(os.path.join(REPO, OUT, "logs"), exist_ok=True)
    with open(os.path.join(REPO, OUT, "logs", "{}_dx{}_{}.log".format(form, dx, shares)), "w") as fh:
        fh.write(" ".join(cmd) + "\n" + r.stdout + "\n--- stderr\n" + r.stderr)
    if r.returncode != 0:
        print("  FAILED ({}): {}".format(r.returncode, r.stderr.strip().splitlines()[-1:] if r.stderr else ""),
              flush=True)
        return "failed"
    return "ran in {:.0f} s".format(time.perf_counter() - t)


# ---------------------------------------------------------------------------------------- finite-element reference
class Reference:
    """A saved P2 solution of Step 1 and a point locator on its mesh."""

    def __init__(self, path):
        import spheral_m2_elastic_fe as fe6
        from scipy.spatial import cKDTree
        from spheral_frag.geometry import TetLocator
        self.basis, self.u, self.lam, self.mu, self.meta = fe6.load_solution(path)
        p, t = self.basis.mesh.p, self.basis.mesh.t
        self.locator = TetLocator(p.T, t.T)
        self.tree = cKDTree(p[:, t].mean(axis=1).T)

    def cells(self, points):
        """The containing element of each point; outside the faceted ball the element whose centroid is nearest
        (its polynomial extrapolated) -- and the number of such points."""
        c, _ = self.locator.locate(points)
        out = c < 0
        if out.any():
            c[out] = self.tree.query(points[out])[1]
        return c, int(out.sum())

    def stress(self, points, chunk=20000):
        """sigma (n, 3, 3) at points (n, 3) (m): the exact gradient of the P2 field in the element (as Step 1's
        `stress_at`, with this class's element search), and the number of points outside the mesh."""
        points = np.asarray(points, dtype=float)
        out = np.zeros((len(points), 3, 3))
        n_out = 0
        b = self.basis
        for lo in range(0, len(points), chunk):
            x = points[lo:lo + chunk].T
            cells, n = self.cells(x.T)
            n_out += n
            pts = b.mapping.invF(x[:, :, None], tind=cells)
            grad = np.zeros((3, 3, x.shape[1]))
            dofs = b.element_dofs[:, cells]
            for k in range(b.Nbfun):
                g = b.elem.gbasis(b.mapping, pts, k, tind=cells)[0].grad[:, :, :, 0]
                grad += g*self.u[dofs[k]][None, None, :]
            eps = 0.5*(grad + grad.transpose(1, 0, 2))
            tr = eps[0, 0] + eps[1, 1] + eps[2, 2]
            s = 2.0*self.mu*eps
            for i in range(3):
                s[i, i] += self.lam*tr
            out[lo:lo + chunk] = s.transpose(2, 0, 1)
        return out, n_out


_REFS = {}


def reference(path):
    if path not in _REFS:
        _REFS[path] = Reference(os.path.join(REPO, path))
    return _REFS[path]


def fe_summary():
    with open(os.path.join(REPO, FE_DIR, "summary.json")) as fh:
        return json.load(fh)


# ----------------------------------------------------------------------------------------------------- measuring
def points3(form, x_ref):
    """Reference positions as 3D points: 3D as is; RZ (x, r) -> (x, r, 0)."""
    x_ref = np.asarray(x_ref, dtype=float)
    if form == "3d":
        return x_ref
    return np.stack([x_ref[:, 0], x_ref[:, 1], np.zeros(len(x_ref))], axis=1)


def wrms(a, w):
    return math.sqrt(float(np.sum(w*a*a)/np.sum(w))) if len(a) else math.nan


def stress_errors(sig_p, sig_fe, r, V):
    """Max and volume-weighted RMS of vM(sigma_p - sigma_FE) and |dP| over VM0, interior and shell."""
    d = sig_p - sig_fe
    vm = run.von_mises(d)/VM0
    dps = np.trace(d, axis1=1, axis2=2)/3.0/VM0            # signed: tr(sigma_p - sigma_FE) / 3 (= -dP)
    dp = np.abs(dps)
    out = {}
    core, inner = r <= CORE, r <= INTERIOR
    out["core_vm_max"] = float(vm[core].max())
    out["core_vm_rms"] = wrms(vm[core], V[core])
    out["core_p_mean"] = float(np.sum(V[core]*dps[core])/np.sum(V[core]))
    out["int_p_mean"] = float(np.sum(V[inner]*dps[inner])/np.sum(V[inner]))
    for name, sel in (("int", r <= INTERIOR), ("shell", r > INTERIOR)):
        out[name + "_vm_max"] = float(vm[sel].max())
        out[name + "_vm_rms"] = wrms(vm[sel], V[sel])
        out[name + "_p_max"] = float(dp[sel].max())
        out[name + "_p_rms"] = wrms(dp[sel], V[sel])
    return out


def state(d, k, table, form):
    ck = checkpoint.read(os.path.join(d, "checkpoints", checkpoint.filename(k)))
    st = run.stress_state(ck, table)
    order = np.argsort(st["pid"], kind="stable")
    st = {key: (v[order] if isinstance(v, np.ndarray) and len(v) == len(order) else v) for key, v in st.items()}
    st["p3"] = points3(form, st["x_ref"])
    st["r"] = np.linalg.norm(st["p3"], axis=1)
    return ck, st


def measure(d, k, table, form, fine, coarse=None, prev=None):
    """The stress measures at frame k of run directory d (module docstring). `prev` the state at k - 1."""
    ck, st = state(d, k, table, form)
    sig_fe, n_out = fine.stress(st["p3"])
    m = stress_errors(st["sigma"], sig_fe, st["r"], st["V"])
    m.update(n_particles=len(st["pid"]), n_interior=int(np.sum(st["r"] <= INTERIOR)),
             n_shell=int(np.sum(st["r"] > INTERIOR)), n_outside_fe=n_out,
             particle_vm_max=float(st["von_mises"].max())/VM0)
    if coarse is not None:
        sig_c, _ = coarse.stress(st["p3"])
        vmc = run.von_mises(sig_fe - sig_c)/VM0
        m["fe_conv_int_vm_max"] = float(vmc[st["r"] <= INTERIOR].max())
        m["fe_conv_shell_vm_max"] = float(vmc[st["r"] > INTERIOR].max())
    if prev is not None:
        assert np.array_equal(prev["pid"], st["pid"])
        dd = st["sigma"] - prev["sigma"]
        m["settle_vm"] = float(run.von_mises(dd).max())/VM0
        m["settle_p"] = float(np.abs(np.trace(dd, axis1=1, axis2=2)).max())/3.0/VM0
    m["_state"] = st
    m["_fe"] = sig_fe
    return m


def surface(d, k, form, dx, shares="volume", prepared=None):
    """The surface data the loads package builds at the start of interval k -> k + 1, recomputed from checkpoint k
    with the package's functions (module docstring), and the drags of its loads with its normals and with exact
    radial ones."""
    from spheral_frag.runner import nodes as rn
    from spheral_frag.runner.packages import loads as pl
    prepared = os.path.join(REPO, prepared or PREPARED)
    ck = checkpoint.read(os.path.join(d, "checkpoints", checkpoint.filename(k)))
    x = np.asarray(ck["x"], dtype=float)
    V = np.asarray(ck["m"], dtype=float)/np.asarray(ck["rho"], dtype=float)
    h = rn.NPERH*dx
    g, surf, normal = pl.detect_surface(form, x, V, x, h)
    exposed = pl.exposure(form, x, surf, dx)
    area_k = pl.frame_surface_area(prepared, k)
    if shares == "patch":
        A = np.zeros(len(x))
        A[surf] = pl.patch_shares(form, np.asarray(ck["x_ref"], dtype=float)[surf], *pl.frame_patches(prepared, k))
    else:
        raw = pl.raw_shares(form, x, V)
        A = np.where(surf, raw*area_k/float(np.sum(raw[surf])), 0.0)
    table = pl.lee_table(pl.read_tables(prepared, [LOADS_FRAME])[0], None)
    idx = np.flatnonzero(surf)
    th = pl.inclination_deg(normal[idx])
    p, tau = pl.particle_loads(table, th, exposed[idx])
    F = pl.surface_forces(normal[idx], A[idx], p, tau)
    p3 = points3(form, ck["x_ref"])[idx]
    radial = p3/np.linalg.norm(p3, axis=1)[:, None]
    n3 = points3(form, normal[idx]) if form == "rz" else normal[idx]
    ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", n3, radial), -1.0, 1.0)))
    rad = radial if form == "3d" else radial[:, :2]
    pr, taur = pl.particle_loads(table, pl.inclination_deg(rad), exposed[idx])
    Fr = pl.surface_forces(rad, A[idx], pr, taur)
    An = A[idx][:, None]*normal[idx]
    closure = np.sum(An, axis=0)
    if form == "rz":
        closure = closure[:1]                                 # rings: only the axial sum is a net vector
    nx = normal[idx, 0]
    if form == "3d":                         # the loads' moment about the body's mass centre (RZ: zero by symmetry)
        m_ = np.asarray(ck["m"], dtype=float)
        xc = (m_[:, None]*x).sum(axis=0)/m_.sum()
        moment = float(np.linalg.norm(np.cross(x[idx] - xc, F).sum(axis=0)))/(abs(float(F[:, 0].sum()))*R)
    else:
        moment = 0.0
    return dict(moment_rel=moment, n_surface=int(surf.sum()), n_exposed=int(exposed.sum()), area_frame_m2=area_k,
                n_zero_area=int(np.count_nonzero(A[idx] == 0.0)),
                area_cover=float(A.sum())/(4.0*math.pi*R*R),
                closure=float(np.linalg.norm(closure))/(4.0*math.pi*R*R),
                projected_windward=float(np.sum(A[idx][nx > 0]*nx[nx > 0]))/(math.pi*R*R),
                normal_err_max_deg=float(ang.max()), normal_err_rms_deg=float(np.sqrt(np.mean(ang*ang))),
                drag_surface_N=-float(F[:, 0].sum()), drag_radial_normals_N=-float(Fr[:, 0].sum()),
                _surface=dict(x=x[idx], normal=normal[idx], theta=th, area=A[idx], p=p, tau=tau))


def section(st, form, dx, sf):
    """(the particles' axial force through the mid-plane x = 0, the load the nose half x > 0 transmits there): the
    first sum sigma_xx over the particles of the lattice plane x_ref = 0 times their section (3D dx^2, RZ the ring's
    2 pi r dx), the second the surface forces on x > 0 (half on x = 0) less the deceleration's inertial force on the
    same mass, -a_frame sum m; equilibrium makes them equal (the finite element's section integral is -15.48 N)."""
    from spheral_frag.runner.packages import loads as pl
    row = np.abs(st["x_ref"][:, 0]) < 1e-6*dx
    if form == "3d":
        fsec = float(np.sum(st["sigma"][row, 0, 0]))*dx*dx
    else:
        fsec = float(np.sum(st["sigma"][row, 0, 0]*2.0*math.pi*st["x_ref"][row, 1]))*dx
    F = pl.surface_forces(sf["normal"], sf["area"], sf["p"], sf["tau"])[:, 0]
    xs = sf["x"][:, 0]
    w = np.where(xs > 1e-6*dx, 1.0, np.where(np.abs(xs) <= 1e-6*dx, 0.5, 0.0))
    xm = st["x_ref"][:, 0]
    wm = np.where(xm > 1e-6*dx, 1.0, np.where(row, 0.5, 0.0))
    a = float(F.sum())/float(np.sum(st["m"]))
    return fsec, float(np.sum(w*F)) - a*float(np.sum(wm*st["m"]))


def analyse(only=None):
    table = material.MaterialTable.load(os.path.join(REPO, PREPARED))
    fes = fe_summary()
    drag_ref = float(fes["tables"]["lee"]["continuum_drag_N"])
    fe_decel = -float(fes["meshes"]["1.6667"]["cases"]["lee"]["a_ms2"][0])
    fine, coarse = reference(FE_FINE), reference(FE_COARSE)
    rows, profiles = [], {}
    for shares in SHARES:
        for form, dx_mm, procs in FORMS:
            if only and (form, dx_mm, shares) not in only:
                continue
            d = run_dir(form, dx_mm, shares)
            if not os.path.isfile(os.path.join(d, "summary.json")):
                print("  missing: {}".format(naming.run_name_from_dir(d)))
                continue
            with open(os.path.join(d, "summary.json")) as fh:
                summ = json.load(fh)
            rows += measure_run(d, form, dx_mm, shares, table, fine, coarse, drag_ref, fe_decel, summ["processes"],
                                profiles)
    path = os.path.join(REPO, CSV)
    if not only:
        write_csv(path, rows)
        print("wrote {} ({} rows)".format(CSV, len(rows)))
    report(rows)
    plots(profiles, fine)
    return rows


def measure_run(d, form, dx_mm, shares, table, fine, coarse, drag_ref, fe_decel, processes, profiles=None):
    iv = {int(r["k1"]): r for r in run.read_intervals(os.path.join(d, "intervals.csv"))}
    sf = surface(d, FRAMES[0] + 1, form, dx_mm*1e-3, shares)
    rows = []
    _, prev = state(d, FRAMES[0], table, form)
    for k in range(FRAMES[0] + 1, FRAMES[1] + 1):
        m = measure(d, k, table, form, fine, coarse if k == COMPARE_K else None, prev)
        prev = m.pop("_state")
        sig_fe = m.pop("_fe")
        r = iv[k]

        def num(c):
            v = r.get(c, "")
            return float(v) if v not in ("", None) else math.nan
        mass = float(np.sum(prev["m"]))
        hist = np.load(os.path.join(d, "history", "interval_{:05d}.npz".format(k - 1)))
        s = float(hist["scale"])
        h = hist["history"]
        decel = -s*float(h[-1, 1 + (2 if form == "rz" else 3)])
        drag = -float(h[-1, 1])
        row = dict(form=form, dx_mm=dx_mm, area_shares=shares, k=k, processes=processes, scale=num("scale"),
                   ke_ratio=num("ke_ratio"), ke_internal_ratio=num("ke_internal_ratio"),
                   drag_N=drag, drag_ref_N=drag_ref, drag_rel=drag/drag_ref - 1.0,
                   frame_table_drag_N=num("table_drag_N"), decel_ms2=decel, decel_ref_ms2=drag_ref/mass,
                   decel_rel=decel/(drag_ref/mass) - 1.0, fe_decel_ms2=fe_decel,
                   flight_decel_ms2=num("flight_decel_ms2"), mass_kg=mass,
                   momentum_axial_Ns=num("momentum_axial_Ns"), momentum_over_impulse=num("momentum_over_impulse"),
                   v_com_removed_ms=num("v_com_removed_ms"), recentre_shift_dx=num("recentre_shift_dx"),
                   max_disp_ref_dx=num("max_disp_ref_dx"),
                   monitored=int(num("monitored")), crossings_err_est=num("crossings_err_est"),
                   crossings_err_pressure=num("crossings_err_pressure"),
                   crossings_err_limit=num("crossings_err_limit"), run_crossings=int(num("run_crossings")),
                   steps=int(num("steps")), wall_s=num("wall_s"), wall_per_step_ms=1e3*num("wall_per_step_s"),
                   run=naming.run_name_from_dir(d))
        row.update({c: m.get(c, math.nan) for c in COLUMNS if c in m})
        row.update({c: sf[c] for c in sf if not c.startswith("_")})
        fsec, fapp = section(prev, form, dx_mm*1e-3, sf["_surface"])
        row.update(section_force_N=fsec, section_applied_N=fapp, section_ratio=fsec/fapp)
        if k == FRAMES[0] + 2 and abs(drag/sf["drag_surface_N"] - 1.0) > 1e-12:     # built from checkpoint k - 1
            raise RuntimeError("{} k {}: the recomputed surface's drag {} differs from the run's {}".format(
                form, k, sf["drag_surface_N"], drag))
        if int(num("n_surface")) != sf["n_surface"] or int(num("n_exposed")) != sf["n_exposed"]:
            raise RuntimeError("{} k {}: the recomputed surface ({}, {}) differs from the run's ({}, {})".format(
                form, k, sf["n_surface"], sf["n_exposed"], r["n_surface"], r["n_exposed"]))
        rows.append(row)
        if profiles is not None and k == COMPARE_K:
            profiles[(form, dx_mm, shares)] = dict(state=prev, fe=sig_fe, surface=sf["_surface"])
    return rows


def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(COLUMNS)
        for row in rows:
            row = dict({c: math.nan for c in COLUMNS}, **row)      # fe_conv_*: COMPARE_K only; settle_*: k > K0
            wr.writerow([repr(float(row[c])) if isinstance(row[c], (float, np.floating)) else row[c]
                         for c in COLUMNS])


def report(rows):
    print("\n{:6s} {:3s} {:>4s} {:>2s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
        "shares", "frm", "dx", "k", "int_vm", "int_rms", "int_p", "sh_vm", "sh_rms", "sh_p", "settle", "KEint/E",
        "drag_rel", "decel"))
    for r in rows:
        print("{:6s} {:3s} {:4.1f} {:2d} {:9.3e} {:9.3e} {:9.3e} {:9.3e} {:9.3e} {:9.3e} {:9.2e} {:9.2e} {:+9.4f} {:9.4f}"
              .format(r["area_shares"], r["form"], r["dx_mm"], r["k"], r["int_vm_max"], r["int_vm_rms"], r["int_p_max"],
                      r["shell_vm_max"], r["shell_vm_rms"], r["shell_p_max"], r.get("settle_vm", math.nan),
                      r["ke_internal_ratio"], r["drag_rel"], r["decel_ms2"]))


# --------------------------------------------------------------------------------------------------------- plots
INK, MUTED = "#1f2328", "#6e7781"
COMP = {"xx": (0, 0, "#0969da"), "yy": (1, 1, "#bf3989"), "zz": (2, 2, "#1a7f37")}


def line_particles(form, p3, dx, line):
    """Particles on the axis (3D: y = z = 0; RZ: the row nearest the axis, r = dx/2) or the 90-degree meridian
    (3D: x = z = 0, y > 0; RZ: the row x = 0), sorted along it; the coordinate along it (signed x, or y)."""
    tol = 1e-3*dx
    if line == "axis":
        if form == "3d":
            sel = (np.abs(p3[:, 1]) < tol) & (np.abs(p3[:, 2]) < tol)
        else:
            sel = np.abs(p3[:, 1] - p3[:, 1].min()) < tol
        s = p3[:, 0]
    else:
        sel = (np.abs(p3[:, 0]) < tol) & (np.abs(p3[:, 2]) < tol) & (p3[:, 1] > 0.0)
        s = p3[:, 1]
    idx = np.flatnonzero(sel)
    return idx[np.argsort(s[idx])], s


def plots(profiles, fine):
    if not profiles:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = [(f, dx) for f, dx, _ in FORMS if any((f, dx, sh) in profiles for sh in SHARES)]
    rows_ = [sh for sh in SHARES if any((f, dx, sh) in profiles for f, dx in cols)]
    os.makedirs(os.path.join(REPO, OUT), exist_ok=True)
    for line, xlabel in (("axis", "x (mm): tail -50 ... nose +50"), ("meridian", "y (mm), x = 0")):
        fig, axes = plt.subplots(len(rows_), len(cols), figsize=(4.2*len(cols), 3.8*len(rows_)), sharey=True,
                                 squeeze=False)
        for a_row, sh in zip(axes, rows_):
            for ax, (form, dx_mm) in zip(a_row, cols):
                key = (form, dx_mm, sh)
                if key not in profiles:
                    ax.set_visible(False)
                    continue
                st = profiles[key]["state"]
                idx, s = line_particles(form, st["p3"], dx_mm*1e-3, line)
                off = st["p3"][idx[0]].copy() if len(idx) else np.zeros(3)
                t = np.linspace(-R + 1e-4, R - 1e-4, 401) if line == "axis" else np.linspace(0.0, R - 1e-4, 201)
                pts = np.zeros((len(t), 3))
                if line == "axis":
                    pts[:, 0] = t
                    pts[:, 1] = off[1]
                else:
                    pts[:, 1] = t
                sig_line, _ = fine.stress(pts)
                for name, (i, j, color) in COMP.items():
                    ax.plot(t*1e3, sig_line[:, i, j], color=color, lw=1.4, label="FE " + name)
                    ax.plot(s[idx]*1e3, st["sigma"][idx, i, j], "o", color=color, ms=3.5, mfc="none",
                            label="particles " + name)
                ax.set_title("{} {} mm, {} shares".format(form.upper(), dx_mm, sh), fontsize=10)
                ax.set_xlabel(xlabel)
                ax.grid(color="#d0d7de", lw=0.5)
            a_row[0].set_ylabel("stress (Pa); RZ zz = hoop")
        axes[0][0].legend(fontsize=7, frameon=False, ncol=2)
        fig.suptitle("Check 6: frame 182's loads on the elastic sphere, {} (k {})".format(
            "along the flight axis" if line == "axis" else "along the 90-degree meridian", COMPARE_K), fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(REPO, OUT, "profile_{}.png".format(line)), dpi=110)
        plt.close(fig)
    # von Mises of the difference against the reference radius, per run
    shades = {("rz", 2.2): "#c6dbef", ("rz", 1.1): "#2171b5", ("3d", 3.0): "#a1d99b", ("3d", 2.2): "#238b45"}
    fig, axes = plt.subplots(1, len(rows_), figsize=(6.0*len(rows_), 4.2), sharey=True, squeeze=False)
    for ax, sh in zip(axes[0], rows_):
        for form, dx_mm in cols:
            key = (form, dx_mm, sh)
            if key not in profiles:
                continue
            st, fe_p = profiles[key]["state"], profiles[key]["fe"]
            ax.plot(st["r"]*1e3, run.von_mises(st["sigma"] - fe_p)/VM0, ".", ms=2, color=shades[(form, dx_mm)],
                    label="{} {} mm".format(form.upper(), dx_mm))
        ax.axvline(INTERIOR*1e3, color=MUTED, lw=0.8, ls="--")
        ax.set_xlabel("|x_ref| (mm)")
        ax.set_title("{} shares".format(sh), fontsize=10)
        ax.set_yscale("log")
        ax.grid(color="#d0d7de", lw=0.5)
        ax.legend(fontsize=8, frameon=False, markerscale=4)
    axes[0][0].set_ylabel("vM(sigma_particle - sigma_FE) / 3,510.5 Pa")
    fig.tight_layout()
    fig.savefig(os.path.join(REPO, OUT, "error_vs_radius.png"), dpi=110)
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--only", choices=("run", "analyse"))
    p.add_argument("--max-load", type=float, default=2.0)
    p.add_argument("--force", action="store_true", help="start runs even above --max-load")
    a = p.parse_args(argv)
    os.chdir(REPO)
    if a.only in (None, "run"):
        build_bodies()
        for shares in SHARES:
            for form, dx, procs in FORMS:
                print("{} {} mm, {} shares: ".format(form, dx, shares), end="", flush=True)
                print(launch(form, dx, shares, procs, a.max_load, a.force), flush=True)
    if a.only in (None, "analyse"):
        analyse()
    return 0


if __name__ == "__main__":
    sys.exit(main())
