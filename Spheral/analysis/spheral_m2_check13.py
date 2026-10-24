"""Spheral M2 check 13: RZ against 3D on the real flight's first interval range (M2 plan, Task 10; spec §12 check 13).

    "$PY" analysis/spheral_m2_check13.py [--only run|probe|analyse] [--forms rz:2.2,...] [--max-load 2.0] [--force]

**The case** (plan fact 8, decisions 1, 7, 8, 9). The MVP flight's 100 mm sphere from frame 0, replayed over frames
0-47 (0-23.5 s: heating and thermal stress before the first film at k 48; no particle leaves the body), every
setting at its default: `--strength flight` (Chen's law, the provisional extension below 623 K), the loads and the
deceleration on with the patch area shares, SolidSPH, 10 crossings per interval with the crossings monitor at its
defaults (the first interval and every 10th after it, j = 9, 19, 29, 39; limit 5 % on the larger of the von Mises
and pressure estimates), damping 0.05. RZ at 2.2 mm (813 particles, 1 process) and 1.1 mm (3,243, 4 processes), 3D
at 2.2 mm (49,173, 18 processes); the bodies are Task 2's (rebuilt byte-identical on 2026-10-11). The checkpoints
of every frame stay in the run directories (`checkpoints/checkpoint_<k>.npz`): Task 11 switches from k 47.

**Matching** (`match`). RZ and 3D are compared at the RZ particles' reference positions (x, r). In 3D every particle's
stress is resolved into (axial, radial, hoop) at its own azimuth (`cylindrical`); the 3D particles nearest an RZ
ring in the (x, rho) half-plane -- all particles at the nearest (x, rho), i.e. every azimuth of that lattice ring, 4
to 16 of them -- are averaged with volume weights ("averaged over azimuth"), at most 0.5 spacings from the ring. RZ
2.2 mm against RZ 1.1 mm the same way (the four nearest 1.1 mm rings). Compared: the von Mises stress of the
averaged tensor and the pressure, as target minus source, over the 3D run's largest von Mises at k (`vm_peak_3d_Pa`,
the same scale for all pairs: 94.6, 73.3, 63.1, 58.9 MPa at k 12, 24, 36, 47); the maximum, the volume-weighted RMS
and the signed mean of each per region; `vmt_*` the von Mises of the tensor difference; `dT_max_K` the largest
difference of the imposed temperatures at the matched positions (the inputs' own mismatch: the ring positions differ
by up to half a spacing).

**Regions** (`regions`, by the target's reference radius, with dx_c the coarser spacing of the pair): `axis` the RZ
row nearest the axis (r = dx/2) and `axis_interior` its part with |x| <= R - 2 dx_c; `surface` R - dx_c < |x_ref|
(the outermost spacing), `subsurface` the next spacing, `shell` both (the outer one-two spacings that carry SPH's
free-surface error, checks 4 and 6), `interior` |x_ref| <= R - 2 dx_c, `core` |x_ref| <= 25 mm, `all`.

**The floor** (Task 2): the body stage samples each RZ particle's temperature at azimuth 0 and records the largest
difference to the samples at 90, 180 and 270 degrees per frame (the frames' departure from axisymmetry): 2.3-2.6 K
(RZ 2.2 mm) and 2.5-2.7 K (RZ 1.1 mm) at the compared frames. As a stress, E alpha dT / (1 - nu) (fully restrained;
E 71.7 GPa, nu 0.33, alpha_eff 2.34e-5) it is 5.7-6.7 MPa, 7-11 % of the peak (`floor_stress`): an RZ run cannot be
expected to agree with the 3D one closer than that where the temperature field is not axisymmetric. Measured beside
it: the 3D run's own spread of von Mises over the azimuths of a ring (`spread_*`; 3-6 % of the peak in the core, up
to 17-19 % in the interior next to the hot nose's surface).

**Measured** (2026-10-11; the CSVs; thresholds in tests/test_spheral_frag_runner_check13.py). Fractions of the 3D
peak von Mises at k.

- *RZ 2.2 mm against 3D 2.2 mm agrees to the floor in the core and along the axis.* Core max |dvM| 5.0-12.8 % and
  |dP| 3.2-6.6 %, RMS 2.1-5.9 % and 1.3-2.5 %; along the axis (interior part) |dvM| 6.6-11.8 %, |dP| 5.2-6.2 % at
  k 12 and 24, rising to 19-30 % at k 36 and 47 at the nose end of the axis row (x 40-46 mm, 750-775 K: the eutectic
  melting at the solidus, SolidSPHRZ's axis row, and the crossings error there, below). Interior max |dvM| 11-16 %,
  |dP| 10-30 % (the same nose particles), RMS 2.4-3.9 % and 2.3-3.3 %. The surface (outer 2.2 mm) |dvM| 28-33 %,
  |dP| 18-28 %, RMS 10-13 % and 7.5-9.2 %; the subsurface |dvM| 17-20 %, |dP| 19-33 %: both forms carry SPH's
  free-surface error in their outer one or two spacings (check 4: 5.5-9.8 % of the exact peak in each form, of
  different sign and shape; check 6: up to 1.4-3 times the peak under loads), and here the peak itself sits in
  that layer, so the surface comparison measures the two forms' surface errors against each other, not either
  against the truth.
- *RZ's known offset* (check 6: SolidSPHRZ carries 5-8 % more axial force through the mid-plane than equilibrium
  requires, 3D 2 % less; check 4: a near-uniform isotropic core offset, RZ +1.4 %, 3D -1.6 to -2.4 %). Here the core's
  mean axial stress is 4.0-4.7 % larger in RZ than in 3D (`core_sxx_ratio` 1.040-1.047; 49.8 against 47.8 MPa
  at k 12) and the core's mean pressure difference is -0.8 to -1.3 % of the peak (RZ more tensile): the same sign as
  both checks, within the floor. Check 6's offset scaled to these loads (its 8 % of the load stress, 3.5 kPa x
  drag / 30.7 N) is at most 1.5e-6 of the peak (`rz_load_offset`): the loads' share is negligible; the offset is
  SolidSPHRZ's under the thermal stress itself.
- *RZ 1.1 mm is not usable in the core after k 5: SolidSPHRZ's tensile pairing.* The core is in axial tension (the
  heated shell expands; core sigma_xx 28 MPa at k 3, 45-50 MPa from k 9 in RZ 2.2 mm and 3D), and at 1.1 mm the
  particles along every row within r < 33 mm pair up axially: the second difference of the axial displacement along a
  row (`pairing`, in spacings) grows from 0.0017 (k 3) by about 3 per interval to 0.050 at k 7, 0.067 at k 9 and 0.10
  at k 47 (no two particles closer than 0.91 dx: not clumping as check 7 counted it), and the core loses its axial
  tension: 20.0 MPa at k 9 against RZ 2.2 mm's 47.8 and 3D's 45.9, 7.2 against 31.0 and 29.6 at k 47 (`core_sxx_ratio`
  against 3D 0.24-0.34 at k 12-47). Its core pressure is off by +22 to +31 % of the peak from 3D's and -32 to -46 %
  from RZ 2.2 mm's. The pairing grows per step, not per unit of flight time: rerunning interval 4 -> 5 from the run's
  checkpoint (probe `pairing_rz11_k4`, 1 process) multiplies it by 3.2 at 10 crossings, 7.8 at 20 and 17 at 40, so the
  monitor's switches to more crossings feed it. At 2.2 mm the same zigzag is weak: RZ 0.0047 and 3D 0.0083 at most up
  to k 36 (on the lines next to the axis from k 24), rising to 0.013 (RZ, k 44) and 0.015 (3D, k 42) after the
  switches to more crossings (RZ to 20 and 40, 3D to 20; 3D's k 39 probe: 0.0086 -> 0.010, 0.014, 0.017 at 20, 40,
  80), with the core's axial stress still agreeing to 4-5 %: the same instability, slower at the coarser spacing
  (fewer steps per interval), and a risk for long replays at more crossings. Spheral's Monaghan artificial stress
  (SPH's `epsTensile`, tried at 0.3 on the same interval) slowed a single interval beyond 20 minutes on one process
  (60 s without it) and was not pursued. The RZ 1.1 mm rows are in the tables, flagged by the tests, not held to the
  RZ-against-3D thresholds.
- *RZ 2.2 mm against RZ 1.1 mm* therefore differs by the pairing in the core (|dP| up to 54-87 %, mean -32 to -46 %)
  and by 32-90 % of the peak in von Mises at its largest; it is no spacing study of the axisymmetric form.

**The crossings monitor** (Asha's rule: the crossings error below 5 % over the flight). Estimates, the larger of
von Mises and pressure, on the monitored intervals 0, 9, 19, 29, 39:
RZ 2.2 mm 1.09, 0.85, 0.54, 6.43 (pressure; switched to 20), 9.48 % (pressure, 20 against 40; switched to 40);
RZ 1.1 mm 1.04, 7.96 (von Mises; switched to 20), 2.85, 7.86 (von Mises; switched to 40), 2.80 %;
3D 2.2 mm 0.66, 0.25, 0.42, 3.66, 26.1 % (pressure; von Mises 16.9 %; switched to 20).
Below 5 % while the body is elastic (k < 9); above it once it yields (RZ 1.1 mm's outer layer at 475 K from k 9;
the nose at 730-790 K from k 29, the eutectic's 0.6 % volume change at the 750 K solidus and the semi-solid flow
law above it). The probes (`PROBES`, `analysis/spheral_m2_check13_probe.py`: the interval rerun from the run's own
checkpoint at more crossings) show the estimates are real and converge: RZ 2.2 mm k 29 -> 30 against 160
crossings 5.1 / 5.8 % (von Mises / pressure) at 10, 3.0 / 2.6 % at 20, 1.4 / 1.2 % at 40; k 39 -> 40 against 320
9.0 / 14.3 % at 10, 5.9 / 5.7 % at 20, 3.8 / 2.9 % at 40, 0.8 / 0.6 % at 160; RZ 1.1 mm k 9 -> 10 against 80
7.0 / 5.5 % at 10, 3.1 / 2.4 % at 20; 3D 2.2 mm k 39 -> 40 against 80 11.0 / 15.0 % at 20 and 4.5 / 5.1 % at 40
(its rerun at 20 equals the run's accepted state bitwise; the monitor's estimate at 20 against 40 is 14.5 / 19.7 %),
largest on the axis at x 41.8 mm, 751 K -- the eutectic, in 3D as in RZ, so not SolidSPHRZ's axis row. The rule as
built leaves gaps: an interval between two
monitored ones is not checked (RZ 2.2 mm ran 30 -> 39 at 20 crossings, whose error at k 39 is 5.9 %; 3D ran 30 -> 39
at 10, whose error at k 39 is 17-26 %), and after the last switch nothing confirms the doubled crossings are below
5 % (3D's last eight intervals ran at 20 crossings, whose error at k 39 is 11-15 % by the probe).

**Cost per interval** against spec §5's gate table (10 crossings, the median step; the gate's "with the Python
classes" column per interval = its hot-phase hours / 399). At 10 crossings the steps match the gate (RZ 2.2 mm
1,649-1,661 against 1,660; RZ 1.1 mm 3,332-3,346 against 3,330; 3D 1,608-1,625 against 1,620) and so does the cost:
RZ 2.2 mm 4.1 ms/step, 6.8 s per interval (gate 9.0 s); RZ 1.1 mm on 4 processes 6.1 ms/step, 20.4 s (gate 25.3 s on
6); 3D 2.2 mm on 18 128 ms/step, 207 s (gate 219 s). The monitor and its switches multiply it in the plastic phase:
with them, the mean over the run is 14.0 s (RZ 2.2 mm), 69.9 s (RZ 1.1 mm) and 281 s (3D) per interval, 1.56, 2.77
and 1.28 times the gate's -- at that rate the 399-interval hot phase would take 1.6, 7.7 and 31 h against the gate's
1.0, 2.8 and 24.3 h, and more where the hot phase needs 40 crossings throughout (3D at 40: about 820 s per interval,
91 h). The 3D run's
first three intervals ran beside this study's own 1-process probes and slowed to 196-252 ms/step; they are flagged
(`disturbed`) and left out of the cost measures (their stresses are unaffected).

**Probes** (`--only probe`): the reruns above write to `spheral_output/m2_probe/task10/` (git-ignored); the table
`data/m2_check13_probe.csv` has one row per probe and crossings number: the monitor's estimate against 2 N,
the error against the most crossings and where it is largest, the plastic particles, the pairing amplitude at the
start and end and the core's axial stress.

Outputs: `data/m2_check13.csv` (one row per pair, compared frame and region), `m2_check13_runs.csv` (one row
per run and interval), `m2_check13_probe.csv`; the runs in `spheral_output/m2/check13/runs/`, their logs in
`spheral_output/m2/check13/logs/`, the profile plots `spheral_output/m2/check13/profile_{axis,surface}.png`.
Resumable: bodies and runs whose outputs exist are skipped (run.py's summary.json); an interrupted run resumes from
its last checkpoint. Run under `caffeinate -i -s` (the host's sleep suspends Docker's VM), one long run at a time,
and nothing else on the VM while the 3D run is timed (a single extra process slows its 18 ranks by up to 2 times).
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
# the hourglass stiffness (Task 10b): Task 10's runs were made before it, without it (_hourglass-0 in their names);
# Task 10b's reruns, after the controller's review, set it to the default
HOURGLASS = 0.0
FRAMES = (0, 47)                            # before the first film (k 48; plan fact 8)
COMPARE_K = (12, 24, 36, 47)
OUT = os.path.join("spheral_output", "m2", "check13")
RUNS_OUT = os.path.join(OUT, "runs")
CSV = os.path.join("data", "m2_check13.csv")
LAUNCHER = os.path.join("spheral_frag", "container", "spheral")
R = 0.05                                    # m

FORMS = [("rz", 2.2, 1), ("rz", 1.1, 4), ("3d", 2.2, 18)]      # form, dx mm, processes (as checks 4-6)


def body_dir(form, dx):
    return naming.body_dir(os.path.join("spheral_output", "body"),
                           naming.body_name(os.path.basename(PREPARED), *FRAMES, form, dx))


def run_name(form, dx):
    """Every setting at its default (decisions 1, 7, 8, 9) but the hourglass (`HOURGLASS`)."""
    return naming.run_name(os.path.basename(PREPARED), "replay", FRAMES, form, dx,
                           naming.canonical_brackets({"hourglass": HOURGLASS}), 1)


def run_dir(form, dx, out=None):
    return naming.run_dir(out or os.path.join(REPO, RUNS_OUT), run_name(form, dx))


# ------------------------------------------------------------------------------------------------------ measuring
PAIRS = [(("rz", 2.2), ("3d", 2.2)),         # (target, source): the source mapped onto the target's particles
         (("rz", 1.1), ("3d", 2.2)),
         (("rz", 2.2), ("rz", 1.1))]
CORE = 0.025                                # m: |x_ref| <= this is the core (as check 6)
REGIONS = ("axis", "axis_interior", "surface", "subsurface", "shell", "interior", "core", "all")
VM0_CHECK6 = 3510.5                         # Pa: check 6's reference peak von Mises under frame 182's loads
DRAG_CHECK6 = 30.717                        # N: frame 182's table drag over the exact sphere (check 6)
RZ_SECTION_EXCESS = 0.08                    # check 6: SolidSPHRZ carries 5.0-7.7 % more axial force (upper, rounded)
# spec §5's gate table (the 100 mm hot phase, 399 intervals of 0.5 s, 10 crossings; "with the Python classes")
GATE = {("3d", 2.2): dict(particles=49173, processes=18, step_ms=113.0, steps=1620, wall_h=20.3, wall_py_h=24.3),
        ("rz", 2.2): dict(particles=813, processes=1, step_ms=3.0, steps=1660, wall_h=0.6, wall_py_h=1.0),
        ("rz", 1.1): dict(particles=3243, processes=6, step_ms=4.2, steps=3330, wall_h=1.6, wall_py_h=2.8)}
GATE_INTERVALS = 399
# 3D intervals whose wall time the study's own probes slowed (1-process probes beside the 18-rank run, 2026-10-11:
# 196-252 ms/step against 125-133); their stresses are unaffected, their cost is left out of the cost measures
DISTURBED = {("3d", 2.2): (0, 1, 2)}

COLUMNS = ["target", "source", "k", "t_s", "region", "n", "vm_peak_3d_Pa", "vm_max", "vm_rms", "vm_mean",
           "p_max", "p_rms", "p_mean", "vmt_max", "vmt_rms", "dT_max_K", "spread_max", "spread_rms",
           "match_dist_max_dx", "n_ties_min", "floor_azimuth_K", "floor_stress", "rz_load_offset",
           "target_vm_max", "target_p_absmax", "core_sxx_target_Pa", "core_sxx_source_Pa", "core_sxx_ratio",
           "run_target", "run_source"]
RUN_COLUMNS = ["form", "dx_mm", "processes", "k0", "k1", "t1_s", "n_particles", "scale", "crossings", "run_crossings",
               "repeats", "steps", "wall_s", "wall_per_step_ms", "ke_internal_ratio", "ke_ratio", "monitored",
               "crossings_err_est", "crossings_err_pressure", "crossings_err_limit", "crossings_err_limit_basis",
               "crossings_switch", "monitor_wall_s", "drag_mean_N", "table_drag_N", "decel_mean_ms2",
               "flight_decel_ms2", "max_disp_ref_dx", "recentre_shift_dx", "vm_peak_Pa", "p_absmax_Pa", "T_max_K",
               "pairing_dx", "core_sxx_Pa", "disturbed",
               "gate_steps", "gate_step_ms", "gate_interval_s", "gate_interval_py_s", "run"]
RUNS_CSV = os.path.join("data", "m2_check13_runs.csv")


def cylindrical(xr, sigma):
    """(x, rho, sigma in (axial, radial, hoop) rows/columns) of 3D particles at reference positions xr (n, 3); on the
    axis the radial direction is taken as +y (any is the same for an axisymmetric field)."""
    y, z = xr[:, 1], xr[:, 2]
    rho = np.hypot(y, z)
    c = np.where(rho > 0.0, y/np.where(rho > 0.0, rho, 1.0), 1.0)
    s = np.where(rho > 0.0, z/np.where(rho > 0.0, rho, 1.0), 0.0)
    Q = np.zeros((len(xr), 3, 3))
    Q[:, 0, 0] = 1.0
    Q[:, 1, 1], Q[:, 1, 2] = c, s
    Q[:, 2, 1], Q[:, 2, 2] = -s, c
    return xr[:, 0], rho, np.einsum("nij,njk,nlk->nil", Q, sigma, Q)


def load_state(form, dx, k, table, out=None):
    """The run's state at frame k (pid order): x, rho (reference axial and radial positions), sigma in (axial,
    radial, hoop), von Mises, P, T, V."""
    d = run_dir(form, dx, out)
    ck = checkpoint.read(os.path.join(d, "checkpoints", checkpoint.filename(k)))
    st = run.stress_state(ck, table)
    o = np.argsort(st["pid"], kind="stable")
    xr = np.asarray(st["x_ref"], dtype=float)[o]
    sig = np.asarray(st["sigma"], dtype=float)[o]
    if form == "3d":
        x, rho, sig = cylindrical(xr, sig)
    else:                                     # RZ rows/columns are already (axial, radial, hoop)
        x, rho = xr[:, 0], xr[:, 1]
    return dict(form=form, dx=dx*1e-3, pid=np.asarray(st["pid"])[o], x=x, rho=rho, sigma=sig, x_ref=xr,
                x_cur=np.asarray(st["x"], dtype=float)[o],
                vm=run.von_mises(sig), P=np.asarray(st["P"], dtype=float)[o], T=np.asarray(st["T"], dtype=float)[o],
                V=np.asarray(st["V"], dtype=float)[o], run=naming.run_name_from_dir(d))


def match(tgt, src):
    """Map the source onto the target's particles: for each target particle (x, r) the source particles nearest it
    in the (axial, radial) half-plane -- in 3D every particle of the nearest (x, rho) ring position, i.e. all its
    azimuths; in RZ the nearest rings, ties included -- averaged with volume weights. Returns the averaged tensor,
    P, T, the spread of the source's von Mises over the tied particles (3D: over azimuth), the match distance and
    the number of tied particles."""
    from scipy.spatial import cKDTree
    pts = np.stack([src["x"], src["rho"]], axis=1)
    tree = cKDTree(pts)
    q = np.stack([tgt["x"], tgt["rho"]], axis=1)
    tol = 1e-6*src["dx"]
    kq = 16
    while True:
        dist, idx = tree.query(q, k=min(kq, len(pts)))
        if kq >= len(pts) or np.all(dist[:, -1] > dist[:, 0] + tol):
            break
        kq *= 2
    tie = dist <= dist[:, :1] + tol
    w = np.where(tie, src["V"][idx], 0.0)
    w /= w.sum(axis=1, keepdims=True)
    sig = np.einsum("nk,nkij->nij", w, src["sigma"][idx])
    P = np.einsum("nk,nk->n", w, src["P"][idx])
    T = np.einsum("nk,nk->n", w, src["T"][idx])
    vm = np.where(tie, src["vm"][idx], np.nan)
    spread = np.nanmax(vm, axis=1) - np.nanmin(vm, axis=1)
    return dict(sigma=sig, P=P, T=T, vm=run.von_mises(sig), spread=spread, dist=dist[:, 0], n_ties=tie.sum(axis=1))


def regions(tgt, dx_c):
    """Boolean masks over the target's particles (module docstring, "Regions")."""
    rr = np.hypot(tgt["x"], tgt["rho"])
    axis = np.abs(tgt["rho"] - tgt["rho"].min()) < 1e-6*tgt["dx"]
    return {"axis": axis, "axis_interior": axis & (np.abs(tgt["x"]) <= R - 2.0*dx_c), "surface": rr > R - dx_c,
            "subsurface": (rr > R - 2.0*dx_c) & (rr <= R - dx_c), "shell": rr > R - 2.0*dx_c,
            "interior": rr <= R - 2.0*dx_c, "core": rr <= CORE,
            "all": np.ones(len(rr), dtype=bool)}


def wrms(a, w):
    return math.sqrt(float(np.sum(w*a*a)/np.sum(w))) if len(a) else math.nan


def body_azimuth_spread(form, dx):
    """{k: the body stage's largest temperature difference between azimuths 0, 90, 180, 270 deg} (Task 2)."""
    with open(os.path.join(REPO, body_dir(form, dx), "body.json")) as fh:
        doc = json.load(fh)
    return {int(f["k"]): float(f["azimuth_spread_K"]) for f in doc["per_frame"]}


def compare(tgt, src, k, t_s, scale, floors, thermo, load_stress):
    """One row per region: the target's von Mises and pressure against the source's matched values, over the 3D
    run's largest von Mises at k (`scale`)."""
    m = match(tgt, src)
    dvm = tgt["vm"] - m["vm"]
    dP = tgt["P"] - m["P"]
    vmt = run.von_mises(tgt["sigma"] - m["sigma"])
    dT = np.abs(tgt["T"] - m["T"])
    dx_c = max(tgt["dx"], src["dx"])
    E, nu, alpha = thermo
    sxx_t, sxx_s = core_sxx(tgt), core_sxx(src)
    floor_K = floors[(tgt["form"], round(tgt["dx"]*1e3, 3))].get(k, math.nan) if tgt["form"] == "rz" else math.nan
    rows = []
    for name, sel in regions(tgt, dx_c).items():
        V = tgt["V"][sel]
        rows.append(dict(
            target="{}{}".format(tgt["form"], round(tgt["dx"]*1e3, 3)), source="{}{}".format(
                src["form"], round(src["dx"]*1e3, 3)), k=k, t_s=t_s, region=name, n=int(sel.sum()),
            vm_peak_3d_Pa=scale,
            vm_max=float(np.abs(dvm[sel]).max())/scale, vm_rms=wrms(dvm[sel], V)/scale,
            vm_mean=float(np.sum(V*dvm[sel])/np.sum(V))/scale,
            p_max=float(np.abs(dP[sel]).max())/scale, p_rms=wrms(dP[sel], V)/scale,
            p_mean=float(np.sum(V*dP[sel])/np.sum(V))/scale,
            vmt_max=float(vmt[sel].max())/scale, vmt_rms=wrms(vmt[sel], V)/scale,
            dT_max_K=float(dT[sel].max()),
            spread_max=float(m["spread"][sel].max())/scale, spread_rms=wrms(m["spread"][sel], V)/scale,
            match_dist_max_dx=float(m["dist"][sel].max())/src["dx"], n_ties_min=int(m["n_ties"][sel].min()),
            floor_azimuth_K=floor_K, floor_stress=E*alpha*floor_K/(1.0 - nu)/scale,
            rz_load_offset=RZ_SECTION_EXCESS*load_stress/scale,
            target_vm_max=float(tgt["vm"][sel].max())/scale, target_p_absmax=float(np.abs(tgt["P"][sel]).max())/scale,
            core_sxx_target_Pa=sxx_t, core_sxx_source_Pa=sxx_s, core_sxx_ratio=sxx_t/sxx_s,
            run_target=tgt["run"], run_source=src["run"]))
    return rows, m


def run_rows(form, dx, table):
    d = run_dir(form, dx)
    with open(os.path.join(d, "summary.json")) as fh:
        summ = json.load(fh)
    g = GATE[(form, dx)]
    rows = []
    for r in run.read_intervals(os.path.join(d, "intervals.csv")):
        def num(c):
            v = r.get(c, "")
            return float(v) if v not in ("", None) else math.nan
        st = load_state(form, dx, int(r["k1"]), table)
        rows.append(dict(form=form, dx_mm=dx, processes=summ["processes"], k0=int(r["k0"]), k1=int(r["k1"]),
                         t1_s=num("t1_s"), n_particles=int(num("n_particles")), scale=num("scale"),
                         crossings=int(num("crossings")), run_crossings=int(num("run_crossings")),
                         repeats=int(num("repeats")), steps=int(num("steps")), wall_s=num("wall_s"),
                         wall_per_step_ms=1e3*num("wall_per_step_s"), ke_internal_ratio=num("ke_internal_ratio"),
                         ke_ratio=num("ke_ratio"), monitored=int(num("monitored")),
                         crossings_err_est=num("crossings_err_est"),
                         crossings_err_pressure=num("crossings_err_pressure"),
                         crossings_err_limit=num("crossings_err_limit"),
                         crossings_err_limit_basis=r.get("crossings_err_limit_basis", ""),
                         crossings_switch=int(num("crossings_switch")), monitor_wall_s=num("monitor_wall_s"),
                         drag_mean_N=num("drag_mean_N"), table_drag_N=num("table_drag_N"),
                         decel_mean_ms2=num("decel_mean_ms2"), flight_decel_ms2=num("flight_decel_ms2"),
                         max_disp_ref_dx=num("max_disp_ref_dx"), recentre_shift_dx=num("recentre_shift_dx"),
                         vm_peak_Pa=float(st["vm"].max()), p_absmax_Pa=float(np.abs(st["P"]).max()),
                         T_max_K=float(st["T"].max()), pairing_dx=pairing(form, st["x_cur"], st["x_ref"], dx*1e-3),
                         core_sxx_Pa=core_sxx(st), disturbed=int(int(r["k0"]) in DISTURBED.get((form, dx), ())),
                         gate_steps=g["steps"], gate_step_ms=g["step_ms"],
                         gate_interval_s=3600.0*g["wall_h"]/GATE_INTERVALS,
                         gate_interval_py_s=3600.0*g["wall_py_h"]/GATE_INTERVALS, run=naming.run_name_from_dir(d)))
    return rows, summ


def write_rows(path, columns, rows):
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(columns)
        for row in rows:
            row = dict({c: math.nan for c in columns}, **row)
            wr.writerow([repr(float(row[c])) if isinstance(row[c], (float, np.floating)) else row[c]
                         for c in columns])


def analyse():
    import spheral_m2_check4 as c4
    table = material.MaterialTable.load(os.path.join(REPO, PREPARED))
    thermo = c4.constants(table)
    floors = {(f, dx): body_azimuth_spread(f, dx) for f, dx, _ in FORMS if f == "rz"}
    rrows, summs = [], {}
    for form, dx, _ in FORMS:
        rr, summs[(form, dx)] = run_rows(form, dx, table)
        rrows += rr
    drag = {(r["form"], r["dx_mm"], r["k1"]): r["drag_mean_N"] for r in rrows}
    rows, profiles = [], {}
    for k in COMPARE_K:
        states = {(f, dx): load_state(f, dx, k, table) for f, dx, _ in FORMS}
        scale = float(states[("3d", 2.2)]["vm"].max())
        load_stress = VM0_CHECK6*drag[("3d", 2.2, k)]/DRAG_CHECK6
        t_s = 0.5*k
        for tkey, skey in PAIRS:
            r_, m = compare(states[tkey], states[skey], k, t_s, scale, floors, thermo, load_stress)
            rows += r_
            profiles[(tkey, skey, k)] = (states[tkey], m, scale)
    write_rows(os.path.join(REPO, CSV), COLUMNS, rows)
    write_rows(os.path.join(REPO, RUNS_CSV), RUN_COLUMNS, rrows)
    prows = probe_rows(table)
    write_rows(os.path.join(REPO, PROBE_CSV), PROBE_COLUMNS, prows)
    print("wrote {} ({} rows), {} ({} rows) and {} ({} rows)".format(CSV, len(rows), RUNS_CSV, len(rrows),
                                                                     PROBE_CSV, len(prows)))
    report(rows, rrows, summs)
    print("\nprobes (errors over the probe's largest von Mises; ref = the most crossings)")
    for r in prows:
        g = lambda c: r.get(c, math.nan)        # noqa: E731
        print("  {:22s} N {:3d}: monitor est {:.4f} pressure {:.4f} | against N {}: vM {:.4f} P {:.4f} (rms vM {:.4f}) "
              "at x {:.1f} r {:.1f} mm, T {:.0f} K, {:.1f} dx deep; plastic {} | pairing {:.4f} -> {:.4f} dx, core "
              "sxx {:.1f} MPa".format(
                  r["probe"], r["crossings"], g("est_2n"), g("pressure_2n"), r["ref_crossings"], g("err_ref_vm"),
                  g("err_ref_p"), g("err_ref_vm_rms"), g("at_x_mm"), g("at_r_mm"), g("at_T_K"), g("at_depth_dx"),
                  r["n_plastic"], r["pairing_start_dx"], r["pairing_end_dx"], r["core_sxx_Pa"]/1e6))
    plots(profiles)
    return rows, rrows


def report(rows, rrows, summs):
    print("\nRZ against 3D (fractions of the 3D run's largest von Mises at k)")
    print(("{:7s} {:7s} {:>2s} {:14s} {:>5s}" + " {:>8s}"*8 + " {:>7s} {:>8s}").format(
        "target", "source", "k", "region", "n", "vm_max", "vm_rms", "vm_mean", "p_max", "p_rms", "p_mean", "spread",
        "floor", "dT_K", "peak_MPa"))
    for r in rows:
        print("{:7s} {:7s} {:2d} {:14s} {:5d} {:8.4f} {:8.4f} {:+8.4f} {:8.4f} {:8.4f} {:+8.4f} {:8.4f} {:8.4f} "
              "{:7.2f} "
              "{:8.2f}".format(r["target"], r["source"], r["k"], r["region"], r["n"], r["vm_max"], r["vm_rms"],
                               r["vm_mean"], r["p_max"], r["p_rms"], r["p_mean"], r["spread_max"], r["floor_stress"],
                               r["dT_max_K"], r["vm_peak_3d_Pa"]/1e6))
    print("\nruns")
    for (form, dx), s in summs.items():
        rr = [r for r in rrows if (r["form"], r["dx_mm"]) == (form, dx)]
        clean = [r for r in rr if not r["disturbed"]]
        wall = math.fsum(r["wall_s"] for r in rr)
        mon = math.fsum(r["monitor_wall_s"] for r in rr if r["monitored"])
        cwall = math.fsum(r["wall_s"] + (r["monitor_wall_s"] if r["monitored"] else 0.0) for r in clean)/len(clean)
        g = GATE[(form, dx)]
        at10 = [r for r in clean if r["crossings"] == 10]
        print("  {} {} mm on {} processes: {} intervals, {} steps, {:.1f} s of intervals + {:.1f} s of monitoring; "
              "max internal KE/E {:.3e}; repeats {}; final crossings {}".format(
                  form, dx, s["processes"], len(rr), sum(r["steps"] for r in rr), wall, mon,
                  max(r["ke_internal_ratio"] for r in rr), sum(r["repeats"] for r in rr),
                  s["crossings_monitor"]["final_crossings"]))
        for n in sorted({r["crossings"] for r in clean}):
            q = [r for r in clean if r["crossings"] == n]
            print("    at {:2d} crossings ({} intervals): {:.0f}-{:.0f} steps, {:.2f} ms/step (median), {:.1f} "
                  "s/interval (median)".format(n, len(q), min(r["steps"] for r in q), max(r["steps"] for r in q),
                                    float(np.median([r["wall_per_step_ms"] for r in q])),
                                    float(np.median([r["wall_s"] for r in q]))))
        if at10:
            m10 = float(np.median([r["wall_s"] for r in at10]))
            print("    gate (spec §5, 10 crossings): {} steps, {} ms/step built-in, {:.1f} s/interval ({:.1f} s with "
                  "the Python classes); measured at 10 crossings {:.1f} s/interval = {:.2f} x the Python gate; this "
                  "run's mean with its switches and monitoring (undisturbed intervals) {:.1f} s/interval = {:.2f} x "
                  "(hot phase {:.1f} h against {} h)"
                  .format(g["steps"], g["step_ms"], 3600*g["wall_h"]/GATE_INTERVALS,
                          3600*g["wall_py_h"]/GATE_INTERVALS, m10, m10/(3600*g["wall_py_h"]/GATE_INTERVALS),
                          cwall, cwall/(3600*g["wall_py_h"]/GATE_INTERVALS), cwall*GATE_INTERVALS/3600.0,
                          g["wall_py_h"]))
        for r in rr:
            if r["monitored"]:
                print("    monitor k {:2d}->{:2d} at {:2d}: von Mises {:.4f}  pressure {:.4f}  limit {:.4f} "
                      "({}){}".format(
                    r["k0"], r["k1"], r["crossings"] if not r["crossings_switch"] else r["run_crossings"]//2,
                    r["crossings_err_est"], r["crossings_err_pressure"], r["crossings_err_limit"],
                    r["crossings_err_limit_basis"], "  -> switched to {}".format(r["run_crossings"])
                    if r["crossings_switch"] else ""))


# ------------------------------------------------------------------------------------------------- crossings probe
PROBE_OUT = os.path.join("spheral_output", "m2_probe", "task10")
PROBE_SCRIPT = os.path.join("analysis", "spheral_m2_check13_probe.py")
PROBE_CSV = os.path.join("data", "m2_check13_probe.csv")
# The probes (module docstring, "Probes"): reruns of one interval from the run's checkpoint at several crossings;
# eps > 0 with Spheral's Monaghan artificial stress (epsTensile) in the hydro, a diagnostic only
PROBES = [dict(name="monitor_rz22_k29", form="rz", dx=2.2, procs=1, k=29, Ns=(10, 20, 40, 80, 160), eps=0.0,
               dir="rz22"),
          dict(name="monitor_rz22_k39", form="rz", dx=2.2, procs=1, k=39, Ns=(10, 20, 40, 80, 160, 320), eps=0.0,
               dir="rz22"),
          dict(name="monitor_rz11_k9", form="rz", dx=1.1, procs=4, k=9, Ns=(10, 20, 40, 80), eps=0.0, dir="rz11"),
          dict(name="pairing_rz11_k4", form="rz", dx=1.1, procs=1, k=4, Ns=(10, 20, 40), eps=0.0, dir="rz11_n1"),
          dict(name="monitor_3d22_k39", form="3d", dx=2.2, procs=18, k=39, Ns=(20, 40, 80), eps=0.0, dir="3d22")]
PROBE_COLUMNS = ["probe", "form", "dx_mm", "processes", "eps_tensile", "k0", "k1", "crossings", "ref_crossings",
                 "steps", "wall_s", "ke_internal_ratio", "est_2n", "pressure_2n", "limit_2n", "err_ref_vm",
                 "err_ref_p", "err_ref_vm_rms", "at_x_mm", "at_r_mm", "at_T_K", "at_depth_dx", "vm_peak_Pa",
                 "n_plastic", "plastic_strain_max", "pairing_start_dx", "pairing_end_dx", "core_sxx_Pa"]


def probe_file(pr, n):
    suffix = "_epst{:g}".format(pr["eps"]) if pr["eps"] > 0.0 else ""
    return os.path.join(PROBE_OUT, pr["dir"], "k{}_N{}{}.npz".format(pr["k"], n, suffix))


def run_probes():
    for pr in PROBES:
        if all(os.path.isfile(os.path.join(REPO, probe_file(pr, n))) for n in pr["Ns"]):
            continue
        cmd = [LAUNCHER] + (["-n", str(pr["procs"])] if pr["procs"] > 1 else []) + [
            PROBE_SCRIPT, os.path.relpath(run_dir(pr["form"], pr["dx"]), REPO), str(pr["k"]),
            os.path.join(PROBE_OUT, pr["dir"]), ",".join(str(n) for n in pr["Ns"])] + (
            [repr(pr["eps"])] if pr["eps"] > 0.0 else [])
        print(" ".join(cmd), flush=True)
        subprocess.run(cmd, cwd=REPO, check=True, stdout=subprocess.DEVNULL)


def pairing(form, x, x_ref, dx, rmax=0.03):
    """The axial pairing amplitude: the largest |second difference| of the axial displacement (x - x_ref)[0] / dx
    between consecutive particles of each lattice line parallel to the axis (RZ: a row of fixed r; 3D: fixed
    (y, z)), over the particles with |x_ref| < rmax (30 mm: inside the core and the region under axial tension,
    away from the surface layer, whose plastic flow is not pairing).
    A smooth displacement field gives ~(dx/L)^2 u/dx; SPH's tensile pairing shows as a zigzag of order 0.01-0.1."""
    x, x_ref = np.asarray(x, dtype=float), np.asarray(x_ref, dtype=float)
    key = np.round(x_ref[:, 1:]/dx*2.0).astype(np.int64)                   # half-spacing units: exact per line
    key = key[:, 0] if key.shape[1] == 1 else key[:, 0]*100003 + key[:, 1]
    u = (x[:, 0] - x_ref[:, 0])/dx
    best = 0.0
    order = np.lexsort((x_ref[:, 0], key))
    ks, us = key[order], u[order]
    xs = np.sqrt(np.einsum("ij,ij->i", x_ref, x_ref))[order]
    starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    for lo, hi in zip(starts, np.r_[starts[1:], len(ks)]):
        m = xs[lo:hi] < rmax
        uu = us[lo:hi][m]
        if len(uu) > 2:
            best = max(best, float(np.abs(np.diff(uu, 2)).max()))
    return best


def core_sxx(st):
    """The volume-weighted mean axial stress over the core |x_ref| <= CORE (Pa)."""
    xr = np.asarray(st["x_ref"], dtype=float)
    c = np.sqrt(np.einsum("ij,ij->i", xr, xr)) <= CORE
    V = np.asarray(st["V"])[c]
    return float(np.sum(V*np.asarray(st["sigma"])[c, 0, 0])/np.sum(V))


def probe_rows(table):
    rows = []
    for pr in PROBES:
        form, dx, k, Ns = pr["form"], pr["dx"], pr["k"], pr["Ns"]
        ck0 = checkpoint.read(os.path.join(run_dir(form, dx), "checkpoints", checkpoint.filename(k)))
        p0 = pairing(form, ck0["x"], ck0["x_ref"], dx*1e-3)
        st = {}
        for n in Ns:
            z = dict(np.load(os.path.join(REPO, probe_file(pr, n))))
            s = run.stress_state(dict(z, scale=float(z["scale"])), table)
            s.update(steps=int(z["steps"]), wall=float(z["wall"]), ke=float(z["ke_internal_ratio"]),
                     eps=np.asarray(z["plastic_strain"]), pairing=pairing(form, z["x"], z["x_ref"], dx*1e-3))
            st[n] = s
        ref = st[Ns[-1]]
        top = float(ref["von_mises"].max())
        for n in Ns:
            row = dict(probe=pr["name"], form=form, dx_mm=dx, processes=pr["procs"], eps_tensile=pr["eps"], k0=k,
                       k1=k + 1, crossings=n, ref_crossings=Ns[-1], steps=st[n]["steps"], wall_s=st[n]["wall"],
                       ke_internal_ratio=st[n]["ke"], vm_peak_Pa=top,
                       n_plastic=int(np.count_nonzero(st[n]["eps"] > 0.0)),
                       plastic_strain_max=float(st[n]["eps"].max()), pairing_start_dx=p0,
                       pairing_end_dx=st[n]["pairing"], core_sxx_Pa=core_sxx(st[n]))
            if 2*n in st:
                est = run.crossings_error(st[n], st[2*n])
                pe = run.pressure_error(st[n], st[2*n])
                row.update(est_2n=est[0], pressure_2n=pe, limit_2n=max(est[0], pe))
            if n != Ns[-1]:
                assert np.array_equal(st[n]["pid"], ref["pid"])
                dvm = run.von_mises(st[n]["sigma"] - ref["sigma"])
                dP = np.abs(st[n]["P"] - ref["P"])
                i = int(np.argmax(np.maximum(dvm, dP)))
                xr = np.asarray(ref["x_ref"])[i]
                rr = math.hypot(*xr)
                row.update(err_ref_vm=float(dvm.max())/top, err_ref_p=float(dP.max())/top,
                           err_ref_vm_rms=math.sqrt(float(np.mean(dvm*dvm)))/top, at_x_mm=float(xr[0])*1e3,
                           at_r_mm=float(math.hypot(*xr[1:]))*1e3, at_T_K=float(ref["T"][i]),
                           at_depth_dx=(R - rr)/(dx*1e-3))
            rows.append(row)
    return rows


# --------------------------------------------------------------------------------------------------------- plots
def plots(profiles):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(REPO, OUT), exist_ok=True)
    colors = {"rz2.2": "#0969da", "rz1.1": "#bf3989", "3d2.2": "#1a7f37"}
    for line in ("axis", "surface"):
        fig, axes = plt.subplots(2, len(COMPARE_K), figsize=(4.0*len(COMPARE_K), 7.0), squeeze=False)
        for j, k in enumerate(COMPARE_K):
            for tkey, skey in PAIRS[:2]:
                tgt, m, scale = profiles[(tkey, skey, k)]
                dx_c = max(tgt["dx"], 2.2e-3)
                sel = regions(tgt, dx_c)[line]
                if line == "axis":
                    s = tgt["x"][sel]*1e3
                else:
                    s = np.degrees(np.arctan2(tgt["rho"][sel], tgt["x"][sel]))
                o = np.argsort(s)
                tname = "{}{}".format(tkey[0], tkey[1])
                axes[0][j].plot(s[o], tgt["vm"][sel][o]/1e6, "o-", ms=3, lw=0.8, color=colors[tname], label=tname)
                axes[1][j].plot(s[o], tgt["P"][sel][o]/1e6, "o-", ms=3, lw=0.8, color=colors[tname], label=tname)
                if tkey == ("rz", 2.2):
                    axes[0][j].plot(s[o], m["vm"][sel][o]/1e6, "s", ms=3, mfc="none", color=colors["3d2.2"],
                                    label="3d2.2 at rz2.2's rings")
                    axes[1][j].plot(s[o], m["P"][sel][o]/1e6, "s", ms=3, mfc="none", color=colors["3d2.2"],
                                    label="3d2.2 at rz2.2's rings")
            axes[0][j].set_title("k {} (t {:.1f} s)".format(k, 0.5*k), fontsize=10)
            for i in range(2):
                axes[i][j].grid(color="#d0d7de", lw=0.5)
                axes[i][j].set_xlabel("x (mm), tail -50 ... nose +50" if line == "axis" else
                                      "polar angle from the nose (deg)")
        axes[0][0].set_ylabel("von Mises (MPa)")
        axes[1][0].set_ylabel("pressure (MPa)")
        axes[0][0].legend(fontsize=7, frameon=False)
        fig.suptitle("Check 13: RZ against 3D along the {}".format("flight axis (first row)" if line == "axis" else
                                                                  "surface (outer 2.2 mm)"), fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(REPO, OUT, "profile_{}.png".format(line)), dpi=110)
        plt.close(fig)


# ------------------------------------------------------------------------------------------------------- running
def build_bodies():
    """The replay bodies of frame 0 for k 0-47 (Task 2's; already built there, so normally skipped)."""
    for form, dx, _ in FORMS:
        if os.path.isfile(os.path.join(REPO, body_dir(form, dx), "body.json")):
            continue
        subprocess.run([sys.executable, "-m", "spheral_frag", "body", "--prepared", PREPARED, "--frames",
                        "{}-{}".format(*FRAMES), "--form", form, "--dx-mm", str(dx), "--quiet"], cwd=REPO, check=True)


def launch(form, dx, procs, max_load, force):
    import spheral_m2_check4 as c4
    d = run_dir(form, dx)
    if os.path.isfile(os.path.join(d, "summary.json")):
        return "done"
    if procs > 4 and not force:
        c4.wait_for_load(0.5*max_load)
    cmd = [LAUNCHER] + (["-n", str(procs)] if procs > 1 else []) + [
        "spheral_frag/runner/run.py", "replay", "--body", body_dir(form, dx), "--prepared", PREPARED,
        "--frames", "{}-{}".format(*FRAMES), "--hourglass", repr(HOURGLASS), "--out", RUNS_OUT,
        "--max-load", repr(max_load)]
    if os.path.isdir(d):
        cmd.append("--resume")
    if force:
        cmd.append("--force")
    print(" ".join(cmd), flush=True)
    t = time.perf_counter()
    os.makedirs(os.path.join(REPO, OUT, "logs"), exist_ok=True)
    log = os.path.join(REPO, OUT, "logs", "{}_dx{}.log".format(form, dx))
    with open(log, "a") as fh:
        fh.write(" ".join(cmd) + "\n")
        fh.flush()
        r = subprocess.run(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        print("  FAILED ({}); see {}".format(r.returncode, log), flush=True)
        return "failed"
    return "ran in {:.0f} s".format(time.perf_counter() - t)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--only", choices=("run", "probe", "analyse"))
    p.add_argument("--forms", default=None, help="comma-separated form:dx to run (default all), e.g. rz:2.2")
    p.add_argument("--max-load", type=float, default=2.0)
    p.add_argument("--force", action="store_true", help="start runs even above --max-load")
    a = p.parse_args(argv)
    os.chdir(REPO)
    if a.only in (None, "run"):
        build_bodies()
        sel = None if not a.forms else {(f.split(":")[0], float(f.split(":")[1])) for f in a.forms.split(",")}
        for form, dx, procs in FORMS:
            if sel is not None and (form, dx) not in sel:
                continue
            print("{} {} mm: ".format(form, dx), flush=True)
            print(launch(form, dx, procs, a.max_load, a.force), flush=True)
    if a.only in (None, "probe"):
        run_probes()
    if a.only in (None, "analyse"):
        analyse()
    return 0


if __name__ == "__main__":
    sys.exit(main())
