"""The Spheral runner's command line: the mass-scaled replay (Spheral M2 plan, Task 5; spec §9.1, §11.3).

    spheral_frag/container/spheral [-n N] spheral_frag/runner/run.py replay --body DIR --prepared DIR --frames K0-K1
        [--crossings 10] [--max-crossings 80] [--damping 0.05] [--scheme sph] [--strength flight] [--rate-min 1e-6]
        [--every-check 1] [--seed 1] [--out spheral_output/runs] [--resume] [--force] [--max-load 2.0]
        [--stop-at K] [--loads on|off] [--mass-scaling on|off]

**The interval loop.** For each interval k -> k+1 of the body's frames between K0 and K1 (spec §9.1):

1. Start from the checkpoint of frame k (`checkpoints/checkpoint_<k>.npz`; at K0 the body stage's particles, stress
   free at rho_free(T(K0)), at rest). Every interval starts from the file, read back, in the uninterrupted run as in
   a resumed one: a frame boundary is a restart (runner/nodes.py), which is what makes resuming bitwise.
2. The scale factor s = (c_l dt_k / (N_c D))^2: c_l the cold (300 K) longitudinal sound speed at s = 1
   (`damping.longitudinal_speed(table, 1)`, 6,146.9 m/s for AA7075), dt_k the interval's flight time, N_c the sound
   crossings per interval (`--crossings`), D the body's current axial extent (the particles' axial span plus one
   spacing, `damping.axial_extent`'s definition, from the checkpoint's positions). Masses and densities (and in RZ
   massRZ and massDensityRZ) are multiplied by s / s_previous, the EOS's scale set to s, the node list's density
   range to (0.2, 4) s rho_free(300 K); u is never scaled (D2), and the stored time step by sqrt(s / s_previous).
3. A fresh body, packages (ImposedTemperature, Damping with omega_1 = pi c_l,s / D, SurfaceLoads from the prepared
   tables of frames k and k+1 and frame k's surface area, Deceleration), SolidSPH (`S.SPH` on the solid node list,
   compatible energy, integrated density, as M0) and CheapSynchronousRK2 with the restored step count, cycle and last
   time step; derivatives primed; advanced to the frame time t_{k+1} (Spheral's time is the flight time).
4. The slowness test (spec §9.1): the kinetic energy KE = sum m v^2 / 2 with Spheral's (scaled) masses and velocities,
   against the stored elastic energy E = sum V (s:s / (4 G) + P^2 / (2 K)) with V = m / rho (the physical volume: both
   scale with s), s the deviatoric stress (in RZ with its hoop element S_zz), P the EOS's pressure from (rho, u) at
   the run's s, G(T) and K(T) the table's. Both are energies of the scaled problem (the stresses do not depend on s).
   The test keys on the *internal* kinetic energy, KE - KE_com, the KE relative to the centre-of-mass motion (RZ:
   the axial centre-of-mass velocity; 3D: the full one) (controller decision at the review of Task 5, 2026-10-10):
   KE_int < `KE_FRACTION` (5 %) E passes; otherwise the interval is repeated from its checkpoint with 2 N_c, up to
   `--max-crossings`, then the run fails (exit 1, `failed.json`). Where E = 0 the ratio is 0 if KE_int = 0, else
   inf. The whole KE, KE_com and both ratios (`ke_ratio`, `ke_internal_ratio`) are recorded in intervals.csv.
   Reason: spec §9.1's test is about ringing in the quasi-static body, and SolidSPHRZ does not conserve axial
   momentum by itself (Task 4), so in RZ a rigid axial drift, which carries no stress, makes up most of the whole KE
   (measured on the real flight, RZ 2.2 mm, k 0-4: 87-97 %). The drift is not a missing axis boundary:
   SpheralController inserts AxisBoundaryRZ itself (`insertAxisBoundary`, at construction, prepended to every
   package then appended -- all of the replay's are).
5. At the boundary the particles whose `remove_k` is k+1 (D1: the next frame's body no longer contains their
   reference position) are removed into the film account (`record.Removal`, reason centre crossing) with their
   physical mass (body.npz's m) and enthalpy h(T) of the body stage's held temperature (the last sampled state).
   Then, in this order: `removed.npz`, a record check if due (`checks/check_<n>.npz`), the interval's row in
   `intervals.csv`, the deceleration history (`history/interval_<k>.npz`), and last the checkpoint of frame k+1.
   Before the checkpoint the velocities are put back into the main body's frame: the mass-weighted mean velocity of
   the remaining particles (RZ: its axial component) is subtracted (`v_com_removed_ms`), because the removed
   particles carry momentum away and SolidSPHRZ does not conserve axial momentum (Task 4); the interval's own net
   momentum before this correction is recorded (`momentum_axial_Ns`, `momentum_norm_Ns`, and over the interval's
   drag impulse |mean sum F| dt, `momentum_over_impulse`). Then, in RZ, every particle's axial position is shifted
   by the same amount so that the remaining body's mass-weighted mean axial displacement from x_ref is zero
   (`recentre_shift_dx`, the shift in spacings; controller decision at the review of Task 5): SolidSPHRZ's drift
   otherwise accumulates in position (0.014 dx per interval on the real RZ 2.2 mm body). In 3D the mean displacement
   is measured and recorded the same way (its norm in spacings) and removed only above `RECENTRE_3D_DX` (1e-9 dx;
   `recentre_applied` says whether it was). 3D conserves momentum, and before any removal the mean displacement is
   round-off (1.3e-16 dx on the synthetic 10 mm sphere); a removal leaves the remaining particles with a mean
   displacement of their own (the removed nose had expanded forward; 1.6e-4 to 1.6e-3 dx there), which is then
   removed as in RZ (measured 2026-10-10, tests/test_spheral_frag_runner_replay.py).

**Check 4's switches** (Spheral M2 plan, Task 6). `--loads off` leaves out SurfaceLoads and Deceleration (a pure
thermal-stress problem; the load and drag columns are NaN, n_surface -1, and the displacement within an interval is
taken from the interval's start positions). `--mass-scaling off` (needs `--loads off`) keeps s = 1 and runs each
interval in (t1 - t0) / sqrt(s_N) = N_c D / c_l seconds of Spheral time, s_N being the scale factor the scaled run
would use (the same D from the checkpoint's positions), starting where the previous interval's Spheral time ended:
ImposedTemperature, the damping (omega_1 at s = 1) and the step all see the frame interval mapped affinely onto that
shorter one. The equivalence: with t = sqrt(s) t' the scaled momentum equation s rho dv/dt = div sigma becomes the
unscaled one in t', and the damping a = -2 zeta omega_1 v with omega_1 proportional to 1/sqrt(s) maps likewise, so
both runs solve the same dimensionless problem (N_c crossings per ramp, the same zeta) and their stresses must agree
up to round-off and the step sequence; KE and the elastic energy are equal (the velocities scale as 1/sqrt(s), the
masses as s). A disagreement would mean that something divides by an unscaled mass or reads an unscaled density.
`time_scale` (1/sqrt(s_N) unscaled, 1 scaled) and `spheral_t1_s` are recorded per interval; the record's times stay
the flight times. Both switches are run-name brackets (`loads`, `mass_scaling`, default on).

**Run directory** `<out>/<run name>/` (`naming.run_name(prepared, "replay", (K0, K1), form, dx_mm, brackets, seed)`;
the brackets crossings, damping, scheme, strength, rate_min and the body's sampler and azimuth are suffixes when they
differ from their defaults):

  meta.json                  record.make_meta (the analyse contract) + the run's settings and provenance
  checks/check_<n>.npz       record.Particles at K0 (n = 0) and every --every-check intervals (always at K1):
                             physical masses and densities (body.npz's m; rho / s), T the imposed temperature, h = u,
                             x and v Spheral's (RZ: x[:, :2] = (r, z), v likewise; v is the scaled problem's), group 0
  removed.npz                every removal so far (record.Removal); `check` is the first check without the mass
  checkpoints/checkpoint_<k>.npz   the state at frame k after its removals (runner/checkpoint.py; the scale in
                             force during the interval that ended there, the driver's counters in `extra`)
  history/interval_<k>.npz   per accepted step: t_mid, sum F, a_frame (Deceleration.history), and the scale
  intervals.csv              one row per accepted interval (INTERVAL_COLUMNS)
  failed.json                only when the slowness test failed at the most crossings allowed
  summary.json               last; its existence makes the run skip

**Stresses per particle** are in the checkpoints: `S` (deviatoric, (n, dim, dim); RZ `S_zz` the hoop element),
`rho` and `u` (scaled density, unscaled energy) give the pressure through the EOS at the checkpoint's `scale`;
`stress_state(ckpt, table)` below returns P, the full deviatoric tensor, von Mises, T and the physical masses (numpy
only, usable from drama_env).

Exit codes: 0 ok (also: run skipped, or stopped at --stop-at); 1 a model failure (Spheral error, the slowness test
failed at --max-crossings, accounts not closing); 2 bad arguments, a missing or inconsistent input, Spheral not
available, the VM's load above --max-load (more than 4 processes) without --force.

Spheral is imported inside the functions (the import-boundary test inspects this module without it)."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import sys
import time

import numpy as np

if __name__ == "__main__":
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from spheral_frag import loads as core_loads, material, naming, record  # noqa: E402
from spheral_frag.runner import checkpoint, nodes as rn  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PIN = os.path.join(REPO_ROOT, "spheral_frag", "container", "pin.json")
DEFAULT_OUT = os.path.join("spheral_output", "runs")
SUMMARY_SCHEMA = "spheral_frag.run_summary"
CFL = 0.25                     # Spheral's Courant number (M0's bench.py)
KE_FRACTION = 0.05             # spec §9.1's slowness test: internal KE below 5 % of the stored elastic energy
RECENTRE_3D_DX = 1e-9          # 3D: the body's mean displacement from x_ref is removed above this (in spacings)
MAX_CROSSINGS = 80             # plan Task 5: the most crossings per interval before the run fails
MAX_LOAD = 2.0                 # the Linux VM's 1-min load above which a run on more than 4 processes refuses
TIME_RTOL = 1e-12              # an interval ends when Spheral's time is within this of the frame time
SCHEMES = ("sph", "crksph")
STATUS_STOPPED = "stopped"

INTERVAL_COLUMNS = [
    "k0", "k1", "t0_s", "t1_s", "scale", "crossings", "repeats", "D_m", "omega1", "steps", "wall_s",
    "wall_per_step_s", "dt_mean_s", "ke_J", "elastic_J", "ke_ratio", "ke_com_J", "ke_internal_ratio", "n_particles", "n_removed",
    "removed_mass_kg", "removed_enthalpy_J", "decel_mean_ms2", "decel_end_ms2", "flight_decel_ms2", "drag_mean_N", "table_drag_N",
    "n_surface", "n_exposed", "max_disp_interval_dx", "max_disp_ref_dx", "momentum_axial_Ns", "momentum_norm_Ns",
    "momentum_over_impulse", "v_com_removed_ms", "recentre_shift_dx", "recentre_applied", "attempts_ke_ratio",
    "time_scale", "spheral_t1_s",
]


class Failure(Exception):
    """A model failure: exit 1 with the message. `collective`: raised on every rank alike (no MPI abort needed)."""

    def __init__(self, message, collective=False):
        super().__init__(message)
        self.collective = collective


class BadInput(Exception):
    """Bad arguments or inputs: exit 2 with the message."""


def comm():
    from mpi4py import MPI
    return MPI.COMM_WORLD


def say(*a):
    if comm().rank == 0:
        print(*a, flush=True)


# ---------------------------------------------------------------------------------------------------------------
# Numpy helpers (no Spheral; usable from drama_env)
# ---------------------------------------------------------------------------------------------------------------

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit(root=REPO_ROOT):
    """The repository's HEAD commit read from .git (the container has no git), or None."""
    try:
        gd = os.path.join(root, ".git")
        with open(os.path.join(gd, "HEAD")) as fh:
            head = fh.read().strip()
        if not head.startswith("ref: "):
            return head
        ref = head[5:]
        p = os.path.join(gd, ref)
        if os.path.isfile(p):
            with open(p) as fh:
                return fh.read().strip()
        with open(os.path.join(gd, "packed-refs")) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) == 2 and parts[1] == ref:
                    return parts[0]
    except OSError:
        pass
    return None


def scale_factor(c_l, dt_interval, crossings, extent):
    """s = (c_l dt / (N_c D))^2 (spec §9.1): N_c sound crossings of D in the interval's flight time dt."""
    return (float(c_l)*float(dt_interval)/(float(crossings)*float(extent)))**2


def axial_extent_arrays(x, dx):
    """D: the axial span of the positions plus one spacing (damping.axial_extent on gathered arrays)."""
    x = np.asarray(x)
    return float(x[:, 0].max() - x[:, 0].min()) + float(dx)


def deviatoric_full(S, S_zz=None):
    """(n, 3, 3) deviatoric stress: 3D as given; RZ the (z, r) block with the hoop element S_zz in [2, 2]."""
    S = np.asarray(S, dtype=np.float64)
    if S.shape[1] == 3:
        return S
    out = np.zeros((len(S), 3, 3))
    out[:, :2, :2] = S
    out[:, 2, 2] = 0.0 if S_zz is None else np.asarray(S_zz, dtype=np.float64)
    return out


def stress_state(ck, table):
    """Per-particle stress state of a read checkpoint (numpy only): dict of pid, T (K), P (Pa, the EOS's at the
    checkpoint's scale), S (n, 3, 3) deviatoric (RZ: rows/columns (z, r, theta)), sigma = S - P I, von_mises
    sqrt(3/2 S:S), m and rho (physical: divided by the scale), V = m / rho, x, x_ref."""
    s = float(ck["scale"])
    T = table.temperature(ck["u"])
    P = material.thermal_eos(table, ck["rho"], ck["u"], s)["P"]
    Sd = deviatoric_full(ck["S"], ck.get("S_zz"))
    sigma = Sd - P[:, None, None]*np.eye(3)[None]
    vm = np.sqrt(1.5*np.einsum("nij,nij->n", Sd, Sd))
    return dict(pid=ck["pid"], T=T, P=P, S=Sd, sigma=sigma, von_mises=vm, m=ck["m"]/s, rho=ck["rho"]/s,
                V=ck["m"]/ck["rho"], x=ck["x"], x_ref=ck["x_ref"])


def energies(state, table, scale):
    """(KE, elastic energy, KE about the centre of mass) of a gathered state of the scaled problem (module
    docstring, step 4). RZ: masses are the rings', velocities (v_z, v_r); only the axial centre-of-mass motion is
    a net motion."""
    m = state["m"]
    v = state["v"]
    ke = 0.5*math.fsum(m*np.einsum("ij,ij->i", v, v))
    M = math.fsum(m)
    p = np.array([math.fsum(m*v[:, c]) for c in range(v.shape[1])])
    if v.shape[1] == 2:
        p[1] = 0.0
    ke_com = 0.5*float(p @ p)/M
    T = table.temperature(state["u"])
    P = material.thermal_eos(table, state["rho"], state["u"], scale)["P"]
    G = material.thermal_shear_modulus(table, T)
    K = table.bulk_modulus(T)
    Sd = deviatoric_full(state["S"], state.get("S_zz"))
    ss = np.einsum("nij,nij->n", Sd, Sd)
    V = m/state["rho"]
    dev = np.where(G > 0.0, ss/(4.0*np.where(G > 0.0, G, 1.0)), 0.0)
    el = math.fsum(V*(dev + P*P/(2.0*K)))
    return ke, el, ke_com


def slowness(ke, ke_com, el, fraction=KE_FRACTION):
    """(passes, internal ratio, whole ratio): the test keys on the internal KE, KE - KE_com (module docstring)."""
    internal = ke_ratio(max(ke - ke_com, 0.0), el)
    return internal < fraction, internal, ke_ratio(ke, el)


def recentre(form, m, x, x_ref, dx, threshold_dx=RECENTRE_3D_DX):
    """(x shifted, shift in spacings, applied): RZ: the axial positions shifted by one amount so that the
    mass-weighted mean axial displacement from x_ref is zero (shift signed); 3D: the mean displacement vector removed
    only when its norm exceeds threshold_dx spacings (its norm returned either way)."""
    M = math.fsum(m)
    d = np.array([math.fsum(m*(x[:, j] - x_ref[:, j])) for j in range(x.shape[1])])/M
    if form == "rz":
        x = x.copy()
        x[:, 0] = x[:, 0] - d[0]
        return x, -float(d[0])/dx, True
    norm = float(np.linalg.norm(d))/dx
    if norm > threshold_dx:
        return x - d, norm, True
    return x, norm, False


def ke_ratio(ke, el):
    if el > 0.0:
        return ke/el
    return 0.0 if ke == 0.0 else math.inf


def smoothing_length(H):
    """h = det(H)^(-1/dim) per particle (the round H's 1/H_xx)."""
    H = np.asarray(H, dtype=np.float64)
    return np.linalg.det(H)**(-1.0/H.shape[1])


def to_record_xv(form, a):
    """Spheral (n, dim) -> the record's (n, 3): 3D as is; RZ (x, r) -> (r, z, 0)."""
    a = np.asarray(a, dtype=np.float64)
    out = np.zeros((len(a), 3))
    if form == "3d":
        out[:] = a
    else:
        out[:, 0] = a[:, 1]
        out[:, 1] = a[:, 0]
    return out


def read_intervals(path):
    if not os.path.isfile(path):
        return []
    with open(path, newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh)]


def write_intervals(path, rows):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(INTERVAL_COLUMNS)
        for r in rows:
            w.writerow([_cell(r[c]) for c in INTERVAL_COLUMNS])
    os.replace(tmp, path)


def _cell(v):
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    return str(v)


def write_json(path, doc):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2, default=_json_default, allow_nan=True)
        fh.write("\n")
    os.replace(tmp, path)


def _json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


# ---------------------------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------------------------

class Inputs:
    """The body stage's body, the prepared input and the run's settings, checked against each other."""

    def __init__(self, a):
        self.body_dir = os.path.abspath(a.body)
        self.prepared_dir = os.path.abspath(a.prepared)
        for d, need in ((self.body_dir, ("body.json", "body.npz", "temperature.npy")),
                        (self.prepared_dir, ("prepare.json", "material_table.npz", "loads.npz", "flight.npz"))):
            for f in need:
                if not os.path.isfile(os.path.join(d, f)):
                    raise BadInput(f"{d}: missing {f}")
        with open(os.path.join(self.body_dir, "body.json"), encoding="utf-8") as fh:
            self.body_doc = json.load(fh)
        self.prepared = os.path.basename(self.prepared_dir.rstrip(os.sep))
        if self.body_doc["prepared"]["name"] != self.prepared:
            raise BadInput(f"the body was built from {self.body_doc['prepared']['name']}, not {self.prepared}")
        tsha = sha256_file(os.path.join(self.prepared_dir, "material_table.npz"))
        if self.body_doc["prepared"]["material_table_sha256"] != tsha:
            raise BadInput("the prepared material table differs from the one the body was built with")
        self.table = material.MaterialTable.load(self.prepared_dir)
        self.form = self.body_doc["form"]
        self.dx_mm = float(self.body_doc["dx_mm"])
        self.dx = float(self.body_doc["dx_m"])
        with np.load(os.path.join(self.body_dir, "body.npz"), allow_pickle=False) as z:
            self.body = {k: z[k] for k in z.files}
        self.frames_k = [int(k) for k in self.body["frames_k"]]
        self.times = [float(t) for t in self.body["frames_time_s"]]
        k0, k1 = parse_frames(a.frames)
        if k0 not in self.frames_k or k1 not in self.frames_k or k1 <= k0:
            raise BadInput(f"--frames {a.frames}: need two of the body's frames {self.frames_k[0]}-"
                           f"{self.frames_k[-1]} with K1 > K0")
        self.k0, self.k1 = k0, k1
        self.ks = [k for k in self.frames_k if k0 <= k <= k1]
        with np.load(os.path.join(self.prepared_dir, "loads.npz"), allow_pickle=False) as z:
            have = set(int(k) for k in z["k"])
        lack = [k for k in self.ks if k not in have]
        if lack:
            raise BadInput(f"the prepared loads lack frames {lack[:5]}")
        with np.load(os.path.join(self.prepared_dir, "flight.npz"), allow_pickle=False) as z:
            self.flight = {k: z[k] for k in z.files}
        with open(os.path.join(self.prepared_dir, "prepare.json"), encoding="utf-8") as fh:
            prep_frames = json.load(fh).get("frames", [])
        # the frames' own integrated drag, the load tables' (with the lee extension; prepare's table_lee_N): the
        # reference of the particles' drag (M1 decision 15), the history's drag a diagnostic
        self.table_drag = {int(f["k"]): (f.get("drag") or {}).get("table_lee_N") for f in prep_frames}
        try:
            material.parse_flow_law(a.strength)
        except ValueError as e:
            raise BadInput(f"--strength: {e}") from None
        if a.scheme not in SCHEMES:
            raise BadInput(f"--scheme must be one of {SCHEMES}")
        if a.scheme == "crksph":
            raise BadInput("--scheme crksph is not implemented yet (Spheral M2 plan, Task 9)")
        if not (a.crossings > 0 and a.max_crossings >= a.crossings):
            raise BadInput("need 0 < --crossings <= --max-crossings")
        if a.damping < 0.0 or not a.rate_min > 0.0 or a.every_check < 1:
            raise BadInput("need --damping >= 0, --rate-min > 0, --every-check >= 1")
        if a.mass_scaling == "off" and a.loads != "off":
            raise BadInput("--mass-scaling off runs each interval in N_c D / c_l seconds, not in flight time: it "
                           "needs --loads off (the load tables are in flight time)")
        brackets = {"crossings": a.crossings, "damping": a.damping, "scheme": a.scheme, "strength": a.strength,
                    "rate_min": a.rate_min, "sampler": self.body_doc.get("sampler", "p1"),
                    "azimuth_deg": float(self.body_doc.get("azimuth_deg", 0.0)), "loads": a.loads,
                    "mass_scaling": a.mass_scaling}
        try:
            self.brackets = naming.canonical_brackets(brackets)
            self.name = naming.run_name(self.prepared, "replay", (k0, k1), self.form, self.dx_mm, self.brackets,
                                        a.seed)
        except ValueError as e:
            raise BadInput(str(e)) from None
        self.lee = {"base_fraction": self.brackets["base_pressure_pct"]/100.0,
                    "separation_deg": self.brackets["separation_deg"],
                    "shear_factor": self.brackets["lee_shear"]}

    def column(self, k):
        return self.frames_k.index(k)

    def table_drag_mean(self, k, k_next):
        """Mean of frames k and k+1's table drag (N); a frame without loads takes the other's, as SurfaceLoads."""
        a, b = self.table_drag.get(k), self.table_drag.get(k_next)
        a = b if a is None else a
        b = a if b is None else b
        return math.nan if a is None else 0.5*(float(a) + float(b))

    def flight_decel(self, t):
        return float(np.interp(t, self.flight["time_s"], self.flight["deceleration_ms2"]))


def parse_frames(text):
    parts = str(text).replace(":", "-").split("-")
    if len(parts) != 2 or not all(p.strip().isdigit() for p in parts):
        raise BadInput(f"--frames must be K0-K1, got {text!r}")
    return int(parts[0]), int(parts[1])


# ---------------------------------------------------------------------------------------------------------------
# The replay
# ---------------------------------------------------------------------------------------------------------------

class Replay:
    def __init__(self, a, inp, run_dir):
        self.a, self.inp, self.run_dir = a, inp, run_dir
        self.form, self.dx = inp.form, inp.dx
        self.S, _ = rn.spheral_module(self.form)
        from spheral_frag.runner.eos import make_eos
        from spheral_frag.runner.strength import make_strength
        from spheral_frag.runner.packages import damping as pd
        self.pd = pd
        self.table = inp.table
        self.eos = make_eos(self.S, self.table, 1.0)
        self.strength = make_strength(self.S, self.eos, a.strength, a.rate_min)
        self.kernel = self.S.TableKernel(self.S.WendlandC4Kernel(), 1000)
        self.rho_ref = float(self.table.rho_free(material.T_ZERO))
        self.c_l = float(pd.longitudinal_speed(self.table, 1.0))
        b = inp.body
        self.body_pid = b["pid"].astype(np.int64)
        self.body_order = np.argsort(self.body_pid, kind="stable")
        self.body_m = b["m"]
        self.remove_k = b["remove_k"]
        self.T_table = np.load(os.path.join(inp.body_dir, "temperature.npy"), mmap_mode="r")
        self.ck_dir = os.path.join(run_dir, "checkpoints")
        self.hist_dir = os.path.join(run_dir, "history")
        self.removed_path = os.path.join(run_dir, "removed.npz")
        self.intervals_path = os.path.join(run_dir, "intervals.csv")

    # -- lookups by pid
    def rows(self, pid):
        i = np.searchsorted(self.body_pid, pid, sorter=self.body_order)
        r = self.body_order[np.clip(i, 0, len(self.body_pid) - 1)]
        if not np.array_equal(self.body_pid[r], pid):
            raise Failure("a particle's pid is not in body.npz")
        return r

    def T_at(self, pid, k):
        return np.asarray(self.T_table[:, self.inp.column(k)], dtype=np.float64)[self.rows(pid)]

    # -- the start
    def initial_state(self):
        """The body stage's particles present at K0 (a particle the body stage removed at or before K0 is not),
        stress free and at rest, at s = 1."""
        b = self.inp.body
        keep = ~((b["remove_k"] >= 0) & (b["remove_k"] <= self.inp.k0))
        pid = b["pid"][keep].astype(np.int64)
        order = np.argsort(pid, kind="stable")
        idx = np.flatnonzero(keep)[order]
        T = np.asarray(self.T_table[:, self.inp.column(self.inp.k0)], dtype=np.float64)[idx]
        if not np.all(np.isfinite(T)):
            raise Failure("a particle present at K0 has no temperature")
        n, dim = len(idx), rn.DIM[self.form]
        x = np.ascontiguousarray(b["x_ref"][idx], dtype=np.float64)
        rho = np.asarray(self.table.rho_free(T), dtype=np.float64)
        m = np.asarray(b["m"][idx], dtype=np.float64)
        st = dict(pid=pid[order], x_ref=x, x=x.copy(), v=np.zeros((n, dim)),
                  H=rn.h_tensor_rows(dim, n, self.dx, rn.NPERH), m=m, rho=rho,
                  u=np.asarray(self.table.enthalpy(T), dtype=np.float64), S=np.zeros((n, dim, dim)),
                  plastic_strain=np.zeros(n), plastic_strain_rate=np.zeros(n))
        if self.form == "rz":
            st.update(m_rz=m/(2.0*math.pi*x[:, 1]), rho_rz=rho.copy(), S_zz=np.zeros(n), H_zz=np.zeros(n))
        return st, T

    # -- one attempt of one interval
    def attempt(self, k, k_next, ck, crossings, attempt_no):
        S, a, inp = self.S, self.a, self.inp
        from SpheralController import SpheralController
        from spheral_frag.runner.packages import loads as pl, temperature as pt
        pd = self.pd
        t0, t1 = inp.times[inp.column(k)], inp.times[inp.column(k_next)]
        arrays = checkpoint.arrays(ck)
        s_prev = float(ck["scale"])
        D = axial_extent_arrays(arrays["x"], self.dx)
        s = scale_factor(self.c_l, t1 - t0, crossings, D)
        # Spheral's time: the flight time when mass-scaled; unscaled (s = 1), the interval runs in
        # (t1 - t0) / sqrt(s) = N_c D / c_l seconds from the previous interval's Spheral end time, the same
        # dimensionless problem (module docstring, "Without mass scaling")
        if a.mass_scaling == "off":
            time_scale = 1.0/math.sqrt(s)
            s = 1.0
            tau0 = float(ck["extra"]["spheral_time"])
        else:
            time_scale = 1.0
            tau0 = t0
        tau1 = tau0 + (t1 - t0)*time_scale if time_scale != 1.0 else t1
        times = inp.times if time_scale == 1.0 else [tau0 + (t - t0)*time_scale for t in inp.times]
        ratio = s/s_prev
        if ratio != 1.0:
            for key in ("m", "rho") + (("m_rz", "rho_rz") if self.form == "rz" else ()):
                arrays[key] = arrays[key]*ratio
        self.eos.scale = s
        rho_range = (rn.RHO_RANGE_FACTORS[0]*s*self.rho_ref, rn.RHO_RANGE_FACTORS[1]*s*self.rho_ref)
        body = rn.make_body(S, self.form, arrays, self.eos, self.strength, dx=self.dx, kernel=self.kernel,
                            name=f"body_k{k:05d}_a{attempt_no}", rho_range=rho_range)
        hist = pt.TemperatureHistory(self.T_table, times, self.body_pid)
        imposed = pt.make_imposed_temperature(S, body, hist, self.table, tau0)
        D_body = pd.axial_extent(body)
        omega1 = pd.breathing_omega(self.table, s, D_body)
        packages = [imposed]
        if a.damping > 0.0:
            packages.append(pd.make_damping(S, body, a.damping, omega1))
        loads = decel = None
        if a.loads == "on":
            tables = pl.read_tables(inp.prepared_dir, [k, k_next])
            area = pl.frame_surface_area(inp.prepared_dir, k)
            loads = pl.make_surface_loads(S, body, tables, (t0, t1), area, lee=inp.lee)
            decel = pl.make_deceleration(S, loads, s)
            packages += [loads, decel]

        hydro = S.SPH(dataBase=body.db, W=body.kernel, cfl=CFL, compatibleEnergyEvolution=True,
                      densityUpdate=S.IntegrateDensity, XSPH=False, ASPH=False)
        integ = S.CheapSynchronousRK2Integrator(body.db)
        integ.appendPhysicsPackage(hydro)
        for p in packages:
            integ.appendPhysicsPackage(p)
        integ.dtMin, integ.dtMax, integ.verbose = 1e-15, 1.0, False
        ex = ck["extra"]
        steps0, cycle0 = int(ex["steps"]), int(ex["cycle"])
        ctl = SpheralController(integ, body.kernel, statsStep=10**9, printStep=10**9, redistributeStep=None,
                                restartStep=None, vizStep=None, vizTime=None, SPH=True, iterateInitialH=False,
                                initialTime=tau0)
        dt0 = float(ex["dt"])
        integ.lastDt = (0.1*CFL*self.dx/float(pd.longitudinal_speed(self.table, s))) if dt0 <= 0.0 \
            else dt0*math.sqrt(ratio)
        integ.currentCycle = cycle0
        ctl.totalSteps = steps0
        rn.prime_derivatives(body, integ)

        c = comm()
        c.Barrier()
        wall0 = time.perf_counter()
        t_end = tau1 - TIME_RTOL*max(1.0, abs(tau1))
        while integ.currentTime < t_end:
            integ.step(tau1)
            ctl.totalSteps += 1
            ctl.doPeriodicWork()
        c.Barrier()
        wall = time.perf_counter() - wall0
        steps = ctl.totalSteps - steps0
        disp_interval = loads.max_displacement() if loads is not None else None
        state = body.gather()
        info = None
        if c.rank == 0:
            ke, el, ke_com = energies(state, self.table, s)
            h = decel.history() if decel is not None else np.zeros((0, 1 + 2*body.dim))
            if disp_interval is None:                    # loads off: from the start positions, by pid
                d0 = state["x"] - arrays["x"][np.searchsorted(arrays["pid"], state["pid"])]
                disp_interval = float(np.sqrt(np.einsum("ij,ij->i", d0, d0)).max())/self.dx
            dim = body.dim
            F = h[:, 1:1 + dim] if len(h) else np.zeros((0, dim))
            a_fr = h[:, 1 + dim:1 + 2*dim] if len(h) else np.zeros((0, dim))
            m, v = state["m"], state["v"]
            pm = np.array([math.fsum(m*v[:, j]) for j in range(dim)])
            p_axial = float(pm[0])
            p_norm = abs(p_axial) if self.form == "rz" else float(np.linalg.norm(pm))
            F_mean = F.mean(axis=0) if len(F) else np.zeros(dim)
            F_scale = abs(float(F_mean[0])) if self.form == "rz" else float(np.linalg.norm(F_mean))
            impulse = F_scale*(t1 - t0)
            dref = state["x"] - state["x_ref"]
            info = dict(scale=s, crossings=crossings, D_m=D, D_body_m=D_body, omega1=omega1, steps=steps,
                        wall_s=wall, ke_J=ke, elastic_J=el, ke_ratio=ke_ratio(ke, el), ke_com_J=ke_com,
                        ke_internal_ratio=ke_ratio(max(ke - ke_com, 0.0), el),
                        n_particles=len(state["pid"]),
                        decel_mean_ms2=float(np.mean(-s*a_fr[:, 0])) if len(a_fr) else math.nan,
                        decel_end_ms2=float(-s*a_fr[-1, 0]) if len(a_fr) else math.nan,
                        flight_decel_ms2=inp.flight_decel(0.5*(t0 + t1)),
                        drag_mean_N=float(np.mean(-F[:, 0])) if len(F) else math.nan,
                        table_drag_N=inp.table_drag_mean(k, k_next),
                        n_surface=loads.n_surface if loads is not None else -1,
                        n_exposed=loads.n_exposed if loads is not None else -1, max_disp_interval_dx=disp_interval,
                        max_disp_ref_dx=float(np.sqrt(np.einsum("ij,ij->i", dref, dref)).max())/self.dx,
                        momentum_axial_Ns=p_axial, momentum_norm_Ns=p_norm,
                        momentum_over_impulse=(p_norm/impulse) if impulse > 0.0 else math.nan,
                        time=float(integ.currentTime), dt=float(integ.lastDt), steps_total=int(ctl.totalSteps),
                        cycle=int(integ.currentCycle), damping_max_fraction=float(
                            packages[1].max_fraction) if a.damping > 0.0 else 0.0,
                        history=h, wall_per_step_s=wall/max(steps, 1), dt_mean_s=(tau1 - tau0)/max(steps, 1),
                        time_scale=time_scale, spheral_t1_s=float(integ.currentTime))
        del ctl, integ, hydro, packages, loads, decel, imposed, body
        info = c.bcast(info if c.rank == 0 else None, root=0)
        return state, info

    # -- the boundary
    def boundary(self, k_next, state, info, removals, next_check, interval_index, rows, row):
        """Remove, book, write (rank 0). Returns (next_check, the new checkpoint path)."""
        inp = self.inp
        t1 = inp.times[inp.column(k_next)]
        pid = state["pid"]
        gone = self.remove_k[self.rows(pid)] == k_next
        n_gone = int(np.count_nonzero(gone))
        due = (interval_index % self.a.every_check == 0) or k_next == inp.k1
        if n_gone:
            gp = pid[gone]
            T_held = self.T_at(gp, k_next)
            if not np.all(np.isfinite(T_held)):
                raise Failure("a removed particle has no held temperature")
            m_phys = self.body_m[self.rows(gp)]
            h = np.asarray(self.table.enthalpy(T_held), dtype=np.float64)
            removals.append(record.Removal(id=gp, t=np.full(n_gone, t1), check=np.full(n_gone, next_check),
                                           mass=m_phys, h=h,
                                           reason=np.full(n_gone, record.REASON_CENTRE_CROSSING)))
            row["removed_mass_kg"] = math.fsum(m_phys)
            row["removed_enthalpy_J"] = math.fsum(m_phys*h)
        else:
            row["removed_mass_kg"] = 0.0
            row["removed_enthalpy_J"] = 0.0
        row["n_removed"] = n_gone
        keep = ~gone
        new = {key: val[keep] for key, val in state.items()}
        # back into the main body's frame: the removed particles carry momentum away (and SolidSPHRZ does not
        # conserve axial momentum), so the rest is given zero net momentum (RZ: the axial component only)
        m = new["m"]
        v_com = np.array([math.fsum(m*new["v"][:, j]) for j in range(new["v"].shape[1])])/math.fsum(m)
        if self.form == "rz":
            v_com[1] = 0.0
        new["v"] = new["v"] - v_com
        row["v_com_removed_ms"] = float(v_com[0]) if self.form == "rz" else float(np.linalg.norm(v_com))
        new["x"], shift, applied = recentre(self.form, m, new["x"], new["x_ref"], self.dx)
        row["recentre_shift_dx"] = shift
        row["recentre_applied"] = int(applied)
        record.write_removed(self.removed_path, record.Removal.concat(removals))
        if due:
            self.write_check(new, info["scale"], k_next, t1, next_check)
            next_check += 1
        rows.append(row)
        write_intervals(self.intervals_path, rows)
        os.makedirs(self.hist_dir, exist_ok=True)
        np.savez(os.path.join(self.hist_dir, f"interval_{inp.ks[interval_index - 1]:05d}.npz"),
                 history=info["history"], scale=info["scale"], k0=inp.ks[interval_index - 1], k1=k_next)
        extra = dict(time=t1, dt=info["dt"], steps=info["steps_total"], cycle=info["cycle"],
                     next_check=next_check, crossings=info["crossings"], spheral_time=info["time"])
        path = os.path.join(self.ck_dir, checkpoint.filename(k_next))
        checkpoint.write_state(path, new, form=self.form, dx=self.dx, nperh=rn.NPERH, run=inp.name, k=k_next,
                               scale=info["scale"], processes=comm().size, extra=extra)
        return next_check, path

    def write_check(self, st, scale, k, t, n):
        pid = st["pid"]
        T = self.T_at(pid, k)
        p = record.Particles(t=t, id=pid, mass=self.body_m[self.rows(pid)], x=to_record_xv(self.form, st["x"]),
                             v=to_record_xv(self.form, st["v"]), T=T, f_l=self.table.liquid_fraction(T),
                             h=st["u"], rho=st["rho"]/scale, damage=np.zeros(len(pid)),
                             h_smooth=smoothing_length(st["H"]), group=np.zeros(len(pid), dtype=np.int64),
                             frame=k, check=n)
        record.write_check(self.run_dir, p)


def read_ck(path):
    c = comm()
    ck = checkpoint.read(path) if c.rank == 0 else None
    return c.bcast(ck, root=0)


def latest_checkpoint(ck_dir, ks):
    best = None
    for k in ks:
        p = os.path.join(ck_dir, checkpoint.filename(k))
        if os.path.isfile(p):
            best = (k, p)
    return best


def replay(a):
    c = comm()
    inp = Inputs(a)
    run_dir = os.path.join(os.path.abspath(a.out), inp.name)
    summary_path = os.path.join(run_dir, "summary.json")
    if os.path.isfile(summary_path):
        say(f"replay: done already: {run_dir}")
        return 0
    load = c.bcast(os.getloadavg() if c.rank == 0 else None, root=0)
    if c.size > 4 and load[0] > a.max_load and not a.force:
        raise BadInput(f"load average {load[0]:.2f} > {a.max_load} (the Linux VM's); not starting "
                       "(--force overrides)")
    rp = Replay(a, inp, run_dir)
    ck_dir = rp.ck_dir
    exists = os.path.isdir(run_dir) and os.listdir(run_dir)
    if exists and not a.resume:
        raise BadInput(f"{run_dir} holds an unfinished run: --resume continues it (or remove it)")
    found = c.bcast(latest_checkpoint(ck_dir, inp.ks[:-1] + [inp.ks[-1]]) if c.rank == 0 else None, root=0)
    if a.resume and found is None and exists and c.rank == 0:
        shutil.rmtree(run_dir)                          # nothing to resume from: start over
    c.Barrier()

    pin = json.load(open(PIN))
    say(f"replay {inp.name}")
    say(f"  body {inp.body_dir}\n  {inp.form}, dx {inp.dx_mm} mm, frames {inp.k0}-{inp.k1} ({len(inp.ks) - 1} "
        f"intervals), {c.size} processes, load {load[0]:.2f}")
    if rp.strength.provisional:
        say("  PROVISIONAL: " + "; ".join(rp.strength.provisional))
    t_start = time.perf_counter()
    if found is None:
        state0, T0 = rp.initial_state()
        m0 = math.fsum(state0["m"])
        h0 = math.fsum(state0["m"]*state0["u"])
        if c.rank == 0:
            os.makedirs(run_dir, exist_ok=True)
            meta = record.make_meta(inp.name, inp.prepared, "replay", inp.form, inp.dx_mm, inp.brackets, a.seed,
                                    m0, h0, len(state0["pid"]), body_dir=inp.body_dir,
                                    body=inp.body_doc.get("name"), prepared_dir=inp.prepared_dir,
                                    frames=[inp.k0, inp.k1], settings=settings(a, rp), pin=pin_info(pin),
                                    git=git_commit(), provisional=list(rp.strength.provisional))
            record.write_meta(os.path.join(run_dir, "meta.json"), meta)
            record.write_removed(rp.removed_path, record.Removal.empty())
            rp.write_check(state0, 1.0, inp.k0, inp.times[inp.column(inp.k0)], 0)
            write_intervals(rp.intervals_path, [])
            checkpoint.write_state(os.path.join(ck_dir, checkpoint.filename(inp.k0)), state0, form=inp.form,
                                   dx=inp.dx, nperh=rn.NPERH, run=inp.name, k=inp.k0, scale=1.0, processes=c.size,
                                   extra=dict(time=inp.times[inp.column(inp.k0)], dt=0.0, steps=0, cycle=0,
                                              next_check=1, crossings=0, spheral_time=0.0))
        c.Barrier()
        k_start = inp.k0
    else:
        k_start = found[0]
        say(f"  resuming from {found[1]}")
    ck = read_ck(os.path.join(ck_dir, checkpoint.filename(k_start)))
    if ck["run"] != inp.name:
        raise BadInput(f"checkpoint of {ck['run']}, not {inp.name}")
    next_check = int(ck["extra"]["next_check"])
    meta = c.bcast(record.read_meta(os.path.join(run_dir, "meta.json")) if c.rank == 0 else None, root=0)
    rows, removals = [], []
    if c.rank == 0:
        t_k = inp.times[inp.column(k_start)]
        rows = [r for r in read_intervals(rp.intervals_path) if int(r["k1"]) <= k_start]
        old = record.read_removed(rp.removed_path)
        removals = [old.subset(old.t <= t_k)] if len(old) else []
        for n, p in record.list_checks(run_dir):       # checks written after the checkpoint: rewritten below
            if n >= next_check:
                os.remove(p)

    i_start = inp.ks.index(k_start)
    for i in range(i_start, len(inp.ks) - 1):
        k, k_next = inp.ks[i], inp.ks[i + 1]
        crossings = a.crossings
        tried = []
        while True:
            state, info = rp.attempt(k, k_next, ck, crossings, len(tried))
            ok, internal, _ = slowness(info["ke_J"], info["ke_com_J"], info["elastic_J"], a.ke_fraction)
            tried.append((crossings, internal, info["ke_ratio"]))
            say(f"  k {k:5d}->{k_next:5d}  s {info['scale']:.4e}  N_c {crossings:3d}  steps {info['steps']:6d}  "
                f"wall {info['wall_s']:8.2f} s ({info['wall_per_step_s']*1e3:.2f} ms/step)  KE/E "
                f"{info['ke_ratio']:.3e} (internal {internal:.3e})  decel {info['decel_mean_ms2']:.4g} (flight {info['flight_decel_ms2']:.4g}) "
                f"m/s2  drag {info['drag_mean_N']:.4g} (tables {info['table_drag_N']:.4g}) N"
                f"  disp {info['max_disp_interval_dx']:.3e} dx  p/impulse {info['momentum_over_impulse']:.2e}")
            if ok:
                break
            if 2*crossings > a.max_crossings:
                if c.rank == 0:
                    write_json(os.path.join(run_dir, "failed.json"),
                               {"run": inp.name, "k0": k, "k1": k_next,
                                "attempts": [{"crossings": n, "ke_internal_ratio": ri, "ke_ratio": rw}
                                             for n, ri, rw in tried],
                                "keyed_on": "ke_internal_ratio", "ke_fraction": a.ke_fraction,
                                "max_crossings": a.max_crossings})
                raise Failure(f"interval {k}->{k_next}: internal KE/elastic energy {internal:.3e} >= {a.ke_fraction} "
                              f"at {crossings} crossings (the most allowed: {a.max_crossings}); attempts {tried}", collective=True)
            crossings *= 2
        row = None
        if c.rank == 0:
            row = {col: info.get(col, math.nan) for col in INTERVAL_COLUMNS}
            row.update(k0=k, k1=k_next, t0_s=inp.times[inp.column(k)], t1_s=inp.times[inp.column(k_next)],
                       repeats=len(tried) - 1, attempts_ke_ratio=";".join(f"{n}:{ri:.6e}/{rw:.6e}" for n, ri, rw in tried))
            next_check, path = rp.boundary(k_next, state, info, removals, next_check, i + 1, rows, row)
        c.Barrier()
        ck = read_ck(os.path.join(ck_dir, checkpoint.filename(k_next)))
        next_check = int(ck["extra"]["next_check"])
        if a.stop_at is not None and k_next == a.stop_at and k_next != inp.k1:
            say(f"replay: stopped at frame {k_next} (--stop-at); --resume continues")
            return 0

    # summary (last)
    if c.rank == 0:
        rows = read_intervals(rp.intervals_path)
        rem = record.read_removed(rp.removed_path)
        final_mass = math.fsum(rp.body_m[rp.rows(ck["pid"])])
        residual = math.fsum([final_mass, math.fsum(rem.mass), -meta["initial_mass_kg"]])
        fl = [float(r["ke_internal_ratio"]) for r in rows]
        summary = {
            "schema": SUMMARY_SCHEMA, "schema_version": 1, "run": inp.name, "status": "ok",
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "frames": [inp.k0, inp.k1], "n_intervals": len(rows), "processes": c.size,
            "n_particles_initial": meta["n_particles"], "n_particles_final": int(len(ck["pid"])),
            "initial_mass_kg": meta["initial_mass_kg"], "removed_mass_kg": math.fsum(rem.mass),
            "n_removed": len(rem), "mass_residual_kg": residual,
            "mass_rel_residual": residual/meta["initial_mass_kg"],
            "max_ke_internal_ratio": max(fl) if fl else None,
            "max_ke_ratio": max(float(r["ke_ratio"]) for r in rows) if rows else None,
            "max_abs_recentre_shift_dx": max(abs(float(r["recentre_shift_dx"])) for r in rows) if rows else None, "repeats": sum(int(r["repeats"]) for r in rows),
            "steps": sum(int(r["steps"]) for r in rows), "wall_s_intervals": sum(float(r["wall_s"]) for r in rows),
            "wall_s_this_invocation": time.perf_counter() - t_start,
            "max_disp_interval_dx": max(float(r["max_disp_interval_dx"]) for r in rows),
            "max_abs_momentum_over_impulse": max(abs(float(r["momentum_over_impulse"])) for r in rows),
            "settings": settings(a, rp), "brackets": inp.brackets, "pin": pin_info(pin), "git": git_commit(),
            "provisional": list(rp.strength.provisional),
            "files": {"meta": "meta.json", "checks": "checks/", "removed": "removed.npz",
                      "checkpoints": "checkpoints/", "history": "history/", "intervals": "intervals.csv"},
        }
        write_json(summary_path, summary)
        say(f"replay done: {len(rows)} intervals, {summary['steps']} steps, {summary['wall_s_intervals']:.1f} s; "
            f"removed {summary['n_removed']} ({summary['removed_mass_kg']:.6g} kg); mass residual "
            f"{summary['mass_rel_residual']:.2e} m0; max internal KE/E {summary['max_ke_internal_ratio']:.3e}")
        if abs(summary["mass_rel_residual"]) > 1e-12:
            raise Failure(f"accounts do not close: {summary['mass_rel_residual']:.3e} m0")
    c.Barrier()
    return 0


def settings(a, rp):
    return {"crossings": a.crossings, "max_crossings": a.max_crossings, "damping": a.damping, "scheme": a.scheme,
            "strength": a.strength, "rate_min": a.rate_min, "every_check": a.every_check, "seed": a.seed,
            "cfl": CFL, "ke_fraction": a.ke_fraction, "surface_loads": a.loads, "mass_scaling": a.mass_scaling, "ke_test": "internal (KE - KE_com)", "c_l_cold_ms": rp.c_l, "nperh": rn.NPERH,
            "loads": {"threshold": 0.5, "exposure": "first_body", "gap": 1.5, "normal_smoothing": 2.0,
                      "lee": rp.inp.lee}}


def pin_info(pin):
    return {k: pin.get(k) for k in ("spheral_commit", "image_tag", "image_id", "image_platform")}


def build_parser():
    p = argparse.ArgumentParser(prog="run.py", description="Spheral runner (spheral_frag M2)")
    sub = p.add_subparsers(dest="mode", required=True)
    r = sub.add_parser("replay", help="the mass-scaled replay (spec §9.1)")
    r.add_argument("--body", required=True, help="the body stage's directory (spheral_output/body/<body>)")
    r.add_argument("--prepared", required=True, help="the prepared input (spheral_output/prepare/<prepared>)")
    r.add_argument("--frames", required=True, help="K0-K1: the replay's first and last frame")
    r.add_argument("--crossings", type=int, default=naming.BRACKET_DEFAULTS["crossings"])
    r.add_argument("--max-crossings", type=int, default=MAX_CROSSINGS)
    r.add_argument("--damping", type=float, default=naming.BRACKET_DEFAULTS["damping"])
    r.add_argument("--scheme", default=naming.BRACKET_DEFAULTS["scheme"])
    r.add_argument("--strength", default=naming.BRACKET_DEFAULTS["strength"])
    r.add_argument("--rate-min", type=float, default=naming.BRACKET_DEFAULTS["rate_min"])
    r.add_argument("--loads", choices=("on", "off"), default=naming.BRACKET_DEFAULTS["loads"],
                   help="the surface loads and deceleration packages (default on; check 4 runs them off)")
    r.add_argument("--mass-scaling", choices=("on", "off"), default=naming.BRACKET_DEFAULTS["mass_scaling"],
                   help="off: s = 1 and each interval runs in N_c D / c_l seconds of Spheral time (check 4's "
                   "unscaled twin; needs --loads off)")
    r.add_argument("--every-check", type=int, default=1)
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--out", default=DEFAULT_OUT)
    r.add_argument("--resume", action="store_true", help="continue an unfinished run from its last checkpoint")
    r.add_argument("--force", action="store_true", help="start even above --max-load")
    r.add_argument("--max-load", type=float, default=MAX_LOAD)
    r.add_argument("--ke-fraction", type=float, default=KE_FRACTION, help=argparse.SUPPRESS)   # tests only
    r.add_argument("--stop-at", type=int, default=None, help="stop after the checkpoint of this frame (tests, "
                   "splitting a long run)")
    return p


def main(argv):
    try:
        a = build_parser().parse_args(argv)
    except SystemExit as e:
        return 2 if e.code else 0
    try:
        import Spheral  # noqa: F401
    except ImportError as e:
        print(f"run.py: Spheral is not available: {e}", file=sys.stderr)
        return 2
    c = comm()
    try:
        return replay(a)
    except BadInput as e:
        if c.rank == 0:
            print(f"run.py: {e}", file=sys.stderr, flush=True)
        return 2
    except Exception as e:                  # a model failure (Failure, a Spheral or package error)
        import traceback
        collective = getattr(e, "collective", False)
        if c.rank == 0 or not collective:
            print(f"run.py: FAILED: {e}", file=sys.stderr, flush=True)
            if not isinstance(e, Failure):
                traceback.print_exc()
        if c.size > 1 and not collective:
            c.Abort(1)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
