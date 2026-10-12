"""The Spheral runner's command line: the mass-scaled replay (Spheral M2 plan, Task 5; spec §9.1, §11.3).

    spheral_frag/container/spheral [-n N] spheral_frag/runner/run.py replay --body DIR --prepared DIR --frames K0-K1
        [--crossings 10] [--max-crossings 80] [--damping 0.05] [--scheme sph|crksph] [--strength flight] [--rate-min 1e-6]
        [--every-check 1] [--seed 1] [--out spheral_output/runs] [--resume] [--force] [--max-load 2.0]
        [--stop-at K] [--loads on|off] [--mass-scaling on|off] [--monitor-every 10] [--monitor-limit 0.05]
        [--loads-frame K] [--area-shares patch|volume] [--hourglass 0.01]

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

**Held particles** (Spheral M2 plan, Task 7, check 5). A body whose body.npz carries `hold` (the body stage's
`--geometry bar`) has those particles held at rest by the `Hold` package (packages/hold.py), built on each attempt's
fresh node list for the local held particles (by pid) and appended after every other package: it sets their rows of
the hydro acceleration to zero, so their velocity stays zero (the start is at rest and the checkpoints keep it so)
through both stages of CheapSynchronousRK2 and their position never changes: exact, not approximate
(intervals.csv's `held_max_speed_ms` and `held_max_disp_m`, the largest |v| and |x - x_ref| over the held particles
at the interval's end, are 0.0 in every check-5 run). Spheral's `ConstantVelocityBoundary` is not used: it resets
the velocity only at the end of a CheapSynchronousRK2 step, so a held particle drifts by a dt^2/2 per step
(packages/hold.py has the measurement). A held body is fixed in space, so the boundary's
frame corrections (the net-momentum removal and the re-centring, step 5) are not applied: they would move the held
particles; `v_com_removed_ms` is 0 and `recentre_shift_dx` the measured mean axial (RZ) or norm (3D) displacement,
`recentre_applied` 0. A body that is not the frames' (body.json `geometry` other than "frame") needs `--loads off`;
the geometry is a run-name bracket (`geometry`, default frame).

**One frame's loads held** (Spheral M2 plan, Task 8, check 6). `--loads-frame K` (default -1: off) gives
SurfaceLoads frame K's prepared load table (with the lee extension at the run's brackets) at every interval, held
constant in time (`make_surface_loads` with one table), instead of the tables of the interval's own frames k and
k + 1. Everything else stays the interval's: the body, its frames' times and temperatures, the surface detection and
exposure on the current positions, and the area calibration (the shares sum to frame k's surface area -- the body's
own frame, since the body is k's, not K's). Over the first interval K0 -> K0 + 1 the table ramps in linearly
from zero (`SurfaceLoads.ramp_from_zero`), and is held from K0 + 1 on. Reason (measured on check 6, RZ 2.2 mm,
2026-10-10): applied at full strength from the stress-free state at rest, the loads set the body ringing (a dynamic
overshoot towards twice the static stress), which the damping (zeta 0.05: amplitude e^(-zeta pi N_c) = 0.21 per
interval at 10 crossings) has not settled by the end of the first interval, so the crossings monitor compared two
different ringing phases, estimated 32 % and switched the run to 20 crossings; the static state does not depend on
the path (linear elasticity), and on the flight the loads also change gradually from frame to frame. The interval's
`table_drag_N` is frame K's own integrated drag (prepare.json's `table_lee_N`, over frame K's finite-element body,
which is not the run's body when K is later than the body's frames: check 6 compares with the table's integral over
the exact sphere instead) and `flight_decel_ms2` the history's at frame K's time (a diagnostic). Needs `--loads on`;
a run-name bracket (`loads_frame`, default -1).

**Area shares** (Task 8, check 6). `--area-shares patch` (default since decision 9, 2026-10-11: each of frame k's
finite-element surface patches gives its area to the surface particle whose reference position is nearest) or
`volume` (Task 4's V^(2/3) shares with one calibration factor, which on a lattice misplace 10-13 % of the load by
inclination; packages/loads.py, "Area shares", has the reason and the measurement). A run-name bracket
(`area_shares`, default patch; the runs made with loads on before decision 9 carry `_areashares-volume`). Needs
`--loads on` when not the default.

**Schemes** (Spheral M2 plan, Task 9, check 7). `--scheme sph` (default): SolidSPH (`S.SPH` on the solid node
list; SolidSPHRZ in RZ). `--scheme crksph`: SolidCRKSPH (`S.CRKSPH`; SolidCRKSPHRZ in RZ, which Spheral 116c71f has),
linear-order reproducing kernels, with integrated density and mass-over-density volumes (spec §14.1 conditions 1
and 2): the controller's `volumeType=VolumeType.MassOverDensity` makes its VolumeUpdate package set V = m / rho
before RKCorrections, and SolidCRKSPH then evolves V by its continuity policy (dV/dt = -m/rho^2 drho/dt, the same
law while m is constant). Both with XSPH off, round H, the factories' own artificial viscosity (SPH's
MonaghanGingold, CRKSPH's LimitedMonaghanGingold) and compatible energy -- except SolidCRKSPHRZ, which cannot take a
step with compatible energy at 116c71f (it enrols its self-accelerations as a reference_wrapper that
StateDerivatives::Zero cannot visit; measured 2026-10-11) and runs without it; the replay imposes u every step, so
the energy equation is unused either way (`make_hydro`). Every package, the restart at each frame boundary
(`prime_derivatives`, no initial H iteration -- SpheralController's iterateIdealH segfaults under SolidCRKSPHRZ
before the corrections exist), the checkpoints and the monitor are the same for both. A run-name bracket (`scheme`,
default sph). Measured in check 7 (analysis/spheral_m2_check7.py).

**The crossings-error monitor** (Spheral M2 plan, decision 7, Asha 2026-10-10; Task 6b). Check 4 found 10 crossings
per interval within the plan's rule at every held state but lagging 20 crossings by 0.6 % of the peak stress right
after a ramp, so the replay estimates its own crossings error as it runs. On a monitored interval, once the interval
has been accepted at N crossings (after any slowness repeats), the same interval is run again from the same
checkpoint (the one read from disk at the interval's start; the rerun writes nothing to the run's checkpoints) at 2 N,
and the two end states are compared particle by particle (by pid) in their full stress sigma = S - P I
(`stress_state`, physical units):

    err = 2 max_i vM(sigma_N,i - sigma_2N,i) / max_i vM(sigma_N,i)      (`crossings_error`)

vM the von Mises stress sqrt(3/2 dev:dev), so the pressure part of the difference does not enter. The factor 2: the
quasi-static lag falls as 1/N, so sigma_N - sigma_inf = 2 (sigma_N - sigma_2N). Recorded with it, for context, the same
with the RMS over particles (`crossings_err_rms`) and with the 99th percentile (`crossings_err_p99`) in place of the
maximum, and the hydrostatic part the von Mises stress leaves out, 2 max |P_N - P_2N| / max vM(sigma_N)
(`crossings_err_pressure`; on check 4's ramp it is 1-5 times the von Mises estimate). **The limit applies to the
larger of the two** (controller, after Task 6b found the ramp's lag mostly hydrostatic; Asha's intent is the error,
not only its deviatoric part): err_limit = max(crossings_err_est, crossings_err_pressure) (`limit_error`), recorded
as `crossings_err_limit` with `crossings_err_limit_basis` ("von_mises" or "pressure", whichever decided; von Mises on
a tie); "err" below means err_limit.

*When it runs* (Task 10b, `MonitorRule`, `monitor_interval`; Asha's requirement: the crossings error below 5 % over
the flight). Check 13 found the first rule (the first interval and every M-th, one doubling on a switch) leaving
errors above the limit between monitored intervals (3D 30 -> 39 at 10 crossings: 17-26 % at k 39) and after a switch
(3D's last eight intervals at 20 crossings: 11-15 %). Now: (1) the first interval (K0 -> K0+1) and every
`--monitor-every` M-th after it (interval index j from K0 with j = M - 1 mod M) are monitored, M = 0 turning the
monitor off altogether; (2) once any estimate exceeds half the limit (F/2, 2.5 %) every later interval is monitored
(`monitor_watch` in intervals.csv and the checkpoints: kept for the rest of the run, also across a resume); (3) within
a monitored interval, while err > F the run takes the 2 N state (it passes the slowness test as any attempt does, else
it is repeated at 4 N), continues at 2 N, and estimates the same interval again against 4 N, doubling until the
estimate is at or below F, so every accepted monitored interval has its own estimate below the limit
(`monitor_chain`: "N:err:basis" per estimate, `crossings_switch` the number of doublings); (4) the interval after a
switch is monitored too (`monitor_force` in the checkpoints); (5) the crossings never step back down. When err <= F
at the first estimate the N run's state stays the accepted one, so a monitored run's accepted states equal an
unmonitored run's bitwise until a switch. The cap: an estimate that would need more than `--max-crossings` (2 N above
it on an interval the rule requires, or err > F where 4 N exceeds it) fails the run loudly (exit 1, `failed.json`
with the chain and reason "crossings_error"); an interval due only by the schedule whose 2 N exceeds the cap is
logged and not monitored. The per-particle von Mises stresses of the last comparison go to
`monitor/interval_<k>.npz` (scratch, not read back). Cost: one interval at 2 N per estimate -- about 2/M of the run's
time while the estimates stay below F/2, three times the run's time (N + 2 N) once the watch is on, more on the
intervals that double.

**The hourglass stiffness** (Task 10b, `--hourglass KAPPA`, packages/hourglass.py). SolidSPH's tensile pairing under
the replay's core tension (check 13: RZ 1.1 mm's lattice planes shifting alternately along the axis, x3.2 per
interval at 10 crossings, the core's axial stress lost) is removed by a stiffness kappa G V / dx^2 against the
non-affine part of each particle's displacement from its reference position (the affine least-squares fit of its
neighbours' displacements), with a viscous part on the same mode; both vanish on every affine field and conserve
momentum. kappa 0.01 (the default since Task 10b; 0 turns it off): on check 13's RZ 1.1 mm interval 4 -> 5 the
pairing falls x0.18, 0.08, 0.08 over the interval at 10, 20, 40 crossings (x3.2, 7.9, 17 without) and the core's
axial stress follows RZ 2.2 mm's to 0.3 %; checks 4-6 are unchanged or slightly closer to their references; the cost
is 1-9 % per step (analysis/spheral_m2_task10b.py, data/spheral/m2_task10b_*.csv). A run-name bracket (`hourglass`,
default 0.01; every run made before carries `_hourglass-0`).

**Run directory** `<out>/<prepared>/<run suffix>/` (`naming.run_dir(out, run name)`, Spheral M2 Task 7b): the run
name `naming.run_name(prepared, "replay", (K0, K1), form, dx_mm, brackets, seed)` split at its `__`, the prepared
input's name one level, the rest (`replay_k00000-00010_rz_dx2.20mm_seed1[_<bracket>...]`) the next, so that no path
component exceeds 255 characters; the brackets crossings, damping, scheme, strength, rate_min, the body's sampler and
azimuth, loads, mass_scaling, monitor_every, monitor_limit, the body's geometry, loads_frame, area_shares and hourglass are suffixes when they differ from
their defaults. The full name is the run's identity (meta.json's and summary.json's `run_name`;
`naming.run_name_from_dir` rebuilds it from the directory):

  meta.json                  record.make_meta (the analyse contract) + the run's settings and provenance
  checks/check_<n>.npz       record.Particles at K0 (n = 0) and every --every-check intervals (always at K1):
                             physical masses and densities (body.npz's m; rho / s), T the imposed temperature, h = u,
                             x and v Spheral's (RZ: x[:, :2] = (r, z), v likewise; v is the scaled problem's), group 0
  removed.npz                every removal so far (record.Removal); `check` is the first check without the mass
  checkpoints/checkpoint_<k>.npz   the state at frame k after its removals (runner/checkpoint.py; the scale in
                             force during the interval that ended there, the driver's counters in `extra`)
  history/interval_<k>.npz   per accepted step: t_mid, sum F, a_frame (Deceleration.history), and the scale
  intervals.csv              one row per accepted interval (INTERVAL_COLUMNS)
  monitor/interval_<k>.npz   per monitored interval: pid, von Mises of sigma_N, sigma_2N and their difference
  failed.json                only when the slowness test failed at the most crossings allowed, or the crossings
                             error exceeded --monitor-limit where the run could not double its crossings again
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
    "time_scale", "spheral_t1_s", "monitored", "crossings_err_est", "crossings_err_rms", "crossings_err_p99",
    "crossings_err_pressure", "monitor_crossings", "monitor_wall_s", "monitor_ke_internal_ratio", "crossings_switch",
    "run_crossings", "crossings_err_limit", "crossings_err_limit_basis", "n_held", "held_max_speed_ms",
    "held_max_disp_m", "n_zero_area", "hourglass_energy_ratio", "monitor_chain", "monitor_watch",
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


def crossings_error(st_n, st_2n):
    """(est, rms, p99): the crossings-error estimate of the state at N crossings against the same interval at 2 N
    (module docstring, "The crossings-error monitor"): 2 max_i vM(sigma_N,i - sigma_2N,i) / max_i vM(sigma_N,i), and
    the same with the RMS and the 99th percentile over particles in place of the maximum. `st_n`, `st_2n`: dicts with
    `pid` and `sigma` (n, 3, 3) (`stress_state`'s); the particles are matched by pid and must be the same set.
    Where sigma_N is zero everywhere the estimate is 0 if the difference is too, else inf."""
    vm_n, vm_d = von_mises_diff(st_n, st_2n)
    top = float(vm_n.max()) if len(vm_n) else 0.0
    if top <= 0.0:
        same = not np.any(vm_d)
        return (0.0, 0.0, 0.0) if same else (math.inf, math.inf, math.inf)
    return (2.0*float(vm_d.max())/top, 2.0*math.sqrt(float(np.mean(vm_d*vm_d)))/top,
            2.0*float(np.percentile(vm_d, 99.0))/top)


def pressure_error(st_n, st_2n):
    """2 max_i |P_N,i - P_2N,i| / max_i vM(sigma_N,i), P = -tr(sigma)/3: the hydrostatic part of the difference,
    which the von Mises estimate does not see; recorded for context (`crossings_err_pressure`)."""
    pa, pb = np.asarray(st_n["pid"]), np.asarray(st_2n["pid"])
    if len(pa) != len(pb) or not np.array_equal(np.sort(pa), np.sort(pb)):
        raise ValueError("the two states hold different particles")
    ob = np.argsort(pb, kind="stable")
    rb = ob[np.searchsorted(pb, pa, sorter=ob)]
    sa = np.asarray(st_n["sigma"], dtype=np.float64)
    dp = np.abs(np.trace(sa - np.asarray(st_2n["sigma"], dtype=np.float64)[rb], axis1=1, axis2=2))/3.0
    top = float(von_mises(sa).max()) if len(sa) else 0.0
    if top <= 0.0:
        return 0.0 if not np.any(dp) else math.inf
    return 2.0*float(dp.max())/top


def limit_error(est, pressure):
    """(err_limit, basis): the larger of the von Mises estimate and the pressure part, and which one it is
    ("von_mises" on a tie, or where the pressure part is missing, NaN); the monitor's limit keys on it (module
    docstring)."""
    if math.isnan(pressure) or est >= pressure:
        return float(est), "von_mises"
    return float(pressure), "pressure"


def von_mises(sigma):
    """sqrt(3/2 dev:dev) of (n, 3, 3) tensors."""
    sigma = np.asarray(sigma, dtype=np.float64)
    dev = sigma - (np.trace(sigma, axis1=1, axis2=2)/3.0)[:, None, None]*np.eye(3)[None]
    return np.sqrt(1.5*np.einsum("nij,nij->n", dev, dev))


def von_mises_diff(st_n, st_2n):
    """(vM(sigma_N), vM(sigma_N - sigma_2N)) per particle, in st_n's pid order (the pid sets must agree)."""
    pa, pb = np.asarray(st_n["pid"]), np.asarray(st_2n["pid"])
    if len(pa) != len(pb) or not np.array_equal(np.sort(pa), np.sort(pb)):
        raise ValueError("the two states hold different particles")
    ob = np.argsort(pb, kind="stable")
    rb = ob[np.searchsorted(pb, pa, sorter=ob)]
    sa = np.asarray(st_n["sigma"], dtype=np.float64)
    sb = np.asarray(st_2n["sigma"], dtype=np.float64)[rb]
    return von_mises(sa), von_mises(sa - sb)


def monitor_due(j, every):
    """Whether interval j (counted from K0, 0-based) is monitored: the first, and every `every`-th after it (j = every
    - 1 mod every); never when every is 0."""
    return every > 0 and (j == 0 or j % every == every - 1)


class MonitorRule:
    """When the crossings-error monitor runs (module docstring, "The crossings-error monitor"; Task 10b, Asha's
    5 % over the flight): the first interval and every `every`-th after it (`monitor_due`); every interval once any
    estimate has exceeded half the limit (`watch`, kept for the rest of the run); and the interval after any switch
    (`force`). Within a monitored interval the driver doubles the crossings until the estimate is at or below the
    limit (never stepping back down). `every` 0 turns the monitor off altogether."""

    def __init__(self, every, limit, watch=False, force=False):
        self.every, self.limit = int(every), float(limit)
        self.watch, self.force = bool(watch), bool(force)

    def due(self, j):
        return self.every > 0 and (monitor_due(j, self.every) or self.watch or self.force)

    def required(self):
        """Whether the monitoring is required by an earlier estimate (not only by the schedule)."""
        return self.every > 0 and (self.watch or self.force)

    def observe(self, err):
        if err > 0.5*self.limit:
            self.watch = True


class MonitorCapReached(Exception):
    """The monitor's rule needs more crossings than --max-crossings allows (`monitor_interval`)."""

    def __init__(self, message, n, last, chain):
        super().__init__(message)
        self.n, self.last, self.chain = n, last, chain


def monitor_interval(rule, n, max_crossings, estimate, switch):
    """One monitored interval accepted at n crossings (Task 10b's rule; `MonitorRule`): estimate(n) compares the
    accepted state with a rerun at 2 n and returns (est, err, basis); while err > rule.limit, switch(n, err, basis)
    accepts the 2 n state (returning the crossings it was accepted at: 2 n, or more after slowness repeats) and the
    interval is estimated again against twice that. Every estimate feeds rule.observe (the every-interval watch above
    half the limit); a switch sets rule.force (the next interval is monitored too). Raises MonitorCapReached when an
    estimate would need more than max_crossings. Returns (chain [(n, err, basis)], accepted n, number of switches)."""
    chain, switches = [], 0
    rule.force = False
    while True:
        if 2*n > max_crossings:
            raise MonitorCapReached(f"the rule needs this interval monitored at {n} against {2*n} crossings, above "
                                    f"--max-crossings {max_crossings}", n, None, chain)
        est, err, basis = estimate(n)
        chain.append((n, err, basis))
        rule.observe(err)
        if err <= rule.limit:
            return chain, n, switches
        if 4*n > max_crossings:
            raise MonitorCapReached("the crossings-error estimate (the larger of the von Mises and pressure "
                                    "estimates) exceeds --monitor-limit and the run cannot continue at twice the "
                                    "crossings and still be monitored (4 N > --max-crossings)", n,
                                    (est, err, basis), chain)
        n = switch(n, err, basis)
        switches += 1
        rule.force = True


def fail_monitor(run_dir, inp, k, k_next, n, last, chain, a, explanation):
    """failed.json and a collective Failure: the monitor's rule needs more crossings than --max-crossings allows."""
    c = comm()
    doc = {"run": inp.name, "k0": k, "k1": k_next, "reason": "crossings_error", "crossings": n,
           "monitor_crossings": 2*n, "monitor_limit": a.monitor_limit, "max_crossings": a.max_crossings,
           "chain": [{"crossings": m, "limit": e, "basis": b} for m, e, b in chain], "explanation": explanation}
    if last is not None:
        est, err, basis = last
        doc.update(crossings_err_est=est[0], crossings_err_rms=est[1], crossings_err_p99=est[2],
                   crossings_err_pressure=est[3], crossings_err_limit=err, crossings_err_limit_basis=basis)
    if c.rank == 0:
        write_json(os.path.join(run_dir, "failed.json"), doc)
    lim = f"crossings error {last[1]:.3e} ({last[2]}) > {a.monitor_limit} at {n} crossings, and " if last else ""
    raise Failure(f"interval {k}->{k_next}: {lim}{explanation} (cap --max-crossings {a.max_crossings} reached)",
                  collective=True)


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
            w.writerow([_cell(r.get(c, "")) for c in INTERVAL_COLUMNS])     # rows of an older run lack columns
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
        # --loads-frame K: frame K's table held at every interval (module docstring, "One frame's loads held")
        self.loads_frame = int(a.loads_frame)
        if a.area_shares != naming.BRACKET_DEFAULTS["area_shares"] and a.loads != "on":
            raise BadInput("--area-shares needs --loads on")
        if self.loads_frame >= 0:
            if a.loads != "on":
                raise BadInput("--loads-frame needs --loads on")
            if self.loads_frame not in have:
                raise BadInput(f"--loads-frame {self.loads_frame}: the prepared loads have no such frame")
        elif self.loads_frame != -1:
            raise BadInput("--loads-frame must be a frame index, or -1 (off)")
        with np.load(os.path.join(self.prepared_dir, "flight.npz"), allow_pickle=False) as z:
            self.flight = {k: z[k] for k in z.files}
        with open(os.path.join(self.prepared_dir, "prepare.json"), encoding="utf-8") as fh:
            prep_frames = json.load(fh).get("frames", [])
        # the frames' own integrated drag, the load tables' (with the lee extension; prepare's table_lee_N): the
        # reference of the particles' drag (M1 decision 15), the history's drag a diagnostic
        self.table_drag = {int(f["k"]): (f.get("drag") or {}).get("table_lee_N") for f in prep_frames}
        self.prep_times = {int(f["k"]): float(f["time_s"]) for f in prep_frames if "time_s" in f}
        if self.loads_frame >= 0 and self.table_drag.get(self.loads_frame) is None:
            raise BadInput(f"--loads-frame {self.loads_frame}: that frame's table carries no loads")
        try:
            material.parse_flow_law(a.strength)
        except ValueError as e:
            raise BadInput(f"--strength: {e}") from None
        if a.scheme not in SCHEMES:
            raise BadInput(f"--scheme must be one of {SCHEMES}")
        if not (a.crossings > 0 and a.max_crossings >= a.crossings):
            raise BadInput("need 0 < --crossings <= --max-crossings")
        if a.damping < 0.0 or not a.rate_min > 0.0 or a.every_check < 1:
            raise BadInput("need --damping >= 0, --rate-min > 0, --every-check >= 1")
        if a.mass_scaling == "off" and a.loads != "off":
            raise BadInput("--mass-scaling off runs each interval in N_c D / c_l seconds, not in flight time: it "
                           "needs --loads off (the load tables are in flight time)")
        if not a.hourglass >= 0.0:
            raise BadInput("need --hourglass >= 0")
        if a.monitor_every < 0 or not a.monitor_limit > 0.0:
            raise BadInput("need --monitor-every >= 0 (0: off) and --monitor-limit > 0")
        self.geometry = self.body_doc.get("geometry", "frame")
        if self.geometry != "frame" and a.loads != "off":
            raise BadInput(f"the body's geometry {self.geometry!r} is not the frames' body: it needs --loads off")
        self.hold = self.body.get("hold")
        if self.hold is not None and self.hold.shape != self.body["pid"].shape:
            raise BadInput("body.npz's hold does not match its particles")
        brackets = {"crossings": a.crossings, "damping": a.damping, "scheme": a.scheme, "strength": a.strength,
                    "rate_min": a.rate_min, "sampler": self.body_doc.get("sampler", "p1"),
                    "azimuth_deg": float(self.body_doc.get("azimuth_deg", 0.0)), "loads": a.loads,
                    "mass_scaling": a.mass_scaling, "monitor_every": a.monitor_every,
                    "monitor_limit": a.monitor_limit, "geometry": self.geometry, "loads_frame": self.loads_frame,
                    "area_shares": a.area_shares, "hourglass": a.hourglass}
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
        """Mean of frames k and k+1's table drag (N); a frame without loads takes the other's, as SurfaceLoads.
        With --loads-frame K, frame K's."""
        if self.loads_frame >= 0:
            k = k_next = self.loads_frame
        a, b = self.table_drag.get(k), self.table_drag.get(k_next)
        a = b if a is None else a
        b = a if b is None else b
        return math.nan if a is None else 0.5*(float(a) + float(b))

    def flight_decel(self, t):
        """The history's deceleration at flight time t (with --loads-frame K: at frame K's time)."""
        if self.loads_frame >= 0:
            t = self.prep_times[self.loads_frame]
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
        # the bar's held particles (check 5; module docstring, "Held particles"), by pid
        self.held_pid = np.sort(b["pid"][b["hold"] != 0].astype(np.int64)) if "hold" in b else np.zeros(0, np.int64)
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
            tables = pl.read_tables(inp.prepared_dir, [inp.loads_frame] if inp.loads_frame >= 0 else [k, k_next])
            area = pl.frame_surface_area(inp.prepared_dir, k)
            patches = pl.frame_patches(inp.prepared_dir, k) if a.area_shares == "patch" else None
            loads = pl.make_surface_loads(S, body, tables, (t0, t1), area, lee=inp.lee, area_shares=a.area_shares,
                                          patches=patches)
            if inp.loads_frame >= 0 and k == inp.k0:
                loads.ramp_from_zero()          # the held table ramps in over the first interval (module docstring)
            decel = pl.make_deceleration(S, loads, s)
            packages += [loads, decel]
        hg = None
        if a.hourglass > 0.0:
            from spheral_frag.runner.packages import hourglass as phg
            hg = phg.make_hourglass(S, body, a.hourglass, self.table)
            packages.append(hg)
        if len(self.held_pid):
            # last: it zeroes the held rows after every other package has added to the acceleration
            from spheral_frag.runner.packages import hold as ph
            packages.append(ph.make_hold(S, body, np.flatnonzero(np.isin(body.local_pids(), self.held_pid))))

        hydro = make_hydro(S, a.scheme, body)
        integ = S.CheapSynchronousRK2Integrator(body.db)
        integ.appendPhysicsPackage(hydro)
        for p in packages:
            integ.appendPhysicsPackage(p)
        integ.dtMin, integ.dtMax, integ.verbose = 1e-15, 1.0, False
        ex = ck["extra"]
        steps0, cycle0 = int(ex["steps"]), int(ex["cycle"])
        ctl = SpheralController(integ, body.kernel, statsStep=10**9, printStep=10**9, redistributeStep=None,
                                restartStep=None, vizStep=None, vizTime=None, SPH=True, iterateInitialH=False,
                                initialTime=tau0, **controller_options(S, a.scheme))
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
            hm = np.isin(state["pid"], self.held_pid)
            held = dict(n_held=int(hm.sum()),
                        held_max_speed_ms=float(np.abs(v[hm]).max()) if hm.any() else math.nan,
                        held_max_disp_m=float(np.abs(dref[hm]).max()) if hm.any() else math.nan)
            info = dict(**held,scale=s, crossings=crossings, D_m=D, D_body_m=D_body, omega1=omega1, steps=steps,
                        wall_s=wall, ke_J=ke, elastic_J=el, ke_ratio=ke_ratio(ke, el), ke_com_J=ke_com,
                        ke_internal_ratio=ke_ratio(max(ke - ke_com, 0.0), el),
                        n_particles=len(state["pid"]),
                        decel_mean_ms2=float(np.mean(-s*a_fr[:, 0])) if len(a_fr) else math.nan,
                        decel_end_ms2=float(-s*a_fr[-1, 0]) if len(a_fr) else math.nan,
                        flight_decel_ms2=inp.flight_decel(0.5*(t0 + t1)),
                        drag_mean_N=float(np.mean(-F[:, 0])) if len(F) else math.nan,
                        table_drag_N=inp.table_drag_mean(k, k_next),
                        n_surface=loads.n_surface if loads is not None else -1,
                        n_exposed=loads.n_exposed if loads is not None else -1,
                        n_zero_area=loads.n_zero_area if loads is not None else -1, max_disp_interval_dx=disp_interval,
                        max_disp_ref_dx=float(np.sqrt(np.einsum("ij,ij->i", dref, dref)).max())/self.dx,
                        momentum_axial_Ns=p_axial, momentum_norm_Ns=p_norm,
                        momentum_over_impulse=(p_norm/impulse) if impulse > 0.0 else math.nan,
                        time=float(integ.currentTime), dt=float(integ.lastDt), steps_total=int(ctl.totalSteps),
                        cycle=int(integ.currentCycle), damping_max_fraction=float(
                            packages[1].max_fraction) if a.damping > 0.0 else 0.0,
                        history=h, wall_per_step_s=wall/max(steps, 1), dt_mean_s=(tau1 - tau0)/max(steps, 1),
                        time_scale=time_scale, spheral_t1_s=float(integ.currentTime),
                        hourglass_energy_ratio=ke_ratio(hg.last_energy, el) if hg is not None else 0.0)
        del ctl, integ, hydro, packages, loads, decel, imposed, body, hg
        info = c.bcast(info if c.rank == 0 else None, root=0)
        return state, info

    # -- one interval: attempts until the slowness test passes
    def run_interval(self, k, k_next, ck, crossings, tag, first=None):
        """Attempts of interval k -> k_next from checkpoint `ck` at `crossings`, doubled until the slowness test
        passes (module docstring, step 4); `first` an attempt already run at `crossings` (the monitor's). Returns
        (state, info, tried [(crossings, internal ratio, whole ratio)]); raises Failure (failed.json) past
        --max-crossings."""
        a, inp, c = self.a, self.inp, comm()
        tried = []
        while True:
            if first is not None:
                state, info = first
                first = None
            else:
                state, info = self.attempt(k, k_next, ck, crossings, f"{tag}{len(tried)}")
            ok, internal, _ = slowness(info["ke_J"], info["ke_com_J"], info["elastic_J"], a.ke_fraction)
            tried.append((crossings, internal, info["ke_ratio"]))
            say(f"  k {k:5d}->{k_next:5d}  s {info['scale']:.4e}  N_c {crossings:3d}  steps {info['steps']:6d}  "
                f"wall {info['wall_s']:8.2f} s ({info['wall_per_step_s']*1e3:.2f} ms/step)  KE/E "
                f"{info['ke_ratio']:.3e} (internal {internal:.3e})  decel {info['decel_mean_ms2']:.4g} (flight "
                f"{info['flight_decel_ms2']:.4g}) m/s2  drag {info['drag_mean_N']:.4g} (tables "
                f"{info['table_drag_N']:.4g}) N  disp {info['max_disp_interval_dx']:.3e} dx  p/impulse "
                f"{info['momentum_over_impulse']:.2e}")
            if ok:
                return state, info, tried
            if 2*crossings > a.max_crossings:
                if c.rank == 0:
                    write_json(os.path.join(self.run_dir, "failed.json"),
                               {"run": inp.name, "k0": k, "k1": k_next, "reason": "slowness",
                                "attempts": [{"crossings": n, "ke_internal_ratio": ri, "ke_ratio": rw}
                                             for n, ri, rw in tried],
                                "keyed_on": "ke_internal_ratio", "ke_fraction": a.ke_fraction,
                                "max_crossings": a.max_crossings})
                raise Failure(f"interval {k}->{k_next}: internal KE/elastic energy {internal:.3e} >= "
                              f"{a.ke_fraction} at {crossings} crossings (the most allowed: {a.max_crossings}); "
                              f"attempts {tried}", collective=True)
            crossings *= 2

    def monitor_compare(self, k, state, info, st2, info2):
        """(est, rms, p99) of the accepted state at N against the rerun at 2 N (rank 0), with the per-particle von
        Mises stresses written to monitor/interval_<k>.npz."""
        sa = stress_state(dict(state, scale=info["scale"]), self.table)
        sb = stress_state(dict(st2, scale=info2["scale"]), self.table)
        est = crossings_error(sa, sb) + (pressure_error(sa, sb),)
        vm_n, vm_d = von_mises_diff(sa, sb)
        vm_2n, _ = von_mises_diff(sb, sa)
        order = np.argsort(sb["pid"], kind="stable")
        vm_2n = vm_2n[order[np.searchsorted(sb["pid"], sa["pid"], sorter=order)]]
        d = os.path.join(self.run_dir, "monitor")
        os.makedirs(d, exist_ok=True)
        np.savez(os.path.join(d, f"interval_{k:05d}.npz"), pid=sa["pid"], vm_n=vm_n, vm_2n=vm_2n, vm_diff=vm_d,
                 crossings=int(info["crossings"]), monitor_crossings=int(info2["crossings"]),
                 est=est[0], rms=est[1], p99=est[2], pressure=est[3])
        return est

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
        m = new["m"]
        if len(self.held_pid):
            # a held body (check 5) is fixed in space by its holds: no frame correction (it would move the held
            # particles); the net momentum (the holds' reaction) and the mean displacement are recorded only
            row["v_com_removed_ms"] = 0.0
            row["recentre_shift_dx"] = recentre(self.form, m, new["x"], new["x_ref"], self.dx, math.inf)[1] \
                if self.form == "3d" else -float(math.fsum(m*(new["x"][:, 0] - new["x_ref"][:, 0]))/math.fsum(m))/self.dx
            row["recentre_applied"] = 0
        else:
            # back into the main body's frame: the removed particles carry momentum away (and SolidSPHRZ does not
            # conserve axial momentum), so the rest is given zero net momentum (RZ: the axial component only)
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
                     next_check=next_check, crossings=info["crossings"], spheral_time=info["time"],
                     run_crossings=info["run_crossings"], monitor_watch=info.get("monitor_watch", 0),
                     monitor_force=info.get("monitor_force", 0))
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


def make_hydro(S, scheme, body):
    """The hydro package of the scheme (module docstring, "Schemes"): SolidSPH or SolidCRKSPH (SolidSPHRZ or
    SolidCRKSPHRZ in RZ, chosen by Spheral's factories from the coordinate registrar and the solid node list), both
    with compatible energy, integrated density, no XSPH, round H and the factories' own artificial viscosity."""
    if scheme == "sph":
        return S.SPH(dataBase=body.db, W=body.kernel, cfl=CFL, compatibleEnergyEvolution=True,
                     densityUpdate=S.IntegrateDensity, XSPH=False, ASPH=False)
    if scheme == "crksph":
        # SolidCRKSPHRZ with compatible energy cannot step at 116c71f: it enrolls its self-accelerations in the
        # derivatives as a reference_wrapper of a FieldList, which StateDerivatives::Zero (called by every
        # integrator step) cannot visit ("AnyVisitor ERROR: unable to process unknown data type", measured
        # 2026-10-11). The replay imposes u every step (ImposedTemperature), so the energy equation is unused:
        # RZ CRKSPH runs without compatible energy, 3D keeps it as SolidSPH does.
        return S.CRKSPH(dataBase=body.db, W=body.kernel, cfl=CFL, compatibleEnergyEvolution=body.form != "rz",
                        densityUpdate=S.IntegrateDensity, XSPH=False, ASPH=False)
    raise BadInput(f"--scheme must be one of {SCHEMES}")


def controller_options(S, scheme):
    """SpheralController's keyword arguments that depend on the scheme: under CRKSPH the volumes are mass over
    density (spec §14.1 condition 2), which the controller's VolumeUpdate package computes before RKCorrections."""
    if scheme == "crksph":
        return {"volumeType": S.VolumeType.MassOverDensity}
    return {}


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
    try:
        run_dir = naming.run_dir(os.path.abspath(a.out), inp.name)
    except ValueError as e:
        raise BadInput(str(e)) from None
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
                                              next_check=1, crossings=0, spheral_time=0.0,
                                              run_crossings=a.crossings, monitor_watch=0, monitor_force=0))
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
    # the run's crossings per interval: --crossings, or twice it after each switch by the monitor (carried in the
    # checkpoints, so a resumed run continues at them)
    base = int(ck["extra"].get("run_crossings", a.crossings))
    rule = MonitorRule(a.monitor_every, a.monitor_limit, watch=bool(int(ck["extra"].get("monitor_watch", 0))),
                       force=bool(int(ck["extra"].get("monitor_force", 0))))
    for i in range(i_start, len(inp.ks) - 1):
        k, k_next = inp.ks[i], inp.ks[i + 1]
        state, info, tried = rp.run_interval(k, k_next, ck, base, "")
        repeats = len(tried) - 1
        mon = dict(monitored=0, crossings_err_est=math.nan, crossings_err_rms=math.nan, crossings_err_p99=math.nan,
                   crossings_err_pressure=math.nan, monitor_crossings=0, monitor_wall_s=0.0,
                   monitor_ke_internal_ratio=math.nan, crossings_switch=0, crossings_err_limit=math.nan,
                   crossings_err_limit_basis="", monitor_chain="", monitor_watch=int(rule.watch))
        if rule.due(i):
            n_acc = int(info["crossings"])
            if 2*n_acc > a.max_crossings and not rule.required():
                say(f"  monitor: k {k}->{k_next} not monitored: 2 x {n_acc} crossings > --max-crossings "
                    f"{a.max_crossings}")
            else:
                box = dict(wall=0.0, state=state, info=info)

                def estimate(n):
                    st2, info2 = rp.attempt(k, k_next, ck, 2*n, f"m{2*n}")     # unique per estimate
                    est = None
                    if c.rank == 0:
                        est = rp.monitor_compare(k, box["state"], box["info"], st2, info2)
                    est = c.bcast(est, root=0)
                    _, internal2, _ = slowness(info2["ke_J"], info2["ke_com_J"], info2["elastic_J"], a.ke_fraction)
                    err, basis = limit_error(est[0], est[3])
                    box["wall"] += info2["wall_s"]
                    box["pending"] = (st2, info2)
                    mon.update(monitored=1, crossings_err_est=est[0], crossings_err_rms=est[1],
                               crossings_err_p99=est[2], crossings_err_pressure=est[3], monitor_crossings=2*n,
                               monitor_wall_s=box["wall"], monitor_ke_internal_ratio=internal2,
                               crossings_err_limit=err, crossings_err_limit_basis=basis)
                    say(f"  monitor: k {k}->{k_next}  N_c {n} against {2*n}: crossings error {est[0]:.3e} "
                        f"(rms {est[1]:.3e}, p99 {est[2]:.3e}, pressure {est[3]:.3e}; limit {a.monitor_limit})  "
                        f"wall {info2['wall_s']:.2f} s")
                    return est, err, basis

                def switch(n, err, basis):
                    say(f"  monitor: crossings error {err:.3e} ({basis}) > {a.monitor_limit}: k {k}->{k_next} takes "
                        f"the {2*n}-crossing state, the run continues at {2*n} crossings, and the interval is "
                        f"monitored again against {4*n}")
                    st, inf, tried2 = rp.run_interval(k, k_next, ck, 2*n, f"s{2*n}_", first=box["pending"])
                    box.update(state=st, info=inf, tried2=box.get("tried2", []) + [tried2])
                    return int(inf["crossings"])

                try:
                    chain, n_acc, n_switch = monitor_interval(rule, n_acc, a.max_crossings, estimate, switch)
                except MonitorCapReached as e:
                    fail_monitor(run_dir, inp, k, k_next, e.n, e.last, e.chain, a, str(e))
                state, info = box["state"], box["info"]
                for t2 in box.get("tried2", []):
                    tried += t2
                    repeats += len(t2) - 1
                base = max(base, n_acc)                     # never back down
                mon["crossings_switch"] = n_switch
                mon["monitor_chain"] = ";".join(f"{n}:{e:.6e}:{b}" for n, e, b in chain)
                mon["monitor_watch"] = int(rule.watch)
        info["run_crossings"] = base
        info["monitor_watch"] = int(rule.watch)
        info["monitor_force"] = int(rule.force)
        row = None
        if c.rank == 0:
            row = {col: info.get(col, math.nan) for col in INTERVAL_COLUMNS}
            row.update(k0=k, k1=k_next, t0_s=inp.times[inp.column(k)], t1_s=inp.times[inp.column(k_next)],
                       repeats=repeats, attempts_ke_ratio=";".join(f"{n}:{ri:.6e}/{rw:.6e}" for n, ri, rw in tried),
                       run_crossings=base, **mon)
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
            "holds": {"n_held": len(rp.held_pid), "method": "Hold package (zero acceleration on the held rows; packages/hold.py)",
                      "max_speed_ms": max((float(r["held_max_speed_ms"]) for r in rows), default=None)
                      if len(rp.held_pid) else None,
                      "max_disp_m": max((float(r["held_max_disp_m"]) for r in rows), default=None)
                      if len(rp.held_pid) else None},
            "crossings_monitor": monitor_summary(rows, a),
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


def monitor_summary(rows, a):
    """summary.json's account of the crossings-error monitor over the run's intervals.csv rows."""
    def num(r, c):
        v = r.get(c, "")
        return float(v) if v not in ("", None) else math.nan
    mon = [r for r in rows if str(r.get("monitored", "0")) == "1"]
    def lim(r):        # rows written before the limit took the larger of the two lack the column
        v = num(r, "crossings_err_limit")
        return v if not math.isnan(v) else limit_error(num(r, "crossings_err_est"), num(r, "crossings_err_pressure"))[0]

    def basis(r):
        b = r.get("crossings_err_limit_basis", "")
        return b if b else limit_error(num(r, "crossings_err_est"), num(r, "crossings_err_pressure"))[1]
    ests = [(lim(r), r) for r in mon]
    best = max(ests, key=lambda e: e[0]) if ests else None
    return {
        "monitor_every": a.monitor_every, "monitor_limit": a.monitor_limit,
        "limit_basis": "max(von_mises, pressure)",
        "intervals_monitored": [[int(r["k0"]), int(r["k1"])] for r in mon],
        "max_crossings_err_limit": best[0] if best else None,
        "max_crossings_err_limit_basis": basis(best[1]) if best else None,
        "max_at": [int(best[1]["k0"]), int(best[1]["k1"])] if best else None,
        "max_crossings_err_est": max(num(r, "crossings_err_est") for r in mon) if mon else None,
        "max_crossings_err_rms": max(num(r, "crossings_err_rms") for r in mon) if mon else None,
        "max_crossings_err_p99": max(num(r, "crossings_err_p99") for r in mon) if mon else None,
        "max_crossings_err_pressure": max(num(r, "crossings_err_pressure") for r in mon) if mon else None,
        "estimates": [{"k0": int(r["k0"]), "k1": int(r["k1"]), "crossings": int(r["crossings"]),
                       "monitor_crossings": int(r["monitor_crossings"]), "est": num(r, "crossings_err_est"),
                       "rms": num(r, "crossings_err_rms"), "p99": num(r, "crossings_err_p99"),
                       "pressure": num(r, "crossings_err_pressure"), "limit": lim(r), "basis": basis(r)}
                      for r in mon],
        "switches": [{"k0": int(r["k0"]), "k1": int(r["k1"]), "to_crossings": int(r["run_crossings"]),
                      "est": num(r, "crossings_err_est"), "limit": lim(r), "basis": basis(r),
                      "doublings": int(r.get("crossings_switch", "0") or 0), "chain": r.get("monitor_chain", "")}
                     for r in rows if int(r.get("crossings_switch", "0") or 0) > 0],
        "rule": "first interval and every monitor_every-th; every interval once an estimate exceeds half the limit; "
                "the interval after a switch; within a monitored interval the crossings double until the estimate "
                "is at or below the limit (Task 10b)",
        "watch_from": next(([int(r["k0"]), int(r["k1"])] for r in rows if str(r.get("monitor_watch", "0")) == "1"),
                           None),
        "intervals_unmonitored": sum(1 for r in rows if str(r.get("monitored", "0")) != "1"),
        "final_crossings": int(rows[-1]["run_crossings"]) if rows and rows[-1].get("run_crossings") not in
        ("", None) else a.crossings,
        "monitor_wall_s": math.fsum(num(r, "monitor_wall_s") for r in mon),
    }


def settings(a, rp):
    return {"crossings": a.crossings, "max_crossings": a.max_crossings, "damping": a.damping, "scheme": a.scheme,
            "strength": a.strength, "rate_min": a.rate_min, "every_check": a.every_check, "seed": a.seed,
            "cfl": CFL, "ke_fraction": a.ke_fraction, "surface_loads": a.loads, "mass_scaling": a.mass_scaling,
            "monitor_every": a.monitor_every, "monitor_limit": a.monitor_limit, "loads_frame": a.loads_frame,
            "area_shares": a.area_shares, "hourglass": a.hourglass,
            "hydro": {"sph": "SolidSPH, compatible energy, IntegrateDensity, no XSPH",
                      "crksph": "SolidCRKSPH (linear), IntegrateDensity, MassOverDensity volumes, no XSPH, compatible "
                                "energy " + ("off (SolidCRKSPHRZ)" if rp.form == "rz" else "on")}[a.scheme],
            "ke_test": "internal (KE - KE_com)", "c_l_cold_ms": rp.c_l, "nperh": rn.NPERH,
            "loads": {"threshold": 0.5, "exposure": "first_body", "gap": 1.5, "normal_smoothing": 2.0,
                      "lee": rp.inp.lee}}


def pin_info(pin):
    return {k: pin.get(k) for k in ("spheral_commit", "image_tag", "image_id", "image_platform")}


def build_parser():
    p = argparse.ArgumentParser(prog="run.py", description="Spheral runner (spheral_frag M2)")
    sub = p.add_subparsers(dest="mode", required=True)
    r = sub.add_parser("replay", help="the mass-scaled replay (spec §9.1)")
    r.add_argument("--body", required=True,
                   help="the body stage's directory (spheral_output/body/<prepared>/<body suffix>)")
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
    r.add_argument("--monitor-every", type=int, default=naming.BRACKET_DEFAULTS["monitor_every"],
                   help="estimate the crossings error on the first interval and every M-th after it by a rerun at "
                   "twice the crossings (decision 7; 0: off)")
    r.add_argument("--monitor-limit", type=float, default=naming.BRACKET_DEFAULTS["monitor_limit"],
                   help="above this crossings-error estimate the run continues at twice the crossings")
    r.add_argument("--loads-frame", type=int, default=naming.BRACKET_DEFAULTS["loads_frame"],
                   help="hold this frame's load table at every interval (check 6; default -1: each interval's own "
                   "frames)")
    r.add_argument("--area-shares", choices=("patch", "volume"), default=naming.BRACKET_DEFAULTS["area_shares"],
                   help="the surface particles' area shares: patch (default: the frame's patch areas to the nearest "
                   "surface particle) or volume (V^(2/3), calibrated)")
    r.add_argument("--hourglass", type=float, default=naming.BRACKET_DEFAULTS["hourglass"],
                   help="the hourglass stiffness kappa against SolidSPH's tensile pairing (packages/hourglass.py; "
                   "default 0.01; 0: off)")
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
