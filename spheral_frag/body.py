"""`python -m spheral_frag body`: the particle body built from a prepared frame, and each particle's temperature at
every frame of the replay range (spec §6.2-6.3, §12 check 3; M2 plan, Task 2 and decisions D1, D4).

    "$PY" -m spheral_frag body --prepared spheral_output/prepare/<prepared> --frames K0-K1 --form rz|3d --dx-mm DX
          [--azimuth-deg 0] [--out spheral_output/body] [--force] [--quiet]

**The lattice, clipped by the finite-element body (D4).** 3D: the cubic lattice of spacing dx with a node at the
origin of the finite-element body's coordinates (the sphere's initial centre, as M0's `bench.lattice_generator`
puts one at the body's centre). RZ: rows at r = (j + 1/2) dx and axial positions i dx (M0's `RZGenerator` lattice)
in the half-plane through the flight axis, which is the plane z = 0, y >= 0 of the body coordinates
(`--azimuth-deg` rotates it about the axis: y = r cos(phi), z = r sin(phi)). The flight axis is the body frame's
x axis (`contract.V_HAT`; a prepared input with another v_hat is refused). A lattice point belongs to the body when
frame k0's `geometry.TetLocator` contains it: the same test as removal. The lattice is fixed by dx and the origin,
not by the frame, so bodies built at different k share their points.

**Masses and volumes** (M1 decision 3: match the mass, report the volume). V_i = c dx^3 (3D) or c 2 pi r_i dx^2 (RZ:
an RZ particle stands for its whole ring, plan fact 2), m_i = rho_free(T_i) V_i with T_i the particle's temperature
at k0 and rho_free the material table's free density; c makes sum m_i equal the finite-element solid mass
sum phi_e rho V_e of frame k0 (the film and deep accounts are not particles). The density is rho_free(T_i), the
stress-free state. The volume difference sum V_i against sum phi_e V_e is reported, not removed.

**Temperatures and removal at the reference position (D1).** Each particle keeps the position it was built at, in
the body's coordinates (`x_sample`, (n, 3); for RZ the half-plane point at the body's azimuth). For every frame
k0 ... k1 of the prepared input the particle's temperature is the P1 interpolant of the frame's nodal T there
(`TetLocator.locate`, float64), and its removal frame `remove_k` is the first frame whose body does not contain it
(`geometry.centre_crossed`'s test), -1 if none. A removed particle stays removed. `temperature.npy` holds, per
particle and frame column: the sampled value while the particle is in the body; at its removal frame the value of
the frame before (held, so that the runner's interpolation over the last interval k -> k + 1 of a particle it
removes at k + 1 stays finite, and the booked enthalpy is that of the last sampled state); NaN after it. RZ
diagnostic: the same points at azimuths +90, +180 and +270 degrees, the largest |T(phi) - T(azimuth)| over the
particles present at each frame (the frames' departure from axisymmetry); points that fall outside the body at
another azimuth are counted, not compared.

**Check 3 numbers** (spec §12), per frame, in body.json: the finite-element solid mass (sum m_i equals frame k0's
by c), the volume difference at k0, and the heat content sum m_i h(T_i(k)) of the particles present at k against
the finite element's, `fe.fe_heat_content` of the frame through the material table (the number prepare records as
`fe_heat_content_J`; the table's enthalpy is the finite element's less h_FE(300 K), bitwise). The relative
difference is NaN where the finite element's heat content is zero (a frame at a uniform 300 K); the absolute one is
recorded beside it. Measured on the MVP flight (data/spheral/m2_check3.csv, 2026-10-09): the heat-content difference
is about 1 % at worst (k = 182) at every spacing, because it comes from the mass weighting, not from sampling --
m_i = rho_free(T_i) V_i makes hot particles lighter than the finite element's constant-density mass, and a particle
in a partly consumed element (phi < 1) carries that element's whole volume.

**Outputs** in `<out>/<body name>/` (`naming.body_name`):

  body.npz         pid (int64, 0..n-1 in lattice order), x_ref ((n, 3) in 3D; (n, 2) in RZ: axial x and r),
                   x_sample (n, 3), V, m, rho, T0 (float64), remove_k (int64, -1 = never), frames_k (int64),
                   frames_time_s, dx_m, c, azimuth_deg (scalars)
  temperature.npy  (n, n_frames) float64, Fortran order: each frame's column is contiguous, so a memory-mapped
                   read of one frame touches one block (the runner maps it, `np.load(..., mmap_mode="r")`)
  body.json        last: provenance (the prepared input, its prepare.json and material table sha256), the settings,
                   c, the check-3 numbers per frame, timings, file sizes and sha256

**Resume.** A body whose body.json exists with the same signature (prepare.json and material table sha256, form,
dx, azimuth, frames, sampler) and whose files match their recorded sha256 is skipped; `--force` rebuilds. body.npz
and temperature.npy are byte-identical when rebuilt (fixed zip timestamps; np.lib.format headers).

**Samplers.** The temperature comes from a sampler `(k, frame, locator, x_sample) -> T` (`P1_SAMPLER` by default);
removal always follows the frame's body. `--synthetic radial` (Task 6, check 4) replaces the frames' temperatures by
T = 300 K + dT min(k / n_ramp, 1) (|x| / R)^2 with R = 50 mm and |x| the distance from the body coordinates' origin
(the 100 mm sphere's centre), k the absolute frame index: zero at frame 0, ramped in linearly to frame n_ramp, held
after. dT (`--radial-delta-k`, default 50 K) and n_ramp (`--radial-ramp-frames`, default 1) are part of the
sampler's name ("radial" at the defaults, else "radialdt<K>r<n>"), which the body's name carries
(`naming.body_name(..., synthetic=...)`) and the run name's `sampler` bracket repeats. Built on the MVP flight's
frames 0-3 (the intact sphere, no removal), it gives check 4's case: a one-interval ramp held for two.

Exit codes: 0 ok; 1 a model failure (an empty body, a mass that c does not close, a frame that cannot be read);
2 bad arguments or a missing prepared input or frame.

Numpy and the standard library at module level (the core's import boundary)."""
from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import sys
import time
import zipfile
from dataclasses import dataclass, field

import numpy as np

from . import __version__, contract, fe, frames, material, naming
from .geometry import TetLocator

BODY_SCHEMA = "spheral_frag.body"
BODY_SCHEMA_VERSION = 1
FORMS = naming.FORMS
RZ_AZIMUTHS_DEG = (90.0, 180.0, 270.0)       # the RZ diagnostic's other half-planes, relative to the body's
MASS_RTOL = 1e-12                            # sum m_i against the FE solid mass: exact by c up to round-off
DEFAULT_SAMPLER = "p1"


# ------------------------------------------------------------------------------------------------------- helpers
class _TableMaterial:
    """What `fe.fe_heat_content` needs, served by a MaterialTable (as prepare does with --material-table)."""

    def __init__(self, table):
        self.rho = table.rho_solid
        self.enthalpy = table.enthalpy


def fe_heat_content(frame, table) -> float:
    """J, the finite element's heat content of `frame` above 300 K through the material table (prepare's number)."""
    return fe.fe_heat_content(frame, _TableMaterial(table))


def fe_solid_mass(frame, rho) -> float:
    """kg, sum phi_e rho V_e over the active tetrahedra (frames.mass_checks' solid mass)."""
    phi = frame.tet.get("phi")
    vol = frames.tet_volumes(frame.points, frame.tets)
    return math.fsum(rho * vol * (phi if phi is not None else 1.0))


def fe_volumes(frame) -> tuple[float, float]:
    """(sum V_e, sum phi_e V_e), m^3."""
    phi = frame.tet.get("phi")
    vol = frames.tet_volumes(frame.points, frame.tets)
    return math.fsum(vol), math.fsum(vol * (phi if phi is not None else 1.0))


def rotate_about_axis(x, deg) -> np.ndarray:
    """Points (n, 3) rotated by `deg` about the x axis (y -> z)."""
    if deg == 0.0:
        return np.array(x, dtype=np.float64)
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    out = np.array(x, dtype=np.float64)
    out[:, 1] = c * x[:, 1] - s * x[:, 2]
    out[:, 2] = s * x[:, 1] + c * x[:, 2]
    return out


# ------------------------------------------------------------------------------------------------------- lattice
def lattice(form, dx, lo, hi, azimuth_deg=0.0) -> tuple[np.ndarray, np.ndarray]:
    """(x_ref, x_sample) of every lattice point covering the box [lo, hi] (body coordinates, m).

    3D: i dx on each axis (a node at the origin); x_ref = x_sample (n, 3), x fastest, then y, then z.
    RZ: axial i dx, radius (j + 1/2) dx up to the box's largest distance from the x axis; x_ref (n, 2) = (x, r),
    x fastest, then r; x_sample (n, 3) = (x, r cos(phi), r sin(phi)) at the azimuth phi."""
    dx = float(dx)
    if not (math.isfinite(dx) and dx > 0.0):
        raise ValueError("dx must be positive, got {!r}".format(dx))
    lo = np.asarray(lo, dtype=np.float64)
    hi = np.asarray(hi, dtype=np.float64)
    i0 = np.floor(lo / dx).astype(np.int64)
    i1 = np.ceil(hi / dx).astype(np.int64)
    if form == "3d":
        ax = [np.arange(i0[d], i1[d] + 1, dtype=np.int64) * dx for d in range(3)]
        z, y, x = np.meshgrid(ax[2], ax[1], ax[0], indexing="ij")
        pts = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
        return pts, pts.copy()
    if form == "rz":
        rmax = float(np.sqrt(max(lo[1] ** 2, hi[1] ** 2) + max(lo[2] ** 2, hi[2] ** 2)))
        nr = int(math.ceil(rmax / dx - 0.5)) + 1
        ax = np.arange(i0[0], i1[0] + 1, dtype=np.int64) * dx
        r = (np.arange(nr, dtype=np.float64) + 0.5) * dx
        rr, xx = np.meshgrid(r, ax, indexing="ij")
        x_ref = np.stack([xx.ravel(), rr.ravel()], axis=1)
        x_sample = np.zeros((len(x_ref), 3))
        x_sample[:, 0] = x_ref[:, 0]
        x_sample[:, 1] = x_ref[:, 1]
        return x_ref, rotate_about_axis(x_sample, float(azimuth_deg))
    raise ValueError("form must be one of {}, got {!r}".format(FORMS, form))


def reference_volumes(form, x_ref, dx, c) -> np.ndarray:
    """V_i = c dx^3 (3D) or c 2 pi r_i dx^2 (RZ, the whole ring)."""
    dx = float(dx)
    if form == "3d":
        return np.full(len(x_ref), c * dx ** 3)
    return c * 2.0 * math.pi * x_ref[:, 1] * dx ** 2


def masses(form, x_ref, dx, T0, table, fe_mass) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """(c, V, m, rho): rho = rho_free(T0), c such that sum rho V = fe_mass (module docstring)."""
    rho = np.asarray(table.rho_free(T0), dtype=np.float64)
    unit = reference_volumes(form, x_ref, dx, 1.0)
    c = float(fe_mass) / math.fsum(rho * unit)
    V = reference_volumes(form, x_ref, dx, c)
    return c, V, rho * V, rho


# ------------------------------------------------------------------------------------------------------ samplers
def p1_sampler(k, frame, locator, x, located=None) -> np.ndarray:
    """The P1 interpolant of the frame's nodal T at x (NaN outside); `located` reuses a (tet, weights) pair."""
    idx, w = located if located is not None else locator.locate(x)
    out = np.full(len(idx), np.nan)
    ok = idx >= 0
    out[ok] = np.einsum("pi,pi->p", w[ok], np.asarray(frame.node["T"], dtype=np.float64)[frame.tets[idx[ok]]])
    return out


P1_SAMPLER = p1_sampler

# Task 6's synthetic radial profile (check 4: the heated sphere against Timoshenko and Goodier)
RADIAL_T_BASE_K = 300.0      # K, the uniform start (stress free)
RADIAL_DELTA_K = 50.0        # K, the profile's rise from the centre to r = RADIAL_R_M at full strength
RADIAL_RAMP_FRAMES = 1       # frames over which the profile ramps in linearly from frame 0; held after
RADIAL_R_M = 0.05            # m, the 100 mm sphere's radius (the profile's r is |x| in body coordinates)


def radial_temperature(k, x, delta_K=RADIAL_DELTA_K, ramp_frames=RADIAL_RAMP_FRAMES, R=RADIAL_R_M,
                       T_base=RADIAL_T_BASE_K) -> np.ndarray:
    """T = T_base + delta_K min(k / ramp_frames, 1) (|x| / R)^2 at body-coordinate points x (n, 3), frame index k:
    zero at frame 0, full at frame ramp_frames, held after (deterministic from the frame index alone)."""
    f = min(max(float(k), 0.0)/float(ramp_frames), 1.0)
    x = np.asarray(x, dtype=np.float64)
    r2 = np.einsum("ij,ij->i", x, x)/(float(R)*float(R))
    return float(T_base) + (float(delta_K)*f)*r2


def radial_sampler_name(delta_K=RADIAL_DELTA_K, ramp_frames=RADIAL_RAMP_FRAMES) -> str:
    """"radial" at the defaults, else "radialdt<K>r<frames>" (integer kelvin and frames: names hold [a-z0-9])."""
    if float(delta_K) == RADIAL_DELTA_K and int(ramp_frames) == RADIAL_RAMP_FRAMES:
        return "radial"
    if not (float(delta_K).is_integer() and float(delta_K) > 0.0 and int(ramp_frames) >= 1
            and float(ramp_frames) == int(ramp_frames)):
        raise ValueError("--radial-delta-k must be a positive whole number of kelvin and --radial-ramp-frames a "
                         "positive integer, got {!r}, {!r}".format(delta_K, ramp_frames))
    return "radialdt{}r{}".format(int(delta_K), int(ramp_frames))


def make_radial_sampler(delta_K=RADIAL_DELTA_K, ramp_frames=RADIAL_RAMP_FRAMES):
    """A sampler (k, frame, locator, x) -> T of `radial_temperature`: the frame serves only removal (module
    docstring), so a point outside the frame's body is NaN as with P1."""
    def sampler(k, frame, locator, x):
        T = radial_temperature(k, x, delta_K, ramp_frames)
        T[locator.locate(x)[0] < 0] = np.nan
        return T
    sampler.params = {"delta_K": float(delta_K), "ramp_frames": int(ramp_frames), "R_m": RADIAL_R_M,
                      "T_base_K": RADIAL_T_BASE_K}
    return sampler


SAMPLERS = {"p1": p1_sampler, "radial": make_radial_sampler()}


# ------------------------------------------------------------------------------------------------------- the body
@dataclass
class Body:
    form: str
    dx: float
    azimuth_deg: float
    pid: np.ndarray
    x_ref: np.ndarray
    x_sample: np.ndarray
    V: np.ndarray
    m: np.ndarray
    rho: np.ndarray
    T0: np.ndarray
    remove_k: np.ndarray
    frames_k: np.ndarray
    frames_time_s: np.ndarray
    c: float
    temperature: np.ndarray                  # (n, n_frames), Fortran order
    check: dict = field(default_factory=dict)
    per_frame: list = field(default_factory=list)
    timing_s: dict = field(default_factory=dict)

    def arrays(self) -> dict:
        return {"pid": self.pid, "x_ref": self.x_ref, "x_sample": self.x_sample, "V": self.V, "m": self.m,
                "rho": self.rho, "T0": self.T0, "remove_k": self.remove_k, "frames_k": self.frames_k,
                "frames_time_s": self.frames_time_s, "dx_m": np.float64(self.dx), "c": np.float64(self.c),
                "azimuth_deg": np.float64(self.azimuth_deg)}


def clip(locator, x_sample) -> np.ndarray:
    """True where frame k0's body contains the lattice point (D4)."""
    return locator.locate(x_sample)[0] >= 0


def build_body(load_frame, ks, form, dx, table, azimuth_deg=0.0, sampler=P1_SAMPLER, temperature=None,
               progress=None, points=None) -> Body:
    """The body of frame ks[0], its temperatures at every frame of `ks` and the check-3 numbers.

    load_frame: k -> a compact Frame (`frames.load_prepared_frame` of the prepared input, or one built in memory).
    temperature: None (an in-memory Fortran-order array) or a writable (n, len(ks)) array to fill, made by the
    caller once n is known: a callable n -> array (the CLI passes a memory-mapped .npy).
    progress: None or a callable (i, k, seconds) for each frame.
    points: None (the lattice) or candidate points (x_ref, x_sample) to clip and sample instead (tests)."""
    if form not in FORMS:
        raise ValueError("form must be one of {}, got {!r}".format(FORMS, form))
    if form == "3d" and azimuth_deg != 0.0:
        raise ValueError("the azimuth applies to the RZ form only")
    ks = [int(k) for k in ks]
    if not ks:
        raise ValueError("no frames")
    t_all = time.perf_counter()
    timing = {}
    t = time.perf_counter()
    f0 = load_frame(ks[0])
    loc0 = TetLocator.from_frame(f0)
    timing["frame0_load_locator"] = time.perf_counter() - t

    t = time.perf_counter()
    if points is None:
        x_ref, x_sample = lattice(form, dx, f0.points.min(axis=0), f0.points.max(axis=0), azimuth_deg)
    else:
        x_ref, x_sample = (np.asarray(a, dtype=np.float64) for a in points)
    n_candidates = len(x_ref)
    located = loc0.locate(x_sample)
    keep = located[0] >= 0
    x_ref, x_sample = np.ascontiguousarray(x_ref[keep]), np.ascontiguousarray(x_sample[keep])
    located = (located[0][keep], located[1][keep])
    n = len(x_ref)
    timing["lattice_clip"] = time.perf_counter() - t
    if n == 0:
        raise ValueError("frame {}'s body contains no lattice point at dx = {} m".format(ks[0], dx))

    T0 = np.asarray(sampler(ks[0], f0, loc0, x_sample, located), dtype=np.float64) if sampler is P1_SAMPLER \
        else np.asarray(sampler(ks[0], f0, loc0, x_sample), dtype=np.float64)
    if not np.all(np.isfinite(T0)):
        raise ValueError("frame {}: a temperature at a particle inside the body is not finite".format(ks[0]))
    t = time.perf_counter()
    fe_mass0 = fe_solid_mass(f0, table.rho_solid)
    c, V, m, rho = masses(form, x_ref, dx, T0, table, fe_mass0)
    vol_tets0, vol_phi0 = fe_volumes(f0)
    timing["masses"] = time.perf_counter() - t
    sum_m, sum_V = math.fsum(m), math.fsum(V)

    T = temperature(n) if callable(temperature) else np.empty((n, len(ks)), order="F")
    remove_k = np.full(n, -1, dtype=np.int64)
    present = np.ones(n, dtype=bool)
    per_frame = []
    t_frames = time.perf_counter()
    rz_alt = [rotate_about_axis(x_sample, a) for a in RZ_AZIMUTHS_DEG] if form == "rz" else []
    for j, k in enumerate(ks):
        tf = time.perf_counter()
        if j == 0:
            fr, loc = f0, loc0
            Tk = T0
            removed_here = np.zeros(0, dtype=np.int64)
        else:
            fr = load_frame(k)
            loc = TetLocator.from_frame(fr)
            Tk = np.full(n, np.nan)
            idx = np.flatnonzero(present)
            lk = loc.locate(x_sample[idx])
            inside = lk[0] >= 0
            if sampler is P1_SAMPLER:
                Tk[idx[inside]] = p1_sampler(k, fr, loc, x_sample[idx[inside]], (lk[0][inside], lk[1][inside]))
            else:
                Tk[idx[inside]] = np.asarray(sampler(k, fr, loc, x_sample[idx[inside]]), dtype=np.float64)
            removed_here = idx[~inside]
            remove_k[removed_here] = k
            present[removed_here] = False
            Tk[removed_here] = T[removed_here, j - 1]            # held at the removal frame (module docstring)
        T[:, j] = Tk
        Th = Tk[present]
        heat_p = math.fsum(m[present] * table.enthalpy(Th)) if len(Th) else 0.0
        heat_fe = fe_heat_content(fr, table)
        row = {"k": k, "time_s": float(fr.time_s), "n_present": int(np.count_nonzero(present)),
               "n_removed": int(len(removed_here)), "fe_solid_mass_kg": fe_solid_mass(fr, table.rho_solid),
               "particle_mass_kg": math.fsum(m[present]), "heat_particles_J": heat_p, "heat_fe_J": heat_fe,
               "heat_abs_diff_J": heat_p - heat_fe, "heat_rel_diff": _rel(heat_p, heat_fe, zero=math.nan),
               "T_min_K": float(Th.min()) if len(Th) else None, "T_max_K": float(Th.max()) if len(Th) else None}
        if rz_alt:
            spread, n_out = 0.0, 0
            for xa in rz_alt:
                Ta = np.full(n, np.nan)
                Ta[present] = p1_sampler(k, fr, loc, xa[present]) if sampler is P1_SAMPLER else \
                    np.asarray(sampler(k, fr, loc, xa[present]), dtype=np.float64)
                ok = present & np.isfinite(Ta)
                n_out += int(np.count_nonzero(present & ~np.isfinite(Ta)))
                if ok.any():
                    spread = max(spread, float(np.abs(Ta[ok] - Tk[ok]).max()))
            row["azimuth_spread_K"] = spread
            row["n_outside_other_azimuths"] = n_out
        row["seconds"] = time.perf_counter() - tf
        per_frame.append(row)
        if progress is not None:
            progress(j, k, row["seconds"])
    timing["frames"] = time.perf_counter() - t_frames
    timing["total"] = time.perf_counter() - t_all

    check = {"k0": ks[0], "n_candidates": n_candidates, "n": n, "c": c, "fe_solid_mass_k0_kg": fe_mass0,
             "sum_m_kg": sum_m, "mass_rel_diff": _rel(sum_m, fe_mass0), "sum_V_m3": sum_V,
             "fe_volume_tets_k0_m3": vol_tets0, "fe_volume_phi_k0_m3": vol_phi0,
             "volume_rel_diff": _rel(sum_V, vol_phi0), "heat_rel_diff_k0": per_frame[0]["heat_rel_diff"],
             "max_abs_heat_rel_diff": max((abs(r["heat_rel_diff"]) for r in per_frame
                                           if math.isfinite(r["heat_rel_diff"])), default=math.nan),
             "n_removed": int(np.count_nonzero(remove_k >= 0))}
    if form == "rz":
        check["max_azimuth_spread_K"] = max(r["azimuth_spread_K"] for r in per_frame)
    return Body(form=form, dx=float(dx), azimuth_deg=float(azimuth_deg), pid=np.arange(n, dtype=np.int64),
                x_ref=x_ref, x_sample=x_sample, V=V, m=m, rho=rho, T0=T0, remove_k=remove_k,
                frames_k=np.array(ks, dtype=np.int64),
                frames_time_s=np.array([r["time_s"] for r in per_frame], dtype=np.float64), c=c, temperature=T,
                check=check, per_frame=per_frame, timing_s=timing)


def _rel(a, b, zero=None) -> float:
    """a / b - 1; where b = 0: 0 if a = 0 too, else `zero` (inf by default; the heat content passes NaN, since at a
    uniform 300 K both heat contents are zero up to the interpolation's round-off and no relative number exists)."""
    if b != 0.0:
        return a / b - 1.0
    if a == 0.0:
        return 0.0
    return math.inf if zero is None else zero


# ---------------------------------------------------------------------------------------------- prepared inputs
def prepared_frames(prepared_dir) -> tuple[dict, list[int]]:
    """(prepare.json, the prepared frame indices)."""
    with open(os.path.join(prepared_dir, "prepare.json"), encoding="utf-8") as fh:
        doc = json.load(fh)
    return doc, [int(f["k"]) for f in doc["frames"] if "file" in f]


def frame_loader(prepared_dir):
    def load(k):
        return frames.load_prepared_frame(os.path.join(prepared_dir, "frames", frames.prepared_frame_name(k)))
    return load


def load_body(body_dir, mmap=True) -> tuple[dict, np.ndarray]:
    """(body.npz's arrays, temperature (n, n_frames), memory-mapped read-only by default). Numpy only."""
    with np.load(os.path.join(body_dir, "body.npz"), allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    T = np.load(os.path.join(body_dir, "temperature.npy"), mmap_mode="r" if mmap else None, allow_pickle=False)
    return arrays, T


def parse_frames(text) -> tuple[int, int]:
    parts = text.replace(":", "-").split("-")
    if len(parts) not in (1, 2) or not all(p.strip().isdigit() for p in parts):
        raise ValueError("--frames must be K0-K1 (or K0:K1, or K) with non-negative integers, got {!r}".format(text))
    k0, k1 = int(parts[0]), int(parts[-1])
    if k1 < k0:
        raise ValueError("--frames {} is reversed".format(text))
    return k0, k1


# ---------------------------------------------------------------------------------------------------- the command
def _say(quiet, *a):
    if not quiet:
        print(*a, flush=True)


def body(args) -> int:
    from .prepare import git_commit, jsonable, sha256_file, write_json
    t_start = time.perf_counter()
    quiet = args.quiet
    prepared = os.path.abspath(args.prepared)
    if not os.path.isfile(os.path.join(prepared, "prepare.json")):
        print("--prepared: no prepare.json in {}".format(prepared), file=sys.stderr)
        return 2
    try:
        k0, k1 = parse_frames(args.frames)
        dx = float(args.dx_mm) * 1e-3
        if not (math.isfinite(dx) and dx > 0.0):
            raise ValueError("--dx-mm must be positive, got {!r}".format(args.dx_mm))
        azimuth = float(args.azimuth_deg)
        if not math.isfinite(azimuth):
            raise ValueError("--azimuth-deg must be finite")
        if args.form == "3d" and azimuth != 0.0:
            raise ValueError("--azimuth-deg applies to the RZ form only")
        sampler_name = getattr(args, "synthetic", None) or DEFAULT_SAMPLER
        if sampler_name == "radial":
            dk = getattr(args, "radial_delta_k", RADIAL_DELTA_K)
            rf = getattr(args, "radial_ramp_frames", RADIAL_RAMP_FRAMES)
            sampler_name = radial_sampler_name(dk, rf)
            SAMPLERS.setdefault(sampler_name, make_radial_sampler(dk, rf))
        if sampler_name not in SAMPLERS:
            raise ValueError("unknown sampler {!r}".format(sampler_name))
        prep_doc, available = prepared_frames(prepared)
    except (ValueError, OSError, KeyError) as e:
        print(e, file=sys.stderr)
        return 2
    ks = [k for k in available if k0 <= k <= k1]
    absent = [k for k in (k0, k1) if k not in available]
    if absent:
        print("--frames: the prepared input has no frame {} (frames {}-{})".format(
            absent, available[0] if available else None, available[-1] if available else None), file=sys.stderr)
        return 2
    v_hat = prep_doc.get("flight", {}).get("v_hat", list(contract.V_HAT))
    if [float(v) for v in v_hat] != [1.0, 0.0, 0.0]:
        print("the prepared input's v_hat {} is not the body frame's x axis; the body stage needs v_hat = "
              "(1, 0, 0)".format(v_hat), file=sys.stderr)
        return 2
    prep_name = os.path.basename(prepared.rstrip(os.sep))
    try:
        name = naming.body_name(prep_name, k0, k1, args.form, float(args.dx_mm), azimuth,
                                None if sampler_name == DEFAULT_SAMPLER else sampler_name)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 2
    out = os.path.join(os.path.abspath(args.out), name)
    json_path = os.path.join(out, "body.json")
    table_path = os.path.join(prepared, "material_table.npz")
    signature = {"schema_version": BODY_SCHEMA_VERSION,
                 "prepare_json_sha256": sha256_file(os.path.join(prepared, "prepare.json")),
                 "material_table_sha256": sha256_file(table_path), "form": args.form, "dx_m": dx,
                 "azimuth_deg": azimuth, "frames": ks, "sampler": sampler_name}
    files = ("body.npz", "temperature.npy")
    if os.path.isfile(json_path) and not args.force:
        try:
            with open(json_path, encoding="utf-8") as fh:
                old = json.load(fh)
            if old.get("signature") == jsonable(signature) and all(
                    os.path.isfile(os.path.join(out, f)) and sha256_file(os.path.join(out, f)) == old["files"][f]
                    ["sha256"] for f in files):
                print("built already: {}".format(out))
                return 0
        except (OSError, ValueError, KeyError):
            pass
    try:
        table = material.MaterialTable.load(table_path)
    except (OSError, ValueError) as e:
        print("material table: {}".format(e), file=sys.stderr)
        return 2
    os.makedirs(out, exist_ok=True)
    if os.path.exists(json_path):
        os.remove(json_path)                          # body.json last: a partial rebuild never looks complete
    tpath = os.path.join(out, "temperature.npy")
    ttmp = tpath + ".tmp"
    holder = {}

    def temperature(n):
        holder["T"] = np.lib.format.open_memmap(ttmp, mode="w+", dtype=np.float64, shape=(n, len(ks)),
                                                fortran_order=True, version=(1, 0))
        return holder["T"]

    def progress(i, k, s):
        if not quiet and (i % 10 == 0 or i == len(ks) - 1):
            print("  frame {:5d}  {:.2f} s".format(k, s), flush=True)

    _say(quiet, "body {} ({} frames) -> {}".format(name, len(ks), out))
    try:
        b = build_body(frame_loader(prepared), ks, args.form, dx, table, azimuth, SAMPLERS[sampler_name],
                       temperature=temperature, progress=progress)
    except (ValueError, OSError, KeyError, EOFError, zipfile.BadZipFile) as e:
        print("body failed: {}".format(e), file=sys.stderr)
        if os.path.exists(ttmp):
            del holder["T"]
            os.remove(ttmp)
        return 1
    t = time.perf_counter()
    b.temperature.flush()
    del b.temperature, holder["T"]
    os.replace(ttmp, tpath)
    frames.write_npz(os.path.join(out, "body.npz"), b.arrays())
    t_write = time.perf_counter() - t
    sizes = {f: os.path.getsize(os.path.join(out, f)) for f in files}
    passed = abs(b.check["mass_rel_diff"]) <= MASS_RTOL
    doc = {
        "schema": BODY_SCHEMA, "schema_version": BODY_SCHEMA_VERSION, "name": name,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "spheral_frag_version": __version__, "git": git_commit(),
        "command": " ".join(["python -m spheral_frag"] + list(getattr(args, "argv", None) or ["body"])),
        "signature": signature,
        "prepared": {"dir": prepared, "name": prep_name, "prepare_json_sha256": signature["prepare_json_sha256"],
                     "material_table_sha256": signature["material_table_sha256"],
                     "fe_run": prep_doc.get("fe_run", {}).get("name"), "material": table.name},
        "form": args.form, "dx_m": dx, "dx_mm": float(args.dx_mm), "azimuth_deg": azimuth, "sampler": sampler_name,
        "sampler_params": getattr(SAMPLERS[sampler_name], "params", None),
        "origin_m": [0.0, 0.0, 0.0], "axis": [1.0, 0.0, 0.0],
        "frames": {"k0": k0, "k1": k1, "n": len(ks), "k": ks},
        "n_particles": int(len(b.pid)), "c": b.c,
        "check3": b.check, "per_frame": b.per_frame, "passed": passed,
        "temperature": {"file": "temperature.npy", "shape": [int(len(b.pid)), len(ks)], "dtype": "float64",
                        "order": "F", "after_removal": "held at the removal frame, NaN after"},
        "files": {f: {"size_bytes": sizes[f], "sha256": sha256_file(os.path.join(out, f))} for f in files},
        "timing_s": {**b.timing_s, "write": t_write, "command_total": time.perf_counter() - t_start},
    }
    write_json(json_path, doc)
    c3 = b.check
    _say(quiet, "  {} particles, c {:.6f}, mass {:+.1e}, volume {:+.3e}, heat at k0 {:+.3e} (max {:.3e}){}; "
               "{:.1f} s, {:.1f} MB".format(c3["n"], b.c, c3["mass_rel_diff"], c3["volume_rel_diff"],
                                            c3["heat_rel_diff_k0"], c3["max_abs_heat_rel_diff"],
                                            ", azimuth spread {:.3f} K".format(c3["max_azimuth_spread_K"])
                                            if args.form == "rz" else "", doc["timing_s"]["command_total"],
                                            sum(sizes.values()) / 1e6))
    if not passed:
        print("mass check failed: sum m_i / FE mass - 1 = {:.3e} (threshold {:.0e})".format(c3["mass_rel_diff"],
                                                                                            MASS_RTOL),
              file=sys.stderr)
        return 1
    return 0


DEFAULT_OUT = os.path.join(fe.REPO_ROOT, "spheral_output", "body")


def add_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument("--prepared", required=True, help="prepared input directory (spheral_output/prepare/<prepared>)")
    p.add_argument("--frames", required=True, help="K0-K1 (or K0:K1, or K), inclusive: the body is built at K0 and "
                   "its temperatures sampled at every prepared frame up to K1")
    p.add_argument("--form", required=True, choices=FORMS, help="rz (half-plane rings) or 3d (cubic lattice)")
    p.add_argument("--dx-mm", required=True, type=float, help="particle spacing, mm")
    p.add_argument("--azimuth-deg", type=float, default=0.0, help="RZ: the half-plane's azimuth about the flight "
                   "axis from +y towards +z (default 0; RZ only)")
    p.add_argument("--out", default=DEFAULT_OUT, help="parent of the body directory (default spheral_output/body)")
    p.add_argument("--synthetic", choices=("radial",), default=None,
                   help="replace the frames' temperatures by a synthetic profile (removal still follows the frames' "
                   "body): radial = T_base + dT min(k / ramp, 1) (|x| / R)^2, R = 50 mm, T_base 300 K (check 4)")
    p.add_argument("--radial-delta-k", type=float, default=RADIAL_DELTA_K,
                   help="radial: the rise dT from the centre to r = R, whole kelvin (default 50)")
    p.add_argument("--radial-ramp-frames", type=int, default=RADIAL_RAMP_FRAMES,
                   help="radial: frames over which the profile ramps in from frame 0, held after (default 1)")
    p.add_argument("--force", action="store_true", help="rebuild even if built already")
    p.add_argument("--quiet", action="store_true", help="print only the summary and errors")
