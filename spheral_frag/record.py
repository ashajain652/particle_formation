"""The fragment record, the unresolved-debris log and the mass accounts (spec §11; M1 plan, Task 8).

Three parts:

1. **The runner -> analyse contract** (`Particles`, `Removal`, the run's `meta.json`). The runner (M3's
   `runner/fragments.py`; until then `tests/spheral_frag_synthetic.py::fake_run`) writes into its run directory

     checks/check_<n:06d>.npz   one `Particles` per fragment check n = 0, 1, ... (CHECK_ARRAYS)
     removed.npz                every particle mass that left the particles for the film account (REMOVED_ARRAYS)
     meta.json                  form, spacing, prepared name, brackets, seed, initial mass and enthalpy (META_KEYS)

   and `analyse` (Task 9) reads them. Plain numeric npz without pickles and JSON, because the writer runs under
   Spheral's Python (numpy 1.26) and the reader under drama_env (numpy 2; plan fact 10).

2. **Rows**: `fragment_row` (FRAGMENT_COLUMNS, spec §11.2 in its group order) and `debris_row` (DEBRIS_COLUMNS),
   from the particles of one group at one check. Definitions chosen where the spec leaves them open (plan Task 8,
   decision 7) are stated at each function.

3. **Accounts** (spec §11.3): at every check, starting mass = main body + resolved fragments (with attached dust) +
   unresolved debris + film account, summed exactly (`math.fsum`), each removal carrying its own enthalpy.

The material enters through the small `RecordMaterial` protocol (the viscosity and the fully liquid enthalpy), which
Task 3's `material.MaterialTable` satisfies; this module never imports `material`.

Numpy and the standard library only (the runner imports this module under Spheral's Python)."""
from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass, fields as dc_fields
from typing import Protocol, runtime_checkable

import numpy as np

from . import naming
from .contract import V_HAT

SCHEMA_VERSION = 1

# ---------------------------------------------------------------------------------------------------------------
# Constants (each with its source)
# ---------------------------------------------------------------------------------------------------------------

SIGMA_LIQUID = 0.80        # N/m, fully liquid 7075 (Bainbridge & Taylor 2013; spec §11.2; AA7075_scheil, decision 1)
F_L_FLUID = 0.5            # a particle is fluid when more than half liquid (spec §11.2; 895.1 K on the Scheil curve)
DUST_DAMAGE = 0.99         # dust: largest principal damage at least 0.99 (spec §11.1)
LINK_H = naming.BRACKET_DEFAULTS["link_h"]                # §11.1 linking distance in smoothing lengths
MIN_PARTICLES = naming.BRACKET_DEFAULTS["min_particles"]  # §11.1 resolution floor
CLEAR_H = 3.0              # §9.3 "clear": nearest particle more than 3 smoothing lengths from the main body
PILCH_ERDMAN = (12.0, 1.077, 1.6)   # We_crit = 12 (1 + 1.077 Oh^1.6) (Pilch and Erdman 1987; spec §11.2)
THOMSEN_P = 1.6075         # Knud Thomsen's ellipsoid-area approximation, largest relative error 1.061 %

REASON_CENTRE_CROSSING, REASON_SINK, REASON_FLOOR = 0, 1, 2      # Removal.reason (spec §9.1, §10, §14.1)
REASONS = {REASON_CENTRE_CROSSING: "centre_crossing", REASON_SINK: "sink", REASON_FLOOR: "floor"}

FORMS = naming.FORMS                         # ("rz", "3d")
ROUTES = ("separation", "slurry", "detachment", "synthetic")    # §9.2 routes 1-3; "synthetic" for fake runs
MECHANISMS = ("tearing", "slurry_breakup", "neck_failure", "ring_release")   # spec §11.2 Origin
PHASES = ("fluid", "mixed", "solid_or_mush")                    # spec §11.2 phase state
NOT_APPLICABLE = "not applicable"


@runtime_checkable
class RecordMaterial(Protocol):
    """What the record needs from the material (Task 3's `material.MaterialTable` provides both).

    `viscosity(T, shear_rate=0.0)`: Pa s, spec §8.3's slurry and liquid rows (Li et al. with the bridge, 1.3 mPa s
    above 908 K), inf below 50 % liquid. `h_full_liquid()`: J/kg from 300 K, the enthalpy at the liquidus, fully
    liquid (the zero of `Particles.h`)."""

    def viscosity(self, T, shear_rate=0.0): ...

    def h_full_liquid(self) -> float: ...


# ---------------------------------------------------------------------------------------------------------------
# 1. The runner -> analyse contract
# ---------------------------------------------------------------------------------------------------------------

# check_<n:06d>.npz: name -> (dtype, shape, units). n = number of particles at the check; () is a scalar.
CHECK_ARRAYS: dict[str, tuple[str, tuple, str]] = {
    "check": ("i8", (), "check index n (the file's own number)"),
    "frame": ("i8", (), "prepared frame k of the interval the check lies in (its start frame)"),
    "t": ("f8", (), "s, flight time of the check"),
    "id": ("i8", ("n",), "permanent particle ID (spec §6.2)"),
    "mass": ("f8", ("n",), "kg (rz form: the whole ring's)"),
    "x": ("f8", ("n", 3), "m, body frame, flight toward +x; rz form: (r, z) in x[:, :2], x[:, 2] = 0"),
    "v": ("f8", ("n", 3), "m/s relative to the main body; rz form: (v_r, v_z, 0)"),
    "T": ("f8", ("n",), "K"),
    "f_l": ("f8", ("n",), "-"),
    "h": ("f8", ("n",), "J/kg from 300 K"),
    "rho": ("f8", ("n",), "kg/m^3"),
    "damage": ("f8", ("n",), "-, largest principal damage"),
    "h_smooth": ("f8", ("n",), "m, smoothing length"),
    "group": ("i8", ("n",), "the fragment finder's per-call label"),
}

# removed.npz: name -> (dtype, shape, units). r = number of removals (a particle may appear several times: sink
# steps, then its floor removal).
REMOVED_ARRAYS: dict[str, tuple[str, tuple, str]] = {
    "id": ("i8", ("r",), "permanent particle ID"),
    "t": ("f8", ("r",), "s, flight time of the removal"),
    "check": ("i8", ("r",), "first check whose particles no longer hold this mass"),
    "mass": ("f8", ("r",), "kg removed (sink: m_old - m_new as booked; centre crossing, floor: all that was left)"),
    "h": ("f8", ("r",), "J/kg from 300 K carried by the removed mass (sink: the sprayed liquid's, §14.1)"),
    "reason": ("i1", ("r",), "0 centre crossing, 1 sink, 2 floor"),
}

# meta.json keys (all required)
META_KEYS = ("schema", "schema_version", "run_name", "prepared", "mode", "form", "dx_mm", "brackets", "seed",
             "initial_mass_kg", "initial_enthalpy_J", "n_particles")
META_SCHEMA = "spheral_frag.run"


def _as(a, dtype, ndim, name):
    out = np.ascontiguousarray(a, dtype=dtype)
    if out.ndim != ndim:
        raise ValueError("{} must have {} dimension(s), got shape {}".format(name, ndim, out.shape))
    return out


@dataclass
class Particles:                  # one fragment check
    t: float
    id: np.ndarray                # (n,) int64, permanent (spec §6.2)
    mass: np.ndarray              # kg (in rz form: the whole ring's)
    x: np.ndarray                 # (n, 3) m, body frame; rz form: (r, z) in x[:, :2], x[:, 2] = 0
    v: np.ndarray                 # (n, 3) m/s relative to the main body
    T: np.ndarray                 # K
    f_l: np.ndarray               # -
    h: np.ndarray                 # J/kg from 300 K
    rho: np.ndarray               # kg/m3
    damage: np.ndarray            # largest principal damage
    h_smooth: np.ndarray          # m, smoothing length
    group: np.ndarray             # int64, the finder's per-call label
    frame: int = -1               # prepared frame k of the check's interval (-1: unknown)
    check: int = -1               # check index n

    def __post_init__(self):
        self.t = float(self.t)
        self.frame, self.check = int(self.frame), int(self.check)
        self.id = _as(self.id, np.int64, 1, "id")
        self.group = _as(self.group, np.int64, 1, "group")
        n = self.id.shape[0]
        for name in ("mass", "T", "f_l", "h", "rho", "damage", "h_smooth"):
            setattr(self, name, _as(getattr(self, name), np.float64, 1, name))
        for name in ("x", "v"):
            setattr(self, name, _as(getattr(self, name), np.float64, 2, name))
        for name in ("mass", "T", "f_l", "h", "rho", "damage", "h_smooth", "group", "x", "v"):
            a = getattr(self, name)
            if a.shape[0] != n or (a.ndim == 2 and a.shape[1] != 3):
                raise ValueError("{} has shape {}, expected ({}{})".format(name, a.shape, n,
                                                                          ", 3" if a.ndim == 2 else ""))
        if len(np.unique(self.id)) != n:
            raise ValueError("permanent particle IDs are not unique")

    @property
    def n(self) -> int:
        return int(self.id.shape[0])


@dataclass
class Removal:                    # particles leaving to the film account
    id: np.ndarray                # (r,) int64
    t: np.ndarray                 # s
    check: np.ndarray             # int64: first check whose particles no longer hold this mass
    mass: np.ndarray              # kg
    h: np.ndarray                 # J/kg from 300 K
    reason: np.ndarray            # int8: 0 centre crossing, 1 sink, 2 floor

    def __post_init__(self):
        self.id = _as(self.id, np.int64, 1, "id")
        self.check = _as(self.check, np.int64, 1, "check")
        self.reason = _as(self.reason, np.int8, 1, "reason")
        for name in ("t", "mass", "h"):
            setattr(self, name, _as(getattr(self, name), np.float64, 1, name))
        r = self.id.shape[0]
        if any(getattr(self, f.name).shape[0] != r for f in dc_fields(self)):
            raise ValueError("Removal arrays differ in length")
        bad = set(np.unique(self.reason).tolist()) - set(REASONS)
        if bad:
            raise ValueError("unknown removal reason(s) {}".format(sorted(bad)))

    @classmethod
    def empty(cls) -> "Removal":
        return cls(*([np.zeros(0)] * 6))

    def __len__(self):
        return int(self.id.shape[0])

    def until(self, check: int) -> "Removal":
        """The removals booked into the film account at check `check` (their `check` <= it)."""
        return self.subset(self.check <= check)

    def subset(self, sel) -> "Removal":
        return Removal(*(getattr(self, f.name)[sel] for f in dc_fields(self)))

    @classmethod
    def concat(cls, parts) -> "Removal":
        parts = [p for p in parts if len(p)] or [cls.empty()]
        return cls(*(np.concatenate([getattr(p, f.name) for p in parts]) for f in dc_fields(cls)))


def check_path(run_dir, n: int) -> str:
    return os.path.join(run_dir, "checks", "check_{:06d}.npz".format(n))


def _savez(path, arrays):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp.npz"
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def _load_npz(path, schema, what):
    with np.load(path, allow_pickle=False) as z:
        missing = sorted(set(schema) - set(z.files))
        if missing:
            raise ValueError("{} {} lacks {}".format(what, path, ", ".join(missing)))
        out = {}
        for name, (dtype, shape, _) in schema.items():
            a = z[name]
            if a.dtype != np.dtype(dtype) or a.ndim != len(shape):
                raise ValueError("{} {}: {} is {} with shape {}, expected {} of {} dimension(s)".format(
                    what, path, name, a.dtype, a.shape, dtype, len(shape)))
            out[name] = a
        return out


def write_check(run_dir, p: Particles, n: int | None = None) -> str:
    """Write `p` as checks/check_<n:06d>.npz (n defaults to p.check) and return the path."""
    n = p.check if n is None else int(n)
    if n < 0:
        raise ValueError("a check needs a non-negative index, got {}".format(n))
    arrays = {"check": np.int64(n), "frame": np.int64(p.frame), "t": np.float64(p.t)}
    for name in CHECK_ARRAYS:
        if name not in arrays:
            arrays[name] = np.asarray(getattr(p, name), dtype=CHECK_ARRAYS[name][0])
    path = check_path(run_dir, n)
    _savez(path, arrays)
    return path


def read_check(path) -> Particles:
    a = _load_npz(path, CHECK_ARRAYS, "check file")
    return Particles(t=float(a["t"]), frame=int(a["frame"]), check=int(a["check"]),
                     **{k: a[k] for k in CHECK_ARRAYS if k not in ("t", "frame", "check")})


def list_checks(run_dir) -> list[tuple[int, str]]:
    """(n, path) of every check file in the run directory, in order of n."""
    d = os.path.join(run_dir, "checks")
    if not os.path.isdir(d):
        return []
    out = []
    for fn in os.listdir(d):
        if fn.startswith("check_") and fn.endswith(".npz") and fn[6:-4].isdigit():
            out.append((int(fn[6:-4]), os.path.join(d, fn)))
    return sorted(out)


def write_removed(path, removed: Removal) -> str:
    _savez(path, {name: np.asarray(getattr(removed, name), dtype=spec[0]) for name, spec in REMOVED_ARRAYS.items()})
    return path


def read_removed(path) -> Removal:
    a = _load_npz(path, REMOVED_ARRAYS, "removal file")
    return Removal(**a)


def make_meta(run_name, prepared, mode, form, dx_mm, brackets, seed, initial_mass_kg, initial_enthalpy_J,
              n_particles, **extra) -> dict:
    """The run's meta.json content (META_KEYS plus any `extra` provenance, e.g. the Spheral pin)."""
    meta = {"schema": META_SCHEMA, "schema_version": SCHEMA_VERSION, "run_name": run_name, "prepared": prepared,
            "mode": mode, "form": form, "dx_mm": float(dx_mm), "brackets": dict(brackets), "seed": int(seed),
            "initial_mass_kg": float(initial_mass_kg), "initial_enthalpy_J": float(initial_enthalpy_J),
            "n_particles": int(n_particles)}
    meta.update(extra)
    validate_meta(meta)
    return meta


def validate_meta(meta: dict) -> None:
    """ValueError unless `meta` has every META_KEYS entry and agrees with its own run name."""
    missing = [k for k in META_KEYS if k not in meta]
    if missing:
        raise ValueError("meta.json lacks {}".format(", ".join(missing)))
    if meta["schema"] != META_SCHEMA or meta["schema_version"] != SCHEMA_VERSION:
        raise ValueError("meta.json schema {!r} v{} is not {!r} v{}".format(
            meta["schema"], meta["schema_version"], META_SCHEMA, SCHEMA_VERSION))
    if meta["form"] not in FORMS:
        raise ValueError("form must be one of {}, got {!r}".format(FORMS, meta["form"]))
    if not (meta["initial_mass_kg"] > 0.0):
        raise ValueError("initial_mass_kg must be positive")
    parsed = naming.parse_run_name(meta["run_name"])
    brackets = dict(naming.BRACKET_DEFAULTS, **meta["brackets"])
    for key, want in (("prepared", meta["prepared"]), ("mode", meta["mode"]), ("form", meta["form"]),
                      ("dx_mm", float(meta["dx_mm"])), ("seed", meta["seed"]), ("brackets", brackets)):
        if parsed[key] != want:
            raise ValueError("meta.json {} = {!r} disagrees with the run name's {!r}".format(key, want, parsed[key]))


def write_meta(path, meta: dict) -> str:
    validate_meta(meta)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2, sort_keys=False, allow_nan=False)
        fh.write("\n")
    return path


def read_meta(path) -> dict:
    with open(path, encoding="utf-8") as fh:
        meta = json.load(fh)
    validate_meta(meta)
    return meta


# ---------------------------------------------------------------------------------------------------------------
# The flight state at release (spec §11.2 "When")
# ---------------------------------------------------------------------------------------------------------------

# Names in the prepared flight table (Task 9's flight.npz): the history's contract keys plus the two columns
# prepare derives (T_air_K from fe.air_temperature, deceleration_ms2 = load_factor_g * g0).
FLIGHT_TIME = "time_s"
FLIGHT_STATE = ("altitude_km", "velocity_kms", "flight_path_deg", "density_kgm3", "T_air_K", "deceleration_ms2")


class FlightTable:
    """The prepared flight table, interpolated linearly in time (spec §11.2; plan Task 8)."""

    def __init__(self, columns: dict):
        missing = [k for k in (FLIGHT_TIME,) + FLIGHT_STATE if k not in columns]
        if missing:
            raise ValueError("flight table lacks {}".format(", ".join(missing)))
        self.time_s = np.asarray(columns[FLIGHT_TIME], dtype=np.float64)
        if self.time_s.ndim != 1 or len(self.time_s) < 1 or np.any(np.diff(self.time_s) <= 0):
            raise ValueError("flight table times must be one-dimensional and increasing")
        self.columns = {k: np.asarray(columns[k], dtype=np.float64) for k in FLIGHT_STATE}

    @classmethod
    def load(cls, path) -> "FlightTable":
        with np.load(path, allow_pickle=False) as z:
            return cls({k: z[k] for k in z.files})

    def at(self, t: float) -> dict:
        t = float(t)
        tol = 1e-9 * max(1.0, abs(t))
        if not (self.time_s[0] - tol <= t <= self.time_s[-1] + tol):
            raise ValueError("t = {} s is outside the flight table's {}-{} s".format(t, self.time_s[0],
                                                                                  self.time_s[-1]))
        return {k: float(np.interp(t, self.time_s, v)) for k, v in self.columns.items()}


def _flight_state(flight, t):
    if hasattr(flight, "at"):
        return flight.at(t)
    missing = [k for k in FLIGHT_STATE if k not in flight]
    if missing:
        raise ValueError("flight state lacks {}".format(", ".join(missing)))
    return {k: float(flight[k]) for k in FLIGHT_STATE}


# ---------------------------------------------------------------------------------------------------------------
# 2. Groups: classification, following, clearance
# ---------------------------------------------------------------------------------------------------------------

def _nearest(xa, ha, xb, hb, chunk=2048):
    """For each point of a: (index into b of the nearest in units of the pair's mean smoothing length, that
    normalised distance). Brute force in chunks (numpy only; the runner uses Spheral's own neighbour search)."""
    idx = np.empty(len(xa), dtype=np.int64)
    dist = np.empty(len(xa))
    for s in range(0, len(xa), chunk):
        d = np.linalg.norm(xa[s:s + chunk, None, :] - xb[None, :, :], axis=2)
        q = d / (0.5 * (ha[s:s + chunk, None] + hb[None, :]))
        j = np.argmin(q, axis=1)
        idx[s:s + chunk] = j
        dist[s:s + chunk] = q[np.arange(len(j)), j]
    return idx, dist


def classify_groups(p: Particles, min_particles=MIN_PARTICLES, dust_damage=DUST_DAMAGE, link_h=LINK_H) -> dict:
    """Main body, fragments and debris at one check (spec §11.1).

    Dust is every particle with damage >= `dust_damage`, whatever its finder label (dust links nothing). The intact
    particles keep the finder's labels; the intact group with the most mass is the main body (ties: the lowest
    label). A dust particle attaches to the group of its nearest intact particle when it lies within
    `link_h` times the pair's mean smoothing length 0.5 (h_i + h_j) of it; otherwise it is isolated dust. An intact
    group other than the main body with at least `min_particles` intact particles is a fragment; one with fewer is
    debris (with its attached dust). All masks are boolean over the particles of `p`.

    Returns {"main_label", "main" (mask with attached dust), "fragments" {label: mask with attached dust},
    "fragment_dust" {label: mask of its attached dust}, "debris_groups" {label: mask with attached dust},
    "debris_dust" (mask of isolated dust), "dust" (mask of all dust)}."""
    dust = p.damage >= dust_damage
    intact = ~dust
    if not intact.any():
        raise ValueError("no intact particle: there is no main body")
    labels = p.group[intact]
    ulab, inv = np.unique(labels, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    bounds = np.searchsorted(inv[order], np.arange(len(ulab) + 1))
    m_sorted = p.mass[intact][order]
    gmass = np.array([math.fsum(m_sorted[bounds[i]:bounds[i + 1]]) for i in range(len(ulab))])
    main_label = int(ulab[np.flatnonzero(gmass == gmass.max())[0]])

    owner = np.full(p.n, np.iinfo(np.int64).min, dtype=np.int64)   # group each particle counts with
    owner[intact] = p.group[intact]
    di = np.flatnonzero(dust)
    ii = np.flatnonzero(intact)
    attached = np.zeros(p.n, dtype=bool)
    if len(di):
        j, q = _nearest(p.x[di], p.h_smooth[di], p.x[ii], p.h_smooth[ii])
        ok = q <= link_h
        owner[di[ok]] = p.group[ii[j[ok]]]
        attached[di[ok]] = True

    out = {"main_label": main_label, "main": owner == main_label, "fragments": {}, "fragment_dust": {},
           "debris_groups": {}, "debris_dust": dust & ~attached, "dust": dust}
    for lab in ulab.tolist():
        if lab == main_label:
            continue
        sel = owner == lab
        if np.count_nonzero(intact & sel) >= min_particles:
            out["fragments"][lab] = sel
            out["fragment_dust"][lab] = sel & dust
        else:
            out["debris_groups"][lab] = sel
    return out


def group_ids(p: Particles, groups: dict[int, np.ndarray]) -> dict[int, np.ndarray]:
    """{label: sorted permanent IDs} of a {label: mask} dict (the input of `match_groups`)."""
    return {lab: np.sort(p.id[sel]) for lab, sel in groups.items()}


def match_groups(prev: dict[int, np.ndarray], curr: dict[int, np.ndarray]) -> dict[int, int]:
    """Follow groups from one check to the next by their permanent IDs (spec §11.1).

    `prev` and `curr` map each check's per-call labels to the permanent IDs of the group's particles. Returns
    {prev label: curr label}: the current group holding the largest share of the previous group's particles that
    still exist (ties: the lowest current label). A previous group none of whose particles remain is absent. A group
    that splits maps to the piece holding most of it; two groups that merge map to the same label."""
    if not curr:
        return {}
    clabels = np.array(sorted(curr), dtype=np.int64)
    cid = np.concatenate([np.asarray(curr[c], dtype=np.int64) for c in clabels])
    clab = np.concatenate([np.full(len(curr[c]), c, dtype=np.int64) for c in clabels])
    if len(np.unique(cid)) != len(cid):
        raise ValueError("a particle ID belongs to two current groups")
    order = np.argsort(cid, kind="stable")
    cid, clab = cid[order], clab[order]
    out = {}
    for plabel in sorted(prev):
        ids = np.asarray(prev[plabel], dtype=np.int64)
        pos = np.searchsorted(cid, ids)
        pos = np.minimum(pos, len(cid) - 1)
        found = cid[pos] == ids
        if not found.any():
            continue
        labs, counts = np.unique(clab[pos[found]], return_counts=True)
        out[int(plabel)] = int(labs[np.flatnonzero(counts == counts.max())[0]])
    return out


def is_clear(p: Particles, sel, main_sel, clear_h=CLEAR_H) -> bool:
    """A group is clear (spec §9.3) when its nearest particle is more than `clear_h` mean smoothing lengths
    0.5 (h_i + h_j) from the main body's nearest particle and it moves away: its mass-weighted velocity relative to
    the main body's, projected on the direction from that nearest main particle to that nearest group particle, is
    positive. In rz form the distance is in the meridional half-plane, where two bodies of revolution are nearest."""
    sel, main_sel = np.asarray(sel, bool), np.asarray(main_sel, bool)
    if not sel.any() or not main_sel.any() or (sel & main_sel).any():
        raise ValueError("is_clear needs two non-empty disjoint selections")
    gi, mi = np.flatnonzero(sel), np.flatnonzero(main_sel)
    j, q = _nearest(p.x[gi], p.h_smooth[gi], p.x[mi], p.h_smooth[mi])
    a = int(np.argmin(q))
    if q[a] <= clear_h:
        return False
    direction = p.x[gi[a]] - p.x[mi[j[a]]]
    v_rel = _mean(p.v[gi], p.mass[gi]) - _mean(p.v[mi], p.mass[mi])
    return bool(np.dot(v_rel, direction) > 0.0)


# ---------------------------------------------------------------------------------------------------------------
# Numbers for the rows
# ---------------------------------------------------------------------------------------------------------------

def _mean(a, w):
    """Weighted mean along axis 0, the weights summed exactly."""
    a, w = np.asarray(a, float), np.asarray(w, float)
    W = math.fsum(w)
    if a.ndim == 1:
        return math.fsum(a * w) / W
    return np.array([math.fsum(a[:, k] * w) for k in range(a.shape[1])]) / W


def principal_lengths(x, w) -> np.ndarray:
    """2 sqrt(5 lambda_i) from the `w`-weighted covariance of the points `x` (n, 3), in descending order: exact
    for a uniform solid ellipsoid (variance a^2/5 along a semi-axis a)."""
    lam = np.linalg.eigvalsh(_cov(x, w))[::-1]
    return 2.0 * np.sqrt(5.0 * np.clip(lam, 0.0, None))


def _cov(x, w):
    x, w = np.asarray(x, float), np.asarray(w, float)
    d = x - _mean(x, w)
    return (d * w[:, None]).T @ d / math.fsum(w)


def ellipsoid_area(a, b, c, p=THOMSEN_P) -> float:
    """Surface area of the ellipsoid of semi-axes a, b, c by Thomsen's approximation (approximate: at most
    1.061 % from the exact area)."""
    ap, bp, cp = a ** p, b ** p, c ** p
    return float(4.0 * math.pi * ((ap * bp + ap * cp + bp * cp) / 3.0) ** (1.0 / p))


def ellipse_perimeter(a, b) -> float:
    """Ramanujan's second approximation (relative error below 4e-5 up to a/b = 5)."""
    hh = ((a - b) / (a + b)) ** 2 if a + b > 0 else 0.0
    return float(math.pi * (a + b) * (1.0 + 3.0 * hh / (10.0 + math.sqrt(4.0 - 3.0 * hh))))


def weber(rho_inf, V, d, sigma=SIGMA_LIQUID):
    """We = rho_inf V^2 d / sigma (spec §11.2; Step 4's convention with d = 2r)."""
    return rho_inf * V * V * d / sigma


def ohnesorge(mu, rho, d, sigma=SIGMA_LIQUID):
    """Oh = mu / sqrt(rho sigma d)."""
    return mu / np.sqrt(rho * sigma * d)


def breakup_threshold(oh):
    """We_crit = 12 (1 + 1.077 Oh^1.6) (Pilch and Erdman 1987); 12 for an inviscid liquid (Step 3's threshold)."""
    c, k, e = PILCH_ERDMAN
    return c * (1.0 + k * np.power(oh, e))


def breakup_flag(we, oh) -> bool:
    """We > 12 (1 + 1.077 Oh^1.6)."""
    return bool(we > breakup_threshold(oh))


def phase_state(f_l, mass) -> tuple[str, float]:
    """("fluid" if every particle is more than half liquid, "mixed" if some are, "solid_or_mush" otherwise;
    the fluid mass fraction)."""
    f_l, mass = np.asarray(f_l, float), np.asarray(mass, float)
    fluid = f_l > F_L_FLUID
    frac = math.fsum(mass[fluid]) / math.fsum(mass)
    state = "fluid" if fluid.all() else ("mixed" if fluid.any() else "solid_or_mush")
    return state, frac


def _check_form(form):
    if form not in FORMS:
        raise ValueError("form must be one of {}, got {!r}".format(FORMS, form))


def _position_and_velocity(p, sel, form, ring):
    """Body-frame centroid and mean velocity (3D: mass-weighted). rz form: the meridional half-plane is placed at
    phi = 0 of the body frame, (r, z) -> (x = z, y = r, z = 0); a ring's position is its cross-section's
    mass-weighted (z, r), a cap's lies on the axis (z, 0); the velocity (v_r, v_z) -> (v_z, v_r, 0)."""
    m = p.mass[sel]
    if form == "3d":
        return _mean(p.x[sel], m), _mean(p.v[sel], m)
    r_bar, z_bar = _mean(p.x[sel][:, :2], m)
    vr, vz = _mean(p.v[sel][:, :2], m)
    return np.array([z_bar, r_bar if ring else 0.0, 0.0]), np.array([vz, vr, 0.0])


def _angles_deg(pos, form):
    """(angle from the nose: between the position from the body-frame origin and V_HAT; angle around the axis:
    atan2(z, y), NaN in rz form, where every angle is the same)."""
    v = np.asarray(V_HAT, float)
    along = float(np.dot(pos, v))
    across = float(np.linalg.norm(pos - along * v))
    theta = math.degrees(math.atan2(across, along)) if (along or across) else float("nan")
    phi = math.degrees(math.atan2(pos[2], pos[1])) if form == "3d" else float("nan")
    return theta, phi


def _thermal(p, sel, table):
    m = p.mass[sel]
    out = {"h_mean_Jkg": _mean(p.h[sel], m), "T_mean_K": _mean(p.T[sel], m), "T_min_K": float(p.T[sel].min()),
           "T_max_K": float(p.T[sel].max()), "f_l_mean": _mean(p.f_l[sel], m)}
    if table is None:
        out["heat_to_melt_J"] = float("nan")
    else:
        out["heat_to_melt_J"] = math.fsum(m * np.maximum(float(table.h_full_liquid()) - p.h[sel], 0.0))
    return out


def _shape(p, sel, form):
    """Volume, equivalent diameter, principal lengths, approximate surface area, ring flag and radius.

    3D: lengths 2 sqrt(5 lambda) of the mass-weighted covariance; area of the ellipsoid of those semi-axes
    (Thomsen). rz form: a ring when no member touches the axis (every r > its smoothing length). A ring's lengths are
    its cross-section's, 2 sqrt(4 lambda) of the (r, z) covariance weighted by cross-section area m / (rho r) (exact
    for a uniform ellipse), plus its circumference 2 pi R with R the mass-weighted r (plan Task 8); its area by Pappus,
    the cross-section's perimeter (Ramanujan) times 2 pi r_c with r_c the cross-section's area centroid. A cap's
    lengths are those of the body of revolution, from its mass-weighted moments (axial variance var(z), transverse
    variance <r^2> / 2 on each of the two transverse axes), and its area the ellipsoid's of those semi-axes."""
    m, rho = p.mass[sel], p.rho[sel]
    volume = math.fsum(m / rho)
    d_eq = (6.0 * volume / math.pi) ** (1.0 / 3.0)
    ring, ring_radius = False, float("nan")
    if form == "3d":
        L = principal_lengths(p.x[sel], m)
        area = ellipsoid_area(*(L / 2.0))
    else:
        r, z = p.x[sel][:, 0], p.x[sel][:, 1]
        ring = bool(np.all(r > p.h_smooth[sel]))
        if ring:
            ring_radius = _mean(r, m)
            w = m / (rho * r)                      # cross-section area per particle (uniform in the plane)
            lam = np.linalg.eigvalsh(_cov(np.column_stack([r, z]), w))[::-1]
            ab = 2.0 * np.sqrt(4.0 * np.clip(lam, 0.0, None)) / 2.0       # cross-section semi-axes
            L = np.sort(np.array([2.0 * ab[0], 2.0 * ab[1], 2.0 * math.pi * ring_radius]))[::-1]
            area = ellipse_perimeter(*ab) * 2.0 * math.pi * _mean(r, w)
        else:
            var_axial = _mean((z - _mean(z, m)) ** 2, m)
            var_trans = _mean(r * r, m) / 2.0
            L = np.sort(2.0 * np.sqrt(5.0 * np.array([var_axial, var_trans, var_trans])))[::-1]
            area = ellipsoid_area(*(L / 2.0))
    return {"volume_m3": volume, "d_eq_m": d_eq, "L1_m": float(L[0]), "L2_m": float(L[1]), "L3_m": float(L[2]),
            "area_m2": area}, ring, ring_radius


# ---------------------------------------------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------------------------------------------

BRACKET_COLUMNS = ["bracket_" + k for k in naming.BRACKET_DEFAULTS]

FRAGMENT_GROUPS: dict[str, list[str]] = {        # spec §11.2, in its group order
    "When": ["t_release_s", "frame", *FLIGHT_STATE],
    "Where and how fast": ["x_m", "y_m", "z_m", "theta_deg", "phi_deg", "vx_ms", "vy_ms", "vz_ms", "v_rel_ms"],
    "Size and shape": ["mass_kg", "count", "n_particles", "dust_mass_kg", "volume_m3", "d_eq_m", "L1_m", "L2_m",
                       "L3_m", "area_m2"],
    "Thermal state": ["h_mean_Jkg", "T_mean_K", "T_min_K", "T_max_K", "f_l_mean", "heat_to_melt_J"],
    "Origin": ["route", "mechanism", "resolved", "ring", "ring_radius_m"],
    "Phase and breakup": ["phase", "fluid_mass_fraction", "weber", "ohnesorge", "mu_Pas", "sigma_Nm",
                          "breakup_threshold", "breakup"],
    "Provenance": ["run", "window", "form", "dx_mm", "seed", *BRACKET_COLUMNS],
}
FRAGMENT_COLUMNS: list[str] = ["fragment"] + [c for cols in FRAGMENT_GROUPS.values() for c in cols]

DEBRIS_COLUMNS: list[str] = [
    "debris", "kind", "t_s", "frame",
    "x_m", "y_m", "z_m", "theta_deg", "phi_deg", "vx_ms", "vy_ms", "vz_ms", "v_rel_ms",
    "h_mean_Jkg", "T_mean_K", "T_min_K", "T_max_K", "f_l_mean", "heat_to_melt_J",
    "mass_kg", "n_particles", "size_upper_m",
]
DEBRIS_KINDS = ("dust", "small_group")

PROVENANCE_KEYS = ("run", "window", "dx_mm", "seed", "brackets")


def fragment_row(p: Particles, sel, dust_sel, table: RecordMaterial, flight, route, mechanism, provenance: dict,
                 form="3d", number=0, sigma=SIGMA_LIQUID) -> dict:
    """One row of the fragment record (FRAGMENT_COLUMNS) for the particles `sel` of `p` (attached dust included;
    `dust_sel` marks the dust among them), at the check `p` holds (the release check).

    `flight`: a FlightTable (evaluated at p.t) or a mapping of FLIGHT_STATE. `provenance`: run, window, dx_mm, seed
    and brackets (missing brackets take naming.BRACKET_DEFAULTS). `count` is 1 (Step 4's record counts the drops a
    row stands for; a resolved fragment is one). Weber and Ohnesorge numbers only for fluid fragments, NaN and
    breakup "not applicable" otherwise: d the equivalent diameter, V the flight speed and rho_inf the air density at
    release, rho the fragment's own density M / V, mu = table.viscosity(T_mw, shear_rate=0) at the mass-weighted
    temperature (Li et al. in the low-shear limit; decision 7), sigma 0.80 N/m."""
    _check_form(form)
    if route not in ROUTES:
        raise ValueError("route must be one of {}, got {!r}".format(ROUTES, route))
    if mechanism not in MECHANISMS:
        raise ValueError("mechanism must be one of {}, got {!r}".format(MECHANISMS, mechanism))
    if not isinstance(table, RecordMaterial):
        raise TypeError("table must provide viscosity(T, shear_rate) and h_full_liquid() (RecordMaterial)")
    missing = [k for k in PROVENANCE_KEYS if k not in provenance]
    if missing:
        raise ValueError("provenance lacks {}".format(", ".join(missing)))
    sel = np.asarray(sel, bool)
    dust_sel = np.asarray(dust_sel, bool) & sel
    if not sel.any():
        raise ValueError("empty fragment")
    m = p.mass[sel]
    mass = math.fsum(m)

    shape, ring, ring_radius = _shape(p, sel, form)
    pos, vel = _position_and_velocity(p, sel, form, ring)
    theta, phi = _angles_deg(pos, form)
    thermal = _thermal(p, sel, table)
    state = _flight_state(flight, p.t)
    phase, fluid_frac = phase_state(p.f_l[sel], m)

    we = oh = mu = thr = float("nan")
    breakup = NOT_APPLICABLE
    if phase == "fluid":
        rho_frag = mass / shape["volume_m3"]
        mu = float(np.asarray(table.viscosity(thermal["T_mean_K"], shear_rate=0.0)))
        we = float(weber(state["density_kgm3"], state["velocity_kms"] * 1e3, shape["d_eq_m"], sigma))
        oh = float(ohnesorge(mu, rho_frag, shape["d_eq_m"], sigma))
        thr = float(breakup_threshold(oh))
        breakup = "yes" if breakup_flag(we, oh) else "no"

    brackets = dict(naming.BRACKET_DEFAULTS, **provenance["brackets"])
    row = {"fragment": int(number), "t_release_s": p.t, "frame": p.frame, **state,
           "x_m": pos[0], "y_m": pos[1], "z_m": pos[2], "theta_deg": theta, "phi_deg": phi,
           "vx_ms": vel[0], "vy_ms": vel[1], "vz_ms": vel[2], "v_rel_ms": float(np.linalg.norm(vel)),
           "mass_kg": mass, "count": 1, "n_particles": int(np.count_nonzero(sel)),
           "dust_mass_kg": math.fsum(p.mass[dust_sel]), **shape, **thermal,
           "route": route, "mechanism": mechanism, "resolved": True, "ring": ring, "ring_radius_m": ring_radius,
           "phase": phase, "fluid_mass_fraction": fluid_frac, "weber": we, "ohnesorge": oh, "mu_Pas": mu,
           "sigma_Nm": float(sigma), "breakup_threshold": thr, "breakup": breakup,
           "run": provenance["run"], "window": provenance["window"], "form": form,
           "dx_mm": float(provenance["dx_mm"]), "seed": int(provenance["seed"]),
           **{"bracket_" + k: brackets[k] for k in naming.BRACKET_DEFAULTS}}
    return {c: row[c] for c in FRAGMENT_COLUMNS}


def _extent(x, n_exact=2000):
    """Largest distance between two of the points (exact up to n_exact points, else the bounding box diagonal, an
    upper bound)."""
    if len(x) < 2:
        return 0.0
    if len(x) > n_exact:
        return float(np.linalg.norm(x.max(axis=0) - x.min(axis=0)))
    best = 0.0
    for i in range(len(x) - 1):
        best = max(best, float(np.max(np.linalg.norm(x[i + 1:] - x[i], axis=1))))
    return best


def debris_row(p: Particles, sel, flight, dx, table: RecordMaterial | None = None, form="3d", number=0,
               kind="small_group") -> dict:
    """One row of the unresolved-debris log (DEBRIS_COLUMNS): the fragment record's time, position, velocity and
    thermal columns, the mass, and an upper bound on size: the group's largest extent plus one spacing `dx` [m]
    (rz form: the extent of the cross-section in the meridional plane). `flight` is accepted for the time check
    (the log carries no flight state); `table` gives the heat to melt (NaN without it)."""
    _check_form(form)
    if kind not in DEBRIS_KINDS:
        raise ValueError("kind must be one of {}, got {!r}".format(DEBRIS_KINDS, kind))
    sel = np.asarray(sel, bool)
    if not sel.any():
        raise ValueError("empty debris selection")
    if flight is not None:
        _flight_state(flight, p.t)
    ring = form == "rz" and bool(np.all(p.x[sel][:, 0] > p.h_smooth[sel]))
    pos, vel = _position_and_velocity(p, sel, form, ring)
    theta, phi = _angles_deg(pos, form)
    row = {"debris": int(number), "kind": kind, "t_s": p.t, "frame": p.frame,
           "x_m": pos[0], "y_m": pos[1], "z_m": pos[2], "theta_deg": theta, "phi_deg": phi,
           "vx_ms": vel[0], "vy_ms": vel[1], "vz_ms": vel[2], "v_rel_ms": float(np.linalg.norm(vel)),
           **_thermal(p, sel, table), "mass_kg": math.fsum(p.mass[sel]), "n_particles": int(np.count_nonzero(sel)),
           "size_upper_m": _extent(p.x[sel]) + float(dx)}
    return {c: row[c] for c in DEBRIS_COLUMNS}


def _cell(v):
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, (float, np.floating)):
        return repr(float(v))                 # shortest text that round-trips; nan as 'nan'
    return str(v)


def write_rows(path, rows: list[dict], columns: list[str]) -> str:
    """CSV with `columns` as header; floats written so that they read back bitwise, booleans as true/false."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(columns)
        for row in rows:
            if list(row) != list(columns):
                raise ValueError("row keys differ from the columns")
            w.writerow([_cell(row[c]) for c in columns])
    return path


# ---------------------------------------------------------------------------------------------------------------
# 3. Mass accounts (spec §11.3)
# ---------------------------------------------------------------------------------------------------------------

ACCOUNTS = ("main", "fragments", "debris", "film")


def accounts(m0, main, fragments, debris, film: Removal) -> dict:
    """Mass and enthalpy of every account at one check, and the mass residual.

    `main`, `fragments`, `debris`: (mass, h) arrays of the particles in each account (kg, J/kg); `film`: the
    removals booked by this check (`Removal.until(n)`). Every sum is exact (`math.fsum` of the stored floats, then
    one rounding), so the residual is the booked floats' own: residual_kg = (sum of every particle and removal mass)
    - m0, computed as one exact sum. The enthalpy of each account is fsum(m h); the film's is split by reason, and
    carries every removal's own h, not the removed particle's average. Enthalpy is not conserved (the replay
    imposes the frame temperatures), so it is reported, not balanced."""
    parts = {"main": main, "fragments": fragments, "debris": debris, "film": (film.mass, film.h)}
    out = {}
    masses = []
    for name, (m, h) in parts.items():
        m, h = np.asarray(m, float), np.asarray(h, float)
        masses.append(m)
        out["mass_{}_kg".format(name)] = math.fsum(m)
        out["enthalpy_{}_J".format(name)] = math.fsum(m * h)
    for code, reason in REASONS.items():
        s = film.reason == code
        out["mass_film_{}_kg".format(reason)] = math.fsum(film.mass[s])
        out["enthalpy_film_{}_J".format(reason)] = math.fsum(film.mass[s] * film.h[s])
    out["m0_kg"] = float(m0)
    out["mass_total_kg"] = math.fsum(np.concatenate(masses))
    out["residual_kg"] = math.fsum(np.concatenate(masses + [np.array([-float(m0)])]))
    out["rel_residual"] = out["residual_kg"] / float(m0)
    out["enthalpy_total_J"] = math.fsum([out["enthalpy_{}_J".format(a)] for a in ACCOUNTS])
    return out


def accounts_at(m0, p: Particles, groups: dict, removed: Removal, check: int | None = None) -> dict:
    """`accounts` for the particles of one check classified by `classify_groups`: main body, fragments with their
    attached dust, debris (sub-floor groups and isolated dust), and the removals booked by that check."""
    check = p.check if check is None else check
    n = p.n
    frag = np.zeros(n, bool)
    for s in groups["fragments"].values():
        frag |= s
    deb = groups["debris_dust"].copy()
    for s in groups["debris_groups"].values():
        deb |= s
    main = groups["main"]
    if np.count_nonzero(main) + np.count_nonzero(frag) + np.count_nonzero(deb) != n or (main & (frag | deb)).any() \
            or (frag & deb).any():
        raise ValueError("the groups do not partition the particles")
    out = accounts(m0, (p.mass[main], p.h[main]), (p.mass[frag], p.h[frag]), (p.mass[deb], p.h[deb]),
                   removed.until(check))
    out.update({"check": int(check), "t": p.t, "n_main": int(np.count_nonzero(main)),
                "n_fragments": len(groups["fragments"]), "n_debris_groups": len(groups["debris_groups"]),
                "n_debris_dust": int(np.count_nonzero(groups["debris_dust"])), "n_removals": len(removed.until(check))})
    return out
