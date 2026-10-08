"""The frame contract: what `spheral_frag` reads from a finite-element run directory (M1 plan, "The frame contract").

`FE_FIELDS` is the one table in the package that spells a finite-element name. Every other module asks
`fe_name(key)` for it, so a field renamed by the finite-element export is fixed here and nowhere else (M1 Task 10
edits only this table: names and the `confirmed` flags).

Where the items live in a run directory `<outdir>/<run>` (plan, facts 1-5):

  node   point data of `<run>/vtk/field_<k>.vtu` (all mesh nodes, dead ones included, unreferenced)
  tet    cell data of the same vtu (the active tetrahedra only)
  patch  cell data of `<run>/vtk/surface_<k>.vtp` (the outward surface triangles, on the same node array)
  frame  per frame, from the `field.pvd` collection (time printed with 6 significant digits)
  history  columns of `<outdir>/<run>.csv` (one row per macro step, 9 significant digits)
  run    keys of `<outdir>/<run>.json`; a dotted name is a path into the JSON (`settings.seed`)

Units of the per-step fields (review focus 3): `release_rate` is kg/m^2 **per macro step**, not per second;
`deep_thickness` is a mass per area expressed as a thickness, m_d / (rho_l A); `p_w` = 0 means "not evaluated"
(the export writes 0.0 where the surface flow was not evaluated, e.g. before the first step or on a face exposed
by that step's deaths; under decision 4 of the plan it writes the loads on every step that evaluated the flow,
with or without film).

Numpy and the standard library only (the runner imports this module under Spheral's Python)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

CONTRACT_VERSION = 1
V_HAT = (1.0, 0.0, 0.0)        # body-frame flight direction: the default of every body class; exported as v_hat (decision 4, confirm)
LOCATIONS = ("node", "tet", "patch", "frame", "history", "run")


class ContractError(Exception):
    """The finite-element run does not meet the contract (a required item missing, frame times unmatched, ...):
    `prepare` exits 1 with the message."""


@dataclass(frozen=True)
class FieldSpec:
    key: str            # spheral_frag's name
    name: str           # the finite-element file's name (the one place to change it)
    location: str       # "node" | "tet" | "patch" | "frame" | "history" | "run"
    units: str
    required: bool
    nan_allowed: bool
    confirmed: bool     # checked against real frames (Task 10)
    consumer: str       # spec section that reads it
    note: str = ""      # status as read from the code (plan, contract table); what Task 10 must confirm


def _f(key, name, location, units, required, consumer, note="", nan_allowed=False, confirmed=False):
    return FieldSpec(key, name, location, units, required, nan_allowed, confirmed, consumer, note)


_HISTORY = (
    ("time_s", "s", "§7.1, §9.1", "on main"),
    ("altitude_km", "km", "§11.2", "on main"),
    ("velocity_kms", "km/s", "§7.1, §11.2", "on main"),
    ("mass_kg", "kg", "§7.1, mass check", "on main"),
    ("flight_path_deg", "deg", "§11.2", "on main"),
    ("lat_deg", "deg", "§11.2, air temperature", "on main"),
    ("lon_deg", "deg", "§11.2, air temperature", "on main"),
    ("density_kgm3", "kg/m^3", "§11.2 (Weber number)", "on main"),
    ("dynamic_pressure_Pa", "Pa", "§7", "on main"),
    ("load_factor_g", "g0 (|a_drag| / g0)", "§7.1 body force; drag = mass_kg * load_factor_g * g0",
     "on main; there is no deceleration or drag column (fact 4)"),
    ("p_w_stag_Pa", "Pa", "§7.2 lee base pressure reference", "Step 3 melt column"),
    ("film_mass_kg", "kg", "mass check, §11.3", "Step 3 melt column"),
    ("deep_mass_kg", "kg", "mass check, §11.3", "Step 3 melt column"),
    ("sprayed_mass_kg", "kg (cumulative)", "release_rate check, §11.3", "Step 3 melt column"),
)

_RUN = (
    ("run_name", "run_name", "-", "§4.4 names, provenance", "on main"),
    ("diameter_mm", "inputs.diameter_mm", "mm", "§6.2 body, provenance", "on main"),
    ("initial_mass_kg", "inputs.mass_kg", "kg", "§11.3 accounts", "on main"),
    ("atmosphere", "settings.atmosphere", "-", "air temperature (§11.2), names", "on main"),
    ("wind", "settings.wind", "-", "air temperature (§11.2), names", "on main"),
    ("macro_step_s", "settings.macro_step_s", "s", "§9.1, release_rate per step", "on main (thermal runs)"),
    ("frames_every", "settings.frames_every", "macro steps", "§9.1 frame times", "on main (thermal runs)"),
    ("seed", "settings.seed", "-", "names, provenance", "Step 3 amendment of 2026-10-05"),
    ("material", "settings.material", "-", "§8.1 material table",
     "Step 3's CLI; key name read from the Step 3 plan (confirm)"),
)

FE_FIELDS: tuple[FieldSpec, ...] = (
    # geometry: arrays rather than named fields; the "name" says where the reader finds them
    _f("points", "points", "node", "m, body frame, flight toward +x", True, "§6.2 body, all geometry",
       "on main; vtu and vtp share one node array (confirm)"),
    _f("tets", "TETRA", "tet", "node indices", True, "§6.2-6.3",
       "on main; the active tetrahedra only"),
    _f("faces", "triangles", "patch", "node indices, outward", True, "§6.2, §7, §10",
       "on main; orientation (confirm)"),
    # node (vtu point data)
    _f("T", "T", "node", "K", True, "§6.3, §9.1", "on main"),
    _f("f_l", "liquid_fraction", "node", "-", True, "§6.3, §10",
       "Step 3; equals the material's f_l(T) (confirm)"),
    # tet (vtu cell data)
    _f("phi", "phi", "tet", "-", True, "mass, §6.2", "Step 3; the element's remaining fraction"),
    # patch (vtp cell data)
    _f("q_conv", "q_conv", "patch", "W/m^2", True, "§9.4 heating package", "on main"),
    _f("q_rad", "q_rad", "patch", "W/m^2", False, "informational", "on main"),
    _f("T_patch", "T_patch", "patch", "K", False, "informational", "on main"),
    _f("film_thickness", "film_thickness", "patch", "m", True, "§10 zones", "Step 3"),
    _f("film_T", "film_T", "patch", "K", True, "§14.1 sink enthalpy", "Step 3"),
    _f("p_w", "p_w", "patch", "Pa; 0 = not evaluated", True, "§7.2 loads",
       "Step 3; written on every evaluated step (decision 4, confirm)"),
    _f("tau", "tau", "patch", "Pa, magnitude along the flow-direction tangent", True, "§7.2 loads",
       "Step 3; sign convention (confirm); 0 where p_w is 0"),
    _f("release_rate", "release_rate", "patch", "kg/m^2 per macro step", True, "§10, §14.1 sink",
       "Step 3; sum rate * A = the step's sprayed mass (confirm)"),
    _f("delta_m", "delta_m", "patch", "m; NaN = undefined", True, "§10 zones",
       "Step 3 amendment of 2026-10-02; NaN where the closure is not Girin's", nan_allowed=True),
    _f("deep_thickness", "deep_thickness", "patch", "m (mass per area, m_d / (rho_l A))", False,
       "§10, §2 backlog", "Step 3 amendment of 2026-10-02; read, not physical at 0.5 s"),
    _f("closure", "closure", "patch", "0 Girin / 1 Couette", False, "consistency with delta_m",
       "Step 3; 1 where not evaluated or not Girin's (decision 4, confirm)"),
    _f("r_droplet", "r_droplet", "patch", "m; NaN where nothing released", False, "informational", "Step 3",
       nan_allowed=True),
    _f("we_s", "we_s", "patch", "-", False, "informational", "Step 3"),
    _f("kn_local", "kn_local", "patch", "-", False, "informational", "Step 3; NaN where not evaluated (decision 4, confirm)",
       nan_allowed=True),
    # frame
    _f("frame_time_s", "field.pvd:timestep", "frame", "s, 6 significant digits", True, "§9.1",
       "on main; matched to the history row within 1e-6 s * max(1, t)"),
    # history (<run>.csv)
    *(_f(col, col, "history", units, True, consumer, note) for col, units, consumer, note in _HISTORY),
    # run (<run>.json)
    *(_f(key, name, "run", units, True, consumer, note) for key, name, units, consumer, note in _RUN),
    _f("v_hat", "v_hat", "run", "-, unit vector in the body frame", True, "§7.2 inclination",
       "exported into the run JSON (decision 4); key path and value (confirm)"),
)

FRAME_TIME_RTOL = 1e-6          # frame time matches a history row within FRAME_TIME_RTOL * max(1, t) seconds

_BY_KEY = {f.key: f for f in FE_FIELDS}


def field_spec(key) -> FieldSpec:
    try:
        return _BY_KEY[key]
    except KeyError:
        raise KeyError("no contract item {!r}".format(key)) from None


def fe_name(key) -> str:
    """The finite-element file's name of contract item `key`."""
    return field_spec(key).name


def key_of(location, name) -> str:
    """Inverse of `fe_name` within a location."""
    for f in FE_FIELDS:
        if f.location == location and f.name == name:
            return f.key
    raise KeyError("no contract item named {!r} at {!r}".format(name, location))


def fields(location) -> list[FieldSpec]:
    if location not in LOCATIONS:
        raise ValueError("location must be one of {}, got {!r}".format(LOCATIONS, location))
    return [f for f in FE_FIELDS if f.location == location]


def required(location) -> list[FieldSpec]:
    return [f for f in fields(location) if f.required]


def unconfirmed() -> list[str]:
    return [f.key for f in FE_FIELDS if not f.confirmed]


# the connectivity items are arrays of the Frame, not entries of its node/tet/patch dicts
GEOMETRY_KEYS = ("points", "tets", "faces")


def data_fields(location) -> list[FieldSpec]:
    """The named data fields at a mesh location (node, tet or patch), connectivity excluded."""
    return [f for f in fields(location) if f.key not in GEOMETRY_KEYS]


@dataclass
class Frame:                               # one frame in memory, finite-element numbering
    k: int
    time_s: float
    points: np.ndarray                     # (n_nodes, 3) float64, m
    tets: np.ndarray                       # (n_tets, 4) int64, active tetrahedra
    faces: np.ndarray                      # (n_faces, 3) int64, outward surface triangles
    node: dict[str, np.ndarray] = field(default_factory=dict)    # key -> (n_nodes,)
    tet: dict[str, np.ndarray] = field(default_factory=dict)     # key -> (n_tets,)
    patch: dict[str, np.ndarray] = field(default_factory=dict)   # key -> (n_faces,)


# The compact per-frame file `frame_<k:05d>.npz` (Task 4): published name -> (dtype, shape, units). Shapes use the
# symbols n (nodes kept: those used by an active tetrahedron), m (active tetrahedra), f (surface patches); () is a
# scalar. `patch_<key>` exists for every patch data field of FE_FIELDS that the frame carries.
PREPARED_FRAME_ARRAYS: dict[str, tuple[str, tuple, str]] = {
    "k": ("i8", (), "frame index"),
    "time_s": ("f8", (), "s"),
    "points": ("f8", ("n", 3), "m, body frame"),
    "node_ids": ("i8", ("n",), "finite-element node index of each compact node"),
    "tets": ("i4", ("m", 4), "compact node indices"),
    "faces": ("i4", ("f", 3), "compact node indices, outward"),
    "node_T": ("f8", ("n",), "K"),
    "node_f_l": ("f8", ("n",), "-"),
    "tet_phi": ("f8", ("m",), "-"),
    **{"patch_" + f.key: ("f8", ("f",), f.units) for f in data_fields("patch")},
    "patch_area": ("f8", ("f",), "m^2"),
    "patch_normal": ("f8", ("f", 3), "-, outward unit normal"),
    "patch_centroid": ("f8", ("f", 3), "m"),
    "patch_theta": ("f8", ("f",), "rad, angle between the outward normal and v_hat"),
    "patch_thickness": ("f8", ("f",), "m, inward ray to the opposite surface (Task 5)"),
    "patch_slurry_depth": ("f8", ("f",), "m, contiguous f_l > 0.5 below the surface (Task 6)"),
    "patch_liquid_depth": ("f8", ("f",), "m, contiguous f_l >= 1 below the surface (Task 6)"),
}

# arrays a prepared frame may lack: optional finite-element fields, and the depth arrays under --no-thickness
PREPARED_FRAME_OPTIONAL = frozenset(
    ["patch_" + f.key for f in data_fields("patch") if not f.required]
    + ["patch_thickness", "patch_slurry_depth", "patch_liquid_depth"])
