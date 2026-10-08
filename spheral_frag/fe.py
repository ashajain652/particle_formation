"""The one place `spheral_frag` touches the finite-element package (M1 plan, Global constraints and Task 3).

Prepare side only: every function imports `reentry_model` inside itself, and Spheral's Python never imports this
module's functions. Which `reentry_model`: `import_fe_package(path)` puts `path` first on `sys.path` (default: this
repository), so `prepare --fe-package prototype/proto3` reads Step 3's material until Step 3 lands on main; tests
read that path from `SPHERAL_FRAG_FE_PACKAGE`. One process holds one `reentry_model`: asking for a second one from
another location raises instead of mixing two packages' modules.

`material_arrays` reads the finite-element `Material`'s private node tables (`_T_h`, `_h_nodes`, `_latent_slope`,
`_fl_T`, `_fl`). That is deliberate and confined to this function: they are what makes the core's evaluation
round-off-exact at any temperature (`material.py`)."""
from __future__ import annotations

import hashlib
import importlib
import math
import os
import sys

import numpy as np

from . import material as _material

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MECHANICAL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "aa7075_mechanical.json")
FE_PACKAGE_ENV = "SPHERAL_FRAG_FE_PACKAGE"
_PKG = "reentry_model"


def _under(file, root):
    return os.path.realpath(file).startswith(os.path.realpath(root) + os.sep)


def import_fe_package(path=None):
    """`reentry_model` from `path` (the directory holding the package; default the repository). Raises
    FileNotFoundError when `path` holds no package, RuntimeError when a `reentry_model` from elsewhere is already
    imported in this process."""
    root = os.path.abspath(os.fspath(path)) if path else REPO_ROOT
    if not os.path.isfile(os.path.join(root, _PKG, "__init__.py")):
        raise FileNotFoundError("no {} package in {}".format(_PKG, root))
    mod = sys.modules.get(_PKG)
    if mod is not None:
        if not _under(mod.__file__, root):
            raise RuntimeError("{} is already imported from {}; one process holds one finite-element package"
                               .format(_PKG, os.path.dirname(os.path.dirname(mod.__file__))))
        return mod
    if not sys.path or os.path.abspath(sys.path[0] or os.getcwd()) != root:
        sys.path.insert(0, root)
    mod = importlib.import_module(_PKG)
    if not _under(mod.__file__, root):
        raise RuntimeError("{} resolved to {}, not {}".format(_PKG, mod.__file__, root))
    return mod


def package_provenance(module) -> dict:
    """path (the directory holding the package), the package's __file__ and a SHA-256 over its .py files (sorted by
    relative path; each contributes its path, a NUL, its bytes and a NUL)."""
    pkg_dir = os.path.dirname(os.path.abspath(module.__file__))
    h = hashlib.sha256()
    for dirpath, dirnames, filenames in os.walk(pkg_dir):
        dirnames[:] = sorted(d for d in dirnames if d != "__pycache__")
        for fn in sorted(filenames):
            if fn.endswith(".py"):
                p = os.path.join(dirpath, fn)
                h.update(os.path.relpath(p, pkg_dir).replace(os.sep, "/").encode() + b"\0")
                with open(p, "rb") as fh:
                    h.update(fh.read() + b"\0")
    return {"path": os.path.dirname(pkg_dir), "file": os.path.abspath(module.__file__), "sha256": h.hexdigest()}


def material_file(name, package=None):
    """The DRAMA material JSON `name` (a path, or a name under the package's data/materials/)."""
    if os.path.isfile(os.fspath(name)):
        return os.path.abspath(name)
    pkg = import_fe_package(package)
    p = os.path.join(pkg.DATA_DIR, "materials", "{}.json".format(name))
    if not os.path.isfile(p):
        raise FileNotFoundError("no material {!r} in {}".format(name, os.path.dirname(p)))
    return p


def fe_material(name, package=None):
    """`Material.from_drama_json` of that package for the material `name` (a path or a name)."""
    path = material_file(name, package)
    return importlib.import_module(_PKG + ".material").Material.from_drama_json(path)


def _fe_T_feed(mat):
    try:
        return float(mat.T_feed)
    except AttributeError:
        return math.nan


def material_arrays(mat, mechanical: dict) -> tuple[dict, tuple]:
    """The material table's arrays (material.NODE_ARRAYS, SCALARS, GRID_COLUMNS) for the finite-element material
    `mat`, and the provisional columns. Reads the FE's private node tables (module docstring)."""
    T_h = np.asarray(mat._T_h, dtype=float)
    T_cp, cp_table = np.asarray(mat.T_cp, dtype=float), np.asarray(mat.cp_table, dtype=float)
    slope = getattr(mat, "_latent_slope", None)
    slope = np.zeros(len(T_h) - 1) if slope is None else np.asarray(slope, dtype=float)
    melts = bool(getattr(mat, "melts", mat.latent_heat > 0.0 and np.isfinite(mat.T_liquidus)))
    fl_T_fe = getattr(mat, "_fl_T", None)
    if not melts:
        kind, fl_T, fl = _material.FL_NONE, np.zeros(0), np.zeros(0)
    elif fl_T_fe is not None:
        kind, fl_T, fl = _material.FL_INTERP, np.asarray(fl_T_fe, dtype=float), np.asarray(mat._fl, dtype=float)
    else:
        kind, fl_T, fl = _material.FL_LINEAR, np.array([mat.T_solidus, mat.T_liquidus], float), np.array([0.0, 1.0])
    if kind == _material.FL_INTERP:
        T_half, T_bridge = float(np.interp(0.5, fl, fl_T)), float(np.interp(0.9, fl, fl_T))
    elif kind == _material.FL_LINEAR:
        T_half = float(mat.T_solidus + 0.5 * (mat.T_liquidus - mat.T_solidus))
        T_bridge = float(mat.T_solidus + 0.9 * (mat.T_liquidus - mat.T_solidus))
    else:
        T_half = T_bridge = math.inf
    liq = getattr(mat, "liquid", None)
    rho_l, mu_l, sig_l = (float(liq.rho), float(liq.mu), float(liq.sigma)) if liq is not None else (math.nan,) * 3
    rho_s_json = float(mechanical["entries"]["rho_solid"]["value"])
    if rho_s_json != float(mat.rho):
        raise ValueError("mechanical rho_solid {} differs from the FE density {}".format(rho_s_json, mat.rho))
    a = {
        "h_T": T_h, "h_nodes": np.asarray(mat._h_nodes, dtype=float),
        "h_cp0": np.interp(T_h[:-1], T_cp, cp_table), "h_cp1": np.interp(T_h[1:], T_cp, cp_table),
        "h_slope": slope, "cp_end": np.array([cp_table[0], cp_table[-1]]), "T_cp": T_cp, "cp_table": cp_table,
        "fl_T": fl_T, "fl": fl, "fl_kind": kind,
        "T_ref_fe": float(importlib.import_module(_PKG + ".material").T_REF),
        "h_offset_300": float(mat.enthalpy(_material.T_ZERO)),
        "T_solidus": float(mat.T_solidus), "T_liquidus": float(mat.T_liquidus), "T_feed": _fe_T_feed(mat),
        "T_half": T_half, "T_bridge": T_bridge, "latent_heat": float(mat.latent_heat),
        "rho_solid": float(mat.rho), "rho_liquid": rho_l, "mu_liquid": mu_l, "sigma_liquid": sig_l,
    }
    lo, hi, step = _material.T_GRID
    T_grid = np.arange(lo, hi + 0.5 * step, step)
    f_l = np.asarray(mat.liquid_fraction(T_grid), dtype=float)
    cols, prov = _material.mechanical_columns(T_grid, f_l, a["T_solidus"], T_half, a["rho_solid"], rho_l, mechanical)
    a.update(T_grid=T_grid, h300=np.asarray(mat.enthalpy(T_grid), dtype=float) - a["h_offset_300"], f_l=f_l, **cols)
    return a, prov


def build_material_table(name, out_path, package=None, mechanical_path=DEFAULT_MECHANICAL, command=None):
    """Build and write the material table of the FE material `name` (prepare side). Returns the MaterialTable."""
    pkg = import_fe_package(package)
    mat = fe_material(name, package)
    mechanical = _material.load_mechanical(mechanical_path)
    arrays, prov = material_arrays(mat, mechanical)
    header = {"fe_material": mat.name, "fe_material_file": material_file(name, package),
              "fe_package": package_provenance(pkg), "mechanical_file": os.path.abspath(mechanical_path),
              "mechanical": mechanical, "provisional": list(prov), "command": command,
              "interp_fma": _material.interp_fuses(), "numpy": np.__version__,
              "arrays": {**_material.NODE_ARRAYS, **_material.SCALARS, **_material.GRID_COLUMNS}}
    _material.write_table(out_path, arrays, header)
    return _material.MaterialTable.load(out_path)


def sample_temperatures(table, n_random=100_000, seed=0, lo=250.0, hi=1500.0):
    """The check's temperatures: n_random uniform in [lo, hi] plus every enthalpy and liquid-fraction node and its
    neighbours +-1e-9 K."""
    rng = np.random.default_rng(seed)
    nodes = np.union1d(table.h_T, table.fl_T)
    return np.concatenate([rng.uniform(lo, hi, n_random), nodes, nodes - 1e-9, nodes + 1e-9])


def fe_samples(mat, T) -> dict:
    """The FE material's own values at T (for comparing a table against it in another process)."""
    T = np.asarray(T, dtype=float)
    h300 = float(mat.enthalpy(_material.T_ZERO))
    h_fe = np.asarray(mat.enthalpy(T), dtype=float)
    el = getattr(mat, "enthalpy_liquid", None)
    T_l = float(mat.T_liquidus)
    return {"T": T, "h_fe": h_fe, "h_fe_300": np.float64(h300), "h300": h_fe - h300,
            "f_l": np.asarray(mat.liquid_fraction(T), dtype=float), "cp": np.asarray(mat.cp(T), dtype=float),
            "T_from_h_fe": np.asarray(mat.temperature_from_enthalpy(h_fe), dtype=float),
            "h_full_liquid": np.float64((float(el(T_l)) if el is not None else float(mat.enthalpy(T_l))) - h300)}


def compare_samples(table, s: dict) -> dict:
    """Measured agreement of the core table with the FE values `s` (from fe_samples)."""
    T = s["T"]
    h = table.enthalpy(T)
    dh = np.abs(h - s["h300"])
    fl = table.liquid_fraction(T)
    cp = table.cp(T)
    T_back = table.temperature(h)
    T_inv = table.temperature(s["h300"])
    hfl = table.h_full_liquid()
    return {
        "n": int(len(T)),
        "max_abs_dh_J_per_kg": float(dh.max()), "max_h_J_per_kg": float(np.abs(s["h300"]).max()),
        "h_bitwise": bool(np.array_equal(h, s["h300"])),
        "max_abs_dfl": float(np.abs(fl - s["f_l"]).max()), "fl_bitwise": bool(np.array_equal(fl, s["f_l"])),
        "cp_bitwise": bool(np.array_equal(cp, s["cp"])),
        "max_abs_roundtrip_K": float(np.abs(T_back - T).max()),
        "max_abs_dT_inverse_vs_fe_K": float(np.abs(T_inv - s["T_from_h_fe"]).max()),
        "h_full_liquid_J_per_kg": hfl, "dh_full_liquid": float(abs(hfl - float(s["h_full_liquid"])))
        if math.isfinite(hfl) else (0.0 if hfl == float(s["h_full_liquid"]) else math.inf),
    }


def check_table(mat, table, n_random=100_000, seed=0) -> dict:
    """compare_samples against the FE material in this process."""
    return compare_samples(table, fe_samples(mat, sample_temperatures(table, n_random, seed)))


def air_temperature(run_json: dict, history: dict, package=None) -> np.ndarray:
    """K, the run's air temperature at every history row: the run's atmosphere rebuilt with
    `reentry_model.cli.make_atmosphere` (settings.atmosphere, settings.wind, inputs.epoch_utc) and its `state(t, h,
    lat, lon).T` at the row's time, altitude (clamped at the ground, as the trajectory clamps it) and position."""
    import datetime
    import_fe_package(package)
    cli = importlib.import_module(_PKG + ".cli")
    s, inp = run_json["settings"], run_json["inputs"]
    epoch = datetime.datetime.strptime(inp.get("epoch_utc", "2024-08-01T12:53:07"), "%Y-%m-%dT%H:%M:%S")
    atm = cli.make_atmosphere(s["atmosphere"], epoch, s.get("wind", "none"))[0]
    t, h = np.asarray(history["time_s"], float), np.asarray(history["altitude_km"], float) * 1e3
    lat, lon = np.radians(np.asarray(history["lat_deg"], float)), np.radians(np.asarray(history["lon_deg"], float))
    return np.array([atm.state(float(ti), max(float(hi), 0.0), float(la), float(lo)).T
                     for ti, hi, la, lo in zip(t, h, lat, lon)])


def fe_heat_content(frame, mat) -> float:
    """J: sum over nodes of the lumped mass (sum_e phi_e rho V_e / 4 over the active tetrahedra) times
    h_FE(T_node) - h_FE(300 K). Recorded for M2's check 3; not checked in M1."""
    p, tets = frame.points, frame.tets
    a, b, c, d = (p[tets[:, i]] for i in range(4))
    V = np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)) / 6.0
    phi = frame.tet.get("phi", np.ones(len(tets)))
    m = np.zeros(len(p))
    np.add.at(m, tets.ravel(), np.repeat(phi * float(mat.rho) * V / 4.0, 4))
    h = np.asarray(mat.enthalpy(frame.node["T"]), dtype=float) - float(mat.enthalpy(_material.T_ZERO))
    return math.fsum(m * h)

