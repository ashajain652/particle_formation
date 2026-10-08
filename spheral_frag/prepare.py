"""`python -m spheral_frag prepare`: a finite-element run directory -> Spheral's input files (spec §4.2, §6.1; M1
plan, Task 9).

    "$PY" -m spheral_frag prepare --fe-run <outdir>/<run> [--frames K0:K1] [--every N] [--fe-package PATH]
          [--mechanical spheral_frag/data/aa7075_mechanical.json] [--material-table PATH]
          [--outdir spheral_output/prepare] [--no-thickness] [--force] [--quiet]

Reads the run through the frame contract (`frames.read_fe_run`), builds the material table from the finite-element
material the run names (`fe.build_material_table`, in the finite-element package of `--fe-package`), and for every
selected frame k0, k0 + N, ... <= k1: reads it, checks it (the surface as written, then reoriented outward; mass
against the history; nodal f_l against the table), compacts it, derives the patch geometry, measures the thickness
and the layer depths and classifies the zones (unless `--no-thickness`), and bins its loads. Writes into
`<outdir>/<prepare_name>/`:

  frames/frame_<k:05d>.npz   the compact frame (contract.PREPARED_FRAME_ARRAYS)
  frames/frame_<k:05d>.json  that frame's entry of prepare.json and the settings it was built with (resume sidecar)
  material_table.npz/.json   the material table (material.MaterialTable)
  loads.npz, loads.json      the windward load tables per frame (loads.stack; the runner applies `loads.with_lee`
                             with its own brackets) and their header
  flight.npz                 the history's contract columns plus T_air_K, deceleration_ms2 and drag_N
  prepare.json               last: the summary (PREPARE_SCHEMA below)

**Layer depths (decision 9).** The depths march along the derived surface's outward normals, the patch field
`n_derived` (contract), when the frame carries it; the march is then bounded by the body itself (it ends where the
ray leaves the active tetrahedra), not by the thickness map, which is measured along the facet normal. The facet-
normal depths (capped by the thickness map, as Task 6 built them) are kept beside them as a diagnostic
(`patch_*_depth_facet`). Without `n_derived` the depths are the facet-normal ones and the frame records
`"normals": "facet"`.

**Resume.** An existing `prepare.json` with the same finite-element run JSON (sha256), frame selection and options
is skipped (one line; exit 0, or 1 if it recorded a failed check). Otherwise every frame whose npz exists with the
sha256 and settings its sidecar records is reused, and only the others are rebuilt; `--force` rebuilds everything.
The frame files, the loads and the flight table are byte-identical when rebuilt (fixed zip timestamps).

**Exit codes:** 0 ok; 1 a contract or check failure (prepare.json is still written, `passed: false`, except when
the run itself misses contract items: then no prepared name can be formed and only the message is printed); 2 bad
arguments, a missing run directory, frame, package or material, or a missing pyvista/vtk.

Prepare side only: imports pyvista and vtk (through `frames` and `geometry`) and the finite-element package (through `fe`) inside
its functions; module level is numpy and the standard library."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import time

import numpy as np

from . import __version__, contract, fe, frames, geometry, loads, material, naming, record
from .contract import ContractError

REPO_ROOT = fe.REPO_ROOT
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "spheral_output", "prepare")
PREPARE_SCHEMA = "spheral_frag.prepare"
PREPARE_SCHEMA_VERSION = 1

# Thresholds of the checks prepare gates on (each with its source; Task 10 tightens or confirms them on the flight).
VOLUME_RTOL = 1e-12            # plan Task 4: divergence volume of the oriented surface vs the tetrahedra (round-off)
MASS_RTOL = 5e-9               # plan Task 4 as measured: the history's 9 significant digits (half a unit, digit 1)
MATERIAL_H_RTOL, MATERIAL_H_ATOL = 1e-12, 1e-9    # plan Task 3: |h_core - h_FE| <= 1e-12 max|h| + 1e-9 J/kg
FILM_LIMITS_M = (2e-3, 3e-3)   # spec §10's film limit and its bracket (naming.BRACKET_DEFAULTS film_limit_mm, 3 mm)


# ---------------------------------------------------------------------------------------------------------- helpers
def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def arrays_digest(path) -> str:
    """SHA-256 over an npz's arrays (sorted names, each with its dtype, shape and bytes): equal for equal contents
    even where the file's bytes differ (np.savez stamps the time; the material table is written that way)."""
    h = hashlib.sha256()
    with np.load(path, allow_pickle=False) as z:
        for name in sorted(z.files):
            a = np.ascontiguousarray(z[name])
            h.update("{}\0{}\0{}\0".format(name, a.dtype.str, a.shape).encode() + a.tobytes())
    return h.hexdigest()


def jsonable(x):
    """x with numpy scalars and arrays as Python values and every non-finite float as None (strict JSON)."""
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        x = float(x)
        return x if math.isfinite(x) else None
    return x


def write_json(path, doc) -> None:
    """Strict JSON (NaN and inf as null), written to a temporary file and renamed into place."""
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(jsonable(doc), fh, indent=1, allow_nan=False)
        fh.write("\n")
    os.replace(tmp, path)


def git_commit() -> dict:
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
                                         stderr=subprocess.DEVNULL).decode().strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"],
                                             cwd=REPO_ROOT, stderr=subprocess.DEVNULL).strip())
        return {"commit": commit, "dirty": dirty}
    except Exception:                                 # not a checkout, or no git: provenance says so
        return {"commit": "unknown", "dirty": None}


def directory_size(path) -> int:
    return sum(os.path.getsize(os.path.join(d, f)) for d, _, files in os.walk(path) for f in files)


def parse_frames(text):
    """'K0:K1' (inclusive) or 'K' -> (k0, k1); ValueError otherwise."""
    parts = text.split(":")
    if len(parts) not in (1, 2) or not all(p.strip().isdigit() for p in parts):
        raise ValueError("--frames must be K0:K1 or K with non-negative integers, got {!r}".format(text))
    k0, k1 = int(parts[0]), int(parts[-1])
    if k1 < k0:
        raise ValueError("--frames {} is reversed".format(text))
    return k0, k1


class _TableMaterial:
    """What `fe.fe_heat_content` needs from a material, served by a MaterialTable (whose enthalpy is the finite
    element's less h_FE(300 K), bitwise; enthalpy(300 K) is exactly 0): used with --material-table."""

    def __init__(self, table):
        self.rho = table.rho_solid
        self.enthalpy = table.enthalpy


def say(quiet, *args):
    if not quiet:
        print(*args, flush=True)


# ------------------------------------------------------------------------------------------------------ the frame
def _surface_entry(sc, oriented_faces, oriented_points, n_flipped):
    """The prepare JSON's surface numbers: as written (sc, from surface_checks on the file's winding) and after
    orient_outward (open directed edges, divergence volume)."""
    n_open, _ = frames._edge_counts(oriented_faces)
    nvec, cen = frames.face_vectors(oriented_points, oriented_faces)
    vol_div = math.fsum(np.einsum("ij,ij->i", nvec, cen)) / 6.0
    vol_tets = sc["volume_tets"]
    return {"n_open_directed_edges": int(n_open), "n_nonmanifold_edges": sc["n_nonmanifold_edges"],
            "n_nonmanifold_vertices": sc["n_nonmanifold_vertices"],
            "volume_rel_diff": abs(vol_div / vol_tets - 1.0) if vol_tets > 0 else math.nan,
            "volume_div_m3": vol_div,
            "n_inward_faces_as_written": sc["n_inward_faces"], "n_flipped": int(n_flipped),
            "n_open_directed_edges_as_written": sc["n_open_directed_edges"],
            "n_boundary_edges": sc["n_boundary_edges"], "n_faces_without_owner": sc["n_faces_without_owner"],
            "n_boundary_faces_missing": sc["n_boundary_faces_missing"],
            "faces_on_active_nodes": sc["faces_on_active_nodes"]}


def _zones_entry(cf, liquid, slurry, area):
    film, dm = cf.patch["film_thickness"], cf.patch["delta_m"]
    deep = cf.patch.get("deep_thickness")
    out = {}
    for limit in FILM_LIMITS_M:
        tag = "{:g}mm".format(limit * 1e3)
        out[tag] = geometry.zone_summary(geometry.classify_zones(film, liquid, slurry, dm, film_limit=limit), area)
        if deep is not None:
            zd = geometry.classify_zones(film, liquid, slurry, dm, film_limit=limit, deep_thickness=deep,
                                         include_deep=True)
            out[tag + "_with_deep"] = {"bulk_area_m2": geometry.bulk_slurry_area(zd, area),
                                       "max_layer_m": geometry.zone_summary(zd, area)["max_layer_m"]}
    head = out["{:g}mm".format(FILM_LIMITS_M[0] * 1e3)]
    out["bulk_area_m2_at_2mm"] = head["bulk_area_m2"]
    out["max_layer_m"] = head["max_layer_m"]
    return out


def _depth_info(info, liquid, slurry):
    s = slurry[np.isfinite(slurry)]
    l_ = liquid[np.isfinite(liquid)]
    return {"n_marched": info["n_marched"], "n_wrong_way": info["n_wrong_way"], "n_left_body": info["n_left_body"],
            "max_slurry_depth_m": float(s.max()) if len(s) else math.nan,
            "max_liquid_depth_m": float(l_.max()) if len(l_) else math.nan}


def process_frame(run, k, table, mat, v_hat, out_dir, thickness=True) -> tuple[dict, loads.LoadTable | None]:
    """Read, check, compact, measure and write frame k; returns (its prepare.json entry, its windward load table).
    A frame that misses required contract items or cannot be oriented gets an entry with `missing` or `error` and
    no table (and no npz)."""
    timing = {}
    t0 = time.perf_counter()
    fr = frames.read_frame(run, k)
    timing["read"] = time.perf_counter() - t0
    entry = {"k": int(k), "time_s": fr.time_s, "history_row": run.frame_row(k)}
    miss = frames.missing(fr)
    if miss:
        entry["missing"] = miss
        return entry, None
    hrow = frames.history_row(run, entry["history_row"])

    t = time.perf_counter()
    sc = frames.surface_checks(fr)
    try:
        ofr, flipped = frames.orient_outward(fr)
    except ContractError as e:
        entry["error"] = str(e)
        entry["surface"] = {k_: sc[k_] for k_ in ("n_faces_without_owner", "n_boundary_faces_missing",
                                                  "n_open_directed_edges", "n_inward_faces")}
        return entry, None
    surface = _surface_entry(sc, ofr.faces, ofr.points, len(flipped))
    mc = frames.mass_checks(fr, table.rho_solid, table.rho_liquid, hrow)
    cc = frames.consistency_checks(fr, table)
    cf, node_ids = frames.compact(ofr)
    der = frames.derived_patch_arrays(cf, v_hat)
    heat = fe.fe_heat_content(cf, mat if mat is not None else _TableMaterial(table))
    timing["checks"] = time.perf_counter() - t

    t = time.perf_counter()
    p_w, tau = cf.patch["p_w"], cf.patch["tau"]
    table_k = loads.build_table(der["theta"], der["area"], p_w, tau, loads.history_p_stag(hrow))
    dc = loads.drag_comparison(der["theta"], der["area"], p_w, tau, hrow)
    timing["loads"] = time.perf_counter() - t

    extra = {}
    n_derived = cf.patch.get("n_derived")
    if n_derived is not None:
        norm = np.linalg.norm(n_derived, axis=1)
        entry["n_derived"] = {"max_norm_error": float(np.nanmax(np.abs(norm - 1.0))) if len(norm) else 0.0,
                              "n_against_facet": int(np.count_nonzero(
                                  np.einsum("ij,ij->i", n_derived, der["normal"]) <= 0.0)),
                              "n_nonfinite": int(np.count_nonzero(~np.isfinite(n_derived).all(axis=1)))}
    if thickness:
        t = time.perf_counter()
        thick = geometry.thickness_map(cf.points, cf.faces, der["normal"], der["centroid"])
        timing["thickness"] = time.perf_counter() - t
        t = time.perf_counter()
        locator = geometry.TetLocator.from_frame(cf)
        timing["locator"] = time.perf_counter() - t
        t = time.perf_counter()
        liq_f, slu_f, info_f = geometry.layer_depths(cf, locator, der["normal"], der["centroid"], max_depth=thick,
                                                     return_info=True)
        timing["depths_facet"] = time.perf_counter() - t
        if n_derived is not None:
            t = time.perf_counter()
            liq, slu, info = geometry.layer_depths(cf, locator, n_derived, der["centroid"], return_info=True)
            timing["depths"] = time.perf_counter() - t
            extra.update(slurry_depth_facet=slu_f, liquid_depth_facet=liq_f)
            depths = {"normals": "derived", **_depth_info(info, liq, slu), "facet": _depth_info(info_f, liq_f, slu_f)}
        else:
            liq, slu, info = liq_f, slu_f, info_f
            timing["depths"] = timing.pop("depths_facet")
            depths = {"normals": "facet", **_depth_info(info, liq, slu)}
        extra.update(thickness=thick, slurry_depth=slu, liquid_depth=liq)
        entry["thickness"] = geometry.thickness_summary(thick)
        entry["thickness"]["n_thin_at_2.2mm"] = int(np.count_nonzero(geometry.thin_patches(thick, 2.2e-3)))
        entry["depths"] = depths
        entry["zones"] = _zones_entry(cf, liq, slu, der["area"])

    t = time.perf_counter()
    path = os.path.join(out_dir, "frames", frames.prepared_frame_name(k))
    frames.write_prepared_frame(path, cf, node_ids, {**der, **extra})
    timing["write"] = time.perf_counter() - t
    timing["total"] = time.perf_counter() - t0

    entry.update({
        "file": os.path.join("frames", frames.prepared_frame_name(k)), "sha256": sha256_file(path),
        "n_nodes": len(cf.points), "n_tets": len(cf.tets), "n_faces": len(cf.faces),
        "volume_m3": sc["volume_tets"], "fe_mass_kg": mc["fe_mass"], "mass_rel_diff": mc["mass_rel_diff"],
        "mass": {k_: mc[k_] for k_ in mc if k_ not in ("fe_mass", "mass_rel_diff")},
        "fe_heat_content_J": heat, "surface": surface,
        "consistency": {"max_abs_dfl": cc["max_abs_dfl"], "fl_bitwise": cc["fl_bitwise"],
                        "n_delta_m_mismatch": cc["n_delta_m_mismatch"], "n_girin": cc["n_girin"],
                        "nan_forbidden": cc["nan_forbidden"]},
        "absent_optional": frames.absent_optional(fr),
        "loads_valid": dc["loads_valid"],
        "drag": {"patch_N": dc["D_patch"], "table_N": dc["D_table"], "history_N": dc["D_hist"],
                 "table_lee_N": dc["D_table_lee"], "lee_share": dc["lee_share"], "n_valid": dc["n_valid"],
                 "theta_last_deg": dc["theta_last_deg"], "rel_table_patch": dc["rel_table_patch"],
                 "rel_patch_hist": dc["rel_patch_hist"]},
        "timing_s": timing,
    })
    return entry, table_k


def table_from_prepared(path, run, k) -> loads.LoadTable:
    """The windward load table of a prepared frame file (a reused frame: identical to the one built from the FE
    frame, since the npz holds the same theta, area, p_w and tau)."""
    a = frames.load_prepared_arrays(path)
    hrow = frames.history_row(run, run.frame_row(k))
    return loads.build_table(a["patch_theta"], a["patch_area"], a["patch_p_w"], a["patch_tau"],
                             loads.history_p_stag(hrow))


# ------------------------------------------------------------------------------------------------- flight table
def flight_arrays(run, T_air) -> dict[str, np.ndarray]:
    """The prepared flight table: every history column of the contract under its key, T_air_K, deceleration_ms2
    (= load_factor_g g0) and drag_N (= mass_kg load_factor_g g0) -- the history has neither (plan fact 4)."""
    out = {f.key: np.asarray(run.history[f.name], dtype=np.float64) for f in contract.fields("history")}
    lf = out["load_factor_g"]
    out["T_air_K"] = np.asarray(T_air, dtype=np.float64)
    out["deceleration_ms2"] = lf * loads.G0
    out["drag_N"] = out["mass_kg"] * lf * loads.G0
    missing = [c for c in (record.FLIGHT_TIME,) + record.FLIGHT_STATE if c not in out]
    if missing:                                       # the record's flight state (record.FlightTable) must be here
        raise ContractError("the flight table lacks {}".format(", ".join(missing)))
    return out


# ------------------------------------------------------------------------------------------------------- checks
def _check(value, threshold, source, passed):
    return {"value": value, "threshold": threshold, "source": source, "passed": passed}


def aggregate_checks(entries, material_check, material_name, table_name, missing):
    """The prepare JSON's `checks`: gating checks (passed True/False) with a priori thresholds, and measured ones
    (threshold and passed null) whose thresholds Task 10 sets on the real flight."""
    ok = [e for e in entries if "sha256" in e]
    c = {}
    c["contract_items"] = _check(len(missing), 0, "contract.FE_FIELDS: every required item present", not missing)
    bad = [e["k"] for e in entries if "error" in e]
    c["frames_oriented"] = _check(len(bad), 0, "frames.orient_outward: every face has one active owner", not bad)

    def worst(key, sub=None, fn=max, default=0):
        vals = [(e[sub] if sub else e)[key] for e in ok]
        vals = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
        return fn(vals) if vals else default

    v = worst("n_open_directed_edges", "surface")
    c["open_directed_edges"] = _check(v, 0, "a priori: the active tetrahedra's boundary, wound outward, is closed as an "
                                      "oriented surface (plan fact 7)", v == 0)
    v = max([e["surface"]["n_faces_without_owner"] + e["surface"]["n_boundary_faces_missing"] for e in ok] or [0])
    c["surface_matches_tets"] = _check(v, 0, "a priori: the surface is the active tetrahedra's boundary (plan fact 7)",
                                       v == 0)
    v = worst("volume_rel_diff", "surface", default=0.0)
    c["volume_rel_diff"] = _check(v, VOLUME_RTOL, "plan Task 4: round-off", v <= VOLUME_RTOL)
    v = worst("mass_rel_diff", default=0.0)
    c["mass_rel_diff"] = _check(v, MASS_RTOL, "plan Task 4 as measured: the history's 9 significant digits",
                                v <= MASS_RTOL)
    v = worst("max_abs_dfl", "consistency", default=0.0)
    bitwise = all(e["consistency"]["fl_bitwise"] for e in ok)
    c["f_l_bitwise"] = _check(v, 0.0, "plan Task 3: nodal liquid_fraction equals the table's f_l(T) bitwise", bitwise)
    v = sum(e["consistency"]["n_delta_m_mismatch"] or 0 for e in ok)
    c["delta_m_where_girin"] = _check(v, 0, "contract: delta_m finite exactly where closure is Girin's", v == 0)
    v = sum(sum(e["consistency"]["nan_forbidden"].values()) for e in ok)
    c["nan_forbidden"] = _check(v, 0, "contract: NaN only in fields with nan_allowed", v == 0)
    c["material_name"] = _check(table_name, material_name, "the table is the run's own material (decision 1)",
                                table_name == material_name)
    if material_check is None:
        c["material_table_h"] = _check(None, None, "not run: --material-table given (no finite-element material)", None)
        c["material_table_f_l"] = _check(None, None, "not run: --material-table given", None)
    else:
        thr = MATERIAL_H_RTOL * material_check["max_h_J_per_kg"] + MATERIAL_H_ATOL
        c["material_table_h"] = _check(material_check["max_abs_dh_J_per_kg"], thr, "plan Task 3: 1e-12 max|h| + 1e-9",
                                       material_check["max_abs_dh_J_per_kg"] <= thr)
        c["material_table_f_l"] = _check(material_check["max_abs_dfl"], 0.0, "plan Task 3: bitwise",
                                         material_check["fl_bitwise"])
    measured = "measured; threshold set on the real flight by Task 10"
    with_loads = [e for e in ok if e["loads_valid"]]
    c["frames_without_loads"] = _check(len(ok) - len(with_loads), None, measured + " (plan fact 3)", None)
    c["drag_binning"] = _check(max([abs(e["drag"]["rel_table_patch"]) for e in with_loads
                                    if e["drag"]["rel_table_patch"] is not None
                                    and math.isfinite(e["drag"]["rel_table_patch"])] or [0.0]),
                               None, measured + ": max |D_table / D_patch - 1|", None)
    c["drag_vs_history"] = _check(max([abs(e["drag"]["rel_patch_hist"]) for e in with_loads
                                       if e["drag"]["rel_patch_hist"] is not None
                                       and math.isfinite(e["drag"]["rel_patch_hist"])] or [0.0]),
                                  None, measured + ": max |D_patch / D_hist - 1| (spec §7.3)", None)
    with_depths = [e for e in ok if "depths" in e]
    if with_depths:
        c["thickness_nan"] = _check(sum(e["thickness"]["n_nan"] for e in with_depths), None,
                                    measured + ": inward rays that escape", None)
        c["depth_wrong_way"] = _check(sum(e["depths"]["n_wrong_way"] for e in with_depths), None,
                                      measured + ": depth normals pointing out of the body", None)
        c["depth_left_body"] = _check(sum(e["depths"]["n_left_body"] for e in with_depths), None,
                                      measured + ": slurry layers that reach the body's far side", None)
    return c


# -------------------------------------------------------------------------------------------------------- prepare
def _signature(run_json_sha, csv_sha, pkg, table_sha, thickness, v_hat):
    """What a prepared frame depends on beyond its own FE files (resume: a sidecar must match it)."""
    return {"fe_json_sha256": run_json_sha, "fe_csv_sha256": csv_sha, "fe_package_sha256": pkg["sha256"],
            "material_table_arrays_sha256": table_sha, "thickness": bool(thickness), "v_hat": list(map(float, v_hat)),
            "contract_version": contract.CONTRACT_VERSION, "spheral_frag_version": __version__,
            "depth_ds_m": geometry.DEPTH_DS_M}


def _material_table(args, out, run, material_name, command):
    """(MaterialTable, FE material or None, the round-off check or None, its source). Raises FileNotFoundError
    (exit 2) for a material the package or the given file does not hold."""
    dst = os.path.join(out, "material_table.npz")
    if args.material_table:
        src = os.path.abspath(args.material_table)
        if os.path.isdir(src):
            src = os.path.join(src, "material_table.npz")
        if not (os.path.isfile(src) and os.path.isfile(material.header_path(src))):
            raise FileNotFoundError("no material table {} (with its JSON header)".format(src))
        if os.path.abspath(src) != os.path.abspath(dst):
            shutil.copyfile(src, dst)
            shutil.copyfile(material.header_path(src), material.header_path(dst))
        return material.MaterialTable.load(dst), None, None, {"from_file": src}
    mat = fe.fe_material(material_name, args.fe_package)
    table = fe.build_material_table(material_name, dst, args.fe_package, args.mechanical, command)
    return table, mat, fe.check_table(mat, table), {"from_fe_package": True}


def prepare(args) -> int:
    t_start = time.perf_counter()
    quiet = args.quiet
    try:
        import pyvista  # noqa: F401
        import vtkmodules.vtkCommonDataModel  # noqa: F401
    except ImportError as e:
        print("prepare needs pyvista and vtk: {}".format(e), file=sys.stderr)
        return 2
    if args.every < 1:
        print("--every must be a positive integer", file=sys.stderr)
        return 2
    try:
        sel_range = parse_frames(args.frames) if args.frames else None
    except ValueError as e:
        print(e, file=sys.stderr)
        return 2
    try:
        run = frames.read_fe_run(args.fe_run)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        return 2
    except ContractError as e:
        print("contract failure: {}".format(e), file=sys.stderr)
        return 1

    ks_all = [k for k, _, _ in run.frames]
    k0, k1 = sel_range if sel_range else (ks_all[0], ks_all[-1])
    absent = [k for k in (k0, k1) if k not in run.frame_files]
    if absent:
        print("--frames: run {} has no frame {} (frames {}-{})".format(run.name, absent, ks_all[0], ks_all[-1]),
              file=sys.stderr)
        return 2
    selected = [k for k in ks_all if k0 <= k <= k1 and (k - k0) % args.every == 0]
    name = naming.prepare_name(run.name, k0, k1, args.every)
    out = os.path.join(os.path.abspath(args.outdir), name)
    json_path = os.path.join(out, "prepare.json")
    json_sha, csv_sha = sha256_file(run.json_path), sha256_file(run.csv_path)
    options = {"thickness": not args.no_thickness, "material_table": os.path.abspath(args.material_table)
               if args.material_table else None, "mechanical": os.path.abspath(args.mechanical),
               "fe_package_arg": os.path.abspath(args.fe_package) if args.fe_package else None}

    if os.path.isfile(json_path) and not args.force:
        try:
            with open(json_path, encoding="utf-8") as fh:
                old = json.load(fh)
        except (OSError, ValueError):
            old = None
        if old and old.get("fe_run", {}).get("json_sha256") == json_sha and \
                [f["k"] for f in old.get("frames", [])] == selected and old.get("options") == jsonable(options):
            print("prepared already ({}): {}{}".format("passed" if old.get("passed") else "FAILED checks", out,
                                                      "" if old.get("passed") else "; --force to rebuild"))
            return 0 if old.get("passed") else 1

    try:
        pkg = fe.import_fe_package(args.fe_package)
    except (FileNotFoundError, RuntimeError) as e:
        print("--fe-package: {}".format(e), file=sys.stderr)
        return 2
    pkg_info = fe.package_provenance(pkg)
    try:
        material_name = frames.run_value(run, "material")
        v_hat = frames.run_v_hat(run)
    except ContractError as e:
        print("contract failure: {}".format(e), file=sys.stderr)
        return 1
    os.makedirs(os.path.join(out, "frames"), exist_ok=True)
    command = " ".join(["python -m spheral_frag"] + list(getattr(args, "argv", None) or ["prepare"]))
    try:
        table, mat, material_check, table_source = _material_table(args, out, run, material_name, command)
    except FileNotFoundError as e:
        print("material: {} (the run's material is {!r}; pass the --fe-package that produced the run, or "
              "--material-table)".format(e, material_name), file=sys.stderr)
        return 2
    say(quiet, "prepare {} -> {}".format(run.name, out))
    say(quiet, "  material table {} ({}); provisional columns: {}".format(
        table.name, "built from the FE package" if mat is not None else "given", ", ".join(table.provisional)
        or "none"))

    T_air = fe.air_temperature(run.run_json, run.history, args.fe_package)
    flight = flight_arrays(run, T_air)
    frames.write_npz(os.path.join(out, "flight.npz"), flight)

    signature = _signature(json_sha, csv_sha, pkg_info, arrays_digest(os.path.join(out, "material_table.npz")),
                           options["thickness"], v_hat)
    entries, tables, reused = [], [], 0
    missing = []
    for i, k in enumerate(selected):
        npz = os.path.join(out, "frames", frames.prepared_frame_name(k))
        side = os.path.splitext(npz)[0] + ".json"
        entry, was_reused = None, False
        if not args.force and os.path.isfile(npz) and os.path.isfile(side):
            try:
                with open(side, encoding="utf-8") as fh:
                    doc = json.load(fh)
                if doc.get("signature") == jsonable(signature) and doc["entry"]["sha256"] == sha256_file(npz):
                    entry = doc["entry"]
                    table_k = table_from_prepared(npz, run, k)
                    reused += 1
                    was_reused = True
            except (OSError, ValueError, KeyError):
                entry = None
        if entry is None:
            entry, table_k = process_frame(run, k, table, mat, v_hat, out, thickness=options["thickness"])
            if "sha256" in entry:
                write_json(side, {"signature": signature, "entry": entry})
        entries.append(jsonable(entry))
        if "missing" in entry:
            missing = sorted(set(missing) | set(entry["missing"]))
            print("frame {} misses required contract items: {}".format(k, ", ".join(entry["missing"])),
                  file=sys.stderr)
            if i == 0:
                break                                  # the first frame already fails the contract: stop there
            continue
        if "error" in entry:
            print("frame {}: {}".format(k, entry["error"]), file=sys.stderr)
            continue
        tables.append((k, entry["time_s"], table_k))
        if not quiet and (i % 50 == 0 or i == len(selected) - 1):
            tt = entry.get("timing_s", {})
            say(quiet, "  frame {:5d}  t {:8.2f} s  faces {:6d}  {}".format(
                k, entry["time_s"], entry["n_faces"],
                "reused" if was_reused else "{:.2f} s (thickness {:.2f}, depths {:.2f})".format(
                    tt.get("total", 0.0), tt.get("thickness", 0.0), tt.get("depths", 0.0))))

    if tables:
        arrays = loads.stack([t for _, _, t in tables])
        arrays["k"] = np.array([k for k, _, _ in tables], dtype=np.int64)
        arrays["time_s"] = np.array([t for _, t, _ in tables], dtype=np.float64)
        frames.write_npz(os.path.join(out, "loads.npz"), arrays)
        header = loads.table_header(tables[0][2])
        header["lee"] = ("not applied: the tables are windward, as the frames carry them; the runner applies "
                         "loads.with_lee with its run's brackets (base_pressure_pct, separation_deg, lee_shear)")
        header["arrays"] = {"k": "frame index", "time_s": "s", "p": "Pa (n_frames, n_bins)",
                            "tau": "Pa (n_frames, n_bins)", "area": "m^2 (n_frames, n_bins)",
                            "lee": "int8, bins set by the lee extension (none here)",
                            "theta_last_deg": "deg, NaN on a frame without loads", "p_stag": "Pa",
                            "theta_deg": "bin centres", "edges_deg": "bin edges"}
        write_json(os.path.join(out, "loads.json"), header)

    checks = aggregate_checks(entries, material_check, material_name, table.name, missing)
    passed = all(c["passed"] is not False for c in checks.values())
    first = next((e for e in entries if "sha256" in e), None)
    present = [f.key for f in contract.FE_FIELDS if f.location in ("history", "run", "frame") or
               (first is not None and f.key not in first["absent_optional"] and f.key not in missing)]
    absent_opt = sorted({key for e in entries for key in e.get("absent_optional", [])})
    depth_kinds = sorted({e["depths"]["normals"] for e in entries if "depths" in e})
    per_frame = [e["timing_s"] for e in entries if "timing_s" in e]

    def median(key):
        vals = [t[key] for t in per_frame if key in t]
        return statistics.median(vals) if vals else None

    doc = {
        "schema": PREPARE_SCHEMA, "schema_version": PREPARE_SCHEMA_VERSION,
        "contract_version": contract.CONTRACT_VERSION, "name": name,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "spheral_frag_version": __version__, "git_commit": git_commit()["commit"], "git": git_commit(),
        "command": command, "options": options,
        "fe_run": {"name": run.name, "dir": run.directory, "json": run.json_path, "csv": run.csv_path,
                   "json_sha256": json_sha, "csv_sha256": csv_sha, "package": pkg_info,
                   "settings": run.run_json.get("settings", {})},
        "flight": {"v_hat": v_hat, "diameter_m": float(frames.run_value(run, "diameter_mm")) * 1e-3,
                   "initial_mass_kg": float(frames.run_value(run, "initial_mass_kg")), "material": material_name,
                   "macro_step_s": float(frames.run_value(run, "macro_step_s")),
                   "frames_every": int(frames.run_value(run, "frames_every")), "seed": frames.run_value(run, "seed"),
                   "file": "flight.npz", "columns": sorted(flight), "n_rows": int(len(flight["time_s"]))},
        "selection": {"k0": k0, "k1": k1, "every": args.every, "n_frames": len(selected)},
        "contract": {"present": present, "absent_optional": absent_opt, "unconfirmed": contract.unconfirmed(),
                     "missing": missing},
        "material_table": {"file": "material_table.npz", "fe_material": table.name, "run_material": material_name,
                           "provisional": list(table.provisional), "source": table_source,
                           "npz_sha256": table.header.get("npz_sha256"),
                           "max_abs_dh_J_per_kg": material_check["max_abs_dh_J_per_kg"] if material_check else None,
                           "fl_bitwise": material_check["fl_bitwise"] if material_check else None,
                           "check": material_check},
        "loads": {"file": "loads.npz" if tables else None, "header": "loads.json" if tables else None,
                  "n_frames": len(tables), "lee_applied": False},
        "depth_normals": depth_kinds[0] if len(depth_kinds) == 1 else ("mixed" if depth_kinds else None),
        "frames": entries,
        "checks": checks, "passed": passed,
        "timing_s": {"total": time.perf_counter() - t_start, "per_frame_median": median("total"),
                     "read_median": median("read"), "checks_median": median("checks"),
                     "thickness_median": median("thickness"), "locator_median": median("locator"),
                     "depths_median": median("depths"),
                     "depths_facet_median": median("depths_facet"), "write_median": median("write"),
                     "n_frames_built": len(per_frame) - reused, "n_frames_reused": reused},
        "size_bytes": directory_size(out),
    }
    write_json(json_path, doc)
    failed = [n for n, c in checks.items() if c["passed"] is False]
    say(quiet, "  {} frames ({} reused), {:.1f} s, {:.1f} MB; checks {}".format(
        len(entries), reused, doc["timing_s"]["total"], doc["size_bytes"] / 1e6,
        "passed" if passed else "FAILED: " + ", ".join(failed)))
    if not passed:
        for n in failed:
            c = checks[n]
            print("check {} failed: value {} against {} ({})".format(n, c["value"], c["threshold"], c["source"]),
                  file=sys.stderr)
        return 1
    return 0


def add_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument("--fe-run", required=True, help="finite-element run: <outdir>/<run>, <outdir>/<run>.json, or "
                   "an <outdir> holding one run")
    p.add_argument("--frames", default=None, help="K0:K1, inclusive (default: every frame of the run)")
    p.add_argument("--every", type=int, default=1, help="every N-th frame from K0 (default 1)")
    p.add_argument("--fe-package", default=None, help="directory holding the finite-element package that produced the run "
                   "(default: this repository's); e.g. prototype/proto3 until Step 3 lands")
    p.add_argument("--mechanical", default=fe.DEFAULT_MECHANICAL, help="mechanical input file of the material table")
    p.add_argument("--material-table", default=None, help="use this material table (npz with its JSON header) "
                   "instead of building one from the finite-element material (the synthetic fixtures' material)")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR, help="parent of the prepared directory (default "
                   "spheral_output/prepare)")
    p.add_argument("--no-thickness", action="store_true", help="skip the thickness map, the depths and the zones")
    p.add_argument("--force", action="store_true", help="rebuild every file even if prepared already")
    p.add_argument("--quiet", action="store_true", help="print only the summary and errors")
