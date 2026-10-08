"""Frame import and the compact prepared frame (spec §6.1-6.2, §12 check "Frame import"; M1 plan, Task 4).

Three parts:

1. **Reading** a finite-element run directory through the frame contract (`contract.FE_FIELDS`): `read_fe_run`
   (the run JSON, the history CSV and the pvd collections, each frame matched to its history row) and `read_frame`
   (one frame's vtu and vtp, in finite-element numbering, every array exactly as the files hold it). pyvista is
   imported inside `read_frame` (prepare side only).
2. **Checks** that return numbers: `surface_checks` (closure as an oriented surface, non-manifold edges and
   vertices, the surface against the active tetrahedra, divergence volume, orientation face by face),
   `mass_checks` (sum phi rho V + film + deep against the history), `consistency_checks` (nodal f_l against the
   material table, delta_m against the closure, NaN where the contract forbids it).
3. **The compact frame**: `orient_outward` (each surface triangle wound away from its owning tetrahedron),
   `compact` (unreferenced nodes dropped, renumbered), `derived_patch_arrays`, `write_prepared_frame` and the
   runner's numpy-only reader `load_prepared_frame`.

Orientation (decision taken here). `read_frame` keeps the file's winding, so that the import can be compared with
the files bitwise and `surface_checks` can report what the export wrote (`n_inward_faces`). The prepared frame is
always wound outward: `orient_outward` swaps two vertices of every face whose right-hand normal points into its
owning tetrahedron. The fix is exact, not a heuristic -- every surface face has exactly one active owner, and the
owner's opposite vertex fixes the inside -- and it changes nothing the finite element computed (patch fields are per
face; Step 3's `VolumeMesh.surface()` flips the normals it uses but not the node order, so the solver's normals were
outward and only the stored winding is off). Flagging alone would make `prepare` fail a whole flight on a storage
artefact; reorienting silently would hide it. So `prepare` does both: it records the count as written and prepares
the reoriented frame. A face without an active owner cannot be oriented and is a contract failure.

Numpy and the standard library at module level (the runner imports `load_prepared_frame` under Spheral's Python)."""
from __future__ import annotations

import csv
import json
import math
import os
import re
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import dataclass, field, replace

import numpy as np

from . import contract
from .contract import ContractError, Frame, V_HAT, fe_name

PVD_FIELD = "field.pvd"
PVD_SURFACE = "surface.pvd"
_FRAME_FILE = {"field": re.compile(r"^field_(\d+)\.vtu$"), "surface": re.compile(r"^surface_(\d+)\.vtp$")}
_NPZ_DATE = (1980, 1, 1, 0, 0, 0)            # fixed zip timestamps: a rebuilt prepared frame is byte-identical

# the four faces of a tetrahedron (a, b, c, d) and the local index of each face's opposite vertex
_TET_FACES = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]], dtype=np.int64)


# ---------------------------------------------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class FERun:
    name: str
    directory: str                            # <outdir>/<run> (holds vtk/)
    run_json: dict
    history: dict[str, np.ndarray]            # every CSV column (the contract's and any others), float64
    frames: list[tuple[int, float, int]]      # (k, time from the pvd, matching history row), ascending k
    csv_path: str = ""
    json_path: str = ""
    frame_files: dict[int, tuple[str, str]] = field(default_factory=dict)   # k -> (vtu path, vtp path)

    def frame_row(self, k) -> int:
        for kk, _, row in self.frames:
            if kk == k:
                return row
        raise KeyError("run {} has no frame {}".format(self.name, k))

    def frame_time(self, k) -> float:
        for kk, t, _ in self.frames:
            if kk == k:
                return t
        raise KeyError("run {} has no frame {}".format(self.name, k))


def _resolve_run(path):
    """(outdir, run name) from `<outdir>/<run>`, `<outdir>/<run>.json`, or an `<outdir>` holding exactly one run."""
    path = os.path.abspath(os.fspath(path)).rstrip(os.sep)
    if path.endswith(".json"):
        if not os.path.isfile(path):
            raise FileNotFoundError("no run JSON {}".format(path))
        return os.path.dirname(path), os.path.basename(path)[:-5]
    if not os.path.isdir(path):
        raise FileNotFoundError("no finite-element run directory {}".format(path))
    if os.path.isdir(os.path.join(path, "vtk")):
        return os.path.dirname(path), os.path.basename(path)
    runs = sorted(fn[:-5] for fn in os.listdir(path)
                  if fn.endswith(".json") and os.path.isdir(os.path.join(path, fn[:-5], "vtk")))
    if len(runs) == 1:
        return path, runs[0]
    raise FileNotFoundError("{} is neither a run directory (<outdir>/<run> with vtk/) nor an output directory "
                            "holding exactly one run (found {})".format(path, runs or "none"))


def json_value(doc, dotted):
    """The value at a dotted path of the run JSON (KeyError if absent)."""
    for part in dotted.split("."):
        if not isinstance(doc, dict) or part not in doc:
            raise KeyError(dotted)
        doc = doc[part]
    return doc


def run_value(run: FERun, key):
    """Contract run item `key` (e.g. "seed", "material", "v_hat") of the run."""
    return json_value(run.run_json, fe_name(key))


def run_v_hat(run: FERun) -> np.ndarray:
    v = np.asarray(run_value(run, "v_hat"), dtype=np.float64)
    if v.shape != (3,) or not np.isfinite(v).all() or abs(np.linalg.norm(v) - 1.0) > 1e-12:
        raise ContractError("v_hat must be a unit 3-vector, got {}".format(v.tolist()))
    return v


def read_history(csv_path) -> dict[str, np.ndarray]:
    with open(csv_path, newline="") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        raise ContractError("empty history {}".format(csv_path))
    header, body = rows[0], [r for r in rows[1:] if r]
    if len(set(header)) != len(header):
        raise ContractError("duplicate history columns in {}".format(csv_path))
    try:
        data = np.array([[float(v) for v in r] for r in body], dtype=np.float64).reshape(len(body), len(header))
    except ValueError as e:
        raise ContractError("history {} is not numeric: {}".format(csv_path, e)) from None
    return {name: np.ascontiguousarray(data[:, i]) for i, name in enumerate(header)}


def read_pvd(path) -> list[tuple[str, str]]:
    """(timestep as printed, file) per DataSet of a pvd collection, in file order."""
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as e:
        raise ContractError("unreadable pvd {}: {}".format(path, e)) from None
    return [(d.get("timestep"), d.get("file")) for d in root.iter("DataSet")]


def match_row(times, t, rtol=contract.FRAME_TIME_RTOL) -> int:
    """The history row whose time matches `t` within rtol * max(1, |t|) s; ContractError if none."""
    tol = rtol * max(1.0, abs(t))
    d = np.abs(np.asarray(times) - t)
    i = int(np.argmin(d)) if len(d) else -1
    if i < 0 or not d[i] <= tol:
        raise ContractError("no history row at t = {!r} s (within {:.3g} s); nearest {}".format(
            t, tol, "none" if i < 0 else "{!r} s".format(float(times[i]))))
    if np.count_nonzero(d <= tol) > 1:
        raise ContractError("several history rows match t = {!r} s".format(t))
    return i


def read_fe_run(path) -> FERun:
    """The run at `path` (`<outdir>/<run>`, `<outdir>/<run>.json`, or an `<outdir>` with one run).

    Raises FileNotFoundError for a missing directory or file (exit 2), ContractError for a contract miss: a history
    column or run item of the contract absent, the two pvd collections disagreeing, a listed frame file absent, or a
    frame time without a history row (exit 1)."""
    outdir, name = _resolve_run(path)
    json_path = os.path.join(outdir, name + ".json")
    csv_path = os.path.join(outdir, name + ".csv")
    run_dir = os.path.join(outdir, name)
    vtk_dir = os.path.join(run_dir, "vtk")
    for p in (json_path, csv_path, os.path.join(vtk_dir, PVD_FIELD), os.path.join(vtk_dir, PVD_SURFACE)):
        if not os.path.isfile(p):
            raise FileNotFoundError("finite-element run {}: missing {}".format(name, p))
    with open(json_path) as fh:
        doc = json.load(fh)
    history = read_history(csv_path)

    problems = ["history column {} ({})".format(f.name, f.key) for f in contract.fields("history")
                if f.name not in history]
    for f in contract.fields("run"):
        try:
            json_value(doc, f.name)
        except KeyError:
            if f.required:
                problems.append("run JSON item {} ({})".format(f.name, f.key))
    if problems:
        raise ContractError("finite-element run {} misses contract items: {}".format(name, "; ".join(problems)))
    if json_value(doc, fe_name("run_name")) != name:
        raise ContractError("run JSON {} names the run {!r}".format(json_path, json_value(doc, fe_name("run_name"))))

    entries = {}
    for kind, pvd in (("field", PVD_FIELD), ("surface", PVD_SURFACE)):
        ks = {}
        for ts, fn in read_pvd(os.path.join(vtk_dir, pvd)):
            m = _FRAME_FILE[kind].match(fn or "")
            if m is None or ts is None:
                raise ContractError("{}: unexpected DataSet timestep={!r} file={!r}".format(pvd, ts, fn))
            k = int(m.group(1))
            if k in ks:
                raise ContractError("{}: frame {} listed twice".format(pvd, k))
            ks[k] = (ts, os.path.join(vtk_dir, fn))
        entries[kind] = ks
    if set(entries["field"]) != set(entries["surface"]):
        raise ContractError("field.pvd and surface.pvd list different frames: {} vs {}".format(
            sorted(entries["field"]), sorted(entries["surface"])))
    frames, files = [], {}
    times = history[fe_name("time_s")]
    for k in sorted(entries["field"]):
        (ts, vtu), (ts2, vtp) = entries["field"][k], entries["surface"][k]
        if ts != ts2:
            raise ContractError("frame {}: field.pvd time {} but surface.pvd time {}".format(k, ts, ts2))
        for p in (vtu, vtp):
            if not os.path.isfile(p):
                raise ContractError("frame {} is listed but {} is missing".format(k, p))
        t = float(ts)
        frames.append((k, t, match_row(times, t)))
        files[k] = (vtu, vtp)
    if not frames:
        raise ContractError("finite-element run {} has no frames".format(name))
    return FERun(name=name, directory=run_dir, run_json=doc, history=history, frames=frames, csv_path=csv_path,
                 json_path=json_path, frame_files=files)


def history_row(run: FERun, row) -> dict[str, float]:
    return {col: float(v[row]) for col, v in run.history.items()}


def read_frame(run: FERun, k) -> Frame:
    """Frame k in finite-element numbering: every node of the file (dead ones included), the active tetrahedra, the
    surface triangles in the file's winding, and every contract data field the files carry, as the files hold them
    (prepare side only: pyvista). Absent fields are simply absent; `missing(frame)` lists the required ones."""
    import pyvista as pv
    vtu, vtp = run.frame_files[k]
    grid, poly = pv.read(vtu), pv.read(vtp)
    points = np.asarray(grid.points)
    if points.dtype != np.float64:
        points = points.astype(np.float64)                    # exact for float32
    if np.asarray(poly.points).shape != points.shape or \
            np.asarray(poly.points, dtype=np.float64).tobytes() != points.tobytes():
        raise ContractError("frame {}: the vtu and the vtp do not share one node array".format(k))
    celltypes = np.asarray(grid.celltypes)
    if len(celltypes) and not (celltypes == pv.CellType.TETRA).all():
        raise ContractError("frame {}: the vtu holds cells other than tetrahedra: {}".format(
            k, sorted(set(int(c) for c in celltypes) - {int(pv.CellType.TETRA)})))
    tets = np.asarray(grid.cells_dict.get(pv.CellType.TETRA, np.zeros((0, 4))), dtype=np.int64).reshape(-1, 4)
    if poly.n_lines or poly.n_verts or poly.n_strips:
        raise ContractError("frame {}: the vtp holds cells other than polygons".format(k))
    flat = np.asarray(poly.faces, dtype=np.int64)
    if flat.size != 4 * poly.n_cells or (flat.size and not (flat.reshape(-1, 4)[:, 0] == 3).all()):
        raise ContractError("frame {}: the vtp's polygons are not all triangles".format(k))
    faces = flat.reshape(-1, 4)[:, 1:].copy()

    def take(data, location, n):
        out = {}
        for f in contract.data_fields(location):
            if f.name in data:
                a = np.asarray(data[f.name])
                if a.shape != (n,):
                    raise ContractError("frame {}: {} {} has shape {}, expected ({},)".format(k, location, f.name,
                                                                                             a.shape, n))
                out[f.key] = a if a.dtype == np.float64 else a.astype(np.float64)
        return out

    return Frame(k=k, time_s=run.frame_time(k), points=points, tets=tets, faces=faces,
                 node=take(grid.point_data, "node", len(points)), tet=take(grid.cell_data, "tet", len(tets)),
                 patch=take(poly.cell_data, "patch", len(faces)))


def _present(frame, location):
    return {"node": frame.node, "tet": frame.tet, "patch": frame.patch}[location]


def missing(frame) -> list[str]:
    """Keys of the required node, tet and patch contract items the frame lacks (empty connectivity counts)."""
    out = [key for key, arr in (("points", frame.points), ("tets", frame.tets), ("faces", frame.faces))
           if arr is None or len(arr) == 0]
    for loc in ("node", "tet", "patch"):
        out += [f.key for f in contract.data_fields(loc) if f.required and f.key not in _present(frame, loc)]
    return out


def absent_optional(frame) -> list[str]:
    return [f.key for loc in ("node", "tet", "patch") for f in contract.data_fields(loc)
            if not f.required and f.key not in _present(frame, loc)]


# ---------------------------------------------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------------------------------------------

def tet_volumes(points, tets) -> np.ndarray:
    """Unsigned volumes of the tetrahedra, m^3."""
    x = points[tets]
    return np.abs(np.einsum("ij,ij->i", x[:, 1] - x[:, 0], np.cross(x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]))) / 6.0


def face_vectors(points, faces):
    """(right-hand normal times twice the area, centroid) per face."""
    a, b, c = (points[faces[:, i]] for i in range(3))
    return np.cross(b - a, c - a), (a + b + c) / 3.0


def _triple_keys(tri, n):
    s = np.sort(np.asarray(tri, dtype=np.int64), axis=1)
    return (s[:, 0] * n + s[:, 1]) * n + s[:, 2]


def surface_owners(points, tets, faces):
    """For each face: its owning active tetrahedron (-1 if the face is no face of exactly one active tetrahedron)
    and the owner's opposite vertex (-1 likewise); plus how many faces used by exactly one active tetrahedron (the
    true boundary) are absent from `faces`, and how many `faces` are used by two active tetrahedra (interior)."""
    n = max(len(points), 1)
    owner = np.full(len(faces), -1, np.int64)
    opposite = np.full(len(faces), -1, np.int64)
    if len(tets) == 0:
        return owner, opposite, 0, 0
    keys = _triple_keys(tets[:, _TET_FACES].reshape(-1, 3), n)  # face j of tet i at slot 4 i + j, opposite vertex j
    order = np.argsort(keys, kind="stable")
    uniq, first, counts = np.unique(keys[order], return_index=True, return_counts=True)
    fkeys = _triple_keys(faces, n) if len(faces) else np.zeros(0, np.int64)
    pos = np.minimum(np.searchsorted(uniq, fkeys), len(uniq) - 1)
    found = uniq[pos] == fkeys
    once = found & (counts[pos] == 1)
    slot = order[first[pos[once]]]
    owner[once] = slot // 4
    opposite[once] = tets[slot // 4, slot % 4]
    n_interior = int(np.count_nonzero(found & ~once))
    n_boundary_missing = int(np.count_nonzero(~np.isin(uniq[counts == 1], fkeys)))
    return owner, opposite, n_boundary_missing, n_interior


def _edge_counts(faces):
    """(directed edge instances without a matching reverse instance, face count of every undirected edge)."""
    if len(faces) == 0:
        return 0, np.zeros(0, np.int64)
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]).astype(np.int64)
    n = int(e.max()) + 1
    uf, cf = np.unique(e[:, 0] * n + e[:, 1], return_counts=True)
    ur, cr = np.unique(e[:, 1] * n + e[:, 0], return_counts=True)   # the reverses
    pos = np.minimum(np.searchsorted(ur, uf), len(ur) - 1)
    m_rev = np.where(ur[pos] == uf, cr[pos], 0)
    und = np.sort(e, axis=1)
    _, ucount = np.unique(und[:, 0] * n + und[:, 1], return_counts=True)
    return int(np.maximum(cf - m_rev, 0).sum()), ucount


def _nonmanifold_vertices(faces):
    """Vertices whose incident triangles do not form one edge-connected fan (a bowtie: two pieces of surface that
    touch only at the vertex). Returns their node indices."""
    nf = len(faces)
    if nf == 0:
        return np.zeros(0, np.int64)
    corner_v = faces.ravel().astype(np.int64)                    # corner 3 f + i = (vertex faces[f, i], face f)
    n = int(corner_v.max()) + 1
    # two corners at the same vertex are linked when their faces share an edge through that vertex
    pairs = []
    for i in range(3):
        for j in (1, 2):
            a = 3 * np.arange(nf) + i                            # corner at faces[:, i]
            b = 3 * np.arange(nf) + (i + j) % 3                  # the edge (faces[:, i], faces[:, (i+j)%3])
            pairs.append(np.stack([a, b], axis=1))
    half = np.concatenate(pairs)                              # (corner at v, corner at w) along edge v-w of a face
    v = corner_v[half[:, 0]]
    w = corner_v[half[:, 1]]
    key = v * n + w                                              # directed by the corner's own vertex
    order = np.argsort(key, kind="stable")
    key, half = key[order], half[order]
    # consecutive entries with the same (v, w) key are corners at v of faces sharing edge v-w
    same = np.flatnonzero(key[1:] == key[:-1])
    pi = half[same, 0]
    pj = half[same + 1, 0]
    label = np.arange(3 * nf)
    while True:
        m = np.minimum(label[pi], label[pj])
        new = label.copy()
        np.minimum.at(new, pi, m)
        np.minimum.at(new, pj, m)
        new = new[new]                                           # pointer jumping
        if np.array_equal(new, label):
            break
        label = new
    comps = np.unique(np.stack([corner_v, label], axis=1), axis=0)
    verts, nfans = np.unique(comps[:, 0], return_counts=True)
    return verts[nfans > 1]


# ---------------------------------------------------------------------------------------------------------------
# Checks (each returns numbers)
# ---------------------------------------------------------------------------------------------------------------

def surface_checks(frame) -> dict:
    """The surface of `frame` as given (in the file's winding when it comes from `read_frame`):

    n_open_directed_edges   directed edge instances without a matching reverse instance; 0 = closed as an
                            oriented surface. A removed triangle gives 3; an inward-wound triangle among outward
                            neighbours gives 6 (its three edges duplicate the neighbours' and lack reverses)
    n_boundary_edges        undirected edges in exactly one triangle
    n_nonmanifold_edges     undirected edges in more than two triangles (plan fact 7)
    n_nonmanifold_vertices  vertices whose triangles form more than one edge-connected fan (plan fact 7)
    faces_on_active_nodes   every face vertex is used by an active tetrahedron (n_faces_off_active_nodes)
    n_faces_without_owner   faces that are no face of exactly one active tetrahedron
    n_boundary_faces_missing  faces used by exactly one active tetrahedron that the surface lacks
    surface_matches_tets    both of the above are zero (the surface is the active tetrahedra's boundary)
    volume_div, volume_tets, rel_diff   (1/3) sum A n.c over the faces; sum of active tet volumes; |div/tets - 1|
    n_inward_faces          faces whose right-hand normal points into their owning tetrahedron
    outward                 volume_div > 0, every face owned, and no inward face"""
    faces, tets, points = frame.faces, frame.tets, frame.points
    n_open, ucount = _edge_counts(faces)
    nm_vertices = _nonmanifold_vertices(faces)
    used = np.zeros(len(points), bool)
    used[tets.ravel()] = True
    off = int(np.count_nonzero(~used[faces].all(axis=1))) if len(faces) else 0
    owner, opposite, n_missing, n_interior = surface_owners(points, tets, faces)
    nvec, cen = face_vectors(points, faces)
    has_owner = owner >= 0
    side = np.einsum("ij,ij->i", nvec[has_owner], cen[has_owner] - points[opposite[has_owner]])
    inward = np.zeros(len(faces), bool)
    inward[np.flatnonzero(has_owner)[side < 0.0]] = True
    vol_div = float(math.fsum(np.einsum("ij,ij->i", nvec, cen)) / 6.0)   # (1/3) sum A n.c, nvec = 2 A n
    vol_tets = float(math.fsum(tet_volumes(points, tets)))
    n_without = int(np.count_nonzero(~has_owner))
    return {
        "n_faces": int(len(faces)), "n_tets": int(len(tets)), "n_nodes": int(len(points)),
        "n_nodes_active": int(np.count_nonzero(used)),
        "n_open_directed_edges": n_open,
        "n_boundary_edges": int(np.count_nonzero(ucount == 1)),
        "n_nonmanifold_edges": int(np.count_nonzero(ucount > 2)),
        "n_nonmanifold_vertices": int(len(nm_vertices)),
        "nonmanifold_vertices": nm_vertices.tolist(),
        "faces_on_active_nodes": off == 0, "n_faces_off_active_nodes": off,
        "n_faces_without_owner": n_without, "n_faces_interior": n_interior,
        "n_boundary_faces_missing": n_missing,
        "surface_matches_tets": n_without == 0 and n_missing == 0,
        "volume_div": vol_div, "volume_tets": vol_tets,
        "rel_diff": abs(vol_div / vol_tets - 1.0) if vol_tets > 0 else float("nan"),
        "n_inward_faces": int(np.count_nonzero(inward)),
        "inward_faces": np.flatnonzero(inward).tolist(),
        "outward": bool(vol_div > 0.0 and n_without == 0 and not inward.any()),
    }


def orient_outward(frame) -> tuple[Frame, np.ndarray]:
    """The frame with every surface triangle wound so its right-hand normal points away from its owning active
    tetrahedron (vertices 1 and 2 swapped on the others), and the indices of the faces that were flipped. Face order
    and every patch field are unchanged. ContractError if a face has no active owner (it cannot be oriented)."""
    owner, opposite, _, _ = surface_owners(frame.points, frame.tets, frame.faces)
    if (owner < 0).any():
        raise ContractError("frame {}: {} surface faces are no face of exactly one active tetrahedron; cannot orient "
                            "them".format(frame.k, int(np.count_nonzero(owner < 0))))
    nvec, cen = face_vectors(frame.points, frame.faces)
    flip = np.flatnonzero(np.einsum("ij,ij->i", nvec, cen - frame.points[opposite]) < 0.0)
    if len(flip) == 0:
        return frame, flip
    faces = frame.faces.copy()
    faces[flip] = faces[flip][:, [0, 2, 1]]
    return replace(frame, faces=faces), flip


def mass_checks(frame, rho, rho_liquid, history_row) -> dict:
    """fe_mass = sum phi rho V (active tets) + sum film_thickness rho_l A + sum deep_thickness rho_l A, summed
    exactly, against the history row's mass_kg; film and deep against film_mass_kg and deep_mass_kg where the
    history carries them. `history_row`: column -> value (`history_row(run, row)`). A frame without phi counts
    phi = 1; without film or deep fields those terms are 0 (and reported absent)."""
    area = 0.5 * np.linalg.norm(face_vectors(frame.points, frame.faces)[0], axis=1)
    phi = frame.tet.get("phi")
    vol = tet_volumes(frame.points, frame.tets)
    solid = math.fsum(rho * vol * (phi if phi is not None else 1.0))
    out = {"solid_mass": solid, "volume_tets": math.fsum(vol),
           "volume_phi": math.fsum(vol * (phi if phi is not None else 1.0)),
           "n_partial_tets": int(np.count_nonzero(phi < 1.0)) if phi is not None else 0,
           "phi_min": float(phi.min()) if phi is not None and len(phi) else 1.0}
    total = [solid]
    for key, col in (("film_thickness", "film_mass_kg"), ("deep_thickness", "deep_mass_kg")):
        name = key.split("_")[0]
        if key in frame.patch:
            m = math.fsum(frame.patch[key] * rho_liquid * area)
            total.append(m)
            out[name + "_mass"] = m
            hist = history_row.get(fe_name(col))
            if hist is not None:
                out[name + "_mass_history"] = hist
                out[name + "_abs_diff"] = abs(m - hist)
                out[name + "_rel_diff"] = abs(m / hist - 1.0) if hist != 0 else (0.0 if m == 0 else float("inf"))
        else:
            out[name + "_mass"] = None
    fe_mass = math.fsum(total)
    hist_mass = history_row[fe_name("mass_kg")]
    out.update({"fe_mass": fe_mass, "mass_history": hist_mass, "mass_abs_diff": abs(fe_mass - hist_mass),
                "mass_rel_diff": abs(fe_mass / hist_mass - 1.0)})
    return out


def consistency_checks(frame, table) -> dict:
    """max |f_l(node) - table.liquid_fraction(T(node))| and whether they agree bitwise (table: any object with
    `liquid_fraction(T)`, e.g. `material.MaterialTable`; None skips it); delta_m finite exactly where closure is
    Girin's (0), counted as mismatches; NaN counts of every contract field whose spec forbids NaN."""
    out = {}
    if table is not None and "f_l" in frame.node and "T" in frame.node:
        want = np.asarray(table.liquid_fraction(frame.node["T"]), dtype=np.float64)
        d = np.abs(frame.node["f_l"] - want)
        out["max_abs_dfl"] = float(np.nanmax(d)) if len(d) else 0.0
        out["fl_bitwise"] = bool(np.array_equal(frame.node["f_l"], want))
    else:
        out["max_abs_dfl"], out["fl_bitwise"] = None, None
    if "delta_m" in frame.patch and "closure" in frame.patch:
        girin = frame.patch["closure"] == 0.0
        finite = np.isfinite(frame.patch["delta_m"])
        out["n_delta_m_mismatch"] = int(np.count_nonzero(girin != finite))
        out["n_girin"] = int(np.count_nonzero(girin))
    else:
        out["n_delta_m_mismatch"], out["n_girin"] = None, None
    nan = {}
    for loc in ("node", "tet", "patch"):
        for f in contract.data_fields(loc):
            a = _present(frame, loc).get(f.key)
            if a is not None and not f.nan_allowed:
                c = int(np.count_nonzero(~np.isfinite(a)))
                if c:
                    nan[f.key] = c
    out["nan_forbidden"] = nan
    return out


# ---------------------------------------------------------------------------------------------------------------
# The compact prepared frame
# ---------------------------------------------------------------------------------------------------------------

def compact(frame) -> tuple[Frame, np.ndarray]:
    """The frame on the nodes its active tetrahedra use (ascending finite-element index), renumbered; `node_ids`
    maps compact to finite-element numbering. ContractError if a face uses a node no active tetrahedron uses."""
    node_ids = np.unique(frame.tets.ravel()).astype(np.int64)
    inv = np.full(len(frame.points), -1, np.int64)
    inv[node_ids] = np.arange(len(node_ids))
    faces = inv[frame.faces]
    if (faces < 0).any():
        raise ContractError("frame {}: surface faces use nodes no active tetrahedron uses".format(frame.k))
    return Frame(k=frame.k, time_s=frame.time_s, points=np.ascontiguousarray(frame.points[node_ids]),
                 tets=inv[frame.tets], faces=faces,
                 node={key: np.ascontiguousarray(v[node_ids]) for key, v in frame.node.items()},
                 tet=dict(frame.tet), patch=dict(frame.patch)), node_ids


def derived_patch_arrays(frame, v_hat=V_HAT) -> dict[str, np.ndarray]:
    """Per patch: area (m^2), unit normal of the face's winding (outward on an oriented frame), centroid (m) and theta
    (rad, angle between that normal and v_hat)."""
    v = np.asarray(v_hat, dtype=np.float64)
    v = v / np.linalg.norm(v)
    nvec, centroid = face_vectors(frame.points, frame.faces)
    twice = np.linalg.norm(nvec, axis=1)
    normal = nvec / twice[:, None]
    theta = np.arccos(np.clip(normal @ v, -1.0, 1.0))
    return {"area": 0.5 * twice, "normal": normal, "centroid": centroid, "theta": theta}


def _sym(shape, sizes):
    return tuple(sizes[s] if isinstance(s, str) else s for s in shape)


def prepared_arrays(frame, node_ids, derived) -> dict[str, np.ndarray]:
    """The arrays of `contract.PREPARED_FRAME_ARRAYS` for a compact frame, cast and shape-checked. `derived` holds
    patch arrays by their key without the prefix (area, normal, centroid, theta, and later thickness,
    slurry_depth, liquid_depth). ValueError for an array not in the schema, a wrong shape or a missing required one."""
    sizes = {"n": len(frame.points), "m": len(frame.tets), "f": len(frame.faces)}
    arrays = {"k": frame.k, "time_s": frame.time_s, "points": frame.points, "node_ids": node_ids,
              "tets": frame.tets, "faces": frame.faces}
    arrays.update({"node_" + key: v for key, v in frame.node.items() if key != "node_ids"})
    arrays.update({"tet_" + key: v for key, v in frame.tet.items()})
    arrays.update({"patch_" + key: v for key, v in frame.patch.items()})
    for key, v in derived.items():
        name = "patch_" + key
        if name in arrays and key in frame.patch:
            raise ValueError("derived array {} would replace the frame's patch field".format(key))
        arrays[name] = v
    schema = contract.PREPARED_FRAME_ARRAYS
    unknown = sorted(set(arrays) - set(schema))
    if unknown:
        raise ValueError("arrays not in contract.PREPARED_FRAME_ARRAYS: {}".format(unknown))
    absent = sorted(set(schema) - set(arrays) - contract.PREPARED_FRAME_OPTIONAL)
    if absent:
        raise ValueError("prepared frame lacks required arrays: {}".format(absent))
    out = {}
    for name, value in arrays.items():
        dtype, shape, _ = schema[name]
        a = np.asarray(value)
        if dtype.startswith("i") and a.size and (a.min() < np.iinfo(dtype).min or a.max() > np.iinfo(dtype).max):
            raise ValueError("{} does not fit {}".format(name, dtype))
        a = np.asarray(a, dtype=dtype)               # (ascontiguousarray would turn a scalar into shape (1,))
        if not a.flags.c_contiguous:
            a = a.copy(order="C")
        if a.shape != _sym(shape, sizes):
            raise ValueError("{} has shape {}, expected {}".format(name, a.shape, _sym(shape, sizes)))
        out[name] = a
    return out


def write_npz(path, arrays: dict) -> None:
    """An uncompressed npz without pickles, byte-identical for identical arrays (fixed zip timestamps, sorted names),
    written to a temporary file and renamed into place (an interrupted write leaves no partial file)."""
    tmp = path + ".tmp"
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as zf:
        for name in sorted(arrays):
            info = zipfile.ZipInfo(name + ".npy", date_time=_NPZ_DATE)
            info.compress_type = zipfile.ZIP_STORED
            info.external_attr = 0o644 << 16
            with zf.open(info, "w", force_zip64=True) as fh:
                np.lib.format.write_array(fh, np.asarray(arrays[name]), version=(1, 0), allow_pickle=False)
    os.replace(tmp, path)


def write_prepared_frame(path, frame, node_ids, derived) -> dict[str, np.ndarray]:
    """Write `frame_<k:05d>.npz` (or `path` if it names a file) for a compact, outward frame; returns the arrays.
    Refuses (ContractError) a frame whose surface is not wound outward: the prepared surface is always outward."""
    path = os.fspath(path)
    if os.path.isdir(path):
        path = os.path.join(path, prepared_frame_name(frame.k))
    owner, opposite, _, _ = surface_owners(frame.points, frame.tets, frame.faces)
    nvec, cen = face_vectors(frame.points, frame.faces)
    if (owner < 0).any() or (np.einsum("ij,ij->i", nvec, cen - frame.points[opposite]) < 0.0).any():
        raise ContractError("frame {}: not wound outward or not the active tetrahedra's boundary; "
                            "orient_outward it first".format(frame.k))
    arrays = prepared_arrays(frame, node_ids, derived)
    write_npz(path, arrays)
    return arrays


def prepared_frame_name(k) -> str:
    return "frame_{:05d}.npz".format(int(k))


def load_prepared_arrays(path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as z:
        return {name: z[name] for name in z.files}


def load_prepared_frame(path) -> Frame:
    """The runner's reader (numpy only): a compact Frame whose `node` holds T, f_l and `node_ids` (finite-element
    index per compact node), `tet` phi, and `patch` every `patch_<key>` array under `<key>` -- the frame's patch
    fields and the derived ones (area, normal (f, 3), centroid (f, 3), theta, and the depths when written)."""
    a = load_prepared_arrays(path)
    group = {"node": {}, "tet": {}, "patch": {}}
    for name, v in a.items():
        for prefix in group:
            if name.startswith(prefix + "_") and name != "node_ids":
                group[prefix][name[len(prefix) + 1:]] = v
    group["node"]["node_ids"] = a["node_ids"]
    return Frame(k=int(a["k"]), time_s=float(a["time_s"]), points=a["points"], tets=a["tets"].astype(np.int64),
                 faces=a["faces"].astype(np.int64), node=group["node"], tet=group["tet"], patch=group["patch"])
