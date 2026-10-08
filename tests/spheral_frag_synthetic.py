"""Synthetic finite-element frames and run directories in the frame contract's format (Spheral M1 plan, Task 2).

Small finite-element-like inputs whose answers are known in closed form, so that the core (frame import, thickness,
zones, loads, prepare) can be tested before the real prototype flight exists. Three devices:

  sphere_run      a coarse 100 mm tetrahedral sphere, three frames at 0, 0.5 and 1.0 s (frame 2 with ten nose
                  elements dead and two partly consumed), a history row at every 0.5 s and a run JSON
  dumbbell_frame  two 20 mm spheres joined by a 4 mm neck: thickness 2 r_neck on the neck, 2 R on the spheres
  slab_frame      a 20 x 4 x 12 mm box whose top face carries five columns of known slurry depth (f_l linear in
                  depth, so the f_l = 0.5 crossing is exact under P1 interpolation)

Every name written to a file comes from `spheral_frag.contract.fe_name`; `write_fe_run` writes a run directory exactly
as the finite-element export does (Step 3's `coupled.write_vtk_frame` and `CoupledRun._collection`: pyvista XML, vtu
and vtp on one node array holding every mesh node, the vtu's cells the active tetrahedra only, pvd collections with
the time printed `{:.6g}`; the history CSV at 9 significant digits; the run JSON). The exact answers of each device
are the module constants and the `*_answers` functions below; `tests/fixtures/spheral_frag/README.md` lists them.

Synthetic conventions, stated so that later tasks test against them rather than against the real material:

- The material is `SYNTHETIC_MATERIAL`: solid 2,813 kg/m^3 (the finite element's constant density), liquid
  2,400 kg/m^3, and a liquid fraction **linear in T between 750 and 908 K** (`synthetic_liquid_fraction`). It is not
  AA7075_scheil; nodal `liquid_fraction` equals this law exactly.
- The flow fields follow decision 4 of the M1 plan's review: `p_w`, `tau` and `closure` are written on every frame
  whose step evaluated the surface flow (with or without film). Frame 0 of `sphere_run` is written before any step,
  so nothing was evaluated: `p_w` = `tau` = 0 (the contract's "0 = not evaluated"), `closure` = 1, `delta_m`,
  `kn_local` and `r_droplet` NaN, no film and no release -- while the history row at t = 0 still carries the drag
  the analytic loads would give (a frame without loads, as on a real flight's first frame).
- Evaluated frames carry modified-Newtonian loads in the patch's own inclination theta (outward facet normal against
  v_hat = +x): p_w = p_stag cos^2(theta) for theta < 90 deg, the base pressure BASE_FRACTION p_stag for
  90 <= theta < THETA_EVAL_DEG, and 0 ("not evaluated") beyond, so a load table has bins beyond its last valid one;
  tau = TAU_FRACTION p_stag sin(2 theta) for theta < 90 deg, 0 elsewhere. p_stag is the history's p_w_stag_Pa as
  written to the CSV, so p_w at theta = 0 equals it exactly.
- The surface triangles are wound outward (the contract's "faces: outward"), ordered as `reentry_model.mesh`
  orders boundary faces (ascending sorted node triple).

Prepare side: needs gmsh (the meshes) and pyvista (the writer), so it runs in drama_env, never under Spheral."""
from __future__ import annotations

import json
import math
import os
import shutil
import tempfile

import numpy as np

from spheral_frag import contract
from spheral_frag.contract import Frame, fe_name

G0 = 9.80665                      # m/s^2, reentry_model.constants.G0: load_factor_g = |a_drag| / G0
SIGMA_SB = 5.670374419e-8         # W/(m^2 K^4)
RHO = 2813.0                      # kg/m^3, the finite-element solid (constant; the mesh never expands)
RHO_LIQUID = 2400.0               # kg/m^3, the finite-element liquid
T_SOLIDUS = 750.0                 # K, synthetic linear f_l(T): 0 at and below
T_LIQUIDUS = 908.0                # K, 1 at and above
EMISSIVITY = 0.1
SYNTHETIC_MATERIAL = "synthetic_linear_fl"
CP_MAX = 1.84                     # modified Newtonian stagnation coefficient: p_stag = CP_MAX * q
BASE_FRACTION = 0.01              # lee base pressure on 90 <= theta < THETA_EVAL_DEG, of p_stag
TAU_FRACTION = 0.02               # tau = TAU_FRACTION * p_stag * sin(2 theta) on the windward side
THETA_EVAL_DEG = 150.0            # beyond it the frame carries no loads (p_w = tau = 0)
THETA_GIRIN_DEG = 60.0            # Girin's closure (and a finite delta_m, release) on theta < 60 deg
FILM0 = 5.0e-5                    # m, film_thickness = FILM0 cos(theta) on the windward side (sphere, dumbbell off)
DELTA_M0 = 2.5e-4                 # m, delta_m on the Girin patches
RELEASE0 = 1.0e-3                 # kg/m^2 per macro step, release_rate = RELEASE0 cos(theta) on the Girin patches
R_DROPLET = 2.0e-5                # m, where something is released
KN_LOCAL = 1.0e-3
WE_S0 = 50.0
Q_CONV0 = 2.0e6                   # W/m^2, q_conv = Q_CONV0 (0.1 + 0.9 max(cos theta, 0))
V_HAT = np.array(contract.V_HAT)
MACRO_STEP_S = 0.5
SEED = 12345

# sphere_run
SPHERE_R = 0.05
SPHERE_H = 8.0e-3                 # m, surface element size
SPHERE_H_CORE = 20.0e-3           # m: 753 nodes, 2,497 tetrahedra, 1,194 patches (gmsh 4, HXT, band = 0)
SPHERE_TIMES = (0.0, 0.5, 1.0)
SPHERE_T_WALL = (300.0, 900.0, 960.0)   # K at r = R; T = 300 + (T_wall - 300) (r / R)^2 in every frame
SPHERE_N_DEAD = 10                # frame 2: the owners of the ten patches nearest the nose (first ten distinct)
SPHERE_PHI = (0.5, 0.25)          # frame 2: phi of the owners of the two patches nearest the equator
SPHERE_RUN_NAME = "synthetic_sphere_d100.00mm_seed12345"

# dumbbell_frame
DUMBBELL = dict(R=0.02, r_neck=4.0e-3, L_neck=0.02, h=1.5e-3, h_far=5.0e-3)
DUMBBELL_RUN_NAME = "synthetic_dumbbell_R20.0mm_neck4.0mm"
DUMBBELL_T = 600.0

# slab_frame
SLAB = dict(L=0.02, W=0.004, H=0.012, h=1.0e-3)
SLAB_COLUMN_CELLS = 4                                       # cells per column along x
SLAB_SLURRY_DEPTH = (0.15e-3, 1.0e-3, 2.5e-3, 4.0e-3, 1.0e-3)   # m, depth of f_l = 0.5 per column
SLAB_DELTA_M = (3.0e-4, 3.0e-4, 3.0e-4, 3.0e-4, np.nan)         # m; column 4: delta_m undefined (closure Couette)
SLAB_FILM = 1.0e-4                # m, film_thickness on every top patch
SLAB_DEEP = 2.0e-4                # m (mass per area), deep_thickness on every top patch
SLAB_GRADIENT = 24.0e-3           # m: f_l = 0.5 + (D_c - d) / SLAB_GRADIENT at depth d, strictly inside (0, 1)
SLAB_RUN_NAME = "synthetic_slab_L20.0mm_H12.0mm"


# ------------------------------------------------------------------------------------------------------------ helpers
def r9(x):
    """The value as the history CSV holds it (9 significant digits)."""
    return float("{:.9g}".format(x))


def synthetic_liquid_fraction(T):
    return np.clip((np.asarray(T, dtype=float) - T_SOLIDUS) / (T_LIQUIDUS - T_SOLIDUS), 0.0, 1.0)


def synthetic_temperature(f_l):
    """Inverse of the linear law on 0 < f_l < 1."""
    return T_SOLIDUS + np.asarray(f_l, dtype=float) * (T_LIQUIDUS - T_SOLIDUS)


def tet_volumes(points, tets):
    x = points[tets]
    J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
    return np.abs(np.linalg.det(J)) / 6.0


def outward_surface(points, tets):
    """Boundary triangles of `tets` (faces used by exactly one), wound so the right-hand normal points away from the
    owning tetrahedron's opposite vertex; in `reentry_model.mesh.boundary_faces`' order."""
    from reentry_model import mesh
    faces, opposite = mesh.boundary_faces(np.asarray(tets, dtype=np.int64))
    a, b, c = (points[faces[:, i]] for i in range(3))
    n = np.cross(b - a, c - a)
    inward = np.einsum("ij,ij->i", n, (a + b + c) / 3.0 - points[opposite]) < 0.0
    faces = faces.copy()
    faces[inward] = faces[inward][:, [0, 2, 1]]
    return faces.astype(np.int64), _owners(tets, faces)


def _owners(tets, faces):
    """Index of the tetrahedron owning each boundary face."""
    from reentry_model.mesh import _FACE_OF_VERTEX
    all_faces = np.sort(np.asarray(tets)[:, _FACE_OF_VERTEX].reshape(-1, 3), axis=1)
    lookup = {tuple(f): i // 4 for i, f in enumerate(all_faces)}
    return np.array([lookup[tuple(f)] for f in np.sort(faces, axis=1)], dtype=np.int64)


def patch_geometry(points, faces, v_hat=V_HAT):
    """(area, outward unit normal, centroid, theta [rad]) per patch."""
    a, b, c = (points[faces[:, i]] for i in range(3))
    n = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(n, axis=1)
    normal = n / (2.0 * area[:, None])
    theta = np.arccos(np.clip(normal @ np.asarray(v_hat, dtype=float), -1.0, 1.0))
    return area, normal, (a + b + c) / 3.0, theta


def newtonian_loads(theta, p_stag):
    """(p_w, tau) of an evaluated frame: modified Newtonian windward, base pressure to THETA_EVAL_DEG, 0 beyond."""
    deg = np.degrees(theta)
    windward = deg < 90.0
    p = np.where(windward, p_stag * np.cos(theta) ** 2, np.where(deg < THETA_EVAL_DEG, BASE_FRACTION * p_stag, 0.0))
    tau = np.where(windward, TAU_FRACTION * p_stag * np.sin(2.0 * theta), 0.0)
    return p, tau


def drag(theta, area, p_w, tau):
    """D = sum A (p cos theta + tau sin theta), N along -v_hat: the frame's own loads summed over its patches."""
    return float(np.sum(area * (p_w * np.cos(theta) + tau * np.sin(theta))))


def patch_fields(theta, T_patch, p_stag, evaluated=True, melt=True):
    """Every patch field of the contract for an analytic frame (see the module docstring)."""
    n = len(theta)
    cos = np.cos(theta)
    deg = np.degrees(theta)
    girin = (deg < THETA_GIRIN_DEG) & evaluated & melt
    f = {
        "q_conv": Q_CONV0 * (0.1 + 0.9 * np.maximum(cos, 0.0)),
        "q_rad": EMISSIVITY * SIGMA_SB * T_patch ** 4,
        "T_patch": T_patch.copy(),
        "film_thickness": FILM0 * np.maximum(cos, 0.0) if (evaluated and melt) else np.zeros(n),
        "film_T": T_patch.copy(),
        "we_s": np.where(girin, WE_S0 * cos, 0.0),
        "closure": np.where(girin, 0.0, 1.0),
        "kn_local": np.full(n, KN_LOCAL) if evaluated else np.full(n, np.nan),
        "r_droplet": np.where(girin, R_DROPLET, np.nan),
        "release_rate": np.where(girin, RELEASE0 * cos, 0.0),
        "delta_m": np.where(girin, DELTA_M0, np.nan),
        "deep_thickness": np.zeros(n),
    }
    if evaluated:
        f["p_w"], f["tau"] = newtonian_loads(theta, p_stag)
    else:
        f["p_w"], f["tau"] = np.zeros(n), np.zeros(n)
    return f


def make_frame(k, t, points, tets, T, phi, patch_overrides=None, p_stag=0.0, evaluated=True, melt=True):
    """A contract.Frame on (all) `points` and the active `tets`, its surface and analytic patch fields."""
    points = np.ascontiguousarray(points, dtype=np.float64)
    tets = np.ascontiguousarray(tets, dtype=np.int64)
    T = np.ascontiguousarray(T, dtype=np.float64)
    faces, _ = outward_surface(points, tets)
    area, normal, centroid, theta = patch_geometry(points, faces)
    T_patch = T[faces].mean(axis=1)
    patch = patch_fields(theta, T_patch, p_stag, evaluated, melt)
    patch.update(patch_overrides or {})
    patch = {key: np.ascontiguousarray(v, dtype=np.float64) for key, v in patch.items()}
    return Frame(k=k, time_s=float(t), points=points, tets=tets, faces=faces,
                 node={"T": T, "f_l": synthetic_liquid_fraction(T)},
                 tet={"phi": np.ascontiguousarray(phi, dtype=np.float64)}, patch=patch)


def frame_geometry(frame):
    return patch_geometry(frame.points, frame.faces)


def frame_drag(frame):
    area, _, _, theta = frame_geometry(frame)
    return drag(theta, area, frame.patch["p_w"], frame.patch["tau"])


def frame_masses(frame):
    """(element mass sum phi rho V, film mass, deep mass) in kg."""
    area = frame_geometry(frame)[0]
    solid = math.fsum(frame.tet["phi"] * RHO * tet_volumes(frame.points, frame.tets))
    film = math.fsum(frame.patch["film_thickness"] * RHO_LIQUID * area)
    deep = math.fsum(frame.patch["deep_thickness"] * RHO_LIQUID * area)
    return solid, film, deep


def frame_mass(frame):
    return math.fsum(frame_masses(frame))


def sprayed_this_step(frame):
    """sum release_rate * A: the step's sprayed mass (release_rate is kg/m^2 per macro step)."""
    return math.fsum(frame.patch["release_rate"] * frame_geometry(frame)[0])


# ------------------------------------------------------------------------------------------------------------ history
HISTORY_MAIN = ("time_s", "altitude_km", "velocity_kms", "mass_kg", "flight_path_deg", "lat_deg", "lon_deg",
                "density_kgm3", "dynamic_pressure_Pa", "load_factor_g")
HISTORY_EXTRA = ("surface_T_max_K", "n_active_elements")     # not in the contract: readers must ignore them
HISTORY_MELT = ("p_w_stag_Pa", "film_mass_kg", "deep_mass_kg", "sprayed_mass_kg")


def flight_state(t):
    """A smooth synthetic flight state at time t, each value as the CSV holds it (9 digits)."""
    velocity = r9(7.5 - 0.002 * t)                 # km/s
    density = r9(3.0e-5 * math.exp(t / 5.0))       # kg/m^3
    q = r9(0.5 * density * (velocity * 1e3) ** 2)
    return {"time_s": r9(t), "altitude_km": r9(77.5 - 0.06 * t), "velocity_kms": velocity,
            "flight_path_deg": r9(-0.959331 - 0.01 * t), "lat_deg": r9(29.546 + 0.004 * t),
            "lon_deg": r9(-82.134 + 0.003 * t), "density_kgm3": density, "dynamic_pressure_Pa": q,
            "p_w_stag_Pa": r9(CP_MAX * q)}


def history_rows(frames, states, drags, sprayed_increments):
    """One history row per frame (frames_every = 1, one frame per macro step). `drags` is the drag the history
    records at the row (N); mass and the melt accounts are the frame's own."""
    rows, sprayed = [], 0.0
    for frame, state, D, ds in zip(frames, states, drags, sprayed_increments):
        solid, film, deep = frame_masses(frame)
        mass = r9(math.fsum([solid, film, deep]))
        sprayed += ds
        row = dict(state)
        row.update({"mass_kg": mass, "load_factor_g": r9(D / (mass * G0)),
                    "surface_T_max_K": r9(float(frame.patch["T_patch"].max())),
                    "n_active_elements": float(len(frame.tets)),
                    "film_mass_kg": r9(film), "deep_mass_kg": r9(deep), "sprayed_mass_kg": r9(sprayed)})
        rows.append(row)
    return {col: np.array([r[col] for r in rows]) for col in HISTORY_MAIN + HISTORY_EXTRA + HISTORY_MELT}


# ------------------------------------------------------------------------------------------------------------ run json
def _set_path(doc, dotted, value):
    keys = dotted.split(".")
    for key in keys[:-1]:
        doc = doc.setdefault(key, {})
    doc[keys[-1]] = value


def run_json(run_name, diameter_mm, mass_kg, initial_fe_mass_kg, n_frames, synthetic):
    """The run JSON: the contract's run items at their dotted paths, plus `results` and a `synthetic` block that
    records the device's parameters (provenance; not a contract item)."""
    doc = {"schema_version": 1, "run_name": None, "inputs": {}, "settings": {"thermal": "fem", "heating": "physics", "melt": "on"},
           "results": {"initial_mass_kg": initial_fe_mass_kg, "n_frames": n_frames}}
    values = {"run_name": run_name, "diameter_mm": diameter_mm, "initial_mass_kg": mass_kg, "atmosphere": "us76",
              "wind": "none", "macro_step_s": MACRO_STEP_S, "frames_every": 1, "seed": SEED,
              "material": SYNTHETIC_MATERIAL, "v_hat": [float(v) for v in V_HAT]}
    for f in contract.fields("run"):
        _set_path(doc, f.name, values[f.key])
    doc["synthetic"] = synthetic
    return doc


# ------------------------------------------------------------------------------------------------------------ meshes
_MESHES = {}


def _cached(key, build):
    if key not in _MESHES:
        _MESHES[key] = build()
    points, tets = _MESHES[key]
    return points.copy(), tets.copy()


def sphere_mesh(R=SPHERE_R, h=SPHERE_H, h_core=SPHERE_H_CORE):
    """(points, tets) of `reentry_model.mesh.sphere_mesh` with band = 0 (as conftest's coarse sphere), meshed in a
    temporary directory so that nothing is cached under reentry_model_output/."""
    def build():
        from reentry_model import mesh
        d = tempfile.mkdtemp(prefix="spheral_frag_synthetic_")
        try:
            m = mesh.sphere_mesh(R, h, h_core, d, band=0.0)
        finally:
            shutil.rmtree(d, ignore_errors=True)
        return np.asarray(m.points, dtype=np.float64), np.asarray(m.tets, dtype=np.int64)
    return _cached(("sphere", R, h, h_core), build)


def dumbbell_mesh(R=0.02, r_neck=4e-3, L_neck=0.02, h=1.5e-3, h_far=5e-3):
    """Two spheres of radius R centred at x = -+(L_neck/2 + R), fused (gmsh OCC) with a cylinder of radius r_neck
    from centre to centre; size h within 3 mm of the neck region, h_far elsewhere (Delaunay, single thread)."""
    def build():
        import gmsh
        a = 0.5 * L_neck + R
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.option.setNumber("General.NumThreads", 1)
            gmsh.model.add("dumbbell")
            occ = gmsh.model.occ
            s1, s2 = occ.addSphere(-a, 0.0, 0.0, R), occ.addSphere(a, 0.0, 0.0, R)
            cyl = occ.addCylinder(-a, 0.0, 0.0, 2.0 * a, 0.0, 0.0, r_neck)
            occ.fuse([(3, s1)], [(3, s2), (3, cyl)])
            occ.synchronize()
            box = gmsh.model.mesh.field.add("Box")
            reach = 3.0e-3
            for name, value in (("VIn", h), ("VOut", h_far), ("XMin", -0.5 * L_neck - reach),
                                ("XMax", 0.5 * L_neck + reach), ("YMin", -r_neck - reach), ("YMax", r_neck + reach),
                                ("ZMin", -r_neck - reach), ("ZMax", r_neck + reach), ("Thickness", 4.0e-3)):
                gmsh.model.mesh.field.setNumber(box, name, value)
            gmsh.model.mesh.field.setAsBackgroundMesh(box)
            for option in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
                gmsh.option.setNumber(option, 0)
            gmsh.option.setNumber("Mesh.Algorithm3D", 1)
            gmsh.model.mesh.generate(3)
            tags, coords, _ = gmsh.model.mesh.getNodes()
            types, _, conn = gmsh.model.mesh.getElements(3)
            tet_conn = np.asarray(conn[list(types).index(4)], dtype=np.int64).reshape(-1, 4)
        finally:
            gmsh.finalize()
        tags = np.asarray(tags, dtype=np.int64)
        coords = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
        used = np.unique(tet_conn)                          # drop nodes no tetrahedron uses (curve/point nodes kept)
        index = np.full(tags.max() + 1, -1, dtype=np.int64)
        index[tags] = np.arange(len(tags))
        keep = index[used]
        renumber = np.full(tags.max() + 1, -1, dtype=np.int64)
        renumber[used] = np.arange(len(used))
        return coords[keep], renumber[tet_conn]
    return _cached(("dumbbell", R, r_neck, L_neck, h, h_far), build)


def slab_mesh(L=0.02, W=0.004, H=0.012, h=1e-3):
    from reentry_model import mesh
    m = mesh.box_mesh(L, W, H, h)
    return np.asarray(m.points, dtype=np.float64), np.asarray(m.tets, dtype=np.int64)


# ------------------------------------------------------------------------------------------------------------ sphere
def sphere_temperature(points, T_wall, R=SPHERE_R, T_core=300.0):
    r2 = np.einsum("ij,ij->i", points, points) / R ** 2
    return T_core + (T_wall - T_core) * r2


def sphere_frame(R=SPHERE_R, h=SPHERE_H, k=0, t=0.0, T=None, fields=None, h_core=SPHERE_H_CORE, p_stag=None,
                 evaluated=True, melt=True):
    """One frame of the intact coarse sphere. `T`: None (the frame-1 radial law, T_wall = 900 K), a scalar, an array
    over the nodes or a callable of the points. `fields`: contract key -> array overriding any node, tet or patch
    field (f_l, phi, p_w, ...). `p_stag` defaults to the flight state's at t."""
    points, tets = sphere_mesh(R, h, h_core)
    if T is None:
        T = sphere_temperature(points, SPHERE_T_WALL[1], R)
    elif callable(T):
        T = T(points)
    T = np.broadcast_to(np.asarray(T, dtype=float), (len(points),)).copy()
    if p_stag is None:
        p_stag = flight_state(t)["p_w_stag_Pa"]
    frame = make_frame(k, t, points, tets, T, np.ones(len(tets)), p_stag=p_stag, evaluated=evaluated, melt=melt)
    for key, value in (fields or {}).items():
        where = {"node": frame.node, "tet": frame.tet, "patch": frame.patch}[contract.field_spec(key).location]
        where[key] = np.ascontiguousarray(value, dtype=np.float64)
    return frame


def sphere_run_parts(R=SPHERE_R, h=SPHERE_H, h_core=SPHERE_H_CORE):
    """(frames, history, run JSON, info) of the three-frame sphere run; `info` holds the dead and partly consumed
    elements (indices into the full tetrahedron list) and the drags."""
    points, tets = sphere_mesh(R, h, h_core)
    faces0, owners0 = outward_surface(points, tets)
    theta0 = patch_geometry(points, faces0)[3]
    dead = []
    for i in np.argsort(theta0, kind="stable"):            # the nose: owners of the patches nearest theta = 0
        if owners0[i] not in dead:
            dead.append(int(owners0[i]))
        if len(dead) == SPHERE_N_DEAD:
            break
    active = np.ones(len(tets), dtype=bool)
    active[dead] = False
    tets2 = tets[active]
    faces2, owners2 = outward_surface(points, tets2)       # owners2 index tets2
    theta2 = patch_geometry(points, faces2)[3]
    partial = []
    for i in np.argsort(np.abs(theta2 - 0.5 * np.pi), kind="stable"):   # near the equator, far from the dead nose
        if owners2[i] not in partial:
            partial.append(int(owners2[i]))
        if len(partial) == len(SPHERE_PHI):
            break
    phi2 = np.ones(len(tets2))
    phi2[partial] = SPHERE_PHI

    states = [flight_state(t) for t in SPHERE_TIMES]
    frames, drags = [], []
    for k, (t, T_wall, state) in enumerate(zip(SPHERE_TIMES, SPHERE_T_WALL, states)):
        T = sphere_temperature(points, T_wall, R)
        tk, phik = (tets2, phi2) if k == 2 else (tets, np.ones(len(tets)))
        frame = make_frame(k, t, points, tk, T, phik, p_stag=state["p_w_stag_Pa"], evaluated=k > 0)
        frames.append(frame)
        if k == 0:     # the history's drag at t = 0 is the one the analytic loads would give: the frame has none
            ev = make_frame(k, t, points, tk, T, phik, p_stag=state["p_w_stag_Pa"], evaluated=True)
            drags.append(frame_drag(ev))
        else:
            drags.append(frame_drag(frame))
    sprayed = [sprayed_this_step(f) for f in frames]
    history = history_rows(frames, states, drags, sprayed)
    mass0 = frame_mass(frames[0])
    synthetic = {"device": "sphere_run", "R_m": R, "h_surface_m": h, "h_core_m": h_core, "rho_kgm3": RHO,
                 "rho_liquid_kgm3": RHO_LIQUID, "fl_law": "linear in T", "T_solidus_K": T_SOLIDUS,
                 "T_liquidus_K": T_LIQUIDUS, "T_wall_K": list(SPHERE_T_WALL),
                 "dead_elements_frame2": dead, "partial_elements_frame2": partial, "phi_frame2": list(SPHERE_PHI)}
    doc = run_json(SPHERE_RUN_NAME, 2e3 * R, RHO * 4.0 / 3.0 * math.pi * R ** 3, mass0, len(frames), synthetic)
    info = {"dead": dead, "partial": partial, "tets_full": tets, "drags": drags, "sprayed": sprayed}
    return frames, history, doc, info


def sphere_run(directory, **kw):
    frames, history, doc, _ = sphere_run_parts(**kw)
    return write_fe_run(directory, frames, history, doc)


def sphere_answers(R=SPHERE_R):
    """Closed forms on the sphere (continuum; the faceted fixture differs by its faceting error)."""
    p_stag = flight_state(SPHERE_TIMES[1])["p_w_stag_Pa"]
    return {
        "thickness_m": 2.0 * R,
        "windward_pressure_drag_over_p_stag": 0.5 * math.pi * R ** 2,           # integral of p_stag cos^2 cos dA
        "drag_over_p_stag": math.pi * R ** 2 * (0.5 + TAU_FRACTION - 0.75 * BASE_FRACTION),
        "p_stag_frame1_Pa": p_stag,
        # T = 300 + (T_wall - 300)(r/R)^2 with the linear f_l: f_l = 0.5 at 829 K, f_l = 1 at 908 K
        "slurry_depth_m": {k: R * (1.0 - math.sqrt((0.5 * (T_SOLIDUS + T_LIQUIDUS) - 300.0) / (Tw - 300.0)))
                           for k, Tw in enumerate(SPHERE_T_WALL) if Tw > 0.5 * (T_SOLIDUS + T_LIQUIDUS)},
        "liquid_depth_m": {k: R * (1.0 - math.sqrt((T_LIQUIDUS - 300.0) / (Tw - 300.0)))
                           for k, Tw in enumerate(SPHERE_T_WALL) if Tw > T_LIQUIDUS},
    }


# ------------------------------------------------------------------------------------------------------------ dumbbell
def dumbbell_volume(R=0.02, r_neck=4e-3, L_neck=0.02):
    """Exact volume of the fused solid: two spheres, the neck between them and, at each end, the part of the
    cylinder (centre to centre) that lies outside its sphere."""
    end = math.pi * r_neck ** 2 * R - 2.0 * math.pi / 3.0 * (R ** 3 - (R ** 2 - r_neck ** 2) ** 1.5)
    return 2.0 * 4.0 / 3.0 * math.pi * R ** 3 + math.pi * r_neck ** 2 * L_neck + 2.0 * end


def dumbbell_parts(R=0.02, r_neck=4e-3, L_neck=0.02, h=1.5e-3, h_far=5e-3):
    points, tets = dumbbell_mesh(R, r_neck, L_neck, h, h_far)
    state = flight_state(0.0)
    frame = make_frame(0, 0.0, points, tets, np.full(len(points), DUMBBELL_T), np.ones(len(tets)),
                       p_stag=state["p_w_stag_Pa"], evaluated=True, melt=False)
    history = history_rows([frame], [state], [frame_drag(frame)], [0.0])
    synthetic = {"device": "dumbbell_frame", "R_m": R, "r_neck_m": r_neck, "L_neck_m": L_neck, "h_neck_m": h,
                 "h_far_m": h_far, "sphere_centres_x_m": [-(0.5 * L_neck + R), 0.5 * L_neck + R],
                 "rho_kgm3": RHO, "T_K": DUMBBELL_T}
    doc = run_json(DUMBBELL_RUN_NAME, 2e3 * R, RHO * dumbbell_volume(R, r_neck, L_neck), frame_mass(frame), 1,
                   synthetic)
    return [frame], history, doc


def dumbbell_frame(R=0.02, r_neck=4e-3, L_neck=0.02, h=1.5e-3, h_far=5e-3):
    return dumbbell_parts(R, r_neck, L_neck, h, h_far)[0][0]


def dumbbell_neck_patches(frame, R=0.02, r_neck=4e-3, L_neck=0.02, tol=1e-9):
    """Lateral patches of the neck: every vertex on the cylinder r = r_neck (within `tol`). The free lateral surface
    runs between the junction circles at |x| = 0.5 L_neck + R - sqrt(R^2 - r_neck^2) (10.404 mm), so it is 20.808 mm
    long at r = r_neck, not L_neck (Task 5's measurement: 567 patches, all thinner than 8.8 mm)."""
    radial = np.linalg.norm(frame.points[frame.faces][:, :, 1:], axis=2)       # (f, 3) vertex distance from the axis
    return (np.abs(radial - r_neck) < tol).all(axis=1)


def dumbbell_sphere_patches(frame, R=0.02, r_neck=4e-3, L_neck=0.02, margin_deg=10.0):
    """Sphere patches whose inward ray through the sphere's centre exits on the same sphere (thickness 2R): the
    antipode lies more than `margin_deg` outside the cap the neck occupies (half-angle asin(r_neck/R)); patches of
    the neck's lateral surface (whose centroids can lie within 5 % of R of the sphere) are excluded."""
    a = 0.5 * L_neck + R
    c = frame.points[frame.faces].mean(axis=1)
    side = np.sign(c[:, 0])
    u = (c - np.column_stack([side * a, np.zeros(len(c)), np.zeros(len(c))])) / R
    on_sphere = np.abs(np.linalg.norm(u, axis=1) - 1.0) < 0.05
    cap = math.cos(math.asin(r_neck / R) + math.radians(margin_deg))
    toward_neck_of_antipode = -u[:, 0] * (-side)       # cos(angle between the antipode and the neck direction)
    return (on_sphere & (toward_neck_of_antipode < cap) & (np.abs(c[:, 0]) > 0.5 * L_neck)
            & ~dumbbell_neck_patches(frame, R, r_neck, L_neck))


# ------------------------------------------------------------------------------------------------------------ slab
def slab_columns(points, L=0.02, h=1e-3):
    """Column of each node plane along x: plane i (x = i h) belongs to column min(i // 4, 4)."""
    i = np.rint(points[:, 0] / h).astype(np.int64)
    return np.minimum(i // SLAB_COLUMN_CELLS, len(SLAB_SLURRY_DEPTH) - 1)


def slab_parts(L=0.02, W=0.004, H=0.012, h=1e-3):
    points, tets = slab_mesh(L, W, H, h)
    col = slab_columns(points, L, h)
    depth = H - points[:, 2]
    D = np.asarray(SLAB_SLURRY_DEPTH)[col]
    f_l = 0.5 + (D - depth) / SLAB_GRADIENT
    T = synthetic_temperature(f_l)
    state = flight_state(0.0)
    faces, _ = outward_surface(points, tets)
    area, normal, centroid, theta = patch_geometry(points, faces)
    top = normal[:, 2] > 1.0 - 1e-12
    cube = np.floor(centroid[:, 0] / h).astype(np.int64)
    patch_col = np.minimum(cube // SLAB_COLUMN_CELLS, len(SLAB_SLURRY_DEPTH) - 1)
    delta_m = np.where(top, np.asarray(SLAB_DELTA_M)[patch_col], np.nan)
    overrides = {"film_thickness": np.where(top, SLAB_FILM, 0.0), "deep_thickness": np.where(top, SLAB_DEEP, 0.0),
                 "delta_m": delta_m, "closure": np.where(np.isfinite(delta_m), 0.0, 1.0),
                 "release_rate": np.zeros(len(faces)), "r_droplet": np.full(len(faces), np.nan),
                 "we_s": np.zeros(len(faces))}
    frame = make_frame(0, 0.0, points, tets, T, np.ones(len(tets)), overrides, p_stag=state["p_w_stag_Pa"],
                       evaluated=True)
    history = history_rows([frame], [state], [frame_drag(frame)], [0.0])
    synthetic = {"device": "slab_frame", "L_m": L, "W_m": W, "H_m": H, "h_m": h, "exposed": "top face z = H",
                 "column_cells": SLAB_COLUMN_CELLS, "slurry_depth_m": list(SLAB_SLURRY_DEPTH),
                 "delta_m_m": [None if math.isnan(v) else v for v in SLAB_DELTA_M], "film_m": SLAB_FILM,
                 "deep_m": SLAB_DEEP, "fl_gradient_m": SLAB_GRADIENT, "T_solidus_K": T_SOLIDUS,
                 "T_liquidus_K": T_LIQUIDUS}
    doc = run_json(SLAB_RUN_NAME, 1e3 * L, RHO * L * W * H, frame_mass(frame), 1, synthetic)
    return [frame], history, doc


def slab_frame(L=0.02, H=0.012, h=1e-3, W=0.004):
    return slab_parts(L, W, H, h)[0][0]


def slab_exact_patches(frame, h=1e-3):
    """(top patch mask, column per patch, exact mask): a top patch is exact when both node planes of its cube lie in
    one column, so the vertical line below its centroid sees one column's linear f_l only."""
    area, normal, centroid, theta = patch_geometry(frame.points, frame.faces)
    top = normal[:, 2] > 1.0 - 1e-12
    cube = np.floor(centroid[:, 0] / h).astype(np.int64)
    n_col = len(SLAB_SLURRY_DEPTH)
    col = np.minimum(cube // SLAB_COLUMN_CELLS, n_col - 1)
    exact = top & (np.minimum(cube // SLAB_COLUMN_CELLS, n_col - 1) == np.minimum((cube + 1) // SLAB_COLUMN_CELLS,
                                                                                    n_col - 1))
    return top, np.where(top, col, -1), exact


def slab_answers():
    """Per column: slurry depth, liquid depth (0: the surface f_l stays below 1), film, delta_m, deep, and the layer
    film + slurry (+ deep with include_deep)."""
    out = []
    for c, (D, dm) in enumerate(zip(SLAB_SLURRY_DEPTH, SLAB_DELTA_M)):
        out.append({"column": c, "slurry_depth_m": D, "liquid_depth_m": 0.0, "film_m": SLAB_FILM, "delta_m_m": dm,
                    "deep_m": SLAB_DEEP, "surface_f_l": 0.5 + D / SLAB_GRADIENT, "layer_m": SLAB_FILM + D,
                    "layer_with_deep_m": SLAB_FILM + D + SLAB_DEEP})
    return out


# ------------------------------------------------------------------------------------------------------------ writer
PVD_TEMPLATE = ('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
                '<Collection>\n{}</Collection>\n</VTKFile>\n')     # reentry_model.coupled.PVD_TEMPLATE


def write_history_csv(path, history):
    import csv
    names = list(history)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(names)
        for i in range(len(history[names[0]])):
            w.writerow(["{:.9g}".format(history[k][i]) for k in names])


def write_fe_run(directory, frames, history, run_doc):
    """Write `<directory>/<run>.csv`, `<run>.json` and `<run>/vtk/{field,surface}_<k>.{vtu,vtp}` with the two pvd
    collections, as the finite-element export does. Returns the run directory `<directory>/<run>`."""
    import pyvista as pv
    name = run_doc
    for part in fe_name("run_name").split("."):
        name = name[part]
    run_dir = os.path.join(directory, name)
    vtk_dir = os.path.join(run_dir, "vtk")
    os.makedirs(vtk_dir, exist_ok=True)
    for frame in frames:
        tets = frame.tets
        cells = np.hstack([np.full((len(tets), 1), 4), tets]).ravel()
        grid = pv.UnstructuredGrid(cells, np.full(len(tets), pv.CellType.TETRA), frame.points)
        for f in contract.data_fields("node"):
            if f.key in frame.node:
                grid.point_data[f.name] = frame.node[f.key]
        for f in contract.data_fields("tet"):
            if f.key in frame.tet:
                grid.cell_data[f.name] = frame.tet[f.key]
        grid.save(os.path.join(vtk_dir, "field_{}.vtu".format(frame.k)))
        poly = pv.PolyData(frame.points, np.hstack([np.full((len(frame.faces), 1), 3), frame.faces]).ravel())
        poly.point_data[fe_name("T")] = frame.node["T"]
        for f in contract.data_fields("patch"):
            if f.key in frame.patch:
                poly.cell_data[f.name] = frame.patch[f.key]
        poly.save(os.path.join(vtk_dir, "surface_{}.vtp".format(frame.k)))
    for kind, ext in (("field", "vtu"), ("surface", "vtp")):
        entries = "".join('<DataSet timestep="{:.6g}" file="{}_{}.{}"/>\n'.format(f.time_s, kind, f.k, ext)
                          for f in frames)
        with open(os.path.join(vtk_dir, kind + ".pvd"), "w") as fh:
            fh.write(PVD_TEMPLATE.format(entries))
    write_history_csv(os.path.join(directory, name + ".csv"), history)
    with open(os.path.join(directory, name + ".json"), "w") as fh:
        json.dump(run_doc, fh, indent=2, default=str)
        fh.write("\n")
    return run_dir


def read_frame_pyvista(run_dir, k):
    """A Frame read back with pyvista directly through the contract's names (the tests' independent reader)."""
    import pyvista as pv
    grid = pv.read(os.path.join(run_dir, "vtk", "field_{}.vtu".format(k)))
    poly = pv.read(os.path.join(run_dir, "vtk", "surface_{}.vtp".format(k)))
    tets = np.asarray(grid.cells_dict[pv.CellType.TETRA], dtype=np.int64)
    faces = np.asarray(poly.faces, dtype=np.int64).reshape(-1, 4)[:, 1:]
    node = {f.key: np.asarray(grid.point_data[f.name]) for f in contract.data_fields("node")
            if f.name in grid.point_data}
    tet = {f.key: np.asarray(grid.cell_data[f.name]) for f in contract.data_fields("tet") if f.name in grid.cell_data}
    patch = {f.key: np.asarray(poly.cell_data[f.name]) for f in contract.data_fields("patch")
             if f.name in poly.cell_data}
    return Frame(k=k, time_s=float("nan"), points=np.asarray(grid.points), tets=tets, faces=faces, node=node, tet=tet,
                 patch=patch), np.asarray(poly.points)
