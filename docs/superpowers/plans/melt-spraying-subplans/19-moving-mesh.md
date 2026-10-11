# Sub-plan 19 — the moving 3D mesh (stage 1 of the layered-skin design) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first (required reading for every Step 3 sub-plan), then the spec.

**Goal:** Make the 3D model's mesh move with its surface — corners moved by a prescribed displacement, a deforming zone
beneath them, the enthalpy remapped with exact books, and the mesh rebuilt inside the same surface before any element
degrades — on both thermal backends, with no skins yet.

**Architecture:** A new module `reentry_model/movemesh.py` holds the geometry helpers, the deforming zone (`ZoneMap`),
and `MovingMesh`, which owns a `VolumeMesh`, its thermal solver and its material: `move(u_corners)` moves the corners and
the zone, carries the nodal enthalpy to the moved nodes by its gradient and closes the books with one uniform offset,
reporting the volume and enthalpy each corner swept out of the 3D model; `rebuild()` regenerates the tetrahedra with
gmsh inside the unchanged surface triangulation and carries the enthalpy over; `move_with_policy` decides when.
Each thermal backend gains `update_geometry(elements=None)`. Nothing in the model's run path uses any of this until
stage 2 (sub-plan 20): stage 1 is verified on its own, with prescribed motion.

**Tech Stack:** Python 3.12 (`drama_env`), NumPy, SciPy, gmsh 4.15.2, scikit-fem backend, FEniCSx 0.11 backend
(`fenicsx_env`), pytest.

**Spec:** `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` (committed cc15b3f), §5.1 (the corners'
tangential smoothing), §5.3 (the moving 3D mesh) and §11 (level 1's moving-mesh checks, build stage 1). Sub-plans 20
(skins, heat only), 21 (flow) and 22 (spraying and the flight) build on this one and are written in full after it is
built and measured, as the spec gates each stage on the one before.

**How this plan's code was checked (2026-10-11).** Every code and diff block below was applied to a copy of
`prototype/work-2026-10-08-skin/code/` at its head 49d88e0 and tested there: the new module's 26 tests, the three new
FEniCSx tests, and the whole unit tier, which goes from 303 to 329 passed with 1 skipped and the same five known
failures (two failed and three errors, all needing the melting SESAM references of Task 11 of the original plan); the
FEniCSx tier goes from 14 to 17 passed. The tested copy is committed at `prototype/work-2026-10-11-moving-mesh-plan/code/`,
and this document is generated from it by `prototype/work-2026-10-11-moving-mesh-plan/make_plan19.py`, which also checks
that the blocks, applied in order to the base, reproduce the tested files exactly. Measurements quoted below were each
made once, on the cases named; they are recorded in Task 8 as facts.

## Global Constraints

- Work in a fresh work copy, `prototype/work-2026-10-12-moving-mesh/code/`, a git clone of
  `prototype/work-2026-10-08-skin/code/` checked out at 49d88e0 (Task 0). It is a git repository of its own, ignored by
  the outer one; commit there after every task. Never edit `prototype/proto3/`, the package `reentry_model/` at the repo
  root, the specs or the monolithic plan.
- Interpreters: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`;
  `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python`, run with
  `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang`.
- Unit tier, run in the work copy: `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py`.
  Baseline at 49d88e0: 303 passed, 1 skipped, and five known failures (`test_reentry_model_compare.py`'s melt plots and
  the melting reference fixtures — two failed, three errors). Record the counts after every task. No other failure is
  acceptable.
- FEniCSx tier: `FI_PROVIDER=tcp CC=... "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q`. Baseline 14 passed.
- Nothing in the existing run path changes behaviour: no CLI flag, no history column, no default. A melting or
  non-melting flight from the work copy must be bit-identical to the base's (Task 8 Step 1 checks one).
- Books: in every test that moves or rebuilds, total enthalpy plus what the corners swept equals the enthalpy before, to
  relative 1e-12; the rebuilt volume equals the old volume to relative 1e-12.
- gmsh runs on one thread in a rebuild: fenicsx_env's gmsh 4.15.2 segfaults meshing a discrete surface's volume with
  HXT on four threads (found while planning, 2026-10-11).
- Exit codes: a move no rebuild can absorb raises `RuntimeError` (exit code 1 when it reaches the CLI in stage 2); a
  surface the validity gate rejects raises `movemesh.RemeshError`, a `RuntimeError`.
- Long runs under `caffeinate -i`, on mains power, the power source logged at start and end (`pmset -g batt`).

## Review Focus

1. **Corners moving outward** — rigid metal joining the 3D model, as freezing will in stage 2: the swept volume is
   negative, the books still close and a uniform field stays uniform. Pinned by
   `test_an_outward_move_books_metal_coming_in` (Task 4).
2. **Corners moving along the surface** — the tangential smoothing every move adds: it sweeps no volume and keeps the
   books. Pinned by `test_a_tangential_move_sweeps_nothing_and_keeps_the_books` (Task 4).
3. **A solver set up again on a rebuilt mesh** — the scikit-fem backend kept its AMG hierarchy across `setup` and
   crashed on the first solve that reused it ("matrix and preconditioner have different shapes", found while
   measuring); FEniCSx must keep exact books across a rebuild too. Pinned by
   `test_a_solver_set_up_on_a_new_mesh_drops_its_old_amg_hierarchy` (Task 2) and
   `test_a_rebuild_keeps_the_energy_on_fenicsx` (Task 5).
4. **Corner identity across repeated rebuilds** — stage 2's columns live on the corners in a fixed order, and gmsh
   renumbers nodes. Pinned by `test_corners_keep_their_order_through_repeated_rebuilds` (Task 5).
5. **A move no rebuild can absorb** — a fold in the surface sharper than one triangle (found while measuring, Task 8):
   a clear model failure, not a hang or an inverted element. Pinned by
   `test_a_move_refused_on_a_fresh_mesh_is_a_model_failure` (Task 6).

---

### Task 0: The work copy

**Files:** none in the code; one line in the outer repository's `.gitignore`.

- [ ] **Step 1: Clone the stage-1 skin work copy at its reviewed head**

```bash
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
git clone work-2026-10-08-skin/code work-2026-10-12-moving-mesh/code
git -C work-2026-10-12-moving-mesh/code checkout -b moving-mesh 49d88e0
```

- [ ] **Step 2: Ignore it in the outer repository** — append to `.gitignore`, under the existing line for
  `prototype/work-2026-10-08-skin/code/`:

```
prototype/work-2026-10-12-moving-mesh/code/
```

- [ ] **Step 3: Record the baseline.** Run the unit tier and the FEniCSx tier in the work copy (commands in Global
  Constraints). Expected: 303 passed, 1 skipped, two failed and three errors (the five known); FEniCSx 14 passed.

### Task 1: Element and surface geometry

**Files:**
- Create: `reentry_model/movemesh.py` (its header and the geometry section)
- Create: `tests/test_reentry_model_movemesh.py` (its header, helpers and the Task 1 tests)

**Interfaces:**
- Produces: `RemeshError(RuntimeError)`; `signed_volumes(points, tets) -> (ne,)`; `radius_ratio(points, tets) -> (ne,)`
  (3 r_in / R_circ, 1 regular); `triangle_quality(points, faces) -> (nf,)` (2 r_in / R_circ); `max_dihedral(points,
  faces) -> float` degrees; `element_gradients(points, tets) -> ((ne, 4, 3), (ne,))`; `nodal_gradient(points, tets, f,
  elements=None) -> (n, 3)`; `lumped_masses(points, tets, rho) -> (n,)`; `orient_outward(points, faces, normals) ->
  faces`; `corner_frame(points, faces, corners) -> (normals (nc, 3), dual_areas (nc,))`; `surface_edges(faces) -> (ne,
  2)`; `tangential_smoothing(points, faces, normal, lam=SMOOTH_LAMBDA, clamp=SMOOTH_CLAMP, edges=None) -> (n, 3)`;
  constants `DEFORM_DEPTH = 15e-3`, `RATIO_MIN = 0.05`, `MIN_GAP = 20`, `ZONE_NEIGHBOURS = 12`, `SMOOTH_LAMBDA = 0.5`,
  `SMOOTH_CLAMP = 0.35`.

- [ ] **Step 1: Write the test file's header, its helpers and the Task 1 tests** (`tests/test_reentry_model_movemesh.py`).
  The helper `inward` uses `MovingMesh.frame` from Task 4; it is defined now so that the file grows only by appending.

```python
"""The moving 3D mesh (layered-skin design §5.3, sub-plan 19): element quality, the deforming zone, the remap and its
books, the rebuild inside a fixed surface, the validity gate, the rebuild policy, and a moving front against an exact
solution."""
import numpy as np
import pytest
from scipy.special import erfc

from reentry_model import material as material_mod
from reentry_model import mesh as mesh_mod
from reentry_model import movemesh
from reentry_model.thermal import thermal_solver


# -- element and surface geometry (Task 1) -------------------------------------------------------------------------
def constant_material(k=30.0, rho=3000.0, cp=1000.0):
    """A non-melting material with constant properties (diffusivity k / (rho cp) = 1e-5 m2/s by default)."""
    return material_mod.Material("constant", rho, 0.0, np.array([0.0, 5000.0]), np.array([cp, cp]),
                                 np.array([0.0, 5000.0]), np.array([k, k]))


def solver_on(mesh, mat, T=300.0, name="skfem"):
    solver = thermal_solver(name, linear_solver="direct")
    solver.setup(mesh, mat, 0.0)
    solver.set_temperature(T)
    return solver


def inward(moving, depth):
    """Corner displacements moving every corner `depth` metres along its inward normal."""
    normal, _ = moving.frame()
    return -depth * normal


def test_radius_ratio_is_one_for_a_regular_tet_and_zero_for_a_flat_one():
    regular = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    flat = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], dtype=float)
    tets = np.array([[0, 1, 2, 3]])
    assert movemesh.radius_ratio(regular, tets)[0] == pytest.approx(1.0, abs=1e-12)
    assert movemesh.radius_ratio(flat, tets)[0] == pytest.approx(0.0, abs=1e-12)


def test_triangle_quality_and_fold_of_a_flat_and_a_folded_pair():
    flat = np.array([[0, 0, 0], [1, 0, 0], [0.5, np.sqrt(3) / 2, 0], [0.5, -np.sqrt(3) / 2, 0]], dtype=float)
    faces = np.array([[0, 1, 2], [1, 0, 3]])
    assert np.allclose(movemesh.triangle_quality(flat, faces), 1.0)
    assert movemesh.max_dihedral(flat, faces) == pytest.approx(0.0, abs=1e-6)
    folded = flat.copy()
    folded[3] = [0.5, 0.0, -np.sqrt(3) / 2]                   # the second triangle turned down by 90 degrees
    assert movemesh.max_dihedral(folded, faces) == pytest.approx(90.0, abs=1e-6)


def test_nodal_gradient_is_exact_for_a_linear_field():
    m = mesh_mod.box_mesh(4e-3, 2e-3, 2e-3, 1e-3)
    f = 3.0 + m.points @ np.array([2.0, -1.0, 0.5])
    g = movemesh.nodal_gradient(m.points, m.tets, f)
    assert np.allclose(g, [2.0, -1.0, 0.5], atol=1e-9)


def test_corner_frame_gives_outward_normals_and_dual_areas_summing_to_the_surface():
    m = mesh_mod.box_mesh(4e-3, 2e-3, 2e-3, 1e-3)
    s = m.surface()
    corners = np.unique(s.faces)
    faces = movemesh.orient_outward(m.points, s.faces, s.normals)
    normal, area = movemesh.corner_frame(m.points, faces, corners)
    assert area.sum() == pytest.approx(s.areas.sum(), rel=1e-12)
    centre = m.points.mean(axis=0)
    assert np.all(np.einsum("ij,ij->i", normal, m.points[corners] - centre) > 0.0)


def test_tangential_smoothing_moves_along_the_surface_and_evens_it_out(coarse_sphere_mesh):
    m = coarse_sphere_mesh
    s = m.surface()
    corners = np.unique(s.faces)
    local = np.full(m.n_nodes, -1)
    local[corners] = np.arange(len(corners))
    faces = local[movemesh.orient_outward(m.points, s.faces, s.normals)]
    pts = m.points[corners].copy()
    rng = np.random.default_rng(5)
    r = np.linalg.norm(pts, axis=1, keepdims=True)
    jitter = rng.normal(0.0, 0.5e-3, pts.shape)
    jitter -= np.einsum("ij,ij->i", jitter, pts / r)[:, None] * pts / r      # along the sphere only
    pts = (pts + jitter) * (0.05 / np.linalg.norm(pts + jitter, axis=1, keepdims=True))
    normal = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    q0 = movemesh.triangle_quality(pts, faces).min()
    for _ in range(10):
        step = movemesh.tangential_smoothing(pts, faces, normal)
        assert np.abs(np.einsum("ij,ij->i", step, normal)).max() < 1e-15
        pts = pts + step
        pts *= 0.05 / np.linalg.norm(pts, axis=1, keepdims=True)
        normal = pts / 0.05
    assert movemesh.triangle_quality(pts, faces).min() > q0


```

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q`
Expected: collection error, `ImportError: cannot import name 'movemesh'`.

- [ ] **Step 3: Create `reentry_model/movemesh.py` with its header and the geometry section**

```python
"""The moving 3D mesh (layered-skin design, docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md §5.3).

The 3D mesh's surface nodes are the corners whose columns carry the skins; they move by the metal crossing the rigid
front, and the 3D nodes within a depth D below the surface (the deforming zone) move with them, each by a weighted
average of the movements of the surface nodes near it, the weight falling linearly from one at the surface to zero at
depth D. After every move the temperature field is remapped onto the moved nodes and the energy books are closed
exactly; before any element degrades the mesh is rebuilt inside the same surface triangulation.

The heat equation's moving-mesh term is applied as a remap at each move (the operator-split, "Lagrangian plus remap"
form of the arbitrary Lagrangian-Eulerian method): the nodal enthalpy is carried to each moved node along its
displacement by the nodal gradient, h_i' = h_i + grad(h)_i . u_i, which is exact for uniform and linear fields; then
the total enthalpy is set to its old value less what left through the moving surface, by one uniform enthalpy offset on
the moved nodes, so the discrete geometric conservation law holds (a uniform field stays exactly uniform) and the books
close to round-off. The conduction step that follows runs on the new geometry with the backends' fixed-mesh scheme.

The corners also move along the surface: `tangential_smoothing` nudges each toward the mean of its neighbours in its
tangent plane (design section 5.1), which keeps the fixed surface triangulation even. Measured while planning
(sub-plan 19, 2026-10-11) on the 100 mm sphere receding at its nose by 0.5 mm per step in proportion to how directly
each corner faces the flow: without it the run failed at 20.5 mm of nose recession, with it at 27.0 mm. What stopped it
then was the shape, not the triangles: the angle between neighbouring surface triangles (`max_dihedral`) grew from
4 to 82 degrees as that motion sharpened the shoulder, and elements whose four nodes all lie on a fold that sharp are
flattened by their own corners' motion, which no rebuild inside the same surface can undo. The fold is reported.

What left through the surface is booked per corner: the volume each corner sweeps is its dual area times its normal
displacement (positive inward), rescaled so the corners' volumes add up exactly to the change of the mesh's volume, and
it carries the corner's own nodal enthalpy. The layered skins take exactly that mass and enthalpy (stage 2).

Rebuild: gmsh meshes the volume inside the current surface triangulation without touching it (a discrete surface in
a geo surface loop, generate(3) only; checked 2026-10-11: every surface triangle and corner is kept to the last bit,
though gmsh renumbers nodes, so the corners are matched back by their coordinates). The enthalpy is interpolated onto
the new nodes in the old mesh and its total made exact by a uniform offset, as the 2026-09-27 remeshing design does;
the pre-correction residual is reported with every rebuild. The validity gate stops on a non-manifold surface, more
than one connected body, or a non-positive enclosed volume."""
import time
from dataclasses import dataclass

import numpy as np

from . import mesh as mesh_mod

DEFORM_DEPTH = 15.0e-3       # m, depth D of the deforming zone: the dense band's depth (mesh.DEFAULT_BAND)
RATIO_MIN = 0.05             # radius ratio below which a move is refused and the mesh rebuilt (2026-09-27 design §7.1)
MIN_GAP = 20                 # macro steps between rebuilds triggered by compression (a refused move ignores it)
ZONE_NEIGHBOURS = 12         # surface nodes averaged for each interior node of the zone
SMOOTH_LAMBDA = 0.5          # tangential smoothing: share of the way to the neighbours' mean per move
SMOOTH_CLAMP = 0.35          # of the mean surface edge: the most a corner is nudged per move (2026-09-27 design)
_REF_GRADS = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])


class RemeshError(RuntimeError):
    """The surface cannot be rebuilt: non-manifold, more than one body, or a non-positive enclosed volume."""


# -- element and surface geometry ----------------------------------------------------------------------------------
def signed_volumes(points, tets):
    x = points[tets]
    return np.einsum("ij,ij->i", x[:, 1] - x[:, 0], np.cross(x[:, 2] - x[:, 0], x[:, 3] - x[:, 0])) / 6.0


def radius_ratio(points, tets):
    """3 r_in / R_circ per tetrahedron: 1 for a regular one, 0 for a flat one."""
    x = points[tets]
    u, v, w = x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]
    triple = np.einsum("ij,ij->i", u, np.cross(v, w))
    vol = np.abs(triple) / 6.0
    faces = (np.cross(v - u, w - u), np.cross(w, v), np.cross(u, w), np.cross(v, u))
    area = 0.5 * sum(np.linalg.norm(f, axis=1) for f in faces)
    r_in = 3.0 * vol / np.where(area > 0.0, area, 1.0)
    c = ((u * u).sum(1)[:, None] * np.cross(v, w) + (v * v).sum(1)[:, None] * np.cross(w, u)
         + (w * w).sum(1)[:, None] * np.cross(u, v))
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.linalg.norm(c, axis=1) / np.abs(2.0 * triple)
        ratio = np.where(vol > 0.0, 3.0 * r_in / R, 0.0)
    return np.nan_to_num(ratio, nan=0.0)


def triangle_quality(points, faces):
    """2 r_in / R_circ per triangle: 1 for an equilateral one, 0 for a flat one."""
    a, b, c = (points[faces[:, i]] for i in range(3))
    la, lb, lc = np.linalg.norm(b - c, axis=1), np.linalg.norm(c - a, axis=1), np.linalg.norm(a - b, axis=1)
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    s = 0.5 * (la + lb + lc)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = 2.0 * (area / s) / (la * lb * lc / (4.0 * area))
    return np.nan_to_num(q, nan=0.0)


def max_dihedral(points, faces):
    """The largest angle [degrees] between the normals of two surface triangles sharing an edge: how sharply the
    surface folds at the scale of one triangle."""
    a, b, c = (points[faces[:, i]] for i in range(3))
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1)[:, None]
    edges = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(faces)), 3)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    e, o = edges[order], owner[order]
    pair = np.all(e[1:] == e[:-1], axis=1)
    cos = np.einsum("ij,ij->i", n[o[1:][pair]], n[o[:-1][pair]])
    return float(np.degrees(np.arccos(np.clip(cos.min(), -1.0, 1.0)))) if cos.size else 0.0


def element_gradients(points, tets):
    """Gradients of the four P1 basis functions in every tetrahedron, (ne, 4, 3), and the volumes (ne,)."""
    x = points[tets]
    J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
    vol = np.abs(np.linalg.det(J)) / 6.0
    return np.einsum("eij,aj->eai", np.linalg.inv(J).transpose(0, 2, 1), _REF_GRADS), vol


def nodal_gradient(points, tets, f, elements=None):
    """Volume-weighted average over each node's elements of grad f, f a nodal P1 field: exact for a linear f. With
    `elements`, only those elements are used: the gradient is then right at every node all of whose elements are in
    the list (the moving nodes, when the list is the elements they touch)."""
    if elements is not None:
        tets = tets[elements]
    grads, vol = element_gradients(points, tets)
    g = np.einsum("ea,eai->ei", f[tets], grads) * vol[:, None]
    out = np.zeros((len(points), 3))
    for a in range(4):
        np.add.at(out, tets[:, a], g)
    w = np.bincount(tets.ravel(), weights=np.repeat(vol, 4), minlength=len(points))
    return out / np.where(w > 0.0, w, 1.0)[:, None]


def lumped_masses(points, tets, rho):
    """rho times a quarter of each incident element's volume, per node: the backends' lumped capacity weights."""
    vol = np.abs(signed_volumes(points, tets))
    return rho * np.bincount(tets.ravel(), weights=np.repeat(vol / 4.0, 4), minlength=len(points))


def orient_outward(points, faces, normals):
    """The triangles reordered where needed so that (b - a) x (c - a) points along their outward normals."""
    faces = np.array(faces, dtype=np.int64)
    a, b, c = (points[faces[:, i]] for i in range(3))
    flip = np.einsum("ij,ij->i", np.cross(b - a, c - a), normals) < 0.0
    faces[flip] = faces[flip][:, [0, 2, 1]]
    return faces


def corner_frame(points, faces, corners):
    """Area-weighted outward vertex normals and dual areas (a third of each touching triangle) at `corners`."""
    a, b, c = (points[faces[:, i]] for i in range(3))
    n2 = np.cross(b - a, c - a)                         # twice the area times the unit normal (outward for mesh faces)
    normal = np.zeros((len(points), 3))
    area = np.zeros(len(points))
    for i in range(3):
        np.add.at(normal, faces[:, i], n2)
        np.add.at(area, faces[:, i], np.linalg.norm(n2, axis=1) / 6.0)
    normal = normal[corners]
    return normal / np.linalg.norm(normal, axis=1)[:, None], area[corners]


def surface_edges(faces):
    """The unique edges of a triangulation, (ne, 2), each with its smaller node id first."""
    return np.unique(np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1), axis=0)


def tangential_smoothing(points, faces, normal, lam=SMOOTH_LAMBDA, clamp=SMOOTH_CLAMP, edges=None):
    """Displacements that even out the surface triangles without moving the surface (layered-skin design section 5.1):
    each corner moves `lam` of the way toward the mean of its neighbours, projected onto its tangent plane, the move
    clamped to `clamp` mean edge lengths. Tangential moves sweep no volume to first order. `edges` (from
    `surface_edges(faces)`) saves recomputing them every move."""
    if edges is None:
        edges = surface_edges(faces)
    n = len(points)
    total = np.zeros((n, 3))
    np.add.at(total, edges[:, 0], points[edges[:, 1]])
    np.add.at(total, edges[:, 1], points[edges[:, 0]])
    degree = np.bincount(edges.ravel(), minlength=n)
    L = total / np.maximum(degree, 1)[:, None] - points
    L -= np.einsum("ij,ij->i", L, normal)[:, None] * normal
    step = lam * L
    limit = clamp * np.linalg.norm(points[edges[:, 0]] - points[edges[:, 1]], axis=1).mean()
    size = np.linalg.norm(step, axis=1)
    return step * np.minimum(1.0, limit / np.where(size > 0.0, size, 1.0))[:, None]



```

- [ ] **Step 4: Run them to see them pass**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q`
Expected: 5 passed.

- [ ] **Step 5: Unit tier and commit.** Expected 308 passed, 1 skipped, the five known.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py
git commit -m "Moving mesh: element and surface geometry helpers"
```

### Task 2: The backends on moved geometry

**Files:**
- Modify: `reentry_model/thermal/skfem_backend.py` (`setup`, new `update_geometry`)
- Modify: `reentry_model/thermal/fenicsx_backend.py` (`setup`, new `update_geometry`)
- Test: `tests/test_reentry_model_movemesh.py` (append), `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Produces: `solver.update_geometry(elements=None)` on both backends — the mesh's `points` have moved in place;
  scikit-fem recomputes the element matrices of `elements` (all if None) and the surface; FEniCSx copies the coordinates
  into dolfinx's geometry and refreshes the surface (it recomputes every Jacobian at assembly anyway). `setup` on a new
  mesh now discards the scikit-fem AMG hierarchy.

- [ ] **Step 1: Append the Task 2 tests** to `tests/test_reentry_model_movemesh.py`:

```python
# -- the backends on moved geometry (Task 2) -----------------------------------------------------------------------
def test_moved_geometry_gives_the_operators_of_a_fresh_solver(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat, 450.0)
    r = np.linalg.norm(m.points, axis=1)
    m.points *= (1.0 - 0.04 * (r / 0.05) ** 4)[:, None]       # an uneven squeeze toward the centre, in place
    outer = np.flatnonzero((r[m.tets] > 0.02).any(axis=1))    # the elements whose nodes moved noticeably...
    solver.update_geometry(np.arange(m.n_elements))           # ...but every node moved, so recompute all
    fresh = solver_on(mesh_mod.VolumeMesh(m.points.copy(), m.tets.copy()), mat, 450.0)
    T = solver.temperature()
    K1, M1 = solver.operators(T)
    K2, M2 = fresh.operators(T)
    assert outer.size and abs(K1 - K2).max() == 0.0 and abs(M1 - M2).max() == 0.0
    assert np.array_equal(solver.areas, fresh.areas)


def test_a_solver_set_up_on_a_new_mesh_drops_its_old_amg_hierarchy(coarse_sphere_mesh, uniform_test_mesh):
    mat = constant_material()
    solver = thermal_solver("skfem")                          # AMG, the default
    solver.setup(coarse_sphere_mesh, mat, 0.0)
    solver.set_temperature(300.0)
    solver.step(0.5, np.full(len(solver.faces), 1e5), 300.0)
    solver.setup(uniform_test_mesh, mat, 0.0)                 # what a rebuild does
    solver.set_temperature(300.0)
    solver.step(0.5, np.full(len(solver.faces), 1e5), 300.0)  # raised "different shapes" before the fix


```

  and to the end of `tests/test_reentry_model_fenicsx.py`:

```python


def test_moved_geometry_steps_like_a_fresh_solver(coarse_sphere_mesh):
    """Moving mesh (sub-plan 19): after update_geometry, a step equals a fresh solver's on the moved mesh."""
    from reentry_model import mesh as mesh_mod
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    s1 = solver("fenicsx", m, mat, linear_solver="direct")
    s1.set_temperature(500.0 + 1000.0 * m.points[:, 0])
    r = np.linalg.norm(m.points, axis=1)
    m.points *= (1.0 - 0.04 * (r / R) ** 4)[:, None]
    s1.update_geometry()
    s2 = solver("fenicsx", mesh_mod.VolumeMesh(m.points.copy(), m.tets.copy()), mat, linear_solver="direct")
    s2.set_temperature(s1.temperature())
    q = np.full(len(s1.faces), 2e5)
    r1, r2 = s1.step(0.5, q, 300.0), s2.step(0.5, q, 300.0)
    assert np.abs(r1.T - r2.T).max() < 1e-9


```

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q -k "fresh_solver or amg"`
Expected: `AttributeError: 'SkfemThermalSolver' object has no attribute 'update_geometry'`, and the AMG test fails with
`ValueError: matrix and preconditioner have different shapes`.

- [ ] **Step 3: Apply the scikit-fem change**

```diff
--- a/reentry_model/thermal/skfem_backend.py
+++ b/reentry_model/thermal/skfem_backend.py
@@ -92,9 +92,21 @@
         n = len(self.points)
         self.pattern = _Pattern(np.repeat(self.tets, 4, axis=1).ravel(), np.tile(self.tets, (1, 4)).ravel(), n)
         self.T, self._T_prev = np.full(n, 300.0), None
+        self._ml = None                           # an AMG hierarchy belongs to one mesh: a rebuilt mesh needs a new one
         self.phi = np.where(mesh.active, 1.0, 0.0)
         self.film_mass = np.zeros(n)              # kg per node: the melt film riding on the boundary (Step 3)
         self.pinned = np.zeros(n, dtype=bool)
+        self._refresh_surface()
+
+    def update_geometry(self, elements=None):
+        """The mesh's node coordinates have moved (the moving 3D mesh, layered-skin design section 5.3): recompute the
+        element matrices of `elements` (all if None) and the surface. `self.points` is the mesh's own array, so the
+        coordinates are already here."""
+        if elements is None:
+            self.vol, self.Ke, self.Mk = element_matrices(self.points, self.tets)
+        else:
+            elements = np.asarray(elements, dtype=np.int64)
+            self.vol[elements], self.Ke[elements], self.Mk[elements] = element_matrices(self.points, self.tets[elements])
         self._refresh_surface()
 
     def _refresh_surface(self):
```

- [ ] **Step 4: Apply the FEniCSx change**

```diff
--- a/reentry_model/thermal/fenicsx_backend.py
+++ b/reentry_model/thermal/fenicsx_backend.py
@@ -66,6 +66,8 @@
         self.cell_nodes = self.node_of_dof[self.cell_dofs]                   # mesh node ids of each dolfinx cell
         centroids = mesh.points[mesh.tets].mean(axis=1)
         self.element_of_cell = cKDTree(centroids).query(mesh.points[self.cell_nodes].mean(axis=1))[1]
+        # the VolumeMesh node behind each dolfinx geometry node, so that moved coordinates can be copied in
+        self._geometry_nodes = np.asarray(self.msh.geometry.input_global_indices, dtype=np.int64)
         self.coef = fem.Function(self.V0)                                   # phi_e k(T_e) per cell
         u, v = ufl.TrialFunction(self.V), ufl.TestFunction(self.V)
         self.a = fem.form(self.coef * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx)
@@ -88,6 +90,13 @@
         self._refresh_surface()
 
     # -- geometry / fractions -----------------------------------------------------------------------------------
+    def update_geometry(self, elements=None):
+        """The mesh's node coordinates have moved (the moving 3D mesh, layered-skin design section 5.3): copy them into
+        dolfinx's geometry, which every assembly reads, and refresh the surface. The nodal pieces read the VolumeMesh.
+        `elements` is accepted for the protocol's sake: dolfinx recomputes every Jacobian at assembly anyway."""
+        self.msh.geometry.x[:, :3] = self.mesh.points[self._geometry_nodes]
+        self._refresh_surface()
+
     def _refresh_surface(self):
         surface = self.mesh.surface()
         self.faces, self.areas = surface.faces, surface.areas
```

- [ ] **Step 5: Run them to see them pass**

Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q` — expected 7 passed — and
`FI_PROVIDER=tcp CC=... "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q -k moved_geometry` — expected 1 passed.

- [ ] **Step 6: Both tiers and commit.** Expected unit tier 310 passed, 1 skipped, the five known; FEniCSx 15 passed.

```bash
git add reentry_model/thermal/skfem_backend.py reentry_model/thermal/fenicsx_backend.py tests/test_reentry_model_movemesh.py tests/test_reentry_model_fenicsx.py
git commit -m "Thermal backends follow moved node coordinates, and a solver set up on a new mesh drops its old AMG hierarchy"
```

### Task 3: The deforming zone

**Files:**
- Modify: `reentry_model/movemesh.py` (append the zone section)
- Test: `tests/test_reentry_model_movemesh.py` (append)

**Interfaces:**
- Consumes: Task 1's `corner_frame`, `orient_outward` (tests).
- Produces: `ZoneMap(points, corners, depth=DEFORM_DEPTH, neighbours=ZONE_NEIGHBOURS)` with `.displacement(u_corners
  (nc, 3)) -> (n, 3)`, `.interior` (zone node ids that are not corners), `.corners`, `.depth`.

- [ ] **Step 1: Append the Task 3 tests**

```python
# -- the deforming zone (Task 3) -----------------------------------------------------------------------------------
def test_zone_moves_corners_exactly_and_leaves_deep_nodes_still(coarse_sphere_mesh):
    m = coarse_sphere_mesh
    corners = np.unique(m.surface().faces)
    zone = movemesh.ZoneMap(m.points, corners, depth=10e-3)
    u = np.tile([1e-4, 0.0, 0.0], (len(corners), 1))
    U = zone.displacement(u)
    assert np.array_equal(U[corners], u)
    deep = np.linalg.norm(m.points, axis=1) < 0.05 - 10.5e-3
    assert np.all(U[deep] == 0.0)
    shallow = np.setdiff1d(zone.interior, corners)
    assert np.all(U[shallow, 0] > 0.0) and np.all(U[shallow, 0] <= 1e-4 + 1e-15)


def test_uniform_recession_compresses_the_zone_uniformly(coarse_sphere_mesh):
    m = coarse_sphere_mesh
    s = m.surface()
    corners = np.unique(s.faces)
    normal, _ = movemesh.corner_frame(m.points, movemesh.orient_outward(m.points, s.faces, s.normals), corners)
    zone = movemesh.ZoneMap(m.points, corners, depth=20e-3)
    moved = m.points + zone.displacement(-2e-3 * normal)
    r0, r1 = np.linalg.norm(m.points, axis=1), np.linalg.norm(moved, axis=1)
    near = (r0 > 0.05 - 8e-3) & (r0 < 0.05 - 1e-3)      # well inside the zone, away from the surface
    # each node's depth below the surface shrinks by about the zone's factor (D - r) / D = 0.9
    assert np.median(((0.05 - 2e-3) - r1[near]) / (0.05 - r0[near])) == pytest.approx(0.9, abs=0.03)


```

- [ ] **Step 2: Run them to see them fail** — Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q -k zone`.
  Expected: `AttributeError: module 'reentry_model.movemesh' has no attribute 'ZoneMap'`.

- [ ] **Step 3: Append the zone section to `reentry_model/movemesh.py`**

```python
# -- the deforming zone ----------------------------------------------------------------------------------------------
class ZoneMap:
    """Which 3D nodes move with the surface, and how: surface nodes by their own displacement; every other node within
    `depth` of the surface by w_i times an inverse-square-distance average of its ZONE_NEIGHBOURS nearest surface
    nodes' displacements, w_i = 1 - d_i / depth (d_i its distance to the nearest surface node). Built once per mesh."""

    def __init__(self, points, corners, depth=DEFORM_DEPTH, neighbours=ZONE_NEIGHBOURS):
        from scipy.spatial import cKDTree
        if depth <= 0.0:
            raise ValueError("the deforming zone needs a positive depth")
        self.n_points, self.corners, self.depth = len(points), np.asarray(corners, dtype=np.int64), float(depth)
        tree = cKDTree(points[self.corners])
        d = tree.query(points)[0]
        w = np.clip(1.0 - d / self.depth, 0.0, 1.0)
        is_corner = np.zeros(len(points), dtype=bool)
        is_corner[self.corners] = True
        self.interior = np.flatnonzero((w > 0.0) & ~is_corner)
        k = min(neighbours, len(self.corners))
        dist, idx = tree.query(points[self.interior], k=k)
        dist, idx = dist.reshape(len(self.interior), k), idx.reshape(len(self.interior), k)
        K = 1.0 / (dist + 1e-3 * self.depth) ** 2
        self.idx, self.K = idx, K / K.sum(axis=1)[:, None]
        self.w = w[self.interior]

    def displacement(self, u_corners):
        """(n_corners, 3) corner displacements -> (n_points, 3) node displacements."""
        u = np.asarray(u_corners, dtype=float)
        U = np.zeros((self.n_points, 3))
        U[self.corners] = u
        U[self.interior] = self.w[:, None] * np.einsum("ik,ikj->ij", self.K, u[self.idx])
        return U



```

- [ ] **Step 4: Run them to see them pass** — expected 9 passed in the file.

- [ ] **Step 5: Unit tier and commit.** Expected 312 passed, 1 skipped, the five known.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py
git commit -m "Moving mesh: the deforming zone moves the nodes under the corners in proportion to their depth"
```

### Task 4: Moving the mesh — the remap and its books

**Files:**
- Modify: `reentry_model/movemesh.py` (append `MoveResult` and the class `MovingMesh` up to and including `move`)
- Test: `tests/test_reentry_model_movemesh.py` (append)

**Interfaces:**
- Consumes: Tasks 1–3; the solver's `temperature()`, `set_temperature()`, `energy()`, `update_geometry(elements)`,
  `emissivity`; the material's `enthalpy`, `temperature_from_enthalpy`, `rho`.
- Produces: `MoveResult(refused, swept_volume (nc,), swept_enthalpy (nc,), correction, worst_ratio, rebuilt=False)`;
  `MovingMesh(mesh, solver, material, depth=DEFORM_DEPTH, ratio_min=RATIO_MIN, min_gap=MIN_GAP, sizes=None)` with
  `.mesh`, `.solver`, `.material`, `.corners` (current node id of each corner, fixed order), `.faces_corner` (surface
  triangles in corner numbering, outward), `.edges_corner`, `.zone`, `.compression` (inward motion of each corner since
  the last rebuild), `.steps_since_rebuild`, `.rebuilds`, `.last_rebuild`; methods `frame() -> (normals, dual_areas)`,
  `smoothing(normal=None, lam=None, clamp=None) -> (nc, 3)`, `surface_quality() -> float`, `surface_fold() -> float`
  degrees, `move(u_corners) -> MoveResult`. Stage 2 credits each column with `swept_volume` and `swept_enthalpy`.

- [ ] **Step 1: Append the Task 4 tests**

```python
# -- moving the mesh: the remap and its books (Task 4) -------------------------------------------------------------
def test_a_uniform_field_stays_exactly_uniform_while_the_mesh_moves(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat, 700.0)
    moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3)
    rng = np.random.default_rng(3)
    for _ in range(5):
        r = moving.move(inward(moving, 2e-4) * rng.uniform(0.0, 1.0, (len(moving.corners), 1)))
        assert not r.refused
    assert np.abs(solver.temperature() - 700.0).max() < 1e-9


def test_the_books_close_what_the_corners_swept_left_the_body(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat)
    r_node = np.linalg.norm(m.points, axis=1)
    solver.set_temperature(300.0 + 4000.0 * (r_node - 0.03).clip(0.0))           # warm near the surface
    moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3)
    E0, V0 = solver.energy(), m.volume()
    r = moving.move(inward(moving, 3e-4))
    assert r.swept_volume.sum() == pytest.approx(V0 - moving.mesh.volume(), rel=1e-12)
    assert solver.energy() + r.swept_enthalpy.sum() == pytest.approx(E0, rel=1e-12)
    # the offset that closed the books is reported; on this kinked field over 4 mm cells it is 0.7 % of what the
    # corners swept (measured while planning), the lumped quadrature's error on a field that is not linear
    assert 0.0 < abs(r.correction) < 0.02 * r.swept_enthalpy.sum()


def test_a_linear_field_is_carried_to_the_moved_nodes(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat)
    T_of = lambda p: 500.0 + 2000.0 * p[:, 0]
    solver.set_temperature(T_of(m.points))
    moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3)
    r = moving.move(inward(moving, 2e-4))
    # the remap is exact for a linear field; only the books' uniform offset moves it: 3.3e-6 K here, from a 0.0035 J
    # offset against 3877 J swept (measured while planning)
    assert np.abs(solver.temperature() - T_of(moving.mesh.points)).max() < 1e-4
    assert abs(r.correction) < 1e-5 * r.swept_enthalpy.sum()


def test_a_move_that_would_invert_an_element_is_refused_and_changes_nothing(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat, 400.0)
    moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3)
    before, T_before = m.points.copy(), solver.temperature()
    u = np.zeros((len(moving.corners), 3))
    u[0] = -8e-3 * moving.frame()[0][0]                                           # one corner 8 mm in: past its elements
    r = moving.move(u)
    assert r.refused
    assert np.array_equal(moving.mesh.points, before) and np.array_equal(solver.temperature(), T_before)


def test_an_outward_move_books_metal_coming_in(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat, 500.0)
    moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3)
    E0, V0 = solver.energy(), m.volume()
    r = moving.move(inward(moving, -2e-4))                    # outward: rigid metal joining the 3D model (stage 2)
    assert not r.refused and r.swept_volume.sum() < 0.0
    assert r.swept_volume.sum() == pytest.approx(V0 - moving.mesh.volume(), rel=1e-12)
    assert solver.energy() + r.swept_enthalpy.sum() == pytest.approx(E0, rel=1e-12)
    assert np.abs(solver.temperature() - 500.0).max() < 1e-9


def test_a_tangential_move_sweeps_nothing_and_keeps_the_books(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat)
    solver.set_temperature(300.0 + 2000.0 * m.points[:, 2])
    moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3)
    normal, _ = moving.frame()
    tangent = np.cross(normal, [0.0, 0.0, 1.0])               # along the surface, around the z axis
    E0 = solver.energy()
    r = moving.move(1e-4 * tangent)
    assert abs(r.swept_volume.sum()) < 1e-3 * 1e-4 * moving.frame()[1].sum()
    assert solver.energy() + r.swept_enthalpy.sum() == pytest.approx(E0, rel=1e-12)


def test_the_moving_mesh_reports_its_surface_quality_and_fold(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    moving = movemesh.MovingMesh(m, solver_on(m, mat), mat)
    # the coarse 4 mm sphere as gmsh made it: worst triangle 0.14, largest fold 14.4 degrees (measured while planning)
    assert 0.1 < moving.surface_quality() <= 1.0
    assert 5.0 < moving.surface_fold() < 20.0
    assert np.linalg.norm(moving.smoothing(), axis=1).max() > 0.0


```

- [ ] **Step 2: Run them to see them fail** — Run: `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q`.
  Expected: the Task 4 tests fail with `AttributeError: module 'reentry_model.movemesh' has no attribute 'MovingMesh'`.

- [ ] **Step 3: Append `MoveResult` and `MovingMesh` (through `move`) to `reentry_model/movemesh.py`**

```python
# -- moving the mesh -------------------------------------------------------------------------------------------------
@dataclass
class MoveResult:
    refused: bool                       # the move would invert or degrade an element: nothing was changed
    swept_volume: np.ndarray            # (n_corners,) m3 each corner's motion took out of the 3D model (negative: added)
    swept_enthalpy: np.ndarray          # (n_corners,) J that volume carried
    correction: float                   # J, the uniform offset that closed the books (pre-correction residual)
    worst_ratio: float                  # smallest radius ratio among the elements the move touched
    rebuilt: bool = False               # move_with_policy rebuilt the mesh before (or to make) this move


class MovingMesh:
    """A VolumeMesh, its thermal solver and material, moved by its corners (the surface nodes, in a fixed order that
    survives rebuilds) and rebuilt when needed. `mesh` and the solver's mesh are replaced at every rebuild; callers
    read `self.mesh` afresh. `corners` holds the current node id of each corner."""

    def __init__(self, mesh, solver, material, depth=DEFORM_DEPTH, ratio_min=RATIO_MIN, min_gap=MIN_GAP,
                 sizes=None):
        self.mesh, self.solver, self.material = mesh, solver, material
        self.depth, self.ratio_min, self.min_gap = float(depth), float(ratio_min), int(min_gap)
        self.sizes = sizes if sizes is not None else default_sizes(mesh)
        s = mesh.surface()
        self.corners = np.unique(s.faces)
        self._set_surface(orient_outward(mesh.points, s.faces, s.normals))
        self.compression = np.zeros(len(self.corners))      # inward normal motion of each corner since the last rebuild
        self.steps_since_rebuild, self.rebuilds, self.last_rebuild = 0, 0, None

    def _set_surface(self, faces_nodes):
        local = np.full(self.mesh.n_nodes, -1, dtype=np.int64)
        local[self.corners] = np.arange(len(self.corners))
        self.faces_corner = local[faces_nodes]              # surface triangles in corner numbering (rebuild-proof)
        self.edges_corner = surface_edges(self.faces_corner)
        self.zone = ZoneMap(self.mesh.points, self.corners, self.depth)

    def frame(self):
        """Outward normals and dual areas at the corners, on the current geometry."""
        return corner_frame(self.mesh.points, self.corners[self.faces_corner], self.corners)

    def smoothing(self, normal=None, lam=None, clamp=None):
        """Tangential displacements of the corners that even out the surface triangles (add them to a move)."""
        if normal is None:
            normal = self.frame()[0]
        return tangential_smoothing(self.mesh.points[self.corners], self.faces_corner, normal,
                                    SMOOTH_LAMBDA if lam is None else lam, SMOOTH_CLAMP if clamp is None else clamp,
                                    self.edges_corner)

    def surface_quality(self):
        """The worst surface triangle's quality, 2 r_in / R_circ."""
        return float(triangle_quality(self.mesh.points[self.corners], self.faces_corner).min())

    def surface_fold(self):
        """The largest angle between neighbouring surface triangles' normals [degrees]."""
        return max_dihedral(self.mesh.points[self.corners], self.faces_corner)

    def move(self, u_corners):
        """Move the corners by `u_corners` (n_corners, 3), the zone with them; remap the enthalpy; close the books.
        Refuses, changing nothing, if any element the move touches would invert or fall below the radius-ratio floor."""
        mesh, solver, mat = self.mesh, self.solver, self.material
        u = np.asarray(u_corners, dtype=float)
        U = self.zone.displacement(u)
        is_moving = np.any(U != 0.0, axis=1)
        moving = np.flatnonzero(is_moving)
        empty = np.zeros(len(self.corners))
        if moving.size == 0:
            return MoveResult(False, empty, empty.copy(), 0.0, 1.0)
        touched = np.flatnonzero(is_moving[mesh.tets].any(axis=1))
        new_points = mesh.points + U
        vol_new = signed_volumes(new_points, mesh.tets[touched])
        sign = np.sign(signed_volumes(mesh.points, mesh.tets[touched]))
        ratio = radius_ratio(new_points, mesh.tets[touched])
        worst = float(ratio.min())
        if np.any(vol_new * sign <= 0.0) or worst < self.ratio_min:
            return MoveResult(True, empty, empty.copy(), 0.0, worst)
        h = mat.enthalpy(solver.temperature())
        h_new = h + np.einsum("ij,ij->i", nodal_gradient(mesh.points, mesh.tets, h, touched), U)
        E_old = solver.energy()
        normal, area = self.frame()
        dV = -area * np.einsum("ij,ij->i", normal, u)       # volume each corner sweeps, positive inward
        V_old = float(np.abs(signed_volumes(mesh.points, mesh.tets[touched])).sum())     # only touched elements change
        mesh.points[:] = new_points
        mesh.derived_points[:] = new_points
        dV_total = V_old - float(np.abs(vol_new).sum())
        if dV.sum() != 0.0:
            dV = dV * (dV_total / dV.sum())
        swept_h = mat.rho * dV * h[self.corners]
        solver.update_geometry(touched)
        solver.set_temperature(mat.temperature_from_enthalpy(h_new))
        correction = (E_old - swept_h.sum()) - solver.energy()
        m = lumped_masses(mesh.points, mesh.tets[touched], mat.rho)      # complete at the moving nodes
        h_new[moving] += correction / m[moving].sum()
        solver.set_temperature(mat.temperature_from_enthalpy(h_new))
        self.compression += np.maximum(-np.einsum("ij,ij->i", normal, u), 0.0)
        self.steps_since_rebuild += 1
        return MoveResult(False, dV, swept_h, float(correction), worst)

```

- [ ] **Step 4: Run them to see them pass** — expected 16 passed in the file. The measured sizes behind the tolerances,
  each from one run while planning: on the steep, kinked field of the books test over 4 mm cells the offset that closes
  the books is 0.7 % of what the corners swept; on a linear field it is 0.0035 J against 3877 J swept and moves no node
  by more than 3.3e-6 K.

- [ ] **Step 5: Unit tier and commit.** Expected 319 passed, 1 skipped, the five known.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py
git commit -m "Moving mesh: corners and zone move, the enthalpy is remapped by its gradient and the books close exactly"
```

### Task 5: The rebuild

**Files:**
- Modify: `reentry_model/movemesh.py` (add the method `rebuild` to `MovingMesh`; append the rebuild section)
- Test: `tests/test_reentry_model_movemesh.py` (append), `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Consumes: Tasks 1–4; `mesh.VolumeMesh`; the solver's `setup(mesh, material, emissivity)`.
- Produces: `MovingMesh.rebuild(reason="requested") -> RebuildResult(residual, relative_residual, n_nodes, n_elements,
  seconds, reason)`; module functions `default_sizes(mesh) -> dict`, `check_surface(points, faces) -> volume` (raises
  `RemeshError`), `rebuild_volume(corner_points, faces_corner, h_surface, h_core, band, ramp) -> (points, tets,
  node_of_corner)`, `locate(points, tets, query, candidates=16) -> (elements, barycentric (nq, 4))`. After a rebuild
  `moving.mesh` and `moving.solver.mesh` are the new mesh; `moving.corners` maps each corner to its new node.

- [ ] **Step 1: Append the Task 5 tests**

```python
# -- the rebuild (Task 5) ------------------------------------------------------------------------------------------
def test_a_rebuild_keeps_the_surface_exactly_and_the_energy_exactly(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat)
    solver.set_temperature(300.0 + 4000.0 * (np.linalg.norm(m.points, axis=1) - 0.03).clip(0.0))
    moving = movemesh.MovingMesh(m, solver, mat, depth=15e-3, sizes={"h_surface": 4e-3, "h_core": 20e-3, "band": 0.0, "ramp": 0.05})
    moving.move(inward(moving, 1e-3))
    corner_points = moving.mesh.points[moving.corners].copy()
    faces_before = moving.faces_corner.copy()
    V, E = moving.mesh.volume(), solver.energy()
    result = moving.rebuild()
    assert np.array_equal(moving.mesh.points[moving.corners], corner_points)
    assert np.array_equal(moving.faces_corner, faces_before)
    new_faces = np.sort(moving.mesh.surface().faces, axis=1)
    expected = np.sort(moving.corners[moving.faces_corner], axis=1)
    assert {tuple(f) for f in new_faces} == {tuple(f) for f in expected}
    assert moving.mesh.volume() == pytest.approx(V, rel=1e-12)
    assert solver.energy() == pytest.approx(E, rel=1e-12)
    assert np.all(movemesh.signed_volumes(moving.mesh.points, moving.mesh.tets) != 0.0)
    assert result.relative_residual < 1e-3 and moving.rebuilds == 1 and solver.mesh is moving.mesh


def test_a_rebuild_keeps_a_uniform_field_uniform(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat, 650.0)
    moving = movemesh.MovingMesh(m, solver, mat, sizes={"h_surface": 4e-3, "h_core": 20e-3, "band": 0.0, "ramp": 0.05})
    moving.rebuild()
    assert np.abs(solver.temperature() - 650.0).max() < 1e-9


def test_corners_keep_their_order_through_repeated_rebuilds(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    moving = movemesh.MovingMesh(m, solver_on(m, mat), mat, sizes={"h_surface": 4e-3, "h_core": 20e-3, "band": 0.0, "ramp": 0.05})
    first = moving.mesh.points[moving.corners].copy()
    for _ in range(2):
        moving.move(inward(moving, 3e-4))
        expected = moving.mesh.points[moving.corners].copy()
        moving.rebuild()
        assert np.array_equal(moving.mesh.points[moving.corners], expected)
    assert np.linalg.norm(moving.mesh.points[moving.corners], axis=1).max() < np.linalg.norm(first, axis=1).max()
    assert moving.rebuilds == 2


@pytest.mark.parametrize("case", ["two bodies", "non-manifold", "inside out"])
def test_the_validity_gate_stops_on_a_surface_it_cannot_rebuild(case):
    tet = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])            # outward
    if case == "two bodies":
        points, f = np.vstack([tet, tet + 5.0]), np.vstack([faces, faces + 4])
    elif case == "non-manifold":
        points, f = np.vstack([tet, tet[[0, 1]], [[0, -1, 0], [0, 0, -1]]]), np.vstack([faces, [[0, 1, 4], [0, 4, 5], [1, 5, 4], [0, 5, 1]]])
    else:
        points, f = tet, faces[:, [0, 2, 1]]
    with pytest.raises(movemesh.RemeshError):
        movemesh.check_surface(points, f)


```

  and to the end of `tests/test_reentry_model_fenicsx.py`:

```python
def test_a_rebuild_keeps_the_energy_on_fenicsx(coarse_sphere_mesh):
    """Moving mesh (sub-plan 19): a rebuild re-runs setup on the new mesh; the books stay exact and a step follows."""
    from reentry_model import mesh as mesh_mod
    from reentry_model import movemesh
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    s1 = solver("fenicsx", m, mat)
    s1.set_temperature(400.0 + 2000.0 * m.points[:, 0])
    moving = movemesh.MovingMesh(m, s1, mat, sizes={"h_surface": 4e-3, "h_core": 20e-3, "band": 0.0, "ramp": 0.05})
    E = s1.energy()
    moving.rebuild()
    assert s1.energy() == pytest.approx(E, rel=1e-12) and s1.mesh is moving.mesh
    s1.step(0.5, np.full(len(s1.faces), 1e5), 300.0)


```

- [ ] **Step 2: Run them to see them fail** — expected `AttributeError: 'MovingMesh' object has no attribute 'rebuild'`
  and `module 'reentry_model.movemesh' has no attribute 'check_surface'`.

- [ ] **Step 3: Add the method `rebuild` at the end of the class `MovingMesh`** (after `move`, indented as a method):

```python
    def rebuild(self, reason="requested"):
        """Regenerate the volume mesh inside the current surface triangulation and carry the enthalpy over."""
        t0 = time.perf_counter()
        mesh, solver, mat = self.mesh, self.solver, self.material
        faces = self.corners[self.faces_corner]
        check_surface(mesh.points, faces)
        points, tets, node_of_corner = rebuild_volume(mesh.points[self.corners], self.faces_corner, **self.sizes)
        h_old = mat.enthalpy(solver.temperature())
        E_old = solver.energy()
        elem, bary = locate(mesh.points, mesh.tets, points)
        h_new = np.einsum("ij,ij->i", bary, h_old[mesh.tets[elem]])
        h_new[node_of_corner] = h_old[self.corners]          # the surface is the same: corners keep their enthalpy
        params = dict(mesh.params, rebuilds=self.rebuilds + 1)
        new_mesh = mesh_mod.VolumeMesh(points, tets, params)
        solver.setup(new_mesh, mat, solver.emissivity)
        solver.set_temperature(mat.temperature_from_enthalpy(h_new))
        residual = E_old - solver.energy()
        h_new += residual / lumped_masses(points, tets, mat.rho).sum()
        solver.set_temperature(mat.temperature_from_enthalpy(h_new))
        self.mesh, self.corners = new_mesh, node_of_corner
        self._set_surface(node_of_corner[self.faces_corner])
        self.compression[:] = 0.0
        self.steps_since_rebuild = 0
        self.rebuilds += 1
        self.last_rebuild = RebuildResult(float(residual), float(residual / E_old) if E_old else 0.0, new_mesh.n_nodes,
                                          new_mesh.n_elements, time.perf_counter() - t0, reason)
        return self.last_rebuild

```

- [ ] **Step 4: Append the rebuild section at the end of `reentry_model/movemesh.py`**

```python
# -- rebuilding the mesh ---------------------------------------------------------------------------------------------
@dataclass
class RebuildResult:
    residual: float                     # J, total enthalpy the interpolation lost or gained before the offset
    relative_residual: float            # that over the total enthalpy
    n_nodes: int
    n_elements: int
    seconds: float
    reason: str


def default_sizes(mesh):
    """gmsh sizes for a rebuild, from the generator's parameters: a sphere's dense band, or a box's uniform h."""
    p = mesh.params
    if "h_surface_m" in p:
        return {"h_surface": p["h_surface_m"], "h_core": p["h_core_m"], "band": p.get("band_m", mesh_mod.DEFAULT_BAND),
                "ramp": max(p.get("band_ramp_m", mesh_mod.DEFAULT_BAND_RAMP), 1e-9)}
    h = p.get("h", mesh_mod.DEFAULT_H_SURFACE)
    return {"h_surface": h, "h_core": h, "band": 0.0, "ramp": 1.0}


def check_surface(points, faces):
    """The validity gate: every edge used by exactly two triangles, one connected surface, positive enclosed volume."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    faces = np.asarray(faces, dtype=np.int64)
    edges = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    if np.any(counts != 2):
        raise RemeshError("the surface is not manifold: {} edges are not shared by exactly two triangles".format(
            int(np.sum(counts != 2))))
    nodes = np.unique(faces)
    local = np.full(points.shape[0], -1, dtype=np.int64)
    local[nodes] = np.arange(len(nodes))
    e = local[edges]
    graph = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(len(nodes), len(nodes)))
    n_parts = connected_components(graph, directed=False)[0]
    if n_parts != 1:
        raise RemeshError("the body has come apart into {} pieces; fragmentation is out of scope".format(n_parts))
    a, b, c = (points[faces[:, i]] for i in range(3))
    volume = float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)
    if volume <= 0.0:
        raise RemeshError("the surface encloses a non-positive volume ({:.3g} m3)".format(volume))
    return volume


def rebuild_volume(corner_points, faces_corner, h_surface, h_core, band, ramp):
    """Tetrahedra inside the closed triangulation (corner_points, faces_corner), keeping it exactly. Returns the new
    points, tets and the new node id of every corner (matched by coordinates: gmsh renumbers the nodes)."""
    import gmsh
    from scipy.spatial import cKDTree
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("rebuild")
        surface = gmsh.model.addDiscreteEntity(2)
        gmsh.model.mesh.addNodes(2, surface, list(range(1, len(corner_points) + 1)), corner_points.ravel().tolist())
        gmsh.model.mesh.addElementsByType(surface, 2, [], (np.asarray(faces_corner) + 1).ravel().tolist())
        loop = gmsh.model.geo.addSurfaceLoop([surface])
        gmsh.model.geo.addVolume([loop])
        gmsh.model.geo.synchronize()
        dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", [surface])
        thr = gmsh.model.mesh.field.add("Threshold")
        d_min = band
        d_max = band + ramp if band > 0.0 else max(ramp, 1e-9)
        for key, value in (("InField", dist), ("SizeMin", h_surface), ("SizeMax", h_core), ("DistMin", d_min),
                           ("DistMax", d_max)):
            gmsh.model.mesh.field.setNumber(thr, key, value)
        gmsh.model.mesh.field.setAsBackgroundMesh(thr)
        for option in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
            gmsh.option.setNumber(option, 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)          # HXT, as mesh.generate_sphere_mesh
        # one thread: fenicsx_env's gmsh 4.15.2 segfaults meshing a discrete surface's volume with HXT on four threads
        # (checked 2026-10-11; drama_env's does not), and a rebuild takes under two seconds on one
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.option.setNumber("Mesh.MaxNumThreads3D", 1)
        gmsh.model.mesh.generate(3)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        types, _, nodes = gmsh.model.mesh.getElements(3)
    finally:
        gmsh.finalize()
    if 4 not in list(types):
        raise RemeshError("gmsh produced no tetrahedra inside the surface")
    coords = np.asarray(coords).reshape(-1, 3)
    tags = np.asarray(tags, dtype=np.int64)
    index = np.full(tags.max() + 1, -1, dtype=np.int64)
    index[tags] = np.arange(len(tags))
    tets = index[np.asarray(nodes[list(types).index(4)], dtype=np.int64).reshape(-1, 4)]
    used = np.unique(tets)                                   # drop nodes no tetrahedron uses (none expected)
    renumber = np.full(len(coords), -1, dtype=np.int64)
    renumber[used] = np.arange(len(used))
    coords, tets = coords[used], renumber[tets]
    d, node_of_corner = cKDTree(coords).query(corner_points)
    if d.max() > 0.0:
        raise RemeshError("gmsh moved a surface node by {:.3g} m".format(float(d.max())))
    return coords, tets, node_of_corner.astype(np.int64)


def locate(points, tets, query, candidates=16):
    """For each query point, an element of (points, tets) holding it and its barycentric coordinates there; a point
    just outside every candidate gets the candidate it is least outside of (linear extrapolation)."""
    from scipy.spatial import cKDTree
    centroids = points[tets].mean(axis=1)
    k = min(candidates, len(tets))
    _, cand = cKDTree(centroids).query(query, k=k)
    cand = cand.reshape(len(query), k)
    x = points[tets[cand]]                                   # (nq, k, 4, 3)
    J = np.stack([x[:, :, 1] - x[:, :, 0], x[:, :, 2] - x[:, :, 0], x[:, :, 3] - x[:, :, 0]], axis=3)
    lam = np.linalg.solve(J, (query[:, None, :] - x[:, :, 0])[..., None])[..., 0]
    bary = np.concatenate([1.0 - lam.sum(axis=2, keepdims=True), lam], axis=2)    # (nq, k, 4)
    best = np.argmax(bary.min(axis=2), axis=1)
    rows = np.arange(len(query))
    return cand[rows, best], bary[rows, best]

```

- [ ] **Step 5: Run them to see them pass** — expected 22 passed in the file and the FEniCSx rebuild test passing.
  A rebuild of the coarse test sphere takes well under a second.

- [ ] **Step 6: Both tiers and commit.** Expected unit tier 325 passed, 1 skipped, the five known; FEniCSx 16 passed.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py tests/test_reentry_model_fenicsx.py
git commit -m "Moving mesh: rebuild inside the unchanged surface triangulation, enthalpy carried over with exact totals"
```

### Task 6: The rebuild policy

**Files:**
- Modify: `reentry_model/movemesh.py` (add `needs_rebuild` and `move_with_policy` to `MovingMesh`, after `rebuild`)
- Test: `tests/test_reentry_model_movemesh.py` (append)

**Interfaces:**
- Produces: `MovingMesh.needs_rebuild() -> "compression" | None` (the zone pressed to half its depth, after at least
  `min_gap` moves); `MovingMesh.move_with_policy(u_corners) -> MoveResult` (rebuild first if compression calls for it;
  if the move is refused, rebuild with reason "quality" and move again; refused again raises `RuntimeError` containing
  "even after a rebuild"). Stage 2 calls `move_with_policy` once per macro step.

- [ ] **Step 1: Append the Task 6 tests**

```python
# -- the rebuild policy (Task 6) -----------------------------------------------------------------------------------
def test_compression_triggers_a_rebuild_only_after_the_minimum_gap(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    moving = movemesh.MovingMesh(m, solver_on(m, mat), mat, depth=4e-3, min_gap=3,
                                 sizes={"h_surface": 4e-3, "h_core": 20e-3, "band": 0.0, "ramp": 0.05})
    rebuilt_at = []
    for step in range(8):
        r = moving.move_with_policy(inward(moving, 0.8e-3))
        if r.rebuilt:
            rebuilt_at.append(step)
    # compression passes half the 4 mm zone (2 mm) after the third move; the gap of 3 is met then too
    assert rebuilt_at and rebuilt_at[0] == 3 and moving.rebuilds >= 1


def test_a_refused_move_rebuilds_and_then_moves(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    moving = movemesh.MovingMesh(m, solver_on(m, mat), mat, depth=8e-3, ratio_min=0.1, min_gap=1000,
                                 sizes={"h_surface": 4e-3, "h_core": 20e-3, "band": 0.0, "ramp": 0.05})
    for _ in range(40):                                  # press the zone until a move is refused (no compression rebuild)
        u = inward(moving, 0.4e-3)
        if moving.move(u).refused:
            break
    else:
        pytest.fail("the zone never degraded")
    r = moving.move_with_policy(u)
    assert r.rebuilt and not r.refused and moving.last_rebuild.reason == "quality"


def test_a_move_refused_on_a_fresh_mesh_is_a_model_failure(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    moving = movemesh.MovingMesh(m, solver_on(m, mat), mat, depth=4e-3,
                                 sizes={"h_surface": 4e-3, "h_core": 20e-3, "band": 0.0, "ramp": 0.05})
    with pytest.raises(RuntimeError, match="even after a rebuild"):
        moving.move_with_policy(inward(moving, 3.5e-3))   # most of the 4 mm zone in one step


```

- [ ] **Step 2: Run them to see them fail** — expected `AttributeError: 'MovingMesh' object has no attribute
  'move_with_policy'`.

- [ ] **Step 3: Add the two methods to `MovingMesh`, after `rebuild` and before the `# -- rebuilding the mesh`
  section**

```python
    def needs_rebuild(self):
        """The reason a rebuild is due by compression (the zone pressed to half its depth), or None."""
        if self.steps_since_rebuild >= self.min_gap and self.compression.max() > 0.5 * self.depth:
            return "compression"
        return None

    def move_with_policy(self, u_corners):
        """Rebuild first if compression calls for it; move; if the move is refused, rebuild and move again. A move
        refused on a fresh mesh is a model failure (RuntimeError, exit code 1)."""
        rebuilt = False
        reason = self.needs_rebuild()
        if reason:
            self.rebuild(reason)
            rebuilt = True
        result = self.move(u_corners)
        if result.refused:
            self.rebuild("quality")
            rebuilt = True
            result = self.move(u_corners)
            if result.refused:
                raise RuntimeError("a corner move of up to {:.3g} m degrades the mesh even after a rebuild".format(
                    float(np.linalg.norm(u_corners, axis=1).max())))
        result.rebuilt = rebuilt
        return result



```

- [ ] **Step 4: Run them to see them pass** — expected 25 passed in the file.

- [ ] **Step 5: Unit tier and commit.** Expected 328 passed, 1 skipped, the five known.

```bash
git add reentry_model/movemesh.py tests/test_reentry_model_movemesh.py
git commit -m "Moving mesh: rebuild when the zone is half compressed or a move would degrade an element"
```

### Task 7: A moving front against an exact solution

**Files:**
- Test: `tests/test_reentry_model_movemesh.py` (append), `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Consumes: Tasks 1–6 and each backend's `step(dt, q_conv, T_amb, dirichlet=...)`.
- Produces: `run_moving_front(name)` in the test module, which the FEniCSx test imports.

This is level 1's check of the moving mesh and the conduction together. The face x = 0 of a 20 mm bar of 0.1 mm cells
is held at 600 K and moves into the bar at 10 mm/s for 0.5 s; the exact temperature ahead of a front at a fixed
temperature moving at constant speed into a solid at T0 is (Carslaw and Jaeger) θ = ½ [erfc((ξ − v t)/(2√(α t)))
e^(−v ξ/α) + erfc((ξ + v t)/(2√(α t)))], with θ = (T − T0)/(Ts − T0), ξ the distance ahead of the front and α = 1e-5 m²/s.
It was checked while planning against a 4000-cell finite-difference solution in the front's frame: largest difference
3.6e-6. The bar's sides slide with the zone, so they stay on their planes.

- [ ] **Step 1: Append the Task 7 tests**

```python
# -- a moving front against an exact solution (Task 7) -------------------------------------------------------------
def exact_front(xi, t, alpha, v):
    """theta = (T - T0) / (Ts - T0) ahead of a front at fixed Ts moving at v into a solid at T0 (Carslaw & Jaeger)."""
    s = 2.0 * np.sqrt(alpha * t)
    return 0.5 * (erfc((xi - v * t) / s) * np.exp(-v * xi / alpha) + erfc((xi + v * t) / s))


def run_moving_front(name):
    """A 20 mm bar of 0.1 mm cells whose x = 0 face is held at 600 K and moves into the bar at 10 mm/s for 0.5 s."""
    mat = constant_material()
    m = mesh_mod.box_mesh(20e-3, 0.2e-3, 0.2e-3, 0.1e-3)
    solver = solver_on(m, mat, 300.0, name)
    D, v, dt = 15e-3, 1e-2, 5e-3
    moving = movemesh.MovingMesh(m, solver, mat, depth=D)
    x_front = 0.0
    for _ in range(100):
        cx = moving.mesh.points[moving.corners, 0]
        u = np.zeros((len(moving.corners), 3))
        u[:, 0] = v * dt * np.clip(1.0 - (cx - x_front) / D, 0.0, 1.0)              # sides slide with the zone
        r = moving.move(u)
        assert not r.refused
        x_front += v * dt
        front = np.flatnonzero(np.abs(moving.mesh.points[:, 0] - x_front) < 1e-9)
        solver.step(dt, np.zeros(len(solver.faces)), 300.0, dirichlet=(front, np.full(len(front), 600.0)))
    return moving, solver, x_front


def test_a_moving_front_follows_the_exact_solution():
    moving, solver, x_front = run_moving_front("skfem")
    x = moving.mesh.points[:, 0]
    theta = (solver.temperature() - 300.0) / 300.0
    exact = exact_front(x - x_front, 0.5, 1e-5, 1e-2)
    assert np.abs(theta - exact).max() < 0.005          # 0.0010 measured while planning: 0.1 % of the rise
```

  and to the end of `tests/test_reentry_model_fenicsx.py`:

```python
def test_a_moving_front_matches_the_skfem_backend():
    """Moving mesh (sub-plan 19): the moving-front case of test_reentry_model_movemesh on both backends."""
    import sys, os
    sys.path.insert(0, os.path.dirname(__file__))
    from test_reentry_model_movemesh import run_moving_front
    _, s_fx, _ = run_moving_front("fenicsx")
    _, s_sk, _ = run_moving_front("skfem")
    assert np.abs(s_fx.temperature() - s_sk.temperature()).max() < 1e-6
```

- [ ] **Step 2: Run them** — `"$PY" -m pytest tests/test_reentry_model_movemesh.py -q -k front` (expected 1 passed: the
  largest error was 0.0010 of the temperature rise while planning, against a tolerance of 0.005) and the FEniCSx tier
  (expected 17 passed; the two backends agree to 1e-6 K).

- [ ] **Step 3: Both tiers and commit.** Expected unit tier 329 passed, 1 skipped, the five known; FEniCSx 17 passed.

```bash
git add tests/test_reentry_model_movemesh.py tests/test_reentry_model_fenicsx.py
git commit -m "Moving mesh: a heated front moving into a bar follows the exact solution on both backends"
```

### Task 8: Measurement on the 100 mm sphere, and the record

**Files:**
- Create: `analysis/movemesh_recession.py`
- Modify (outside the work copy, documents only): `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`
  (facts), this sub-plan (a results section)

- [ ] **Step 1: The run path is unchanged.** Nothing in the model's run path calls the new code, so a flight from the
  work copy must be bit-identical to one from the base. Run the 100 mm non-melting flight to 30 s from each, with the
  same seed, and compare:

```bash
P="/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
for copy in work-2026-10-08-skin work-2026-10-12-moving-mesh; do
  (cd "$P/$copy/code" && caffeinate -i "$PY" -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 \
     --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --heating physics --t-max 30 --seed 12345 \
     --outdir "$P/work-2026-10-12-moving-mesh/runs/bitid_$copy")
done
"$PY" - "$P/work-2026-10-12-moving-mesh/runs/bitid_work-2026-10-08-skin" "$P/work-2026-10-12-moving-mesh/runs/bitid_work-2026-10-12-moving-mesh" <<'PYEOF'
import glob, json, sys
import pandas as pd
a, b = (pd.read_csv(glob.glob(d + "/*.csv")[0]) for d in sys.argv[1:3])
assert list(a.columns) == list(b.columns) and a.equals(b), "history differs"
ra, rb = (json.load(open(glob.glob(d + "/*.json")[0]))["results"] for d in sys.argv[1:3])
differ = [k for k in ra if k != "runtime_s" and ra[k] != rb.get(k)]
assert not differ, differ
print("bit-identical:", len(a), "rows,", len(a.columns), "columns,", len(ra) - 1, "result fields")
PYEOF
```

  Expected: `bit-identical: 61 rows, 31 columns, ...` (checked while planning on the plan's tested copy). Only the
  `results` block is compared, without `runtime_s`; the rest of the JSON holds file paths and the git commit, which
  differ by construction.

- [ ] **Step 2: Create the measurement script**

```python
"""Sub-plan 19 measurement: the moving 3D mesh on the default 100 mm sphere through a prescribed recession of the
front, with heating, to measure rebuild count, rebuild time, energy residuals, element quality and the cost of a
macro step on moving geometry.

Every corner recedes along its inward normal by NOSE_STEP x max(cos theta, 0) per 0.5 s macro step (theta its angle
from the flight direction, -x), so the nose goes 0.5 mm per step, a speed of 1 mm/s, as on the 100 mm flight at
50-60 s (fact: 16.6 mm at 50 s, 29.0 mm at 60 s). The front is heated at Q_STAG x max(cos theta, 0) (no melting
material, `AA7075_nomelt`) so the field the remap and the rebuilds carry has the flight's steep near-surface gradient.
The run stops at NOSE_TOTAL of nose recession.

    python analysis/movemesh_recession.py [--nose-total-mm 30] [--solver skfem|fenicsx] [--out PATH.json]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reentry_model import material, mesh, movemesh  # noqa: E402
from reentry_model.thermal import thermal_solver  # noqa: E402

DT, NOSE_STEP, Q_STAG = 0.5, 0.5e-3, 1.0e6

parser = argparse.ArgumentParser()
parser.add_argument("--nose-total-mm", type=float, default=30.0)
parser.add_argument("--solver", default="skfem")
parser.add_argument("--out", default="movemesh_recession.json")
parser.add_argument("--no-smoothing", action="store_true", help="corners move along their normals only")
parser.add_argument("--smooth-lambda", type=float, default=None)
parser.add_argument("--smooth-clamp", type=float, default=None)
args = parser.parse_args()

mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_nomelt"])
m = mesh.sphere_mesh(0.05)
solver = thermal_solver(args.solver)
solver.setup(m, mat, mat.emissivity)
solver.set_temperature(300.0)
moving = movemesh.MovingMesh(m, solver, mat)
flight = np.array([-1.0, 0.0, 0.0])
rows, t_total, failure = [], time.perf_counter(), None
steps = int(round(args.nose_total_mm * 1e-3 / NOSE_STEP))
for k in range(steps):
    t0 = time.perf_counter()
    normal, _ = moving.frame()
    cos = np.clip(normal @ flight, 0.0, None)
    u = -NOSE_STEP * cos[:, None] * normal
    if not args.no_smoothing:
        u = u + moving.smoothing(normal, args.smooth_lambda, args.smooth_clamp)
    try:
        result = moving.move_with_policy(u)
    except RuntimeError as exc:                       # a move no rebuild can absorb: record where, and stop
        failure = {"failed_at_mm": (k + 1) * NOSE_STEP * 1e3, "error": str(exc)}
        print(json.dumps(failure), flush=True)
        break
    t_move = time.perf_counter() - t0
    surface = moving.mesh.surface()
    q = Q_STAG * np.clip(surface.normals @ flight, 0.0, None)
    t1 = time.perf_counter()
    step = solver.step(DT, q, 300.0)
    t_step = time.perf_counter() - t1
    ratio = movemesh.radius_ratio(moving.mesh.points, moving.mesh.tets)
    row = {"step": k + 1, "nose_mm": (k + 1) * NOSE_STEP * 1e3, "rebuilt": result.rebuilt,
           "rebuilds": moving.rebuilds, "correction_J": result.correction,
           "swept_J": float(result.swept_enthalpy.sum()), "worst_ratio": float(ratio.min()),
           "worst_ratio_moved": result.worst_ratio,
           "surface_quality": moving.surface_quality(), "surface_fold_deg": moving.surface_fold(), "move_s": t_move, "step_s": t_step,
           "T_max": float(step.T.max()), "n_nodes": moving.mesh.n_nodes, "n_elements": moving.mesh.n_elements}
    if result.rebuilt:
        rb = moving.last_rebuild
        row.update({"rebuild_s": rb.seconds, "rebuild_residual_J": rb.residual,
                    "rebuild_relative_residual": rb.relative_residual, "rebuild_reason": rb.reason})
    rows.append(row)
    if (k + 1) % 10 == 0 or result.rebuilt:
        print(json.dumps({key: (round(v, 6) if isinstance(v, float) else v) for key, v in row.items()}), flush=True)
summary = {"failure": failure, "completed_mm": rows[-1]["nose_mm"] if rows else 0.0, "solver": args.solver, "steps": steps, "nose_total_mm": args.nose_total_mm, "rebuilds": moving.rebuilds,
           "wall_s": time.perf_counter() - t_total,
           "rebuild_s": [r["rebuild_s"] for r in rows if "rebuild_s" in r],
           "rebuild_relative_residuals": [r["rebuild_relative_residual"] for r in rows if "rebuild_s" in r],
           "max_abs_correction_over_swept": max(abs(r["correction_J"]) / r["swept_J"] for r in rows if r["swept_J"] > 0),
           "worst_ratio_min": min(r["worst_ratio"] for r in rows),
           "surface_quality_min": min(r["surface_quality"] for r in rows), "smoothing": not args.no_smoothing,
           "move_s_median": float(np.median([r["move_s"] for r in rows])),
           "step_s_median": float(np.median([r["step_s"] for r in rows]))}
print(json.dumps(summary, indent=1))
json.dump({"summary": summary, "rows": rows}, open(args.out, "w"), indent=1)
```

- [ ] **Step 3: Run it with and without the tangential smoothing**, one after the other (not side by side: the run time
  is measured), under `caffeinate -i` on mains power:

```bash
"$PY" analysis/movemesh_recession.py --out ../measure/recession_smooth.json > ../measure/recession_smooth.log 2>&1
"$PY" analysis/movemesh_recession.py --no-smoothing --out ../measure/recession_nosmooth.json > ../measure/recession_nosmooth.log 2>&1
```

  While planning (one run each, 2026-10-11): with smoothing the nose receded 29.5 mm before a move failed even after a
  rebuild, 20 rebuilds, the first at 10.5 mm by compression and every later one by quality — every macro step from
  21 mm; without smoothing 25.0 mm and 15 rebuilds. A rebuild took 1.6–2.3 s on one gmsh thread (130 000–177 000
  tetrahedra), against the spec's estimate of about 42 s; its pre-correction residual was 0.07–0.27 % of the total
  enthalpy, and a move's offset at most 1.4 % of what the corners swept. The worst surface triangle stayed at 0.71 with
  smoothing (0.36 without). What ends the run is the shape: the largest angle between neighbouring surface triangles grew
  from 19.5 degrees at 10 mm and 51.9 at 20 mm to 95 degrees, as recession in proportion to how directly a corner faces
  the flow sharpened the shoulder; an element whose four nodes lie on a fold sharper than one triangle is flattened by
  its own corners' motion, and a rebuild inside the same surface cannot undo it. The moment of failure depends on the
  meshes gmsh happens to make: with four gmsh threads the same runs failed at 27.0 and 20.5 mm. A macro step's mesh work
  took a median 0.45 s with smoothing and 0.27 s without, against 0.19–0.20 s for the conduction step.

- [ ] **Step 4: Record the facts** in `00-shared-context.md` (the next numbers after 101): the rebuild cost and residual,
  the policy's rebuild record, the smoothing's effect, the fold at failure and the dependence on gmsh's meshes, the mesh
  work against the conduction step, the AMG fix, and the one-thread gmsh setting with its reason. Mark each measured
  once.

- [ ] **Step 5: Add a results section to this sub-plan** with the counts of every tier after each task, the
  measurements, and anything that differed from the numbers above. The fold limit is recorded, not resolved, here:
  Asha decided on 2026-10-11 to measure it under heating-driven recession in stage 2 (sub-plan 20) before building any
  surface refinement.
