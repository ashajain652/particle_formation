# Sub-plan: Task 1 — Surface mesh (dense band), active set, box mesh

> Body extracted verbatim from `2026-09-20-melt-spraying.md` (lines 300–768) **and amended on 2026-09-27 by the
> section immediately below, which takes precedence wherever the two disagree.** Read `00-shared-context.md` first
> (Global Constraints, the measured facts — now 45 — and the file structure), including the 2026-09-27 amendment
> block with facts 38–45, before refining this into an implementation plan.

**Depends on:** Step 2's existing mesh module only. Nothing else in this Step 3 plan needs to land first.
**Produces, for later tasks:** the dense-band surface mesh, the active/dead element set (with an O(1) face table), the box mesh, the `derived_points` array and `surface(derived=...)`, and the surface helper geometry (tangents, projected areas, edges) that Tasks 3, 5, 6, 9, 16 and 17 all read from.
**Character:** numerics/mesh geometry.
**Read before implementing:** the amendment below **first**, then measured facts 38–43 in the shared context. Fact 6
(gmsh's boundary-layer extrusion does not work on the OCC sphere, so prism layers are built by radial projection)
still holds and still explains the `add_prism_layers` code below — but prism layers are no longer the default and no
longer the mechanism that keeps the surface resolved. Do not spend time improving them.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan — file list, function signatures, test plan, and acceptance criteria — precise enough that another agent could implement it without reading the master document.

---

## Amendment of 2026-09-27 — the dense band replaces the prism stack

Source: `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md` §6 and §8, approved 2026-09-27.
Measured facts 38–43. **This section overrides the extracted body wherever they conflict.** Two new tasks, 16
(`recession.py`) and 17 (`remesh.py`), consume what this one produces; build this first.

### Why the extracted design is superseded

The prism stack keeps the surface resolved for 3.75 mm. Fact 14 has the nose receding **77 mm**, so the skin covers
**4.9 %** of the recession; for the rest, fact 39 measures the recession quantum at **about 6 mm** because the graded
core *is* the surface. Refining the mesh instead is unaffordable: 1.0 mm cells over a 15 mm band cost 3.07 M
elements and a five-hour flight, against 310 k and fifty-one minutes at 2.16 mm. And the stack cannot survive
remeshing at all — a prism layer is an **offset**, which self-intersects once it exceeds the local concave radius of
curvature, and the windward face erodes into a dish.

The stack is also not load-bearing. Fact 12's sensitivity varied the layer count 2/4/6 — total fine radial depth
0.75 mm to 15.75 mm — and moved the sprayed mass by **under 0.2 %** and the demise altitude by **under 0.3 %**, while
costing 15 430 × 3 × 4 = **185 160 tetrahedra**, 73 % of the present mesh. Fact 45's near-isothermal body explains
why: there is no thin thermal boundary layer inside the solid to resolve, so the resolution requirement is
geometric, not thermal.

### What changes

**1. A distance-field band, not prism layers.** The existing `Distance` + `Threshold` field in
`generate_sphere_mesh` gains a non-zero `DistMin`:

```
SizeMin = h_surface        (2.16 mm)     SizeMax = h_core   (8 mm)
DistMin = band_thickness   (15 mm)       DistMax = band_thickness + ramp
```

which reads: keep cells at `h_surface` everywhere within 15 mm of the surface, then grow them. A distance query is
well defined at every interior point of **any** shape, concave included, because nothing is being offset — which is
precisely why it survives remeshing where prisms do not.

New constants: `DEFAULT_BAND = 15e-3`, `DEFAULT_BAND_RAMP = 8e-3`. **`DEFAULT_LAYERS` becomes 0.** `sphere_mesh`
gains `band=DEFAULT_BAND`; `generate_sphere_mesh` gains `band` and `ramp` and must include them in
`mesh_file_name`, or a banded mesh will collide in the cache with an unbanded one of the same `h`.

**2. Prism layers stay, demoted.** `layers`, `layer_thicknesses`, `split_prisms` and `add_prism_layers` are kept
verbatim and tested exactly as the body below specifies — fact 12's pinned verification devices already run with
`--prism-layers 0`, and Task 14 needs both settings to re-measure. They are simply no longer the default, and no
remeshed mesh ever has them.

**3. `derived_points` and `surface(derived=...)`.** `VolumeMesh` gains `derived_points`, a `(n, 3)` array
initialised to a copy of `points`, which Task 16 writes and **nothing else ever modifies**. `surface(derived=False)`
takes a keyword: with `derived=True` the returned `SurfaceMesh` has its centroids, normals and areas computed from
`derived_points` instead. `points` itself is never touched, and that is what keeps the thermal solver — which reads
element geometry volumetrically — completely unaffected (spec §5, fact 45).

`SurfaceMesh` gains `smoothed_normals(passes=8, lam=0.53, mu=-0.55, clamp=0.35)`, a Taubin-smoothed normal field,
because the raw staircase normals deviate from smooth by up to 45° (a step is one cell tall and one cell wide) and θ
is what the Lees distribution reads. This is one sparse matrix-vector product over 15 430 patches. **Applied to
today's mesh it is already an improvement, independently of everything else here.**

**4. Coarseness diagnostics.** `SurfaceMesh.edge_lengths()` and, on `VolumeMesh`, `element_depths()` — each
element's distance to the boundary, computed once per mesh with a `cKDTree` over the boundary triangle centroids.
Task 17's remesh triggers read both; Task 10 writes `surface_edge_max_mm` and `surface_edge_mean_mm` to the history
so every run records how coarse its surface actually was.

### Interfaces, superseding the `Interfaces:` line above

- `sphere_mesh(radius, h_surface, h_core, mesh_dir, layers=0, layer_thickness=..., growth=..., band=DEFAULT_BAND, ramp=DEFAULT_BAND_RAMP)`
- `generate_sphere_mesh(radius, h_surface, h_core, mesh_dir, band=DEFAULT_BAND, ramp=DEFAULT_BAND_RAMP)`
- `mesh_file_name(radius, h_surface, h_core, band, ramp)`
- `VolumeMesh.derived_points`, `.surface(derived=False)`, `.element_depths()`
- `SurfaceMesh.smoothed_normals(...)`, `.edge_lengths()`
- `DEFAULT_LAYERS = 0`, `DEFAULT_BAND = 15e-3`, `DEFAULT_BAND_RAMP = 8e-3`

Everything else in the extracted `Interfaces:` line is unchanged: the active set, the face table, `deactivate`,
`element_faces`, `active_nodes`, `volume`, `surface`, `box_mesh`, `tangent_from`, `projected_area`, `edges`.

### Tests to add, alongside the five in the body below

```python
def test_band_keeps_cells_fine_through_its_depth(tmp_path_factory):
    """The whole 15 mm band is at h_surface; only beyond it does the size grow (fact 39 measures what the
    ungraded mesh does instead: 2.74 mm at the surface rising to 8.32 mm at the centre)."""
    m = mesh.sphere_mesh(R, 2.16e-3, 8e-3, str(tmp_path_factory.mktemp("m")), layers=0, band=15e-3)
    d, a = m.element_depths(), (m.element_volumes() / 0.117851) ** (1 / 3.0)
    assert a[d < 10e-3].mean() == pytest.approx(2.16e-3, rel=0.35)         # inside the band
    assert a[d > 25e-3].mean() > 1.8 * a[d < 10e-3].mean()                 # coarsens beyond it
    assert m.params["band_m"] == 15e-3 and m.params["layers"] == 0

def test_mesh_file_name_separates_banded_from_unbanded():
    assert mesh.mesh_file_name(R, 2e-3, 8e-3, 15e-3, 8e-3) != mesh.mesh_file_name(R, 2e-3, 8e-3, 0.0, 8e-3)

def test_derived_points_never_touch_the_solver_geometry(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    v0, p0 = m.volume(), m.points.copy()
    bn = m.boundary_nodes()
    m.derived_points[bn] *= 0.98                                            # a 2 % recession, derived only
    assert np.allclose(m.points, p0) and m.volume() == pytest.approx(v0)    # element geometry untouched
    assert m.surface(derived=True).area < 0.999 * m.surface().area          # but the derived surface shrank

def test_smoothed_normals_remove_the_staircase(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    s = m.surface()
    rng = np.random.default_rng(0)
    m.derived_points[m.boundary_nodes()] += rng.normal(0.0, 2.0e-4, (len(m.boundary_nodes()), 3))   # cell-scale noise
    sd = m.surface(derived=True)
    exact = s.normals / np.linalg.norm(s.normals, axis=1)[:, None]
    raw = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", sd.normals, exact), -1, 1)))
    sm = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", sd.smoothed_normals(), exact), -1, 1)))
    assert sm.mean() < 0.5 * raw.mean()                                     # smoothing halves the normal error

def test_element_depths_are_zero_at_the_surface_and_maximal_at_the_centre(layered_mesh):
    d = layered_mesh.element_depths()
    assert d.min() < 1e-3 and d.max() == pytest.approx(R, rel=0.15)
```

### Acceptance criteria for this task

1. A banded mesh of the 100 mm sphere at `h_surface` 2.16 mm, `h_core` 8 mm, band 15 mm, `layers=0` has **310 000
   ± 15 %** tetrahedra (spec §6), and its cache file name differs from the unbanded one.
2. Cells stay within 35 % of `h_surface` throughout the band and coarsen beyond it.
3. `points` is bit-identical before and after any write to `derived_points`; `volume()` is unchanged.
4. `smoothed_normals` at least halves the mean normal error against an analytic sphere under cell-scale noise.
5. The five existing tests in the body below still pass with `layers=2`, unchanged — the prism path is demoted, not
   broken.
6. `"$PY" -m pytest tests/test_reentry_model_mesh.py -q` green, and the session-scoped fixtures in `conftest.py`
   are reused rather than re-meshing (CLAUDE.md: the thermal tests dominate the unit-tier runtime).

### Corrections to the extracted body

- Element **death** is split across three tasks and the body's `Interfaces:` line misattributes one of them. This
  task owns the **mechanism**: `VolumeMesh.active`, `deactivate(elements) -> (gone_face_ids, new_face_ids)`,
  `element_faces`, `n_active`, `active_nodes`, and the face table that makes a deactivation cost O(dead elements)
  rather than a re-sort, with `volume()` and `surface()` taken over the active set. It owns no policy at all. The
  **decision** — a patch owner dies at `φ_e ≤ PHI_DEATH`, interior elements keep `PHI_MIN = 1e-3` so no cavity can
  open, deaths cascade within the step, and the dying element's remainder joins the film — is **Task 9**, not
  Task 10 as the body's description of `element_faces` says. Task 3 owns the solver's response (`set_fractions`,
  with 0 meaning dead, and the pinning of nodes that have lost all material).
- The body's code block still reads `PHI_DEATH = 0.05` and `DEFAULT_LAYERS = 4`. Both are superseded: 0.50 by
  fact 44 (in Task 9, where the constant lives) and 0 by the amendment above.


### Not in this task

Moving the nodes (Task 16), remeshing (Task 17), and anything touching `thermal/` — spec §5 adopts Design B, in
which the solver's mesh never moves, so **no backend change belongs here**.

---

### Task 1 as originally extracted (superseded in part by the amendment above)


**Files:**
- Modify (replace): `reentry_model/mesh.py`
- Test: `tests/test_reentry_model_mesh.py` (append)

**Interfaces:**
- Consumes: gmsh generation and `VolumeMesh`/`SurfaceMesh` of Step 2 (kept: `generate_sphere_mesh`, `load_mesh`, `boundary_faces`, `mesh_file_name`, `centre_node`, `boundary_nodes`, `angles_to`, `facet_mean`, `patch_toward`).
- Produces: `sphere_mesh(radius, h_surface, h_core, mesh_dir, layers=0, layer_thickness=0.25e-3, growth=2.0)`; `layer_thicknesses(layers, layer_thickness, growth)`; `split_prisms(bottom, top)`; `add_prism_layers(inner, radius, layers, layer_thickness, growth)`; `box_mesh(lx, ly, lz, h)`; `VolumeMesh.active` (bool per element), `.element_layer` (−1 core, 0 outermost layer), `.n_active`, `.deactivate(elements) -> (gone_face_ids, new_face_ids)`, `.active_nodes()`, `.element_faces(elements) -> (n, 4) face ids` (the stable ids of an element's own four faces, which the death hand-over of Task 10 uses to find the patches a vanishing element exposes), `.volume()` (active), `.surface()` over the active set with `SurfaceMesh.owner` (element per patch) and `.face_ids` (stable ids into the face table `_face_nodes`); `SurfaceMesh.tangent_from(v_hat)`, `.projected_area(v_hat)`, `.edges() -> (edge node pairs, patch i, patch j)`; constants `DEFAULT_LAYERS = 4`, `DEFAULT_LAYER_THICKNESS = 0.25e-3`, `DEFAULT_LAYER_GROWTH = 2.0`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_mesh.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: prism layers, the active set, the box mesh and the surface helpers

@pytest.fixture(scope="module")
def layered_mesh(tmp_path_factory):
    """100 mm sphere, 4 mm / 20 mm inner mesh with two 0.25 mm / 0.5 mm prism layers (~5 k nodes)."""
    return mesh.sphere_mesh(R, 4e-3, 20e-3, str(tmp_path_factory.mktemp("meshes")), layers=2)


def test_layer_thicknesses_and_defaults():
    assert np.allclose(mesh.layer_thicknesses(4, 0.25e-3, 2.0), [0.25e-3, 0.5e-3, 1e-3, 2e-3])
    assert mesh.DEFAULT_LAYERS == 4 and mesh.DEFAULT_LAYER_THICKNESS == 0.25e-3 and mesh.DEFAULT_LAYER_GROWTH == 2.0
    with pytest.raises(ValueError):
        mesh.sphere_mesh(R, 4e-3, 20e-3, layers=-1)
    with pytest.raises(ValueError):
        mesh.sphere_mesh(0.001, 4e-3, 20e-3, layers=4)                      # layers thicker than half the radius


def test_split_prisms_is_conforming_and_fills_the_prism():
    bottom = np.array([[0, 1, 2], [1, 3, 2]])                                   # two triangles sharing the edge 1-2
    tets = mesh.split_prisms(bottom, bottom + 4)
    assert tets.shape == (6, 4)
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=float)
    m = mesh.VolumeMesh(pts, tets)
    assert m.volume() == pytest.approx(1.0) and (m.element_volumes() > 0).all()
    faces, _ = mesh.boundary_faces(tets)
    assert len(faces) == 12                                                     # the cube's 6 faces, 2 triangles each: no internal single-use face


def test_layered_sphere_geometry(layered_mesh):
    m = layered_mesh
    s = m.surface()
    assert set(np.unique(m.element_layer)) == {-1, 0, 1}
    outer = np.bincount(m.element_layer + 1)
    assert outer[1] == outer[2] == 3 * s.n_patches                            # three tets per prism, one prism per patch
    assert s.area == pytest.approx(4.0 * math.pi * R ** 2, rel=2e-3) and m.volume() == pytest.approx(4.0 / 3.0 * math.pi * R ** 3, rel=3e-3)
    radii = np.linalg.norm(m.points, axis=1)
    shells = np.unique(np.round(radii[radii > R - 1e-3] * 1e3, 3))
    assert np.allclose(shells, [49.25, 49.75, 50.0])                         # inner mesh at R - 0.75 mm, layers 0.5 and 0.25 mm
    assert m.element_layer[s.owner].max() == 0 and m.element_layer[s.owner].min() == 0
    e, i, j = s.edges()
    assert len(e) == 3 * s.n_patches // 2                                    # closed manifold surface
    assert np.abs((s.areas[:, None] * s.normals).sum(axis=0)).max() < 1e-12 * s.area
    assert s.projected_area([1.0, 0.0, 0.0]) == pytest.approx(math.pi * R ** 2, rel=2e-3)
    t = s.tangent_from([1.0, 0.0, 0.0])
    theta = s.angles_to([1.0, 0.0, 0.0])
    ok = (theta > 0.1) & (theta < 3.0)
    assert np.allclose(np.einsum("ij,ij->i", t[ok], s.normals[ok]), 0.0, atol=1e-9)      # in the tangent plane
    assert (t[ok] @ np.array([1.0, 0.0, 0.0]) < 0.0).all()                              # away from the stagnation point
    assert m.params["layers"] == 2 and m.params["layer_thickness_m"] == 0.25e-3


def test_deactivation_updates_the_boundary(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points, layered_mesh.tets, element_layer=layered_mesh.element_layer.copy())
    s0 = m.surface()
    volume0 = m.volume()
    owners = s0.owner[s0.normals[:, 0] > 0.9]
    gone, new = m.deactivate(owners)
    assert m.n_active == m.n_elements - owners.size and gone.size == owners.size and new.size > 0
    s1 = m.surface()
    assert s1.n_patches == s0.n_patches - gone.size + new.size
    assert m.volume() == pytest.approx(volume0 - layered_mesh.element_volumes()[owners].sum())
    assert np.abs((s1.areas[:, None] * s1.normals).sum(axis=0)).max() < 1e-12 * s1.area      # still closed
    assert set(s1.face_ids[np.isin(s1.face_ids, s0.face_ids)]) == set(s0.face_ids) - set(gone)
    assert m.active[s1.owner].all()
    again = m.deactivate(owners)                                               # already dead: no change
    assert again[0].size == 0 and m.n_active == m.n_elements - owners.size
    assert np.isin(m.active_nodes(), np.unique(m.tets[m.active])).all()


def test_box_mesh():
    b = mesh.box_mesh(0.03, 0.004, 0.004, 0.001)
    assert b.n_nodes == 31 * 5 * 5 and b.n_elements == 30 * 4 * 4 * 6
    assert b.volume() == pytest.approx(0.03 * 0.004 * 0.004) and (b.element_volumes() > 0).all()
    s = b.surface()
    assert s.area == pytest.approx(2 * (0.03 * 0.004 * 2 + 0.004 * 0.004)) and s.n_patches == 2 * (30 * 4 * 2 * 2 + 4 * 4 * 2)
    assert (np.abs(np.abs(s.normals).max(axis=1) - 1.0) < 1e-12).all()      # axis-aligned faces
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: the five new tests fail (`layer_thicknesses`, `split_prisms`, `box_mesh` missing; `sphere_mesh` has no `layers`).

- [ ] **Step 3: Replace `reentry_model/mesh.py`**

```python
"""Tetrahedral sphere meshes (gmsh, graded, optional prism layers) and their surface geometry. Spec sections 5
(Step 2) and 5 (Step 3).

A VolumeMesh is node coordinates, tetrahedra, an active mask (Step 3: elements that have melted away are
deactivated and drop out of the surface, the operators and the mass) and the boundary triangles of the active set
(faces used by exactly one active element) with outward unit normals, centroids, areas and owner elements -- so any
gmsh volume mesh loads, tagged or not. A face table (every face of the mesh with its one or two elements) is built
once so that deactivation updates the boundary in O(faces) without re-sorting.

`sphere_mesh(R, h_surface, h_core)` is the Step 2 generator: element size growing linearly from h_surface at the
surface to h_core at the centre (gmsh Distance/Threshold field), a node embedded at the centre, cached by
(R, h_surface, h_core). With `layers = n > 0` (Step 3 default for melting runs) the graded gmsh mesh is generated for
the inner sphere R - T (T the total layer thickness) and n prism layers are built on it by radial projection of its
boundary triangulation: shells at R - T + t_(n-1) + ... , thicknesses t0 growth^k from the outside in (the outermost
layer is t0 thick), each prism split into three tetrahedra by the smallest-node-id diagonal rule (conforming across
prisms). The construction is exact on a sphere and needs no gmsh boundary-layer machinery (spike of 2026-09-20:
gmsh's extrudeBoundaryLayer is a geo-kernel function that cannot be attached to the OCC sphere without
re-parametrising; the radial construction gives the same layers deterministically). `box_mesh` is a structured
Kuhn-split box for the analytic (Stefan) tests.
"""
import os
from dataclasses import dataclass, field

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MESH_DIR = os.path.join(REPO_ROOT, "reentry_model_output", "meshes")
DEFAULT_H_SURFACE = 2.0e-3          # m; 2 mm / 8 mm gives 18.9 k nodes on the 100 mm sphere (measured 2026-09-18)
DEFAULT_H_CORE = 8.0e-3             # m
DEFAULT_LAYERS = 4                  # prism layers for melting runs (spec Step 3 section 5)
DEFAULT_LAYER_THICKNESS = 0.25e-3   # m, outermost layer
DEFAULT_LAYER_GROWTH = 2.0
_FACE_OF_VERTEX = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])     # face opposite each tet vertex


@dataclass
class SurfaceMesh:
    faces: np.ndarray            # (nf, 3) node ids of the boundary triangles ("patches")
    centroids: np.ndarray        # (nf, 3) m
    normals: np.ndarray          # (nf, 3) outward unit normals
    areas: np.ndarray            # (nf,) m2
    owner: np.ndarray = None     # (nf,) element owning each patch
    face_ids: np.ndarray = None  # (nf,) ids in the mesh's face table (stable across deactivations)

    @property
    def n_patches(self):
        return len(self.faces)

    @property
    def area(self):
        return float(self.areas.sum())

    def angles_to(self, v_hat):
        """theta per patch [rad]: the angle between the outward normal and v_hat, the direction the body moves in.
        The stagnation patch has its normal along v_hat (theta = 0); theta > pi/2 is leeward."""
        v = np.asarray(v_hat, dtype=float)
        v = v / np.linalg.norm(v)
        return np.arccos(np.clip(self.normals @ v, -1.0, 1.0))

    def facet_mean(self, nodal):
        """Mean of a nodal field over each patch's three nodes."""
        return np.asarray(nodal)[self.faces].mean(axis=1)

    def patch_toward(self, direction):
        """Index of the patch whose outward normal is closest to `direction`."""
        d = np.asarray(direction, dtype=float)
        return int(np.argmax(self.normals @ (d / np.linalg.norm(d))))

    def tangent_from(self, v_hat):
        """Unit surface direction away from the stagnation point (increasing theta) per patch; zero where undefined."""
        v = np.asarray(v_hat, dtype=float)
        v = v / np.linalg.norm(v)
        t = -v[None, :] + (self.normals @ v)[:, None] * self.normals
        n = np.linalg.norm(t, axis=1)
        return np.where(n[:, None] > 1e-12, t / np.maximum(n, 1e-300)[:, None], 0.0)

    def projected_area(self, v_hat):
        """Area projected on the plane normal to v_hat: sum of A max(0, n.v) (pi R^2 for a sphere)."""
        v = np.asarray(v_hat, dtype=float)
        return float((self.areas * np.maximum(0.0, self.normals @ (v / np.linalg.norm(v)))).sum())

    def edges(self):
        """Patch adjacency: (edge node pairs (ne, 2), patch i, patch j, edge lengths) for edges shared by two patches."""
        e = np.sort(self.faces[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
        patch = np.repeat(np.arange(self.n_patches), 3)
        order = np.lexsort((e[:, 1], e[:, 0]))
        e, patch = e[order], patch[order]
        same = np.all(e[1:] == e[:-1], axis=1)
        k = np.flatnonzero(same)
        return e[k], patch[k], patch[k + 1]


@dataclass
class VolumeMesh:
    points: np.ndarray                                   # (n, 3) m
    tets: np.ndarray                                     # (ne, 4) node ids
    params: dict = field(default_factory=dict)           # generator parameters / source file
    element_layer: np.ndarray = None                     # (ne,) prism layer index (0 = outermost) or -1 (graded core)
    active: np.ndarray = field(init=False, repr=False, default=None)
    _face_nodes: np.ndarray = field(init=False, repr=False, default=None)    # (nfaces, 3) node ids, one entry per unique face
    _face_elements: np.ndarray = field(init=False, repr=False, default=None) # (nfaces, 2) elements sharing the face, -1 if one
    _face_opposite: np.ndarray = field(init=False, repr=False, default=None) # (nfaces, 2) opposite vertex in each of those elements
    _face_count: np.ndarray = field(init=False, repr=False, default=None)    # active elements per face
    _element_faces: np.ndarray = field(init=False, repr=False, default=None) # (ne, 4) face ids of each element
    _boundary: np.ndarray = field(init=False, repr=False, default=None)      # face ids of the current boundary

    def __post_init__(self):
        self.points = np.array(self.points, dtype=float)
        self.tets = np.asarray(self.tets, dtype=np.int64)
        self.active = np.ones(len(self.tets), dtype=bool)
        if self.element_layer is None:
            self.element_layer = np.full(len(self.tets), -1, dtype=np.int64)
        self._build_face_table()

    def _build_face_table(self):
        faces = self.tets[:, _FACE_OF_VERTEX].reshape(-1, 3)                    # 4 per element, element-major
        opposite = self.tets.ravel()                                            # vertex opposite face 4e + v is tets[e, v]
        element = np.repeat(np.arange(len(self.tets)), 4)
        key = np.sort(faces, axis=1)
        _, first, inverse, counts = np.unique(key, axis=0, return_index=True, return_inverse=True, return_counts=True)
        inverse = inverse.ravel()
        nfaces = len(first)
        self._face_nodes = faces[first]
        self._face_elements = np.full((nfaces, 2), -1, dtype=np.int64)
        self._face_opposite = np.full((nfaces, 2), -1, dtype=np.int64)
        order = np.argsort(inverse, kind="stable")
        slot = np.zeros(nfaces, dtype=np.int64)
        for k in order:                                                         # <= 2 entries per face
            f = inverse[k]
            self._face_elements[f, slot[f]], self._face_opposite[f, slot[f]] = element[k], opposite[k]
            slot[f] += 1
        self._face_count = counts.astype(np.int64)
        self._element_faces = inverse.reshape(len(self.tets), 4)              # face id of each element's four faces
        self._boundary = np.flatnonzero(self._face_count == 1)

    @property
    def n_nodes(self):
        return len(self.points)

    @property
    def n_elements(self):
        return len(self.tets)

    @property
    def n_active(self):
        return int(self.active.sum())

    def element_volumes(self):
        x = self.points[self.tets]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        return np.abs(np.linalg.det(J)) / 6.0

    def volume(self):
        """Volume of the active elements."""
        return float(self.element_volumes()[self.active].sum())

    def element_faces(self, elements):
        """Face-table ids of each given element's four faces, (n, 4). Stable across deactivations, so a dying element
        can be asked which of its own faces have just become boundary (Step 3: where its melt film goes)."""
        return self._element_faces[np.asarray(elements, dtype=np.int64)]

    def deactivate(self, elements):
        """Remove elements from the active set; the boundary (and `surface()`) follows. Returns the face ids that
        stopped being boundary faces (their patches vanished) and those that became boundary faces."""
        elements = np.unique(np.asarray(elements, dtype=np.int64))
        elements = elements[self.active[elements]]
        if elements.size == 0:
            return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
        self.active[elements] = False
        faces = self._element_faces[elements].ravel()
        before = self._face_count[faces].copy()
        np.subtract.at(self._face_count, faces, 1)
        after = self._face_count[faces]
        gone, new = np.unique(faces[(before == 1) & (after == 0)]), np.unique(faces[(before == 2) & (after == 1)])
        self._boundary = np.flatnonzero(self._face_count == 1)
        return gone, new

    def surface(self):
        """The boundary patches of the active set with their current geometry."""
        f = self._boundary
        faces = self._face_nodes[f]
        active_slot = np.where(self.active[self._face_elements[f, 0]], 0, 1)
        owner = self._face_elements[f, active_slot]
        opposite = self._face_opposite[f, active_slot]
        a, b, c = (self.points[faces[:, i]] for i in range(3))
        n = np.cross(b - a, c - a)
        areas = 0.5 * np.linalg.norm(n, axis=1)
        n = n / (2.0 * areas)[:, None]
        centroids = (a + b + c) / 3.0
        inward = np.einsum("ij,ij->i", n, centroids - self.points[opposite]) < 0.0
        n[inward] *= -1.0
        return SurfaceMesh(faces, centroids, n, areas, owner, f.copy())

    def boundary_nodes(self):
        return np.unique(self._face_nodes[self._boundary])

    def active_nodes(self):
        """Nodes belonging to at least one active element."""
        return np.unique(self.tets[self.active])

    def centre_node(self):
        """Node nearest the origin (the embedded centre node of a generated sphere)."""
        return int(np.argmin(np.linalg.norm(self.points, axis=1)))


def boundary_faces(tets):
    """Faces used by exactly one tetrahedron, as (nf, 3) node ids, and the opposite vertex of that tetrahedron."""
    faces = tets[:, _FACE_OF_VERTEX].reshape(-1, 3)
    opposite = np.repeat(tets, 4, axis=0)[np.arange(4 * len(tets)), np.tile(np.arange(4), len(tets))]
    _, first, counts = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_counts=True)
    keep = first[counts == 1]
    return faces[keep], opposite[keep]


def mesh_file_name(radius, h_surface, h_core):
    return "sphere_R{:.3f}mm_hs{:.3f}mm_hc{:.3f}mm.msh".format(radius * 1e3, h_surface * 1e3, h_core * 1e3)


def generate_sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR):
    """gmsh sphere of `radius` graded from h_surface (surface) to h_core (centre); returns the cached .msh path."""
    os.makedirs(mesh_dir, exist_ok=True)
    path = os.path.join(mesh_dir, mesh_file_name(radius, h_surface, h_core))
    if os.path.isfile(path):
        return path
    import gmsh                                   # imported here: only mesh generation needs gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("sphere")
        vol = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, radius)
        centre = gmsh.model.occ.addPoint(0.0, 0.0, 0.0, h_core)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.embed(0, [centre], 3, vol)
        surfaces = [s[1] for s in gmsh.model.getBoundary([(3, vol)], oriented=False)]
        gmsh.model.addPhysicalGroup(3, [vol], tag=1, name="body")
        gmsh.model.addPhysicalGroup(2, surfaces, tag=2, name="surface")
        dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", surfaces)
        gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)
        thr = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thr, "InField", dist)
        gmsh.model.mesh.field.setNumber(thr, "SizeMin", h_surface)
        gmsh.model.mesh.field.setNumber(thr, "SizeMax", h_core)
        gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(thr, "DistMax", radius)
        gmsh.model.mesh.field.setAsBackgroundMesh(thr)
        for option in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
            gmsh.option.setNumber(option, 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)          # HXT
        gmsh.model.mesh.generate(3)
        tmp = path[:-4] + ".tmp.msh"           # gmsh.write needs the .msh extension to pick the format
        gmsh.write(tmp)
    finally:
        gmsh.finalize()
    os.replace(tmp, path)                      # atomic: an interrupted generation never leaves a reusable partial .msh
    return path


def load_mesh(path):
    """Any gmsh/meshio volume mesh with tetrahedra (other cell types are ignored)."""
    import meshio
    m = meshio.read(path)
    tets = [c.data for c in m.cells if c.type == "tetra"]
    if not tets:
        raise ValueError("no tetrahedra in {}".format(path))
    return VolumeMesh(np.asarray(m.points, dtype=float), np.vstack(tets), {"path": os.path.abspath(path)})


def layer_thicknesses(layers, layer_thickness, growth):
    """Thickness of each prism layer from the outside in: t0, t0 g, t0 g^2, ..."""
    return layer_thickness * growth ** np.arange(layers)


def split_prisms(bottom, top):
    """Three tetrahedra per prism (bottom (np, 3) and top (np, 3) node ids, top[k] above bottom[k]) with the diagonal
    of every quad face drawn from its smallest node id -- so two prisms sharing a quad face split it identically."""
    bottom, top = np.asarray(bottom, dtype=np.int64), np.asarray(top, dtype=np.int64)
    # rotate each triangle so that its smallest bottom id comes first (the top follows the same rotation)
    shift = np.argmin(bottom, axis=1)
    idx = (shift[:, None] + np.arange(3)[None, :]) % 3
    b = np.take_along_axis(bottom, idx, axis=1)
    t = np.take_along_axis(top, idx, axis=1)
    b0, b1, b2, t0, t1, t2 = b[:, 0], b[:, 1], b[:, 2], t[:, 0], t[:, 1], t[:, 2]
    # diagonals from b0 on the faces (b0 b1 t1 t0) and (b2 b0 t0 t2); on (b1 b2 t2 t1) from min(b1, b2)
    case = b1 < b2
    tets = np.where(case[:, None, None],
                    np.stack([np.stack([b0, b1, b2, t2], 1), np.stack([b0, b1, t2, t1], 1), np.stack([b0, t1, t2, t0], 1)], 1),
                    np.stack([np.stack([b0, b1, b2, t1], 1), np.stack([b0, b2, t1, t2], 1), np.stack([b0, t1, t2, t0], 1)], 1))
    return tets.reshape(-1, 4)


def add_prism_layers(inner, radius, layers, layer_thickness, growth):
    """Prism layers by radial projection of the inner sphere mesh's boundary triangulation out to `radius`
    (thicknesses t0 growth^k from the outside in). Returns the layered VolumeMesh (inner elements first, then the
    layers from the inside out; `element_layer` is 0 for the outermost layer, -1 for the graded core)."""
    t = layer_thicknesses(layers, layer_thickness, growth)
    r_inner = radius - t.sum()
    faces, _ = boundary_faces(inner.tets)
    shell_nodes = np.unique(faces)
    ns = len(shell_nodes)
    local = -np.ones(inner.n_nodes, dtype=np.int64)
    local[shell_nodes] = np.arange(ns)
    lf = local[faces]                                        # faces in shell-local ids
    direction = inner.points[shell_nodes]
    direction = direction / np.linalg.norm(direction, axis=1)[:, None]
    radii = r_inner + np.cumsum(t[::-1])                     # shell radii from the inside out (outermost = radius)
    shell_ids = [shell_nodes] + [inner.n_nodes + k * ns + np.arange(ns) for k in range(layers)]
    points = np.vstack([inner.points] + [direction * r for r in radii])
    tets, layer_index = [inner.tets], [np.full(inner.n_elements, -1, dtype=np.int64)]
    for k in range(layers):
        new = split_prisms(shell_ids[k][lf], shell_ids[k + 1][lf])
        tets.append(new)
        layer_index.append(np.full(len(new), layers - 1 - k, dtype=np.int64))
    return VolumeMesh(points, np.vstack(tets), dict(inner.params), np.concatenate(layer_index))


def sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR, layers=0,
                layer_thickness=DEFAULT_LAYER_THICKNESS, growth=DEFAULT_LAYER_GROWTH):
    """Graded gmsh sphere (layers = 0, Step 2) or the graded inner sphere with `layers` prism layers on it."""
    if layers < 0 or layer_thickness <= 0.0 or growth <= 0.0:
        raise ValueError("layers must be >= 0 and layer_thickness, growth > 0")
    total = float(layer_thicknesses(layers, layer_thickness, growth).sum()) if layers else 0.0
    if total >= 0.5 * radius:
        raise ValueError("the prism layers ({:.3g} m) must be thinner than half the radius".format(total))
    mesh = load_mesh(generate_sphere_mesh(radius - total, h_surface, h_core, mesh_dir))
    if layers:
        mesh = add_prism_layers(mesh, radius, layers, layer_thickness, growth)
    mesh.params.update({"radius_m": radius, "h_surface_m": h_surface, "h_core_m": h_core, "layers": layers,
                        "layer_thickness_m": layer_thickness, "layer_growth": growth})
    return mesh


def box_mesh(lx, ly, lz, h):
    """Structured box [0, lx] x [0, ly] x [0, lz] of cubes of side ~h, each split into six tetrahedra (Kuhn
    decomposition, conforming). For the analytic tests (Stefan front along x with the x = 0 face heated)."""
    n = [max(1, int(round(L / h))) for L in (lx, ly, lz)]
    xs = [np.linspace(0.0, L, k + 1) for L, k in zip((lx, ly, lz), n)]
    X, Y, Z = np.meshgrid(*xs, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    nid = np.arange(points.shape[0]).reshape(n[0] + 1, n[1] + 1, n[2] + 1)
    c = nid[:-1, :-1, :-1].ravel()
    dx, dy, dz = nid[1, 0, 0] - nid[0, 0, 0], nid[0, 1, 0] - nid[0, 0, 0], 1
    v = lambda i, j, k: c + i * dx + j * dy + k * dz
    corners = [v(0, 0, 0), v(1, 0, 0), v(1, 1, 0), v(0, 1, 0), v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)]
    kuhn = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
    tets = np.vstack([np.column_stack([corners[a], corners[b], corners[cc], corners[d]]) for a, b, cc, d in kuhn])
    return VolumeMesh(points, tets, {"box": (lx, ly, lz), "h": h, "path": None})
```


- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: 10 passed (the five Step 2 tests unchanged; the layered fixture builds a 4 mm / 20 mm inner sphere with two layers, ~5 k nodes). Then the unit tier: `"$PY" -m pytest -m "not drama and not reference" -q` — everything that used `mesh.surface()` still passes (the boundary now comes from the face table; the coarse-mesh fixtures are unchanged).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/mesh.py tests/test_reentry_model_mesh.py
git commit -m "Add the dense surface band, the active element set, the derived surface and the box mesh (Step 3 Task 1)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

