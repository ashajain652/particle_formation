# Mesh: dense band and derived surface (Step 3 Task 1 amendment) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give `mesh.py` a dense surface band set by a gmsh size field, a recession-displaced *derived* surface that the surface physics can read while the solver's own geometry stands still, and the two coarseness diagnostics the remesh triggers need.

**Architecture:** Four additive changes to one module. `SurfaceMesh` gains the node coordinates it was built from, so it can measure its own edges and smooth its own normals. `VolumeMesh` gains `derived_points` — a displaced copy that **nothing in `thermal/` ever sees** — and `surface(derived=True)` to read it. The gmsh `Threshold` field gains a non-zero `DistMin`, which is the whole of the dense band. Nothing is removed: prism layers keep working, just no longer by default.

**Tech Stack:** Python 3.12, numpy, scipy (`cKDTree`, `scipy.sparse`), gmsh 4.15.2, meshio, pytest. scipy and gmsh are imported lazily inside functions, matching the module's existing style.

**Spec:** `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md` (§6 the band, §8 the derived surface). Sub-plan: `docs/superpowers/plans/melt-spraying-subplans/01-mesh-active-set.md`, amendment section. Measured facts 38–45 in `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`.

## Global Constraints

- **Work in the prototype, not the package.** All edits go to `prototype/proto3/reentry_model/mesh.py` and `prototype/proto3/tests/test_reentry_model_mesh.py`. `reentry_model/mesh.py` at the repo root is the Step 2 version and is **not** touched — Step 3 has never landed there. `prototype/README.md` is the authority: edit `proto3/`, run its tests there, regenerate the plan later.
- **Interpreter:** `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`. There is no virtualenv and no `pip install -e .`; run pytest from inside `prototype/proto3` so `reentry_model` resolves.
- **Run tests as:** `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest tests/test_reentry_model_mesh.py -q`
- **Never modify `VolumeMesh.points`.** `derived_points` is a separate array. This is the single invariant the whole design rests on (spec §5, Design B): the thermal solver reads element geometry from `points`, so if `points` moves, every precomputed element matrix is invalidated and skfem's assembly goes from 0.004 s to 0.074 s (measured fact 45).
- **Additive only.** `layers`, `layer_thicknesses`, `split_prisms`, `add_prism_layers`, `deactivate`, the face table and `box_mesh` all keep working unchanged. The five existing Step 3 mesh tests must still pass untouched.
- **Comment density and naming must match the surrounding file** — terse trailing comments carrying a measured number or a reason, docstrings that say *why*, not *what*.
- **No new module-level imports.** `scipy` and `gmsh` are imported inside the function that needs them.
- **Commit messages** are a scope prefix and a sentence describing the behaviour that changed, ending with the attribution line shown in each task's commit step.

## Review Focus

Five things the spec implies, that a reasonable person would expect to work, and that no task's main tests would otherwise exercise. Each has a test assigned to the task that owns the code.

1. **`band=0` must reproduce Step 2's grading exactly.** The pinned verification devices of measured fact 12 depend on the old behaviour; a band that silently changes them invalidates the committed reference numbers. → Task 2.
2. **`band` at or beyond the radius must not produce a degenerate size field.** `DistMin == DistMax` is a division by zero inside gmsh's Threshold; the sensible behaviour is a uniformly fine mesh. → Task 2.
3. **A banded and an unbanded mesh of the same radius and sizes must not collide in the cache.** They are different meshes; sharing one `.msh` would silently serve the wrong one and the error would surface as an inexplicable physics change. → Task 2.
4. **`element_depths` on a fully dead mesh must not raise.** `surface()` returns zero patches, and `cKDTree` on an empty point set raises. Demise is a normal end state, not an error. → Task 1.
5. **A degenerate patch must not silently produce NaN normals.** A derived surface can collapse a patch to zero area; `smoothed_normals` dividing by that length yields NaN, which propagates into θ and then into the heating without any error being raised. → Task 4.

---

### Task 1: Surface coordinates, edge lengths and element depths

**Files:**
- Modify: `prototype/proto3/reentry_model/mesh.py` — `SurfaceMesh` dataclass (line 36–43), `VolumeMesh.surface` (line 178)
- Test: `prototype/proto3/tests/test_reentry_model_mesh.py` (append)

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces: `SurfaceMesh.points` (the `(n, 3)` coordinate array the surface was built from, a reference not a copy); `SurfaceMesh.edge_lengths() -> (nf, 3)`; `VolumeMesh.element_depths() -> (ne,)` in metres, `-1.0` for dead elements. Task 2's band test and Task 4's smoothing both need `points`; Task 17's remesh triggers need both accessors.

- [ ] **Step 1: Write the failing tests**

Append to `prototype/proto3/tests/test_reentry_model_mesh.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Amendment of 2026-09-27: coarseness diagnostics, the dense band and the derived surface


def test_surface_carries_its_coordinates_and_measures_its_edges(layered_mesh):
    s = layered_mesh.surface()
    assert s.points is layered_mesh.points                       # a reference, not a copy: no per-step duplication
    e = s.edge_lengths()
    assert e.shape == (s.n_patches, 3) and (e > 0).all()
    # an equilateral patch of area A has edge sqrt(4A/sqrt(3)); the layered mesh's surface is the 4 mm inner mesh
    assert e.mean() == pytest.approx(math.sqrt(4.0 * s.areas.mean() / math.sqrt(3.0)), rel=0.15)
    assert e.max() < 3.0 * e.mean()                              # no wildly stretched patch on a fresh sphere


def test_element_depths_run_from_zero_at_the_surface_to_the_radius(layered_mesh):
    d = layered_mesh.element_depths()
    assert d.shape == (layered_mesh.n_elements,)
    assert d.min() >= 0.0 and d.min() < 1.0e-3                   # the outermost prism sits on the surface
    assert d.max() == pytest.approx(R, rel=0.2)                  # the deepest element is near the centre
    radii = np.linalg.norm(layered_mesh.points[layered_mesh.tets].mean(axis=1), axis=1)
    assert np.corrcoef(d, R - radii)[0, 1] > 0.99                # depth is the radial complement, on a sphere


def test_element_depths_mark_dead_elements_and_survive_total_demise(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    dead = m.surface().owner[:20]
    m.deactivate(dead)
    d = m.element_depths()
    assert (d[dead] == -1.0).all() and (d[m.active] >= 0.0).all()
    m.deactivate(np.arange(m.n_elements))                        # demise: nothing left, and no patches either
    assert m.n_active == 0 and m.surface().n_patches == 0
    assert (m.element_depths() == -1.0).all()                    # review focus 4: must not raise on an empty tree
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest tests/test_reentry_model_mesh.py -q -k "carries_its_coordinates or element_depths"`
Expected: 3 failures — `SurfaceMesh` has no attribute `points`, no `edge_lengths`; `VolumeMesh` has no `element_depths`.

- [ ] **Step 3: Add the `points` field and `edge_lengths` to `SurfaceMesh`**

Add `points` as the **last** field of the dataclass, after `face_ids`, so every existing positional construction keeps working:

```python
    face_ids: np.ndarray = None  # (nf,) ids in the mesh's face table (stable across deactivations)
    points: np.ndarray = None    # the (n, 3) coordinates these patches index -- `points` or `derived_points`
```

and this method, next to `facet_mean`:

```python
    def edge_lengths(self):
        """The three edge lengths of every patch, (nf, 3) m. `edge_lengths().max()` is the surface coarseness the
        history records and the second remesh trigger watches: on a fresh 100 mm sphere it is 2.16 mm, and once the
        erosion front is into the graded core it reaches 8 mm (measured fact 39)."""
        x = self.points[self.faces]
        return np.linalg.norm(x[:, [1, 2, 0]] - x, axis=2)
```

- [ ] **Step 4: Pass the coordinates from `surface()`**

In `VolumeMesh.surface`, the return becomes:

```python
        return SurfaceMesh(faces, centroids, n, areas, owner, f.copy(), self.points)
```

- [ ] **Step 5: Add `element_depths` to `VolumeMesh`**

Next to `active_nodes`:

```python
    def element_depths(self):
        """Distance from each active element's centroid to the nearest boundary patch [m]; -1 for dead elements.
        The first remesh trigger fires when the deepest element that still owns a patch is more than half the band
        down, i.e. when the erosion front has eaten through half the fine mesh. Measured against boundary *centroids*
        rather than the triangles themselves, which over-reads by at most half a patch edge -- immaterial against a
        15 mm band, and it keeps this a single cKDTree query."""
        from scipy.spatial import cKDTree            # imported here, like gmsh and meshio: not every mesh needs it
        depths = np.full(len(self.tets), -1.0)
        live = np.flatnonzero(self.active)
        patches = self.surface()
        if live.size == 0 or patches.n_patches == 0:             # demised: no elements, or none with a boundary
            return depths
        centres = self.points[self.tets[live]].mean(axis=1)
        depths[live] = cKDTree(patches.centroids).query(centres)[0]
        return depths
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: PASS, including the five pre-existing Step 3 mesh tests, unchanged.

- [ ] **Step 7: Commit**

```bash
git add prototype/proto3/reentry_model/mesh.py prototype/proto3/tests/test_reentry_model_mesh.py
git commit -m "mesh: carry the surface's coordinates, and measure patch edges and element depths

The two coarseness diagnostics the remesh triggers read, and the node
coordinates a SurfaceMesh needs to measure or smooth itself.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: The dense band

**Files:**
- Modify: `prototype/proto3/reentry_model/mesh.py` — constants (line 28–32), `mesh_file_name` (line 215), `generate_sphere_mesh` (line 219), `sphere_mesh` (line 317)
- Test: `prototype/proto3/tests/test_reentry_model_mesh.py` (append)

**Interfaces:**
- Consumes: `VolumeMesh.element_depths()` from Task 1.
- Produces: `DEFAULT_BAND = 15.0e-3`, `DEFAULT_BAND_RAMP = 8.0e-3`, `DEFAULT_LAYERS = 0`; `mesh_file_name(radius, h_surface, h_core, band=0.0, ramp=DEFAULT_BAND_RAMP)`; `generate_sphere_mesh(..., band=DEFAULT_BAND, ramp=DEFAULT_BAND_RAMP)`; `sphere_mesh(..., band=DEFAULT_BAND, ramp=DEFAULT_BAND_RAMP)`; `params["band_m"]` and `params["band_ramp_m"]`.

**Why:** the prism stack keeps the surface resolved for 3.75 mm, but the nose recedes 77 mm (fact 14), so it covers 4.9 % of the recession and the graded core — 2.74 mm at the surface rising to 8.32 mm at the centre (fact 39) — becomes the surface. A distance field is well defined at every interior point of any shape, so unlike a prism offset it survives remeshing into a concave dish.

- [ ] **Step 1: Write the failing tests**

```python
@pytest.fixture(scope="module")
def banded_mesh(tmp_path_factory):
    """50 mm sphere, 4 mm/20 mm sizes with an 8 mm dense band (small numbers to keep the test cheap; production is
    2.16 mm/8 mm with a 15 mm band)."""
    return mesh.sphere_mesh(R, 4e-3, 20e-3, str(tmp_path_factory.mktemp("meshes")), layers=0, band=8e-3, ramp=8e-3)


def test_band_keeps_the_cells_fine_through_its_whole_depth(banded_mesh):
    d = banded_mesh.element_depths()
    edge = (banded_mesh.element_volumes() / 0.117851) ** (1.0 / 3.0)    # edge of the regular tet of equal volume
    inside, outside = edge[(d >= 0) & (d < 6e-3)], edge[d > 20e-3]
    assert inside.mean() == pytest.approx(4e-3, rel=0.4)                # the band is at h_surface throughout
    assert outside.mean() > 1.6 * inside.mean()                         # and only then does it coarsen
    assert banded_mesh.params["band_m"] == 8e-3 and banded_mesh.params["band_ramp_m"] == 8e-3


def test_band_zero_reproduces_the_step_2_grading(tmp_path):
    """Review focus 1: the pinned verification devices of fact 12 depend on the old behaviour exactly, so band = 0
    must still grade from h_surface at the surface to h_core at the centre rather than hold h_surface anywhere."""
    m = mesh.sphere_mesh(0.01, 2e-3, 6e-3, str(tmp_path), band=0.0)
    d = m.element_depths()
    edge = (m.element_volumes() / 0.117851) ** (1.0 / 3.0)
    near, deep = edge[(d >= 0) & (d < 2e-3)], edge[d > 5e-3]
    assert near.size and deep.size
    assert deep.mean() > 1.5 * near.mean()                              # graded, not banded
    assert m.params["band_m"] == 0.0
    assert "_b0.000mm_" in mesh.mesh_file_name(0.01, 2e-3, 6e-3, 0.0, 8e-3)


def test_band_beyond_the_radius_gives_a_uniform_mesh_and_does_not_divide_by_zero(tmp_path):
    """Review focus 2: DistMin == DistMax is degenerate inside gmsh's Threshold field."""
    m = mesh.sphere_mesh(0.01, 2e-3, 8e-3, str(tmp_path), band=0.05, ramp=8e-3)
    edge = (m.element_volumes() / 0.117851) ** (1.0 / 3.0)
    assert np.isfinite(edge).all() and (edge > 0).all()
    assert np.percentile(edge, 90) < 2.0 * np.percentile(edge, 10)      # uniformly fine; percentiles, since gmsh leaves slivers
    with pytest.raises(ValueError):
        mesh.sphere_mesh(0.01, 2e-3, 8e-3, str(tmp_path), band=-1e-3)
    with pytest.raises(ValueError):
        mesh.sphere_mesh(0.01, 2e-3, 8e-3, str(tmp_path), ramp=0.0)


def test_banded_and_unbanded_meshes_do_not_share_a_cache_file():
    """Review focus 3: sharing one .msh would silently serve the wrong mesh."""
    assert mesh.mesh_file_name(R, 2e-3, 8e-3, 15e-3, 8e-3) != mesh.mesh_file_name(R, 2e-3, 8e-3, 0.0, 8e-3)
    assert mesh.mesh_file_name(R, 2e-3, 8e-3, 15e-3, 8e-3) != mesh.mesh_file_name(R, 2e-3, 8e-3, 15e-3, 4e-3)
    assert mesh.DEFAULT_LAYERS == 0 and mesh.DEFAULT_BAND == 15e-3 and mesh.DEFAULT_BAND_RAMP == 8e-3
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest tests/test_reentry_model_mesh.py -q -k band`
Expected: failures — `sphere_mesh() got an unexpected keyword argument 'band'`, and `DEFAULT_BAND` missing.

- [ ] **Step 3: Add the constants**

Replace `DEFAULT_LAYERS = 4` and add two below it:

```python
DEFAULT_LAYERS = 0                  # prism layers: superseded as the default by the dense band (amendment 2026-09-27)
DEFAULT_LAYER_THICKNESS = 0.25e-3   # m, outermost layer
DEFAULT_LAYER_GROWTH = 2.0
DEFAULT_BAND = 15.0e-3              # m; cells stay at h_surface this far down, so the surface survives 7.5 mm of
                                    # recession before the mesh must be rebuilt. 310 k tets on the 100 mm sphere
DEFAULT_BAND_RAMP = 8.0e-3          # m; the depth over which the size then grows from h_surface to h_core
```

- [ ] **Step 4: Put the band in the cache key**

```python
def mesh_file_name(radius, h_surface, h_core, band=0.0, ramp=DEFAULT_BAND_RAMP):
    """Cache key. The band and its ramp are part of it: two meshes of the same radius and sizes but different bands
    are different meshes, and sharing one file would silently serve the wrong one."""
    return "sphere_R{:.3f}mm_hs{:.3f}mm_hc{:.3f}mm_b{:.3f}mm_r{:.3f}mm.msh".format(
        radius * 1e3, h_surface * 1e3, h_core * 1e3, band * 1e3, ramp * 1e3)
```

- [ ] **Step 5: Give the size field a non-zero `DistMin`**

Change the signature and the two field lines in `generate_sphere_mesh`:

```python
def generate_sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR,
                         band=DEFAULT_BAND, ramp=DEFAULT_BAND_RAMP):
    """gmsh sphere of `radius`, at h_surface throughout the outer `band` and then graded to h_core over `ramp`;
    returns the cached .msh path. `band = 0` is the Step 2 grading, which the pinned verification devices use."""
    os.makedirs(mesh_dir, exist_ok=True)
    path = os.path.join(mesh_dir, mesh_file_name(radius, h_surface, h_core, band, ramp))
```

and, replacing the `DistMin`/`DistMax` pair:

```python
        d_min = min(band, radius)                      # a band at or past the centre means uniformly fine
        d_max = max(min(band + ramp, radius), d_min + 1e-9)     # never equal: Threshold divides by the difference
        gmsh.model.mesh.field.setNumber(thr, "DistMin", d_min)
        gmsh.model.mesh.field.setNumber(thr, "DistMax", d_max)
```

- [ ] **Step 6: Thread it through `sphere_mesh`**

```python
def sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR, layers=0,
                layer_thickness=DEFAULT_LAYER_THICKNESS, growth=DEFAULT_LAYER_GROWTH, band=DEFAULT_BAND,
                ramp=DEFAULT_BAND_RAMP):
    """Graded gmsh sphere with a dense outer band, optionally with `layers` prism layers on top of it."""
    if layers < 0 or layer_thickness <= 0.0 or growth <= 0.0:
        raise ValueError("layers must be >= 0 and layer_thickness, growth > 0")
    if band < 0.0 or ramp <= 0.0:
        raise ValueError("band must be >= 0 and ramp > 0")
    total = float(layer_thicknesses(layers, layer_thickness, growth).sum()) if layers else 0.0
    if total >= 0.5 * radius:
        raise ValueError("the prism layers ({:.3g} m) must be thinner than half the radius".format(total))
    mesh = load_mesh(generate_sphere_mesh(radius - total, h_surface, h_core, mesh_dir, band, ramp))
    if layers:
        mesh = add_prism_layers(mesh, radius, layers, layer_thickness, growth)
    mesh.params.update({"radius_m": radius, "h_surface_m": h_surface, "h_core_m": h_core, "layers": layers,
                        "layer_thickness_m": layer_thickness, "layer_growth": growth, "band_m": band,
                        "band_ramp_m": ramp})
    return mesh
```

- [ ] **Step 7: Run the whole mesh module**

Run: `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: PASS.

**Two things will look alarming and are expected.** First, `band` now defaults to 15 mm in `sphere_mesh`, so
**every existing fixture silently becomes a banded mesh** — `coarse_sphere_mesh`, `layered_mesh`, and the fixtures in
`test_reentry_model_melting.py`. Their tests check geometry (closed surface, exact volume and area, shell radii,
three tets per prism), not element counts, so they should still pass; if one fails, read it carefully before
changing it, because a genuine break here is a genuine break. Second, the cache key changed, so **every cached
`.msh` in `reentry_model_output/meshes/` is now a miss and regenerates** — the 100 mm sphere takes a couple of
minutes the first time. Neither is a defect.

The pre-existing `test_generation_is_cached` still passes: it calls `generate_sphere_mesh` positionally with four
arguments and the new ones are keyword-defaulted.

- [ ] **Step 8: Commit**

```bash
git add prototype/proto3/reentry_model/mesh.py prototype/proto3/tests/test_reentry_model_mesh.py
git commit -m "mesh: keep the cells fine through a dense band, not just a prism skin

The prism stack resolves 3.75 mm of a 77 mm recession, so the graded core
becomes the surface and the recession quantum reaches 6 mm. A non-zero DistMin
on the existing Threshold field holds h_surface through the whole band, and
unlike a prism offset a distance field survives remeshing into a concave dish.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: The derived surface

**Files:**
- Modify: `prototype/proto3/reentry_model/mesh.py` — `VolumeMesh` fields (line 94–106), `__post_init__` (line 107), `surface` (line 178)
- Test: `prototype/proto3/tests/test_reentry_model_mesh.py` (append)

**Interfaces:**
- Consumes: `SurfaceMesh.points` from Task 1.
- Produces: `VolumeMesh.derived_points` — an `(n, 3)` array initialised to a copy of `points`, which Task 16 writes and nothing else does; `VolumeMesh.surface(derived=False)`.

**Why:** the shape is currently a binary function of which elements are alive, so the surface jumps by a whole cell. Writing the recession into a *separate* array lets the shape move continuously for everything that reads it, while `points` — and therefore every element matrix the thermal solver precomputes — never changes.

- [ ] **Step 1: Write the failing tests**

```python
def test_derived_points_start_equal_and_are_independent(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    assert m.derived_points is not m.points and np.array_equal(m.derived_points, m.points)
    m.derived_points[0] += 1.0
    assert not np.array_equal(m.derived_points, m.points)          # writing one must not write the other


def test_moving_the_derived_surface_leaves_the_solver_geometry_untouched(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    p0, v0, e0 = m.points.copy(), m.volume(), m.element_volumes().copy()
    a0 = m.surface().area
    bn = m.boundary_nodes()
    m.derived_points[bn] *= 0.98                                    # a 2 % uniform recession, derived only
    assert np.array_equal(m.points, p0)                             # bit-identical: the solver sees nothing
    assert m.volume() == v0 and np.array_equal(m.element_volumes(), e0)
    assert m.surface().area == pytest.approx(a0)                    # the default surface is still the real one
    d = m.surface(derived=True)
    assert d.area == pytest.approx(0.98 ** 2 * a0, rel=1e-6)        # and the derived one shrank by exactly r^2
    assert d.points is m.derived_points
    assert np.array_equal(d.faces, m.surface().faces)               # same patches, different coordinates


def test_derived_surface_normals_stay_outward(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    m.derived_points[m.boundary_nodes()] *= 0.9
    d = m.surface(derived=True)
    outward = np.einsum("ij,ij->i", d.normals, d.centroids / np.linalg.norm(d.centroids, axis=1)[:, None])
    assert (outward > 0.5).all()                                    # still pointing out of the shrunken sphere
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest tests/test_reentry_model_mesh.py -q -k derived`
Expected: 3 failures — `VolumeMesh` has no attribute `derived_points`; `surface()` takes no `derived` argument.

- [ ] **Step 3: Add the field**

After `_boundary` in the dataclass:

```python
    derived_points: np.ndarray = field(init=False, repr=False, default=None)   # recession-displaced copy of points
```

and at the end of `__post_init__`, **after** `self.points` is coerced to float:

```python
        self.derived_points = self.points.copy()
```

- [ ] **Step 4: Let `surface()` read either array**

```python
    def surface(self, derived=False):
        """The boundary patches of the active set with their current geometry.

        With `derived=True` the geometry comes from `derived_points`, the recession-displaced copy, instead of from
        `points`. That is the whole of Design B: the shape every consumer of the surface reads -- the heating angle
        theta, the film tangents, the projected area, the nose fit, the pictures -- moves continuously with the
        recession, while `points` never moves, so the thermal solver's element matrices stay valid and skfem keeps
        rescaling them at 0.004 s instead of recomputing at 0.074 s (measured fact 45)."""
        f = self._boundary
        faces = self._face_nodes[f]
        pts = self.derived_points if derived else self.points
        active_slot = np.where(self.active[self._face_elements[f, 0]], 0, 1)
        owner = self._face_elements[f, active_slot]
        opposite = self._face_opposite[f, active_slot]
        a, b, c = (pts[faces[:, i]] for i in range(3))
        n = np.cross(b - a, c - a)
        areas = 0.5 * np.linalg.norm(n, axis=1)
        n = n / (2.0 * areas)[:, None]
        centroids = (a + b + c) / 3.0
        inward = np.einsum("ij,ij->i", n, centroids - pts[opposite]) < 0.0
        n[inward] *= -1.0
        return SurfaceMesh(faces, centroids, n, areas, owner, f.copy(), pts)
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: PASS, all tests including the pre-existing ones.

- [ ] **Step 6: Commit**

```bash
git add prototype/proto3/reentry_model/mesh.py prototype/proto3/tests/test_reentry_model_mesh.py
git commit -m "mesh: a derived surface that moves while the solver's geometry stands still

derived_points carries the recession for everything that reads the shape;
points never changes, so every precomputed element matrix stays valid.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: Smoothed normals

**Files:**
- Modify: `prototype/proto3/reentry_model/mesh.py` — constants, `SurfaceMesh` (add one method)
- Test: `prototype/proto3/tests/test_reentry_model_mesh.py` (append)

**Interfaces:**
- Consumes: `SurfaceMesh.points` (Task 1), `VolumeMesh.surface(derived=True)` (Task 3).
- Produces: `TAUBIN_PASSES = 8`, `TAUBIN_LAMBDA = 0.53`, `TAUBIN_MU = -0.55`, `TAUBIN_CLAMP = 0.35`; `SurfaceMesh.smoothed_normals(passes=TAUBIN_PASSES, lam=TAUBIN_LAMBDA, mu=TAUBIN_MU, clamp=TAUBIN_CLAMP) -> (nf, 3)` unit vectors, outward, same order as `normals`.

**Why:** a staircase step left by element death is one cell tall and one cell wide, so its patch normal deviates from the true surface by up to 45°, and θ is what the Lees distribution and the film tangents read. Taubin's shrinking-then-inflating pair removes that without the net shrinkage a plain Laplacian causes, and being purely local it assumes nothing about the shape.

- [ ] **Step 1: Write the failing tests**

```python
def test_smoothed_normals_halve_the_error_a_staircase_introduces(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    exact = m.surface()                                              # the analytic sphere's own normals
    rng = np.random.default_rng(0)
    bn = m.boundary_nodes()
    m.derived_points[bn] += rng.normal(0.0, 3.0e-4, (bn.size, 3))    # cell-scale roughness, like element death
    rough = m.surface(derived=True)
    ang = lambda n: np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", n, exact.normals), -1.0, 1.0)))
    assert ang(rough.normals).mean() > 5.0                           # the staircase really is that bad
    assert ang(rough.smoothed_normals()).mean() < 0.5 * ang(rough.normals).mean()


def test_smoothed_normals_are_unit_outward_and_leave_a_clean_sphere_alone(layered_mesh):
    s = layered_mesh.surface()
    n = s.smoothed_normals()
    assert n.shape == s.normals.shape
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0, atol=1e-9)
    assert (np.einsum("ij,ij->i", n, s.normals) > 0.0).all()         # never flipped inward
    ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", n, s.normals), -1.0, 1.0)))
    assert ang.mean() < 2.0                                          # an already-smooth sphere barely moves


def test_smoothed_normals_never_return_nan_on_a_degenerate_patch(layered_mesh):
    """Review focus 5: a collapsed patch divides by zero, and NaN would reach theta and the heating unflagged."""
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    f = m.surface().faces[0]
    m.derived_points[f[1]] = m.derived_points[f[0]]                  # two nodes coincident: zero-area patch
    n = m.surface(derived=True).smoothed_normals()
    assert np.isfinite(n).all()
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest tests/test_reentry_model_mesh.py -q -k smoothed`
Expected: 3 failures — `SurfaceMesh` has no attribute `smoothed_normals`.

- [ ] **Step 3: Add the constants**

Below `DEFAULT_BAND_RAMP`:

```python
TAUBIN_PASSES = 8                   # lambda/mu smoothing of the surface read by the heating (amendment 2026-09-27)
TAUBIN_LAMBDA = 0.53
TAUBIN_MU = -0.55                   # slightly stronger than lambda: the pair passes low frequencies without shrinking
TAUBIN_CLAMP = 0.35                 # of the mean patch edge; a node can then never move past its neighbours
```

- [ ] **Step 4: Add `smoothed_normals` to `SurfaceMesh`**

```python
    def smoothed_normals(self, passes=TAUBIN_PASSES, lam=TAUBIN_LAMBDA, mu=TAUBIN_MU, clamp=TAUBIN_CLAMP):
        """Patch normals of a Taubin-smoothed copy of this surface. Nothing is moved: the copy is thrown away.

        Element death leaves a staircase whose steps are one cell tall and one cell wide, so its normals deviate
        from the true surface by up to 45 degrees -- and theta, the angle between the normal and the flight
        direction, is what the Lees distribution and the film tangents read. Taubin's pair (a shrinking pass at
        +lam, an inflating one at mu slightly more negative) removes the cell-scale roughness without the net
        shrinkage a plain Laplacian would cause. It looks only at a node and its neighbours, so it assumes nothing
        about the shape: any genus, any number of components, concave or convex. Each displacement is clamped to
        `clamp` of the mean patch edge, so a node can never pass its neighbours."""
        import math
        import scipy.sparse as sp                    # imported here, like gmsh and meshio
        node = np.unique(self.faces)
        local = np.full(int(self.faces.max()) + 1, -1, dtype=np.int64)
        local[node] = np.arange(node.size)
        pairs = np.sort(local[self.faces][:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
        pairs = np.unique(pairs, axis=0)
        pairs = pairs[pairs[:, 0] != pairs[:, 1]]                       # a collapsed patch has a self-edge
        w = sp.coo_matrix((np.ones(2 * len(pairs)), (np.r_[pairs[:, 0], pairs[:, 1]],
                                                     np.r_[pairs[:, 1], pairs[:, 0]])),
                          shape=(node.size, node.size)).tocsr()
        degree = np.maximum(np.asarray(w.sum(axis=1)).ravel(), 1.0)     # an isolated node stays where it is
        p = self.points[node].copy()
        step = clamp * math.sqrt(4.0 * self.areas.mean() / math.sqrt(3.0))    # the mean patch's own edge
        for _ in range(passes):
            for weight in (lam, mu):
                delta = weight * (w @ p / degree[:, None] - p)
                length = np.linalg.norm(delta, axis=1)
                over = length > step
                delta[over] *= (step / length[over])[:, None]
                p = p + delta
        moved = self.points.copy()
        moved[node] = p
        a, b, c = (moved[self.faces[:, i]] for i in range(3))
        n = np.cross(b - a, c - a)
        length = np.linalg.norm(n, axis=1)
        flat = length < 1e-30                                           # a patch the recession collapsed
        n = np.where(flat[:, None], self.normals, n / np.where(flat, 1.0, length)[:, None])
        n[np.einsum("ij,ij->i", n, self.normals) < 0.0] *= -1.0         # keep this surface's own orientation
        return n
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: PASS.

- [ ] **Step 6: Run the whole unit tier — nothing else may regress**

Run: `cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/proto3" && "$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py`
Expected: the counts `prototype/README.md` records — 201 passed, with 2 failures and 3 errors that all come from the two melting SESAM reference runs being absent from `proto3/data/reference_runs/`. **Any other failure is a regression and must be fixed before committing.**

- [ ] **Step 7: Commit**

```bash
git add prototype/proto3/reentry_model/mesh.py prototype/proto3/tests/test_reentry_model_mesh.py
git commit -m "mesh: smoothed normals, so the heating stops reading a 45-degree staircase

A death step is one cell tall and one cell wide. Taubin's lambda/mu pair takes
the cell-scale roughness out of the normals the Lees distribution and the film
tangents read, without moving the mesh and without assuming a shape.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Acceptance for the whole plan

1. `"$PY" -m pytest tests/test_reentry_model_mesh.py -q` green in `prototype/proto3`, with the five pre-existing Step 3 mesh tests unmodified.
2. The unit tier matches `prototype/README.md`'s recorded counts; no new failure.
3. `VolumeMesh.points` is bit-identical after any sequence of writes to `derived_points`.
4. A 100 mm sphere at `h_surface` 2.16 mm, `h_core` 8 mm, `band` 15 mm, `layers=0` has **310 000 ± 15 %** tetrahedra. Measure it once and record the number — it is the figure spec §6 predicts and Task 14 will re-measure.
5. `smoothed_normals` at least halves the mean normal error against an analytic sphere under cell-scale noise.

## Not in this plan

Task 16 (`recession.py`, which writes `derived_points`) and Task 17 (`remesh.py`). Any change under `thermal/` — Design B means the solver's mesh never moves, so no backend change belongs here. Any change to `reentry_model/` at the repo root. Regenerating the master plan, which happens once the prototype's tests are green.
