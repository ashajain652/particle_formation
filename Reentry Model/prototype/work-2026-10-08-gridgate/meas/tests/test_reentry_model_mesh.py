"""mesh.py: gmsh sphere generation and cache, boundary faces, outward normals, patch angles."""
import math
import os

import numpy as np
import pytest

from reentry_model import mesh

R = 0.05


def test_coarse_sphere_geometry(coarse_sphere_mesh):
    m = coarse_sphere_mesh
    assert 2000 < m.n_nodes < 5000 and m.n_elements > 4 * m.n_nodes // 2
    assert m.volume() == pytest.approx(4.0 / 3.0 * math.pi * R ** 3, rel=3e-3)
    s = m.surface()
    assert s.area == pytest.approx(4.0 * math.pi * R ** 2, rel=2e-3)
    assert np.allclose(np.linalg.norm(s.centroids, axis=1), R, atol=3e-4)      # centroids sag by ~h^2/(6R) below the sphere
    assert np.all(np.einsum("ij,ij->i", s.normals, s.centroids) > 0.95 * np.linalg.norm(s.centroids, axis=1))   # outward, within 18 deg of radial
    assert np.allclose(np.linalg.norm(s.normals, axis=1), 1.0)
    assert np.linalg.norm(m.points[m.centre_node()]) < 1e-9
    assert m.boundary_nodes().size == np.unique(s.faces).size and m.params["radius_m"] == R


def test_angles_and_patch_lookup(coarse_sphere_mesh):
    s = coarse_sphere_mesh.surface()
    theta = s.angles_to([1.0, 0.0, 0.0])
    assert theta.min() < math.radians(5.0) and theta.max() > math.radians(175.0)
    windward = theta < math.pi / 2
    assert s.areas[windward].sum() == pytest.approx(0.5 * s.area, rel=2e-2)
    assert s.normals[s.patch_toward([1.0, 0.0, 0.0])][0] > 0.99
    assert s.normals[s.patch_toward([-1.0, 0.0, 0.0])][0] < -0.99
    x = coarse_sphere_mesh.points[:, 0]
    assert np.allclose(s.facet_mean(x), s.centroids[:, 0])


def test_boundary_faces_are_single_use_faces():
    tets = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
    faces, opposite = mesh.boundary_faces(tets)
    assert faces.shape == (6, 3) and not any(set(f) == {1, 2, 3} for f in faces)
    assert opposite.shape == (6,)


def test_generation_is_cached(tmp_path):
    p1 = mesh.generate_sphere_mesh(0.01, 2e-3, 4e-3, str(tmp_path))
    mtime = os.path.getmtime(p1)
    p2 = mesh.generate_sphere_mesh(0.01, 2e-3, 4e-3, str(tmp_path))
    assert p1 == p2 and os.path.getmtime(p2) == mtime
    assert os.path.basename(p1) == "sphere_R10.000mm_hs2.000mm_hc4.000mm_b15.000mm_r8.000mm.msh"
    m = mesh.load_mesh(p1)
    assert m.volume() == pytest.approx(4.0 / 3.0 * math.pi * 0.01 ** 3, rel=3e-2)        # crude 2 mm elements on a 10 mm sphere


def test_points_are_mutable_and_surface_follows(coarse_sphere_mesh):
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets)
    area0 = m.surface().area
    m.points *= 0.5
    assert m.surface().area == pytest.approx(0.25 * area0) and m.volume() == pytest.approx(0.125 * coarse_sphere_mesh.volume())


# ---------------------------------------------------------------------------------------------------------------
# Step 3: prism layers, the active set, the box mesh and the surface helpers

@pytest.fixture(scope="module")
def layered_mesh(tmp_path_factory):
    """100 mm sphere, 4 mm / 20 mm inner mesh with two 0.25 mm / 0.5 mm prism layers (~5 k nodes)."""
    return mesh.sphere_mesh(R, 4e-3, 20e-3, str(tmp_path_factory.mktemp("meshes")), layers=2)


def test_layer_thicknesses_and_defaults():
    assert np.allclose(mesh.layer_thicknesses(4, 0.25e-3, 2.0), [0.25e-3, 0.5e-3, 1e-3, 2e-3])
    assert mesh.DEFAULT_LAYERS == 0 and mesh.DEFAULT_LAYER_THICKNESS == 0.25e-3 and mesh.DEFAULT_LAYER_GROWTH == 2.0  # the dense band replaces the prism stack as the default
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
    assert 0.5 * R < d.max() < R      # the deepest element is well inside, but no depth can exceed the radius
    # not approx(R): the depth is measured at element CENTROIDS, and a 20 mm core cell puts the innermost
    # centroid 14 mm from the origin, so the deepest depth this fixture can show is 0.72 R (measured 2026-09-27)
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
    """Review focus 1: the pinned Step 2 verification devices depend on band = 0 being the OLD field exactly --
    DistMin 0, DistMax radius -- not merely 'graded'. Letting DistMax collapse to the ramp still grades, but to
    h_core within 8 mm, which cost 72 % of the elements on the reference geometry (measured 2026-09-27)."""
    assert mesh.band_limits(0.0, 8e-3, 0.05) == (0.0, 0.05)       # the whole radius, not the ramp
    m = mesh.sphere_mesh(0.02, 2e-3, 8e-3, str(tmp_path), band=0.0)
    d = m.element_depths()
    edge = (m.element_volumes() / 0.117851) ** (1.0 / 3.0)
    near, mid = edge[(d >= 0) & (d < 2e-3)], edge[(d > 3e-3) & (d < 5e-3)]
    assert near.size and mid.size
    assert mid.mean() < 0.75 * 8e-3          # a ramp-length DistMax would already be at h_core here
    assert m.params["band_m"] == 0.0


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


def test_params_record_the_band_the_mesh_actually_got(tmp_path):
    """A body smaller than the band clamps it, and the remesh trigger reads band_m to decide when half the band is
    gone -- so it must be the effective value, not the requested one."""
    clamped = mesh.sphere_mesh(0.01, 2e-3, 8e-3, str(tmp_path), band=0.05, ramp=8e-3)
    assert clamped.params["band_m"] == pytest.approx(0.01)      # clamped to the radius, not the 50 mm requested
    plain = mesh.sphere_mesh(0.01, 2e-3, 8e-3, str(tmp_path), band=4e-3, ramp=8e-3)
    assert plain.params["band_m"] == pytest.approx(4e-3)        # unclamped: recorded as asked


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
    m.derived_points[bn] *= 0.998                                   # 0.1 mm inward, well inside the 0.25 mm outermost layer: a displacement
                                                                    # larger than the local element inverts it and the derived surface is then
                                                                    # meaningless (measured 2026-09-27: 0.995, one full layer, inverts 7956)

    def orientation(pts):
        x = pts[m.tets]
        return np.sign(np.linalg.det(np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)))

    assert (orientation(m.derived_points) == orientation(m.points)).all()      # the displacement inverts nothing

    assert np.array_equal(m.points, p0)                             # bit-identical: the solver sees nothing
    assert m.volume() == v0 and np.array_equal(m.element_volumes(), e0)
    assert m.surface().area == pytest.approx(a0)                    # the default surface is still the real one
    d = m.surface(derived=True)
    assert d.area == pytest.approx(0.998 ** 2 * a0, rel=1e-6)       # and the derived one shrank by exactly r^2
    assert d.points is m.derived_points
    assert np.array_equal(d.faces, m.surface().faces)               # same patches, different coordinates


def test_derived_surface_normals_stay_outward(layered_mesh):
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    m.derived_points[m.boundary_nodes()] *= 0.998
    d = m.surface(derived=True)
    outward = np.einsum("ij,ij->i", d.normals, d.centroids / np.linalg.norm(d.centroids, axis=1)[:, None])
    assert (outward > 0.5).all()                                    # still pointing out of the shrunken sphere


def test_smoothed_normals_halve_the_error_a_staircase_introduces(layered_mesh):
    """Amendment 2026-09-27: the noise amplitude here is tuned down from the design's 0.3 mm. On this fixture's
    thinnest prism layer (0.25 mm), a per-node 3-axis Gaussian of standard deviation 0.09 mm is the most this seed
    tolerates without inverting an element (measured: 0.10 mm already inverts 18); much of the 0.25 mm budget is
    eaten by the random direction, unlike the coherent radial scaling the derived-surface tests above use. The
    inversion check below is the same determinant-sign test Task 3 uses, kept here because a silently inverted
    element would make both the rough and the smoothed normal meaningless rather than just noisy."""
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    exact = m.surface()                                              # the analytic sphere's own normals
    rng = np.random.default_rng(0)
    bn = m.boundary_nodes()
    m.derived_points[bn] += rng.normal(0.0, 9.0e-5, (bn.size, 3))    # cell-scale roughness, like element death

    def orientation(pts):
        x = pts[m.tets]
        return np.sign(np.linalg.det(np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)))

    assert (orientation(m.derived_points) == orientation(m.points)).all()      # the noise inverts nothing
    rough = m.surface(derived=True)
    ang = lambda n: np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", n, exact.normals), -1.0, 1.0)))
    assert ang(rough.normals).mean() > 2.0                           # the staircase really is that bad (measured 2.33 deg)
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


def test_a_collapsed_patch_never_produces_a_non_finite_normal(layered_mesh):
    """surface().normals feeds theta, the film tangents and the drag area, so one NaN reaches the physics with
    nothing but a RuntimeWarning to announce it (measured 2026-09-27: 2 of 4598 patches went non-finite)."""
    m = mesh.VolumeMesh(layered_mesh.points.copy(), layered_mesh.tets)
    f = m.surface().faces[0]
    m.derived_points[f[1]] = m.derived_points[f[0]]          # two nodes coincident: a zero-area patch
    s = m.surface(derived=True)
    assert np.isfinite(s.normals).all()
    assert np.allclose(np.linalg.norm(s.normals, axis=1), 1.0)
    assert np.isfinite(s.smoothed_normals()).all()
