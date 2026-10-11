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


# -- the backends on moved geometry (Task 2) -----------------------------------------------------------------------
def test_moved_geometry_gives_the_operators_of_a_fresh_solver(coarse_sphere_mesh):
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    mat = constant_material()
    solver = solver_on(m, mat, 450.0)
    r = np.linalg.norm(m.points, axis=1)
    moved = r > 0.04                                          # only the outer shell (r above 4 cm) moves; the core stays put
    m.points[moved] *= (1.0 - 0.04 * (r[moved] / 0.05) ** 4)[:, None]    # an uneven squeeze toward the centre, in place
    outer = np.flatnonzero(moved[m.tets].any(axis=1))         # the elements with a moved node: the only ones that change
    assert outer.size < m.n_elements                          # a real subset: the core elements keep their matrices
    solver.update_geometry(outer)                             # recompute just those elements, and the surface
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
