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
