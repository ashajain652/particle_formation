"""The corner columns (layered-skin design section 5.1, sub-plan 20 Task 1): dual cells, normals, the edge graph and its
cotangent conduction weights, and the corner surface the coupled run reads."""
import numpy as np
import pytest

from reentry_model import columns, movemesh


def sphere_columns(m):
    s = m.surface()
    corners = np.unique(s.faces)
    local = np.full(m.n_nodes, -1)
    local[corners] = np.arange(len(corners))
    faces = local[movemesh.orient_outward(m.points, s.faces, s.normals)]
    return columns.Columns(m.points[corners], faces), s


def flat_grid(n=8, h=1e-3, jitter=0.25, seed=2):
    """A jittered flat triangulation of the square [0, n h]^2 in z = 0 (open edges), normals +z."""
    rng = np.random.default_rng(seed)
    x, y = np.meshgrid(np.arange(n + 1) * h, np.arange(n + 1) * h, indexing="ij")
    pts = np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])
    inner = (x.ravel() > 0) & (x.ravel() < n * h) & (y.ravel() > 0) & (y.ravel() < n * h)
    pts[inner, :2] += rng.uniform(-jitter, jitter, (inner.sum(), 2)) * h
    idx = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
    faces = np.vstack([np.column_stack([a, b, c]), np.column_stack([a, c, d])])
    return pts, faces, inner


def test_dual_areas_add_up_to_the_surface_and_normals_point_out(coarse_sphere_mesh):
    col, s = sphere_columns(coarse_sphere_mesh)
    assert col.n == len(np.unique(s.faces))
    assert col.area.sum() == pytest.approx(s.areas.sum(), rel=1e-12)
    radial = col.points / np.linalg.norm(col.points, axis=1)[:, None]
    # the area-weighted normal of a 4 mm triangulation of a 50 mm sphere follows the radius: 0.18 degrees off on
    # average and 4.7 at worst, where gmsh left an irregular corner (measured while planning)
    angle = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", col.normal, radial), -1, 1)))
    assert angle.mean() < 0.5 and angle.max() < 6.0
    assert np.all(col.degree >= 3)
    assert len(col.edges) == 3 * len(col.faces) // 2           # closed: every edge in exactly two triangles


def test_cotangent_weights_are_exact_for_a_linear_field_on_a_flat_surface():
    pts, faces, inner = flat_grid()
    col = columns.Columns(pts, faces)
    f = 2.0 + 300.0 * pts[:, 0] - 700.0 * pts[:, 1]
    assert np.abs(col.laplacian(f)[inner]).max() < 1e-12
    assert np.all(col.normal[:, 2] > 0.999999)


def test_on_an_equilateral_grid_the_weight_is_the_dual_face_over_the_edge():
    h = 1e-3
    rows = 6
    pts = np.array([[i * h + 0.5 * h * (j % 2), j * h * np.sqrt(3) / 2, 0.0] for j in range(rows) for i in range(rows)])
    faces = []
    for j in range(rows - 1):
        for i in range(rows - 1):
            a, b, c, d = j * rows + i, j * rows + i + 1, (j + 1) * rows + i, (j + 1) * rows + i + 1
            faces += [[a, b, c], [b, d, c]] if j % 2 == 0 else [[a, b, d], [a, d, c]]
    col = columns.Columns(pts, np.array(faces))
    interior = np.all(col._tri >= 0, axis=1)
    assert np.allclose(col.edge_length[interior], h)
    assert np.allclose(col.cot_weight[interior], 1.0 / np.sqrt(3.0))
    assert np.allclose(col.face_length[interior] / col.edge_length[interior], col.cot_weight[interior])


def test_the_cotangent_laplacian_of_a_sphere_sees_its_curvature(coarse_sphere_mesh):
    col, _ = sphere_columns(coarse_sphere_mesh)
    # sum_j w_ij (x_j - x_i) = -2 A_i n_i / R on a sphere of radius R; with barycentric areas it is right on average
    # (0.5 % high, with a scatter of 7.6 % from corner to corner, measured while planning)
    lap = col.laplacian(col.points)
    curv = -np.einsum("ij,ij->i", lap, col.normal) / (2.0 * col.area)
    assert np.mean(curv) == pytest.approx(1.0 / 0.05, rel=0.02)


def test_update_follows_moved_corners_and_the_top_surface_sits_above_the_base(coarse_sphere_mesh):
    col, _ = sphere_columns(coarse_sphere_mesh)
    area0 = col.area.sum()
    col.update(col.points * 0.9)                                  # a sphere 10 % smaller
    assert col.area.sum() == pytest.approx(0.81 * area0, rel=1e-12)
    top = columns.CornerSurface.from_points(col.top_points(np.full(col.n, 1e-3)), col.faces)
    r = np.linalg.norm(top.points, axis=1)
    lift = r - np.linalg.norm(col.points, axis=1)                 # 1 mm times the cosine of the normal's tilt
    assert np.all(lift <= 1e-3 + 1e-12) and np.all(lift > 0.99e-3)
    assert top.n_patches == col.n and top.area > col.area.sum()
    v = np.array([-1.0, 0.0, 0.0])
    assert top.angles_to(v)[top.patch_toward(v)] < np.radians(5.0)
    assert top.projected_area(v) == pytest.approx(np.pi * (0.045 + 1e-3) ** 2, rel=0.01)


# -- the surface operators (sub-plan 21 Task 2) -------------------------------------------------------------------------
def obtuse_corners(col, points=None):
    p = col.points if points is None else points
    _, _, cot = columns.triangle_geometry(p, col.faces)
    out = np.zeros(col.n, dtype=bool)
    for v in range(3):
        np.logical_or.at(out, col.faces[:, v], (cot < 0.0).any(axis=1))
    return out


def test_the_gradient_is_exact_for_a_field_linear_on_a_flat_surface():
    pts, faces, inner = flat_grid()
    col = columns.Columns(pts, faces)
    g = col.gradient(2.0 + 300.0 * pts[:, 0] - 700.0 * pts[:, 1])
    assert np.abs(g - np.array([300.0, -700.0, 0.0])).max() < 1e-9          # every corner, the edges' included


def test_on_a_sphere_the_gradient_of_a_linear_field_is_its_tangential_part(coarse_sphere_mesh):
    """On the 4 mm test sphere within 0.09 % of a linear field's tangential part on average and 0.8 % at worst, at the
    irregular corner gmsh left (0.04 % at worst on the 2 mm sphere; measured 2026-10-18)."""
    col, _ = sphere_columns(coarse_sphere_mesh)
    a = np.array([3.0, -1.0, 2.0])
    exact = a - (col.normal @ a)[:, None] * col.normal
    err = np.linalg.norm(col.gradient(col.points @ a) - exact, axis=1) / np.linalg.norm(a)
    assert err.mean() < 0.002 and err.max() < 0.01


def test_the_curvature_of_a_sphere_is_two_over_its_radius_at_every_corner(coarse_sphere_mesh):
    """Every corner within 2e-4 of 2/R, those touching an obtuse triangle included (33 of the 2 254, which the mixed
    Voronoi area alone left off by up to 176 %; 8.4e-5 measured at worst, 2026-10-18), and the same on a sphere 2 %
    larger whose corners are given."""
    col, _ = sphere_columns(coarse_sphere_mesh)
    assert obtuse_corners(col).sum() > 0                                     # the fallback is exercised
    R = np.linalg.norm(col.points, axis=1)
    assert np.abs(col.curvature() * R / 2.0 - 1.0).max() < 2e-4
    assert np.abs(col.curvature(1.02 * col.points) * 1.02 * R / 2.0 - 1.0).max() < 2e-4


def test_the_curvature_of_a_plane_is_zero():
    pts, faces, inner = flat_grid()
    assert np.abs(columns.Columns(pts, faces).curvature()[inner]).max() < 1e-9
