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


# -- the dual faces as vectors (sub-plan 21 Task 6a) --------------------------------------------------------------------
def test_the_dual_faces_close_every_column_so_a_uniform_flow_neither_gathers_nor_spreads_metal(coarse_sphere_mesh):
    """`face_vector`: each dual face's two segments, each its normal in its triangle times its length, summed and
    pointing from the edge's first corner to its second -- what a velocity crossing the face is dotted with. Around
    every interior corner of a jittered plane the faces close to round-off, so a uniform flow carries as much into each
    column as out of it; the edge's direction times the face's length does not close where the triangles are not
    symmetric about their edges, and on this plane a uniform 1 m/s would gather or spread a column's metal at up to
    539 1/s (measured 2026-10-21; 2.9e-13 with the faces' vectors). On the test sphere the faces close in each column's
    tangent plane to 3.7 % of the mean face length at worst and 0.03 % on average; what is left is along the normal,
    the surface's curvature (24 % of a face length on average)."""
    pts, faces, inner = flat_grid()
    col = columns.Columns(pts, faces)
    i, j = col.edges[:, 0], col.edges[:, 1]
    S = col.face_vector
    assert np.all(np.einsum("ij,ij->i", S, pts[j] - pts[i]) > 0.0)                    # from the first corner out
    assert np.all(np.linalg.norm(S, axis=1) <= col.face_length * (1.0 + 1e-12))
    u = np.array([1.0, 0.3, 0.0])
    assert np.abs((col.incidence @ (S @ u))[inner] / col.area[inner]).max() < 1e-9
    e = pts[j] - pts[i]
    naive = e / np.linalg.norm(e, axis=1)[:, None] * col.face_length[:, None]
    assert np.abs((col.incidence @ (naive @ u))[inner] / col.area[inner]).max() > 100.0
    col, _ = sphere_columns(coarse_sphere_mesh)
    closure = col.incidence @ col.face_vector
    tangential = closure - np.einsum("ij,ij->i", closure, col.normal)[:, None] * col.normal
    assert np.linalg.norm(tangential, axis=1).max() < 0.075 * col.face_length.mean()
    assert np.linalg.norm(tangential, axis=1).mean() < 5e-4 * col.face_length.mean()


# -- the thin film's curvature (sub-plan 21 Task 6b) --------------------------------------------------------------------
def test_the_film_curvature_is_the_base_s_plus_the_top_s_cotangent_laplacian_of_the_heights(coarse_sphere_mesh):
    """On a plane a cosine's curvature is the compact second difference of the heights (to 4e-7, the cotangents of its
    own tilted top), and on the test sphere a uniform layer's is the base's, Task 2's, exactly. The cotangents are the
    top surface's, but a cotangent depends on angles alone: on a sphere a uniform layer's top has the base's angles but
    for the corners' normals not being quite radial (the two operators on a smooth bump differ by 0.35 % on average at
    3.4 mm, measured 2026-10-21), and the area is the column's, the base's, which the energy argument requires."""
    pts, faces, inner = flat_grid(n=24, jitter=0.0)
    col = columns.Columns(pts, faces)
    H = 1e-6 * np.cos(2.0 * np.pi * pts[:, 0] / 8e-3)
    k2 = (2.0 * np.sin(np.pi * 1e-3 / 8e-3) / 1e-3) ** 2
    np.testing.assert_allclose(col.film_curvature(H)[inner], k2 * H[inner], rtol=0, atol=1e-6 * k2 * 1e-6)
    col, _ = sphere_columns(coarse_sphere_mesh)
    uniform = np.full(col.n, 3.4e-3)
    np.testing.assert_allclose(col.film_curvature(uniform), col.curvature(), rtol=0, atol=1e-9)
    bump = 1e-6 * np.cos(col.points[:, 0] / 0.02)
    top, base = col.film_operator(uniform) @ bump, col.film_operator(np.zeros(col.n)) @ bump
    big = np.abs(base) > 0.5 * np.abs(base).max()
    assert np.mean(top[big] / base[big]) == pytest.approx(1.0, abs=0.01)


def test_a_pile_up_lowers_its_own_curvature_at_an_obtuse_corner_of_the_top_but_raises_the_film_s(coarse_sphere_mesh):
    """Why surface tension takes the thin film's curvature (Task 6b; Asha, 2026-10-21): a corner touching an obtuse
    triangle takes its neighbours' mean curvature (Task 2), so raising it makes its neighbours concave and lowers its
    own -- a pile-up there would draw metal in, surface tension running backwards: raising each of the test sphere's 33
    such corners by 1 um lowers its own curvature at 32 and leaves it unchanged at the other, by up to 0.083 1/m, where
    every other corner's rises by at least 0.21 1/m. The film's curvature rises at every one of them, by 0.19-0.75 1/m
    (measured 2026-10-21)."""
    col, _ = sphere_columns(coarse_sphere_mesh)
    H = np.full(col.n, 3e-4)
    _, _, cot = columns.triangle_geometry(col.top_points(H), col.faces)
    obtuse = np.zeros(col.n, dtype=bool)
    for v in range(3):
        np.logical_or.at(obtuse, col.faces[:, v], (cot < 0.0).any(axis=1))
    corners = np.flatnonzero(obtuse)
    assert len(corners) > 20
    own_top, own_film = [], []
    for c in corners:
        bumped = H.copy()
        bumped[c] += 1e-6
        own_top.append(col.curvature(col.top_points(bumped))[c] - col.curvature(col.top_points(H))[c])
        own_film.append(col.film_curvature(bumped)[c] - col.film_curvature(H)[c])
    own_top, own_film = np.array(own_top), np.array(own_film)
    assert np.all(own_top <= 0.0) and (own_top < 0.0).sum() >= 30 and np.all(own_film > 0.0)


def test_the_faces_exchange_times_the_film_operator_has_no_growing_mode(coarse_sphere_mesh):
    """The film's curvature is (1 / A_i) dE / dH_i for the thin film's surface energy E, A_i the column's area that turns
    its volume into its height (Asha, 2026-10-21): the faces' exchange B diag(w) B^T / A times the film operator is a
    product of two symmetric non-negative matrices once scaled by the areas, so no pattern of heights grows by itself --
    on the test sphere with a rough top (0-3.4 mm) every eigenvalue is real and at least -1e-14 of the largest (measured
    3.3e-5 against 7.7e11 for a uniform top, 2026-10-21; Task 2's curvature of the top gave -1.5e10, a growth of about
    80 1/s for 0.3 mm of liquid)."""
    import scipy.linalg as la
    col, _ = sphere_columns(coarse_sphere_mesh)
    H = np.random.default_rng(3).uniform(0.0, 3.4e-3, col.n)
    B = col.incidence.toarray()
    K = (B * col.conduct_weight[None, :]) @ B.T / col.area[:, None]
    ev = la.eigvals(K @ col.film_operator(H).toarray())
    assert np.abs(ev.imag).max() <= 1e-9 * ev.real.max() and ev.real.min() > -1e-14 * ev.real.max()
