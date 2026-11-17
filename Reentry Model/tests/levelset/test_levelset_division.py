"""Tests of the division of the layer between the gas surface phi and the skin base beta among the markers (design
sections 21.5 and 21.11; plan 03b Task 1, planning ruling 1): each column's cross-section and each pair's shared boundary
at every depth, the same-side rule, metal thinner than twice the skin depth, and each piece's closure factor.

Every phi here is a body's exact signed distance sampled at the nodes of a box mesh refined through the layer, so its
level surfaces down to the skin depth are well defined."""
import functools
import math
import time

import numpy as np

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import inside_volumes, patches
from reentry_model.levelset.division import DIVISION_LEVELS, DIVISION_TOL, divide
from reentry_model.levelset.levelset import pieces
from reentry_model.levelset.markers import OPPOSED_NORMALS, OWN_SHARE, ownership

from cases import DESIGN_SPACING, box_sdf, polygon_sdf, seeded

H = 1.2e-3                 # m: the skin depth (design section 21.3)
OFF = 0.13e-3              # m: the bodies' offset off the grid (the stage 3 box's)
GRID_OFF = 0.07e-3         # m: the marker grids' offset off the mesh's node planes
FINE = 0.35e-3             # m: the cells through the layer of the curved and edged bodies
BOX_LENGTH, BOX_HALF, BOX_H0 = 8e-3, 5e-3, 2.8e-3     # m: the stage 3 box (analysis/levelset/stage3_cases.py)


def _layer_mesh(base, sdf, fine=FINE):
    """`base` refined to `fine` [m] through the layer (-H < sdf < 0) and two fine cells either side of it."""
    def target(x, s):
        f = sdf(x)
        pad = 2.0 * fine + 0.87 * s
        return np.where((f < pad) & (f > -H - pad), fine, np.inf)
    return build(base, target)


def _slabs_sdf(*spans):
    """The exact signed distance to the union of disjoint slabs lo < y < hi [m] (the minimum of theirs)."""
    def sdf(p):
        y = p[:, 1]
        return np.min([np.maximum(lo - y, y - hi) for lo, hi in spans], axis=0)
    return sdf


def _sheet(y, up, lo, hi, shift=0.0):
    """Markers on the plane y [m] on a grid DESIGN_SPACING apart in x and z from lo + half a spacing (+ shift + GRID_OFF)
    to hi, with the normal +y (`up`) or -y. Returns (positions (n, 3), normals (n, 3))."""
    g = np.arange(lo + 0.5 * DESIGN_SPACING + shift, hi, DESIGN_SPACING) + GRID_OFF
    X, Z = np.meshgrid(g, g, indexing="ij")
    x = np.column_stack([X.ravel(), np.full(X.size, y), Z.ravel()])
    return x, np.tile([0.0, 1.0 if up else -1.0, 0.0], (len(x), 1))


def _pair(d, a, b):
    """The index of the boundary between markers a and b in the division d, or None."""
    i, j = min(a, b), max(a, b)
    e = np.flatnonzero((d.i == i) & (d.j == j))
    return int(e[0]) if len(e) else None


def _seeded_box(lo, hi, base):
    """A box's exact signed distance on `base` refined through its layer, and markers seeded on it as the body seeds
    them (one per DESIGN_SPACING voxel)."""
    sdf = box_sdf(lo, hi)
    mesh = _layer_mesh(base, sdf)
    phi = sdf(mesh.points)
    return mesh, phi, seeded(mesh, phi, DESIGN_SPACING, base.lo)


def test_the_columns_sum_to_the_layer_to_round_off():
    """A 6 mm cube (0.13 mm off the grid, 0.35 mm cells through the layer, markers seeded as the body seeds them), whose
    layer is mostly full cells of one size: the columns' summed volume, `d.volume.sum()`, equals the layer's volume summed
    exactly (`math.fsum`) to within 1e-14 of it, a hundredth of the stage 3 cases' mass tolerance. The closure's sum of
    the layer must not be a running sum: on many equal cells that rounds the same way at every addition, and the
    columns carried its error (measured once on this fixture, 1.0e-13 of the layer too small with `np.bincount`, 4.7e-16
    with an exact sum; Task 7d, ruling T7d-a)."""
    side = np.full(3, 6e-3)
    lo = np.full(3, OFF - 3e-3)
    base = BoxBase.around((-3e-3 - BOX_H0,) * 3, (3e-3 + BOX_H0,) * 3, BOX_H0)
    mesh, phi, mk = _seeded_box(lo, lo + side, base)
    exact = math.fsum(inside_volumes(mesh, phi) - inside_volumes(mesh, phi + H))
    d = divide(mesh, phi, H, mk.x, mk.normal)
    assert abs(d.volume.sum() - exact) <= 1e-14 * exact


def test_a_flat_slab_divides_into_prisms():
    """A flat face (phi = y - OFF, the metal below it filling the mesh) with markers on a 2 mm grid: every interior
    column is a 2 mm square prism, each pair of interior neighbours shares a 2 mm boundary at every level, and the
    columns tile the layer exactly. Each interior column's centroid at each level (plan 03b Task 6, the controller's
    ruling T6-b: `Division.centroid`, each level triangle's share at the triangle's centroid) lies on that level, to
    1e-12 m, and beneath its marker, to CENTROID_SLACK."""
    base = BoxBase((-6e-3, -2.5e-3, -6e-3), 0.5e-3, (24, 6, 24))
    mesh = build(base)
    phi = mesh.points[:, 1] - OFF
    x, nrm = _sheet(OFF, True, -6e-3, 6e-3)
    d = divide(mesh, phi, H, x, nrm)
    assert d.area.shape == (len(x), DIVISION_LEVELS + 1) and d.levels[-1] == H
    span = DESIGN_SPACING * 2.0
    inner = np.all(np.abs(x[:, [0, 2]] - GRID_OFF) < span, axis=1)    # the 4 x 4 markers whose squares lie inside
    assert inner.sum() == 16
    assert np.abs(d.area[inner] / DESIGN_SPACING ** 2 - 1.0).max() <= 1e-9
    rows = np.flatnonzero(inner)
    checked = 0
    for a in rows:
        for b in rows[rows > a]:
            if abs(np.linalg.norm(x[a] - x[b]) - DESIGN_SPACING) < 1e-9:
                e = _pair(d, a, b)
                assert e is not None
                assert np.abs(d.length[e] / DESIGN_SPACING - 1.0).max() <= 1e-9
                assert abs(d.normal[e] @ (x[b] - x[a]) / DESIGN_SPACING - 1.0) <= 1e-12
                checked += 1
    assert checked == 24
    assert np.abs(d.factor - 1.0).max() <= 1e-12
    assert not d.opposed.any()
    c = d.centroid[inner]
    off_level = np.abs(c[..., 1] - (OFF - d.levels)[None, :]).max()
    sideways = np.linalg.norm(c[..., [0, 2]] - x[inner][:, None, [0, 2]], axis=2).max()
    print("flat slab: interior centroids {:.1e} m off their levels, {:.4f} mm from beneath their markers".format(
        off_level, sideways * 1e3))
    assert off_level <= 1e-12 and sideways <= CENTROID_SLACK


CENTROID_SLACK = 0.05e-3   # m: a tenth of the slab mesh's 0.5 mm cells (each triangle's share is placed at its centroid;
                           # under 5e-8 m measured once here, the squares' cut triangles placed symmetrically)


def test_the_base_level_s_ownership_divides_beta_s_patches_whole():
    """Plan 03b Task 5's fix round 1 (the controller's ruling T5-j): the division keeps its base level's ownership by the
    plain area (`Division.base_owner`), the groups of beta the body's base law uses. On the flat slab of the prism test
    the base level's patches lie in the cells phi + H cuts, each divided whole among the markers (their columns sum to
    the patches' areas to 1e-12), and each interior marker holds its 2 mm square of the base to 1e-9, which is its
    column's base A(H) (the faces' |grad phi| is one)."""
    base = BoxBase((-6e-3, -2.5e-3, -6e-3), 0.5e-3, (24, 6, 24))
    mesh = build(base)
    phi = mesh.points[:, 1] - OFF
    x, nrm = _sheet(OFF, True, -6e-3, 6e-3)
    d = divide(mesh, phi, H, x, nrm)
    P = patches(mesh, phi + H)
    assert np.array_equal(np.sort(d.base_cell), np.sort(P.cell))
    assert d.base_owner.shape == (len(x), len(d.base_cell))
    np.testing.assert_allclose(np.asarray(d.base_owner.sum(axis=0)).ravel(), d.base_patch_area, rtol=1e-12)
    inner = np.all(np.abs(x[:, [0, 2]] - GRID_OFF) < 2.0 * DESIGN_SPACING, axis=1)
    held = np.asarray(d.base_owner.sum(axis=1)).ravel()
    assert np.abs(held[inner] / DESIGN_SPACING ** 2 - 1.0).max() <= 1e-9
    np.testing.assert_allclose(held[inner], d.base[inner], rtol=1e-9)


def test_a_spherical_shell_divides_by_the_shell_volume():
    """A sphere of radius 5 mm, markers seeded 2 mm apart as the body seeds them: the columns' volumes
    close the shell to within DIVISION_TOL, and every column is a cone, its base (R - H)^2 / R^2 of its top to 2 %."""
    R = 5e-3
    c = np.array([OFF, -0.07e-3, 0.11e-3])
    base = BoxBase.around(c - R - BOX_H0, c + R + BOX_H0, BOX_H0)
    sdf = lambda p: np.linalg.norm(p - c, axis=1) - R
    mesh = _layer_mesh(base, sdf)
    phi = sdf(mesh.points)
    mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
    d = divide(mesh, phi, H, mk.x, mk.normal)
    ratio = d.base / d.area[:, 0] / ((R - H) / R) ** 2
    print("sphere: {} markers, factor - 1 = {:.2e}, base over top over the cone's from {:.4f} to {:.4f}".format(
        mk.n, d.factor[0] - 1.0, ratio.min(), ratio.max()))
    assert len(d.factor) == 1
    assert abs(d.factor[0] - 1.0) <= DIVISION_TOL
    assert np.abs(ratio - 1.0).max() <= 0.02


def test_columns_beside_a_convex_edge_narrow_and_meet_on_the_halfway_plane():
    """The stage 3 box (8 x 10 x 10 mm, 0.13 mm off the grid, markers seeded as the body seeds them) and the edge of its
    10 x 10 face (x = OFF, normal -x) with the +y face: every column of either face owning a strip of phi that touches
    the edge (away from the box's corners) narrows with depth, and shares a boundary of positive length at depth with a
    column of the other face, a pair whose normals do not oppose.

    The narrowing is judged on the column's mean cross-section over bands one cell deep, each band's below the last:
    the depth the P1 mesh resolves. Level by level, the volume per unit depth (the coarea measure, ruling T1-a) of a
    column beside the edge jumps where the level surfaces pass the mesh's planes of nodes, since the nearly flat P1
    cells along the edge hold their volume within a narrow depth there: measured once (96 levels, before Task 1's fix
    round 1 placed such cells' volumes whole at their depths), rises of 12 to 44 % of the top at single levels, one
    cell apart. (The level surfaces' own areas fall smoothly but for rises of 4e-4 of the top near the base in two
    columns; Task 1's report.)"""
    lo = np.array([OFF, OFF - BOX_HALF, OFF - BOX_HALF])
    hi = lo + np.array([BOX_LENGTH, 2.0 * BOX_HALF, 2.0 * BOX_HALF])
    base = BoxBase.around((-BOX_H0, -BOX_HALF - BOX_H0, -BOX_HALF - BOX_H0),
                          (BOX_LENGTH + BOX_H0, BOX_HALF + BOX_H0, BOX_HALF + BOX_H0), BOX_H0)
    mesh, phi, mk = _seeded_box(lo, hi, base)
    t0 = time.perf_counter()
    d = divide(mesh, phi, H, mk.x, mk.normal)
    took = time.perf_counter() - t0
    largest = np.argmax(np.abs(mk.normal), axis=1)        # a marker's face: its normal's largest component (beside
    on_face = (largest == 0) & (mk.normal[:, 0] < 0.0)   # the edge its normal blends the two faces')
    on_top = (largest == 1) & (mk.normal[:, 1] > 0.0)
    along = OFF + np.linspace(-BOX_HALF + 1e-3, BOX_HALF - 1e-3, 161)
    beside = 0.02e-3                                       # m: the strip's points this far from the edge on each face
    touching = set()
    for p, nf in ((np.column_stack([np.full(161, OFF), np.full(161, OFF + BOX_HALF - beside), along]), (-1.0, 0.0, 0.0)),
                  (np.column_stack([np.full(161, OFF + beside), np.full(161, OFF + BOX_HALF), along]), (0.0, 1.0, 0.0))):
        dist = np.linalg.norm(p[:, None, :] - mk.x[None], axis=2)
        dist[:, mk.normal @ np.array(nf) <= OPPOSED_NORMALS] = np.inf        # the same-side rule
        touching |= set(np.argmin(dist, axis=1).tolist())
    touching = [i for i in sorted(touching) if on_face[i] or on_top[i]]
    print("box: {} markers on {} cells, divided in {:.1f} s, factor - 1 = {:.2e}; {} columns touch the edge".format(
        mk.n, len(mesh.tets), took, d.factor[0] - 1.0, len(touching)))
    assert len(touching) >= 6
    B = len(d.levels) - 1
    bands = np.unique(np.r_[np.arange(0, B, int(round(FINE / (H / B)))), B])   # one cell deep each (the last less)
    for i in touching:
        a = d.area[i]
        mean = np.array([np.trapezoid(a[p:q + 1], d.levels[p:q + 1]) / (d.levels[q] - d.levels[p])
                         for p, q in zip(bands[:-1], bands[1:])])
        other = on_top if on_face[i] else on_face
        across = [e for e in np.flatnonzero((d.i == i) | (d.j == i))
                  if other[d.j[e] if d.i[e] == i else d.i[e]] and d.length[e, 1:].max() > 0.0]
        print("  column {} ({}): bands' mean cross-sections {} mm2, {} boundaries across, deepest length {:.2e} m".format(
            i, "face" if on_face[i] else "top", np.round(mean * 1e6, 3).tolist(), len(across),
            max([d.length[e, -1] for e in across], default=0.0)))
        assert np.all(np.diff(mean) < 0.0) and a[-1] < a[0]
        assert across
        assert not d.opposed[across].any()


def test_a_concave_corner_column_widens_with_depth():
    """An L-shaped prism (an 8 mm square less its 4 mm quarter, through the mesh along z): the column on the concave
    edge owns, at depth, the level surface wrapped round the edge, and widens."""
    corners = OFF + np.array([[0.0, 0.0], [8.0, 0.0], [8.0, 4.0], [4.0, 4.0], [4.0, 8.0], [0.0, 8.0]]) * 1e-3
    sdf = polygon_sdf(corners)
    base = BoxBase.around((-BOX_H0, -BOX_H0, -4e-3), (8e-3 + BOX_H0, 8e-3 + BOX_H0, 4e-3), BOX_H0)
    mesh = _layer_mesh(base, sdf)
    phi = sdf(mesh.points)
    mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
    d = divide(mesh, phi, H, mk.x, mk.normal)
    i = int(np.argmin(np.linalg.norm(mk.x - np.array([4e-3 + OFF, 4e-3 + OFF, 0.0]), axis=1)))
    print("L: the concave edge's column {}: top {:.3e} m2, base {:.3e} m2; factor - 1 = {:.2e}".format(
        i, d.area[i, 0], d.base[i], d.factor[0] - 1.0))
    assert d.area[i, -1] > d.area[i, 0]


def test_ownership_never_crosses_a_gap():
    """Two sheets of metal 1.0 mm apart, markers 2 mm apart on each face, the grids on the facing faces a spacing's half
    apart, so that much of each facing face is nearer a marker across the gap than any on its own sheet: no triangle of
    one sheet is owned by a marker of the other, and no shared boundary joins them."""
    spans = ((OFF - 4.0e-3, OFF - 0.5e-3), (OFF + 0.5e-3, OFF + 4.0e-3))
    sdf = _slabs_sdf(*spans)
    mesh = _layer_mesh(BoxBase((-3e-3, -5e-3, -3e-3), 1e-3, (6, 10, 6)), sdf, fine=0.25e-3)
    phi = sdf(mesh.points)
    parts = [_sheet(spans[0][0], False, -3e-3, 3e-3), _sheet(spans[0][1], True, -3e-3, 3e-3),
             _sheet(spans[1][0], False, -3e-3, 3e-3, shift=-0.5 * DESIGN_SPACING),
             _sheet(spans[1][1], True, -3e-3, 3e-3)]
    x = np.concatenate([p[0] for p in parts])
    nrm = np.concatenate([p[1] for p in parts])
    sheet = x[:, 1] > OFF
    P = patches(mesh, phi)
    Q = len(P.tris)
    tri_sheet = P.tris.mean(axis=1)[:, 1] > OFF
    size = OWN_SHARE * DESIGN_SPACING
    across = lambda W: np.any(sheet[W.row] != tri_sheet[W.col])
    assert across(ownership(x, P.tris, np.arange(Q), Q, size).tocoo())          # nearest alone crosses the gap
    W = ownership(x, P.tris, np.arange(Q), Q, size, normals=nrm, tri_normal=P.normal[P.tri_patch]).tocoo()
    assert not across(W)
    d = divide(mesh, phi, H, x, nrm)
    assert len(d.i) > 0 and np.all(sheet[d.i] == sheet[d.j])
    assert len(d.factor) == 2 and len(np.unique(d.piece[sheet])) == 1 and d.piece[sheet][0] != d.piece[~sheet][0]


def test_metal_thinner_than_twice_the_depth_has_opposed_neighbours():
    """A slab 2.0 mm thick, less than twice the skin depth: its two faces' columns meet face to face at its mid-plane,
    neighbours whose normals oppose (design section 21.11). A slab 3.0 mm thick keeps 0.6 mm of the volume between
    them: no opposed neighbours."""
    found = {}
    for thick in (2.0e-3, 3.0e-3):
        span = (OFF - 0.5 * thick, OFF + 0.5 * thick)
        sdf = _slabs_sdf(span)
        mesh = _layer_mesh(BoxBase((-3e-3, -3e-3, -3e-3), 1e-3, (6, 6, 6)), sdf, fine=0.25e-3)
        phi = sdf(mesh.points)
        bottom, top = _sheet(span[0], False, -3e-3, 3e-3), _sheet(span[1], True, -3e-3, 3e-3)
        d = divide(mesh, phi, H, np.concatenate([bottom[0], top[0]]), np.concatenate([bottom[1], top[1]]))
        found[thick] = d.opposed
    print("opposed neighbours: {} (2.0 mm), {} (3.0 mm)".format(found[2.0e-3].sum(), found[3.0e-3].sum()))
    assert found[2.0e-3].any()
    assert not found[3.0e-3].any()


WEDGE_FACE = 8e-3          # m: the wedge's and the square's faces, over 4H long, so that thickness is not the cause


def test_a_wedge_sharper_than_sixty_degrees_has_opposed_neighbours_and_a_right_angle_has_none():
    """Design 21.11's other case: an edge sharper than about 60 degrees. A prism through the mesh whose cross-section is
    a wedge of 45 degrees, its two faces WEDGE_FACE long (its other two corners 67.5 degrees, whose normals do not
    oppose), markers seeded as the body seeds them: columns of its two faces share boundaries whose normals oppose
    (n_i . n_j = cos 135 degrees, -0.71), all of them beside the 45 degree edge. A square prism of the same face length,
    every edge 90 degrees: none."""
    half = np.radians(22.5)
    apex = np.array([OFF, OFF])
    wedge = [apex, apex + WEDGE_FACE * np.array([np.cos(half), np.sin(half)]),
             apex + WEDGE_FACE * np.array([np.cos(half), -np.sin(half)])]
    square = [OFF + np.array(c) * WEDGE_FACE for c in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))]
    upper, lower = np.array([-np.sin(half), np.cos(half), 0.0]), np.array([-np.sin(half), -np.cos(half), 0.0])
    found = {}
    for name, corners in (("wedge", wedge), ("square", square)):
        c = np.array(corners)
        base = BoxBase.around((c[:, 0].min() - BOX_H0, c[:, 1].min() - BOX_H0, -4e-3),
                              (c[:, 0].max() + BOX_H0, c[:, 1].max() + BOX_H0, 4e-3), BOX_H0)
        sdf = polygon_sdf(corners)
        mesh = _layer_mesh(base, sdf)
        phi = sdf(mesh.points)
        mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
        d = divide(mesh, phi, H, mk.x, mk.normal)
        e = np.flatnonzero(d.opposed)
        found[name] = (mk, d, e)
        print("{}: {} markers, {} opposed pair(s), their markers {} mm from the 45 degree edge at most".format(
            name, mk.n, len(e), np.round(np.linalg.norm(mk.x[np.r_[d.i[e], d.j[e]], :2] - apex, axis=1).max() * 1e3, 2)
            if len(e) else "-"))
    mk, d, e = found["wedge"]
    assert len(e) > 0
    a, b = mk.normal[d.i[e]], mk.normal[d.j[e]]
    assert np.all(((a @ upper > 0.9) & (b @ lower > 0.9)) | ((a @ lower > 0.9) & (b @ upper > 0.9)))
    assert np.linalg.norm(mk.x[np.r_[d.i[e], d.j[e]], :2] - apex, axis=1).max() <= 4.0 * H
    assert len(found["square"][2]) == 0


def _face_markers(lo, hi, spacing=DESIGN_SPACING):
    """Markers on the six faces of the box from `lo` to `hi` [m], on each a grid `spacing` apart from half a spacing (+
    GRID_OFF) inside its corner, each with its face's outward normal (no marker anywhere else). Returns (positions (n,
    3), normals (n, 3))."""
    x, nrm = [], []
    for axis in range(3):
        u, v = [k for k in range(3) if k != axis]
        gu = np.arange(lo[u] + 0.5 * spacing, hi[u], spacing) + GRID_OFF
        gv = np.arange(lo[v] + 0.5 * spacing, hi[v], spacing) + GRID_OFF
        U, V = np.meshgrid(gu, gv, indexing="ij")
        for side, at in ((-1.0, lo[axis]), (1.0, hi[axis])):
            p = np.zeros((U.size, 3))
            p[:, axis], p[:, u], p[:, v] = at, U.ravel(), V.ravel()
            n = np.zeros(3)
            n[axis] = side
            x.append(p), nrm.append(np.tile(n, (len(p), 1)))
    return np.concatenate(x), np.concatenate(nrm)


VOID_CENTRE = np.array([-1.8e-3, 1.5e-3, 0.0])    # m, from the cube's centre: 2.6 mm or more of metal round the void
VOID_RADIUS = 0.6e-3                               # m


def test_columns_on_opposite_faces_of_thick_metal_are_not_opposed():
    """Plan 03b Task 7c (ruling T7-k): the stage 3 box's spurious THIN_METAL pairs in small. A 10 mm cube with markers
    on its six faces only holds a void VOID_RADIUS in radius, 2.6 mm or more from every face. The void's wall on the side
    of the -x and +y faces faces back into the body: its normal (out of the metal, into the void) has +x and -y parts
    beyond -OPPOSED_NORMALS, as the undersides of the box's rim slivers had, and no marker stands on it, as none stood on
    them. The same-side rule forbids that wall to the -x and +y faces' markers and gives it to the nearest it allows, on
    the +z and -z faces about 5 mm away; these tie on the cube's mid-plane, which crosses the wall, so their columns
    share a boundary of positive length there although their faces are 10 mm apart through metal. The pair keeps its
    boundary (the division's ownership is not this test's subject) but is not opposed: no metal there is thinner than
    twice the skin depth, and no edge is sharper than 60 degrees."""
    lo, hi = np.full(3, OFF - 5e-3), np.full(3, OFF + 5e-3)
    cube = box_sdf(lo, hi)
    centre = OFF + VOID_CENTRE
    sdf = lambda p: np.maximum(cube(p), VOID_RADIUS - np.linalg.norm(p - centre, axis=1))
    base = BoxBase.around(lo - BOX_H0, hi + BOX_H0, BOX_H0)
    mesh = _layer_mesh(base, sdf)
    phi = sdf(mesh.points)
    x, nrm = _face_markers(lo, hi)
    d = divide(mesh, phi, H, x, nrm)
    across = np.einsum("ek,ek->e", nrm[d.i], nrm[d.j]) < OPPOSED_NORMALS
    facing = np.flatnonzero(across & (d.length.max(axis=1) > 0.0))
    print("void in a cube: {} markers; {} boundaries of positive length between opposed columns, the longest {:.3f} mm "
          "(markers {:.2f} mm apart); {} opposed".format(
              len(x), len(facing), d.length[facing].max() * 1e3 if len(facing) else 0.0,
              np.linalg.norm(x[d.i[facing[0]]] - x[d.j[facing[0]]]) * 1e3 if len(facing) else 0.0, d.opposed.sum()))
    assert len(facing) > 0                                    # the mechanism is reproduced
    assert np.all(np.abs(nrm[d.i[facing], 2]) == 1.0)       # the +z and -z faces' columns
    assert not d.opposed.any()


THIN_STEP = 0.1e-3         # m: the thin test's samples along the normals in the test below (the division samples at its
                           # level spacing, 0.05 mm)


def test_the_thin_test_rejects_metal_that_is_thick_along_the_normals():
    """The final review's M12 (plan 03b Task 7c's deferred minor): the rejecting half of the thin test
    (`division._thin_between` with `_metal_along`), which no other test reaches. Pairs of markers of opposed normals +y
    and -y stand either side of points on a slab's mid-plane. 0.5 mm off each point, the depths of the point below their
    tangent planes sum to 1.0 mm, under twice the skin depth: in a slab 1.0 mm thick the metal along their normals from
    the point is 0.5 mm each way and the pair is thin; in a slab 6 mm thick the metal along each normal outlasts twice
    the skin depth (`_metal_along` gives inf) and the pair is not thin, whatever its tangent planes say. 2 mm off each
    point, in the thin slab, the tangent planes alone reject the pair. From a point in the gas 0.3 mm above a slab the
    metal along -y is 0."""
    from reentry_model.levelset.division import _metal_along, _thin_between
    at = OFF + np.array([[0.0, 0.0, 0.0], [0.5e-3, 0.0, -0.7e-3], [-1.3e-3, 0.0, 0.9e-3]])
    k = len(at)
    up, down = np.tile([0.0, 1.0, 0.0], (k, 1)), np.tile([0.0, -1.0, 0.0], (k, 1))
    i, j = np.arange(k), np.arange(k, 2 * k)
    reach = np.linspace(0.0, 2.0 * H, int(np.ceil(2.0 * H / THIN_STEP)) + 1)
    found = {}
    for thick in (1.0e-3, 6.0e-3):
        sdf = _slabs_sdf((OFF - 0.5 * thick, OFF + 0.5 * thick))
        mesh = _layer_mesh(BoxBase((-3e-3, -4e-3, -3e-3), 1e-3, (6, 8, 6)), sdf, fine=0.25e-3)
        phi = sdf(mesh.points)
        for off in (0.5e-3, 2e-3):
            x = np.concatenate([at + off * up, at + off * down])
            found[thick, off] = _thin_between(mesh, phi, H, x, np.concatenate([up, down]), i, j, at, THIN_STEP)
        found[thick, "along"] = _metal_along(mesh, phi, at, up, reach), _metal_along(mesh, phi, at, down, reach)
        found[thick, "gas"] = _metal_along(mesh, phi, at + [0.0, 0.5 * thick + 0.3e-3, 0.0], down, reach)
    print("thin test: 1.0 mm slab {} (metal {} mm up), 6 mm slab {} (metal {} up)".format(
        found[1.0e-3, 0.5e-3], np.round(found[1.0e-3, "along"][0] * 1e3, 9), found[6.0e-3, 0.5e-3],
        found[6.0e-3, "along"][0]))
    for a in found[1.0e-3, "along"]:
        np.testing.assert_allclose(a, 0.5e-3, rtol=0.0, atol=1e-12)
    assert np.all(np.isinf(found[6.0e-3, "along"][0])) and np.all(np.isinf(found[6.0e-3, "along"][1]))
    assert found[1.0e-3, 0.5e-3].all()                      # thin both ways
    assert not found[6.0e-3, 0.5e-3].any()                  # thick along the normals: rejected by the metal
    assert not found[1.0e-3, 2e-3].any()                    # far tangent planes: rejected before the metal is read
    assert np.all(found[1.0e-3, "gas"] == 0.0) and np.all(found[6.0e-3, "gas"] == 0.0)


def test_the_same_side_rule_needs_both_normals():
    """The final review's M12 (plan 03b Task 1's deferred minor): `ownership` given the markers' normals without the
    triangles' (or the triangles' without the markers') raises ValueError, not an opaque reshape error; given both, a
    marker allowed the triangle owns the whole of it."""
    import pytest
    x = np.array([[0.2e-3, 0.2e-3, 0.0]])
    tris = np.array([[[0.0, 0.0, 0.0], [1e-3, 0.0, 0.0], [0.0, 1e-3, 0.0]]])
    with pytest.raises(ValueError, match="tri_normal"):
        ownership(x, tris, [0], 1, 1e-3, normals=[[0.0, 0.0, 1.0]])
    with pytest.raises(ValueError, match="tri_normal"):
        ownership(x, tris, [0], 1, 1e-3, tri_normal=[[0.0, 0.0, 1.0]])
    W = ownership(x, tris, [0], 1, 1e-3, normals=[[0.0, 0.0, 1.0]], tri_normal=[[0.0, 0.0, 1.0]])
    assert W.toarray()[0, 0] == pytest.approx(0.5e-6, rel=1e-12)


def test_a_point_halfway_between_two_markers_is_owned_by_one():
    """Two markers 2 mm apart on a flat face, symmetric about a plane of the mesh's nodes, so that the corners of every
    triangle on that plane are exactly as near one marker as the other: each patch's owned areas, with and without the
    same-side rule, sum to its area to round-off, and so do the columns' areas at every level."""
    base = BoxBase((-6e-3, -2.5e-3, -6e-3), 0.5e-3, (24, 6, 24))
    mesh = build(base)
    phi = mesh.points[:, 1] - OFF
    x = np.array([[-1e-3, OFF, 0.0], [1e-3, OFF, 0.0]])
    nrm = np.tile([0.0, 1.0, 0.0], (2, 1))
    P = patches(mesh, phi)
    assert np.any(P.tris[:, :, 0] == 0.0)                                  # corners on the halfway plane
    size = OWN_SHARE * DESIGN_SPACING
    for W in (ownership(x, P.tris, P.tri_patch, len(P.area), size),
              ownership(x, P.tris, P.tri_patch, len(P.area), size, normals=nrm, tri_normal=P.normal[P.tri_patch])):
        assert np.abs(np.asarray(W.sum(axis=0)).ravel() / P.area - 1.0).max() <= 1e-12
    d = divide(mesh, phi, H, x, nrm)
    assert np.abs(d.area.sum(axis=0) / (12e-3 * 12e-3) - 1.0).max() <= 1e-12


def test_two_pieces_close_their_layers_separately():
    """Two 6 mm cubes 4 mm apart, off the grid, markers seeded on both: a closure factor for each piece, each within
    DIVISION_TOL of one, and each marker's piece the cube it stands on."""
    a_lo, b_lo = np.array([OFF - 8e-3, OFF - 3e-3, OFF - 3e-3]), np.array([OFF + 2e-3, OFF - 3e-3, OFF - 3e-3])
    side = np.full(3, 6e-3)
    fa, fb = box_sdf(a_lo, a_lo + side), box_sdf(b_lo, b_lo + side)
    sdf = lambda p: np.minimum(fa(p), fb(p))
    base = BoxBase.around((-8e-3 - BOX_H0, -3e-3 - BOX_H0, -3e-3 - BOX_H0), (8e-3 + BOX_H0, 3e-3 + BOX_H0,
                                                                              3e-3 + BOX_H0), BOX_H0)
    mesh = _layer_mesh(base, sdf)
    phi = sdf(mesh.points)
    mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
    d = divide(mesh, phi, H, mk.x, mk.normal)
    label, _ = pieces(mesh, phi)
    cube = [label[np.argmin(np.linalg.norm(mesh.centroids - (lo + 0.5 * side), axis=1))] for lo in (a_lo, b_lo)]
    print("two cubes: factors - 1 = {}".format(d.factor - 1.0))
    assert len(d.factor) == 2
    assert np.abs(d.factor - 1.0).max() <= DIVISION_TOL
    assert np.array_equal(d.piece, np.where(mk.x[:, 0] < 0.0, cube[0], cube[1]))


@functools.lru_cache(maxsize=1)
def _close_pair():
    """Two 6 mm cubes 0.7 mm apart across gas, off the grid (the review's case for Task 1's fix round 1): the pair's
    signed distance on a base mesh refined through both layers at FINE, markers seeded on both as the body seeds them.
    Returns (mesh, the pair's phi, cube B's own phi, the markers, which of them stand on B, B's lower corner)."""
    gap, side = 0.7e-3, np.full(3, 6e-3)
    a_lo = np.array([OFF - 6e-3 - 0.5 * gap, OFF - 3e-3, OFF - 3e-3])
    b_lo = np.array([OFF + 0.5 * gap, OFF - 3e-3, OFF - 3e-3])
    fa, fb = box_sdf(a_lo, a_lo + side), box_sdf(b_lo, b_lo + side)
    sdf = lambda p: np.minimum(fa(p), fb(p))
    base = BoxBase.around((-7e-3 - BOX_H0, -3e-3 - BOX_H0, -3e-3 - BOX_H0), (7e-3 + BOX_H0, 3e-3 + BOX_H0,
                                                                              3e-3 + BOX_H0), BOX_H0)
    mesh = _layer_mesh(base, sdf)
    phi = sdf(mesh.points)
    mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
    return mesh, phi, fb(mesh.points), mk, mk.x[:, 0] >= OFF, b_lo


def test_cells_too_thin_for_the_levels_keep_their_volume():
    """Task 1's fix round 1 (the review's Important 1): a 6 mm cube from x = 0.48 to 6.48 mm, y and z from -2.87 to 3.13
    mm (both of its faces at an edge equally off their planes of nodes), divided alone with its own markers. Beside each
    edge, at the depths where the two faces' planes of nodes lie equally deep, the P1 distance is constant over whole
    cells, which no level surface cuts: here they hold more than 1 % of the layer, and the levels alone left the columns
    2.1 % short. Placed whole (`division.THIN_SHARE`), the columns close the layer to within DIVISION_TOL."""
    mesh, _, phi, mk, on_b, _ = _close_pair()
    f = phi[mesh.tets]
    layer_volume = inside_volumes(mesh, phi) - inside_volumes(mesh, phi + H)
    flat = (f.max(axis=1) - f.min(axis=1) < 1e-7) & (layer_volume > 0.0)
    d = divide(mesh, phi, H, mk.x[on_b], mk.normal[on_b])
    print("cube B alone: {} markers; {} flat cells hold {:.2%} of the layer; factor - 1 = {:.2e}".format(
        on_b.sum(), flat.sum(), layer_volume[flat].sum() / layer_volume.sum(), d.factor[0] - 1.0))
    assert layer_volume[flat].sum() > 0.01 * layer_volume.sum()
    assert len(d.factor) == 1 and abs(d.factor[0] - 1.0) <= DIVISION_TOL


def test_no_column_crosses_from_one_piece_to_another():
    """Task 1's fix round 1 (the review's Important 2; the controller's ruling T1-i): the same two cubes 0.7 mm apart,
    divided together. Faces of the two cubes that lie side by side across the gap are parallel, so the same-side rule
    alone lets markers of one cube own the other's surface; the same-piece rule does not. Each marker's piece is the
    cube it stands on, no shared boundary joins the two cubes' markers, and each piece's factor is within
    DIVISION_TOL."""
    mesh, phi, _, mk, on_b, b_lo = _close_pair()
    label, _ = pieces(mesh, phi)
    cube_b = label[np.argmin(np.linalg.norm(mesh.centroids - (b_lo + 3e-3), axis=1))]
    P = patches(mesh, phi)
    W = ownership(mk.x, P.tris, P.tri_patch, len(P.area), OWN_SHARE * DESIGN_SPACING, normals=mk.normal,
                  tri_normal=P.normal[P.tri_patch]).tocoo()
    assert np.any(on_b[W.row] != (label[P.cell[W.col]] == cube_b))     # the same-side rule alone crosses the gap
    d = divide(mesh, phi, H, mk.x, mk.normal)
    print("the pair together: factors - 1 = {}".format(d.factor - 1.0))
    assert np.array_equal(d.piece == cube_b, on_b)
    assert len(d.i) > 0 and np.all(on_b[d.i] == on_b[d.j])
    assert len(d.factor) == 2 and np.abs(d.factor - 1.0).max() <= DIVISION_TOL
