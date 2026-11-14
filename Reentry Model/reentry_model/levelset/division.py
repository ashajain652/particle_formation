"""The division of the layer between the gas surface and the skin base among the markers (design section 21.5; plan 03b
Task 1, planning ruling 1).

phi, the gas surface, is a signed distance (negative in the metal) reaching at least the skin depth H into the metal;
the skin base beta is its level surface phi = -H. A point of the layer between them belongs to the marker that owns the
point of phi nearest it (C1's `markers.ownership`, with the same-side rule: never a marker whose normal opposes the
surface's there), and its depth is its distance from phi. Column i is the part marker i owns; its cross-section A_i(z)
is that part's volume per unit depth at depth z (its area, where phi is a distance in every cell): A_i(0) at phi and
A_i(H) its base. One division serves everything that follows it: the skins' cross-sections and bases, and the boundaries
across which runoff and sideways heat pass.

How it is integrated (ruling 1). The layer is cut by `levels` + 1 uniform level surfaces phi = -z_k from 0 to H
(`cutcell.patches` of phi + z_k). Each level surface's triangles are carried to their feet on phi, each corner p to
p + z_k g(p), g phi's recovered gradient (the volume-weighted nodal average of the cells' gradients, interpolated), and
the feet are divided among the markers by `ownership`. A triangle's share for a marker is its foot's share (the map
between the two is affine, so it keeps area fractions) of the triangle's measure on the level surface: its area over
the length of its cell's gradient of phi, |grad phi| (below). A shared boundary's length at a level is the length of
the segments dividing the feet, carried back to the level surface by the same maps and measured there (`ownership`'s
`measure`), summed over the pair's segments. A foot carried to a segment or a point, of no area, goes whole to the
marker owning its centroid. Each column's volume is the trapezoid integral of its areas over the levels. The areas of
each piece of metal (`levelset.pieces`) are then scaled by one factor so that its columns' volumes sum to its layer's
volume, inside_volumes(phi) - inside_volumes(phi + H). The factor's departure from one is the closure measure of design
section 21.8 item 5; more than DIVISION_TOL fails a step's check, which is the caller's (the body records it; the tests
raise on it). A column is counted with the piece holding most of its volume.

Why each triangle's area is divided by |grad phi| (plan 03b Task 1, the controller's ruling T1-a; ruling 1 took the
area itself). The volume of {phi < -z} on the P1 mesh falls with z at the rate sum(area / |grad phi|) over the level
surface's triangles (the coarea formula), so that is the cross-section whose integral over z is the layer's volume. phi
is a distance, |grad phi| = 1, in every cell where it is linear, which is every cell of a flat face and nearly every
cell of a smooth one. But the P1 |grad phi| is not 1 in the cells along a sharp edge: there the interpolant of the
distance (the larger of the two faces' distances) has gradients from 0 to sqrt(2) long, and those cells hold more
volume than their level surfaces' area shows. With the level areas themselves, the columns of edged bodies fell short
of their layers by 0.7 to 5.2 % (two 6 mm cubes and the stage 3 box, on cells of 0.175 to 0.7 mm); with the areas over
|grad phi|, the shortfall is gone, and what remains is the levels' sampling of those edge cells (Task 1's report).
That sampling sets DIVISION_LEVELS (ruling T1-b; measured once each, on the stage 3 box with markers seeded 2 mm
apart): with 24, 48 and 96 levels the factor's departure from one was -2.0e-3, +3.3e-3 and +5.6e-4 at 0.7 mm cells
(divide 2.9, 5.6 and 11.4 s) and -6.6e-4, -8.1e-5 and -1.6e-5 at 0.175 mm (32, 54 and 92 s). It is not monotone in the
levels, and it is not a bound for other bodies: two 6 mm cubes at 0.7 mm cells closed to only 2.5e-3 at 96 levels.

Why the recovered gradient, not each cell's own. Carried by its own cell's gradient, each level triangle would move
rigidly, and the feet of neighbouring triangles would part where the gradient turns. Across a convex edge the level
surface's triangles beside the edge carry to the faces and those of the edge itself to the edge; the boundary between
a column of one face and a column of the other then falls on such a parting, in no foot, and is lost. With the
recovered gradient, neighbouring triangles' feet share their corners and tile phi without gaps, and that boundary is a
dividing segment like any other, in the triangles beside the edge.

Metal thinner than twice H (design section 21.11). Where two faces whose normals oppose lie closer than 2H, their
columns meet face to face on the surface halfway between them. That surface is parallel to the level surfaces, so it has
no length on any of them and the levels alone do not find it. So each node of the layer (-H <= phi < 0) is also given the
marker owning its foot, node - phi g (the same-side rule judged with g itself, which is short where the normal is
undefined, between the faces, and then allows every marker). Two columns whose normals oppose (n_i . n_j <
`markers.OPPOSED_NORMALS`) and whose nodes an edge of the mesh joins share a boundary too: of zero length at every level
where the levels found none, and marked `opposed`, for the caller to stop the step. It reports metal thinner than 2H as
the P1 mesh resolves it, to within about one cell where the cells are coarser than beta's interior.

Numbers here are numerical settings; no physics."""
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from .criteria import nodal_average
from .cutcell import cell_gradients, inside_volumes, patches
from .levelset import pieces
from .markers import OPPOSED_NORMALS, OWN_SHARE, _nearest, ownership
from .transfer import barycentric

DIVISION_LEVELS = 96       # level surfaces below phi down to H, less the one at phi itself (ruling 1 named 24; the
                           # controller's ruling T1-b: the fewest of 24, 48 and 96 that closed the stage 3 box to
                           # DIVISION_TOL at both 0.7 and 0.175 mm cells, its divide at most four times as long as at 24;
                           # module docstring)
DIVISION_TOL = 1e-3        # a piece's closure factor at most this from one passes the step's check (ruling 1)

_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))     # a tetrahedron's edges, by its corners


@dataclass(frozen=True)
class Division:
    """The layer divided among n markers at B + 1 levels (module docstring), with E shared boundaries.

    `area` is each column's volume per unit depth at each level (the coarea measure, module docstring), the
    cross-section whose integral is the column's volume: what the skins fill. On a flat face it is the column's area;
    beside an edge it differs from the area of phi the marker owns (C1's `markers.ownership`), which the gas loads and
    the heat keep (plan 03b Task 6; ruling T1-e)."""
    levels: np.ndarray     # (B + 1,) the levels' depths below phi [m], 0 to H
    area: np.ndarray       # (n, B + 1) each column's volume per unit depth at each level [m2], its piece's factor applied
    volume: np.ndarray     # (n,) each column's volume [m3], the trapezoid integral of its areas
    i: np.ndarray          # (E,) the first marker of each shared boundary
    j: np.ndarray          # (E,) its second, i < j
    length: np.ndarray     # (E, B + 1) each boundary's length at each level [m]
    normal: np.ndarray     # (E, 3) unit, from i toward j along the surface
    opposed: np.ndarray    # (E,) bool: the pair's normals oppose (design section 21.11)
    piece: np.ndarray      # (n,) each column's piece of metal (`levelset.pieces`), -1 where there is none
    factor: np.ndarray     # (P,) each piece's factor on its columns' areas (its closure measure)

    @property
    def base(self):
        """Each column's base [m2], its volume per unit depth at H (n,)."""
        return self.area[:, -1]

    def cross_section(self, z):
        """Each column's cross-section, its volume per unit depth [m2], at depths z [m] below phi, z (n, ...) (or one
        depth for all): linear in z between the levels, the top's above phi and zero below H."""
        n = len(self.area)
        z = np.asarray(z, dtype=float)
        z = np.broadcast_to(z, (n,)) if z.ndim == 0 else z
        flat = z.reshape(n, -1)
        k = np.clip(np.searchsorted(self.levels, flat, side="right") - 1, 0, len(self.levels) - 2)
        z0, z1 = self.levels[k], self.levels[k + 1]
        s = np.clip((flat - z0) / (z1 - z0), 0.0, 1.0)
        row = np.arange(n)[:, None]
        out = (1.0 - s) * self.area[row, k] + s * self.area[row, k + 1]
        return np.where(flat > self.levels[-1], 0.0, out).reshape(z.shape)


def _boundary_normals(x, normal, i, j):
    """Unit vectors from marker i toward marker j along the surface (E, 3): (t_ij - t_ji) / |t_ij - t_ji|, t_ij the
    unit direction from i toward j in i's tangent plane (the faces of `Markers.tessellate`), so that a boundary across an
    edge points round it; the straight direction where that vanishes (normals opposed along the line between)."""
    d = x[j] - x[i]
    unit = lambda v: v / np.maximum(np.linalg.norm(v, axis=1), 1e-300)[:, None]
    ti = unit(d - np.einsum("ek,ek->e", d, normal[i])[:, None] * normal[i])
    tj = unit(-d + np.einsum("ek,ek->e", d, normal[j])[:, None] * normal[j])
    v = ti - tj
    flat = np.linalg.norm(v, axis=1) < 1e-6
    v[flat] = d[flat]
    return unit(v)


def divide(mesh, phi, depth, x, normal, levels=DIVISION_LEVELS, size=None):
    """The layer between phi and phi = -depth divided among the markers at `x` (n, 3) with normals `normal` (n, 3), at
    `levels` + 1 uniform depths (module docstring). phi must be a signed distance at least `depth` [m] deep, and a cell
    beyond. `size` [m] is `ownership`'s finest piece; by default OWN_SHARE of the markers' spacing, the square root of
    phi's area per marker. Returns a `Division`."""
    phi = np.asarray(phi, dtype=float)
    x = np.asarray(x, dtype=float).reshape(-1, 3)
    normal = np.asarray(normal, dtype=float).reshape(-1, 3)
    n = len(x)
    z = depth * np.arange(levels + 1) / levels
    w = np.full(levels + 1, depth / levels)
    w[[0, -1]] *= 0.5                                         # the trapezoid rule's weights
    f = phi[mesh.tets]
    layer = np.flatnonzero((f.min(axis=1) < 0.0) & (f.max(axis=1) >= -depth))
    near = mesh.ring(layer, 1) if len(layer) else layer
    grad = cell_gradients(mesh, phi, near)
    G = nodal_average(mesh, grad, near)
    slope = np.zeros(len(mesh.tets))
    slope[near] = np.linalg.norm(grad, axis=1)                # |grad phi| per cell (a cut cell's is never zero)
    label, n_pieces = pieces(mesh, phi)
    tree = cKDTree(x)
    area = np.zeros((n, levels + 1))
    held = np.zeros((n, n_pieces))                            # each column's volume in each piece [m3]
    found = []
    top = None
    for k, zk in enumerate(z):
        P = patches(mesh, phi + zk, layer)
        if not len(P.tris):
            continue
        Q = len(P.tris)
        cell = np.repeat(P.cell[P.tri_patch], 3)
        corners = P.tris.reshape(-1, 3)
        g = np.einsum("ni,nik->nk", barycentric(mesh, cell, corners), G[mesh.tets[cell]])
        feet = (corners + zk * g).reshape(Q, 3, 3)
        tri_area = 0.5 * np.linalg.norm(np.cross(P.tris[:, 1] - P.tris[:, 0], P.tris[:, 2] - P.tris[:, 0]), axis=1)
        own = tri_area / slope[P.cell[P.tri_patch]]           # each triangle's measure (docstring)
        if top is None:
            top = P
            size = OWN_SHARE * np.sqrt(tri_area.sum() / max(n, 1)) if size is None else size
        tn = P.normal[P.tri_patch]
        W, (bi, bj, bl) = ownership(x, feet, np.arange(Q), Q, size, normals=normal, tri_normal=tn, boundaries=True,
                                    measure=P.tris)
        W = W.tocoo()
        foot = np.bincount(W.col, W.data, minlength=Q)
        ok = foot[W.col] > 0.0
        rows, cols = W.row[ok], W.col[ok]
        vals = W.data[ok] / foot[cols] * own[cols]
        flat = np.flatnonzero(foot <= 0.0)                    # carried to a segment or a point (or owned by none)
        if len(flat):
            _, o = _nearest(tree, feet[flat].mean(axis=1), 1, normal, tn[flat])
            mine = o[:, 0] < n
            rows = np.concatenate([rows, o[mine, 0]])
            cols = np.concatenate([cols, flat[mine]])
            vals = np.concatenate([vals, own[flat[mine]]])
        area[:, k] = np.bincount(rows, vals, minlength=n)
        np.add.at(held, (rows, label[P.cell[P.tri_patch]][cols]), w[k] * vals)
        found.append((np.full(len(bi), k), bi, bj, bl))

    # the shared boundaries, at every level, and the pairs whose normals oppose that meet face to face (docstring)
    lev, bi, bj, bl = (np.concatenate(v) for v in zip(*found)) if found else (np.zeros(0, dtype=np.int64),) * 4
    nodes = np.flatnonzero((phi < 0.0) & (phi >= -depth))
    owner = np.full(len(phi), n, dtype=np.int64)
    if len(nodes) and n:
        owner[nodes] = _nearest(tree, mesh.points[nodes] - phi[nodes, None] * G[nodes], 1, normal, G[nodes])[1][:, 0]
    T = mesh.tets[layer]
    ci, cj = [], []
    for a, b in _EDGES:
        oa, ob = owner[T[:, a]], owner[T[:, b]]
        meet = (oa < n) & (ob < n) & (oa != ob)
        oa, ob = oa[meet], ob[meet]
        apart = np.einsum("ek,ek->e", normal[oa], normal[ob]) < OPPOSED_NORMALS
        ci.append(np.minimum(oa, ob)[apart]), cj.append(np.maximum(oa, ob)[apart])
    key = np.unique(np.concatenate([bi * n + bj] + [u * n + v for u, v in zip(ci, cj)]).astype(np.int64))
    length = np.zeros((len(key), levels + 1))
    np.add.at(length, (np.searchsorted(key, bi * n + bj), lev.astype(np.int64)), bl)
    i, j = key // n, key % n

    # pieces and closure (docstring)
    piece = np.where(held.sum(axis=1) > 0.0, np.argmax(held, axis=1) if n_pieces else -1, -1)
    lost = np.flatnonzero(piece < 0)
    if len(lost) and top is not None:
        _, t = cKDTree(top.centroid).query(x[lost])
        piece[lost] = label[top.cell[t]]
    layer_volume = np.bincount(label[label >= 0], (inside_volumes(mesh, phi) - inside_volumes(mesh, phi + depth))
                               [label >= 0], minlength=n_pieces)
    columns = np.bincount(piece[piece >= 0], (area @ w)[piece >= 0], minlength=n_pieces)
    with np.errstate(divide="ignore", invalid="ignore"):
        factor = np.where(columns > 0.0, layer_volume / columns, np.where(layer_volume > 0.0, np.inf, 1.0))
    area = area * np.where(piece >= 0, factor[np.maximum(piece, 0)] if n_pieces else 1.0, 1.0)[:, None]
    return Division(z, area, area @ w, i, j, length, _boundary_normals(x, normal, i, j),
                    np.einsum("ek,ek->e", normal[i], normal[j]) < OPPOSED_NORMALS, piece, factor)
