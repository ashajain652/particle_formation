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
raise on it).

Pieces (the controller's ruling T1-i, Task 1's fix round 1). Two pieces of metal are parted by gas, so ownership never
crosses from one to the other: the level triangles of a piece, and its cells, are divided among that piece's markers
alone (`ownership` is called piece by piece). A marker's piece is the piece of the cell holding it
(`transfer.locate`); a marker on a face of a cell that holds no metal takes the piece of the patch of phi whose
centroid is nearest it. Every column then lies in one piece, and each piece's factor measures that piece's closure
alone. (With the same-side rule alone, faces of two cubes 0.7 mm apart that lie side by side across the gap are
parallel, and markers of each owned 0.39 % of the other's surface, the factors moving about 2e-3 apart; the review's
measurement.)

Cells too thin for the levels (Task 1's fix round 1). Where both faces at an edge stand equally off their planes of
nodes, the P1 distance is constant over whole cells beside the edge, at the depths where the two faces' planes of nodes
lie equally deep: no level cuts such a cell, and nearly flat ones are cut by one or two levels with a weight of
1 / |grad phi| far above one. On a 6 mm cube so placed, 474 cells constant to 0.1 um held 2.0 % of the layer and the
levels alone left the columns 2.1 % short (the review's measurement, 96 levels). So a cell whose nodes' depths spread
less than THIN_SHARE level spacings, and whose |grad phi| is not one (to UNIT_SLOPE: a cell of a flat face never is),
leaves the levels' sums: its exact volume in the layer goes whole to the column owning its centroid's foot (centroid -
phi g, by the nodes' means, under the same-side rule judged with g itself), shared between the two levels about its
mean depth with weights whose trapezoid integral is that volume.

Why each triangle's area is divided by |grad phi| (plan 03b Task 1, the controller's ruling T1-a; ruling 1 took the
area itself). The volume of {phi < -z} on the P1 mesh falls with z at the rate sum(area / |grad phi|) over the level
surface's triangles (the coarea formula), so that is the cross-section whose integral over z is the layer's volume. phi
is a distance, |grad phi| = 1, in every cell where it is linear, which is every cell of a flat face and nearly every
cell of a smooth one. But the P1 |grad phi| is not 1 in the cells along a sharp edge: there the interpolant of the
distance (the larger of the two faces' distances) has gradients from 0 to sqrt(2) long, and those cells hold more
volume than their level surfaces' area shows. With the level areas themselves, the columns of edged bodies fell short
of their layers by 0.7 to 5.2 % (two 6 mm cubes and the stage 3 box, on cells of 0.175 to 0.7 mm); with the areas over
|grad phi|, the shortfall is gone, and what remains is the levels' sampling of those edge cells (Task 1's report),
which the thin cells' exact volumes (below) take most of. That sets DIVISION_LEVELS (ruling T1-b; measured once each,
on the stage 3 box with markers seeded 2 mm apart, after the fix round): with 24, 48 and 96 levels the factor's
departure from one was -3.7e-4, +3.3e-3 and +5.6e-4 at 0.7 mm cells (divide 1.6, 3.1 and 6.1 s) and -3.6e-4, -5.0e-4
and +2.6e-4 at 0.175 mm (18, 29 and 53 s). It is not monotone in the levels, and it is not a bound for other bodies.

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

The base's ownership (plan 03b Task 5's fix round 1, the controller's ruling T5-j). The base level, phi = -H, is the
skin base beta's contour (`surfaces.skin_base` keeps phi + H's values at the nodes of the cells it cuts), so its patches
are beta's, the patches through which the volume's solve takes each skin's base law. The same ownership that gives each
column its base cross-section is kept with the plain area as its measure (`Division.base_owner`: the area of each patch
of the base each marker owns, each triangle's foot's shares of the triangle's own area, a foot of no area whole to the
marker owning its centroid), so that the group of beta beneath a skin is its own column's base, found once per division
(`skinbody.LevelSetSkinBody._beta_owners` reads it).

The columns' centroids (plan 03b Task 6, the controller's ruling T6-b, amending planning ruling 5). The sideways
conductance between two columns' cells divides by the distance between the columns' centroids over the cells' shared
depth span, so the division keeps each column's first moment at each level (`Division.moment`, `Division.centroid`): each
level triangle's share at the triangle's centroid on the level surface (the share's own centroid lies within the
triangle, a mesh cell across; plan 03b Task 6 measured what this leaves against sampled centroids), a thin cell's volume
at the cell's centroid, and the piece's factor applied as to the areas. Linear between the levels as the areas are, the
moment's integral over a depth span over the area's is the column's centroid there.

Numbers here are numerical settings; no physics."""
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree

from .criteria import nodal_average
from .cutcell import cell_gradients, inside_volumes, patches
from .levelset import pieces
from .markers import OPPOSED_NORMALS, OWN_SHARE, _nearest, ownership
from .transfer import barycentric, locate

DIVISION_LEVELS = 24       # level surfaces below phi down to H, less the one at phi itself (ruling 1; the controller's
                           # ruling T1-b: the fewest of 24, 48 and 96 that close the stage 3 box to DIVISION_TOL at
                           # both 0.7 and 0.175 mm cells, re-measured after Task 1's fix round 1; module docstring)
DIVISION_TOL = 1e-3        # a piece's closure factor at most this from one passes the step's check (ruling 1)
THIN_SHARE = 4.0           # of the level spacing: a cell whose nodes' depths spread less than this is too thin for the
                           # levels, and its volume is placed whole (module docstring; Task 1's fix round 1)
UNIT_SLOPE = 1e-6          # a cell whose |grad phi| is within this of one is a piece of a distance, its level surfaces
                           # exact at any depth, and never counts as thin (module docstring)

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
    base_cell: np.ndarray = None    # (Q,) the cell of each patch of the base level phi = -H (beta's patches; docstring)
    base_patch_area: np.ndarray = None  # (Q,) each such patch's area [m2]
    base_owner: object = None       # (n, Q) csr: the area [m2] of each such patch each marker owns, by the plain area
                                    # (module docstring: the base's ownership; not the piece's factor)
    moment: np.ndarray = None       # (n, B + 1, 3) each column's first moment at each level [m3]: `area` times the
                                    # centroid of its part of the level surface (module docstring; plan 03b Task 6, the
                                    # controller's ruling T6-b), its piece's factor applied as to `area`

    @property
    def base(self):
        """Each column's base [m2], its volume per unit depth at H (n,)."""
        return self.area[:, -1]

    @property
    def centroid(self):
        """Each column's centroid at each level [m] (n, B + 1, 3), the first moment over the cross-section; NaN at a
        level where the column has none (beyond a pinch)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(self.area[..., None] > 0.0, self.moment / self.area[..., None], np.nan)

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


def _marker_pieces(mesh, x, label, top):
    """Each marker's piece of metal (n,): the piece of the cell holding it (`transfer.locate`); for a marker on a face of
    a cell that holds no metal, the piece of the patch of phi (`top`) whose centroid is nearest it; -1 where phi has no
    patch."""
    piece = np.full(len(x), -1, dtype=np.int64)
    if not len(x):
        return piece
    piece = label[locate(mesh, x)[0]]
    lost = np.flatnonzero(piece < 0)
    if len(lost) and len(top.centroid):
        piece[lost] = label[top.cell[cKDTree(top.centroid).query(x[lost])[1]]]
    return piece


def divide(mesh, phi, depth, x, normal, levels=DIVISION_LEVELS, size=None):
    """The layer between phi and phi = -depth divided among the markers at `x` (n, 3) with normals `normal` (n, 3), at
    `levels` + 1 uniform depths (module docstring). phi must be a signed distance at least `depth` [m] deep, and a cell
    beyond. `size` [m] is `ownership`'s finest piece; by default OWN_SHARE of the markers' spacing, the square root of
    phi's area per marker. Returns a `Division`."""
    phi = np.asarray(phi, dtype=float)
    x = np.asarray(x, dtype=float).reshape(-1, 3)
    normal = np.asarray(normal, dtype=float).reshape(-1, 3)
    n = len(x)
    dz = depth / levels
    z = dz * np.arange(levels + 1)
    z[-1] = depth                                             # the base level is beta's contour exactly (docstring)
    w = np.full(levels + 1, dz)
    w[[0, -1]] *= 0.5                                         # the trapezoid rule's weights
    f = phi[mesh.tets]
    layer = np.flatnonzero((f.min(axis=1) < 0.0) & (f.max(axis=1) >= -depth))
    near = mesh.ring(layer, 1) if len(layer) else layer
    grad = cell_gradients(mesh, phi, near)
    G = nodal_average(mesh, grad, near)
    slope = np.zeros(len(mesh.tets))
    slope[near] = np.linalg.norm(grad, axis=1)                # |grad phi| per cell (a cut cell's is never zero)
    exact = inside_volumes(mesh, phi) - inside_volumes(mesh, phi + depth)     # each cell's volume in the layer [m3]
    thin = np.zeros(len(mesh.tets), dtype=bool)               # cells too thin in depth for the levels (docstring)
    thin[layer] = (f[layer].max(axis=1) - f[layer].min(axis=1) < THIN_SHARE * dz) & (exact[layer] > 0.0) & \
        (np.abs(slope[layer] - 1.0) > UNIT_SLOPE)
    label, n_pieces = pieces(mesh, phi)
    top = patches(mesh, phi, layer)
    piece = _marker_pieces(mesh, x, label, top)
    members = [np.flatnonzero(piece == p) for p in range(n_pieces)]            # each piece's markers (the same-piece
    trees = [cKDTree(x[mq]) if len(mq) else None for mq in members]           # rule, docstring)
    if size is None:
        size = OWN_SHARE * np.sqrt(top.area.sum() / max(n, 1))
    area = np.zeros((n, levels + 1))
    moment = np.zeros((n, levels + 1, 3))                     # each column's first moment at each level (docstring)
    found = []
    base_rows, base_cols, base_vals = [], [], []               # the base level's ownership by the plain area
    base_cell, base_patch_area = np.zeros(0, dtype=np.int64), np.zeros(0)
    for k, zk in enumerate(z):
        P = top if k == 0 else patches(mesh, phi + zk, layer)
        if not len(P.tris):
            continue
        Q = len(P.tris)
        tri_cell = P.cell[P.tri_patch]
        cell = np.repeat(tri_cell, 3)
        corners = P.tris.reshape(-1, 3)
        g = np.einsum("ni,nik->nk", barycentric(mesh, cell, corners), G[mesh.tets[cell]])
        feet = (corners + zk * g).reshape(Q, 3, 3)
        tri_area = 0.5 * np.linalg.norm(np.cross(P.tris[:, 1] - P.tris[:, 0], P.tris[:, 2] - P.tris[:, 0]), axis=1)
        own = np.where(thin[tri_cell], 0.0, tri_area / slope[tri_cell])       # each triangle's measure (docstring)
        tn = P.normal[P.tri_patch]
        tri_piece = label[tri_cell]
        tri_centroid = P.tris.mean(axis=1)                    # where each level triangle's share lies (docstring)
        rows, vals, moms = [], [], []
        for p, mq in enumerate(members):
            tq = np.flatnonzero(tri_piece == p)
            if not len(tq) or not len(mq):
                continue
            W, (bi, bj, bl) = ownership(x[mq], feet[tq], np.arange(len(tq)), len(tq), size, normals=normal[mq],
                                        tri_normal=tn[tq], boundaries=True, measure=P.tris[tq])
            W = W.tocoo()
            foot = np.bincount(W.col, W.data, minlength=len(tq))
            ok = foot[W.col] > 0.0
            r, c = W.row[ok], W.col[ok]
            share = W.data[ok] / foot[c]                      # each owner's share of its triangle's foot
            v = share * own[tq[c]]
            flat = np.flatnonzero(foot <= 0.0)                # carried to a segment or a point (or owned by none)
            if len(flat):
                o = _nearest(trees[p], feet[tq[flat]].mean(axis=1), 1, normal[mq], tn[tq[flat]])[1][:, 0]
                mine = o < len(mq)
                r, c = np.concatenate([r, o[mine]]), np.concatenate([c, flat[mine]])
                share = np.concatenate([share, np.ones(int(mine.sum()))])
                v = np.concatenate([v, own[tq[flat[mine]]]])
            rows.append(mq[r]), vals.append(v), moms.append(v[:, None] * tri_centroid[tq[c]])
            if k == levels:                                   # the base's ownership, by the plain area (docstring)
                base_rows.append(mq[r]), base_cols.append(P.tri_patch[tq[c]]), base_vals.append(share * tri_area[tq[c]])
            found.append((np.full(len(bi), k), mq[bi], mq[bj], bl))
        if rows:
            at, mom = np.concatenate(rows), np.concatenate(moms)
            area[:, k] = np.bincount(at, np.concatenate(vals), minlength=n)
            moment[:, k] = np.stack([np.bincount(at, mom[:, c], minlength=n) for c in range(3)], axis=1)
        if k == levels:
            base_cell, base_patch_area = P.cell.copy(), P.area.copy()

    # the thin cells' volumes, each whole to the column owning its centroid's foot, at the two levels about its depth
    for p, mq in enumerate(members):
        cells = np.flatnonzero(thin & (label == p))
        if not len(cells) or not len(mq):
            continue
        T = mesh.tets[cells]
        fc, gc = phi[T].mean(axis=1), G[T].mean(axis=1)
        o = _nearest(trees[p], mesh.points[T].mean(axis=1) - fc[:, None] * gc, 1, normal[mq], gc)[1][:, 0]
        mine = o < len(mq)
        s = np.clip(-fc[mine], 0.0, depth) / dz
        kk = np.minimum(np.floor(s).astype(np.int64), levels - 1)
        rows, V = mq[o[mine]], exact[cells[mine]]
        at = mesh.points[T[mine]].mean(axis=1)                # the cell's centroid, where its volume lies
        lower, upper = V * (1.0 - (s - kk)) / w[kk], V * (s - kk) / w[kk + 1]
        np.add.at(area, (rows, kk), lower)
        np.add.at(area, (rows, kk + 1), upper)
        np.add.at(moment, (rows, kk), lower[:, None] * at)
        np.add.at(moment, (rows, kk + 1), upper[:, None] * at)

    # the shared boundaries, at every level, and the pairs whose normals oppose that meet face to face (docstring)
    lev, bi, bj, bl = (np.concatenate(v) for v in zip(*found)) if found else (np.zeros(0, dtype=np.int64),) * 4
    labelled = layer[label[layer] >= 0]
    node_piece = np.full(len(phi), -1, dtype=np.int64)          # each node's piece, from a layer cell holding it
    for corner in range(4):
        node_piece[mesh.tets[labelled, corner]] = label[labelled]
    nodes = np.flatnonzero((phi < 0.0) & (phi >= -depth))
    owner = np.full(len(phi), n, dtype=np.int64)
    for p, mq in enumerate(members):
        sel = nodes[node_piece[nodes] == p]
        if not len(sel) or not len(mq):
            continue
        o = _nearest(trees[p], mesh.points[sel] - phi[sel, None] * G[sel], 1, normal[mq], G[sel])[1][:, 0]
        mine = o < len(mq)
        owner[sel[mine]] = mq[o[mine]]
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

    # closure, piece by piece (docstring)
    inside = label >= 0
    layer_volume = np.bincount(label[inside], exact[inside], minlength=n_pieces)
    counted = piece >= 0
    columns = np.bincount(piece[counted], (area @ w)[counted], minlength=n_pieces)
    with np.errstate(divide="ignore", invalid="ignore"):
        factor = np.where(columns > 0.0, layer_volume / columns, np.where(layer_volume > 0.0, np.inf, 1.0))
    scale = np.where(counted, factor[np.maximum(piece, 0)] if n_pieces else 1.0, 1.0)
    area, moment = area * scale[:, None], moment * scale[:, None, None]
    cat = lambda parts, dtype: np.concatenate(parts) if parts else np.zeros(0, dtype=dtype)
    base_owner = sp.csr_matrix((cat(base_vals, float), (cat(base_rows, np.int64), cat(base_cols, np.int64))),
                               shape=(n, len(base_cell)))
    return Division(z, area, area @ w, i, j, length, _boundary_normals(x, normal, i, j),
                    np.einsum("ek,ek->e", normal[i], normal[j]) < OPPOSED_NORMALS, piece, factor, base_cell,
                    base_patch_area, base_owner, moment)
