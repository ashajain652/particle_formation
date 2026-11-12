"""The level set's markers (design sections 9.1, 9.2 and 9.6; level-set stage 3 Task 5): Lagrangian points on phi, the
gas surface, each carrying one skin (a row of a `Skins`) and an area (Ruling 1: a marker sits on phi and carries its
skin's base at x - H n, H the skin's thickness).

Motion. `carry` moves the markers by a velocity field; `snap` returns them to phi's zero contour along the normal,
x <- x - phi(x) n(x), n the unit recovered gradient, applied once and iterated once more (section 9.1); a marker then
within NEWTON_REACH of the contour is finished by Newton steps along the gradient of the cell it lies in,
x <- x - phi g / |g|^2, the exact root of the P1 field within the cell (`snap(reach=...)`: the body's last snap round
takes any reach). Task 9, fix round 4: the recovered step assumes phi rises one metre per metre along n; where phi's
cells are far from a distance (|grad phi| 0.20 and 0.88, the slope along n 0.17 and 0.25, on the full-size cube's
second step) it removes only a sixth to a quarter of the distance a pass, and 40 passes left two markers 4e-8 and
4e-10 m off phi against the body's 1e-10 m. The reach keeps the landing where the recovered path lands: Newton from a
millimetre off moves a marker along a gradient up to 73 degrees off n, and on that step's markers landed up to 1.5 mm
(a tenth of them over 0.1 mm) from where 400 recovered passes land, against 2.6 um with a reach of 1 um. phi is read
as its P1 interpolant; the recovered gradient is the nodal (volume-weighted) average of the cells' gradients,
interpolated, so the normal varies continuously along the surface.

Tessellation (Ruling 2). Each marker's neighbours are the markers within NEIGHBOUR_REACH target spacings whose patches
are joined by a path along phi's patch graph (cells sharing a node, each step the distance between the two patches'
centroids; `patch_paths`) no longer than NEIGHBOUR_PATH reaches, so markers across a thin gap of metal or gas, which
the surface joins only far round, are never neighbours. Task 9b: the rule before allowed at most two hops of the graph,
but the patches are pieces of cut cells, on the body's meshes about a sixth of the 2 mm spacing across, so two hops
reached about one spacing. On the 10 mm test cube 123 of its 219 cells were capped and the areas were 0.38 to 2.67
times the area nearest each marker with nothing moving; at the step ends 64 to 95 per cent of the area changes the
tessellation itself made sat on markers whose neighbours had crossed the hop limit. With the path, no cell is capped
there and the areas are 0.68 to 1.10 times the nearest-marker areas (the least at the rounded edges and corners, where
a tangent plane stands in poorly for the surface); on a 10 mm sphere on the design's 0.7 mm cells with markers at
2 mm, 0.98 to 1.02 and none capped. The bound is two reaches, not less, because the centroids' path zigzags across
sliver patches: with one marker per patch on a 25 mm sphere at 2 mm cells it is 1.04 times the straight distance at the
median and 1.7 at the 95th percentile, and at 1.5 reaches a rotation of those markers lost 306 far pairs and changed
areas by 8 %; at two it loses none. A slot 1.5 mm wide keeps its two walls apart from 1 mm above its floor, where the
path round the floor passes the bound. The plan named 1.5 spacings for the reach; with one marker per patch at creation (the
patches' areas spread over three decades) the true Voronoi cells reach past 0.75 spacings on most markers, so at 1.5
the cap bound 78 % of the cells on a 25 mm sphere at 2 mm cells, a capped cell could not grow with the surface (a
uniform stretch to 26 mm gave 6 % less than (26/25)^2) and a rotation changed areas by 73 % (7 % with the two hops). At
3.0 none is capped there (0.55 % at 0.7 mm cells near the surface, one marker per patch), no pair within the reach is
lost to the path bound, and the rotation changes areas by 0.3 %. In each marker's tangent plane the neighbours' offsets, projected, bound its Voronoi
cell (the half-planes y . d <= |d|^2 / 2), clipped from a CAP_SIDES-gon whose inscribed radius is half the reach (no
neighbour beyond the reach could bound the cell inside it); a cell that the cap still bounds is flagged `capped`. Each
edge is a face to the neighbour that made it. Faces are symmetrised: a pair's face has the mean of its two edge lengths
(an edge present on one side only counts as zero on the other) and the unit normal (t_ij - t_ji)/|t_ij - t_ji|, t_ij the
unit tangent-plane direction from i toward j in i's plane, so a face is shared by both markers with opposite normals.
The cells' areas are scaled by one factor so that they sum to the area of phi's patches the markers cover (those whose
centroid lies within the cap of its nearest marker: all of them on a surface markers cover), and the face lengths by its
square root.

Topology. `split` divides a marker above twice the target area into two half a target spacing apart along the longest
direction of its Voronoi cell, each with the same profile and half the mass and area; `merge` joins markers closer than
a third of the target spacing, their skins mixed by the hand-over rule (`Skins.combine`); `trapped_markers` finds the
markers with no zero contour of phi within half a cell, and `release` removes markers with their skins and returns what
the skins held, for the caller to book as a deferred load in the cell holding each. Every change of the markers changes
the skins' rows with them and is counted (`splits`, `merges`, `trapped`).

Re-partition (Task 9b, the controller's ruling of 2026-11-12). `transfer_matrix` follows the area a tessellation hands
out afresh, for the body's step end (`stepend.StepEnd._repartition`, which moves the skins' metal with it). The surface
before and after is sampled (`surface_samples`). Each sample of the old surface belongs to the marker nearest it at its
old position (`x_tess`); each sample of the new surface to its old owner, the nearest marker at its old position carried
along its old normal by the normal part of its move (so that a receding face keeps its owners), and to its new owner,
the nearest at its new position. A marker's stretch is the new surface's area its old owners hold over the old
surface's, each summed over the marker and its neighbours (the new faces) before the ratio, then scaled, on each piece
of the faces' graph, so that the old areas times the stretch sum to the new areas exactly (a marker whose neighbourhood
holds no sample takes the piece's ratio of new to old area). The transfer starts from the new-surface samples whose old
and new owners are one marker or neighbours, each row scaled to its old area times its stretch (a row with none keeps
that area itself); its columns are then brought to the new areas by a potential flow across the faces (`_fit_columns`),
which moves area only within rows and turns no entry negative. Rows and columns hold to round-off, and when nothing moves
the transfer is the old areas on its diagonal. Measured once each (Task 9b): a cold step end of the test cube, markers
moved 13 nm, passed 3.8e-7 of the area between markers; a 10 mm sphere's markers carried onto 10.4 mm found stretches
of 1.0749 to 1.0879 against (10.4/10)^2 = 1.0816, with 0.35 % of the area passed between markers.

Numbers here are numerical settings of the construction; no physics."""
from collections import namedtuple

import numpy as np
from scipy.sparse import coo_matrix, diags, identity
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree

from .criteria import nodal_average, principal_curvatures
from .cutcell import cell_gradients, patches as cut_patches
from .levelset import point_triangle_distance
from .skin import Skins
from .thin import patch_adjacency
from .transfer import locate

NEIGHBOUR_REACH = 3.0      # target spacings: markers farther apart are never neighbours (plan: 1.5; see the docstring)
NEIGHBOUR_PATH = 2.0       # the path along phi's patch graph between neighbours' patches at most this times the reach
                           # (Task 9b; see the docstring)
PATH_CHUNK = 64            # sources of one Dijkstra in `patch_paths` at most
CAP_SIDES = 16             # sides of the polygon a Voronoi cell is clipped from
SNAP_PASSES = 2            # the snap applied once and iterated once (section 9.1)
NEWTON_ROUNDS = 8          # cell-gradient Newton steps after them at most (`_cell_newton`; Task 9, fix round 4)
NEWTON_REACH = 1.0e-6      # m: only markers this close to phi's zero contour take them (module docstring)
NEWTON_TOL = 1.0e-13       # m: points this close to phi's zero contour take no Newton step (a thousandth of the body's
                           # snap tolerance, stepend.SNAP_TOL)
SPLIT_FACTOR = 2.0         # a marker splits above this times the target area (section 9.2)
MERGE_FRACTION = 1.0 / 3.0  # markers closer than this times the target spacing merge (section 9.2)
SPLIT_GAP = 0.5            # target spacings between the two halves of a split (plan, Task 5)
TRAP_FRACTION = 0.5        # cells: a marker with no zero contour this close is trapped (section 9.6)
TRAP_CANDIDATES = 8        # nearest patch triangles (by centroid) searched first for the trapped test
SAMPLE_LEVEL_MAX = 5       # midpoint subdivisions of a patch triangle at most for the surface samples (`surface_samples`)
FIT_PASSES = 4             # flows at most that bring the re-partition's columns to the new areas (`_fit_columns`)
FIT_RTOL = 1e-14           # of a column's area: the fit stops when no column is farther from it than this
FIT_FAIL = 1e-12           # of a column's area: a fit still farther from it after FIT_PASSES raises

Faces = namedtuple("Faces", "i j length normal")


def _field(mesh, phi, x):
    """phi's P1 value and the unit recovered gradient at the points x (n, 3)."""
    cell, b = locate(mesh, x)
    near = mesh.ring(cell, 1)
    G = nodal_average(mesh, cell_gradients(mesh, phi, near), near)
    T = mesh.tets[cell]
    val = np.einsum("ni,ni->n", b, phi[T])
    g = np.einsum("ni,nik->nk", b, G[T])
    return val, g / np.maximum(np.linalg.norm(g, axis=1), 1e-300)[:, None]


def _cell_newton(mesh, phi, x):
    """The points x (n, 3) moved onto phi's P1 zero contour by Newton steps along the gradient of the cell each lies
    in, x <- x - phi(x) g / |g|^2 with g that cell's (constant) gradient: the exact root of the P1 field within the cell,
    reached in one step when the step stays in the cell; a step that leaves it is repeated from the new cell, at most
    NEWTON_ROUNDS rounds, the points within NEWTON_TOL of the contour (or in a cell where phi is flat) left alone. (Task 9,
    fix round 4: the recovered normal is not the cell's gradient and assumes a unit slope, so where the slope along it
    is far from one its steps converge slowly or, across a kink of the P1 contour, cycle; module docstring.)"""
    x = np.array(x, dtype=float)
    todo = np.arange(len(x))
    for _ in range(NEWTON_ROUNDS):
        if not len(todo):
            break
        cell, b = locate(mesh, x[todo])
        val = np.einsum("ni,ni->n", b, phi[mesh.tets[cell]])
        g = cell_gradients(mesh, phi, cell)
        gg = np.einsum("nk,nk->n", g, g)
        go = (np.abs(val) > NEWTON_TOL) & (gg > 0.0)
        todo, val, g, gg = todo[go], val[go], g[go], gg[go]
        x[todo] -= (val / gg)[:, None] * g
    return x


def _subdivision(level):
    """Barycentric centroids (4**level, 3) of the triangles a triangle is cut into by its edge midpoints `level` times."""
    tri = np.eye(3)[None, :, :]
    for _ in range(level):
        a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
        ab, bc, ca = 0.5 * (a + b), 0.5 * (b + c), 0.5 * (c + a)
        tri = np.concatenate([np.stack(t, axis=1) for t in ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])
    return tri.mean(axis=1)


def surface_samples(mesh, phi, size):
    """Points sampling phi's zero contour, with the area each stands for: every patch triangle cut into 4**L equal
    triangles by its edge midpoints, L the least level (at most SAMPLE_LEVEL_MAX) at which they are at most `size` [m]
    across (the square root of twice their area), each at its centroid. Returns (points (k, 3), areas (k,) [m2]); the
    areas sum to the patches' area to round-off."""
    t = cut_patches(mesh, phi).tris
    a = 0.5 * np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1)
    level = np.clip(np.ceil(np.log2(np.maximum(np.sqrt(2.0 * a), 1e-300) / size)), 0, SAMPLE_LEVEL_MAX).astype(int)
    y, w = [np.zeros((0, 3))], [np.zeros(0)]
    for L in np.unique(level):
        sel = level == L
        bary = _subdivision(int(L))
        y.append(np.einsum("sk,qkd->qsd", bary, t[sel]).reshape(-1, 3))
        w.append(np.repeat(a[sel] / len(bary), len(bary)))
    return np.concatenate(y), np.concatenate(w)


def patch_paths(mesh, P, a, b, limit):
    """The shortest path [m] along phi's patch graph (`thin.patch_adjacency`: patches whose cells share a node, each
    step the distance between the two patches' centroids) from patch a[k] to patch b[k] of `P` (a `cutcell.Patches`),
    np.inf where it is longer than `limit` [m]. Dijkstra from the sources in spatial chunks, each on the patches within
    `limit` of its sources (a path no longer than `limit` never leaves them), so that the work and memory follow the
    sources' neighbourhoods, not the whole surface."""
    a, b = np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64)
    out = np.full(len(a), np.inf)
    if not len(a):
        return out
    A = patch_adjacency(mesh, P).tocoo()
    c = P.centroid
    w = np.maximum(np.linalg.norm(c[A.row] - c[A.col], axis=1), 1e-300)    # an explicit zero would be no edge
    G = coo_matrix((w, (A.row, A.col)), shape=A.shape).tocsr()
    tree = cKDTree(c)
    src = np.unique(a)
    key = np.floor((c[src] - c.min(axis=0)) / limit).astype(np.int64)       # sources in spatial chunks
    src = src[np.lexsort((key[:, 2], key[:, 1], key[:, 0]))]
    for k in range(0, len(src), PATH_CHUNK):
        chunk = np.sort(src[k:k + PATH_CHUNK])
        balls = tree.query_ball_point(c[chunk], limit)
        near = np.unique(np.concatenate([np.asarray(q, dtype=np.int64) for q in balls] + [chunk]))
        D = dijkstra(G[near][:, near], indices=np.searchsorted(near, chunk), limit=limit)
        sel = np.flatnonzero(np.isin(a, chunk))
        row = np.searchsorted(chunk, a[sel])
        col = np.minimum(np.searchsorted(near, b[sel]), len(near) - 1)
        inside = near[col] == b[sel]
        out[sel[inside]] = D[row[inside], col[inside]]
    return out


def transfer_matrix(old_x, old_normal, old_area, old_samples, x, area, faces, samples):
    """The step end's re-partition (Task 9b; module docstring, Re-partition): n markers, the same before and after, each
    owning before the area `old_area` [m2] (its skin's) about `old_x` (n, 3) with the normal `old_normal` on the
    surface sampled by `old_samples`, and after the tessellated `area` about `x` on the surface sampled by `samples`
    (each a (points, areas) pair from `surface_samples`), with the tessellation's `faces`. Returns (T, stretch): T
    (n, n) csr [m2], T[j, i] the area of marker j's old column that marker i's new column receives, its rows summing to
    old_area * stretch and its columns to `area`, each to round-off; `stretch` (n,) each old column's stretch."""
    old_x, old_normal, x = (np.asarray(v, dtype=float).reshape(-1, 3) for v in (old_x, old_normal, x))
    old_area, area = np.asarray(old_area, dtype=float), np.asarray(area, dtype=float)
    n = len(x)
    (y0, w0), (y1, w1) = old_samples, samples
    carried = old_x + np.einsum("ij,ij->i", x - old_x, old_normal)[:, None] * old_normal
    _, own0 = cKDTree(old_x).query(y0)                      # the old surface's samples by their old owners
    _, own1 = cKDTree(carried).query(y1)                    # the new surface's by theirs, carried along the old normal
    _, new1 = cKDTree(x).query(y1)                          # and by their new owners
    fi, fj = (np.zeros(0, np.int64),) * 2 if faces is None else (faces.i, faces.j)
    F = coo_matrix((np.ones(len(fi)), (fi, fj)), shape=(n, n))
    W = ((F + F.T + identity(n)) > 0).astype(float).tocsr()
    ncomp, comp = connected_components(W, directed=False)
    with np.errstate(divide="ignore", invalid="ignore"):   # each marker's stretch, smoothed over its neighbours
        s = (W @ np.bincount(own1, w1, minlength=n)) / (W @ np.bincount(own0, w0, minlength=n))
    R = np.zeros(n)
    for c in range(ncomp):                                  # normalised to the covered area's change, piece by piece
        m = comp == c
        A_new, A_old = area[m].sum(), old_area[m].sum()
        sc = np.where(np.isfinite(s[m]) & (s[m] > 0.0), s[m], A_new / A_old if A_old > 0.0 else 1.0)
        R[m] = old_area[m] * sc
        R[m] *= A_new / R[m].sum() if R[m].sum() > 0.0 else 0.0
    stretch = np.divide(R, old_area, out=np.ones(n), where=old_area > 0.0)
    K = coo_matrix((w1, (own1, new1)), shape=(n, n)).tocsr().multiply(W).tocsr()   # between neighbours or itself
    rs = np.asarray(K.sum(axis=1)).ravel()
    T = (diags(np.divide(R, rs, out=np.zeros(n), where=rs > 0.0)) @ K
         + diags(np.where(rs > 0.0, 0.0, R))).tocsr()          # rows to their stretched old areas
    length = np.zeros(0) if faces is None else faces.length
    return _fit_columns(T, area, fi, fj, length, x, comp), stretch


def _fit_columns(T, area, fi, fj, length, x, comp):
    """T (csr) with each column's sum brought to `area`, its rows' sums kept: the columns' shortfalls, r = area less the
    column sums, are met by a flow of area between the columns of neighbours across their faces, the potential flow of
    the faces' graph (conductance the face length over the markers' distance; one potential held per piece `comp`)
    with L p = -r, so that column i passes g = c (p_i - p_k) to column k; the columns are visited from the highest
    potential down, and what one passes on it takes from its own column's metal first (row i's entry T[i, i], moved to
    T[i, k]) and the rest from what it received, in proportion (row j's T[j, i] moved to T[j, k]): every move stays
    within its row and nothing turns negative, since a column passes on less than it holds (it keeps its area).
    Repeated on the round-off left, FIT_PASSES times at most; a fit still FIT_FAIL of an area short raises."""
    n = len(area)
    cols = [dict() for _ in range(n)]
    coo = T.tocoo()
    for j, i, v in zip(coo.row.tolist(), coo.col.tolist(), coo.data.tolist()):
        if v > 0.0:
            cols[i][j] = cols[i].get(j, 0.0) + v

    def shortfall():
        return area - np.array([sum(c.values()) for c in cols])
    if len(fi):
        g = length / np.linalg.norm(x[fi] - x[fj], axis=1)
        C = coo_matrix((g, (fi, fj)), shape=(n, n))
        C = (C + C.T).tocsr()
        L = (diags(np.asarray(C.sum(axis=1)).ravel()) - C).tocsr()
        free = np.ones(n, bool)
        free[np.unique(comp, return_index=True)[1]] = False
        Lf = L[free][:, free].tocsc()
    for _ in range(FIT_PASSES):
        r = shortfall()
        if np.all(np.abs(r) <= FIT_RTOL * area) or not len(fi):
            break
        p = np.zeros(n)
        if free.any():
            p[free] = spsolve(Lf, -r[free])
        out = [[] for _ in range(n)]
        for a, b, q in zip(fi.tolist(), fj.tolist(), (g * (p[fi] - p[fj])).tolist()):
            if q > 0.0:
                out[a].append((b, q))
            elif q < 0.0:
                out[b].append((a, -q))
        for i in np.argsort(-p, kind="stable").tolist():
            if not out[i]:
                continue
            col = cols[i]
            O = sum(q for _, q in out[i])
            own = col.get(i, 0.0)
            take = min(own, O)
            for k, q in out[i]:
                cols[k][i] = cols[k].get(i, 0.0) + q * (take / O)
            col[i] = own - take if take < own else 0.0
            others = [(h, v) for h, v in col.items() if h != i and v > 0.0]
            if O > take and others:
                frac = min((O - take) / sum(v for _, v in others), 1.0)
                for h, v in others:
                    moved = v * frac
                    col[h] = v - moved
                    for k, q in out[i]:
                        cols[k][h] = cols[k].get(h, 0.0) + moved * (q / O)
    r = shortfall()
    if np.any(np.abs(r) > FIT_FAIL * area):
        raise RuntimeError("the re-partition's columns did not meet the new areas: {:.1e} of an area left".format(
            float(np.max(np.abs(r) / area))))
    rows = [j for c in cols for j in c]
    at = [i for i, c in enumerate(cols) for _ in c]
    vals = [v for c in cols for v in c.values()]
    fitted = coo_matrix((vals, (rows, at)), shape=(n, n)).tocsr()
    fitted.eliminate_zeros()
    return fitted


def _tangent_basis(n):
    """Two unit vectors (e1, e2) spanning the plane normal to each n (n, 3)."""
    a = np.eye(3)[np.argmin(np.abs(n), axis=1)]
    e1 = np.cross(n, a)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    return e1, np.cross(n, e1)


def _clip(poly, lab, cnt, sp, j):
    """Clip the convex polygons poly (r, V, 2) with cnt vertices (edge k from vertex k to k+1 labelled lab[:, k]) to
    the side sp <= 0, sp (r, V) the vertices' signed values of the clipping half-plane a . y - c; the new edge is
    labelled j (r,). Returns the new (poly, lab, cnt), of width V + 1."""
    r, V, _ = poly.shape
    rr = np.arange(r)[:, None]
    k = np.arange(V)[None, :]
    real = k < cnt[:, None]
    nxt = np.where(k + 1 < cnt[:, None], k + 1, 0)
    q, sq = poly[rr, nxt], sp[rr, nxt]
    inp, inq = sp <= 0.0, sq <= 0.0
    cross = real & (inp != inq)
    t = np.divide(sp, sp - sq, out=np.zeros_like(sp), where=cross)
    I = poly + t[..., None] * (q - poly)
    pts = np.stack([poly, I], axis=2).reshape(r, 2 * V, 2)
    labs = np.stack([lab, np.where(inp, j[:, None], lab)], axis=2).reshape(r, 2 * V)
    keep = np.stack([real & inp, cross], axis=2).reshape(r, 2 * V)
    dest = np.cumsum(keep, axis=1) - 1
    out = np.zeros((r, V + 1, 2))
    out_lab = np.full((r, V + 1), -1, dtype=np.int64)
    rows = np.broadcast_to(rr, keep.shape)[keep]
    out[rows, dest[keep]] = pts[keep]
    out_lab[rows, dest[keep]] = labs[keep]
    return out, out_lab, keep.sum(axis=1)


class Markers:
    """n markers: positions `x` (n, 3) on phi, unit normals `normal` (n, 3) out of the metal, areas `area` (n,) [m2],
    phi's principal curvatures `kappa` (n, 2) [1/m], the row of each one's skin `skin_index` (n,) in `skins` (a
    `Skins`, or None for markers without skins), the patch of phi holding each `patch` (n,) (-1 for none), the target
    spacing `spacing` [m] and the counters `splits`, `merges`, `trapped`. After `tessellate`: `faces`, `capped` (n,)
    and `axis` (n, 3), the longest direction of each Voronoi cell. `x_tess` and `normal_tess` (n, 3) are each marker's
    position and normal at the last tessellation (at creation before any), the owners of the area it was given there
    (Task 9b, `transfer_matrix`); a split or a merge sets them to the new markers' own."""

    def __init__(self, x, normal, area, kappa, skin_index, patch, spacing, skins=None):
        self.x = np.array(x, dtype=float).reshape(-1, 3)
        n = len(self.x)
        self.normal = np.array(normal, dtype=float).reshape(n, 3)
        self.area = np.array(area, dtype=float).reshape(n)
        self.kappa = np.array(kappa, dtype=float).reshape(n, 2)
        self.skin_index = np.array(skin_index, dtype=np.int64).reshape(n)
        self.patch = np.array(patch, dtype=np.int64).reshape(n)
        self.spacing = float(spacing)
        self.skins = skins
        self.splits = self.merges = self.trapped = 0
        self.faces = None
        self.capped = np.zeros(n, bool)
        self.axis = np.zeros((n, 3))
        self.x_tess, self.normal_tess = self.x.copy(), self.normal.copy()

    @property
    def n(self):
        return len(self.x)

    @classmethod
    def create(cls, mesh, phi, select=None, material=None, settings=None, spacing=None):
        """One marker at the centroid of each selected patch of phi (`select`: a mask or indices over
        `cutcell.patches(mesh, phi)`, None for all), with its patch's area, the recovered normal and phi's principal
        curvatures there. With a `material`, each marker gets an empty skin (row i of a new `Skins` for marker i).
        The target spacing is `spacing` [m], or the square root of the selected patches' mean area."""
        P = cut_patches(mesh, phi)
        idx = np.arange(len(P.cell)) if select is None else np.asarray(select)
        if idx.dtype == bool:
            idx = np.flatnonzero(idx)
        x, area = P.centroid[idx], P.area[idx]
        _, normal = _field(mesh, phi, x)
        kappa = np.column_stack(principal_curvatures(mesh, phi, x)) if len(x) else np.zeros((0, 2))
        spacing = float(np.sqrt(area.mean())) if spacing is None else float(spacing)
        skins = Skins(material, area, kappa, settings) if material is not None else None
        return cls(x, normal, area, kappa, np.arange(len(x)), idx, spacing, skins)

    # -- motion -----------------------------------------------------------------------------------------------------
    def carry(self, velocity, dt):
        """x <- x + velocity(x) dt, velocity a function of the positions (n, 3) [m/s]."""
        self.x = self.x + np.asarray(velocity(self.x), dtype=float) * dt

    def snap(self, mesh, phi, reach=NEWTON_REACH):
        """Back onto phi's zero contour (module docstring): SNAP_PASSES steps along the recovered normal, then Newton
        steps along its own cell's gradient (`_cell_newton`) for each marker within `reach` [m] of the contour (np.inf:
        every marker); the normals are updated. Returns (the markers that took the Newton steps from farther than
        NEWTON_REACH -- the fallback, with a reach beyond it -- and the largest sideways slide [m] of the Newton steps, the
        displacement's part along the surface, normal to the new normal)."""
        for _ in range(SNAP_PASSES):
            val, n = _field(mesh, phi, self.x)
            self.x = self.x - val[:, None] * n
        val, _ = _field(mesh, phi, self.x)
        near = np.abs(val) <= reach
        before = self.x.copy()
        if near.any():
            self.x[near] = _cell_newton(mesh, phi, self.x[near])
        _, self.normal = _field(mesh, phi, self.x)
        d = self.x - before
        slide = np.linalg.norm(d - np.einsum("ij,ij->i", d, self.normal)[:, None] * self.normal, axis=1)
        return int((near & (np.abs(val) > NEWTON_REACH)).sum()), float(slide.max(initial=0.0))

    def base(self):
        """Each marker's skin base [m] (n, 3), x - H n with H its skin's thickness (Ruling 1)."""
        return self.x - self.skins.thickness()[self.skin_index][:, None] * self.normal

    # -- patches and the tessellation -------------------------------------------------------------------------------
    def assign(self, mesh, phi_patches):
        """`patch` <- the patch of `phi_patches` (a `cutcell.Patches`) whose cut cell contains each marker. A marker
        found in an uncut cell (on a face the contour meets, to round-off) takes the patch with the nearest centroid
        within that cell's size; -1 if there is none."""
        P = phi_patches
        of_cell = np.full(len(mesh.tets), -1, dtype=np.int64)
        of_cell[P.cell] = np.arange(len(P.cell))
        cell, _ = locate(mesh, self.x)
        patch = of_cell[cell]
        miss = np.flatnonzero(patch < 0)
        if len(miss) and len(P.cell):
            d, j = cKDTree(P.centroid).query(self.x[miss])
            patch[miss] = np.where(d <= np.sqrt(3.0) * mesh.size[cell[miss]], j, -1)
        self.patch = patch
        return patch

    def _neighbours(self, mesh, P):
        """The neighbour pairs (i < j) of Ruling 2 (module docstring)."""
        reach = NEIGHBOUR_REACH * self.spacing
        pairs = cKDTree(self.x).query_pairs(reach, output_type="ndarray")
        if not len(pairs):
            return pairs.reshape(0, 2)
        pi, pj = self.patch[pairs[:, 0]], self.patch[pairs[:, 1]]
        ok = (pi >= 0) & (pj >= 0)
        pairs, pi, pj = pairs[ok], pi[ok], pj[ok]
        limit = NEIGHBOUR_PATH * reach
        return pairs[patch_paths(mesh, P, pi, pj, limit) <= limit]

    def tessellate(self, mesh, phi):
        """Ruling 2's tessellation (module docstring): assigns the patches, updates the normals and curvatures, sets
        `area`, `capped` and `axis` and returns the faces as `Faces(i, j, length, unit_normal_ij)`, i < j."""
        P = cut_patches(mesh, phi)
        self.assign(mesh, P)
        n = self.n
        _, self.normal = _field(mesh, phi, self.x)
        self.kappa = np.column_stack(principal_curvatures(mesh, phi, self.x))
        pairs = self._neighbours(mesh, P)
        # each marker's neighbours, padded
        a = np.concatenate([pairs[:, 0], pairs[:, 1]])
        b = np.concatenate([pairs[:, 1], pairs[:, 0]])
        order = np.lexsort((np.linalg.norm(self.x[b] - self.x[a], axis=1), a))     # by marker, nearest first
        a, b = a[order], b[order]
        deg = np.bincount(a, minlength=n)
        K = int(deg.max()) if len(a) else 0
        start = np.concatenate([[0], np.cumsum(deg)[:-1]])
        slot = np.arange(len(a)) - start[a]
        nbr = np.full((n, K), -1, dtype=np.int64)
        nbr[a, slot] = b
        e1, e2 = _tangent_basis(self.normal)
        d3 = self.x[np.maximum(nbr, 0)] - self.x[:, None, :]
        D = np.stack([np.einsum("nkj,nj->nk", d3, e1), np.einsum("nkj,nj->nk", d3, e2)], axis=2)
        dd = (D * D).sum(axis=2)
        use = (nbr >= 0) & (dd > (1e-9 * self.spacing) ** 2)
        # the cap and the clipping
        r_in = 0.5 * NEIGHBOUR_REACH * self.spacing
        ang = 2.0 * np.pi * (np.arange(CAP_SIDES) + 0.5) / CAP_SIDES
        cap = (r_in / np.cos(np.pi / CAP_SIDES)) * np.column_stack([np.cos(ang), np.sin(ang)])
        poly = np.broadcast_to(cap, (n, CAP_SIDES, 2)).copy()
        lab = np.full((n, CAP_SIDES), -1, dtype=np.int64)
        cnt = np.full(n, CAP_SIDES, dtype=np.int64)
        for k in range(K):                                # nearest neighbours first: later ones rarely cut
            r = np.flatnonzero(use[:, k])
            if not len(r):
                continue
            V = int(cnt[r].max())
            sp = np.einsum("rvk,rk->rv", poly[r, :V], D[r, k]) - 0.5 * dd[r, k][:, None]
            cuts = ((sp > 0.0) & (np.arange(V)[None, :] < cnt[r][:, None])).any(axis=1)
            r, sp = r[cuts], sp[cuts]
            if not len(r):
                continue
            p2, l2, c2 = _clip(poly[r, :V], lab[r, :V], cnt[r], sp, np.full(len(r), k))
            if poly.shape[1] < V + 1:
                poly = np.pad(poly, ((0, 0), (0, V + 1 - poly.shape[1]), (0, 0)))
                lab = np.pad(lab, ((0, 0), (0, V + 1 - lab.shape[1])), constant_values=-1)
            poly[r, :V + 1], lab[r, :V + 1], cnt[r] = p2, l2, c2
        V = poly.shape[1]
        kk = np.arange(V)[None, :]
        real = kk < cnt[:, None]
        nxt = np.where(kk + 1 < cnt[:, None], kk + 1, 0)
        q = np.take_along_axis(poly, nxt[..., None], axis=1)
        cell_area = 0.5 * np.where(real, poly[..., 0] * q[..., 1] - q[..., 0] * poly[..., 1], 0.0).sum(axis=1)
        edge = np.linalg.norm(q - poly, axis=2)
        self.capped = (real & (lab < 0)).any(axis=1)
        # the longest direction: the principal axis of the cell's vertices
        cnt_safe = np.maximum(cnt, 1)
        cen = np.where(real[..., None], poly, 0.0).sum(axis=1) / cnt_safe[:, None]
        rel = np.where(real[..., None], poly - cen[:, None, :], 0.0)
        cov = np.einsum("nvi,nvj->nij", rel, rel)
        _, vec = np.linalg.eigh(cov)
        v = vec[:, :, -1]
        self.axis = v[:, :1] * e1 + v[:, 1:] * e2
        # faces, symmetrised
        on = real & (lab >= 0)
        rows = np.broadcast_to(np.arange(n)[:, None], on.shape)[on]
        cols = nbr[rows, lab[on]]
        L = coo_matrix((edge[on], (rows, cols)), shape=(n, n)).tocsr()
        S = ((L + L.T) * 0.5).tocoo()
        up = S.row < S.col
        fi, fj, fl = S.row[up].astype(np.int64), S.col[up].astype(np.int64), S.data[up]
        live = fl > 0.0
        fi, fj, fl = fi[live], fj[live], fl[live]

        def tangent(src, dst):
            d = self.x[dst] - self.x[src]
            nn = self.normal[src]
            t = d - np.einsum("ij,ij->i", d, nn)[:, None] * nn
            return t / np.linalg.norm(t, axis=1)[:, None]

        u = tangent(fi, fj) - tangent(fj, fi)
        u /= np.linalg.norm(u, axis=1)[:, None]
        # scale to the area of phi's patches the markers cover
        cell_area = np.where(self.patch >= 0, cell_area, 0.0)
        covered = P.area
        if len(P.cell) and n:
            dist, _ = cKDTree(self.x).query(P.centroid)
            covered = P.area[dist <= r_in / np.cos(np.pi / CAP_SIDES)]
        total = cell_area.sum()
        f = covered.sum() / total if total > 0.0 else 0.0
        self.area = cell_area * f
        self.faces = Faces(fi, fj, fl * np.sqrt(f), u)
        self.x_tess, self.normal_tess = self.x.copy(), self.normal.copy()
        return self.faces

    # -- topology ---------------------------------------------------------------------------------------------------
    def _append(self, src, x):
        """New markers copying the markers `src` at positions x, their skins' rows given later."""
        self.x = np.concatenate([self.x, x])
        for name in ("normal", "area", "kappa", "patch", "capped", "axis", "x_tess", "normal_tess"):
            setattr(self, name, np.concatenate([getattr(self, name), getattr(self, name)[src]]))

    def _drop(self, idx):
        keep = np.ones(self.n, bool)
        keep[idx] = False
        for name in ("x", "normal", "area", "kappa", "patch", "capped", "axis", "skin_index", "x_tess", "normal_tess"):
            setattr(self, name, getattr(self, name)[keep])
        self.faces = None

    def split(self, target_area):
        """Every marker above SPLIT_FACTOR times `target_area` [m2] becomes two, SPLIT_GAP target spacings apart along
        its cell's longest direction (from the last `tessellate`) and centred on it, each with half the area and a skin
        of the same profile with half the mass of every cell (`Skins.divide`). The new markers are appended; the caller
        snaps and tessellates. Returns the indices that split."""
        big = np.flatnonzero(self.area > SPLIT_FACTOR * target_area)
        if not len(big):
            return big
        off = 0.5 * SPLIT_GAP * np.sqrt(target_area) * self.axis[big]
        x = self.x[big]
        self.area[big] *= 0.5
        self._append(big, x - off)
        self.x[big] = x + off
        self.x_tess[big], self.x_tess[-len(big):] = self.x[big], self.x[-len(big):]   # each half owns its own side
        if self.skins is not None:
            new_rows = self.skins.divide(self.skin_index[big])
        else:
            new_rows = self.n - len(big) + np.arange(len(big))
        self.skin_index = np.concatenate([self.skin_index, new_rows])
        self.splits += len(big)
        self.faces = None
        return big

    def merge(self, target_spacing):
        """Markers closer than MERGE_FRACTION of `target_spacing` [m] merge in pairs, the closest first, each marker in
        at most one pair per call: the pair's first keeps the area-weighted position, normal and curvatures and the
        summed area; their skins are mixed by the hand-over rule (`Skins.combine`) and the second leaves. The caller
        snaps and tessellates. Returns the number of merges."""
        pairs = cKDTree(self.x).query_pairs(MERGE_FRACTION * target_spacing, output_type="ndarray")
        if not len(pairs):
            return 0
        dist = np.linalg.norm(self.x[pairs[:, 0]] - self.x[pairs[:, 1]], axis=1)
        taken = np.zeros(self.n, bool)
        keep, gone = [], []
        for p in np.argsort(dist, kind="stable"):
            i, j = pairs[p]
            if not (taken[i] or taken[j]):
                taken[i] = taken[j] = True
                keep.append(i)
                gone.append(j)
        keep, gone = np.array(keep), np.array(gone)
        wa, wb = self.area[keep], self.area[gone]
        A = wa + wb
        fa = np.where(A > 0.0, wa / np.where(A > 0.0, A, 1.0), 0.5)
        fb = 1.0 - fa
        self.x[keep] = fa[:, None] * self.x[keep] + fb[:, None] * self.x[gone]
        nrm = fa[:, None] * self.normal[keep] + fb[:, None] * self.normal[gone]
        self.normal[keep] = nrm / np.maximum(np.linalg.norm(nrm, axis=1), 1e-300)[:, None]
        self.kappa[keep] = fa[:, None] * self.kappa[keep] + fb[:, None] * self.kappa[gone]
        self.x_tess[keep], self.normal_tess[keep] = self.x[keep], self.normal[keep]
        self.area[keep] = A
        if self.skins is not None:
            self.skins.combine(self.skin_index[keep], self.skin_index[gone])
        self._remove(gone)
        self.merges += len(keep)
        return len(keep)

    def trapped_markers(self, mesh, phi, cell_size):
        """Indices of the markers with no zero contour of phi within TRAP_FRACTION of `cell_size` [m] (a number or one
        per marker) (section 9.6)."""
        cell_size = np.broadcast_to(np.asarray(cell_size, dtype=float), (self.n,))
        P = cut_patches(mesh, phi)
        if not len(P.tris):
            return np.arange(self.n)
        tri = P.tris
        cen = tri.mean(axis=1)
        r_max = np.linalg.norm(tri - cen[:, None, :], axis=2).max()
        tree = cKDTree(cen)
        lim = TRAP_FRACTION * cell_size
        kk = min(TRAP_CANDIDATES, len(tri))
        near, js = tree.query(self.x, k=kk)
        near, js = near.reshape(self.n, kk), js.reshape(self.n, kk)
        d = np.full(self.n, np.inf)
        for c in range(kk):
            t = tri[js[:, c]]
            d = np.minimum(d, point_triangle_distance(self.x, t[:, 0], t[:, 1], t[:, 2]))
        # where the nearest candidates are all farther than the limit but some triangle might still reach it
        for i in np.flatnonzero((d > lim) & (near[:, 0] <= lim + r_max)):
            ball = tree.query_ball_point(self.x[i], lim[i] + r_max)
            t = tri[ball]
            p = np.broadcast_to(self.x[i], (len(ball), 3))
            d[i] = min(d[i], point_triangle_distance(p, t[:, 0], t[:, 1], t[:, 2]).min())
        return np.flatnonzero(d > lim)

    def release(self, idx):
        """Remove the markers `idx` (distinct; e.g. the trapped ones) and their skins, counting them as trapped.
        Returns per removed marker, in the order of `idx`, (position (., 3), mass [kg], enthalpy [J], momentum
        [kg m/s] (., 3)) of what its skin held, for the caller to book as a deferred load in the cell holding it."""
        idx = np.asarray(idx, dtype=np.int64).reshape(-1)
        x = self.x[idx].copy()
        mass, H, P = self._remove(idx)
        self.trapped += len(idx)
        return x, mass, H, P

    def _remove(self, idx):
        """Drop markers and their skins' rows, renumbering the remaining rows. Returns the removed skins' totals."""
        idx = np.asarray(idx, dtype=np.int64).reshape(-1)
        if len(np.unique(idx)) != len(idx):
            raise ValueError("markers to remove must be distinct")
        rows = self.skin_index[idx]
        if self.skins is not None:
            out = self.skins.remove(rows)
        else:
            out = (np.zeros(len(idx)), np.zeros(len(idx)), np.zeros((len(idx), 3)))
        self._drop(idx)
        gone = np.zeros(self.n + len(idx), bool)
        gone[rows] = True
        shift = np.cumsum(gone)
        self.skin_index = self.skin_index - shift[self.skin_index]
        return out
