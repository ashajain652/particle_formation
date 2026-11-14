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
centroids; `patch_paths`) no longer than NEIGHBOUR_PATH reaches, less the pairs that face each other across gas: their
normals' dot product below OPPOSED_NORMALS (more than 120 degrees apart) and each in front of the other (n_i . d > 0 and
n_j . d < 0, d from i to j), as the two walls of a slot are. The path bound parts sheets that the surface joins only far
round, whatever lies between them; the normals part facing sheets however near the surface joins them. Neither parts
the two faces of a thin plate of metal (back to back, each behind the other) where the surface joins them within the
bound, nor a wall from the rim of the wall facing it, where the rim's normal blends with the face beyond it (10 short
faces in the slot below).

Task 9b. The rule before allowed at most two hops of the graph, but the patches are pieces of cut cells, on the body's
meshes about a sixth of the 2 mm spacing across, so two hops reached about one spacing. On the 10 mm test cube 123 of its
219 cells were capped and the areas were 0.38 to 2.67 times the area nearest each marker with nothing moving; at the
step ends 64 to 95 per cent of the area changes the tessellation itself made sat on markers whose neighbours had crossed
the hop limit. With the path, no cell is capped there and the areas are 0.68 to 1.10 times the nearest-marker areas (the
least at the rounded edges and corners, where a tangent plane stands in poorly for the surface); on a 10 mm sphere on the
design's 0.7 mm cells with markers at 2 mm, 0.98 to 1.02 and none capped. The bound is two reaches, not less, because
the centroids' path zigzags across sliver patches: with one marker per patch on a 25 mm sphere at 2 mm cells it is 1.04
times the straight distance at the median, 1.7 at the 95th percentile and 2.5 at the 99th, and at 1.5 reaches a rotation
of those markers lost 306 far pairs and changed areas by 8 %; at two it loses none. The path bound alone (fix round 1)
left 23 pairs of a slot's walls neighbours (the slot 1.5 mm wide, markers 1 mm apart, the walls joined round the floor
within the bound up to 2.7 mm above it), and they cut the lowest wall markers' areas to 0.137 of the nearest-marker
area; with the normals none is, and the least is 0.455, by the floor's concave corner. The threshold: walls facing each
other within 60 degrees of parallel (a V narrower than 60 degrees, a slot) are across gas, and a tangent plane cannot
stand for both; at a concave right angle the dot product is about zero (0.34 to 0.47 measured, the normals blending
within a cell of the corner) and the pairs stay. In the slot the facing walls' pairs lie below -0.9 and the pairs of a
wall with the facing rim between -0.9 and -0.5.

The plan named 1.5 spacings for the reach; with one marker per patch at creation (the patches' areas spread over three
decades) the true Voronoi cells reach past 0.75 spacings on most markers, so at 1.5 the cap bound 78 % of the cells on a
25 mm sphere at 2 mm cells, a capped cell could not grow with the surface (a uniform stretch to 26 mm gave 6 % less than
(26/25)^2) and a rotation changed areas by 73 % (7 % with the two hops). At 3.0 none is capped there (0.55 % at 0.7 mm
cells near the surface, one marker per patch), no pair within the reach is lost to the path bound, and the rotation
changes areas by 0.3 %.

In each marker's tangent plane the neighbours' offsets, projected, bound its Voronoi
cell (the half-planes y . d <= |d|^2 / 2), clipped from a CAP_SIDES-gon whose inscribed radius is half the reach (no
neighbour beyond the reach could bound the cell inside it); a cell that the cap still bounds is flagged `capped`. Each
edge is a face to the neighbour that made it. Faces are symmetrised: a pair's face has the mean of its two edge lengths
(an edge present on one side only counts as zero on the other) and the unit normal (t_ij - t_ji)/|t_ij - t_ji|, t_ij the
unit tangent-plane direction from i toward j in i's plane, so a face is shared by both markers with opposite normals.

The areas (stage 3 Task 11's fix round 1, C1, the controller's ruling; it replaces Ruling 2's scaled cells): each
marker's area is the area of phi it owns, the part of phi nearer it than any other marker (`ownership`, exact by planes
between the markers). This is the one partition of phi that the gas loads and the skins' base groups read, so a skin's
heat is its flux over the area it stands for and the beam's total is exact; with the tangent-plane cells' areas the q6p
box's markers took 0.961 of a beam's projected heat at the start and 0.640 after its first melting step, and its face's
central skins took the beam over cells 1.2 to 1.3 times the area beneath them. The tangent-plane cells now give only the
faces (the runoff's exchange, the modes' adjacency, the spray's regions, the re-partition's neighbours; their lengths
scaled by the square root of the factor that brings the cells' areas to phi's; the cells' areas so scaled are kept as
`cell_area`, for diagnosis), the `capped` flag and the cells' longest direction (`axis`, for splits). A marker without a
patch owns nothing (it is released at the step end).

Topology. `split` divides a marker above twice the target area into two half a target spacing apart along the longest
direction of its Voronoi cell, each with the same profile and half the mass and area; `merge` joins markers closer than
a third of the target spacing, their skins mixed by the hand-over rule (`Skins.combine`); `trapped_markers` finds the
markers with no zero contour of phi within half a cell, and `release` removes markers with their skins and returns what
the skins held, for the caller to book as a deferred load in the cell holding each. Every change of the markers changes
the skins' rows with them and is counted (`splits`, `merges`, `trapped`).

The re-partition of the skins' metal after a tessellation is `repartition`'s (Task 9b).

Numbers here are numerical settings of the construction; no physics."""
from collections import namedtuple

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import dijkstra
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
OPPOSED_NORMALS = -0.5     # two markers whose normals' dot product is below this, each in front of the other, face
                           # each other across gas and are never neighbours (Task 9b, fix round 1; see the docstring);
                           # and a marker whose normal's dot product with the surface's is below it never owns that
                           # surface (`ownership`'s same-side rule, design section 21.5)
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
OWN_SHARE = 1.0 / 64.0     # of the target spacing: a piece of phi at most this across is divided among its nearest
                           # markers whatever else lies near (`ownership`)
OWN_CANDIDATES = 4         # markers a piece of phi is divided among at most (`ownership`)
ALLOWED_FIRST = 16         # markers searched first for the nearest that may own a point under the same-side rule
                           # (`_nearest`; four times as many again where too few of them may)

Faces = namedtuple("Faces", "i j length normal")


def owners(x, points):
    """The marker that owns each point on phi (k, 3) [m]: the nearest of the markers at `x` (n, 3). Returns marker indices
    (k,)."""
    if not len(points):
        return np.zeros(0, dtype=np.int64)
    _, m = cKDTree(np.asarray(x, dtype=float)).query(np.asarray(points, dtype=float))
    return np.asarray(m, dtype=np.int64)


def _halfway(y, xu, xw):
    """|y - xu|^2 - |y - xw|^2 = 2 (xw - xu) . (y - (xu + xw) / 2) at the points y (k, V, 3) for the marker pairs xu, xw
    (k, 3) (k, V): negative where y is nearer xu, linear in y (taken from the pair's midpoint, against cancellation)."""
    return 2.0 * np.einsum("kvd,kd->kv", y - 0.5 * (xu + xw)[:, None, :], xw - xu)


def _clip_below(poly, cnt, f):
    """Convex polygons in space (k, V, 3) with `cnt` (k,) corners each, cut to where the linear function with the values
    `f` (k, V) at their corners is not positive. Returns (polygons (k, V + 1, 3), corners (k,))."""
    k, V = f.shape
    idx = np.arange(V)[None, :]
    real = idx < cnt[:, None]
    nxt = np.where(idx + 1 < cnt[:, None], idx + 1, 0)
    fn = np.take_along_axis(f, nxt, axis=1)
    pn = np.take_along_axis(poly, nxt[..., None], axis=1)
    keep = real & (f <= 0.0)
    cross = real & ((f <= 0.0) != (fn <= 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(cross, f / (f - fn), 0.0)
    hit = poly + s[..., None] * (pn - poly)
    cand = np.stack([poly, hit], axis=2).reshape(k, 2 * V, 3)
    use = np.stack([keep, cross], axis=2).reshape(k, 2 * V)
    order = np.argsort(~use, axis=1, kind="stable")[:, :V + 1]
    return np.take_along_axis(cand, order[..., None], axis=1), use.sum(axis=1)


def _polygon_area(poly, cnt):
    """Areas (k,) of convex polygons in space (k, V, 3) with `cnt` (k,) corners each (a fan from the first corner)."""
    out = np.zeros(len(poly))
    for v in range(1, poly.shape[1] - 1):
        tri = np.cross(poly[:, v] - poly[:, 0], poly[:, v + 1] - poly[:, 0])
        out += np.where(v + 1 < cnt, 0.5 * np.linalg.norm(tri, axis=1), 0.0)
    return out


def _below(f):
    """The share of a triangle on which a function linear over it, with the values `f` (k, 3) at its corners, is
    negative (k,): 1 or 0 with every corner on one side; with one corner v alone on its side, the corner of the triangle
    cut off at v, f_v^2 / ((f_v - f_a)(f_v - f_b)), or 1 less it."""
    neg = f < 0.0
    count = neg.sum(axis=1)
    alone = np.where(count == 1, np.argmax(neg, axis=1), np.argmax(~neg, axis=1))
    k = np.arange(len(f))
    fv, fa, fb = f[k, alone], f[k, (alone + 1) % 3], f[k, (alone + 2) % 3]
    with np.errstate(divide="ignore", invalid="ignore"):
        corner = np.clip(np.where(fv == 0.0, 0.0, fv * fv / ((fv - fa) * (fv - fb))), 0.0, 1.0)
    return np.select([count == 0, count == 3, count == 1], [0.0, 1.0, corner], 1.0 - corner)


def _divide(t, a, x, cand, live):
    """The triangles `t` (q, 3, 3) of areas `a` (q,) divided among their first `live` (q,) candidate markers `cand` (q,
    m): each one's share the triangle cut by the planes halfway between it and each other (for two, `_below`'s closed
    form). Returns the areas (q, m) [m2], zero past the live ones, each row summing to its area to round-off (the largest
    share takes what the others' rounding leaves)."""
    q, m = cand.shape
    out = np.zeros((q, m))
    for L in np.unique(live).tolist():
        sel = np.flatnonzero(live == L)
        if L == 1:
            out[sel, 0] = a[sel]
            continue
        if L == 2:
            share = _below(_halfway(t[sel], x[cand[sel, 0]], x[cand[sel, 1]]))
            out[sel, 0], out[sel, 1] = a[sel] * share, a[sel] - a[sel] * share
            continue
        for u in range(L):
            poly, cnt = t[sel], np.full(len(sel), 3)
            for w in range(L):
                if w != u:
                    poly, cnt = _clip_below(poly, cnt, _halfway(poly, x[cand[sel, u]], x[cand[sel, w]]))
            out[sel, u] = _polygon_area(poly, cnt)
        big = np.argmax(out[sel], axis=1)
        rest = out[sel].sum(axis=1) - out[sel, big]
        out[sel, big] = np.maximum(a[sel] - rest, 0.0)
    return out


def _dividing(t, x, cand, live, s):
    """The segments dividing the triangles `t` (q, 3, 3) among their first `live` (q,) candidate markers `cand` (q, m),
    as `_divide` divides them: for each pair (u, w) of a triangle's candidates, the triangle's cut by the plane halfway
    between them, kept where no other candidate is nearer (where u is nearer than each other v, the function `_halfway`
    of (u, v), linear along the cut, is not positive). A segment is measured on `s` (q, 3, 3), triangles matching `t`
    corner for corner, carried there by the affine map between the two (`s` = `t` measures it where it lies). The cut's
    ends follow `_below`'s sides (negative on u's), so a cut along an edge of two triangles counts in one of them; where a
    corner is exactly as near two markers, which of the triangles beside it is whole and which divided depends on the
    nearest-marker query's choice, and such a cut may count in neither or both (of measure zero away from symmetric
    set-ups). Returns (marker i, marker j, length [m]) of each piece, i < j."""
    out_i, out_j, out_l = [], [], []
    for L in np.unique(live).tolist():
        if L < 2:
            continue
        sel = np.flatnonzero(live == L)
        T, S, C = t[sel], s[sel], cand[sel]
        r = np.arange(len(sel))
        for u in range(L):
            for w in range(u + 1, L):
                xu, xw = x[C[:, u]], x[C[:, w]]
                f = _halfway(T, xu, xw)
                neg = f < 0.0
                count = neg.sum(axis=1)
                cut = (count == 1) | (count == 2)
                if not cut.any():
                    continue
                v = np.where(count == 1, np.argmax(neg, axis=1), np.argmax(~neg, axis=1))    # the corner alone
                va, vb = (v + 1) % 3, (v + 2) % 3
                ends = []
                for e in (va, vb):
                    with np.errstate(divide="ignore", invalid="ignore"):
                        tau = np.where(cut, f[r, v] / (f[r, v] - f[r, e]), 0.0)
                    lam = np.zeros((len(sel), 3))
                    lam[r, v] = 1.0 - tau
                    lam[r, e] += tau
                    ends.append(lam)
                lo, hi = np.zeros(len(sel)), np.ones(len(sel))
                p0, p1 = (np.einsum("ki,kid->kd", lam, T) for lam in ends)
                for o in range(L):
                    if o in (u, w):
                        continue
                    xo = x[C[:, o]]
                    g0, g1 = (_halfway(p[:, None, :], xu, xo)[:, 0] for p in (p0, p1))
                    with np.errstate(divide="ignore", invalid="ignore"):
                        at = g0 / (g0 - g1)
                    lo = np.where((g0 > 0.0) & (g1 <= 0.0), np.maximum(lo, at), lo)
                    hi = np.where((g0 <= 0.0) & (g1 > 0.0), np.minimum(hi, at), hi)
                    hi = np.where((g0 > 0.0) & (g1 > 0.0), 0.0, hi)
                q0, q1 = (np.einsum("ki,kid->kd", lam, S) for lam in ends)
                length = np.linalg.norm(q1 - q0, axis=1) * np.clip(hi - lo, 0.0, 1.0)
                keep = cut & (length > 0.0)
                out_i.append(np.minimum(C[keep, u], C[keep, w]))
                out_j.append(np.maximum(C[keep, u], C[keep, w]))
                out_l.append(length[keep])
    if not out_i:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0)
    return np.concatenate(out_i), np.concatenate(out_j), np.concatenate(out_l)


def _nearest(tree, points, k, normals=None, point_normal=None):
    """The k markers of `tree` (a cKDTree of the markers' positions) nearest each of `points` (N, 3): (distances (N, k)
    [m], indices (N, k)), nearest first. With `normals` (n, 3), the markers' normals, and `point_normal` (N, 3), the
    surface's normal at each point, only the markers that may own the point by the same-side rule, normals[m] .
    point_normal > OPPOSED_NORMALS (design section 21.5: ownership never crosses gas; a normal shorter than one loosens
    the rule, and one shorter than -OPPOSED_NORMALS lifts it), searched among the ALLOWED_FIRST nearest and then four
    times as many again until k are found or every marker is searched; a row with fewer is padded with distance inf and
    index n."""
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    if normals is None:
        d, i = tree.query(pts, k=k)
        return np.asarray(d, dtype=float).reshape(len(pts), -1), np.asarray(i, dtype=np.int64).reshape(len(pts), -1)
    n, N = tree.n, len(pts)
    pn = np.asarray(point_normal, dtype=float).reshape(-1, 3)
    dist, idx = np.full((N, k), np.inf), np.full((N, k), n, dtype=np.int64)
    todo, kq = np.arange(N), min(n, max(k, ALLOWED_FIRST))
    while len(todo):
        d, i = tree.query(pts[todo], k=kq)
        d, i = np.asarray(d, dtype=float).reshape(len(todo), kq), np.asarray(i, dtype=np.int64).reshape(len(todo), kq)
        ok = np.einsum("tqd,td->tq", normals[i], pn[todo]) > OPPOSED_NORMALS
        order = np.argsort(~ok, axis=1, kind="stable")[:, :k]                # the allowed first, nearest first
        kk = order.shape[1]
        allowed = np.take_along_axis(ok, order, axis=1)
        dist[todo, :kk] = np.where(allowed, np.take_along_axis(d, order, axis=1), np.inf)
        idx[todo, :kk] = np.where(allowed, np.take_along_axis(i, order, axis=1), n)
        todo = todo[(ok.sum(axis=1) < k) & (kq < n)]
        kq = min(n, 4 * kq)
    return dist, idx


def ownership(x, tris, tri_patch, n_patch, size, normals=None, tri_normal=None, boundaries=False, measure=None):
    """The one partition of phi by the markers at `x` (n, 3) (stage 3 Task 11's fix round 1, C1): every point of phi
    belongs to its nearest marker, so that a marker owns the part of phi inside its Voronoi cell in space. Each marker's
    area is the area it owns (`Markers.tessellate`), the gas loads give it the mean over what it owns
    (`loads.by_patches`), and its skin's base group is the base beneath it (`skinbody.LevelSetSkinBody._beta_owners`).

    Evaluated on the patches' triangles `tris` (Q, 3, 3) [m] (of patch `tri_patch` (Q,) of `n_patch`), exactly: a
    triangle whose three corners have one nearest marker lies wholly in that marker's (convex) cell and is its whole. For
    any other, every point y within r of its centroid c has its nearest marker within d0 + 2 r of c, d0 the distance from
    c to the marker nearest it; when at most OWN_CANDIDATES markers lie that near, the triangle is divided among them by
    the planes halfway between each pair (`_divide`), exactly. Any other is cut into four by its edge midpoints, again and
    again; a piece at most `size` [m] across (the square root of twice its area) is divided among the OWN_CANDIDATES
    markers nearest its centroid whatever else lies near, which could leave it inexact by a part of it another marker's
    cell reaches into. Returns W (n, n_patch), csr [m2]: W[i, p] the area of patch p owned by marker i; its columns sum
    to the patches' areas (the sums of their triangles' areas) to round-off.

    The same-side rule (design section 21.5; plan 03b Task 1): with the markers' `normals` (n, 3) and each triangle's
    unit normal `tri_normal` (Q, 3), marker m may own triangle q only if normals[m] . tri_normal[q] > OPPOSED_NORMALS,
    so that ownership never crosses gas (the facing walls of a slot narrower than the spacing); "nearest" above is then
    nearest among the markers allowed (`_nearest`), and the argument for exactness stands with them. A triangle no
    marker may own is owned by none. With `boundaries`, also the segments dividing the triangles (`_dividing`), measured
    on `measure` (Q, 3, 3), triangles matching `tris` corner for corner (by default `tris` themselves; the division of
    the layer measures on the level surface what it divides on phi, `division.divide`): returns (W, (i, j, length)),
    each pair of markers i < j sharing a boundary once, with its segments' lengths summed [m]. Without normals and
    boundaries the result is the one before these options, bit for bit."""
    x = np.asarray(x, dtype=float).reshape(-1, 3)
    n = len(x)
    t = np.asarray(tris, dtype=float).reshape(-1, 3, 3)
    tp = np.asarray(tri_patch, dtype=np.int64)
    if not n or not len(t):
        W = csr_matrix((n, n_patch))
        return (W, (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0))) if boundaries else W
    tree = cKDTree(x)
    m = min(OWN_CANDIDATES, n)
    a = 0.5 * np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1)
    q = np.arange(len(t))                                       # each piece's triangle, for its normal
    nn = None if normals is None else np.asarray(normals, dtype=float).reshape(-1, 3)
    tn = None if normals is None else np.asarray(tri_normal, dtype=float).reshape(-1, 3)
    s = (t if measure is None else np.asarray(measure, dtype=float).reshape(-1, 3, 3)) if boundaries else None
    near = lambda points, k, tri: _nearest(tree, points, k, nn, None if tn is None else tn[tri])
    _, o = near(t.reshape(-1, 3), 1, np.repeat(q, 3))
    o = o.reshape(-1, 3)
    if nn is not None:
        mine = o[:, 0] < n                                      # a triangle no marker may own stays unowned
        t, tp, o, a, q = t[mine], tp[mine], o[mine], a[mine], q[mine]
        s = None if s is None else s[mine]
    rows, cols, vals = [], [], []
    pairs = []
    while len(t):
        same = (o[:, 0] == o[:, 1]) & (o[:, 0] == o[:, 2])
        rows.append(o[same, 0]), cols.append(tp[same]), vals.append(a[same])
        t, tp, o, a, q = t[~same], tp[~same], o[~same], a[~same], q[~same]
        s = None if s is None else s[~same]
        if not len(t):
            break
        c = t.mean(axis=1)
        r = np.linalg.norm(t - c[:, None, :], axis=2).max(axis=1)
        d, k = near(c, min(m + 1, n) if nn is None else m + 1, q)
        reach = d[:, 0] + 2.0 * r
        live = (d[:, :m] <= reach[:, None]).sum(axis=1)       # the markers within reach, nearest first
        last = np.sqrt(2.0 * a) <= size
        done = last | (d.shape[1] == m) | (d[:, -1] > reach)
        if done.any():
            divided = np.where(last[done], np.isfinite(d[done, :m]).sum(axis=1), live[done])
            share = _divide(t[done], a[done], x, k[done, :m], divided)
            keep = share > 0.0
            rows.append(k[done, :m][keep]), cols.append(np.broadcast_to(tp[done][:, None], share.shape)[keep])
            vals.append(share[keep])
            if boundaries:
                pairs.append(_dividing(t[done], x, k[done, :m], divided, s[done]))
        t, tp, o, a, q = t[~done], tp[~done], o[~done], a[~done], q[~done]
        s = None if s is None else s[~done]
        if not len(t):
            break
        A, B, C = t[:, 0], t[:, 1], t[:, 2]
        AB, BC, CA = 0.5 * (A + B), 0.5 * (B + C), 0.5 * (C + A)
        _, om = near(np.concatenate([AB, BC, CA]), 1, np.tile(q, 3))
        oAB, oBC, oCA = om.reshape(3, -1)
        oA, oB, oC = o[:, 0], o[:, 1], o[:, 2]
        t = np.concatenate([np.stack(v, axis=1) for v in ((A, AB, CA), (AB, B, BC), (CA, BC, C), (AB, BC, CA))])
        o = np.concatenate([np.stack(v, axis=1) for v in ((oA, oAB, oCA), (oAB, oB, oBC), (oCA, oBC, oC),
                                                          (oAB, oBC, oCA))])
        tp, a, q = np.tile(tp, 4), np.tile(0.25 * a, 4), np.tile(q, 4)
        if s is not None:
            A, B, C = s[:, 0], s[:, 1], s[:, 2]
            AB, BC, CA = 0.5 * (A + B), 0.5 * (B + C), 0.5 * (C + A)
            s = np.concatenate([np.stack(v, axis=1) for v in ((A, AB, CA), (AB, B, BC), (CA, BC, C), (AB, BC, CA))])
    rows = np.concatenate(rows) if rows else np.zeros(0, dtype=np.int64)
    W = coo_matrix((np.concatenate(vals) if vals else np.zeros(0), (rows, np.concatenate(cols) if cols else rows)),
                   shape=(n, n_patch)).tocsr()
    W.sum_duplicates()
    if not boundaries:
        return W
    if not pairs:
        return W, (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0))
    i, j, length = (np.concatenate(v) for v in zip(*pairs))
    key, at = np.unique(i * n + j, return_inverse=True)
    total = np.bincount(at.ravel(), length, minlength=len(key))
    return W, (key // n, key % n, total)


def owned(x, P, size):
    """`ownership` of phi's patches `P` (`cutcell.patches`) by the markers at `x` (n, 3), with pieces at most `size` [m]
    across."""
    return ownership(x, P.tris, P.tri_patch, len(P.area), size)


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
    (Task 9b, `repartition.transfer_matrix`); a split or a merge sets them to the new markers' own. `top` (n,) bool
    each marker's mode (design 20.6: False a full-depth skin, True a top-zone skin) and `hold` (n,) the macro steps in a
    row its thin test has given the other answer (`modes`, Task 10); a split's new half takes its parent's mode and
    hold (the same surface: a new half alone could never make the group a change needs), a merge across modes the
    top-zone marker's (the step end hands the volume what the mixed column holds below its fine zone) and otherwise the
    first's, and a release drops them."""

    def __init__(self, x, normal, area, kappa, skin_index, patch, spacing, skins=None):
        self.x = np.array(x, dtype=float).reshape(-1, 3)
        n = len(self.x)
        self.normal = np.array(normal, dtype=float).reshape(n, 3)
        self.area = np.array(area, dtype=float).reshape(n)
        self.cell_area = self.area.copy()
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
        self.top, self.hold = np.zeros(n, bool), np.zeros(n, np.int64)

    def __setstate__(self, state):
        """Unpickling: a state from before Task 9b has no `x_tess`, `normal_tess`; saved between steps, the markers'
        last tessellation was where they are, so they are filled from `x` and `normal`. One from before Task 10 has no
        `top`, `hold`: every skin was full-depth, with no hold counted. One from before Task 11's fix round 1 has no
        `cell_area`: its areas were its tangent-plane cells'."""
        self.__dict__.update(state)
        if "x_tess" not in state:
            self.x_tess = np.array(self.x, dtype=float)
        if "normal_tess" not in state:
            self.normal_tess = np.array(self.normal, dtype=float)
        if "top" not in state:
            self.top = np.zeros(len(self.x), bool)
        if "hold" not in state:
            self.hold = np.zeros(len(self.x), np.int64)
        if "cell_area" not in state:
            self.cell_area = np.array(self.area, dtype=float)

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
        i, j = pairs[:, 0], pairs[:, 1]
        d = self.x[j] - self.x[i]
        ni, nj = self.normal[i], self.normal[j]
        facing = ((np.einsum("ij,ij->i", ni, nj) < OPPOSED_NORMALS) & (np.einsum("ij,ij->i", ni, d) > 0.0)
                  & (np.einsum("ij,ij->i", nj, d) < 0.0))       # across a gap of gas: each in front of the other
        pi, pj = self.patch[i], self.patch[j]
        ok = (pi >= 0) & (pj >= 0) & ~facing
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
        # the areas: the area of phi each marker with a patch owns (`ownership`, the one partition of phi; stage 3 Task
        # 11's fix round 1, C1), a marker without one none; the faces' lengths scaled as the tangent-plane cells would be
        # to phi's area
        cell_area = np.where(self.patch >= 0, cell_area, 0.0)
        total = cell_area.sum()
        f = P.area.sum() / total if total > 0.0 else 0.0
        self.cell_area = cell_area * f                    # the faces' cells, kept for diagnosis (the neighbour tests)
        live = np.flatnonzero(self.patch >= 0)
        self.area = np.zeros(n)
        if len(P.cell) and len(live):
            self.area[live] = np.asarray(owned(self.x[live], P, OWN_SHARE * self.spacing).sum(axis=1)).ravel()
        self.faces = Faces(fi, fj, fl * np.sqrt(f), u)
        self.x_tess, self.normal_tess = self.x.copy(), self.normal.copy()
        return self.faces

    # -- topology ---------------------------------------------------------------------------------------------------
    def _append(self, src, x):
        """New markers copying the markers `src` at positions x, their skins' rows given later."""
        self.x = np.concatenate([self.x, x])
        for name in ("normal", "area", "cell_area", "kappa", "patch", "capped", "axis", "x_tess", "normal_tess", "top",
                     "hold"):
            setattr(self, name, np.concatenate([getattr(self, name), getattr(self, name)[src]]))

    def _drop(self, idx):
        keep = np.ones(self.n, bool)
        keep[idx] = False
        for name in ("x", "normal", "area", "cell_area", "kappa", "patch", "capped", "axis", "skin_index", "x_tess",
                     "normal_tess", "top", "hold"):
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
        snaps and tessellates. Where the two differ in mode the merged marker is top-zone with the top-zone marker's hold
        (`Markers.top`, `Markers.hold`). Returns the number of merges."""
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
        across = self.top[gone] & ~self.top[keep]        # across modes the top-zone marker's mode and hold (Task 10)
        self.top[keep], self.hold[keep] = self.top[keep] | across, np.where(across, self.hold[gone], self.hold[keep])
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
