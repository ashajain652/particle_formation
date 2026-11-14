"""Operations on a P1 level set (negative inside the metal) on an AdaptiveMesh: reset to a signed distance, motion
under a prescribed velocity, connected pieces, the exact volume correction per piece (design §6.2-§6.4), and the move of
the gas surface phi by each marker's exact volume (design §21.4 and §21.6: phi is the evolved level set; plan 03b Task 2,
planning ruling 2; `move_by_volumes`).

How phi moves by volumes. Each node of phi's band takes the marker owning its foot, node - phi g (g phi's recovered
gradient), by `markers.ownership`'s rules: nearest among the markers on its own piece of metal whose normals do not
oppose g (the same-side rule, n . g > `markers.OPPOSED_NORMALS`; the same-piece rule, the controller's ruling T1-i), as
`division.divide` assigns the layer. A marker's group is its nodes, offset by one distance d_i (outward positive). The
group's swept volume is the change of the metal's volume over the part of the mesh its marker owns: each cell near phi
is divided exactly among the markers owning its points' feet, the feet taken by the foot map interpolated linearly from
the cell's corners' feet (so a cell whose four corners' feet one marker owns is wholly that marker's, and any other is cut
by the planes on which two candidates' feet are equally near). Every cell whose volume can change is so divided, so the
groups' swept volumes sum to the change of phi's volume; only cells no marker may own (none of their piece, or all
opposed) are counted in no group. (A cell given whole to the marker owning its centroid's foot, the obvious division,
gave the stage 3 box's groups, in a recession of 1.0 mm, 0.79 to 1.14 of the volume `division.divide` gives their
columns; the exact division gives them within the division's own level sampling, about +-0.8 %, measured once.)

The offsets are found by Newton's method on the groups' swept volumes, with their exact first-order Jacobian (each
owned part's contour, weighted by its cell's hat functions over |grad phi|) and steps halved until the misses fall
(planning ruling 2 named a secant iteration on each group alone: on the stage 3 box, its groups coupled through the
cells their boundaries cross, it left misses of 1e-2 after 8 rounds, and Newton from the offsets d_i = dV_i / A_i(0)
5e-6; from each group's own volume under offsets alike, the first offsets here, 9e-7 after 4, measured once). The
iterations move the unreset level set; the reset follows, once, keeping the contour exactly (`redistance(keep=True)`):
the plain reset moves the contour wherever it bends within a cell, on the stage 3 box by 0.26 mm3 at the first reset
of its exact signed distance and by a little less at each reset after it (0.21, 0.16, 0.14 mm3), always inward at its
convex edges, which would round them and spoil each group's exact volume."""
from collections import namedtuple

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree

from .criteria import nodal_average
from .cutcell import cell_gradients, inside_volumes, patches
from .transfer import interpolate


def _segment_distance(p, a, b):
    ab = b - a
    t = np.einsum("ij,ij->i", p - a, ab) / np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300)
    return np.linalg.norm(p - (a + np.clip(t, 0.0, 1.0)[:, None] * ab), axis=1)


def point_triangle_distance(p, a, b, c):
    """Distance from points p (n,3) to triangles (a, b, c) (n,3 each). The nearest of the three edges, or the
    perpendicular distance where the foot falls inside the triangle; degenerate (zero-area) triangles reduce to
    their edges, so a contour passing through a node gives no NaN."""
    d = np.minimum(np.minimum(_segment_distance(p, a, b), _segment_distance(p, b, c)), _segment_distance(p, c, a))
    n = np.cross(b - a, c - a)
    nn = np.einsum("ij,ij->i", n, n)
    ok = nn > 1e-30 * np.maximum(np.einsum("ij,ij->i", b - a, b - a), 1e-300) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.einsum("ij,ij->i", p - a, n) / nn
        foot = p - h[:, None] * n
        u = np.einsum("ij,ij->i", np.cross(b - foot, c - foot), n) / nn
        v = np.einsum("ij,ij->i", np.cross(c - foot, a - foot), n) / nn
    inside = ok & (u >= 0) & (v >= 0) & (u + v <= 1)
    d[inside] = np.minimum(d[inside], np.abs(h[inside]) * np.sqrt(nn[inside]))
    return d


def redistance(mesh, phi, k=8, band=None, keep=False):
    """Signed distance to the zero contour of phi, keeping phi's sign. Nodal values are exact distances to the
    piecewise-flat contour, searched among the k triangles with the nearest centroids. With `band` (m), only nodes
    within band of a triangle are computed (exactly, by the search below); the rest get +-band, never more than their
    true distance.

    The P1 interpolant of those distances does not in general have the contour as its zero: where the contour bends
    within a cell (at an edge, or between two moves' offsets) it moves, and at a convex edge inward (module docstring).
    With `keep`, the nodes of the cells the contour cuts keep phi's values, so the contour stays exactly where it was and
    phi is a distance at every other node (`move_by_volumes`; plan 03b Task 2)."""
    P = patches(mesh, phi)
    if len(P.tris) == 0:
        return np.where(phi < 0, -np.inf, np.inf)
    tri = P.tris
    kd = cKDTree(tri.mean(axis=1))
    kk = min(k, len(tri))
    d = np.full(len(mesh.points), np.inf if band is None else float(band))
    if band is None:
        sel = np.arange(len(mesh.points))
        _, idx = kd.query(mesh.points, k=kk)
    else:
        r_max = np.linalg.norm(tri - tri.mean(axis=1, keepdims=True), axis=2).max()   # centroid to farthest vertex
        near, _ = kd.query(mesh.points, k=1, distance_upper_bound=band + r_max)
        sel = np.flatnonzero(np.isfinite(near))
        _, idx = kd.query(mesh.points[sel], k=kk)
    idx = idx.reshape(len(sel), kk)
    ds = np.full(len(sel), np.inf)
    for j in range(kk):
        t = tri[idx[:, j]]
        ds = np.minimum(ds, point_triangle_distance(mesh.points[sel], t[:, 0], t[:, 1], t[:, 2]))
    d[sel] = ds if band is None else np.minimum(ds, band)
    out = np.where(phi < 0, -d, d)
    if keep:
        neg = (phi[mesh.tets] < 0).sum(axis=1)
        held = np.unique(mesh.tets[(neg > 0) & (neg < 4)])
        out[held] = phi[held]
    return out


def advect(mesh, phi, velocity, dt, band=None):
    """Semi-Lagrangian step: phi_new(x) = phi(x - u(x) dt), then reset to a distance."""
    x = mesh.points
    return redistance(mesh, interpolate(mesh, phi, x - velocity(x) * dt), band=band)


def _face_pairs(mesh, cells):
    T = mesh.tets[cells]
    f = np.vstack([T[:, [0, 1, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 3]], T[:, [1, 2, 3]]])
    f.sort(axis=1)
    owner = np.tile(np.arange(len(cells)), 4)
    _, inv = np.unique(f, axis=0, return_inverse=True)
    inv = inv.ravel()
    order = np.argsort(inv, kind="stable")
    s_inv, s_own = inv[order], owner[order]
    same = s_inv[1:] == s_inv[:-1]
    return s_own[:-1][same], s_own[1:][same]


def pieces(mesh, phi):
    """Connected pieces of metal (cells holding any metal, joined through shared faces).
    Returns (label per cell, -1 where no metal; number of pieces), pieces ordered by decreasing volume."""
    vol = inside_volumes(mesh, phi)
    cells = np.flatnonzero(vol > 0)
    label = np.full(len(mesh.tets), -1)
    if len(cells) == 0:
        return label, 0
    a, b = _face_pairs(mesh, cells)
    g = coo_matrix((np.ones(len(a)), (a, b)), shape=(len(cells), len(cells)))
    n, lab = connected_components(g, directed=False)
    pv = np.bincount(lab, vol[cells], minlength=n)
    rank = np.empty(n, int)
    rank[np.argsort(-pv, kind="stable")] = np.arange(n)
    label[cells] = rank[lab]
    return label, n


def piece_volumes(mesh, phi, label, n):
    vol = inside_volumes(mesh, phi)
    m = label >= 0
    return np.bincount(label[m], vol[m], minlength=n)


SWEEPS_MAX = 20     # sweeps over the pieces at most in `correct_volumes`
BRENT_XTOL = 1e-18  # m: Brent's method's tolerance on a piece's shift (`correct_volumes`), far below any cell
BRENT_ITER = 200    # its iterations at most


def correct_volumes(mesh, phi, label, targets, rings=2, rtol=1e-14, max_iter=60, misses=None):
    """Shift phi uniformly over each piece (its cells' nodes, dilated by `rings` rings of cells) so that the
    piece's volume equals targets[p] to round-off. Returns (phi, largest relative correction); with `misses` a list,
    each piece's miss before the correction, as a share of its target, is appended to it in label order (stage 3 Task
    11's fix round 2: the body's, piece 0's, beside the worst).

    Pieces whose dilations meet (Task 9, fix round 4): each node is shifted with the first piece (in label order) that
    reaches it, and each cell holding no metal before is counted with the first piece that reaches it, so no metal is
    counted for two pieces; the pieces are corrected in turn, sweep after sweep (at most SWEEPS_MAX), until every piece
    holds its target at once, since a later piece's shift moves the contour in a cell an earlier piece counts. (Counted
    for both, the metal a shift grew in a cell between two pieces 0.25 to 0.4 mm apart on 0.125 mm cells was booked
    twice and the total fell short by up to 3e-5; on the full-size cube's second step, with a 2e-9 m3 piece beside the
    body, a short 1e-14 m3 a move made the step's mass residual -9.4e-9.) One piece is corrected as before. A piece
    that does not converge raises with its target, the volume it held, its last shift and bracket, the iterations, its
    claimed nodes and cells, and the miss before the correction (stage 3 Task 11's review, M6).

    Each piece's shift is found by Newton steps on the contour's area, kept within the bracket of shifts that held too
    little and too much: bisection of the bracket where a Newton step leaves it, and Brent's method on it where a Newton
    step fails to halve the step before the last (rtsafe's test; stage 3 Task 11's fix round 1). A speck a base move makes beside the body,
    asked to grow some 1e5 times, first steps far beyond the nodes it may shift; there its volume no longer responds
    while the cut cells at its claim's edge keep a patch area, and the bare Newton steps crawled for all 60 iterations
    (the q6p box's step to 1.5 s and the shoulder's to 2.0 s stopped so)."""
    phi = phi.copy()
    targets = np.asarray(targets, float)
    if (targets <= 0).any():
        raise ValueError("volume target of piece {} is not positive; remove the piece instead".format(
            int(np.flatnonzero(targets <= 0)[0])))
    claimed_nodes = np.zeros(len(mesh.points), bool)
    claimed_cells = np.zeros(len(mesh.tets), bool)
    sets = []
    for p in range(len(targets)):
        others = np.zeros(len(mesh.points), bool)
        others[mesh.tets[(label >= 0) & (label != p)].ravel()] = True       # nodes of every other piece
        nodes = np.unique(mesh.tets[label == p])
        for _ in range(rings):
            touch = np.isin(mesh.tets, nodes).any(axis=1)
            nodes = np.unique(mesh.tets[touch])
        nodes = nodes[~others[nodes] & ~claimed_nodes[nodes]]
        claimed_nodes[nodes] = True
        cells = np.flatnonzero(np.isin(mesh.tets, nodes).any(axis=1) & ((label < 0) | (label == p)) & ~claimed_cells)
        claimed_cells[cells] = True
        sets.append((nodes, cells, phi[nodes].copy()))
    shift = np.zeros(len(targets))
    worst = 0.0
    for sweep in range(SWEEPS_MAX):
        for p, target in enumerate(targets):
            nodes, cells, base = sets[p]
            lo, hi, c = None, None, shift[p]
            step_before = step = np.inf                    # the last two steps (the bisection safeguard below)
            for _ in range(max_iter):
                phi[nodes] = base - c
                V, c_held = inside_volumes(mesh, phi, cells).sum(), c
                err = target - V
                if sweep == 0 and c == 0.0:
                    worst = max(worst, abs(err) / target)
                    if misses is not None:
                        misses.append(abs(err) / target)
                if abs(err) <= rtol * target:
                    break
                if err > 0:
                    lo = c
                else:
                    hi = c
                A = patches_area(mesh, phi, cells)
                c_new = c + err / max(A, 1e-300)
                bracketed = lo is not None and hi is not None
                if (lo is not None and c_new <= lo) or (hi is not None and c_new >= hi):
                    c_new = 0.5 * (lo + hi) if bracketed else c + 2 * (c_new - c)
                elif bracketed and abs(c_new - c) > 0.5 * abs(step_before):
                    # Newton not halving the step before the last: the contour's area is not the volume's slope here
                    # (cut cells at the claim's edge keep an area that does not respond), so Brent's method on the bracket
                    def held(x, nodes=nodes, base=base, cells=cells, target=target):
                        phi[nodes] = base - x
                        return target - inside_volumes(mesh, phi, cells).sum()
                    c_new = brentq(held, lo, hi, xtol=BRENT_XTOL, rtol=4.0 * np.finfo(float).eps, maxiter=BRENT_ITER)
                step_before, step = step, c_new - c
                c = c_new
            else:
                raise RuntimeError("volume correction of piece {} did not converge: target {:.6e} m3, held {:.6e} m3 "
                                   "at the shift {:.6e} m, bracket [{}, {}] m, {} iterations (sweep {}), {} nodes "
                                   "and {} cells claimed, the largest miss before the correction {:.3e}".format(
                                       p, target, V, c_held, "none" if lo is None else "{:.6e}".format(lo),
                                       "none" if hi is None else "{:.6e}".format(hi), max_iter, sweep, len(nodes),
                                       len(cells), worst))
            shift[p] = c
        if len(targets) == 1 or all(abs(t - inside_volumes(mesh, phi, sets[p][1]).sum()) <= rtol * t
                                    for p, t in enumerate(targets)):
            return phi, worst
    raise RuntimeError("the volume correction did not hold every piece's target at once in {} sweeps".format(SWEEPS_MAX))


def patches_area(mesh, phi, cells):
    return float(patches(mesh, phi, cells).area.sum())


MOVE_ROUNDS = 8            # Newton iterations at most in each round of a move (`move_by_volumes`; planning ruling 2)
MOVE_RTOL = 1e-6           # each group's swept volume within this share of its change (ruling 2)
MOVE_BAND_SHARE = 0.5      # a move whose first estimate is longer than this share of the reset's band is taken in
                           # rounds no longer (ruling 2, as stepend._move_beta takes its moves; this module may not
                           # import stepend, so it keeps its own number)
MOVE_FLOOR = 1e-3          # of the largest change of a move: the least denominator of a group's miss, so that a group
                           # asked for nothing is judged against the move's scale
MOVE_REACH = 4.0           # times the longest first estimate (and within the band): the furthest a group's offset goes
                           # in a round, and the depth of phi's band whose cells are divided among the groups
MOVE_PROFILE = 8           # offsets alike, from none to the reach, at which each group's swept volume is sampled for
                           # its first offset
MOVE_HALVINGS = 6          # halvings of a Newton step at most before the iteration stops
MOVE_CANDIDATES = 8        # markers a cell's feet are divided among at most

Moved = namedtuple("Moved", "phi swept miss rounds")
Moved.__doc__ = """A move of phi by volumes: phi moved, reset to a distance with its contour kept (`redistance(keep=True)`);
each group's swept volume (n,) [m3], outward positive; its miss (n,), |swept - dV| / max(|dV|, MOVE_FLOOR times the
largest |dV|); and the rounds it was taken in (zero when nothing moved)."""

_REF = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])   # the reference tetrahedron
_REST = np.array([[j for j in range(4) if j != lone] for lone in range(4)])             # the corners but one


def _cut(p, q, fp, fq):
    return p + (fp / (fp - fq))[:, None] * (q - p)


def _tvol(a, b, c, d):
    return np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0


def _fraction(f):
    """The share of a tetrahedron on which the linear function with the values f (m, 4) at its corners is negative (m,),
    a corner at exactly zero outside (as `cutcell.inside_volumes`, whose construction on the reference tetrahedron this
    is: the share is the same on any tetrahedron)."""
    neg = f < 0
    nneg = neg.sum(axis=1)
    out = np.where(nneg == 4, 1.0, 0.0)
    for want, flip in ((1, False), (3, True)):
        m = np.flatnonzero(nneg == want)
        if not len(m):
            continue
        ff = f[m]
        lone = np.argmax(~neg[m] if flip else neg[m], axis=1)
        fa = ff[np.arange(len(m)), lone]
        corner = np.ones(len(m))
        for j in range(4):
            with np.errstate(divide="ignore", invalid="ignore"):
                corner = np.where(j != lone, corner * fa / (fa - ff[:, j]), corner)
        out[m] = 1.0 - corner if flip else corner
    m = np.flatnonzero(nneg == 2)
    if len(m):
        order = np.argsort(~neg[m], axis=1, kind="stable")              # the two negative corners first
        fo = np.take_along_axis(f[m], order, axis=1)
        a, b, c, d = (_REF[order][:, i] for i in range(4))
        fa, fb, fc, fd = (fo[:, i] for i in range(4))
        pac, pad, pbc, pbd = _cut(a, c, fa, fc), _cut(a, d, fa, fd), _cut(b, c, fb, fc), _cut(b, d, fb, fd)
        out[m] = 6.0 * (_tvol(a, pac, pad, b) + _tvol(pac, pad, b, pbc) + _tvol(pad, b, pbc, pbd))
    return out


def _below(S, f):
    """The parts of the tetrahedra with corners S (m, 4, D) on which the linear function with the values f (m, 4) at the
    corners is negative, as tetrahedra (k, 4, D), with the index (k,) of the one each came from (the construction of
    `cutcell.inside_subtets`, for corners in any coordinates; here barycentric in a cell of the mesh)."""
    neg = f < 0
    nneg = neg.sum(axis=1)
    m = np.flatnonzero(nneg == 4)
    out, owner = [S[m]], [m]
    for want in (1, 3):
        m = np.flatnonzero(nneg == want)
        if not len(m):
            continue
        r = np.arange(len(m))
        lone = np.argmax(neg[m] if want == 1 else ~neg[m], axis=1)
        o = _REST[lone]
        Sm, fm = S[m], f[m]
        p, fp = Sm[r, lone], fm[r, lone]
        q = [_cut(p, Sm[r, o[:, i]], fp, fm[r, o[:, i]]) for i in range(3)]
        if want == 1:
            tets = ((p, q[0], q[1], q[2]),)
        else:
            a, b, c = (Sm[r, o[:, i]] for i in range(3))
            tets = ((a, b, c, q[0]), (b, c, q[0], q[1]), (c, q[0], q[1], q[2]))
        for tet in tets:
            out.append(np.stack(tet, axis=1)), owner.append(m)
    m = np.flatnonzero(nneg == 2)
    if len(m):
        order = np.argsort(~neg[m], axis=1, kind="stable")
        So = np.take_along_axis(S[m], order[..., None], axis=1)
        fo = np.take_along_axis(f[m], order, axis=1)
        a, b, c, d = (So[:, i] for i in range(4))
        fa, fb, fc, fd = (fo[:, i] for i in range(4))
        pac, pad, pbc, pbd = _cut(a, c, fa, fc), _cut(a, d, fa, fd), _cut(b, c, fb, fc), _cut(b, d, fb, fd)
        for tet in ((a, pac, pad, b), (pac, pad, b, pbc), (pad, b, pbc, pbd)):
            out.append(np.stack(tet, axis=1)), owner.append(m)
    return np.concatenate(out), np.concatenate(owner)


def _contour(S, P, f):
    """The zero contour of the linear function with the values f (m, 4) at the corners S (m, 4, 4) (barycentric in
    cells whose corners are P (m, 4, 3) [m]): its area (m,) [m2] and its centroid's barycentric coordinates (m, 4)."""
    neg = f < 0
    nneg = neg.sum(axis=1)
    area, centre = np.zeros(len(f)), np.zeros((len(f), 4))

    def tri(a, b, c, m):
        pa, pb, pc = (np.einsum("ki,kid->kd", v, P[m]) for v in (a, b, c))
        return 0.5 * np.linalg.norm(np.cross(pb - pa, pc - pa), axis=1)
    for want in (1, 3):
        m = np.flatnonzero(nneg == want)
        if not len(m):
            continue
        r = np.arange(len(m))
        lone = np.argmax(neg[m] if want == 1 else ~neg[m], axis=1)
        o = _REST[lone]
        Sm, fm = S[m], f[m]
        q = [_cut(Sm[r, lone], Sm[r, o[:, i]], fm[r, lone], fm[r, o[:, i]]) for i in range(3)]
        area[m] = tri(q[0], q[1], q[2], m)
        centre[m] = (q[0] + q[1] + q[2]) / 3.0
    m = np.flatnonzero(nneg == 2)
    if len(m):
        order = np.argsort(~neg[m], axis=1, kind="stable")
        So = np.take_along_axis(S[m], order[..., None], axis=1)
        fo = np.take_along_axis(f[m], order, axis=1)
        a, b, c, d = (So[:, i] for i in range(4))
        fa, fb, fc, fd = (fo[:, i] for i in range(4))
        pac, pbc, pbd, pad = _cut(a, c, fa, fc), _cut(b, c, fb, fc), _cut(b, d, fb, fd), _cut(a, d, fa, fd)
        a1, a2 = tri(pac, pbc, pbd, m), tri(pac, pbd, pad, m)
        area[m] = a1 + a2
        whole = np.maximum(a1 + a2, 1e-300)[:, None]
        centre[m] = (a1[:, None] * (pac + pbc + pbd) + a2[:, None] * (pac + pbd + pad)) / (3.0 * whole)
    return area, centre


class _Groups:
    """phi's band divided among the markers at x (n, 3) with normals (n, 3), for one round of a move (module
    docstring): each node's marker (`node`, -1 for none) and each cell's parts, whole (`whole`, `whole_owner`) or as
    tetrahedra barycentric in the cell (`part`, `part_S`, `part_owner`, `part_volume`)."""

    def __init__(self, mesh, phi, x, normal, band=None):
        from .division import _marker_pieces                # (division imports this module)
        self.mesh, self.phi, self.x, self.normal, self.band, self.n = mesh, phi, x, normal, band, len(x)
        self.label, n_pieces = pieces(mesh, phi)
        top = patches(mesh, phi)
        self._piece_of = lambda pts: _marker_pieces(mesh, pts, self.label, top)
        piece = self._piece_of(x) if len(x) else np.zeros(0, dtype=np.int64)
        self.members = [np.flatnonzero(piece == p) for p in range(n_pieces)]
        self.trees = [cKDTree(x[mq]) if len(mq) else None for mq in self.members]
        f = phi[mesh.tets]
        self.fmin, self.fmax = f.min(axis=1), f.max(axis=1)
        self.node_piece = np.full(len(phi), -1, dtype=np.int64)        # a node of a metal cell: that cell's piece
        metal = np.flatnonzero(self.label >= 0)
        for corner in range(4):
            self.node_piece[mesh.tets[metal, corner]] = self.label[metal]
        self.G = np.full((len(phi), 3), np.nan)                         # the recovered gradient, where needed
        self.node = np.full(len(phi), -2, dtype=np.int64)               # -2: not yet asked
        self.covered = np.zeros(len(mesh.tets), dtype=bool)
        self._whole, self._whole_owner, self._part, self._part_S, self._part_owner = [], [], [], [], []
        self._freeze()

    def _nearest(self, p, points, k, point_normal):
        from .markers import _nearest                       # (markers imports this module)
        mq = self.members[p]
        d, o = _nearest(self.trees[p], points, k, self.normal[mq], point_normal)
        return d, np.where(o < len(mq), mq[np.minimum(o, len(mq) - 1)], -1)

    def cover(self, reach=0.0, cells=None):
        """Divide the cells of phi's band |phi| < reach [m] and every cell its zero contour cuts (or the given
        `cells`), those not divided yet."""
        mesh, phi = self.mesh, self.phi
        sel = (self.fmin < reach) & (self.fmax >= -reach) if cells is None else np.isin(np.arange(len(mesh.tets)), cells)
        cells = np.flatnonzero(sel & ~self.covered)
        if not len(cells) or not self.n:
            return
        self.covered[cells] = True
        T = mesh.tets[cells]
        nodes = np.unique(T)
        need = nodes[np.isnan(self.G[nodes, 0])]
        if len(need):
            near = np.unique(mesh.node_cells()[need].indices)
            self.G[need] = nodal_average(mesh, cell_gradients(mesh, phi, near), near)[need]
        foot = lambda v: mesh.points[v] - phi[v, None] * self.G[v]
        # each node's marker: the one owning its foot (module docstring); none beyond the band
        ask = nodes[self.node[nodes] == -2]
        self.node[ask] = -1
        if self.band is not None:
            ask = ask[np.abs(phi[ask]) < self.band]
        on = self.node_piece[ask]                            # each node's piece: its metal cells', or its foot's
        gas = on < 0
        if gas.any():
            on[gas] = self._piece_of(foot(ask[gas]))
        for p in range(len(self.members)):
            at = ask[on == p]
            if len(at) and len(self.members[p]):
                self.node[at] = self._nearest(p, foot(at), 1, self.G[at])[1][:, 0]
        # each cell's candidates: the markers its feet may be nearest (exactly, as `markers.ownership` argues: every
        # point of the feet lies within r of their centre, so its nearest marker lies within d0 + 2 r of the centre)
        F = foot(T.ravel()).reshape(-1, 4, 3)
        Gc = self.G[T].mean(axis=1)
        piece = self.label[cells].copy()
        gas = piece < 0
        if gas.any():
            piece[gas] = self._piece_of(F[gas].mean(axis=1))
        centre = F.mean(axis=1)
        r = np.linalg.norm(F - centre[:, None, :], axis=2).max(axis=1)
        cand = np.full((len(cells), MOVE_CANDIDATES), -1, dtype=np.int64)
        for p in range(len(self.members)):
            at = np.flatnonzero(piece == p)
            if not len(at) or not len(self.members[p]):
                continue
            k = min(MOVE_CANDIDATES, len(self.members[p]))
            d, o = self._nearest(p, centre[at], k, Gc[at])
            cand[at, :k] = np.where((o >= 0) & (d <= d[:, :1] + 2.0 * r[at, None]), o, -1)
        live = (cand >= 0).sum(axis=1)                       # (the candidates come nearest first)
        owner = cand[:, 0].copy()
        many = np.flatnonzero(live > 1)
        if len(many):
            c = cand[many, :int(live[many].max())]
            dist = np.linalg.norm(F[many][:, :, None, :] - self.x[np.maximum(c, 0)][:, None, :, :], axis=3)
            dist[np.broadcast_to((c < 0)[:, None, :], dist.shape)] = np.inf
            near = np.argmin(dist, axis=2)                   # each corner's foot's nearest candidate
            same = (near == near[:, :1]).all(axis=1)         # then one marker's (convex) cell holds the whole cell
            owner[many[same]] = c[same, near[same, 0]]
            many = many[~same]
        whole = (live > 0) & ~np.isin(np.arange(len(cells)), many)
        self._whole.append(cells[whole]), self._whole_owner.append(owner[whole])
        if len(many):
            # (cell, candidate u) pairs: u's part of the cell, where u's foot is nearer than each other candidate's
            L = live[many]
            pc = np.repeat(np.arange(len(many)), L)
            pu = np.concatenate([np.arange(v) for v in L])
            S, src = np.tile(np.eye(4), (len(pc), 1, 1)), np.arange(len(pc))
            for j in range(int(L.max())):
                act = (pu[src] != j) & (j < L[pc[src]])
                if not act.any():
                    continue
                q = src[act]
                xu, xw = self.x[cand[many[pc[q]], pu[q]]], self.x[cand[many[pc[q]], j]]
                h = 2.0 * np.einsum("kd,kvd->kv", xw - xu, F[many[pc[q]]] - 0.5 * (xu + xw)[:, None, :])
                Sa, oa = _below(S[act], np.einsum("kvi,ki->kv", S[act], h))
                S, src = np.concatenate([S[~act], Sa]), np.concatenate([src[~act], q[oa]])
            keep = np.abs(np.linalg.det(S)) > 0.0
            S, src = S[keep], src[keep]
            self._part.append(cells[many[pc[src]]]), self._part_S.append(S)
            self._part_owner.append(cand[many[pc[src]], pu[src]])
        self._freeze()

    def _freeze(self):
        join = lambda v, shape: np.concatenate(v) if v else np.zeros(shape, dtype=np.int64)
        self.whole, self.whole_owner = join(self._whole, 0), join(self._whole_owner, 0)
        self.part, self.part_owner = join(self._part, 0), join(self._part_owner, 0)
        self.part_S = np.concatenate(self._part_S) if self._part_S else np.zeros((0, 4, 4))
        self.part_volume = np.abs(np.linalg.det(self.part_S)) * self.mesh.volumes[self.part]
        self.held = self._held(self.phi)

    def _held(self, f):
        """The metal's volume [m3] in each whole cell and in each part, for the level set f."""
        mesh = self.mesh
        return (mesh.volumes[self.whole] * _fraction(f[mesh.tets[self.whole]]),
                self.part_volume * _fraction(np.einsum("kvi,ki->kv", self.part_S, f[mesh.tets[self.part]])))

    def swept(self, f):
        """Each group's swept volume [m3] (n,) from phi to the level set f, outward positive."""
        (w0, p0), (w1, p1) = self.held, self._held(f)
        return np.bincount(self.whole_owner, w1 - w0, minlength=self.n) + \
            np.bincount(self.part_owner, p1 - p0, minlength=self.n)

    def jacobian(self, f):
        """J (n, n), csr [m2]: J[i, j] the rate at which group i's swept volume grows with group j's outward offset, at
        the level set f (first order: each owned part's contour, weighted by its cell's nodes' hat functions over
        |grad f|)."""
        mesh, rows, cols, vals = self.mesh, [], [], []
        for cells, S, owner in ((self.whole, None, self.whole_owner), (self.part, self.part_S, self.part_owner)):
            T = mesh.tets[cells]
            S = np.broadcast_to(np.eye(4), (len(cells), 4, 4)) if S is None else S
            fv = np.einsum("kvi,ki->kv", S, f[T])
            neg = (fv < 0).sum(axis=1)
            cut = np.flatnonzero((neg > 0) & (neg < 4))
            if not len(cut):
                continue
            a, lam = _contour(np.asarray(S[cut]), mesh.points[T[cut]], fv[cut])
            g = np.linalg.norm(cell_gradients(mesh, f, cells[cut]), axis=1)
            w = a[:, None] * lam / np.maximum(g, 1e-300)[:, None]
            grp = self.node[T[cut]]
            ok = grp >= 0
            rows.append(np.broadcast_to(owner[cut, None], grp.shape)[ok]), cols.append(grp[ok]), vals.append(w[ok])
        if not rows:
            return coo_matrix((self.n, self.n)).tocsr()
        return coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(self.n, self.n)).tocsr()

    def offset(self, d):
        """phi with each grouped node offset by its group's distance d (n,) [m], outward positive, unreset."""
        out = self.phi.copy()
        at = np.flatnonzero(self.node >= 0)
        out[at] -= d[self.node[at]]
        return out


def swept_volumes(mesh, phi, new, x, normal):
    """Each group's swept volume [m3] (n,), outward positive, from phi to the level set `new`, by the division of the
    cells among the markers at x (n, 3) with normals (n, 3) that `move_by_volumes` uses (module docstring): the volume a
    move of phi to `new` would have to ask of each group."""
    x, normal = np.asarray(x, float).reshape(-1, 3), np.asarray(normal, float).reshape(-1, 3)
    groups = _Groups(mesh, np.asarray(phi, float), x, normal)
    groups.cover(cells=np.flatnonzero(inside_volumes(mesh, phi) != inside_volumes(mesh, new)))
    return groups.swept(np.asarray(new, float))


def _first_offsets(groups, target, reach):
    """Each group's first offset [m] (n,): the one at which its own swept volume, with every group offset alike, meets
    its target, between MOVE_PROFILE + 1 offsets from none to +-reach, linearly; none for a group that never responds."""
    out = np.zeros(groups.n)
    steps = np.arange(MOVE_PROFILE + 1) / MOVE_PROFILE
    for sign in (-1.0, 1.0):
        want = np.flatnonzero(np.sign(target) == sign)
        if not len(want):
            continue
        offsets = sign * reach * steps
        rise = np.array([np.zeros(len(want))] + [sign * groups.swept(groups.offset(np.full(groups.n, o)))[want]
                                                 for o in offsets[1:]])
        t = sign * target[want]
        k = np.clip((rise < t).sum(axis=0) - 1, 0, MOVE_PROFILE - 1)
        c = np.arange(len(want))
        lo, hi = rise[k, c], rise[k + 1, c]
        s = np.clip(np.where(hi > lo, (t - lo) / np.where(hi > lo, hi - lo, 1.0), 1.0), 0.0, 1.0)
        out[want] = np.where(rise[-1] > 0.0, offsets[k] + s * (offsets[k + 1] - offsets[k]), 0.0)
    return out


def _round(mesh, groups, target, reach, band, rounds, rtol, floor):
    """One round of a move: the groups' offsets, by Newton's method from `_first_offsets`, at which each group's swept
    volume meets its target [m3] (n,) within rtol of max(|target|, floor); then the reset keeping the contour. Returns
    (phi moved, each group's swept volume)."""
    n = groups.n
    moves = np.bincount(groups.node[groups.node >= 0], minlength=n) > 0    # groups with a node to offset
    scale = np.maximum(np.abs(target), floor)

    def state(d):
        raw = groups.offset(d)
        S = groups.swept(raw)
        r = np.where(moves, (target - S) / scale, 0.0)
        return raw, S, float(r @ r)
    d = _first_offsets(groups, target, reach)
    raw, S, f = state(d)
    for _ in range(rounds):
        if np.all(np.abs(target - S)[moves] <= rtol * scale[moves]):
            break
        J = groups.jacobian(raw)
        act = np.flatnonzero(moves & (J.diagonal() > 0.0))             # groups that respond to their own offset
        step = np.zeros(n)
        try:
            step[act] = splu(J[act][:, act].tocsc()).solve((target - S)[act])
        except RuntimeError:                                 # singular: each group on its own response
            step[act] = (target - S)[act] / J.diagonal()[act]
        for _ in range(MOVE_HALVINGS + 1):
            trial = np.clip(d + step, -reach, reach)
            got = state(trial)
            if got[2] < f:
                break
            step *= 0.5
        else:
            break                                            # no step lowers the misses: stop where they stand
        d, (raw, S, f) = trial, got
    new = redistance(mesh, raw, band=band, keep=True)
    return new, groups.swept(new)


def move_by_volumes(mesh, phi, x, normal, dV, band=None, rounds=MOVE_ROUNDS, rtol=MOVE_RTOL):
    """phi, a signed distance (negative in the metal), moved so that each marker's group (at x (n, 3), normals (n, 3))
    sweeps its volume dV (n,) [m3], outward positive (design section 21.6; planning ruling 2; module docstring). Within
    each round, at most `rounds` Newton iterations, until every group's miss is within rtol. A move whose first estimate,
    each group's dV over its area, is longer than MOVE_BAND_SHARE of the reset's `band` [m] is taken in that many equal
    rounds, the groups divided afresh on each round's phi (the markers left where they are); each round's phi is reset
    with the band, keeping its contour (`redistance(keep=True)`). Nodes beyond the band, and nodes whose foot no marker
    may own, keep their phi; so the band must exceed each round's move by more than the cells the contour crosses, whose
    nodes must hold distances (the body's, 13.6 mm, does). A group sweeps only the metal of the part of the mesh its marker
    owns: asked for more (a corner column narrows to nothing within about a millimetre on the stage 3 box), it misses.
    The per-piece volume correction is the caller's. Returns a `Moved`; with no change asked, phi as it was in no
    round."""
    phi = np.asarray(phi, dtype=float)
    x, normal = np.asarray(x, float).reshape(-1, 3), np.asarray(normal, float).reshape(-1, 3)
    dV = np.asarray(dV, dtype=float)
    n = len(x)
    scale = np.abs(dV).max(initial=0.0)
    if scale == 0.0:
        return Moved(phi.copy(), np.zeros(n), np.zeros(n), 0)
    floor = MOVE_FLOOR * scale
    swept, count, k = np.zeros(n), 1, 0
    while k < count:
        groups = _Groups(mesh, phi, x, normal, band)
        groups.cover(0.0)
        area = np.asarray(groups.jacobian(phi).sum(axis=1)).ravel()
        estimate = np.where(area > 0.0, dV / np.where(area > 0.0, area, 1.0), 0.0)
        if k == 0 and band is not None:
            count = max(1, int(np.ceil(np.abs(estimate).max() / (MOVE_BAND_SHARE * band) - 1e-12)))
        reach = MOVE_REACH * np.abs(estimate).max() / count
        reach = reach if band is None else min(reach, band)
        groups.cover(reach)
        phi, s = _round(mesh, groups, dV / count, reach, band, rounds, rtol, floor / count)
        swept += s
        k += 1
    return Moved(phi, swept, np.abs(swept - dV) / np.maximum(np.abs(dV), floor), count)
