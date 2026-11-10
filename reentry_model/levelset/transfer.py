"""Point location, P1 interpolation and the conservative transfer of P1 fields between two adaptive meshes built
from the same base (design §20.5: exact on refinement, conserving the integral where cells are coarsened)."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg
from scipy.spatial import cKDTree


def _geom(mesh):
    c = mesh._cache
    if "kd" not in c:
        p = mesh.points[mesh.tets]
        J = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=2)   # columns = edges
        c["Jinv"] = np.linalg.inv(J)
        c["x0"] = p[:, 0]
        c["kd"] = cKDTree(p.mean(axis=1))
    return c["kd"], c["Jinv"], c["x0"]


def barycentric(mesh, cells, pts):
    _, Jinv, x0 = _geom(mesh)
    l = np.einsum("nij,nj->ni", Jinv[cells], pts - x0[cells])
    return np.column_stack([1.0 - l.sum(axis=1), l])


def locate(mesh, pts, k=4):
    """(cell, barycentric) for each point; points outside the box get the cell with the least negative coordinate."""
    pts = np.atleast_2d(np.asarray(pts, float))
    kd, _, _ = _geom(mesh)
    cell = np.full(len(pts), -1)
    best = np.full(len(pts), -np.inf)
    todo = np.arange(len(pts))
    for kk in (k, 4 * k, 16 * k, 64 * k):
        if len(todo) == 0:
            break
        kk = min(kk, len(mesh.tets))
        _, cand = kd.query(pts[todo], k=kk)
        cand = cand.reshape(len(todo), kk)
        for j in range(kk):
            b = barycentric(mesh, cand[:, j], pts[todo])
            mn = b.min(axis=1)
            better = mn > best[todo]
            cell[todo[better]], best[todo[better]] = cand[better, j], mn[better]
        todo = todo[best[todo] < -1e-10]
    return cell, barycentric(mesh, cell, pts)


def interpolate(mesh, u, pts):
    cell, b = locate(mesh, pts)
    return np.einsum("ni,ni->n", b, u[mesh.tets[cell]])


# 4-point rule, exact for quadratics on a tetrahedron (barycentric coordinates, equal weights 1/4)
_A, _B = 0.5854101966249685, 0.1381966011250105
_QP = np.array([[_A, _B, _B, _B], [_B, _A, _B, _B], [_B, _B, _A, _B], [_B, _B, _B, _A]])


def mass_matrix(mesh):
    n = len(mesh.points)
    loc = (np.ones((4, 4)) + np.eye(4)) / 20.0
    vals = mesh.volumes[:, None, None] * loc[None]
    rows = np.repeat(mesh.tets, 4, axis=1).ravel()
    cols = np.tile(mesh.tets, (1, 4)).ravel()
    return coo_matrix((vals.ravel(), (rows, cols)), shape=(n, n)).tocsr()


def project(old, new, u_old, rtol=1e-13):
    """L2 projection of the P1 field u_old on `old` onto the P1 space of `new` (both built from one base).

    Each region is integrated on the finer of the two meshes there, so the right-hand side is exact: a new leaf
    at least as deep as the old leaf holding its centroid is integrated on itself; old leaves deeper than the new
    leaf holding their centroid are integrated on themselves. The integral of u is conserved to the solver's
    tolerance, and a field that lies in the new space (pure refinement) is reproduced exactly."""
    n_new = len(new.points)
    b = np.zeros(n_new)
    # case A: new leaves no coarser than the old leaf they lie in
    oc, _ = locate(old, new.centroids)
    A = new.level >= old.level[oc]
    for q in _QP:
        xq = np.einsum("j,njk->nk", q, new.points[new.tets[A]])
        cq, bq = locate(old, xq)
        uq = np.einsum("ni,ni->n", bq, u_old[old.tets[cq]])
        np.add.at(b, new.tets[A], (new.volumes[A] / 4 * uq)[:, None] * q[None, :])
    # case B: old leaves finer than the new leaf they lie in
    nc, _ = locate(new, old.centroids)
    B = old.level > new.level[nc]
    for q in _QP:
        xq = np.einsum("j,njk->nk", q, old.points[old.tets[B]])
        uq = u_old[old.tets[B]] @ q
        cq, bq = locate(new, xq)
        np.add.at(b, new.tets[cq], (old.volumes[B] / 4 * uq)[:, None] * bq)
    M = mass_matrix(new)
    x0 = interpolate(old, u_old, new.points)
    x, info = cg(M, b, x0=x0, rtol=rtol, atol=0.0, maxiter=2000)
    if info != 0:
        raise RuntimeError("projection did not converge (cg info {})".format(info))
    return x


def integral(mesh, u):
    return float((mesh.volumes * u[mesh.tets].mean(axis=1)).sum())
