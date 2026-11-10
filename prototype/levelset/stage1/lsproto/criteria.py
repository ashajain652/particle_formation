"""Curvature of a distance field and the refinement targets of design §20.4-§20.5: the fine region (metal with
psi > -(d + lead)) at the fine size, and cells near either surface no larger than a quarter of its smallest principal
radius of curvature, down to the floor."""
import numpy as np

from .cutcell import cell_gradients
from .transfer import interpolate, locate, barycentric


def nodal_average(mesh, cell_values, cells=None):
    """Volume-weighted average of per-cell values onto the nodes (of all cells, or of the given cells only;
    nodes outside get NaN)."""
    cells = np.arange(len(mesh.tets)) if cells is None else np.asarray(cells)
    w = mesh.volumes[cells]
    flat = cell_values.reshape(len(w), -1)
    T = mesh.tets[cells]
    out = np.zeros((len(mesh.points), flat.shape[1]))
    wsum = np.zeros(len(mesh.points))
    for j in range(4):
        np.add.at(out, T[:, j], flat * w[:, None])
        np.add.at(wsum, T[:, j], w)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (out / wsum[:, None]).reshape((len(mesh.points),) + cell_values.shape[1:])


def hessian_field(mesh, phi, cells=None):
    """Nodal gradient and nodal (symmetrised) Hessian of a P1 field by twice-repeated gradient recovery, over all
    cells or over the given cells (values are exact recoveries only two rings inside that set)."""
    G = nodal_average(mesh, cell_gradients(mesh, phi, cells), cells)
    H = np.stack([cell_gradients(mesh, G[:, k], cells) for k in range(3)], axis=1)     # H[c, k, :] = d(G_k)/dx
    H = 0.5 * (H + np.transpose(H, (0, 2, 1)))
    return G, nodal_average(mesh, H, cells)


def principal_curvatures(mesh, phi, pts, fields=None, rings=3):
    """(k1, k2), k1 >= k2, of the level surfaces of the distance field phi through pts; positive where convex
    (a sphere of radius R has k1 = k2 = 1/R). Without `fields`, the recovery runs only on the cells within
    `rings` rings of the cells holding the points."""
    cell, b = locate(mesh, pts)
    if fields is None:
        fields = hessian_field(mesh, phi, mesh.ring(cell, rings))
    G, H = fields
    g = np.einsum("ni,nik->nk", b, G[mesh.tets[cell]])
    h = np.einsum("ni,nikl->nkl", b, H[mesh.tets[cell]])
    gn = np.linalg.norm(g, axis=1)
    n = g / np.maximum(gn, 1e-300)[:, None]
    P = np.eye(3)[None] - n[:, :, None] * n[:, None, :]
    S = P @ h @ P / np.maximum(gn, 1e-300)[:, None, None]
    w, v = np.linalg.eigh(S)
    along = np.abs(np.einsum("nk,nkj->nj", n, v))           # |n . eigenvector|
    drop = np.argmax(along, axis=1)                          # the normal direction's eigenvalue
    keep = np.array([[j for j in range(3) if j != d] for d in range(3)])[drop]
    kk = np.take_along_axis(w, keep, axis=1)
    return kk.max(axis=1), kk.min(axis=1)


def kappa_max_nodal(mesh, phi, band=None, fields=None):
    """Largest |principal curvature| of the level surface through every node; with `band` (m) only nodes with
    |phi| < band are evaluated (others 0), the recovery running on the cells around them."""
    out = np.zeros(len(mesh.points))
    nodes = np.arange(len(mesh.points)) if band is None else np.flatnonzero(np.abs(phi) < band)
    if len(nodes) == 0:
        return out
    if fields is None and band is not None:
        seed = np.unique(mesh.node_cells()[nodes].indices)
        fields = hessian_field(mesh, phi, mesh.ring(seed, 2))
    k1, k2 = principal_curvatures(mesh, phi, mesh.points[nodes], fields)
    out[nodes] = np.nan_to_num(np.maximum(np.abs(k1), np.abs(k2)))
    return out


def make_target(mesh, phi, psi, kmax_phi, kmax_psi, fine=0.7e-3, margin=4.0e-3, lead=0.6e-3, floor=0.25e-3,
                curvature_factor=0.25):
    """target_size(x, size) for boxmesh.build, evaluated from fields on the current mesh."""
    def target(x, size):
        cell, b = locate(mesh, x)
        T = mesh.tets[cell]
        ph = np.einsum("ni,ni->n", b, phi[T])
        ps = np.einsum("ni,ni->n", b, psi[T])
        reach = 0.94 * np.asarray(size)                   # centroid to farthest vertex of a Kuhn tetrahedron: 0.935 h
        out = np.full(len(x), np.inf)
        in_fine = (ph < reach) & (ps > -(margin + lead + reach))
        out[in_fine] = fine
        for dist, km in ((ph, kmax_phi), (ps, kmax_psi)):
            near = np.abs(dist) < 2 * reach
            k = np.einsum("ni,ni->n", b, km[T])
            want = np.where(k > 0, curvature_factor / np.maximum(k, 1e-300), np.inf)
            out = np.where(near, np.minimum(out, np.maximum(want, floor)), out)
        return out
    return target
