"""The derived surfaces: the skin base beta from the gas surface phi and the skin depth (design §21.4), and the
molten-rigid surface psi from the temperature, keeping only rigid metal connected to the core (crusts float). (The outer
surface phi built from beta and the skins' thickness, design §20.3, was retired by plan 03b Task 5: phi is evolved.)"""
import numpy as np

from .levelset import pieces, redistance


def skin_base(mesh, phi, depth, band=None):
    """beta, the skin base (design §21.4): the surface `depth` [m] inside the gas surface phi, the zero contour of
    phi + depth exactly, reset to a distance with `band` [m] keeping that contour (`redistance(..., keep=True)`: the
    nodes of the cells it cuts keep phi + depth). phi must be a signed distance (negative in the metal) reaching beyond
    `depth` into the metal, so a band no deeper than `depth` raises ValueError (it should reach beyond by a cell).
    Moving a surface inward by a distance keeps its convex edges sharp and rounds its concave corners to radius `depth`.

    Why the contour is kept (plan 03b Task 5, against the controller's ruling T2-b's plain reset). The layer the division
    measures is the metal of phi less that of phi + depth (`division.divide`), and the volume's metal is beta's; the
    books close only if the two agree, cell by cell. The plain reset moves the contour wherever it bends within a cell,
    inward at beta's convex edges: on the stage 3 box it held 0.67 mm3 less than phi + depth (Task 2's report), which
    would leave the skins' columns and the volume overlapping or apart by that much at every acceptance."""
    if band is not None and band <= depth:
        raise ValueError("the skin base needs phi a distance beyond the depth {:.3e} m into the metal, but the band is "
                         "{:.3e} m".format(depth, band))
    return redistance(mesh, np.asarray(phi, dtype=float) + depth, band=band, keep=True)


def molten_rigid_surface(mesh, phi, T, T_rigid, core_point=None, band=None):
    """psi: signed distance to the rigid region, negative inside it; `T_rigid` [K] the material's rigid temperature
    (the shared `Material.T_rigid`, half liquid: 895.1 K on AA7075_scheil, design §8.3).

    The rigid region is the metal (phi < 0) below T_rigid that is connected through rigid metal to the core: the
    piece containing core_point (default: the deepest node of phi), or the largest rigid piece if the core point is
    not rigid. Other rigid pieces are crusts. psi is clipped so that rigid metal never lies outside the metal,
    psi >= phi. With `band` (m) both distance resets are limited to that band. Returns (psi, crust_cells), crust_cells a boolean per cell."""
    g = redistance(mesh, T - T_rigid, band=band)
    if not np.isfinite(g).any() or (g > 0).all():      # nothing below T_rigid anywhere
        g = np.where(T < T_rigid, -np.inf, np.inf)
    q = np.maximum(g, phi)
    label, n = pieces(mesh, q)
    crust = np.zeros(len(mesh.tets), bool)
    if n == 0:
        return np.full(len(mesh.points), np.inf), crust
    if core_point is None:
        core_point = mesh.points[np.argmin(phi)]
    from .transfer import locate
    cell, _ = locate(mesh, np.atleast_2d(core_point))
    keep = label[cell[0]] if label[cell[0]] >= 0 else 0
    crust = (label >= 0) & (label != keep)
    if crust.any():
        core_nodes = np.zeros(len(mesh.points), bool)
        core_nodes[mesh.tets[label == keep].ravel()] = True
        crust_nodes = np.zeros(len(mesh.points), bool)
        crust_nodes[mesh.tets[crust].ravel()] = True
        lift = crust_nodes & ~core_nodes & (q < 0)
        q = np.where(lift, -q, q)
    psi = redistance(mesh, q, band=band)
    return np.maximum(psi, phi), crust
