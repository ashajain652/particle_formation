"""The derived surfaces of design §20.3: the outer surface phi from the skin base beta and the skins' thickness, and
the molten-rigid surface psi from the temperature, keeping only rigid metal connected to the core (crusts float)."""
import numpy as np

from .levelset import pieces, redistance

T_RIGID = 895.1   # K, half liquid on AA7075_scheil (design §8.3)


def outer_surface(mesh, beta, thickness, band=None):
    """phi = the zero contour of beta moved outward by the nodal skin thickness (m), reset to a distance."""
    return redistance(mesh, beta - thickness, band=band)


def molten_rigid_surface(mesh, phi, T, T_rigid=T_RIGID, core_point=None, band=None):
    """psi: signed distance to the rigid region, negative inside it.

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
