"""The moving 3D mesh (layered-skin design, docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md §5.3).

The 3D mesh's surface nodes are the corners whose columns carry the skins; they move by the metal crossing the rigid
front, and the 3D nodes within a depth D below the surface (the deforming zone) move with them, each by a weighted
average of the movements of the surface nodes near it, the weight falling linearly from one at the surface to zero at
depth D. After every move the temperature field is remapped onto the moved nodes and the energy books are closed
exactly; before any element degrades the mesh is rebuilt inside the same surface triangulation.

The heat equation's moving-mesh term is applied as a remap at each move (the operator-split, "Lagrangian plus remap"
form of the arbitrary Lagrangian-Eulerian method): the nodal enthalpy is carried to each moved node along its
displacement by the nodal gradient, h_i' = h_i + grad(h)_i . u_i, which is exact for uniform and linear fields; then
the total enthalpy is set to its old value less what left through the moving surface, by one uniform enthalpy offset on
the moved nodes, so the discrete geometric conservation law holds (a uniform field stays exactly uniform) and the books
close to round-off. The conduction step that follows runs on the new geometry with the backends' fixed-mesh scheme.

The corners also move along the surface: `tangential_smoothing` nudges each toward the mean of its neighbours in its
tangent plane (design section 5.1), which keeps the fixed surface triangulation even. Measured while planning
(sub-plan 19, 2026-10-11) on the 100 mm sphere receding at its nose by 0.5 mm per step in proportion to how directly
each corner faces the flow: without it the run failed at 20.5 mm of nose recession, with it at 27.0 mm. What stopped it
then was the shape, not the triangles: the angle between neighbouring surface triangles (`max_dihedral`) grew from
4 to 82 degrees as that motion sharpened the shoulder, and elements whose four nodes all lie on a fold that sharp are
flattened by their own corners' motion, which no rebuild inside the same surface can undo. The fold is reported.

What left through the surface is booked per corner: the volume each corner sweeps is its dual area times its normal
displacement (positive inward), rescaled so the corners' volumes add up exactly to the change of the mesh's volume, and
it carries the corner's own nodal enthalpy. The layered skins take exactly that mass and enthalpy (stage 2).

Rebuild: gmsh meshes the volume inside the current surface triangulation without touching it (a discrete surface in
a geo surface loop, generate(3) only; checked 2026-10-11: every surface triangle and corner is kept to the last bit,
though gmsh renumbers nodes, so the corners are matched back by their coordinates). The enthalpy is interpolated onto
the new nodes in the old mesh and its total made exact by a uniform offset, as the 2026-09-27 remeshing design does;
the pre-correction residual is reported with every rebuild. The validity gate stops on a non-manifold surface, more
than one connected body, or a non-positive enclosed volume."""
import time
from dataclasses import dataclass

import numpy as np

from . import mesh as mesh_mod

DEFORM_DEPTH = 15.0e-3       # m, depth D of the deforming zone: the dense band's depth (mesh.DEFAULT_BAND)
RATIO_MIN = 0.05             # radius ratio below which a move is refused and the mesh rebuilt (2026-09-27 design §7.1)
MIN_GAP = 20                 # macro steps between rebuilds triggered by compression (a refused move ignores it)
ZONE_NEIGHBOURS = 12         # surface nodes averaged for each interior node of the zone
SMOOTH_LAMBDA = 0.5          # tangential smoothing: share of the way to the neighbours' mean per move
SMOOTH_CLAMP = 0.35          # of the mean surface edge: the most a corner is nudged per move (2026-09-27 design)
_REF_GRADS = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])


class RemeshError(RuntimeError):
    """The surface cannot be rebuilt: non-manifold, more than one body, or a non-positive enclosed volume."""


# -- element and surface geometry ----------------------------------------------------------------------------------
def signed_volumes(points, tets):
    x = points[tets]
    return np.einsum("ij,ij->i", x[:, 1] - x[:, 0], np.cross(x[:, 2] - x[:, 0], x[:, 3] - x[:, 0])) / 6.0


def radius_ratio(points, tets):
    """3 r_in / R_circ per tetrahedron: 1 for a regular one, 0 for a flat one."""
    x = points[tets]
    u, v, w = x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]
    triple = np.einsum("ij,ij->i", u, np.cross(v, w))
    vol = np.abs(triple) / 6.0
    faces = (np.cross(v - u, w - u), np.cross(w, v), np.cross(u, w), np.cross(v, u))
    area = 0.5 * sum(np.linalg.norm(f, axis=1) for f in faces)
    r_in = 3.0 * vol / np.where(area > 0.0, area, 1.0)
    c = ((u * u).sum(1)[:, None] * np.cross(v, w) + (v * v).sum(1)[:, None] * np.cross(w, u)
         + (w * w).sum(1)[:, None] * np.cross(u, v))
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.linalg.norm(c, axis=1) / np.abs(2.0 * triple)
        ratio = np.where(vol > 0.0, 3.0 * r_in / R, 0.0)
    return np.nan_to_num(ratio, nan=0.0)


def triangle_quality(points, faces):
    """2 r_in / R_circ per triangle: 1 for an equilateral one, 0 for a flat one."""
    a, b, c = (points[faces[:, i]] for i in range(3))
    la, lb, lc = np.linalg.norm(b - c, axis=1), np.linalg.norm(c - a, axis=1), np.linalg.norm(a - b, axis=1)
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    s = 0.5 * (la + lb + lc)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = 2.0 * (area / s) / (la * lb * lc / (4.0 * area))
    return np.nan_to_num(q, nan=0.0)


def max_dihedral(points, faces):
    """The largest angle [degrees] between the normals of two surface triangles sharing an edge: how sharply the
    surface folds at the scale of one triangle."""
    a, b, c = (points[faces[:, i]] for i in range(3))
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1)[:, None]
    edges = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(faces)), 3)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    e, o = edges[order], owner[order]
    pair = np.all(e[1:] == e[:-1], axis=1)
    cos = np.einsum("ij,ij->i", n[o[1:][pair]], n[o[:-1][pair]])
    return float(np.degrees(np.arccos(np.clip(cos.min(), -1.0, 1.0)))) if cos.size else 0.0


def element_gradients(points, tets):
    """Gradients of the four P1 basis functions in every tetrahedron, (ne, 4, 3), and the volumes (ne,)."""
    x = points[tets]
    J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
    vol = np.abs(np.linalg.det(J)) / 6.0
    return np.einsum("eij,aj->eai", np.linalg.inv(J).transpose(0, 2, 1), _REF_GRADS), vol


def nodal_gradient(points, tets, f, elements=None):
    """Volume-weighted average over each node's elements of grad f, f a nodal P1 field: exact for a linear f. With
    `elements`, only those elements are used: the gradient is then right at every node all of whose elements are in
    the list (the moving nodes, when the list is the elements they touch)."""
    if elements is not None:
        tets = tets[elements]
    grads, vol = element_gradients(points, tets)
    g = np.einsum("ea,eai->ei", f[tets], grads) * vol[:, None]
    out = np.zeros((len(points), 3))
    for a in range(4):
        np.add.at(out, tets[:, a], g)
    w = np.bincount(tets.ravel(), weights=np.repeat(vol, 4), minlength=len(points))
    return out / np.where(w > 0.0, w, 1.0)[:, None]


def lumped_masses(points, tets, rho):
    """rho times a quarter of each incident element's volume, per node: the backends' lumped capacity weights."""
    vol = np.abs(signed_volumes(points, tets))
    return rho * np.bincount(tets.ravel(), weights=np.repeat(vol / 4.0, 4), minlength=len(points))


def orient_outward(points, faces, normals):
    """The triangles reordered where needed so that (b - a) x (c - a) points along their outward normals."""
    faces = np.array(faces, dtype=np.int64)
    a, b, c = (points[faces[:, i]] for i in range(3))
    flip = np.einsum("ij,ij->i", np.cross(b - a, c - a), normals) < 0.0
    faces[flip] = faces[flip][:, [0, 2, 1]]
    return faces


def corner_frame(points, faces, corners):
    """Area-weighted outward vertex normals and dual areas (a third of each touching triangle) at `corners`."""
    a, b, c = (points[faces[:, i]] for i in range(3))
    n2 = np.cross(b - a, c - a)                         # twice the area times the unit normal (outward for mesh faces)
    normal = np.zeros((len(points), 3))
    area = np.zeros(len(points))
    for i in range(3):
        np.add.at(normal, faces[:, i], n2)
        np.add.at(area, faces[:, i], np.linalg.norm(n2, axis=1) / 6.0)
    normal = normal[corners]
    return normal / np.linalg.norm(normal, axis=1)[:, None], area[corners]


def surface_edges(faces):
    """The unique edges of a triangulation, (ne, 2), each with its smaller node id first."""
    return np.unique(np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1), axis=0)


def tangential_smoothing(points, faces, normal, lam=SMOOTH_LAMBDA, clamp=SMOOTH_CLAMP, edges=None):
    """Displacements that even out the surface triangles without moving the surface (layered-skin design section 5.1):
    each corner moves `lam` of the way toward the mean of its neighbours, projected onto its tangent plane, the move
    clamped to `clamp` mean edge lengths. Tangential moves sweep no volume to first order. `edges` (from
    `surface_edges(faces)`) saves recomputing them every move."""
    if edges is None:
        edges = surface_edges(faces)
    n = len(points)
    total = np.zeros((n, 3))
    np.add.at(total, edges[:, 0], points[edges[:, 1]])
    np.add.at(total, edges[:, 1], points[edges[:, 0]])
    degree = np.bincount(edges.ravel(), minlength=n)
    L = total / np.maximum(degree, 1)[:, None] - points
    L -= np.einsum("ij,ij->i", L, normal)[:, None] * normal
    step = lam * L
    limit = clamp * np.linalg.norm(points[edges[:, 0]] - points[edges[:, 1]], axis=1).mean()
    size = np.linalg.norm(step, axis=1)
    return step * np.minimum(1.0, limit / np.where(size > 0.0, size, 1.0))[:, None]


def swept_volumes(old, new, faces, faces_corner, n_corners):
    """The volume each corner's surface swept, positive inward: every surface triangle's prism between its old and new
    positions, split into three tetrahedra, a third of it to each of the triangle's corners. A corner sliding along a
    flat face sweeps exactly nothing; a corner on a sharp edge books only what the faces it touches swept (Task 4
    review, 2026-10-12, replacing dual area times normal displacement, which mis-shares at edges and, scaled to the
    total, blew up wherever the first-order volumes nearly cancelled). `faces` are the outward triangles' node ids,
    `faces_corner` the same triangles in corner numbering."""
    a, b, c = (old[faces[:, i]] for i in range(3))
    A, B, C = (new[faces[:, i]] for i in range(3))

    def vol(p, q, r, s):
        return np.einsum("ij,ij->i", q - p, np.cross(r - p, s - p)) / 6.0
    prism = -(vol(a, b, c, A) + vol(b, c, A, B) + vol(c, A, B, C))
    out = np.zeros(n_corners)
    for i in range(3):
        np.add.at(out, faces_corner[:, i], prism / 3.0)
    return out


# -- the deforming zone ----------------------------------------------------------------------------------------------
class ZoneMap:
    """Which 3D nodes move with the surface, and how: surface nodes by their own displacement; every other node within
    `depth` of the surface by w_i times an inverse-square-distance average of its ZONE_NEIGHBOURS nearest surface
    nodes' displacements, w_i = 1 - d_i / depth (d_i its distance to the nearest surface node). Built once per mesh."""

    def __init__(self, points, corners, depth=DEFORM_DEPTH, neighbours=ZONE_NEIGHBOURS):
        from scipy.spatial import cKDTree
        if depth <= 0.0:
            raise ValueError("the deforming zone needs a positive depth")
        self.n_points, self.corners, self.depth = len(points), np.asarray(corners, dtype=np.int64), float(depth)
        tree = cKDTree(points[self.corners])
        d = tree.query(points)[0]
        w = np.clip(1.0 - d / self.depth, 0.0, 1.0)
        is_corner = np.zeros(len(points), dtype=bool)
        is_corner[self.corners] = True
        self.interior = np.flatnonzero((w > 0.0) & ~is_corner)
        k = min(neighbours, len(self.corners))
        dist, idx = tree.query(points[self.interior], k=k)
        dist, idx = dist.reshape(len(self.interior), k), idx.reshape(len(self.interior), k)
        K = 1.0 / (dist + 1e-3 * self.depth) ** 2
        self.idx, self.K = idx, K / K.sum(axis=1)[:, None]
        self.w = w[self.interior]

    def displacement(self, u_corners):
        """(n_corners, 3) corner displacements -> (n_points, 3) node displacements."""
        u = np.asarray(u_corners, dtype=float)
        U = np.zeros((self.n_points, 3))
        U[self.corners] = u
        U[self.interior] = self.w[:, None] * np.einsum("ik,ikj->ij", self.K, u[self.idx])
        return U


# -- moving the mesh -------------------------------------------------------------------------------------------------
@dataclass
class MoveResult:
    refused: bool                       # the move would invert or degrade an element: nothing was changed
    swept_volume: np.ndarray            # (n_corners,) m3 each corner's motion took out of the 3D model (negative: added)
    swept_enthalpy: np.ndarray          # (n_corners,) J that volume carried
    correction: float                   # J, the uniform offset that closed the books (pre-correction residual)
    worst_ratio: float                  # smallest radius ratio among the elements the move touched
    rebuilt: bool = False               # move_with_policy rebuilt the mesh before (or to make) this move


def default_sizes(mesh):
    """gmsh sizes for a rebuild, from the generator's parameters: a sphere's dense band, or a box's uniform h."""
    p = mesh.params
    if "h_surface_m" in p:
        return {"h_surface": p["h_surface_m"], "h_core": p["h_core_m"], "band": p.get("band_m", mesh_mod.DEFAULT_BAND),
                "ramp": max(p.get("band_ramp_m", mesh_mod.DEFAULT_BAND_RAMP), 1e-9)}
    h = p.get("h", mesh_mod.DEFAULT_H_SURFACE)
    return {"h_surface": h, "h_core": h, "band": 0.0, "ramp": 1.0}


class MovingMesh:
    """A VolumeMesh, its thermal solver and material, moved by its corners (the surface nodes, in a fixed order that
    survives rebuilds) and rebuilt when needed. `mesh` and the solver's mesh are replaced at every rebuild; callers
    read `self.mesh` afresh. `corners` holds the current node id of each corner."""

    def __init__(self, mesh, solver, material, depth=DEFORM_DEPTH, ratio_min=RATIO_MIN, min_gap=MIN_GAP,
                 sizes=None):
        self.mesh, self.solver, self.material = mesh, solver, material
        self.depth, self.ratio_min, self.min_gap = float(depth), float(ratio_min), int(min_gap)
        self.sizes = sizes if sizes is not None else default_sizes(mesh)
        s = mesh.surface()
        self.corners = np.unique(s.faces)
        self._set_surface(orient_outward(mesh.points, s.faces, s.normals))
        self.compression = np.zeros(len(self.corners))      # inward normal motion of each corner since the last rebuild
        self.steps_since_rebuild, self.rebuilds, self.last_rebuild = 0, 0, None

    def _set_surface(self, faces_nodes):
        local = np.full(self.mesh.n_nodes, -1, dtype=np.int64)
        local[self.corners] = np.arange(len(self.corners))
        self.faces_corner = local[faces_nodes]              # surface triangles in corner numbering (rebuild-proof)
        self.edges_corner = surface_edges(self.faces_corner)
        self.zone = ZoneMap(self.mesh.points, self.corners, self.depth)

    def frame(self):
        """Outward normals and dual areas at the corners, on the current geometry."""
        return corner_frame(self.mesh.points, self.corners[self.faces_corner], self.corners)

    def smoothing(self, normal=None, lam=None, clamp=None):
        """Tangential displacements of the corners that even out the surface triangles (add them to a move)."""
        if normal is None:
            normal = self.frame()[0]
        return tangential_smoothing(self.mesh.points[self.corners], self.faces_corner, normal,
                                    SMOOTH_LAMBDA if lam is None else lam, SMOOTH_CLAMP if clamp is None else clamp,
                                    self.edges_corner)

    def surface_quality(self):
        """The worst surface triangle's quality, 2 r_in / R_circ."""
        return float(triangle_quality(self.mesh.points[self.corners], self.faces_corner).min())

    def surface_fold(self):
        """The largest angle between neighbouring surface triangles' normals [degrees]."""
        return max_dihedral(self.mesh.points[self.corners], self.faces_corner)

    def move(self, u_corners):
        """Move the corners by `u_corners` (n_corners, 3), the zone with them; remap the enthalpy; close the books.
        Refuses, changing nothing, if any element the move touches would invert or fall below the radius-ratio floor."""
        mesh, solver, mat = self.mesh, self.solver, self.material
        u = np.asarray(u_corners, dtype=float)
        U = self.zone.displacement(u)
        is_moving = np.any(U != 0.0, axis=1)
        moving = np.flatnonzero(is_moving)
        empty = np.zeros(len(self.corners))
        if moving.size == 0:
            return MoveResult(False, empty, empty.copy(), 0.0, 1.0)
        touched = np.flatnonzero(is_moving[mesh.tets].any(axis=1))
        new_points = mesh.points + U
        vol_new = signed_volumes(new_points, mesh.tets[touched])
        sign = np.sign(signed_volumes(mesh.points, mesh.tets[touched]))
        ratio = radius_ratio(new_points, mesh.tets[touched])
        worst = float(ratio.min())
        if np.any(vol_new * sign <= 0.0) or worst < self.ratio_min:
            return MoveResult(True, empty, empty.copy(), 0.0, worst)
        h = mat.enthalpy(solver.temperature())
        h_new = h + np.einsum("ij,ij->i", nodal_gradient(mesh.points, mesh.tets, h, touched), U)
        E_old = solver.energy()
        normal, area = self.frame()
        faces = self.corners[self.faces_corner]
        dV = swept_volumes(mesh.points, new_points, faces, self.faces_corner, len(self.corners))
        V_old = float(np.abs(signed_volumes(mesh.points, mesh.tets[touched])).sum())     # only touched elements change
        mesh.points[:] = new_points
        mesh.derived_points[:] = new_points
        # the prisms tile what the surface swept up to how their shared side faces are split; that remainder (third
        # order in the displacement) goes to the corners in proportion to what they swept, so the books are exact
        remainder = (V_old - float(np.abs(vol_new).sum())) - dV.sum()
        if np.abs(dV).sum() > 0.0:
            dV = dV + remainder * np.abs(dV) / np.abs(dV).sum()
        swept_h = mat.rho * dV * h[self.corners]
        solver.update_geometry(touched)
        solver.set_temperature(mat.temperature_from_enthalpy(h_new))
        correction = (E_old - swept_h.sum()) - solver.energy()
        m = lumped_masses(mesh.points, mesh.tets[touched], mat.rho)      # complete at the moving nodes
        h_new[moving] += correction / m[moving].sum()
        solver.set_temperature(mat.temperature_from_enthalpy(h_new))
        self.compression += np.maximum(-np.einsum("ij,ij->i", normal, u), 0.0)
        self.steps_since_rebuild += 1
        return MoveResult(False, dV, swept_h, float(correction), worst)
