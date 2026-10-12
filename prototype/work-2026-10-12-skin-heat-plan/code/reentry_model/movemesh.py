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
(sub-plan 19, Task 8, 2026-10-12, one thread) on the 100 mm sphere receding at its nose by 0.5 mm per step in
proportion to how directly each corner faces the flow: without it the run failed at 25.0 mm of nose recession, with it
at 29.5 mm. What stopped it was the shape, not the triangles: the angle between neighbouring surface triangles
(`max_dihedral`) grew from 19.5 degrees at 10 mm of recession and 51.9 at 20 mm to about 95 degrees as that motion
sharpened the shoulder, and elements whose four nodes all lie on a fold that sharp are flattened by their own corners'
motion, which no rebuild inside the same surface can undo. The fold is reported.

What left through the surface is booked per corner by the exact volume of the prism each surface triangle sweeps
between its old and new positions (`swept_volumes`). Each corner is booked what its own motion swept, to first order
A_t (u_k . n_t) / 3 for the triangle t at corner k (A_t the triangle's area, n_t its outward unit normal, u_k the
corner's displacement), and the second-order rest of the prism, the part of its exact volume that is not the sum of its
three corners' own sweeps, is shared among those corners in proportion to the absolute values of their own sweeps. So a
corner sliding along a flat face books nothing, and a corner that does not move is never credited with what its moving
neighbours swept (a third of each prism to every corner, as sub-plan 19 booked it, credited them with exactly that).
The third-order remainder by which the prisms miss the change of the mesh's volume, which arises in how the prisms'
shared side faces are split, is then spread over the corners in proportion to their own |volume| (in `move`), so the
corners' volumes add up exactly to that change. Each corner's volume carries the corner's own nodal enthalpy before the
move (the surface value, which the exact moving-front solution favours over the mean of the old and remapped values).
The layered skins take exactly that mass and enthalpy (stage 2).

Rebuild: gmsh meshes the volume inside the current surface triangulation without touching it (a discrete surface in
a geo surface loop, generate(3) only; checked 2026-10-11: every surface triangle and corner is kept to the last bit,
though gmsh renumbers nodes, so the corners are matched back by their coordinates, and the rebuilt boundary is checked
against the corner triangles). The enthalpy is carried onto the new mesh by an L2 projection of the old P1 (piecewise
linear) field (`project_enthalpy`): the new nodal values h solve M h = b, M being the new mesh's consistent mass matrix
and b_i the integral of the shape function i times the old field over the new elements by a four-point quadrature rule,
by preconditioned conjugate gradients. The projection reproduces a linear field exactly and keeps the energy where it
was, region by region, to the rule's accuracy; what the quadrature still misses is made exact by a uniform offset over
the new deforming zone only, where the field varies, as the 2026-09-27 remeshing design does for the whole body. The
pre-correction residual is reported with every rebuild and summed in `rebuild_residual_total`. The validity gate stops
on a non-manifold surface, more than one connected body, or a non-positive enclosed volume."""
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
NEUTRAL_PASSES = 3           # normal corrections that make the smoothing volume-neutral per corner (sub-plan 20 Task 0):
                             # with own sweeps each pass removes about 95 % of what is left, the worst corner's sweep
                             # on the 100 mm sphere going from 0.46 um of equivalent depth to 3e-11 m after three
                             # (measured 2026-10-12)
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
    positions (exact, split into three tetrahedra), shared among the triangle's corners by what each corner's own motion
    swept -- to first order A_t (u_k . n_t) / 3 for corner k, n_t the triangle's outward normal -- with the second-order
    rest of the prism shared in proportion to those. A corner sliding along a flat face sweeps exactly nothing, and a
    corner that does not move is never credited with what its moving neighbours swept: booking a third of each prism to
    every corner (sub-plan 19) credited the side walls of the box test with metal its receding front swept, they handed
    it down, bulged out and folded the surface at 1.3 s (sub-plan 20, 2026-10-12). `faces` are the outward
    triangles' node ids, `faces_corner` the same triangles in corner numbering."""
    a, b, c = (old[faces[:, i]] for i in range(3))
    A, B, C = (new[faces[:, i]] for i in range(3))

    def vol(p, q, r, s):
        return np.einsum("ij,ij->i", q - p, np.cross(r - p, s - p)) / 6.0
    prism = -(vol(a, b, c, A) + vol(b, c, A, B) + vol(c, A, B, C))
    n2 = np.cross(b - a, c - a)                                   # twice the area times the outward unit normal
    own = np.stack([-np.einsum("ij,ij->i", Q - P, n2) / 6.0 for P, Q in ((a, A), (b, B), (c, C))], axis=1)
    weight = np.abs(own)
    total = weight.sum(axis=1)
    share = np.where(total[:, None] > 0.0, weight / np.where(total > 0.0, total, 1.0)[:, None], 1.0 / 3.0)
    booked = own + (prism - own.sum(axis=1))[:, None] * share
    out = np.zeros(n_corners)
    for i in range(3):
        np.add.at(out, faces_corner[:, i], booked[:, i])
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
def default_sizes(mesh):
    """gmsh sizes for a rebuild, from the generator's parameters: a sphere's dense band, or a box's uniform h."""
    p = mesh.params
    if "h_surface_m" in p:
        return {"h_surface": p["h_surface_m"], "h_core": p["h_core_m"], "band": p.get("band_m", mesh_mod.DEFAULT_BAND),
                "ramp": max(p.get("band_ramp_m", mesh_mod.DEFAULT_BAND_RAMP), 1e-9)}
    h = p.get("h", mesh_mod.DEFAULT_H_SURFACE)
    return {"h_surface": h, "h_core": h, "band": 0.0, "ramp": 1.0}


@dataclass
class MoveResult:
    refused: bool                       # the move would invert or degrade an element: nothing was changed
    swept_volume: np.ndarray            # (n_corners,) m3 each corner's motion took out of the 3D model (negative: added)
    swept_enthalpy: np.ndarray          # (n_corners,) J that volume carried
    correction: float                   # J, the uniform offset that closed the books (pre-correction residual)
    worst_ratio: float                  # smallest radius ratio among the elements the move touched
    rebuilt: bool = False               # move_with_policy rebuilt the mesh before (or to make) this move


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
        self.rebuild_residual_total = 0.0       # J, summed pre-correction residuals of all rebuilds

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

    def neutral_smoothing(self, normal=None, area=None, passes=NEUTRAL_PASSES):
        """The tangential smoothing made volume-neutral per corner (layered-skin design section 5.1, decided 2026-10-12):
        each pass moves every corner along its normal by the volume the displacement so far sweeps there, divided by
        its dual area, so that a column is not credited with metal the smoothing alone swept. A corner's own normal
        move is booked to itself in full (own sweeps, to first order), so the passes converge quickly; what is left is
        booked to the columns like any other sweep (measured 2026-10-12 on the 100 mm sphere: the worst corner 0.46 um
        of equivalent depth before, 3e-11 m after three passes; the root mean square 9 nm before, 4e-13 m after)."""
        if normal is None or area is None:
            normal, area = self.frame()
        points = self.mesh.points[self.corners]
        u = self.smoothing(normal)
        for _ in range(passes):
            dV = swept_volumes(points, points + u, self.faces_corner, self.faces_corner, len(self.corners))
            u = u + normal * (dV / area)[:, None]
        return u

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
        # the slab is booked at the corners' old enthalpy, not the mean of the old and remapped values: the exact-solution
        # moving-front test is 0.0010 of the rise with the surface value against 0.0093 with the mean, and in stage 2
        # the metal crossing the rigid front is booked at the base node's temperature
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

    def rebuild(self, reason="requested"):
        """Regenerate the volume mesh inside the current surface triangulation and carry the enthalpy over."""
        t0 = time.perf_counter()
        mesh, solver, mat = self.mesh, self.solver, self.material
        faces = self.corners[self.faces_corner]
        check_surface(mesh.points, faces)
        points, tets, node_of_corner = rebuild_volume(mesh.points[self.corners], self.faces_corner, **self.sizes)
        h_old = mat.enthalpy(solver.temperature())
        E_old = solver.energy()
        h_new = project_enthalpy(mesh.points, mesh.tets, h_old, points, tets)
        params = dict(mesh.params, rebuilds=self.rebuilds + 1)
        new_mesh = mesh_mod.VolumeMesh(points, tets, params)
        check_rebuilt_boundary(new_mesh, node_of_corner[self.faces_corner])
        solver.setup(new_mesh, mat, solver.emissivity)
        solver.set_temperature(mat.temperature_from_enthalpy(h_new))
        residual = E_old - solver.energy()
        self.mesh, self.corners = new_mesh, node_of_corner
        self._set_surface(node_of_corner[self.faces_corner])             # builds the new deforming zone
        # what the quadrature misses arises where the field varies, near the surface, so the offset that makes the
        # total exact is spread over the deforming zone only, by lumped mass, and the deep body keeps its values
        in_zone = np.concatenate([self.zone.corners, self.zone.interior])
        m_zone = lumped_masses(points, tets, mat.rho)[in_zone]
        h_new[in_zone] += residual / m_zone.sum()
        solver.set_temperature(mat.temperature_from_enthalpy(h_new))
        self.rebuild_residual_total += float(residual)
        self.compression[:] = 0.0
        self.steps_since_rebuild = 0
        self.rebuilds += 1
        self.last_rebuild = RebuildResult(float(residual), float(residual / E_old) if E_old else 0.0, new_mesh.n_nodes,
                                          new_mesh.n_elements, time.perf_counter() - t0, reason)
        return self.last_rebuild

    def needs_rebuild(self):
        """The reason a rebuild is due by compression (the zone pressed to half its depth), or None."""
        if self.steps_since_rebuild >= self.min_gap and self.compression.max() > 0.5 * self.depth:
            return "compression"
        return None

    def move_with_policy(self, u_corners):
        """Rebuild first if compression calls for it; move; if the move is refused, rebuild and move again. A move
        refused on a fresh mesh is a model failure (RuntimeError, exit code 1)."""
        rebuilt = False
        reason = self.needs_rebuild()
        if reason:
            self.rebuild(reason)
            rebuilt = True
        result = self.move(u_corners)
        if result.refused:
            self.rebuild("quality")
            rebuilt = True
            result = self.move(u_corners)
            if result.refused:
                raise RuntimeError("a corner move of up to {:.3g} m degrades the mesh even after a rebuild".format(
                    float(np.linalg.norm(u_corners, axis=1).max())))
        result.rebuilt = rebuilt
        return result



# -- rebuilding the mesh ---------------------------------------------------------------------------------------------
@dataclass
class RebuildResult:
    residual: float                     # J, enthalpy the L2 projection's quadrature lost or gained before the offset
    relative_residual: float            # that over the total enthalpy
    n_nodes: int
    n_elements: int
    seconds: float
    reason: str


def check_surface(points, faces):
    """The validity gate: every edge used by exactly two triangles, one connected surface, positive enclosed volume."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    faces = np.asarray(faces, dtype=np.int64)
    edges = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    if np.any(counts != 2):
        raise RemeshError("the surface is not manifold: {} edges are not shared by exactly two triangles".format(
            int(np.sum(counts != 2))))
    nodes = np.unique(faces)
    local = np.full(points.shape[0], -1, dtype=np.int64)
    local[nodes] = np.arange(len(nodes))
    e = local[edges]
    graph = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(len(nodes), len(nodes)))
    n_parts = connected_components(graph, directed=False)[0]
    if n_parts != 1:
        raise RemeshError("the body has come apart into {} pieces; fragmentation is out of scope".format(n_parts))
    a, b, c = (points[faces[:, i]] for i in range(3))
    volume = float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)
    if volume <= 0.0:
        raise RemeshError("the surface encloses a non-positive volume ({:.3g} m3)".format(volume))
    return volume


def check_rebuilt_boundary(mesh, expected_faces):
    """Raise RemeshError unless the boundary triangles of `mesh` are exactly `expected_faces` (node ids, any order of
    triangles and of nodes within one): the rebuild must keep the surface it was given, with no face lost or added."""
    def triples(faces):
        return {tuple(f) for f in np.sort(np.asarray(faces, dtype=np.int64), axis=1)}
    got, want = triples(mesh.surface().faces), triples(expected_faces)
    if got != want:
        raise RemeshError("the rebuilt mesh's boundary differs from the surface it was given: {} triangles missing, "
                          "{} extra".format(len(want - got), len(got - want)))


def rebuild_volume(corner_points, faces_corner, h_surface, h_core, band, ramp):
    """Tetrahedra inside the closed triangulation (corner_points, faces_corner), keeping it exactly. Returns the new
    points, tets and the new node id of every corner (matched by coordinates: gmsh renumbers the nodes)."""
    import gmsh
    from scipy.spatial import cKDTree
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("rebuild")
        surface = gmsh.model.addDiscreteEntity(2)
        gmsh.model.mesh.addNodes(2, surface, list(range(1, len(corner_points) + 1)), corner_points.ravel().tolist())
        gmsh.model.mesh.addElementsByType(surface, 2, [], (np.asarray(faces_corner) + 1).ravel().tolist())
        loop = gmsh.model.geo.addSurfaceLoop([surface])
        gmsh.model.geo.addVolume([loop])
        gmsh.model.geo.synchronize()
        dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", [surface])
        thr = gmsh.model.mesh.field.add("Threshold")
        d_min = band
        d_max = band + ramp if band > 0.0 else max(ramp, 1e-9)
        for key, value in (("InField", dist), ("SizeMin", h_surface), ("SizeMax", h_core), ("DistMin", d_min),
                           ("DistMax", d_max)):
            gmsh.model.mesh.field.setNumber(thr, key, value)
        gmsh.model.mesh.field.setAsBackgroundMesh(thr)
        for option in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
            gmsh.option.setNumber(option, 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)          # HXT, as mesh.generate_sphere_mesh
        # one thread: fenicsx_env's gmsh 4.15.2 segfaults meshing a discrete surface's volume with HXT on four threads
        # (checked 2026-10-11; drama_env's does not), and a rebuild takes under two seconds on one
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.option.setNumber("Mesh.MaxNumThreads3D", 1)
        try:
            gmsh.model.mesh.generate(3)
            tags, coords, _ = gmsh.model.mesh.getNodes()
            types, _, nodes = gmsh.model.mesh.getElements(3)
        except Exception as error:
            raise RemeshError("gmsh could not mesh the surface: {}".format(error)) from error
    finally:
        gmsh.finalize()
    if 4 not in list(types):
        raise RemeshError("gmsh produced no tetrahedra inside the surface")
    coords = np.asarray(coords).reshape(-1, 3)
    tags = np.asarray(tags, dtype=np.int64)
    index = np.full(tags.max() + 1, -1, dtype=np.int64)
    index[tags] = np.arange(len(tags))
    tets = index[np.asarray(nodes[list(types).index(4)], dtype=np.int64).reshape(-1, 4)]
    used = np.unique(tets)                                   # drop nodes no tetrahedron uses (none expected)
    renumber = np.full(len(coords), -1, dtype=np.int64)
    renumber[used] = np.arange(len(used))
    coords, tets = coords[used], renumber[tets]
    d, node_of_corner = cKDTree(coords).query(corner_points)
    if d.max() > 0.0:
        raise RemeshError("gmsh moved a surface node by {:.3g} m".format(float(d.max())))
    return coords, tets, node_of_corner.astype(np.int64)


_QUAD_A, _QUAD_B = 0.5854101966249685, 0.1381966011250105      # the 4-point rule exact for quadratics on a tetrahedron
QUADRATURE = np.full((4, 4), _QUAD_B) + np.eye(4) * (_QUAD_A - _QUAD_B)   # barycentric coordinates, weights 1/4 each
PROJECT_CHUNK = 200_000                                      # quadrature points located at a time (bounds the memory)
LOCATE_RETRY = 256                                           # candidates for the points the first 16 missed: 0.41 % of the
                                                             # rebuild's quadrature points on the 4 mm test sphere (worst
                                                             # -0.149) and 0.28 % on the 100 mm sphere (worst -0.270); 256
                                                             # resolve all but 2 of the latter (measured 2026-10-12)


def project_enthalpy(old_points, old_tets, h_old, points, tets):
    """The nodal enthalpy on a new mesh as the L2 projection of the old P1 field: M h = b, M the new mesh's consistent
    mass matrix and b_i the integral of phi_i times the old field over the new elements by the 4-point rule. The
    projection reproduces a linear field exactly and keeps the integral of h -- the energy -- to the rule's accuracy,
    region by region, because the rows of M add up to the lumped masses the solvers' energy uses (sub-plan 20, the spec's
    fallback once the residuals proved large: interpolating at the new nodes and offsetting the total uniformly over the
    zone moved 2 768 J -- 4.9 % of the heat absorbed -- from the cold bulk to the melting front of the level-2 plate in two
    rebuilds, and its centre melted 9.5 % ahead of the exact solution; a lumped projection kept the energy local but
    moved a linear field's surface values by up to 7 K). An old field point just outside the old mesh -- a quadrature
    point across the shared surface -- is extrapolated linearly from the nearest element, as `locate` does."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    x = points[tets]                                         # (ne, 4, 3)
    vol = np.abs(signed_volumes(points, tets))
    qp = np.einsum("qa,eai->eqi", QUADRATURE, x).reshape(-1, 3)
    h_q = np.empty(len(qp))
    for start in range(0, len(qp), PROJECT_CHUNK):
        elem, bary = locate(old_points, old_tets, qp[start:start + PROJECT_CHUNK])
        h_q[start:start + PROJECT_CHUNK] = np.einsum("ij,ij->i", bary, h_old[old_tets[elem]])
    h_q = h_q.reshape(len(tets), 4)
    n = len(points)
    rhs = np.zeros(n)
    for a in range(4):                                       # integral of phi_a h: V/4 sum_q phi_a(q) h(q)
        np.add.at(rhs, tets[:, a], vol / 4.0 * (QUADRATURE[:, a][None, :] * h_q).sum(axis=1))
    local = np.full((4, 4), 1.0 / 20.0) + np.eye(4) / 20.0   # the P1 consistent mass on a tetrahedron, over its volume
    rows, cols = np.repeat(tets, 4, axis=1).ravel(), np.tile(tets, (1, 4)).ravel()
    M = sp.csr_matrix(((vol[:, None, None] * local[None]).ravel(), (rows, cols)), shape=(n, n))
    lumped = np.asarray(M.sum(axis=1)).ravel()
    h, info = spla.cg(M, rhs, x0=rhs / lumped, rtol=1e-13, maxiter=500, M=sp.diags(1.0 / lumped))
    if info != 0:
        raise RuntimeError("the rebuild's projection did not converge ({})".format(info))
    return h


def locate(points, tets, query, candidates=16):
    """For each query point, an element of (points, tets) holding it and its barycentric coordinates there; a point
    just outside every candidate gets the candidate it is least outside of (linear extrapolation). The `candidates`
    elements with the nearest centroids are tried first; the points whose best candidate still has a barycentric
    coordinate below -1e-9 -- in a graded mesh the containing element's centroid can lie beyond sixteen others -- are
    tried again among LOCATE_RETRY candidates from the same tree, keeping the better result per point."""
    from scipy.spatial import cKDTree
    tree = cKDTree(points[tets].mean(axis=1))

    def search(q, count):
        k = min(count, len(tets))
        _, cand = tree.query(q, k=k)
        cand = cand.reshape(len(q), k)
        x = points[tets[cand]]                               # (nq, k, 4, 3)
        J = np.stack([x[:, :, 1] - x[:, :, 0], x[:, :, 2] - x[:, :, 0], x[:, :, 3] - x[:, :, 0]], axis=3)
        lam = np.linalg.solve(J, (q[:, None, :] - x[:, :, 0])[..., None])[..., 0]
        bary = np.concatenate([1.0 - lam.sum(axis=2, keepdims=True), lam], axis=2)    # (nq, k, 4)
        best = np.argmax(bary.min(axis=2), axis=1)
        rows = np.arange(len(q))
        return cand[rows, best], bary[rows, best]

    elem, bary = search(query, candidates)
    missed = np.flatnonzero(bary.min(axis=1) < -1e-9)
    if len(missed) and candidates < LOCATE_RETRY:
        elem2, bary2 = search(query[missed], LOCATE_RETRY)
        better = bary2.min(axis=1) > bary[missed].min(axis=1)
        elem[missed[better]], bary[missed[better]] = elem2[better], bary2[better]
    return elem, bary

