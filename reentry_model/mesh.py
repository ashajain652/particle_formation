"""Tetrahedral sphere meshes (gmsh, graded, optional prism layers) and their surface geometry. Spec sections 5
(Step 2) and 5 (Step 3).

A VolumeMesh is node coordinates, tetrahedra, an active mask (Step 3: elements that have melted away are
deactivated and drop out of the surface, the operators and the mass) and the boundary triangles of the active set
(faces used by exactly one active element) with outward unit normals, centroids, areas and owner elements -- so any
gmsh volume mesh loads, tagged or not. A face table (every face of the mesh with its one or two elements) is built
once so that deactivation updates the boundary in O(faces) without re-sorting.

`sphere_mesh(R, h_surface, h_core)` is the Step 2 generator: element size held at h_surface throughout an outer
`band` (amendment 2026-09-27: the dense band, default 15 mm, that survives a receding surface) and then graded to
h_core over a `ramp` (gmsh Distance/Threshold field), a node embedded at the centre, cached by
(R, h_surface, h_core, band, ramp). With `layers = n > 0` (off by default now that the dense band is the recession
strategy) the graded gmsh mesh is generated for the inner sphere R - T (T the total layer thickness) and n prism
layers are built on it by radial projection of its boundary triangulation: shells at R - T + t_(n-1) + ... ,
thicknesses t0 growth^k from the outside in (the outermost layer is t0 thick), each prism split into three
tetrahedra by the smallest-node-id diagonal rule (conforming across prisms). The construction is exact on a
sphere and needs no gmsh boundary-layer machinery (spike of 2026-09-20:
gmsh's extrudeBoundaryLayer is a geo-kernel function that cannot be attached to the OCC sphere without
re-parametrising; the radial construction gives the same layers deterministically). `box_mesh` is a structured
Kuhn-split box for the analytic (Stefan) tests.
"""
import os
from dataclasses import dataclass, field

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MESH_DIR = os.path.join(REPO_ROOT, "reentry_model_output", "meshes")
DEFAULT_H_SURFACE = 2.0e-3          # m; 2 mm / 8 mm gives 18.9 k nodes on the 100 mm sphere at band = 0
                                    # (the Step 2 field, which every caller here pins); 33.4 k with the
                                    # 15 mm default band (measured 2026-09-18 and 2026-09-27)
                                    # sphere (measured 2026-09-27; band = 0 reproduces the pre-band 18.9 k instead)
DEFAULT_H_CORE = 8.0e-3             # m
DEFAULT_LAYERS = 0                  # prism layers: superseded as the default by the dense band (amendment 2026-09-27)
DEFAULT_LAYER_THICKNESS = 0.25e-3   # m, outermost layer
DEFAULT_LAYER_GROWTH = 2.0
DEFAULT_BAND = 15.0e-3              # m; cells stay at h_surface this far down, so the surface survives 7.5 mm of
                                    # recession before the mesh must be rebuilt. 310 k tets on the 100 mm sphere
DEFAULT_BAND_RAMP = 8.0e-3          # m; the depth over which the size then grows from h_surface to h_core
TAUBIN_PASSES = 8                   # lambda/mu smoothing of the surface read by the heating (amendment 2026-09-27)
TAUBIN_LAMBDA = 0.53
TAUBIN_MU = -0.55                   # slightly stronger than lambda: the pair passes low frequencies without shrinking
TAUBIN_CLAMP = 0.35                 # of the mean patch edge; a node can then never move past its neighbours
_FACE_OF_VERTEX = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])     # face opposite each tet vertex


@dataclass
class SurfaceMesh:
    faces: np.ndarray            # (nf, 3) node ids of the boundary triangles ("patches")
    centroids: np.ndarray        # (nf, 3) m
    normals: np.ndarray          # (nf, 3) outward unit normals
    areas: np.ndarray            # (nf,) m2
    owner: np.ndarray = None     # (nf,) element owning each patch
    face_ids: np.ndarray = None  # (nf,) ids in the mesh's face table (stable across deactivations)
    points: np.ndarray = None    # the (n, 3) coordinates these patches index -- `points` or `derived_points`; this is
                                  # the one field that ALIASES the mesh's array rather than snapshotting it, so a
                                  # write to `derived_points` after this SurfaceMesh was built moves `points` under
                                  # it while `centroids`, `normals` and `areas` stay stale -- call `surface()` again

    @property
    def n_patches(self):
        return len(self.faces)

    @property
    def area(self):
        return float(self.areas.sum())

    def angles_to(self, v_hat):
        """theta per patch [rad]: the angle between the outward normal and v_hat, the direction the body moves in.
        The stagnation patch has its normal along v_hat (theta = 0); theta > pi/2 is leeward."""
        v = np.asarray(v_hat, dtype=float)
        v = v / np.linalg.norm(v)
        return np.arccos(np.clip(self.normals @ v, -1.0, 1.0))

    def facet_mean(self, nodal):
        """Mean of a nodal field over each patch's three nodes."""
        return np.asarray(nodal)[self.faces].mean(axis=1)

    def edge_lengths(self):
        """The three edge lengths of every patch, (nf, 3) m. `edge_lengths().max()` is the surface coarseness the
        history records and the second remesh trigger watches: on a fresh 100 mm sphere it is 2.16 mm, and once the
        erosion front is into the graded core it reaches 8 mm (measured fact 39)."""
        x = self.points[self.faces]
        return np.linalg.norm(x[:, [1, 2, 0]] - x, axis=2)

    def patch_toward(self, direction):
        """Index of the patch whose outward normal is closest to `direction`."""
        d = np.asarray(direction, dtype=float)
        return int(np.argmax(self.normals @ (d / np.linalg.norm(d))))

    def tangent_from(self, v_hat):
        """Unit surface direction away from the stagnation point (increasing theta) per patch; zero where undefined."""
        v = np.asarray(v_hat, dtype=float)
        v = v / np.linalg.norm(v)
        t = -v[None, :] + (self.normals @ v)[:, None] * self.normals
        n = np.linalg.norm(t, axis=1)
        return np.where(n[:, None] > 1e-12, t / np.maximum(n, 1e-300)[:, None], 0.0)

    def projected_area(self, v_hat):
        """Area projected on the plane normal to v_hat: sum of A max(0, n.v) (pi R^2 for a sphere)."""
        v = np.asarray(v_hat, dtype=float)
        return float((self.areas * np.maximum(0.0, self.normals @ (v / np.linalg.norm(v)))).sum())

    def edges(self):
        """Patch adjacency: (edge node pairs (ne, 2), patch i, patch j, edge lengths) for edges shared by two patches."""
        e = np.sort(self.faces[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
        patch = np.repeat(np.arange(self.n_patches), 3)
        order = np.lexsort((e[:, 1], e[:, 0]))
        e, patch = e[order], patch[order]
        same = np.all(e[1:] == e[:-1], axis=1)
        k = np.flatnonzero(same)
        return e[k], patch[k], patch[k + 1]

    def smoothed_normals(self, passes=TAUBIN_PASSES, lam=TAUBIN_LAMBDA, mu=TAUBIN_MU, clamp=TAUBIN_CLAMP):
        """Patch normals of a Taubin-smoothed copy of this surface. Nothing is moved: the copy is thrown away.

        Element death leaves a staircase whose steps are one cell tall and one cell wide, so its normals deviate
        from the true surface by up to 45 degrees -- and theta, the angle between the normal and the flight
        direction, is what the Lees distribution and the film tangents read. Taubin's pair (a shrinking pass at
        +lam, an inflating one at mu slightly more negative) removes the cell-scale roughness without the net
        shrinkage a plain Laplacian would cause. It looks only at a node and its neighbours, so it assumes nothing
        about the shape: any genus, any number of components, concave or convex. Each displacement is clamped to
        `clamp` of the mean patch edge, so a node can never pass its neighbours."""
        import math
        import scipy.sparse as sp                    # imported here, like gmsh and meshio
        node = np.unique(self.faces)
        local = np.full(int(self.faces.max()) + 1, -1, dtype=np.int64)
        local[node] = np.arange(node.size)
        pairs = np.sort(local[self.faces][:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
        pairs = np.unique(pairs, axis=0)
        pairs = pairs[pairs[:, 0] != pairs[:, 1]]                       # a collapsed patch has a self-edge
        w = sp.coo_matrix((np.ones(2 * len(pairs)), (np.r_[pairs[:, 0], pairs[:, 1]],
                                                     np.r_[pairs[:, 1], pairs[:, 0]])),
                          shape=(node.size, node.size)).tocsr()
        raw_degree = np.asarray(w.sum(axis=1)).ravel()
        isolated = raw_degree == 0.0                                    # no neighbours: nothing to average toward
        degree = np.maximum(raw_degree, 1.0)                            # only to keep the division finite below
        p = self.points[node].copy()
        step = clamp * math.sqrt(4.0 * self.areas.mean() / math.sqrt(3.0))    # the mean patch's own edge
        for _ in range(passes):
            for weight in (lam, mu):
                delta = weight * (w @ p / degree[:, None] - p)
                delta[isolated] = 0.0                                   # an isolated node stays exactly where it is
                length = np.linalg.norm(delta, axis=1)
                over = length > step
                delta[over] *= (step / length[over])[:, None]
                p = p + delta
        moved = self.points.copy()
        moved[node] = p
        a, b, c = (moved[self.faces[:, i]] for i in range(3))
        n = np.cross(b - a, c - a)
        length = np.linalg.norm(n, axis=1)
        flat = length < 1e-30                                           # a patch the recession collapsed
        n = np.where(flat[:, None], self.normals, n / np.where(flat, 1.0, length)[:, None])
        n[np.einsum("ij,ij->i", n, self.normals) < 0.0] *= -1.0         # keep this surface's own orientation
        return n


@dataclass
class VolumeMesh:
    points: np.ndarray                                   # (n, 3) m
    tets: np.ndarray                                     # (ne, 4) node ids
    params: dict = field(default_factory=dict)           # generator parameters / source file
    element_layer: np.ndarray = None                     # (ne,) prism layer index (0 = outermost) or -1 (graded core)
    active: np.ndarray = field(init=False, repr=False, default=None)
    _face_nodes: np.ndarray = field(init=False, repr=False, default=None)    # (nfaces, 3) node ids, one entry per unique face
    _face_elements: np.ndarray = field(init=False, repr=False, default=None) # (nfaces, 2) elements sharing the face, -1 if one
    _face_opposite: np.ndarray = field(init=False, repr=False, default=None) # (nfaces, 2) opposite vertex in each of those elements
    _face_count: np.ndarray = field(init=False, repr=False, default=None)    # active elements per face
    _element_faces: np.ndarray = field(init=False, repr=False, default=None) # (ne, 4) face ids of each element
    _boundary: np.ndarray = field(init=False, repr=False, default=None)      # face ids of the current boundary
    derived_points: np.ndarray = field(init=False, repr=False, default=None)   # recession-displaced copy of points

    def __post_init__(self):
        self.points = np.array(self.points, dtype=float)
        self.tets = np.asarray(self.tets, dtype=np.int64)
        self.active = np.ones(len(self.tets), dtype=bool)
        if self.element_layer is None:
            self.element_layer = np.full(len(self.tets), -1, dtype=np.int64)
        self._build_face_table()
        self.derived_points = self.points.copy()

    def _build_face_table(self):
        faces = self.tets[:, _FACE_OF_VERTEX].reshape(-1, 3)                    # 4 per element, element-major
        opposite = self.tets.ravel()                                            # vertex opposite face 4e + v is tets[e, v]
        element = np.repeat(np.arange(len(self.tets)), 4)
        key = np.sort(faces, axis=1)
        _, first, inverse, counts = np.unique(key, axis=0, return_index=True, return_inverse=True, return_counts=True)
        inverse = inverse.ravel()
        nfaces = len(first)
        self._face_nodes = faces[first]
        self._face_elements = np.full((nfaces, 2), -1, dtype=np.int64)
        self._face_opposite = np.full((nfaces, 2), -1, dtype=np.int64)
        order = np.argsort(inverse, kind="stable")
        slot = np.zeros(nfaces, dtype=np.int64)
        for k in order:                                                         # <= 2 entries per face
            f = inverse[k]
            self._face_elements[f, slot[f]], self._face_opposite[f, slot[f]] = element[k], opposite[k]
            slot[f] += 1
        self._face_count = counts.astype(np.int64)
        self._element_faces = inverse.reshape(len(self.tets), 4)              # face id of each element's four faces
        self._boundary = np.flatnonzero(self._face_count == 1)

    @property
    def n_nodes(self):
        return len(self.points)

    @property
    def n_elements(self):
        return len(self.tets)

    @property
    def n_active(self):
        return int(self.active.sum())

    def element_volumes(self):
        x = self.points[self.tets]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        return np.abs(np.linalg.det(J)) / 6.0

    def volume(self):
        """Volume of the active elements."""
        return float(self.element_volumes()[self.active].sum())

    def element_faces(self, elements):
        """Face-table ids of each given element's four faces, (n, 4). Stable across deactivations, so a dying element
        can be asked which of its own faces have just become boundary (Step 3: where its melt film goes)."""
        return self._element_faces[np.asarray(elements, dtype=np.int64)]

    def deactivate(self, elements):
        """Remove elements from the active set; the boundary (and `surface()`) follows. Returns the face ids that
        stopped being boundary faces (their patches vanished) and those that became boundary faces."""
        elements = np.unique(np.asarray(elements, dtype=np.int64))
        elements = elements[self.active[elements]]
        if elements.size == 0:
            return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
        self.active[elements] = False
        faces = self._element_faces[elements].ravel()
        before = self._face_count[faces].copy()
        np.subtract.at(self._face_count, faces, 1)
        after = self._face_count[faces]
        gone, new = np.unique(faces[(before == 1) & (after == 0)]), np.unique(faces[(before == 2) & (after == 1)])
        self._boundary = np.flatnonzero(self._face_count == 1)
        return gone, new

    def surface(self, derived=False):
        """The boundary patches of the active set with their current geometry.

        With `derived=True` the geometry comes from `derived_points`, the recession-displaced copy, instead of from
        `points`. That is the whole of Design B: the shape every consumer of the surface reads -- the heating angle
        theta, the film tangents, the projected area, the nose fit, the pictures -- moves continuously with the
        recession, while `points` never moves, so the thermal solver's element matrices stay valid and skfem keeps
        rescaling them at 0.004 s instead of recomputing at 0.074 s (measured fact 45)."""
        f = self._boundary
        faces = self._face_nodes[f]
        pts = self.derived_points if derived else self.points
        active_slot = np.where(self.active[self._face_elements[f, 0]], 0, 1)
        owner = self._face_elements[f, active_slot]
        opposite = self._face_opposite[f, active_slot]
        a, b, c = (pts[faces[:, i]] for i in range(3))
        n = np.cross(b - a, c - a)
        areas = 0.5 * np.linalg.norm(n, axis=1)
        flat = areas <= 0.0                                  # a derived displacement can collapse a patch
        n = n / np.where(flat, 1.0, 2.0 * areas)[:, None]    # never divide by zero; `flat` rows are fixed below
        centroids = (a + b + c) / 3.0
        away = centroids - pts[opposite]
        if flat.any():                                       # outward by construction, and finite
            n[flat] = away[flat] / np.maximum(np.linalg.norm(away[flat], axis=1), 1e-300)[:, None]
        inward = np.einsum("ij,ij->i", n, away) < 0.0
        n[inward] *= -1.0
        return SurfaceMesh(faces, centroids, n, areas, owner, f.copy(), pts)

    def boundary_nodes(self):
        return np.unique(self._face_nodes[self._boundary])

    def active_nodes(self):
        """Nodes belonging to at least one active element."""
        return np.unique(self.tets[self.active])

    def element_depths(self):
        """Distance from each active element's centroid to the nearest boundary patch [m]; -1 for dead elements.
        The first remesh trigger fires when the deepest element that still owns a patch is more than half the band
        down, i.e. when the erosion front has eaten through half the fine mesh. Measured against boundary *centroids*
        rather than the triangles themselves, which over-reads by at most half a patch edge -- immaterial against a
        15 mm band, and it keeps this a single cKDTree query."""
        from scipy.spatial import cKDTree            # imported here, like gmsh and meshio: not every mesh needs it
        depths = np.full(len(self.tets), -1.0)
        live = np.flatnonzero(self.active)
        patches = self.surface()
        if live.size == 0 or patches.n_patches == 0:             # demised: no elements, or none with a boundary
            return depths
        centres = self.points[self.tets[live]].mean(axis=1)
        depths[live] = cKDTree(patches.centroids).query(centres)[0]
        return depths

    def centre_node(self):
        """Node nearest the origin (the embedded centre node of a generated sphere)."""
        return int(np.argmin(np.linalg.norm(self.points, axis=1)))


def boundary_faces(tets):
    """Faces used by exactly one tetrahedron, as (nf, 3) node ids, and the opposite vertex of that tetrahedron."""
    faces = tets[:, _FACE_OF_VERTEX].reshape(-1, 3)
    opposite = np.repeat(tets, 4, axis=0)[np.arange(4 * len(tets)), np.tile(np.arange(4), len(tets))]
    _, first, counts = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_counts=True)
    keep = first[counts == 1]
    return faces[keep], opposite[keep]


def mesh_file_name(radius, h_surface, h_core, band=DEFAULT_BAND, ramp=DEFAULT_BAND_RAMP):
    """Cache key. The band and its ramp are part of it: two meshes of the same radius and sizes but different bands
    are different meshes, and sharing one file would silently serve the wrong one."""
    return "sphere_R{:.3f}mm_hs{:.3f}mm_hc{:.3f}mm_b{:.3f}mm_r{:.3f}mm.msh".format(
        radius * 1e3, h_surface * 1e3, h_core * 1e3, band * 1e3, ramp * 1e3)


def band_limits(band, ramp, radius):
    """The Threshold field's DistMin and DistMax for a dense band of `band` m on a sphere of `radius`. The band is
    clamped to the radius -- a band at or past the centre just means uniformly fine -- and DistMax is held strictly
    above DistMin, because Threshold divides by the difference. Returned rather than inlined so that `sphere_mesh`
    can record the band the mesh actually got: it differs from the one asked for on a remnant smaller than the band,
    which is exactly where the remesh trigger reads it.

    `band = 0` is the Step 2 field -- graded from the wall to the centre over the whole radius -- and must stay
    exactly that, because the committed verification numbers were measured on it: letting DistMax collapse to the
    ramp instead gives 24 413 tetrahedra where the committed mesh has 87 632 (measured 2026-09-27)."""
    if band <= 0.0:
        return 0.0, radius
    d_min = min(band, radius)
    return d_min, max(min(band + ramp, radius), d_min + 1e-9)


def generate_sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR,
                         band=DEFAULT_BAND, ramp=DEFAULT_BAND_RAMP):
    """gmsh sphere of `radius`, at h_surface throughout the outer `band` and then graded to h_core over `ramp`;
    returns the cached .msh path. `band = 0` is the Step 2 grading, which the pinned verification devices use."""
    os.makedirs(mesh_dir, exist_ok=True)
    path = os.path.join(mesh_dir, mesh_file_name(radius, h_surface, h_core, band, ramp))
    if os.path.isfile(path):
        return path
    import gmsh                                   # imported here: only mesh generation needs gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("sphere")
        vol = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, radius)
        centre = gmsh.model.occ.addPoint(0.0, 0.0, 0.0, h_core)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.embed(0, [centre], 3, vol)
        surfaces = [s[1] for s in gmsh.model.getBoundary([(3, vol)], oriented=False)]
        gmsh.model.addPhysicalGroup(3, [vol], tag=1, name="body")
        gmsh.model.addPhysicalGroup(2, surfaces, tag=2, name="surface")
        dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", surfaces)
        gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)
        thr = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thr, "InField", dist)
        gmsh.model.mesh.field.setNumber(thr, "SizeMin", h_surface)
        gmsh.model.mesh.field.setNumber(thr, "SizeMax", h_core)
        d_min, d_max = band_limits(band, ramp, radius)
        gmsh.model.mesh.field.setNumber(thr, "DistMin", d_min)
        gmsh.model.mesh.field.setNumber(thr, "DistMax", d_max)
        gmsh.model.mesh.field.setAsBackgroundMesh(thr)
        for option in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
            gmsh.option.setNumber(option, 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)          # HXT
        gmsh.model.mesh.generate(3)
        tmp = path[:-4] + ".tmp.msh"           # gmsh.write needs the .msh extension to pick the format
        gmsh.write(tmp)
    finally:
        gmsh.finalize()
    os.replace(tmp, path)                      # atomic: an interrupted generation never leaves a reusable partial .msh
    return path


def load_mesh(path):
    """Any gmsh/meshio volume mesh with tetrahedra (other cell types are ignored)."""
    import meshio
    m = meshio.read(path)
    tets = [c.data for c in m.cells if c.type == "tetra"]
    if not tets:
        raise ValueError("no tetrahedra in {}".format(path))
    return VolumeMesh(np.asarray(m.points, dtype=float), np.vstack(tets), {"path": os.path.abspath(path)})


def layer_thicknesses(layers, layer_thickness, growth):
    """Thickness of each prism layer from the outside in: t0, t0 g, t0 g^2, ..."""
    return layer_thickness * growth ** np.arange(layers)


def split_prisms(bottom, top):
    """Three tetrahedra per prism (bottom (np, 3) and top (np, 3) node ids, top[k] above bottom[k]) with the diagonal
    of every quad face drawn from its smallest node id -- so two prisms sharing a quad face split it identically."""
    bottom, top = np.asarray(bottom, dtype=np.int64), np.asarray(top, dtype=np.int64)
    # rotate each triangle so that its smallest bottom id comes first (the top follows the same rotation)
    shift = np.argmin(bottom, axis=1)
    idx = (shift[:, None] + np.arange(3)[None, :]) % 3
    b = np.take_along_axis(bottom, idx, axis=1)
    t = np.take_along_axis(top, idx, axis=1)
    b0, b1, b2, t0, t1, t2 = b[:, 0], b[:, 1], b[:, 2], t[:, 0], t[:, 1], t[:, 2]
    # diagonals from b0 on the faces (b0 b1 t1 t0) and (b2 b0 t0 t2); on (b1 b2 t2 t1) from min(b1, b2)
    case = b1 < b2
    tets = np.where(case[:, None, None],
                    np.stack([np.stack([b0, b1, b2, t2], 1), np.stack([b0, b1, t2, t1], 1), np.stack([b0, t1, t2, t0], 1)], 1),
                    np.stack([np.stack([b0, b1, b2, t1], 1), np.stack([b0, b2, t1, t2], 1), np.stack([b0, t1, t2, t0], 1)], 1))
    return tets.reshape(-1, 4)


def add_prism_layers(inner, radius, layers, layer_thickness, growth):
    """Prism layers by radial projection of the inner sphere mesh's boundary triangulation out to `radius`
    (thicknesses t0 growth^k from the outside in). Returns the layered VolumeMesh (inner elements first, then the
    layers from the inside out; `element_layer` is 0 for the outermost layer, -1 for the graded core)."""
    t = layer_thicknesses(layers, layer_thickness, growth)
    r_inner = radius - t.sum()
    faces, _ = boundary_faces(inner.tets)
    shell_nodes = np.unique(faces)
    ns = len(shell_nodes)
    local = -np.ones(inner.n_nodes, dtype=np.int64)
    local[shell_nodes] = np.arange(ns)
    lf = local[faces]                                        # faces in shell-local ids
    direction = inner.points[shell_nodes]
    direction = direction / np.linalg.norm(direction, axis=1)[:, None]
    radii = r_inner + np.cumsum(t[::-1])                     # shell radii from the inside out (outermost = radius)
    shell_ids = [shell_nodes] + [inner.n_nodes + k * ns + np.arange(ns) for k in range(layers)]
    points = np.vstack([inner.points] + [direction * r for r in radii])
    tets, layer_index = [inner.tets], [np.full(inner.n_elements, -1, dtype=np.int64)]
    for k in range(layers):
        new = split_prisms(shell_ids[k][lf], shell_ids[k + 1][lf])
        tets.append(new)
        layer_index.append(np.full(len(new), layers - 1 - k, dtype=np.int64))
    return VolumeMesh(points, np.vstack(tets), dict(inner.params), np.concatenate(layer_index))


def sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR, layers=0,
                layer_thickness=DEFAULT_LAYER_THICKNESS, growth=DEFAULT_LAYER_GROWTH, band=DEFAULT_BAND,
                ramp=DEFAULT_BAND_RAMP):
    """Graded gmsh sphere with a dense outer band, optionally with `layers` prism layers on top of it."""
    if layers < 0 or layer_thickness <= 0.0 or growth <= 0.0:
        raise ValueError("layers must be >= 0 and layer_thickness, growth > 0")
    if band < 0.0 or ramp <= 0.0:
        raise ValueError("band must be >= 0 and ramp > 0")
    total = float(layer_thicknesses(layers, layer_thickness, growth).sum()) if layers else 0.0
    if total >= 0.5 * radius:
        raise ValueError("the prism layers ({:.3g} m) must be thinner than half the radius".format(total))
    mesh = load_mesh(generate_sphere_mesh(radius - total, h_surface, h_core, mesh_dir, band, ramp))
    if layers:
        mesh = add_prism_layers(mesh, radius, layers, layer_thickness, growth)
    d_min, d_max = band_limits(band, ramp, radius - total)
    mesh.params.update({"radius_m": radius, "h_surface_m": h_surface, "h_core_m": h_core, "layers": layers,
                        "layer_thickness_m": layer_thickness, "layer_growth": growth, "band_m": d_min,
                        "band_ramp_m": d_max - d_min})   # the effective ramp: a nanometre once the band is clamped to the radius
    return mesh


def box_mesh(lx, ly, lz, h):
    """Structured box [0, lx] x [0, ly] x [0, lz] of cubes of side ~h, each split into six tetrahedra (Kuhn
    decomposition, conforming). For the analytic tests (Stefan front along x with the x = 0 face heated)."""
    n = [max(1, int(round(L / h))) for L in (lx, ly, lz)]
    xs = [np.linspace(0.0, L, k + 1) for L, k in zip((lx, ly, lz), n)]
    X, Y, Z = np.meshgrid(*xs, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    nid = np.arange(points.shape[0]).reshape(n[0] + 1, n[1] + 1, n[2] + 1)
    c = nid[:-1, :-1, :-1].ravel()
    dx, dy, dz = nid[1, 0, 0] - nid[0, 0, 0], nid[0, 1, 0] - nid[0, 0, 0], 1
    v = lambda i, j, k: c + i * dx + j * dy + k * dz
    corners = [v(0, 0, 0), v(1, 0, 0), v(1, 1, 0), v(0, 1, 0), v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)]
    kuhn = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
    tets = np.vstack([np.column_stack([corners[a], corners[b], corners[cc], corners[d]]) for a, b, cc, d in kuhn])
    return VolumeMesh(points, tets, {"box": (lx, ly, lz), "h": h, "path": None})
