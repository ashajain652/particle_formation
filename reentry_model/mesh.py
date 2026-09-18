"""Tetrahedral sphere meshes (gmsh, graded) and their surface geometry. Spec section 5.

A VolumeMesh is node coordinates (a mutable array, so Step 3 can recede the surface), tetrahedra, and the boundary
triangles derived from the tetrahedra themselves (faces used by exactly one element) with outward unit normals,
centroids and areas -- so any gmsh volume mesh loads, tagged or not. The generator makes a sphere whose element
size grows linearly from h_surface at the surface to h_core at the centre (gmsh Distance/Threshold field) with a
node embedded at the centre, and caches the .msh by (R, h_surface, h_core).
"""
import os
from dataclasses import dataclass, field

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MESH_DIR = os.path.join(REPO_ROOT, "reentry_model_output", "meshes")
DEFAULT_H_SURFACE = 2.0e-3          # m; 2 mm / 8 mm gives 18.9 k nodes on the 100 mm sphere (measured 2026-09-18)
DEFAULT_H_CORE = 8.0e-3             # m
_FACE_OF_VERTEX = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])     # face opposite each tet vertex


@dataclass
class SurfaceMesh:
    faces: np.ndarray            # (nf, 3) node ids of the boundary triangles ("patches")
    centroids: np.ndarray        # (nf, 3) m
    normals: np.ndarray          # (nf, 3) outward unit normals
    areas: np.ndarray            # (nf,) m2

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

    def patch_toward(self, direction):
        """Index of the patch whose outward normal is closest to `direction`."""
        d = np.asarray(direction, dtype=float)
        return int(np.argmax(self.normals @ (d / np.linalg.norm(d))))


@dataclass
class VolumeMesh:
    points: np.ndarray                                   # (n, 3) m, mutable
    tets: np.ndarray                                     # (ne, 4) node ids
    params: dict = field(default_factory=dict)           # generator parameters / source file
    _faces: np.ndarray = field(init=False, repr=False, default=None)
    _opposite: np.ndarray = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self.points = np.array(self.points, dtype=float)
        self.tets = np.asarray(self.tets, dtype=np.int64)
        self._faces, self._opposite = boundary_faces(self.tets)

    @property
    def n_nodes(self):
        return len(self.points)

    @property
    def n_elements(self):
        return len(self.tets)

    def element_volumes(self):
        x = self.points[self.tets]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        return np.abs(np.linalg.det(J)) / 6.0

    def volume(self):
        return float(self.element_volumes().sum())

    def surface(self):
        """The boundary patches with their current geometry (recomputed from `points` on every call)."""
        a, b, c = (self.points[self._faces[:, i]] for i in range(3))
        n = np.cross(b - a, c - a)
        areas = 0.5 * np.linalg.norm(n, axis=1)
        n = n / (2.0 * areas)[:, None]
        centroids = (a + b + c) / 3.0
        inward = np.einsum("ij,ij->i", n, centroids - self.points[self._opposite]) < 0.0
        n[inward] *= -1.0
        return SurfaceMesh(self._faces, centroids, n, areas)

    def boundary_nodes(self):
        return np.unique(self._faces)

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


def mesh_file_name(radius, h_surface, h_core):
    return "sphere_R{:.3f}mm_hs{:.3f}mm_hc{:.3f}mm.msh".format(radius * 1e3, h_surface * 1e3, h_core * 1e3)


def generate_sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR):
    """gmsh sphere of `radius` graded from h_surface (surface) to h_core (centre); returns the cached .msh path."""
    os.makedirs(mesh_dir, exist_ok=True)
    path = os.path.join(mesh_dir, mesh_file_name(radius, h_surface, h_core))
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
        gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(thr, "DistMax", radius)
        gmsh.model.mesh.field.setAsBackgroundMesh(thr)
        for option in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
            gmsh.option.setNumber(option, 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)          # HXT
        gmsh.model.mesh.generate(3)
        gmsh.write(path)
    finally:
        gmsh.finalize()
    return path


def load_mesh(path):
    """Any gmsh/meshio volume mesh with tetrahedra (other cell types are ignored)."""
    import meshio
    m = meshio.read(path)
    tets = [c.data for c in m.cells if c.type == "tetra"]
    if not tets:
        raise ValueError("no tetrahedra in {}".format(path))
    return VolumeMesh(np.asarray(m.points, dtype=float), np.vstack(tets), {"path": os.path.abspath(path)})


def sphere_mesh(radius, h_surface=DEFAULT_H_SURFACE, h_core=DEFAULT_H_CORE, mesh_dir=DEFAULT_MESH_DIR):
    mesh = load_mesh(generate_sphere_mesh(radius, h_surface, h_core, mesh_dir))
    mesh.params.update({"radius_m": radius, "h_surface_m": h_surface, "h_core_m": h_core})
    return mesh
