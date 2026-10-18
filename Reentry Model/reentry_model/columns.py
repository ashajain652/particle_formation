"""The corner columns of the layered-skin design (docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md
section 5.1): one column on every corner of the 3D mesh's surface triangulation, owning the corner's dual cell.

The triangulation's connectivity is fixed for the whole flight, so everything here that depends only on it (the
edges, the two triangles on each edge, the corner opposite each edge in each triangle) is built once; `update` recomputes
the geometry when the corners move. The dual cell is barycentric: a third of each touching triangle, bounded by the
segments from each touching edge's midpoint to the triangles' centroids, so the areas add up to the surface's exactly.

Sideways conduction between neighbouring columns uses the cotangent weight of their shared edge, w = (cot a + cot b)/2
with a and b the angles opposite the edge in its two triangles: the flux of a piecewise-linear field through the two
dual segments across the edge is exactly w (T_j - T_i) per unit conductivity and depth, so the discrete operator is
exact for a field linear along the surface. A weight is clipped at zero where the opposite angles are obtuse enough
to make it negative: a negative conductance would break the explicit step's positivity (and, since sub-plan 21 makes the
sideways conduction implicit, the definiteness its conjugate gradients need), and it is rare (3 of 6 756 edges on the
4 mm test sphere, 4 of 28 245 on the default 2 mm sphere, measured 2026-10-12).
The surface gradient and curvature operators the flow needs are added in sub-plan 21, the first stage that uses them."""
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from .movemesh import corner_frame


class Columns:
    """Corner columns on a triangulation (points (nc, 3), faces (nt, 3) in corner numbering, oriented outward)."""

    def __init__(self, points, faces):
        self.faces = np.asarray(faces, dtype=np.int64)
        self.n = int(np.asarray(points).shape[0])
        e = np.vstack([self.faces[:, [0, 1]], self.faces[:, [1, 2]], self.faces[:, [2, 0]]])
        opposite = np.concatenate([self.faces[:, 2], self.faces[:, 0], self.faces[:, 1]])
        triangle = np.tile(np.arange(len(self.faces)), 3)
        key = np.sort(e, axis=1)
        self.edges, inverse = np.unique(key, axis=0, return_inverse=True)
        inverse = inverse.ravel()
        counts = np.bincount(inverse, minlength=len(self.edges))
        if np.any(counts > 2):
            raise ValueError("an edge is shared by more than two triangles: the surface is not a manifold")
        # each edge's (up to) two triangles and the corner opposite it in each; -1 on an open boundary
        order = np.argsort(inverse, kind="stable")
        first = np.searchsorted(inverse[order], np.arange(len(self.edges)))
        self._tri = np.full((len(self.edges), 2), -1, dtype=np.int64)
        self._opp = np.full((len(self.edges), 2), -1, dtype=np.int64)
        self._tri[:, 0], self._opp[:, 0] = triangle[order][first], opposite[order][first]
        two = counts == 2
        self._tri[two, 1], self._opp[two, 1] = triangle[order][first[two] + 1], opposite[order][first[two] + 1]
        self.degree = np.bincount(self.edges.ravel(), minlength=self.n)
        self.update(points)

    def update(self, points):
        """Recompute the geometry for corners at `points` (the connectivity is unchanged)."""
        p = np.asarray(points, dtype=float)
        self.points = p
        self.normal, self.area = corner_frame(p, self.faces, np.arange(self.n))
        self.centroid = p[self.faces].mean(axis=1)
        i, j = self.edges[:, 0], self.edges[:, 1]
        self.edge_length = np.linalg.norm(p[j] - p[i], axis=1)
        mid = 0.5 * (p[i] + p[j])
        cot = np.zeros(len(self.edges))
        face = np.zeros(len(self.edges))
        for s in range(2):
            ok = self._tri[:, s] >= 0
            k = self._opp[ok, s]
            u, v = p[i[ok]] - p[k], p[j[ok]] - p[k]
            cot[ok] += np.einsum("ij,ij->i", u, v) / np.linalg.norm(np.cross(u, v), axis=1)
            face[ok] += np.linalg.norm(mid[ok] - self.centroid[self._tri[ok, s]], axis=1)
        self.cot_weight = 0.5 * cot
        self.conduct_weight = np.maximum(self.cot_weight, 0.0)
        self.face_length = face

    def top_points(self, thickness):
        """The top surface: each corner moved out along its normal by its column's thickness [m]."""
        return self.points + np.asarray(thickness, dtype=float)[:, None] * self.normal

    def laplacian(self, f, weights=None):
        """sum over neighbours j of w_ij (f_j - f_i) per corner (f per corner, any trailing shape)."""
        w = self.cot_weight if weights is None else weights
        f = np.asarray(f, dtype=float)
        i, j = self.edges[:, 0], self.edges[:, 1]
        d = (f[j] - f[i]) * w.reshape((-1,) + (1,) * (f.ndim - 1))
        out = np.zeros_like(f)
        np.add.at(out, i, d)
        np.add.at(out, j, -d)
        return out

    @property
    def incidence(self):
        """B (n columns x E edges), +1 at each edge's first column and -1 at its second, so B @ f gives every column the
        sum over its edges of f, with the first column's sign. Built once: the connectivity is fixed."""
        if getattr(self, "_incidence", None) is None:
            i, j = self.edges[:, 0], self.edges[:, 1]
            e = np.arange(len(i))
            self._incidence = sp.csr_matrix((np.r_[np.ones(len(i)), -np.ones(len(i))], (np.r_[i, j], np.r_[e, e])),
                                            shape=(self.n, len(i)))
        return self._incidence


@dataclass
class CornerSurface:
    """The gas-facing surface of a corner-column body as the coupled run reads it: one 'patch' per corner (its top dual
    cell), with outward unit normals and areas, and the triangles for drawing it."""
    points: np.ndarray            # (nc, 3) the top surface's corners
    faces: np.ndarray             # (nt, 3) in corner numbering
    normals: np.ndarray           # (nc, 3)
    areas: np.ndarray             # (nc,) m2

    @classmethod
    def from_points(cls, points, faces):
        points, faces = np.asarray(points, dtype=float), np.asarray(faces, dtype=np.int64)
        normal, area = corner_frame(points, faces, np.arange(len(points)))
        return cls(points, faces, normal, area)

    @property
    def n_patches(self):
        return len(self.areas)

    @property
    def area(self):
        return float(self.areas.sum())

    def angles_to(self, v_hat):
        v = np.asarray(v_hat, dtype=float)
        return np.arccos(np.clip(self.normals @ (v / np.linalg.norm(v)), -1.0, 1.0))

    def patch_toward(self, direction):
        d = np.asarray(direction, dtype=float)
        return int(np.argmax(self.normals @ (d / np.linalg.norm(d))))

    def projected_area(self, v_hat):
        v = np.asarray(v_hat, dtype=float)
        return float((self.areas * np.clip(self.normals @ (v / np.linalg.norm(v)), 0.0, None)).sum())
