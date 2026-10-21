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

The flow (sub-plan 21) reads two surface operators on the same graph. `gradient` averages the gradients of a field's
linear interpolation over each corner's triangles, by area, and projects the result onto the corner's tangent plane.
`curvature` is the mean curvature (the sum of the principal curvatures) from the cotangent Laplacian of the position
over each corner's mixed Voronoi area (Meyer et al. 2003). That area is exact only where no touching triangle is
obtuse, so a corner touching one takes the mean of its neighbours that touch none (33 of the 4 mm test sphere's 2 254
corners and 48 of the 2 mm sphere's 9 417, which the mixed area alone left off by 15 and 3 % on average and by up to 176
and 23 %). On these 50 mm-radius spheres the curvature is then within 8.4·10^-5 of 2/R at every corner (1.3·10^-5 on
the 2 mm sphere), where the barycentric area had scattered 7.6 and 3.9 %, and the gradient of a linear field within
0.09 % of its tangential part on average (0.02 % on the 2 mm sphere), 0.8 % at the test sphere's one irregular corner
(measured 2026-10-18)."""
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

    def gradient(self, f, points=None):
        """The tangential gradient at every corner of a per-corner field f (n,) [f per m] on the surface through
        `points` (default: the columns' own corners): the gradients of f's linear interpolation on the touching
        triangles, averaged by area, projected onto the corner's tangent plane. Exact for a field linear on a flat
        surface."""
        p = self.points if points is None else np.asarray(points, dtype=float)
        normal = self.normal if points is None else corner_frame(p, self.faces, np.arange(self.n))[0]
        A, unit, _ = triangle_geometry(p, self.faces)
        a, b, c = (p[self.faces[:, v]] for v in range(3))
        fa, fb, fc = (np.asarray(f, dtype=float)[self.faces[:, v]] for v in range(3))
        g = (fa[:, None] * np.cross(unit, c - b) + fb[:, None] * np.cross(unit, a - c)
             + fc[:, None] * np.cross(unit, b - a)) / (2.0 * A[:, None])
        G, W = np.zeros((self.n, 3)), np.zeros(self.n)
        for v in range(3):
            np.add.at(G, self.faces[:, v], g * A[:, None])
            np.add.at(W, self.faces[:, v], A)
        G /= W[:, None]
        return G - np.einsum("ij,ij->i", G, normal)[:, None] * normal

    def curvature(self, points=None):
        """The mean curvature, the sum of the principal curvatures [1/m], at every corner of the surface through
        `points` (default: the columns' own corners), positive where the surface bends away from its outward normal (2/R
        on a sphere): the cotangent Laplacian of the position along the corner's normal over its mixed Voronoi area. A
        corner touching an obtuse triangle takes the mean of its neighbours that touch none, or, where every neighbour
        touches one, of those that took a value in the pass before (module docstring)."""
        p = self.points if points is None else np.asarray(points, dtype=float)
        normal = self.normal if points is None else corner_frame(p, self.faces, np.arange(self.n))[0]
        A, _, cot = triangle_geometry(p, self.faces)
        lap, area = np.zeros((self.n, 3)), np.zeros(self.n)
        obtuse = np.zeros(self.n, dtype=bool)
        blunt = (cot < 0.0).any(axis=1)
        for v in range(3):
            i, j, k = self.faces[:, v], self.faces[:, (v + 1) % 3], self.faces[:, (v + 2) % 3]
            cj, ck = cot[:, (v + 1) % 3], cot[:, (v + 2) % 3]
            dj, dk = p[j] - p[i], p[k] - p[i]
            np.add.at(lap, i, 0.5 * (ck[:, None] * dj + cj[:, None] * dk))
            voronoi = 0.125 * (np.einsum("ij,ij->i", dj, dj) * ck + np.einsum("ij,ij->i", dk, dk) * cj)
            np.add.at(area, i, np.where(cot[:, v] < 0.0, 0.5 * A, np.where(blunt, 0.25 * A, voronoi)))
            np.logical_or.at(obtuse, i, blunt)
        kappa = -np.einsum("ij,ij->i", lap, normal) / area
        i, j = self.edges[:, 0], self.edges[:, 1]
        known = ~obtuse
        while not known.all():                     # outward from the corners whose own value is exact, a ring a pass
            good = known.astype(float)
            total, count = np.zeros(self.n), np.zeros(self.n)
            np.add.at(total, i, kappa[j] * good[j])
            np.add.at(total, j, kappa[i] * good[i])
            np.add.at(count, i, good[j])
            np.add.at(count, j, good[i])
            fill = ~known & (count > 0.0)
            if not fill.any():
                break
            kappa[fill] = total[fill] / count[fill]
            known = known | fill
        return kappa

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


def triangle_geometry(points, faces):
    """Per triangle: its area [m2], its unit normal (outward for outward-oriented faces) and the cotangent of the angle
    at each of its three corners, (nt, 3)."""
    a, b, c = (points[faces[:, v]] for v in range(3))
    n2 = np.cross(b - a, c - a)
    twice = np.linalg.norm(n2, axis=1)
    legs = ((b - a, c - a), (c - b, a - b), (a - c, b - c))
    cot = np.stack([np.einsum("ij,ij->i", u, w) / np.linalg.norm(np.cross(u, w), axis=1) for u, w in legs], axis=1)
    return 0.5 * twice, n2 / twice[:, None], cot


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
