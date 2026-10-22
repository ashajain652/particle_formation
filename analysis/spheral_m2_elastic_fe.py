"""Spheral M2 check 6, Step 1: the finite-element reference -- the 100 mm sphere, elastic, under frame 182's surface
loads and the balancing inertial body force (M2 plan, Task 8 Step 1; spec §7, §12 check 6).

    "$PY" analysis/spheral_m2_elastic_fe.py [--h-mm 5.0 2.5] [--frame 182] [--prepared DIR]

**The problem.** Linear elasticity on a ball of radius R = 50 mm centred at the origin, P2 (quadratic) displacement on
straight-sided gmsh tetrahedra, E = 71.7 GPa and nu = 0.33 (`MaterialTable.young(300)`, `.poisson(300)` of the
prepared flight; constant below the solidus) and rho = rho_free(300 K) = 2,811.62 kg/m^3 from the same table (not the
finite element's 2,813 kg/m^3: the table's free density at 300 K is what the particle body is built with). rho only
sets the deceleration; the stresses do not depend on it (the body force per volume is F / V).

**The tractions -- the same definition as `spheral_frag.runner.packages.loads`.** The body frame's flight direction
is +x: the nose is the +x end and the oncoming flow runs toward -x. At every surface quadrature point with outward
unit normal n (the facet's), theta = arccos(n_x) (`loads.inclination_deg`: 0 at the stagnation point, 180 deg at the
tail), (p, tau) = `spheral_frag.loads.evaluate(table, theta)` (np.interp over the non-empty 1-deg bins, held outside
them), and the traction is

    traction = -p n + tau t,     t = (-x_hat - (-x_hat . n) n) / |...|      (`loads.flow_tangent`; t = 0 on the axis)

so tau > 0 drags the surface toward the tail and the drag is D = integral traction . (-x_hat) dA
= integral (p cos theta + tau sin theta) dA (`spheral_frag.loads.drag`). The table is the prepared windward table of
the frame (`loads.read_tables`) with or without `spheral_frag.loads.with_lee` at its defaults (base 3 % of p_stag,
separation 180 deg, lee shear 0.5) -- the particle package always applies `with_lee` (`loads.lee_table`). Every
surface point is exposed (the undeformed convex sphere: the package's "first body" rule exposes every surface
particle of it too).

**The balancing body force.** a and alpha solve [integral rho (a + alpha x (x - x_c)) dV] = the traction's force
and moment, i.e. the load vector is made orthogonal to the six rigid modes exactly (in the P2 space the rigid
modes are interpolated exactly, so this is the continuous balance of the discrete loads); the body force is
-rho (a + alpha x (x - x_c)), x_c the mesh's centroid. a is the uniform deceleration of the particle runs'
`Deceleration` package (a = sum F / sum m; drag along -x gives a_x < 0 and every point is pushed toward the nose);
alpha, the angular part, only absorbs the faceted mesh's asymmetry and is reported (it vanishes on the exact
sphere, whose pressure has no moment).

**Rigid-body motion, two ways.** (A) "weak": the six rigid modes constrained by Lagrange multipliers,
integral rho u dV = 0 and integral rho (x - x_c) x u dV = 0, solved as (K + beta C C^T) u = f with C = M R (M the
consistent mass matrix, R the rigid modes), exactly equivalent when f is orthogonal to R; preconditioned CG
with pyamg smoothed aggregation (near-null space R). (B) "pinned": the 3-2-1 support -- the vertex nearest the
centre fixed, the one nearest (R/2, 0, 0) held in y and z, the one nearest (0, R/2, 0) held in z -- the condensed
nonsingular system by the same CG. The supports carry no reaction because the load is balanced, so the stresses
of A and B must agree to the solver's tolerance; the run reports their largest difference.

**Why gmsh, not `MeshTet.init_ball`.** gmsh's OCC sphere with a uniform target size gives any resolution (here 2.5
and 1.667 mm, ratio 1.5) with good elements throughout; `init_ball` only refines a fixed coarse ball, halving the size each
step with the projected boundary layer degrading. Straight-sided (subparametric) geometry is kept on purpose: on the
closed polyhedron a uniform pressure has exactly zero net force and moment and gives sigma = -p I to round-off
(the test), and the facets' normals lie within O(h / R) of radial (the surface area of the faceted ball is
reported against 4 pi R^2).

**Sampling.** sigma (xx, yy, zz, xy, yz, xz; Pa) and von Mises at the points of five meridians in the x-y plane,
theta = 0 (the nose half of the flight axis), 45, 90, 135 and 180 deg (the tail half) from the nose, at radii
r = 0, 1, ..., 49 and 49.5 mm (r = 49.5 mm, not 50, lies inside both faceted meshes). The stress is the exact
gradient of the P2 solution in the element containing the point (`stress_at`), not a smoothed field. The
convergence difference is the von Mises of sigma_fine - sigma_coarse at each point over the largest sampled von
Mises of the fine mesh, its maximum taken over the interior (r <= 45 mm) and the near-surface points (r > 45 mm)
separately, the latter also split at 47.5 mm (the outer band is the last one or two elements).

Outputs: `data/spheral/m2_check6_fe.csv` (committed; one row per mesh, load case, meridian and radius) and in
`spheral_output/m2/check6/fe/` the summary `summary.json` (drags, accelerations, residuals, the A-B difference, the
convergence, the cost) and per mesh and case `solution_<case>_h<h>mm.npz` (mesh points, tetrahedra and the P2
dofs; `load_solution` + `stress_at` sample it anywhere, for the comparison with the particles in Step 3).

**Measured** (2026-10-10, frame 182, meshes 2.5 and 1.667 mm: 131,543 and 439,120 tetrahedra, 0.56 and 1.82 M
dofs).

- *The lee extension changes nothing at frame 182*: the frame carries loads to 179.5 deg (`theta_last_deg`), so
  `with_lee` sets no bin and the "lee" and "windward" cases are the same loads and the same stresses (to the
  solver's tolerance). Both are written all the same.
- *Drag and deceleration.* D = integral traction . (-x_hat) dA = 30.710 N on the 1.667 mm mesh against the table's
  continuum integral on the exact sphere 30.717 N (-2.3e-4; -8.5e-3, -2.0e-3, -5.2e-4 at 10, 5, 2.5 mm: O(h^2),
  the faceted area being 2.3e-4 short of 4 pi R^2). a_x = -D / M = -20.869 m/s^2 with the mesh's mass (volume
  -4.1e-4 of the sphere's); with the exact sphere's M = 1.4722 kg the continuum value is -20.865 m/s^2. The
  history at k 182 (diagnostic only: its body has lost mass) has drag 49.58 N, mass 0.572 kg, deceleration
  86.69 m/s^2.
- *Balance.* Net force after the body force 3.5e-15 of D, moment 2.6e-17 of D R (round-off). Before it, the
  traction's moment about the centroid is 1.0e-6 D R and its transverse force 2.3e-5 D (the mesh's asymmetry),
  taken up by alpha ~ 1e-3 rad/s^2 and a_y, a_z ~ 4e-4 m/s^2.
- *Rigid modes.* Weak (Lagrange) and pinned (3-2-1) stresses agree to 1.6e-12 to 4.9e-12 of the largest von Mises;
  the pins' reaction is 2.2e-15 of D.
- *Stresses.* Largest sampled von Mises 3,510.5 Pa (on the axis, 15 mm toward the nose); centre sigma_xx
  -3,459 Pa, sigma_yy = sigma_zz -323 Pa (von Mises 3,137 Pa); nose surface (r 49.5 mm) near hydrostatic,
  sigma_xx -7,675, sigma_yy -7,814 Pa (von Mises 140 Pa); tail surface von Mises 190 Pa.
- *Convergence* (the von Mises of sigma_2.5 - sigma_1.667 over 3,510.5 Pa): interior (r <= 45 mm) max 0.41 %
  (RMS 0.08 %); 45-47.5 mm 0.86 %; the outer band r > 47.5 mm 2.2 % (pointwise stress within the last element or
  two). For comparison 5 -> 2.5 mm gave 1.6 % interior and 1.7 % near the surface, and 10 -> 5 mm 3.3 % and 5.0 %,
  so the 5 mm mesh is not used. Check 6 should compare the particles' interior against this reference at better
  than 0.5 % and treat the outer 2.5 mm as known to about 2 % (the particles' own surface shell, about 2 dx, is
  wider).
- *Cost* (this Mac, 18 cores, threaded BLAS): 2.5 mm assembly 21 s, the weak solve 50 s (204 CG
  iterations), the pinned solve 18 s (56); 1.667 mm assembly 71 s, weak 190 s (234), pinned 62 s (61). The whole
  default run (two meshes, two cases, both solves) takes 12.5 min.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import scipy.sparse as sp                                              # noqa: E402
import scipy.sparse.linalg as spla                                     # noqa: E402
import skfem as fe                                                     # noqa: E402
from skfem.helpers import dot                                          # noqa: E402
from skfem.models.elasticity import lame_parameters, linear_elasticity  # noqa: E402

from spheral_frag import loads as core_loads                          # noqa: E402
from spheral_frag.runner.packages import loads as pkg_loads           # noqa: E402

PREPARED = os.path.join(REPO, "spheral_output", "prepare",
                        "prep_model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_"
                        "material-AA7075_scheil_deeprunoff-off_k00000-01254")
OUT = os.path.join(REPO, "spheral_output", "m2", "check6", "fe")
CSV = os.path.join(REPO, "data", "spheral", "m2_check6_fe.csv")
RADIUS = 0.05                     # m
H_MM = (2.5, 1.6667)              # the two mesh sizes (coarse, fine), mm (5 -> 2.5 mm is not converged; docstring)
FRAME = 182                       # peak dynamic pressure (plan fact 8)
T_REF = 300.0                     # K
MERIDIANS_DEG = (0.0, 45.0, 90.0, 135.0, 180.0)
RADII_MM = tuple(float(r) for r in range(50)) + (49.5,)
INTERIOR_MM = 45.0                # r <= this is "interior" in the convergence measure
OUTER_MM = 47.5                   # r > this is the outer band (the last element or two) within "near surface"
CG_RTOL = 1e-12
FACET_INTORDER = 6                # the tables are piecewise linear in theta on 1-deg bins: a high order on facets
CASES = ("lee", "windward")       # with_lee at its defaults (what the particle package applies), and without
STRESS_KEYS = ("sxx", "syy", "szz", "sxy", "syz", "sxz")


# ---------------------------------------------------------------------------------------------------------------
# Mesh, material, tables
# ---------------------------------------------------------------------------------------------------------------

def ball_mesh(radius=RADIUS, h=0.005):
    """A tetrahedral ball (fe.MeshTet, straight-sided) from gmsh's OCC sphere at the uniform target size h (m)."""
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("ball")
        vol = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, radius)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMin", h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        for option in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature", "Mesh.MeshSizeExtendFromBoundary"):
            gmsh.option.setNumber(option, 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)          # HXT
        gmsh.model.mesh.generate(3)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        etypes, _, enodes = gmsh.model.mesh.getElements(3, vol)
        if list(etypes) != [4]:
            raise RuntimeError(f"gmsh gave element types {list(etypes)}, expected linear tetrahedra (4)")
    finally:
        gmsh.finalize()
    tags = np.asarray(tags, dtype=np.int64)
    p = np.asarray(coords, dtype=float).reshape(-1, 3)
    t_tags = np.asarray(enodes[0], dtype=np.int64).reshape(-1, 4)
    used, inv = np.unique(t_tags, return_inverse=True)
    lookup = np.full(tags.max() + 1, -1, dtype=np.int64)
    lookup[tags] = np.arange(len(tags))
    p = p[lookup[used]]
    t = inv.reshape(-1, 4)
    return fe.MeshTet(p.T.copy(), t.T.copy())


def material(prepared=PREPARED, T=T_REF):
    """(E Pa, nu, rho kg/m^3) at T from the prepared flight's material table."""
    from spheral_frag.material import MaterialTable
    tab = MaterialTable.load(os.path.join(prepared, "material_table.npz"))
    return float(tab.young(T)), float(tab.poisson(T)), float(tab.rho_free(T))


def frame_table(prepared=PREPARED, k=FRAME, lee=True):
    """Frame k's prepared windward table, with `with_lee` at its defaults when `lee` (the particle package's
    `lee_table(table, None)`), else as prepared."""
    table = pkg_loads.read_tables(prepared, [k])[0]
    return pkg_loads.lee_table(table, None) if lee else table


def table_from_arrays(arrays, i, lee=True):
    """Row i of a stacked loads file (e.g. the committed fixture `loads_mvp.npz`) as a table, optionally with the lee
    extension (for the tests, which must not need the git-ignored prepared flight)."""
    sub = {k: (arrays[k][[i]] if k not in ("theta_deg", "edges_deg") else arrays[k])
           for k in ("theta_deg", "edges_deg", "p", "tau", "area", "lee", "theta_last_deg", "p_stag")}
    table = core_loads.unstack(sub)[0]
    return pkg_loads.lee_table(table, None) if lee else table


def continuum_drag(table, radius=RADIUS, n=200001):
    """The table's drag on the exact sphere, 2 pi R^2 integral_0^pi (p cos theta + tau sin theta) sin theta dtheta (N),
    by the trapezoid rule on a fine theta grid (the reference the faceted surface integral converges to)."""
    th = np.linspace(0.0, math.pi, n)
    p, tau = core_loads.evaluate(table, np.degrees(th))
    f = (p*np.cos(th) + tau*np.sin(th))*np.sin(th)
    return float(2.0*math.pi*radius**2*np.sum(0.5*(f[1:] + f[:-1])*np.diff(th)))


# ---------------------------------------------------------------------------------------------------------------
# The discrete problem
# ---------------------------------------------------------------------------------------------------------------

class Problem:
    """P2 linear elasticity on `mesh`: stiffness K, consistent mass M (rho), rigid modes R (N x 6, about the mesh
    centroid x_c), and the boundary facet basis."""

    def __init__(self, mesh, E, nu, rho):
        self.mesh, self.E, self.nu, self.rho = mesh, float(E), float(nu), float(rho)
        self.lam, self.mu = lame_parameters(E, nu)
        self.elem = fe.ElementVector(fe.ElementTetP2())
        self.basis = fe.Basis(mesh, self.elem, intorder=2)
        self.fbasis = fe.FacetBasis(mesh, self.elem, intorder=FACET_INTORDER)
        self.K = fe.asm(linear_elasticity(self.lam, self.mu), self.basis).tocsr()
        rho_ = self.rho

        @fe.BilinearForm
        def mass(u, v, w):
            return rho_*dot(u, v)
        self.M = fe.asm(mass, self.basis).tocsr()
        vols = self.element_volumes()
        self.volume = float(vols.sum())
        cent = mesh.p[:, mesh.t].mean(axis=1)                      # (3, n_elements)
        self.x_c = (cent*vols).sum(axis=1)/self.volume
        self.mass_total = self.rho*self.volume
        self.R = rigid_modes(self.basis, self.x_c)
        self.MR = self.M @ self.R
        self.N = self.basis.N

    def element_volumes(self):
        p = self.mesh.p[:, self.mesh.t]                             # (3, 4, n)
        d = p[:, 1:, :] - p[:, :1, :]
        return np.abs(np.einsum("in,in->n", d[:, 0], np.cross(d[:, 1], d[:, 2], axis=0)))/6.0


def rigid_modes(basis, x_c):
    """The six rigid modes (translations x, y, z; rotations about x, y, z through x_c) at the vector P2 dofs
    (components interleaved: dof i is component i % 3)."""
    comp = np.arange(basis.N) % 3
    X = basis.doflocs - np.asarray(x_c, dtype=float)[:, None]
    R = np.zeros((basis.N, 6))
    for c in range(3):
        R[comp == c, c] = 1.0
    # rotation about axis a: u = e_a x X
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        u = np.cross(e, X.T).T                                      # (3, N)
        R[:, 3 + a] = u[comp, np.arange(basis.N)]
    return R


def traction_function(table):
    """traction(n) -> (3, ...) array: -p(theta) n + tau(theta) t with theta = arccos(n_x) (the package's definition)."""
    def traction(n):
        n = np.asarray(n, dtype=float)
        shape = n.shape[1:]
        nf = n.reshape(3, -1).T
        theta = pkg_loads.inclination_deg(nf)
        p, tau = core_loads.evaluate(table, theta)
        t = pkg_loads.flow_tangent(nf)
        tr = -p[:, None]*nf + tau[:, None]*t
        return tr.T.reshape((3,) + shape)
    return traction


def traction_load(prob, traction):
    """(load vector, integrals) for the surface traction `traction(n)`: integrals of the traction (N), its moment
    about x_c (N m), the drag integral traction . (-x_hat) dA (N) and the surface area (m^2)."""
    @fe.LinearForm
    def form(v, w):
        return dot(traction(w.n), v)
    f = fe.asm(form, prob.fbasis)

    @fe.Functional
    def area(w):
        return 1.0 + 0.0*w.x[0]
    comps = []
    for c in range(3):
        @fe.Functional
        def fc(w, c=c):
            return traction(w.n)[c]
        comps.append(fe.asm(fc, prob.fbasis))
    mom = []
    xc = prob.x_c
    for c in range(3):
        @fe.Functional
        def mc(w, c=c):
            tr = traction(w.n)
            X = [w.x[i] - xc[i] for i in range(3)]
            i, j = (c + 1) % 3, (c + 2) % 3
            return X[i]*tr[j] - X[j]*tr[i]
        mom.append(fe.asm(mc, prob.fbasis))
    force = np.array(comps)
    return f, {"force": force, "moment": np.array(mom), "drag": float(-force[0]),
               "area": float(fe.asm(area, prob.fbasis))}


def balance(prob, f):
    """The load vector with the inertial body force -rho (a + alpha x (x - x_c)) added so that it is orthogonal to
    the rigid modes: returns (f_balanced, a (3,), alpha (3,), residual force and moment before and after)."""
    G = prob.R.T @ prob.MR                                         # 6 x 6: integral rho r_i . r_j
    b = prob.R.T @ f                                               # [force; moment about x_c] of the tractions
    c = np.linalg.solve(G, b)
    fb = f - prob.MR @ c
    return fb, c[:3], c[3:], b, prob.R.T @ fb


def _amg(A, B):
    import pyamg
    ml = pyamg.smoothed_aggregation_solver(A, B=B, symmetry="hermitian", max_coarse=500)
    return ml.aspreconditioner(cycle="V")


def _cg(A, b, M):
    it = [0]

    def count(_):
        it[0] += 1
    x, info = spla.cg(A, b, M=M, rtol=CG_RTOL, atol=0.0, maxiter=5000, callback=count)
    if info != 0:
        raise RuntimeError(f"CG did not converge (info {info})")
    return x, it[0]


def solve_weak(prob, f):
    """(A) u with integral rho u = 0 and integral rho X x u = 0 (Lagrange multipliers): (K + beta C C^T) u = f,
    C = M R; f must be orthogonal to R."""
    C = prob.MR
    beta = float(prob.K.diagonal().mean())/float(np.mean(np.einsum("ij,ij->j", C, C)))*prob.N
    K = prob.K
    op = spla.LinearOperator(K.shape, matvec=lambda x: K @ x + beta*(C @ (C.T @ x)), dtype=float)
    u, it = _cg(op, f, _amg(K, prob.R))
    return u, it


def pinned_dofs(prob):
    """The 3-2-1 support's dofs: the vertex nearest the centre (x, y, z), the one nearest (R/2, 0, 0) (y, z), the one
    nearest (0, R/2, 0) (z). Vertex i's component c is dof basis.nodal_dofs[c, i]."""
    p = prob.mesh.p
    def nearest(q):
        return int(np.argmin(np.sum((p - np.asarray(q)[:, None])**2, axis=0)))
    r = RADIUS
    a, b, c = nearest((0, 0, 0)), nearest((r/2, 0, 0)), nearest((0, r/2, 0))
    nd = prob.basis.nodal_dofs
    return np.array([nd[0, a], nd[1, a], nd[2, a], nd[1, b], nd[2, b], nd[2, c]])


def solve_pinned(prob, f):
    """(B) the 3-2-1 support: the condensed system K_II u_I = f_I."""
    D = pinned_dofs(prob)
    I = np.setdiff1d(np.arange(prob.N), D)
    KI = prob.K[I][:, I].tocsr()
    u = np.zeros(prob.N)
    uI, it = _cg(KI, f[I], _amg(KI, prob.R[I]))
    u[I] = uI
    return u, it, prob.R.T @ (prob.K @ u - f)        # reactions' force and moment (should be ~0: balanced load)


# ---------------------------------------------------------------------------------------------------------------
# Stress sampling
# ---------------------------------------------------------------------------------------------------------------

def stress_at(basis, u, points, lam, mu):
    """sigma at `points` (n, 3), (n, 6) as xx, yy, zz, xy, yz, xz: the exact gradient of the P2 field u in the element
    containing each point."""
    x = np.asarray(points, dtype=float).T
    cells = basis.mesh.element_finder(mapping=basis.mapping)(*x)
    pts = basis.mapping.invF(x[:, :, None], tind=cells)
    grad = np.zeros((3, 3, x.shape[1]))
    dofs = basis.element_dofs[:, cells]                            # (Nbfun, n)
    for k in range(basis.Nbfun):
        g = basis.elem.gbasis(basis.mapping, pts, k, tind=cells)[0].grad[:, :, :, 0]   # (3, 3, n)
        grad += g*u[dofs[k]][None, None, :]
    eps = 0.5*(grad + grad.transpose(1, 0, 2))
    tr = eps[0, 0] + eps[1, 1] + eps[2, 2]
    s = 2.0*mu*eps
    for i in range(3):
        s[i, i] += lam*tr
    return np.stack([s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[1, 2], s[0, 2]], axis=1)


def von_mises(s):
    s = np.asarray(s, dtype=float)
    xx, yy, zz, xy, yz, xz = s.T
    return np.sqrt(0.5*((xx - yy)**2 + (yy - zz)**2 + (zz - xx)**2) + 3.0*(xy**2 + yz**2 + xz**2))


def sample_points(meridians_deg=MERIDIANS_DEG, radii_mm=RADII_MM):
    """[(meridian_deg, r_mm, point (3,) m)]: the meridians in the x-y plane, theta from the nose (+x)."""
    out = []
    for th in meridians_deg:
        c, s = math.cos(math.radians(th)), math.sin(math.radians(th))
        for r in radii_mm:
            out.append((th, r, np.array([r*c, r*s, 0.0])*1e-3))
    return out


def save_solution(path, prob, u, meta):
    np.savez_compressed(path, p=prob.mesh.p, t=prob.mesh.t, u=u, lam=prob.lam, mu=prob.mu,
                        meta=json.dumps(meta))


def load_solution(path):
    """(basis, u, lam, mu, meta) of a saved solution; `stress_at(basis, u, points, lam, mu)` samples it."""
    with np.load(path, allow_pickle=False) as z:
        mesh = fe.MeshTet(z["p"], z["t"])
        basis = fe.Basis(mesh, fe.ElementVector(fe.ElementTetP2()), intorder=2)
        return basis, z["u"], float(z["lam"]), float(z["mu"]), json.loads(str(z["meta"]))


# ---------------------------------------------------------------------------------------------------------------
# The study
# ---------------------------------------------------------------------------------------------------------------

def run(h_mm=H_MM, frame=FRAME, prepared=PREPARED, out=OUT, csv_path=CSV):
    os.makedirs(out, exist_ok=True)
    E, nu, rho = material(prepared)
    tables = {"lee": frame_table(prepared, frame, lee=True), "windward": frame_table(prepared, frame, lee=False)}
    with np.load(os.path.join(prepared, "flight.npz")) as z:        # one row per frame (row k = frame k)
        hist = {k: float(z[k][frame]) for k in ("time_s", "drag_N", "deceleration_ms2", "mass_kg", "p_w_stag_step_Pa")}
    M_sphere = rho*4.0/3.0*math.pi*RADIUS**3
    summary = {"frame": frame, "prepared": os.path.relpath(prepared, REPO), "E_Pa": E, "nu": nu, "rho_kgm3": rho,
               "radius_m": RADIUS, "mass_exact_sphere_kg": M_sphere, "history_diagnostic": hist,
               "tables": {c: {"theta_last_deg": t.theta_last_deg, "p_stag_Pa": t.p_stag,
                              "lee_bins": int(np.count_nonzero(t.lee)),
                              "continuum_drag_N": continuum_drag(t),
                              "continuum_a_x_ms2": -continuum_drag(t)/M_sphere}
                          for c, t in tables.items()},
               "lee_changes_table": bool(not np.array_equal(tables["lee"].p, tables["windward"].p, equal_nan=True)
                                         or not np.array_equal(tables["lee"].tau, tables["windward"].tau,
                                                               equal_nan=True)),
               "meshes": {}}
    pts = sample_points()
    P = np.array([q for _, _, q in pts])
    samples = {}
    for h in h_mm:
        t0 = time.perf_counter()
        mesh = ball_mesh(RADIUS, h*1e-3)
        t_mesh = time.perf_counter() - t0
        t0 = time.perf_counter()
        prob = Problem(mesh, E, nu, rho)
        t_asm = time.perf_counter() - t0
        info = {"h_mm": h, "n_vertices": int(mesh.p.shape[1]), "n_tets": int(mesh.t.shape[1]), "n_dofs": int(prob.N),
                "volume_m3": prob.volume, "volume_rel_sphere": prob.volume/(4.0/3.0*math.pi*RADIUS**3) - 1.0,
                "mass_kg": prob.mass_total, "x_c_m": prob.x_c.tolist(), "t_mesh_s": t_mesh, "t_assemble_s": t_asm,
                "cases": {}}
        print(f"h {h} mm: {info['n_tets']} tets, {prob.N} dofs (mesh {t_mesh:.1f} s, assembly {t_asm:.1f} s)",
              flush=True)
        for case, table in tables.items():
            t0 = time.perf_counter()
            f, I = traction_load(prob, traction_function(table))
            fb, a, alpha, before, after = balance(prob, f)
            t_load = time.perf_counter() - t0
            t0 = time.perf_counter()
            uA, itA = solve_weak(prob, fb)
            t_A = time.perf_counter() - t0
            t0 = time.perf_counter()
            uB, itB, reaction = solve_pinned(prob, fb)
            t_B = time.perf_counter() - t0
            sA = stress_at(prob.basis, uA, P, prob.lam, prob.mu)
            sB = stress_at(prob.basis, uB, P, prob.lam, prob.mu)
            vmA = von_mises(sA)
            FR = abs(I["drag"])*RADIUS
            info["cases"][case] = {
                "surface_area_m2": I["area"], "area_rel_sphere": I["area"]/(4.0*math.pi*RADIUS**2) - 1.0,
                "traction_force_N": I["force"].tolist(), "traction_moment_Nm": I["moment"].tolist(),
                "drag_N": I["drag"], "drag_rel_continuum": I["drag"]/summary["tables"][case]["continuum_drag_N"] - 1,
                "a_ms2": a.tolist(), "alpha_rad_s2": alpha.tolist(),
                "a_x_exact_sphere_mass_ms2": -I["drag"]/M_sphere,
                "residual_force_rel_drag": float(np.linalg.norm(after[:3])/abs(I["drag"])),
                "residual_moment_rel_drag_R": float(np.linalg.norm(after[3:])/FR),
                "unbalanced_moment_rel_drag_R": float(np.linalg.norm(before[3:])/FR),
                "unbalanced_transverse_force_rel_drag": float(np.linalg.norm(before[1:3])/abs(I["drag"])),
                "pinned_reaction_rel_drag": float(np.linalg.norm(reaction[:3])/abs(I["drag"])),
                "weak_vs_pinned_max_rel_vm": float(np.max(von_mises(sA - sB))/np.max(vmA)),
                "max_vm_sampled_Pa": float(np.max(vmA)), "cg_iters_weak": int(itA), "cg_iters_pinned": int(itB),
                "t_load_s": t_load, "t_weak_s": t_A, "t_pinned_s": t_B}
            samples[(h, case)] = sA
            save_solution(os.path.join(out, f"solution_{case}_h{h:g}mm.npz"), prob, uA,
                          {"case": case, "frame": frame, "h_mm": h, "a_ms2": a.tolist(), "alpha": alpha.tolist(),
                           "convention": "x flight axis, nose +x; traction -p n + tau t"})
            c = info["cases"][case]
            print(f"  {case}: drag {I['drag']:.5f} N (continuum {c['drag_rel_continuum']:+.2e}), a_x {a[0]:.4f} m/s^2,"
                  f" residual F {c['residual_force_rel_drag']:.1e} M {c['residual_moment_rel_drag_R']:.1e},"
                  f" weak-pinned {c['weak_vs_pinned_max_rel_vm']:.1e}, max vM {c['max_vm_sampled_Pa']:.1f} Pa,"
                  f" CG {itA}/{itB} its, {t_load:.1f}+{t_A:.1f}+{t_B:.1f} s", flush=True)
        summary["meshes"][f"{h:g}"] = info
    # convergence: coarse against fine
    r = np.array([q[1] for q in pts])
    hc, hf = max(h_mm), min(h_mm)
    conv = {}
    for case in tables:
        sF, sC = samples[(hf, case)], samples[(hc, case)]
        vmax = float(np.max(von_mises(sF)))
        d = von_mises(sF - sC)/vmax
        conv[case] = {"max_vm_fine_Pa": vmax, "interior_max_rel": float(d[r <= INTERIOR_MM].max()),
                      "near_surface_max_rel": float(d[r > INTERIOR_MM].max()),
                      "band_45_47p5_max_rel": float(d[(r > INTERIOR_MM) & (r <= OUTER_MM)].max()),
                      "outer_max_rel": float(d[r > OUTER_MM].max()),
                      "interior_rms_rel": float(np.sqrt(np.mean(d[r <= INTERIOR_MM]**2)))}
        print(f"convergence {case}: interior {conv[case]['interior_max_rel']:.2e}, near-surface "
              f"{conv[case]['near_surface_max_rel']:.2e} of max vM {vmax:.1f} Pa")
    summary["convergence"] = conv
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["h_mm", "case", "meridian_deg", "r_mm", "x_mm", "y_mm", "z_mm", *[k + "_Pa" for k in STRESS_KEYS],
                    "von_mises_Pa", "dvm_coarse_fine_Pa"])
        for h in h_mm:
            for case in tables:
                s = samples[(h, case)]
                vm = von_mises(s)
                dv = von_mises(samples[(hf, case)] - samples[(hc, case)])
                for i, (th, rr, q) in enumerate(pts):
                    w.writerow([f"{h:g}", case, f"{th:g}", f"{rr:g}", *[f"{v*1e3:.4f}" for v in q],
                                *[f"{v:.6g}" for v in s[i]], f"{vm[i]:.6g}", f"{dv[i]:.4g}"])
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--h-mm", type=float, nargs=2, default=list(H_MM), help="coarse and fine mesh sizes (mm)")
    ap.add_argument("--frame", type=int, default=FRAME)
    ap.add_argument("--prepared", default=PREPARED)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--csv", default=CSV)
    args = ap.parse_args(argv)
    if not os.path.isdir(args.prepared):
        print(f"no prepared flight at {args.prepared}", file=sys.stderr)
        return 2
    t0 = time.perf_counter()
    run(tuple(args.h_mm), args.frame, args.prepared, args.out, args.csv)
    print(f"total {time.perf_counter() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
