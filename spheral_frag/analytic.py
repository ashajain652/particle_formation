"""Exact solutions the runner is checked against (Spheral M2 plan, Task 6; spec §12 check 4).

**Thermal stress in a solid elastic sphere** (Timoshenko and Goodier, *Theory of Elasticity*, 3rd ed., §150): a free
sphere of radius R with a radial temperature T(r), constant E, nu and alpha, stress free at a uniform T_ref. With
tau(r) = T(r) - T_ref and the mean J(r) = (1/r^3) int_0^r tau s^2 ds = int_0^1 tau(r u) u^2 du,

    sigma_r     = 2 alpha E / (1 - nu) [ J(R) - J(r) ]
    sigma_theta =   alpha E / (1 - nu) [ 2 J(R) + J(r) - tau(r) ]      (= sigma_phi)

so sigma_r(R) = 0 and, at the centre (J(0) = tau(0)/3), sigma_r = sigma_theta. `sphere_thermal_stress` evaluates J
by composite Gauss-Legendre quadrature on [0, 1] in the scaled variable u (no division by r^3, so r = 0 is regular);
for a smooth profile it is exact to round-off (the default 4 panels of 32 points integrate a polynomial of degree 63
exactly; measured against the closed form below to 1e-15 of the peak, tests/test_spheral_frag_analytic.py).

**The quadratic profile** tau = dT (r/R)^2 (check 4's case): J(r) = dT r^2 / (5 R^2), hence

    sigma_r     = sigma_0 (1 - rho^2),   sigma_theta = sigma_0 (1 - 2 rho^2),   sigma_0 = 2 alpha E dT / (5 (1 - nu))

with rho = r / R: tension sigma_0 at the centre (both components), sigma_theta = -sigma_0 at the surface, so the
largest |sigma| is sigma_0 at both places. For AA7075's table values (E 71.7 GPa, nu 0.33, alpha 23.4e-6 1/K) and
dT = 50 K, sigma_0 = 50.08 MPa. (The M2 plan's "75 MPa, sigma_theta at the surface" is 0.6 alpha E dT / (1 - nu):
not this solution; corrected at Task 6, 2026-10-10.)

**Spherical components** of a stress tensor at a point x relative to the centre: sigma_rr = n . sigma . n with
n = x/|x|, and the tangential mean sigma_tt = (tr sigma - sigma_rr) / 2, which is sigma_theta = sigma_phi for a
spherically symmetric field. In RZ the tensor's rows are (z, r, theta) (axial, cylindrical radius, hoop) and
x = (z, r): n = (z, r, 0)/|x|.

Numpy only (the core's import boundary)."""
from __future__ import annotations

import math

import numpy as np

GAUSS_POINTS = 32
GAUSS_PANELS = 4


def _gauss_unit(points=GAUSS_POINTS, panels=GAUSS_PANELS):
    """Nodes and weights of composite Gauss-Legendre quadrature on [0, 1]."""
    x, w = np.polynomial.legendre.leggauss(points)
    edges = np.linspace(0.0, 1.0, panels + 1)
    nodes, weights = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        nodes.append(0.5*(b - a)*x + 0.5*(a + b))
        weights.append(0.5*(b - a)*w)
    return np.concatenate(nodes), np.concatenate(weights)


def mean_excess(r, T_of_r, T_ref, points=GAUSS_POINTS, panels=GAUSS_PANELS) -> np.ndarray:
    """J(r) = (1/r^3) int_0^r (T(s) - T_ref) s^2 ds = int_0^1 (T(r u) - T_ref) u^2 du (module docstring); J(0) =
    (T(0) - T_ref)/3. `T_of_r` takes an array of radii (m) and returns T (K)."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    u, w = _gauss_unit(points, panels)
    s = r[:, None]*u[None, :]
    tau = np.asarray(T_of_r(s.ravel()), dtype=float).reshape(s.shape) - float(T_ref)
    return (tau*(w*u*u)[None, :]).sum(axis=1)


def sphere_thermal_stress(r, T_of_r, R, E, nu, alpha, T_ref, points=GAUSS_POINTS, panels=GAUSS_PANELS):
    """(sigma_r, sigma_theta) in Pa at radii r (m) of a free elastic sphere of radius R (m) with temperature
    T_of_r(r) (K), stress free at the uniform T_ref (Timoshenko and Goodier; module docstring)."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    if np.any(r < 0.0) or np.any(r > R*(1.0 + 1e-12)):
        raise ValueError("radii must lie in [0, R]")
    c = float(alpha)*float(E)/(1.0 - float(nu))
    JR = float(mean_excess(np.array([R]), T_of_r, T_ref, points, panels)[0])
    J = mean_excess(r, T_of_r, T_ref, points, panels)
    tau = np.asarray(T_of_r(r), dtype=float) - float(T_ref)
    return 2.0*c*(JR - J), c*(2.0*JR + J - tau)


def quadratic_peak(E, nu, alpha, dT) -> float:
    """sigma_0 = 2 alpha E dT / (5 (1 - nu)), Pa: the largest |sigma| of the quadratic profile."""
    return 2.0*float(alpha)*float(E)*float(dT)/(5.0*(1.0 - float(nu)))


def sphere_quadratic_stress(r, R, E, nu, alpha, dT):
    """(sigma_r, sigma_theta), Pa, for T = T_ref + dT (r/R)^2: sigma_0 (1 - rho^2), sigma_0 (1 - 2 rho^2)."""
    rho2 = (np.asarray(r, dtype=float)/float(R))**2
    s0 = quadratic_peak(E, nu, alpha, dT)
    return s0*(1.0 - rho2), s0*(1.0 - 2.0*rho2)


def spherical_components(x, sigma, centre=None):
    """(r, sigma_rr, sigma_tt): the distance from `centre` and the radial and mean tangential normal stresses of the
    (n, 3, 3) tensors `sigma` at positions x ((n, 3) in 3D; (n, 2) = (z, r) in RZ, tensor rows (z, r, theta)).
    At r = 0 the radial direction is taken as the first axis."""
    x = np.asarray(x, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    if x.shape[1] == 2:
        p = np.zeros((len(x), 3))
        p[:, :2] = x
        if centre is not None:
            p[:, 0] -= float(np.asarray(centre, dtype=float).ravel()[0])
    else:
        p = x - (0.0 if centre is None else np.asarray(centre, dtype=float))
    r = np.sqrt(np.einsum("ij,ij->i", p, p))
    n = np.zeros_like(p)
    ok = r > 0.0
    n[ok] = p[ok]/r[ok, None]
    n[~ok, 0] = 1.0
    srr = np.einsum("ni,nij,nj->n", n, sigma, n)
    tr = np.einsum("nii->n", sigma)
    return r, srr, 0.5*(tr - srr)


def expansion_strain_alpha(alpha, T_base, T_ref_expansion):
    """The linear expansion coefficient referred to T_base of a free density rho_s / (1 + alpha (T - T_ref_exp))^3:
    the free length scales as 1 + alpha (T - T_ref_exp), so the strain from T_base is alpha (T - T_base) /
    (1 + alpha (T_base - T_ref_exp)), exactly linear in T (the material table's rho_free, material.mechanical_columns,
    up to its 1 K interpolation)."""
    return float(alpha)/(1.0 + float(alpha)*(float(T_base) - float(T_ref_expansion)))


def binned_profile(r, values, weights, edges):
    """Weighted means of each column of `values` ((n,) or (n, m)) in radial bins: (centres, means (nbins, m),
    counts). Empty bins give NaN."""
    r = np.asarray(r, dtype=float)
    v = np.asarray(values, dtype=float)
    v = v[:, None] if v.ndim == 1 else v
    w = np.asarray(weights, dtype=float)
    edges = np.asarray(edges, dtype=float)
    idx = np.searchsorted(edges, r, side="right") - 1
    nb = len(edges) - 1
    ok = (idx >= 0) & (idx < nb)
    sw = np.bincount(idx[ok], weights=w[ok], minlength=nb)
    means = np.full((nb, v.shape[1]), math.nan)
    for j in range(v.shape[1]):
        s = np.bincount(idx[ok], weights=w[ok]*v[ok, j], minlength=nb)
        with np.errstate(invalid="ignore", divide="ignore"):
            means[:, j] = np.where(sw > 0.0, s/np.where(sw > 0.0, sw, 1.0), math.nan)
    counts = np.bincount(idx[ok], minlength=nb)
    return 0.5*(edges[:-1] + edges[1:]), means, counts
