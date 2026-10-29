"""Sub-plan 21 Task 8: the share of the columns' flow at short wavelengths, for the depth-uniform pressure's limitation
(layered-skin design section 14, assumption 1; Asha, 2026-10-21, option a).

The columns give every layer of a column the same push along the surface (spec 6.4): the pressure varies with depth only
through the deceleration's component into the surface. That is exact for a disturbance whose wavelength is long against
the flowing depth h, kh small (k = 2 pi / lambda). Where kh is not small a disturbance really stirs only the top 1/k or
so of the layer, and the columns' response to it is too fast: its rate -- the growth of a Rayleigh-Taylor mode, the
decay of a capillary ripple, the speed of a wave -- is overstated by

    inertial layer (liquid; its viscous time h^2 / nu is 2.7 s at 1.2 mm):   sqrt(kh / tanh kh)       (spec 14.1)
    viscous layer (slurry; 0.01-0.04 s at 1-3.4 mm and 1-5 Pa s):            (2/3) (kh)^3 (cosh 2kh + 2 (kh)^2 + 1)
                                                                              / (sinh 2kh - 2kh)

the second the ratio of the lubrication rate h^3 k^2 (rho g + sigma k^2) / (3 mu), which is what columns with a viscous
profile and a depth-uniform pressure give, to the Stokes-flow rate of a viscous layer on a rigid base,
(rho g + sigma k^2) (sinh 2kh - 2kh) / (2 mu k (cosh 2kh + 2 (kh)^2 + 1)). At kh = 0.5 they are 1.04 and 1.45, at kh = 1
1.15 and 2.77. Both are 1 for a long wave: the runoff along the meridian, which varies over tens of millimetres, is
where they are a few per cent at most.

The measure (ruling, 2026-10-24): the share of the flowing layers' volume flux q -- each column's sum of its cells'
thicknesses times their velocities [m2/s], what crosses a unit width of the surface -- carried at wavelengths where kh
exceeds KH_SHORT, over the columns whose flowing depth is at least THICK (1 mm). The flux is the measure because it is
what the columns compute wrongly where kh is not small, and what moves the metal: where the runoff takes it and how deep
the flowing layer is that Girin's spraying reads (spec 6.7, 8). A short-wave share of the flux s means that at most the
fraction s of where the metal goes is subject to the overstated rates; the flux-weighted mean excess of the factors
above over the modes, sum P_n (F(k_n h_n) - 1) / sum P_n, says by how much on average. The top surface's height is
measured the same way (its variance at short wavelengths: the shape of the layer, which the rates move). KH_SHORT is
where the inertial factor reaches 1.05 (kh = 0.56); the viscous factor does at kh = 0.17.

Wavelengths on the columns are the modes of their surface Laplacian over the region (`modes`): L phi = k^2 M phi, L the
cotangent Laplacian on the region's edges (each edge's clipped weight, `Columns.conduct_weight`, with free ends at the
region's border), M the columns' areas, so that each mode is a wave of wavenumber k on the surface -- on a flat grid of
squares a sampled cosine of wavenumber kappa is one exactly, with k = (2 / dx) sin(kappa dx / 2). A field f (per column,
or per column and Cartesian component) has the power P_n = sum over components of (phi_n^T M f)^2 in mode n, and each
mode its own depth, h_n = sum_j A_j phi_nj^2 h_j (the modes are M-orthonormal)."""
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

KH_SHORT = 0.56            # kh above which the columns overstate an inertial layer's rates by more than 5 %
THICK = 1.0e-3             # m: the flowing layers measured (spec 14.1: "in layers over about 1 mm")


def inertial_factor(kh):
    """sqrt(kh / tanh kh): how much the columns overstate an inertial layer's rates at kh (1 at kh = 0)."""
    kh = np.asarray(kh, dtype=float)
    safe = np.where(kh > 1e-8, kh, 1.0)
    return np.where(kh > 1e-8, np.sqrt(safe / np.tanh(safe)), 1.0)


def viscous_factor(kh):
    """The lubrication rate over the Stokes rate of a viscous layer on a rigid base at kh (module docstring; 1 at kh = 0):
    how much the columns overstate a viscous layer's rates. Its series 1 + (9/5) (kh)^2 below kh = 0.05, where the
    difference sinh 2kh - 2kh loses its digits."""
    kh = np.asarray(kh, dtype=float)
    x = np.where(kh > 0.05, kh, 1.0)
    exact = (2.0 / 3.0) * x ** 3 * (np.cosh(2 * x) + 2 * x ** 2 + 1.0) / (np.sinh(2 * x) - 2 * x)
    return np.where(kh > 0.05, exact, 1.0 + 1.8 * kh ** 2)


def modes(columns, region=None):
    """The surface Laplacian's modes over the columns `region` ((n,) bool, all if None): (k [1/m, ascending], phi
    (m, m) M-orthonormal by columns, idx) for the m columns idx of the region, from L phi = k^2 M phi -- L the cotangent
    Laplacian on the edges with both ends in the region, its clipped weights, M the columns' areas (module docstring).
    Dense: about 30 s for 4 000 columns."""
    region = np.ones(columns.n, dtype=bool) if region is None else np.asarray(region, dtype=bool)
    idx = np.flatnonzero(region)
    local = np.full(columns.n, -1)
    local[idx] = np.arange(len(idx))
    i, j = columns.edges[:, 0], columns.edges[:, 1]
    keep = region[i] & region[j]
    a, b, w = local[i[keep]], local[j[keep]], columns.conduct_weight[keep]
    m = len(idx)
    L = sp.csr_matrix((np.r_[w, w, -w, -w], (np.r_[a, b, a, b], np.r_[a, b, b, a])), shape=(m, m)).toarray()
    lam, phi = sla.eigh(L, np.diag(columns.area[idx]))
    return np.sqrt(np.clip(lam, 0.0, None)), phi, idx


def wavelength_share(columns, field, depth, region=None, kh_short=KH_SHORT, basis=None):
    """The share of `field` ((n,) or (n, 3), per column) carried at kh above `kh_short` over the columns `region` (all if
    None), and the field-weighted mean excess of the columns' rates over the exact ones (module docstring), each mode at
    its own depth from `depth` [m, (n,)]. `basis`: `modes(columns, region)` made already. Returns a dict: share,
    inertial and viscous (the mean excesses, F - 1), and per mode k [1/m], h [m], kh and power."""
    k, phi, idx = modes(columns, region) if basis is None else basis
    f = np.asarray(field, dtype=float)[idx]
    f = f[:, None] if f.ndim == 1 else f
    A = columns.area[idx]
    power = ((phi.T @ (A[:, None] * f)) ** 2).sum(axis=1)
    h = (A[:, None] * phi ** 2 * np.asarray(depth, dtype=float)[idx][:, None]).sum(axis=0)
    kh = k * h
    total = power.sum()
    if total <= 0.0:
        return {"share": 0.0, "inertial": 0.0, "viscous": 0.0, "k": k, "h": h, "kh": kh, "power": power}
    return {"share": float(power[kh > kh_short].sum() / total),
            "inertial": float((power * (inertial_factor(kh) - 1.0)).sum() / total),
            "viscous": float((power * (viscous_factor(kh) - 1.0)).sum() / total),
            "k": k, "h": h, "kh": kh, "power": power}


def column_flux(sk):
    """Each column's volume flux [m2/s, (n, 3)]: its cells' thicknesses times their velocities, summed."""
    return (sk.dz()[..., None] * sk.U).sum(axis=1)


def kh_at(factor, target, lo=1e-6, hi=10.0):
    """The kh at which `factor` (inertial_factor or viscous_factor, increasing) reaches `target`, by bisection."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if float(factor(mid)) < target else (lo, mid)
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    print("inertial factor 1.05 at kh = %.3f; viscous 1.05 at kh = %.3f" % (kh_at(inertial_factor, 1.05),
                                                                         kh_at(viscous_factor, 1.05)))
    for kh in (0.1, 0.25, 0.5, 1.0, 2.0):
        print("kh %.2f: inertial %.4f, viscous %.4f" % (kh, inertial_factor(kh), viscous_factor(kh)))
