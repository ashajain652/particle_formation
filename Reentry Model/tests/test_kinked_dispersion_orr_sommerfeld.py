"""analysis/kinked_dispersion/orr_sommerfeld.py, the viscous (Orr-Sommerfeld) check of the kinked-profile dispersion
(spec stage 6, Task 10), validated against:

- Lamb's exact relation for a damped capillary wave on deep, quiescent liquid (the free surface's conditions);
- the sign of the long-wave two-layer Couette mode (Charru & Hinch 2000, "stable if the thinner layer is the less
  viscous, and unstable otherwise"), with Yih's (1967) leading-order wave speed (Charru & Hinch, Appendix C);
- Hooper & Boyd's (1983, Eq. 15) short-wave growth of two deep fluids in uniform shear;
- Girin's profile at a Reynolds number V_s delta / nu of 1e4 and 1e5 against the inviscid broken-line growth
  (`reentry_model.shared.dispersion`);
- and resolution doubling, including a flight column of the d050 fixture.

Units: the validations are dimensionless (lengths, speeds and densities of order 1), which the solver accepts as any
consistent set."""
import os

import numpy as np
import pytest

from analysis.kinked_dispersion import orr_sommerfeld as osm
from reentry_model.shared import dispersion

FIXTURE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "reference_values",
                       "d050_layer_columns.npz")


def layered(*layers):
    """Faces (depths, downward), centre speeds, viscosities, densities and face speeds of a column of uniform layers,
    each `(thickness, u_top, u_bottom, mu, rho, cells, grade)`: `cells` cells per layer, the speed linear across it, the
    cells' thicknesses graded geometrically by `grade` (> 1 thins them toward the layer's bottom, < 1 toward its top)."""
    faces, uf, U, mu, rho = [0.0], [layers[0][1]], [], [], []
    for h, u0, u1, m, r, n, grade in layers:
        w = grade ** np.arange(n)[::-1] if grade != 1.0 else np.ones(n)
        edges = np.concatenate([[0.0], np.cumsum(w) / w.sum()]) * h
        top = faces[-1]
        for a, b in zip(edges[:-1], edges[1:]):
            faces.append(top + b)
            uf.append(u0 + (u1 - u0) * b / h)
            U.append(u0 + (u1 - u0) * 0.5 * (a + b) / h)
            mu.append(m)
            rho.append(r)
    return np.array(faces), np.array(U), np.array(mu), np.array(rho), np.array(uf)


# ---- the free surface: a capillary wave on deep, quiescent liquid ---------------------------------------------------

def lamb_root(k, nu, sigma, rho, sign):
    """The complex rate s (perturbations ~ exp(s t)) of a capillary wave of wavenumber k on deep liquid (Lamb 1932,
    section 349): (s + 2 nu k^2)^2 + omega_0^2 = 4 nu^2 k^3 (k^2 + s / nu)^(1/2), omega_0^2 = sigma k^3 / rho, solved by
    Newton's method from the weak-damping root s = -2 nu k^2 + sign i omega_0."""
    w0 = np.sqrt(sigma * k ** 3 / rho)
    s = complex(-2.0 * nu * k ** 2, sign * w0)
    for _ in range(60):
        q = np.sqrt(k ** 2 + s / nu)
        f = (s + 2.0 * nu * k ** 2) ** 2 + w0 ** 2 - 4.0 * nu ** 2 * k ** 3 * q
        df = 2.0 * (s + 2.0 * nu * k ** 2) - 2.0 * nu * k ** 3 / q
        s -= f / df
    return s


def test_capillary_wave_damping_matches_lambs_relation():
    k, nu, sigma, rho = 1.0, 0.1, 1.0, 1.0
    y, U, mu, r, uf = layered((20.0, 0.0, 0.0, nu * rho, rho, 8, 1.6 ** -1))
    spec = osm.spectrum(k, y, U, mu, r, sigma, u_faces=uf)
    for sign in (1.0, -1.0):
        c_exact = 1j * lamb_root(k, nu, sigma, rho, sign) / k        # exp(s t) = exp(-i k c t)
        c = spec.c[np.argmin(np.abs(spec.c - c_exact))]
        assert abs(c - c_exact) < 1e-6 * abs(c_exact)


def test_quiescent_and_couette_columns_are_stable():
    """Plane Couette flow of one fluid is stable at every Reynolds number, and so is a liquid at rest: the solver must
    pick no spurious unstable eigenvalue, at either resolution."""
    for res in (1, 2):
        y, U, mu, r, uf = layered((1.0, 1.0, 0.0, 1e-3, 1.0, 4, 1.0))
        assert osm.growth_rate(1.0, y, U, mu, r, 0.0, top="wall", u_faces=uf, resolution=res).imag < 0.0
        y, U, mu, r, uf = layered((1.0, 0.0, 0.0, 1e-2, 1.0, 4, 1.0))
        assert osm.growth_rate(3.0, y, U, mu, r, 1.0, u_faces=uf, resolution=res).imag < 0.0


def test_shift_invert_and_arpack_agree_with_qz():
    """The eigenvalue taken is the same by the dense shift and invert, by ARPACK and by scipy's QZ on the whole pencil:
    on a Charru & Hinch column (the smallest growth of the set, 1e-6 of the speed scale) and on Girin's profile."""
    (y, U, mu, r, uf), _ = couette(0.8, 0.2, 0.5)
    k = 0.05 / 0.8
    cq = osm.growth_rate(k, y, U, mu, r, 0.0, top="wall", u_faces=uf, method="qz")
    for method in ("dense", "arpack"):
        c = osm.growth_rate(k, y, U, mu, r, 0.0, top="wall", u_faces=uf, method=method)
        assert abs(c - cq) < 1e-9 * abs(U).max()
    (y, U, mu, r, uf), sigma = girin_column(1e3, depth=6.0)
    k = 1.2
    cd = osm.growth_rate(k, y, U, mu, r, sigma, u_faces=uf, method="dense")
    ca = osm.growth_rate(k, y, U, mu, r, sigma, u_faces=uf, method="arpack")
    assert abs(ca - cd) < 1e-9 and cd.imag > 0.0


# ---- the default base profile (no face speeds given) ----------------------------------------------------------------

def test_default_face_speeds_recover_a_two_layer_couette_profile():
    """Stress continuity at the faces (the skins' own rule) recovers a two-layer Couette profile from its cells' centre
    speeds, and the surface speed is extrapolated from the top two centres (the straight-segment fit's rule)."""
    y, U, mu, r, uf = layered((0.3, 1.0, 0.4, 1.0, 1.0, 3, 1.0), (0.7, 0.4, 0.0, 2.0 / (0.4 / 0.7), 1.0, 4, 1.3))
    prof = osm.base_profile(y, U, mu, r)
    assert np.allclose(prof.u_faces, uf, rtol=0.0, atol=1e-12)
    assert not prof.kinked.any()                          # a straight speed in every cell: no kink at any centre


# ---- Charru & Hinch: the sign of the long-wave mode in two-layer Couette flow ----------------------------------------

def couette(h1, h2, m, re=1.0, a1=1.0, rho=1.0):
    """Two-layer Couette flow between walls: a lower layer h1 of viscosity mu1 sheared at a1, an upper one h2 of
    viscosity m mu1 sheared at a1 / m (stress continuity), the bottom wall at rest; mu1 from the lower layer's shear
    Reynolds number a1 h1^2 / nu1 = re."""
    mu1 = rho * a1 * h1 ** 2 / re
    u_i = a1 * h1
    u_w = u_i + a1 / m * h2
    return layered((h2, u_w, u_i, m * mu1, rho, 6, 1.0), (h1, u_i, 0.0, mu1, rho, 6, 1.0)), u_i


def yih_wave_speed(m, d):
    """Charru & Hinch (2000), Appendix C: the long-wave interfacial wave speed over a1 h1 relative to the interface's
    speed, at leading order, for m = mu2 / mu1 and d = h2 / h1 (layer 1 the lower)."""
    D = d ** 4 + 4.0 * d * m + 6.0 * d ** 2 * m + 4.0 * d ** 3 * m + m ** 2
    return 2.0 * (1.0 - m) * d ** 2 * (1.0 + d) / D


@pytest.mark.parametrize("h1, h2, m, unstable", [
    (0.2, 0.8, 0.5, True),     # thin lower layer, the more viscous: unstable
    (0.2, 0.8, 2.0, False),    # thin lower layer, the less viscous: stable
    (0.8, 0.2, 2.0, True),     # thin upper layer, the more viscous: unstable
    (0.8, 0.2, 0.5, False),    # thin upper layer, the less viscous: stable
])
def test_long_wave_couette_sign_follows_charru_hinch(h1, h2, m, unstable):
    (y, U, mu, r, uf), u_i = couette(h1, h2, m)
    k = 0.05 / h1                                          # alpha_1 = k h1 = 0.05, a long wave
    c = osm.growth_rate(k, y, U, mu, r, 0.0, top="wall", u_faces=uf)
    assert (c.imag > 0.0) == unstable
    assert (c.real - u_i) / h1 == pytest.approx(yih_wave_speed(m, h2 / h1), rel=0.05)


# ---- Hooper & Boyd: the short-wave interfacial mode of two deep fluids in uniform shear --------------------------------

def hooper_boyd_growth(alpha, m, a2=1.0):
    """Hooper & Boyd (1983), Eq. 15: the growth rate a2 Im c = a2 m (1 - m)^2 / [2 (1 + m)] alpha^-2 + O(alpha^-6) of
    the short-wave mode (Im c has no alpha^-4 term: their c_2 is real), with alpha = k (nu2 / a2)^(1/2), fluid 2 the
    lower and m = mu2 / mu1 (their Eq. 5c, tangential stress, and 'suppose m < 1, i.e. the lower fluid is the less
    viscous', section 3.3)."""
    return a2 * m * (1.0 - m) ** 2 / (2.0 * (1.0 + m)) / alpha ** 2


@pytest.mark.parametrize("m", [0.5, 2.0])
@pytest.mark.parametrize("alpha", [5.0, 10.0])
def test_short_wave_growth_matches_hooper_boyd(alpha, m):
    a2, nu2, rho = 1.0, 1.0, 1.0                           # the lower fluid sets the units: (nu2 / a2)^(1/2) = 1
    mu2 = nu2 * rho
    mu1, a1 = mu2 / m, m * a2                              # m = mu2 / mu1, and mu1 a1 = mu2 a2
    L = 12.0
    y, U, mu, r, uf = layered((L, a1 * L, 0.0, mu1, rho, 6, 1.0 / 1.5), (L, 0.0, -a2 * L, mu2, rho, 6, 1.5))
    k = alpha
    c = osm.growth_rate(k, y, U, mu, r, 0.0, top="wall", u_faces=uf, method="dense")
    assert k * c.imag == pytest.approx(hooper_boyd_growth(alpha, m, a2), rel=0.05)


# ---- Girin's profile: the viscous growth tends to the inviscid broken-line growth -------------------------------------

def girin_column(re, we=1e4, depth=12.0):
    """Girin's profile in his units (layer depth, surface speed and density 1): speed falling linearly from 1 at the
    surface to 0 at depth 1, at rest below to a wall at `depth`; viscosity 1 / re, surface tension 1 / we."""
    return layered((1.0, 1.0, 0.0, 1.0 / re, 1.0, 8, 1.0), (depth - 1.0, 0.0, 0.0, 1.0 / re, 1.0, 6, 1.0 / 1.6)), 1.0 / we


def broken_line_fastest(we=1e4, depth=12.0):
    prof = dispersion.KinkedProfiles.one([1.0], [1.0, 0.0, 0.0], wall=depth, rho=1.0)
    w = dispersion.fastest_waves(prof, 1.0 / we)
    return float(w.k[0]), float(w.growth[0])


def girin_errors():
    k_bl, g_bl = broken_line_fastest()
    errors = {}
    for re in (1e4, 1e5):
        (y, U, mu, r, uf), sigma = girin_column(re)
        f = osm.fastest(y, U, mu, r, sigma, (0.5 * k_bl, 2.0 * k_bl), u_faces=uf, method="dense")
        assert not f.edge
        errors[re] = (f.growth / g_bl - 1.0, f.k / k_bl - 1.0)
    return errors


@pytest.mark.reference
def test_girin_profile_viscous_growth_tends_to_the_broken_line():
    errors = girin_errors()
    for re, (eg, ek) in errors.items():
        assert abs(eg) < 0.05, (re, eg)
    assert abs(errors[1e5][0]) < abs(errors[1e4][0])


def test_girin_profile_at_re_1e4_at_the_broken_line_wavenumber():
    """The fast tier's Girin check: at the broken line's fastest wavenumber, the viscous growth at Re = 1e4 is within
    5 % of the inviscid one, and doubling the resolution moves it by less than 1e-4 of itself."""
    k_bl, g_bl = broken_line_fastest()
    (y, U, mu, r, uf), sigma = girin_column(1e4)
    g1 = k_bl * osm.growth_rate(k_bl, y, U, mu, r, sigma, u_faces=uf, method="dense").imag
    g2 = k_bl * osm.growth_rate(k_bl, y, U, mu, r, sigma, u_faces=uf, resolution=2, method="dense").imag
    assert g1 == pytest.approx(g_bl, rel=0.05)
    assert g2 == pytest.approx(g1, rel=1e-4)


# ---- the sweep's synthetic columns -----------------------------------------------------------------------------------

def test_synthetic_column_is_stress_continuous_and_its_floor_reynolds_number_is_the_models():
    """floor_sweep's two-layer Couette column: the shear stress is continuous at the kink, its slope ratio is 1 / M,
    and the Re_f the sweep records, Re_1 / (1.225 M), is what the profile path's floor test computes on it."""
    from analysis.kinked_dispersion import floor_sweep as fs
    re1, M, b = 1e3, 100.0, 15.0
    y, U, mu, r, uf = fs.synthetic_column(re1, M, b)
    prof = osm.base_profile(y, U, mu, r, u_faces=uf)
    assert not prof.kinked.any()
    s = -np.diff(uf) / -np.diff(y)                        # dU/dy, y up, per cell
    top = y[:-1] < 1.0
    assert np.allclose(mu * s, mu[0] * s[0], rtol=1e-12)            # one shear stress through the column
    assert s[~top][0] / s[top][0] == pytest.approx(1.0 / M)
    kp = dispersion.KinkedProfiles.one([1.0], [uf[0], uf[np.flatnonzero(y == 1.0)[0]], 0.0], wall=b, rho=1.0)
    re_f = dispersion.floor_reynolds(kp, np.array([[M / re1]]))[0, 0]
    assert re_f == pytest.approx(re1 / (dispersion.FASTEST_KD * M), rel=1e-12)


# ---- a flight column ------------------------------------------------------------------------------------------------

def test_flight_column_is_resolved_by_doubling():
    """A thick column of the d050 fixture, at the wavenumber 1.225 over its first kink's depth (or its top cell's when it
    has none): doubling every domain's degree moves the eigenvalue by less than 1e-3 of the column's top speed."""
    from reentry_model.shared import melt_profile
    data = np.load(FIXTURE)
    cols = melt_profile.LayerColumns(**{f: data[f] for f in ("z", "dz", "u", "T", "m", "rho", "mu", "flowing",
                                                               "liquid", "front")})
    i = int(np.flatnonzero(data["group"] == 0)[0])
    y, U, mu, r = osm.column_cells(cols, i)
    k = 1.225 / (8.0 * (y[1] - y[0]))
    c1 = osm.growth_rate(k, y, U, mu, r, 0.80)
    c2 = osm.growth_rate(k, y, U, mu, r, 0.80, resolution=2)
    assert abs(c2 - c1) < 1e-3 * abs(U[0])
