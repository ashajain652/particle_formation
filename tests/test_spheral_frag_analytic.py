"""spheral_frag.analytic: Timoshenko and Goodier's heated solid sphere (Spheral M2 plan, Task 6 Step 1)."""
import math

import numpy as np
import pytest

from spheral_frag import analytic

E, NU, ALPHA, R = 71.7e9, 0.33, 23.4e-6, 0.05


def test_uniform_temperature_is_stress_free():
    r = np.linspace(0.0, R, 11)
    sr, st = analytic.sphere_thermal_stress(r, lambda x: np.full_like(x, 350.0), R, E, NU, ALPHA, 300.0)
    assert np.abs(sr).max() < 1e-6 and np.abs(st).max() < 1e-6          # Pa; terms of order 1e8 cancel


@pytest.mark.parametrize("profile", [lambda x: 300.0 + 50.0*(x/R)**2, lambda x: 300.0 + 80.0*np.sin(30.0*x)**2,
                                     lambda x: 400.0 - 100.0*np.exp(-(x/0.02)**2)])
def test_free_surface_and_centre_limit(profile):
    sr, st = analytic.sphere_thermal_stress(np.array([0.0, 1e-9, R]), profile, R, E, NU, ALPHA, 300.0)
    scale = ALPHA*E*100.0/(1.0 - NU)
    assert abs(sr[2]) <= 1e-12*scale                                    # sigma_r(R) = 0
    assert abs(sr[0] - st[0]) <= 1e-12*scale                            # sigma_r = sigma_theta at the centre
    assert abs(sr[1] - sr[0]) <= 1e-9*scale                             # and continuous there


def test_quadratic_profile_against_its_closed_form():
    dT = 50.0
    r = np.linspace(0.0, R, 201)
    sr, st = analytic.sphere_thermal_stress(r, lambda x: 300.0 + dT*(x/R)**2, R, E, NU, ALPHA, 300.0)
    er, et = analytic.sphere_quadratic_stress(r, R, E, NU, ALPHA, dT)
    s0 = analytic.quadratic_peak(E, NU, ALPHA, dT)
    assert s0 == pytest.approx(50.08298507e6, rel=1e-9)                # 2 alpha E dT / (5 (1 - nu))
    assert np.abs(sr - er).max() <= 1e-10*s0 and np.abs(st - et).max() <= 1e-10*s0
    assert er[0] == pytest.approx(s0) and et[-1] == pytest.approx(-s0) and er[-1] == 0.0


def test_non_polynomial_profile_converges():
    """A profile no Gauss rule integrates exactly: the default rule agrees with a much finer one to 1e-10."""
    f = lambda x: 300.0 + 80.0*np.sin(30.0*x)**2 + 20.0*np.exp(x/R)
    r = np.linspace(0.0, R, 51)
    a = analytic.sphere_thermal_stress(r, f, R, E, NU, ALPHA, 300.0)
    b = analytic.sphere_thermal_stress(r, f, R, E, NU, ALPHA, 300.0, points=64, panels=32)
    scale = ALPHA*E*100.0/(1.0 - NU)
    assert max(np.abs(a[0] - b[0]).max(), np.abs(a[1] - b[1]).max()) <= 1e-10*scale


def test_spherical_components_3d_and_rz():
    rng = np.random.default_rng(3)
    x = rng.normal(size=(50, 3))
    n = x/np.linalg.norm(x, axis=1)[:, None]
    srr, stt = rng.normal(size=50), rng.normal(size=50)
    sig = stt[:, None, None]*np.eye(3)[None] + (srr - stt)[:, None, None]*np.einsum("ni,nj->nij", n, n)
    r, a, b = analytic.spherical_components(x, sig)
    assert np.allclose(r, np.linalg.norm(x, axis=1)) and np.allclose(a, srr) and np.allclose(b, stt)
    # RZ: (z, r) positions, tensor rows (z, r, theta) of the same field at azimuth 0
    xz = np.stack([x[:, 0], np.hypot(x[:, 1], x[:, 2])], axis=1)
    m = np.zeros((50, 3))
    m[:, :2] = xz/np.linalg.norm(xz, axis=1)[:, None]
    sig2 = stt[:, None, None]*np.eye(3)[None] + (srr - stt)[:, None, None]*np.einsum("ni,nj->nij", m, m)
    r2, a2, b2 = analytic.spherical_components(xz, sig2)
    assert np.allclose(r2, r) and np.allclose(a2, srr) and np.allclose(b2, stt)


def test_expansion_alpha_and_binning():
    a = analytic.expansion_strain_alpha(ALPHA, 300.0, 293.0)
    L = lambda T: 1.0 + ALPHA*(T - 293.0)
    assert L(350.0)/L(300.0) - 1.0 == pytest.approx(a*50.0, rel=1e-14)
    c, m, cnt = analytic.binned_profile([0.1, 0.2, 0.6], [1.0, 3.0, 5.0], [1.0, 1.0, 2.0], [0.0, 0.5, 1.0, 1.5])
    assert np.allclose(c, [0.25, 0.75, 1.25]) and m[0, 0] == 2.0 and m[1, 0] == 5.0 and math.isnan(m[2, 0])
    assert list(cnt) == [2, 1, 0]
