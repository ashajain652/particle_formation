"""surface_flow.py: Ranger's boundary layer, the Thwaites cross-check, the edge state and regimes at 71 km."""
import math
from datetime import datetime

import numpy as np
import pytest

from reentry_model import aero, atmosphere, body, surface_flow as sf
from reentry_model import trajectory as tj

pytest.importorskip("cantera")
EPOCH = datetime(2024, 8, 1, 12, 53, 7)
R100 = tj.InitialState(7500.0, 77500.133, math.radians(-0.959331), math.radians(347.168296), math.radians(29.546067),
                       math.radians(-82.134333), EPOCH)


def potential_flow(R=0.05, V=1000.0, nu=1e-4):
    th = sf.THETA_BINS
    return th, R * th, 1.5 * V * np.sin(th), np.full_like(th, nu), V * 2.0 * R / nu


def test_ranger_form_reproduces_the_closed_form():
    th, s, u, nu, Re = potential_flow()
    d = sf.boundary_layer_thickness(s, u, nu)
    ref = sf.ranger_thickness(0.05, Re, th[1:])
    assert np.abs(d[1:] / ref - 1.0).max() < 3e-3 and np.abs(d[2:] / ref[1:] - 1.0).max() < 1e-3
    assert d[0] == d[1] and sf.ranger_psi(1e-3) == pytest.approx(math.sqrt(3.2), rel=1e-4)      # Psi(0) = sqrt(48/15)


def test_thwaites_shape_agrees_with_ranger_within_five_percent():
    th, s, u, nu, Re = potential_flow()
    ranger = sf.boundary_layer_thickness(s, u, nu)
    thwaites = sf.boundary_layer_thickness(s, u, nu, sf.THWAITES_C, 5)
    m = (th >= math.radians(5.0)) & (th <= math.radians(85.0))
    ratio = ranger[m] / thwaites[m]
    assert 11.5 < ratio.min() and ratio.max() < 13.0 and (ratio.max() - ratio.min()) / ratio.mean() < 0.05    # measured 12.0-12.4


def aero_state_at(t):
    sim = tj.Simulator(R100, body.ConstantBody(1.473), atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(),
                       aero.SesamTable(), tj.Settings(diameter=0.1))
    sim.advance(t)
    return sim.aero_state(sim.t, sim.y[:3], sim.y[3:])


def test_prandtl_meyer_against_the_textbook():
    """nu(M) and its inverse (Anderson, Modern Compressible Flow, Table A.5: nu(2) = 26.38 deg at gamma 1.4)."""
    assert math.degrees(sf.prandtl_meyer(1.0, 1.4)) == pytest.approx(0.0, abs=1e-9)
    assert math.degrees(sf.prandtl_meyer(2.0, 1.4)) == pytest.approx(26.3798, abs=1e-3)
    assert math.degrees(sf.prandtl_meyer(5.0, 1.4)) == pytest.approx(76.9202, abs=1e-3)
    assert sf.mach_from_turn(math.radians(26.3798), 1.4) == pytest.approx(2.0, rel=1e-5)
    assert sf.mach_from_turn(0.0, 1.4) == pytest.approx(1.0)


def test_sonic_point_and_wall_pressure():
    """The Newtonian/Prandtl-Meyer switch: phi* = 43.38 deg at gamma 1.4 and 40.72 deg at 1.15, p_w continuous there,
    monotone, floored at p_inf, and 20-50 x the Newtonian value at 90 degrees (spec section 18)."""
    p_inf, q_inf, p_stag = 5.221, 2141.0, 4165.3                      # 100 mm sphere at 70.0 km, V = 7190 m/s
    for gamma, expect in ((1.4, 43.38), (1.15, 40.72)):
        assert math.degrees(sf.sonic_angle(p_stag, p_inf, q_inf, gamma)) == pytest.approx(expect, abs=0.05)
    theta = np.radians(np.linspace(0.0, 90.0, 181))
    p_w, phi_star = sf.wall_pressure(theta, p_stag, p_inf, q_inf, 1.15)
    assert p_w[0] == pytest.approx(p_stag, rel=1e-6) and np.all(np.diff(p_w) < 0.0) and np.all(p_w >= p_inf)
    newtonian = p_inf + q_inf * (p_stag - p_inf) / q_inf * np.cos(theta) ** 2
    below = theta < phi_star - math.radians(sf.PM_BLEND_DEG)
    assert np.allclose(p_w[below], newtonian[below], rtol=1e-9)        # unchanged below the sonic point
    assert p_w[-1] == pytest.approx(248.8, rel=0.02) and p_w[-1] / newtonian[-1] > 40.0   # 90 deg: 249 Pa vs 5.2 Pa Newtonian (measured C_p: 112-219)
    fine = np.radians(np.linspace(0.0, 90.0, 1801))
    smooth = np.abs(np.diff(sf.wall_pressure(fine, p_stag, p_inf, q_inf, 1.15)[0], 2))
    assert smooth.max() < 15.0 * np.median(smooth[smooth > 0])         # the blend starts at phi*, so no kink (35x if it straddles)
    assert sf.wall_pressure(theta, p_stag, p_inf, q_inf, 1.4)[0][-1] == pytest.approx(144.3, rel=0.02)     # gamma is a factor ~2 on p_w(90)


def test_wall_knudsen_is_the_maxwell_path_over_the_nose_radius():
    """Kn_local = lambda_w/R with the wall gas ideal at (p_w, T_w): lambda_w = mu sqrt(pi R_s T/2)/p_w, and the
    identity Kn = (Ma/Re) sqrt(gamma pi/2) (coefficient 1.4829 at gamma 1.4)."""
    mu_w, T_w, p_w, R = 5.5e-5, 908.0, np.array([3667.0, 231.0]), 0.05
    lam, kn = sf.wall_knudsen(p_w, np.full(2, T_w), mu_w, R)
    assert np.allclose(lam, mu_w * math.sqrt(math.pi * sf.R_SPECIFIC_AIR * T_w / 2.0) / p_w)
    assert np.allclose(kn, lam / R) and kn[0] < 1e-3 < kn[1]           # the rim is the rarefied place, not the nose
    rho, T, V, L, mu = 8.28e-5, 219.6, 7250.0, 0.05, 1.7e-5
    a = math.sqrt(sf.GAMMA_AIR * sf.R_SPECIFIC_AIR * T)
    identity = (V / a) / (rho * V * L / mu) * math.sqrt(sf.GAMMA_AIR * math.pi / 2.0)
    assert sf.mean_free_path_maxwell(mu, rho, T) / L == pytest.approx(identity, rel=1e-12)
    assert math.sqrt(sf.GAMMA_AIR * math.pi / 2.0) == pytest.approx(1.4829, abs=1e-4)


def test_branch_gate_and_closures_at_71_km():
    """The body-scale gate and the wall-scale closure gate (spec section 18). At 71 km the 100 mm sphere has
    Kn_body = 0.0113 -- just above the 0.01 threshold -- so the branch is merged and every closure is flagged, while
    the wall Knudsen number is 1e-4 to 2e-3 (the declared conservatism: lambda_w/lambda_inf ~ 1/170)."""
    a = aero_state_at(43.5)
    assert 70.5e3 < a.h < 71.5e3
    theta = np.radians(np.array([0.5, 10.0, 30.0, 45.0, 60.0, 75.0, 89.0, 100.0, 150.0]))
    T_wall = np.full(theta.size, 908.0)
    flow = sf.SurfaceFlow(gamma_pm=1.15)
    r = flow.evaluate(a, theta, 0.05, 2400.0, T_wall)
    assert r.branch == sf.BRANCH_MERGED and r.kn_body == pytest.approx(0.0113, rel=0.05) and 140.0 < r.re_shock < 210.0
    assert np.all(r.closure == sf.CLOSURE_COUETTE) and np.all(r.flagged)          # merged: nothing is certified
    assert 1e-5 < r.kn_local[0] < 1e-3 and r.kn_local[0] < r.kn_local[6] < 0.01   # compressed cold wall; the rim is the loosest
    assert r.kn_local[0] == pytest.approx(r.kn_body / 84.0, rel=0.2)               # lambda_w/lambda_inf ~ 1/170, halved again by D/R
    assert 1000.0 < r.p_stag < 6000.0 and r.p_w[0] == pytest.approx(r.p_stag, rel=1e-3) and r.p_w[7] == r.p_inf
    assert 150.0 < r.p_w[6] < 350.0                                               # 89 deg: Prandtl-Meyer, not p_inf
    assert 1100.0 < r.u_e[2] < 1500.0 and 1900.0 < r.u_e[3] < 2300.0
    assert 4e-3 < r.delta_a[2] < 9e-3 and np.isnan(r.delta_a[7])
    assert r.tau[3] > 0.0 and r.tau_fm[3] > r.tau[3] and np.all(r.tau[7:] == 0.0)
    assert r.G[3] > 0.0 and np.all(r.G[7:] == 0.0)
    assert np.allclose(r.closure_fractions(np.ones(theta.size), theta <= np.pi / 2), [0.0, 1.0, 0.0])
    with pytest.raises(ValueError):
        sf.SurfaceFlow(rarefied_shear="magic")
    with pytest.raises(ValueError):
        sf.SurfaceFlow(gamma_pm=0.9)


def test_shock_layer_branch_certifies_the_nose_lower_down():
    """Below the gate (60 km, Kn_body 0.0026) the branch is the shock layer and the windward face takes Girin's
    closure; the free-molecular branch has no edge state at all."""
    a = aero_state_at(43.5)
    a.kn = 0.005                                                        # the same state, below the body gate
    theta = np.radians(np.array([0.5, 45.0, 89.0, 120.0]))
    r = sf.SurfaceFlow().evaluate(a, theta, 0.05, 2400.0, np.full(4, 908.0))
    assert r.branch == sf.BRANCH_SHOCK_LAYER and np.all(r.closure[:3] == sf.CLOSURE_GIRIN) and not r.flagged[:3].any()
    assert np.isfinite(r.re_shock) and np.isfinite(r.phi_sonic)
    a.kn = 20.0
    r = sf.SurfaceFlow().evaluate(a, theta, 0.05, 2400.0, np.full(4, 908.0))
    assert r.branch == sf.BRANCH_FREE_MOLECULAR and np.all(np.isnan(r.delta_a)) and np.all(r.u_e == 0.0)
    assert np.isnan(r.p_stag) and np.all(r.p_w == r.p_inf) and np.all(r.closure == sf.CLOSURE_COUETTE)
    assert r.tau[1] == pytest.approx(a.freestream.rho * a.V ** 2 * math.sin(theta[1]) * math.cos(theta[1])) and r.tau[3] == 0.0
