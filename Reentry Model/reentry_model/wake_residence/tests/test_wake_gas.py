"""Tests for the wake gas state (wake_gas.py, Plume.wake_state): the conservation statements the Gaussian profiles are
built to satisfy, the stagnation and speed caps, smoothness at the plume edge, the exactly-ambient cases and the four
wake_flow modes (design spec section 4.3a).
"""
import math, numpy as np, pytest
from scipy.integrate import quad
from reentry_model.wake_residence.wake_gas import gas_state, Ambient, GasState, CP_AIR, R_AIR, WAKE_FLOW_MODES
from reentry_model.wake_residence.laws import EmissionState, sutherland_mu
from reentry_model.wake_residence.overlap import Plume
from reentry_model.wake_residence.track import Track

E = EmissionState(rho=3e-4, T=240.0, p=20.0, mu=sutherland_mu(240.0), mfp=1e-4, U=4000.0, d=2.813, cda=5.5, r0=1.4065)
AMB = Ambient(p=20.0, T=240.0, rho=20.0 / (R_AIR * 240.0)); X = np.array([1.0, 0, 0])

@pytest.mark.parametrize("r_b", [40.0, 5.0])    # 40 m: the core is the energy radius; 5 m: the line-blast radius (25.7 m) is
def test_momentum_and_energy_are_conserved(r_b):    # larger than the core, so the R_E > r_b branch is integrated too
    D = 0.5 * E.rho * E.U**2 * E.cda; R_E = max(r_b, math.sqrt(D / E.p))
    st = lambda r: gas_state(r, r_b, E, X, AMB, "both")
    mom = quad(lambda r: st(r).rho * st(r).u[0] * 2 * math.pi * r, 0, 10 * r_b, limit=200)[0]
    # the enthalpy actually deposited above ambient, from the returned density and temperature (not from q)
    energy = quad(lambda r: st(r).rho * CP_AIR * (st(r).T - AMB.T) * 2 * math.pi * r, 0, 10 * R_E, limit=200)[0]
    assert mom == pytest.approx(D / E.U, rel=1e-3)           # only exact where neither cap binds (r_b chosen so)
    assert energy == pytest.approx(D, rel=1e-3)

def test_caps_hold_close_to_a_fast_parent():
    # a heavy drag on a tiny core in thin air: both the stagnation-temperature cap and the speed cap bind at the centre
    e = EmissionState(**{**E.__dict__, "U": 3000.0, "rho": 1e-2, "p": 1e4, "cda": 50.0})
    amb = Ambient(p=10.0, T=240.0, rho=10.0 / (R_AIR * 240.0))
    s = gas_state(0.0, 0.1, e, X, amb, "both")
    assert s.u[0] == pytest.approx(3000.0, rel=1e-12)
    assert s.T == pytest.approx(240.0 + 3000.0**2 / (2 * CP_AIR), rel=1e-12)

def test_density_is_continuous_across_the_cutoff_for_an_inconsistent_ambient():
    amb = Ambient(p=20.0, T=240.0, rho=1.02 * 20.0 / (R_AIR * 240.0))     # 2 percent off the ideal gas with R_AIR
    r_b = 40.0; R_E = max(r_b, math.sqrt(0.5 * E.rho * E.U**2 * E.cda / E.p)); edge = 6.0 * max(r_b, R_E)
    for mode in ("both", "density", "velocity"):
        inside, outside = gas_state(edge * (1 - 1e-12), r_b, E, X, amb, mode), gas_state(edge * (1 + 1e-12), r_b, E, X, amb, mode)
        assert inside.rho == pytest.approx(outside.rho, rel=1e-12) and outside.rho == amb.rho
        assert inside.T == pytest.approx(outside.T, rel=1e-12) and inside.p == outside.p
        assert inside.u[0] == pytest.approx(outside.u[0], abs=1e-12 * E.U)
    # in the plume the density is the ambient density times (1 - q), q = 1 - T_amb / T
    s = gas_state(5.0, r_b, E, X, amb, "both")
    assert s.rho == pytest.approx(amb.rho * amb.T / s.T, rel=1e-12) and s.rho < amb.rho

def test_state_is_smooth_across_the_plume_edge():
    r_b = 40.0; f = lambda r: gas_state(r, r_b, E, X, AMB, "both")
    for g in (lambda s: s.T, lambda s: s.u[0]):
        d1 = (g(f(r_b + 1e-3)) - g(f(r_b - 1e-3))) / 2e-3
        d2 = (g(f(r_b + 2e-3)) - g(f(r_b - 2e-3))) / 4e-3
        assert d1 == pytest.approx(d2, rel=1e-3)

def test_off_and_far_field_are_exactly_ambient():
    for s in (gas_state(1.0, 40.0, E, X, AMB, "off"), gas_state(1e5, 40.0, E, X, AMB, "both")):
        assert (s.rho, s.T, s.p) == (AMB.rho, AMB.T, AMB.p) and not s.u.any()

def test_modes_switch_one_effect_each():
    v, d = gas_state(5.0, 40.0, E, X, AMB, "velocity"), gas_state(5.0, 40.0, E, X, AMB, "density")
    assert v.T == AMB.T and v.u[0] > 0 and d.T > AMB.T and not d.u.any()

def test_velocity_mode_uses_ambient_density_and_both_uses_the_hot_density():
    v, b = gas_state(5.0, 40.0, E, X, AMB, "velocity"), gas_state(5.0, 40.0, E, X, AMB, "both")
    assert v.rho == AMB.rho and b.rho < AMB.rho and b.T > AMB.T
    assert b.p == AMB.p and b.rho == pytest.approx(AMB.rho * AMB.T / b.T, rel=1e-12)
    assert b.u[0] > v.u[0] > 0                          # same momentum density divided by the thinner gas

def test_velocity_points_along_the_direction_and_direction_is_normalised():
    d = np.array([0.0, 3.0, 4.0])
    s = gas_state(5.0, 40.0, E, d, AMB, "both")
    assert s.u[0] == 0.0 and s.u[2] / s.u[1] == pytest.approx(4.0 / 3.0)
    assert np.linalg.norm(s.u) == pytest.approx(gas_state(5.0, 40.0, E, X, AMB, "both").u[0])

def test_the_energy_radius_grows_with_k_E_and_the_cutoff_scales_with_it():
    # on the line-blast radius sqrt(D/p) = 25.7 m, a core of 1 m: k_E = 2 spreads the energy over four times the area
    a, b = gas_state(0.0, 1.0, E, X, AMB, "density", k_E=1.0), gas_state(0.0, 1.0, E, X, AMB, "density", k_E=2.0)
    q = lambda g: 1 - AMB.T / g.T                       # the excess enthalpy fraction scales with 1 / R_E^2 exactly
    assert q(b) == pytest.approx(q(a) / 4.0, rel=1e-9)
    R_E = math.sqrt(0.5 * E.rho * E.U**2 * E.cda / E.p)
    inside = gas_state(5.9 * R_E, 1.0, E, X, AMB, "density"); outside = gas_state(6.1 * R_E, 1.0, E, X, AMB, "density")
    assert inside.T > AMB.T and outside.T == AMB.T

def test_gas_state_rejects_bad_inputs():
    with pytest.raises(ValueError, match="wake_flow"):
        gas_state(1.0, 40.0, E, X, AMB, "sideways")
    with pytest.raises(ValueError, match="r_b"):
        gas_state(1.0, 0.0, E, X, AMB, "both")
    with pytest.raises(ValueError, match="r_axis"):
        gas_state(-1.0, 40.0, E, X, AMB, "both")
    with pytest.raises(ValueError, match="direction"):
        gas_state(1.0, 40.0, E, np.zeros(3), AMB, "both")
    assert set(WAKE_FLOW_MODES) == {"off", "velocity", "density", "both"}

def test_stagnation_cap_binds_the_temperature():
    # a tiny energy radius with a huge drag would heat past the stagnation temperature without the cap
    e = EmissionState(**{**E.__dict__, "U": 3000.0, "rho": 1e-2, "p": 1e4, "cda": 50.0})
    amb = Ambient(p=10.0, T=240.0, rho=10.0 / (R_AIR * 240.0))
    s = gas_state(0.0, 0.1, e, X, amb, "density")
    assert s.T == pytest.approx(240.0 + 3000.0**2 / (2 * CP_AIR), rel=1e-9)


# ---- Plume.wake_state ----

U = 6000.0
TURB = dict(rho=1e-2, T=250.0, p=300.0, mu=1.6e-5, mfp=1e-5, d=2.813, cda=5.5, r0=1.4065)
AMB_P = Ambient(p=300.0, T=250.0, rho=300.0 / (R_AIR * 250.0))

def line_parent(flow=TURB, T=20.0, n=2001):
    t = np.linspace(0, T, n); r = np.c_[U*t, 0*t, 0*t]
    return Track.from_arrays(t, r, v=np.tile([U, 0, 0], (n, 1)), flow={k: np.full(n, v) for k, v in flow.items()})

def test_wake_state_inside_the_plume_is_hot_forward_and_consistent_with_the_margin():
    p = Plume(line_parent()); x = np.array([U * 10.0 - 100.0, 2.0, 0.0])
    s, m = p.wake_state(x, 10.0, AMB_P, "both")
    assert m.margin > 0 and isinstance(s, GasState)
    assert m.margin == p.margin(x, 10.0).margin
    assert s.T > AMB_P.T and s.rho < AMB_P.rho and s.u[0] > 0 and s.u[1] == 0.0 and s.u[2] == 0.0
    ref = gas_state(m.distance, m.r_b, m.emission, np.array([1.0, 0, 0]), AMB_P, "both")
    assert (s.T, s.rho) == (ref.T, ref.rho) and np.allclose(s.u, ref.u)

def test_wake_state_is_ambient_when_off_but_still_returns_the_margin():
    p = Plume(line_parent()); x = np.array([U * 10.0 - 100.0, 2.0, 0.0])
    s, m = p.wake_state(x, 10.0, AMB_P, "off")
    assert (s.rho, s.T, s.p) == (AMB_P.rho, AMB_P.T, AMB_P.p) and not s.u.any() and m.margin > 0

def test_wake_state_is_ambient_when_no_law_applies():
    rarefied = dict(TURB, mfp=1.0)                       # Knudsen number 0.36: no wake law without a fallback
    p = Plume(line_parent(flow=rarefied))
    s, m = p.wake_state(np.array([U * 10.0 - 100.0, 2.0, 0.0]), 10.0, AMB_P, "both")
    assert math.isnan(m.margin) and (s.rho, s.T, s.p) == (AMB_P.rho, AMB_P.T, AMB_P.p) and not s.u.any()

def test_wake_state_far_from_the_track_is_exactly_ambient():
    p = Plume(line_parent())
    s, m = p.wake_state(np.array([U * 10.0 - 100.0, 5e4, 0.0]), 10.0, AMB_P, "both")
    assert m.margin < 0 and (s.rho, s.T, s.p) == (AMB_P.rho, AMB_P.T, AMB_P.p) and not s.u.any()

def test_wake_state_passes_the_release_time_through():
    p = Plume(line_parent(), tau_pre=2.0); x = np.array([0.0, 1.0, 0.0])      # still at the release point at t = 10 s
    _, sliding = p.wake_state(x, 10.0, AMB_P, "both")
    s, anchored = p.wake_state(x, 10.0, AMB_P, "both", t_release=0.0)
    assert sliding.distance > 4e4 and anchored.distance < 2.0 and anchored.tau_star == pytest.approx(0.0, abs=1e-3)
    assert anchored.margin == p.margin(x, 10.0, t_release=0.0).margin

def test_wake_state_rejects_an_unknown_mode():
    p = Plume(line_parent())
    with pytest.raises(ValueError, match="wake_flow"):
        p.wake_state(np.array([0.0, 0.0, 0.0]), 10.0, AMB_P, "sideways")
