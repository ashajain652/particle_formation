"""The parent's shock-layer field (Task 7b, spec 4.3d): oblique jumps, Lighthill's nose solution, Billig's shape and
the tabulated field with its smoothed edges and pressure gradient."""
import math, numpy as np, pytest
from reentry_model.gas import EquilibriumAir
from reentry_model.wake_residence.shock_layer import (ShockLayerSettings, FreeStream, oblique_jump, lighthill,
    lighthill_standoff_approx, billig_standoff, billig_vertex_radius, ShockShape, ShockLayerField)

R, R0 = 1.4065, 1.4065
FS = FreeStream(rho=1.85e-5, T=199.0, p=1.05, U=7400.0, mfp=4e-3, m_bar=4.81e-26)   # about 80 km
M = FS.U / math.sqrt(1.4 * 287.05 * FS.T)

@pytest.fixture(scope="module")
def air():
    return EquilibriumAir()

@pytest.fixture(scope="module")
def field(air):
    return ShockLayerField.build(FS, R, R0, ShockLayerSettings(), air=air)

def test_billig_formulas():
    assert billig_standoff(R, M) == pytest.approx(R * 0.143 * math.exp(3.24 / M**2))
    assert billig_vertex_radius(R, M) == pytest.approx(R * 1.143 * math.exp(0.54 / (M - 1) ** 1.2))
    s = ShockShape(R, M, billig_standoff(R, M))
    assert s.x_s(0.0) == pytest.approx(R + billig_standoff(R, M))
    for phi in np.radians([0, 45, 89]):
        assert s.inside(R * math.cos(phi), R * math.sin(phi))

def test_normal_jump_then_compression_gives_stagnation(air):
    stag = air.stagnation(FS.rho, FS.T, FS.U)
    ps = oblique_jump(air, FS, math.pi / 2)
    assert ps.w_t == pytest.approx(0.0, abs=1e-6)
    assert ps.p < stag.p and ps.p == pytest.approx(stag.p, rel=0.05)     # compression to rest adds a few percent
    from reentry_model.gas import AIR
    from scipy.optimize import brentq
    h0 = ps.h + 0.5 * ps.w_n**2
    def res(lp):
        air.gas.SPX = ps.s, math.exp(lp), AIR; air.gas.equilibrate("SP")
        return air.gas.enthalpy_mass - h0
    p0 = math.exp(brentq(res, math.log(ps.p), math.log(2 * ps.p)))
    assert p0 == pytest.approx(stag.p, rel=0.005)

def test_oblique_tangential_preserved_and_weak_limit(air):
    beta = math.asin(1 / M)
    ps = oblique_jump(air, FS, math.radians(60))
    assert ps.w_t == pytest.approx(FS.U * math.cos(math.radians(60)))
    # The brief asked for rho_2 / rho_inf within 5 % of 1 at beta + 0.5 degrees; at Mach 26 that is a normal Mach
    # number of 1.23, a real (frozen-air) shock with a density ratio of 1.38 (measured), so it is checked against the
    # perfect-gas Rankine-Hugoniot ratio instead.
    mn = M * math.sin(beta + math.radians(0.5))
    near = oblique_jump(air, FS, beta + math.radians(0.5))
    assert near.rho / FS.rho == pytest.approx(2.4 * mn**2 / (0.4 * mn**2 + 2.0), rel=0.02)

@pytest.mark.parametrize("eps", [0.05, 0.1, 0.2])
def test_lighthill_conditions_and_standoff(eps):
    L = lighthill(eps, 1.0, 1.0)
    g = lambda r: L.A * r**2 + L.B / r + L.K * r**4
    dg = lambda r: 2 * L.A * r - L.B / r**2 + 4 * L.K * r**3
    assert g(1.0) == pytest.approx(0.0, abs=1e-12)
    assert g(L.R_s) == pytest.approx(-0.5 * eps * L.R_s**2, rel=1e-10)
    assert dg(L.R_s) == pytest.approx(-L.R_s, rel=1e-10)
    assert (L.R_s - 1.0) / L.R_s == pytest.approx(lighthill_standoff_approx(eps), rel=0.10)
    ur, uphi = L.velocity(1.0, math.radians(30))
    assert ur == pytest.approx(0.0, abs=1e-12)                       # body is a streamline

def test_field_stagnation_body_shock(field, air):
    stag = air.stagnation(FS.rho, FS.T, FS.U)
    s = field.state(R + 1e-9, 0.0)
    assert s.p == pytest.approx(stag.p, rel=0.01) and s.T == pytest.approx(stag.T, rel=0.01)
    assert abs(s.w_x) < 1.0 and abs(s.w_y) < 1.0 and s.weight == pytest.approx(1.0)
    phi = math.radians(50); ex, ey = math.cos(phi), math.sin(phi)
    b = field.state(R * ex * (1 + 1e-9), R * ey * (1 + 1e-9))
    assert abs(b.w_x * ex + b.w_y * ey) < 1e-3 * FS.U                 # no flow through the body
    stag_p = stag.p
    assert b.p == pytest.approx(FS.p + (stag_p - FS.p) * math.cos(phi)**2, rel=0.01)

def test_mass_flux_through_layer(field):
    assert field.flux_iterations > 0 and field.flux_mismatch < 0.02
    for phi_deg in (30, 60, 89):
        phi = math.radians(phi_deg); rs = field.shape.r_s(phi)
        r = np.linspace(R, rs, 400); flux = 0.0
        for a, b in zip(r[:-1], r[1:]):
            rm = 0.5 * (a + b); st = field.state(rm * math.cos(phi), rm * math.sin(phi))
            u_phi = -st.w_x * math.sin(phi) + st.w_y * math.cos(phi)
            flux += st.rho * abs(u_phi) * 2 * math.pi * rm * math.sin(phi) * (b - a)
        ys = rs * math.sin(phi)
        assert flux == pytest.approx(math.pi * ys**2 * FS.rho * FS.U, rel=0.02)

def test_outside_is_ambient_and_edges_are_continuous(field):
    x_far = field.shape.x_s(0.0) + 20 * FS.mfp
    assert field.state(x_far, 0.0).weight == 0.0
    assert field.state(-1.01 * R0, 0.5 * R).weight == 0.0
    ws = [field.state(x, 1.2 * R).weight for x in np.linspace(0.05, -0.05, 41)]   # across the base plane
    assert max(abs(a - b) for a, b in zip(ws[:-1], ws[1:])) < 0.05
    xs = field.shape.x_s(0.0); d = FS.mfp * 3
    ws = [field.state(x, 0.0).weight for x in np.linspace(xs - 2 * d, xs + 2 * d, 81)]
    assert max(abs(a - b) for a, b in zip(ws[:-1], ws[1:])) < 0.1

def test_pressure_gradient_matches_finite_difference_and_points_to_body(field):
    x, y, h = R + 0.3 * field.standoff, 0.2, 1e-5
    gx, gy = field.excess_pressure_gradient(x, y)
    pe = lambda x, y: (lambda s: s.weight * (s.p - FS.p))(field.state(x, y))
    assert gx == pytest.approx((pe(x + h, y) - pe(x - h, y)) / (2 * h), rel=1e-3)
    assert gy == pytest.approx((pe(x, y + h) - pe(x, y - h)) / (2 * h), rel=1e-3, abs=1e-6)
    assert field.excess_pressure_gradient(R + 0.5 * field.standoff, 0.0)[0] < 0   # pressure rises towards the body

def test_no_altitude_argument():
    import inspect
    assert not any("alt" in p or p in ("h", "h_km") for p in inspect.signature(ShockLayerField.build).parameters)


# ---- the flux-matched shape (Asha's option 2, spec 4.3d) and the review round's checks ----

def test_billig_shape_carries_too_much_gas_and_is_kept_as_the_sensitivity_case(air):
    # Task 7b measured 1.27 times the inflow at 30 degrees and 1.73 at 60 with Billig's shape; that is why the
    # default shape is flux-matched. The Billig case is unchanged: its shape is Billig's and its mismatch is recorded.
    f = ShockLayerField.build(FS, R, R0, ShockLayerSettings(shock_shape="billig"), air=air)
    assert f.flux_iterations == 0 and f.flux_mismatch > 0.5 and f.flags["flux_mismatch_above_tol"]
    assert f.shape.x_s(1.0) == pytest.approx(ShockShape(R, f.mach, f.standoff).x_s(1.0))
    with pytest.raises(ValueError, match="shock_shape"):
        ShockLayerSettings(shock_shape="lighthill")

def test_flux_matched_shape_keeps_the_axis_standoff(field):
    assert field.standoff == pytest.approx(field.lighthill.R_s - R)
    assert field.shape.x_s(0.0) == pytest.approx(R + field.standoff, rel=1e-12)
    assert field.shape.r_s(0.0) == pytest.approx(R + field.standoff, rel=1e-12)
    assert not field.flags["flux_mismatch_above_tol"]

def test_flux_matched_shock_angle_never_below_the_mach_angle(field):
    beta = math.asin(1.0 / field.mach)
    sig = [field.shape.sigma(y) for y in np.linspace(0.0, 50.0 * R, 2001)]
    assert min(sig) >= beta and sig[0] == pytest.approx(math.pi / 2)
    assert np.all(np.diff(sig) <= 1e-9)                               # falls monotonically from the axis
    assert field.shape.sigma(1e4 * R) == pytest.approx(beta, abs=1e-3)   # tends to the Mach angle far out

def test_flux_matched_shock_distance_has_a_continuous_slope(field):
    h = 1e-6
    rs = field.shape.r_s
    for phi_deg in list(range(3, 90, 3)) + [90, 100, 120]:            # the knots, the base-plane ray and behind it
        p = math.radians(phi_deg)
        left, right = (rs(p) - rs(p - h)) / h, (rs(p + h) - rs(p)) / h
        assert right == pytest.approx(left, abs=1e-4 * max(1.0, abs(left)))
    y90 = field.shape.r_s(math.pi / 2)
    assert field.shape.dx_s(y90 * (1 + 1e-9)) == pytest.approx(field.shape.dx_s(y90 * (1 - 1e-9)), abs=1e-6)

@pytest.mark.parametrize("where", ["shock_band", "base_plane_fade"])
def test_pressure_gradient_matches_finite_difference_across_the_edges(field, where):
    pe = lambda x, y: (lambda s: s.weight * (s.p - FS.p))(field.state(x, y))
    if where == "shock_band":                                         # half-way through the smoothing outside the shock
        y = 0.3
        x = field.shape.x_s(y) + 0.5 * 3 * FS.mfp / math.sin(field.shape.sigma(y))
    else:                                                             # a third of the way into the fade behind the base
        x, y = -0.3 * R0, R + 0.5 * (field.shape.r_s(math.pi / 2) - R)
    w = field.state(x, y).weight
    assert 0.05 < w < 0.95
    h = 1e-6
    gx, gy = field.excess_pressure_gradient(x, y)
    assert gx == pytest.approx((pe(x + h, y) - pe(x - h, y)) / (2 * h), rel=1e-3)
    assert gy == pytest.approx((pe(x, y + h) - pe(x, y - h)) / (2 * h), rel=1e-3, abs=1e-3)

def test_surface_state_rejects_angles_outside_0_to_180(field):
    assert field.surface_state(180.0).p == pytest.approx(field.state(R, 0.0).p)
    for bad in (-1.0, 180.5):
        with pytest.raises(ValueError, match="theta_deg"):
            field.surface_state(bad)
