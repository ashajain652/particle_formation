"""thermal/: the scikit-fem conduction solver against scikit-fem's own assembly and the analytic cases of spec
section 10 (lumped limit, radiative cooling, Carslaw-Jaeger, per-step energy balance)."""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from reentry_model import material, thermal
from reentry_model.thermal import SIGMA_SB, skfem_backend

R = 0.05
RHO, CP, K0, EPS = 2813.0, 900.0, 160.0, 0.4


def constant_material(k=K0, cp=CP, rho=RHO, eps=EPS):
    return material.Material("const", rho, eps, [100.0, 20000.0], [cp, cp], [100.0, 20000.0], [k, k])


def solver(mesh, mat, **options):
    s = thermal.thermal_solver("skfem", **options)
    s.setup(mesh, mat, mat.emissivity)
    return s


def lumped_reference(q, A, m, cp, eps, T0, t_end, T_amb=0.0):
    """Exact lumped body: m c_p dT/dt = q A - eps sigma A (T^4 - T_amb^4)."""
    sol = solve_ivp(lambda t, y: [(q * A - eps * SIGMA_SB * A * (y[0] ** 4 - T_amb ** 4)) / (m * cp)], (0.0, t_end), [T0],
                    rtol=1e-12, atol=1e-10)
    return float(sol.y[0, -1])


def test_factory_and_options():
    assert isinstance(thermal.thermal_solver("skfem"), skfem_backend.SkfemThermalSolver)
    with pytest.raises(ValueError):
        thermal.thermal_solver("nope")
    with pytest.raises(ValueError):
        thermal.thermal_solver("skfem", linear_solver="magic")
    with pytest.raises(ValueError):
        thermal.thermal_solver("skfem", max_iterations=0)


def test_operators_match_scikit_fem_assembly(coarse_sphere_mesh):
    """K (element-mean k) against scikit-fem's assembly; the lumped capacity matrix is diag(sum_e V_e/4 rho c_p(T_i))
    (the nodal-enthalpy lumping, not the row sums of the consistent matrix), and the consistent option matches
    scikit-fem's consistent mass with the nodal (P1) c_p."""
    import scipy.sparse as sp
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json())
    T = 300.0 + 400.0 * np.random.default_rng(1).random(coarse_sphere_mesh.n_nodes)
    K, M = s.operators(T)
    K_ref, M_ref = s.reference_operators(T)
    nodal = np.bincount(s.tets.ravel(), np.repeat(s.vol / 4.0, 4) * (RHO * s.material.cp(T))[s.tets].ravel(), minlength=s.points.shape[0])
    assert abs(K - K_ref).max() < 1e-10 * abs(K_ref).max() and abs(M - sp.diags(nodal)).max() < 1e-10 * abs(M_ref).max()
    _, M_c = solver(coarse_sphere_mesh, material.Material.from_drama_json(), lumped_mass=False).operators(T)
    assert abs(M_c - M_ref).max() < 1e-10 * abs(M_ref).max()
    assert abs(np.asarray(K.sum(axis=1))).max() < 1e-9 * abs(K).max()          # rows of K sum to zero
    assert M.sum() == pytest.approx((RHO * s.material.cp(T[s.tets]).mean(axis=1) * s.vol).sum(), rel=1e-12)   # nodal c_p, V_e/4 per node
    u = 2.0 * s.points[:, 0] + 3.0 * s.points[:, 1] - s.points[:, 2]           # linear field: u^T K u = k |grad u|^2 V
    K1, _ = solver(coarse_sphere_mesh, constant_material(k=1.0)).operators(T)
    assert u @ (K1 @ u) == pytest.approx(14.0 * s.vol.sum(), rel=1e-10)


def test_lumped_limit_matches_the_lumped_ode(coarse_sphere_mesh):
    """k x 1e4 (Bi -> 0), uniform flux 1e6 W/m2 with radiation, dt 0.5 s for 50 s: T_eq vs the lumped ODE. The lumped
    body uses the discrete sphere's own area and volume. (k x 1e6 makes CG fail on conditioning: keep 1e4.)"""
    s = solver(coarse_sphere_mesh, constant_material(k=K0 * 1e4), newton_tol=1e-8)
    s.set_temperature(300.0)
    A, m = s.areas.sum(), RHO * s.vol.sum()
    q = np.full(len(s.faces), 1e6)
    for _ in range(100):
        res = s.step(0.5, q, 0.0)
    assert res.Q_conv == pytest.approx(1e6 * A, rel=1e-12) and 1 <= res.iterations <= 6
    T_eq = s.energy() / (m * CP) + material.T_REF
    assert T_eq == pytest.approx(lumped_reference(1e6, A, m, CP, EPS, 300.0, 50.0), rel=1e-3)      # spec: 0.1 %
    assert s.T.max() - s.T.min() < 0.1                                                              # isothermal


def test_radiative_cooling_of_an_isothermal_sphere(coarse_sphere_mesh):
    """No flux, k x 1e4, T0 1500 K, dt 0.5 s for 100 s: T(t) = [T0^-3 + 3 eps sigma A t / (m c_p)]^(-1/3), 0.1 %."""
    s = solver(coarse_sphere_mesh, constant_material(k=K0 * 1e4), newton_tol=1e-8)
    s.set_temperature(1500.0)
    A, m = s.areas.sum(), RHO * s.vol.sum()
    for _ in range(200):
        res = s.step(0.5, np.zeros(len(s.faces)), 0.0)
    exact = (1500.0 ** -3 + 3.0 * EPS * SIGMA_SB * A * 100.0 / (m * CP)) ** (-1.0 / 3.0)
    assert s.energy() / (m * CP) + material.T_REF == pytest.approx(exact, rel=1e-3)
    assert res.Q_rad == pytest.approx(EPS * SIGMA_SB * A * exact ** 4, rel=2e-3) and res.Q_conv == 0.0


def test_carslaw_jaeger_step_change_at_the_surface(uniform_test_mesh):
    """Surface stepped from 300 to 800 K and held (Dirichlet), constant properties, dt 0.05 s: centre temperature vs the
    series 2 sum (-1)^(n+1) exp(-n^2 pi^2 alpha t / R^2) within 1 % of the 500 K rise, volume mean within 0.5 %."""
    s = solver(uniform_test_mesh, constant_material(eps=0.0), newton_tol=1e-10)
    s.set_temperature(300.0)
    alpha = K0 / (RHO * CP)
    boundary = uniform_test_mesh.boundary_nodes()
    centre, V = uniform_test_mesh.centre_node(), s.vol.sum()

    def series(t, coefficient):
        return sum(coefficient(n) * math.exp(-n * n * math.pi ** 2 * alpha * t / R ** 2) for n in range(1, 80))

    for i in range(1, 81):
        s.step(0.05, np.zeros(len(s.faces)), 0.0, dirichlet=(boundary, np.full(boundary.size, 800.0)))
        t = 0.05 * i
        if i in (40, 80):
            T_centre_exact = 800.0 - 500.0 * 2.0 * series(t, lambda n: (-1) ** (n + 1))
            T_mean_exact = 800.0 - 500.0 * (6.0 / math.pi ** 2) * series(t, lambda n: 1.0 / n ** 2)
            T_mean = float((s.vol * s.T[s.tets].mean(axis=1)).sum() / V)
            assert abs(s.T[centre] - T_centre_exact) < 0.01 * 500.0, (t, s.T[centre], T_centre_exact)
            assert abs(T_mean - T_mean_exact) < 0.005 * 500.0, (t, T_mean, T_mean_exact)
    assert np.all(s.T[boundary] == 800.0)


@pytest.mark.parametrize("lumped", [False, True])
def test_energy_balance_every_step(coarse_sphere_mesh, lumped):
    """Constant c_p, cos(theta) flux 2e6 W/m2 on the windward side, radiation on: dE = (Q_conv - Q_rad) dt to 1e-6."""
    s = solver(coarse_sphere_mesh, constant_material(), newton_tol=1e-12, lumped_mass=lumped)
    s.set_temperature(300.0)
    theta = np.arccos(np.clip(coarse_sphere_mesh.surface().normals[:, 0], -1.0, 1.0))
    q = 2e6 * np.where(theta < math.pi / 2, np.cos(theta), 0.0)
    for _ in range(20):
        E0 = s.energy()
        res = s.step(0.5, q, 0.0)
        assert abs((s.energy() - E0) - (res.Q_conv - res.Q_rad) * 0.5) < 1e-6 * abs((res.Q_conv - res.Q_rad) * 0.5)
    assert 700.0 < s.T.max() < 800.0 and 305.0 < s.T.min() < 320.0 and res.Q_rad > 0.0


def test_variable_properties_energy_balance_is_exact_too(coarse_sphere_mesh):
    """With the AA7075 tables the secant heat capacity keeps dE = (Q_conv - Q_rad) dt to the Newton tolerance."""
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json(), newton_tol=1e-12)
    s.set_temperature(300.0)
    q = np.full(len(s.faces), 5e5)
    for _ in range(40):
        E0 = s.energy()
        res = s.step(0.5, q, 0.0)
        assert abs((s.energy() - E0) - (res.Q_conv - res.Q_rad) * 0.5) < 1e-6 * abs((res.Q_conv - res.Q_rad) * 0.5)
    assert s.T.max() > 400.0


def test_direct_and_amg_agree(coarse_sphere_mesh):
    q = np.full(len(coarse_sphere_mesh.surface().faces), 1e6)
    results = []
    for linear in ("direct", "amg"):
        s = solver(coarse_sphere_mesh, material.Material.from_drama_json(), linear_solver=linear, newton_tol=1e-10)
        s.set_temperature(300.0)
        for _ in range(5):
            s.step(0.5, q, 0.0)
        results.append(s.T)
    assert np.abs(results[0] - results[1]).max() < 1e-6


def test_step_result_and_temperature_setters(coarse_sphere_mesh):
    s = solver(coarse_sphere_mesh, constant_material())
    s.set_temperature(np.linspace(300.0, 400.0, coarse_sphere_mesh.n_nodes))
    assert s.temperature()[0] == 300.0 and s.temperature() is not s.T
    res = s.step(0.5, np.zeros(len(s.faces)), 200.0)
    assert isinstance(res, thermal.StepResult) and res.T.shape == (coarse_sphere_mesh.n_nodes,)
    assert res.Q_rad == pytest.approx(s.radiated_power(200.0)) and res.Q_conv == 0.0


def test_newton_non_convergence_raises(coarse_sphere_mesh):
    """Newton loop limited to 1 iteration must fail on a nonlinear step with significant flux."""
    s = solver(coarse_sphere_mesh, constant_material(), max_iterations=1)
    s.set_temperature(300.0)
    with pytest.raises(RuntimeError, match="Newton did not converge"):
        s.step(0.5, np.full(len(s.faces), 1e6), 0.0)
    # Verify that with default max_iterations the same step succeeds and takes >1 iteration
    s2 = solver(coarse_sphere_mesh, constant_material())
    s2.set_temperature(300.0)
    res = s2.step(0.5, np.full(len(s2.faces), 1e6), 0.0)
    assert res.iterations >= 2


# ---------------------------------------------------------------------------------------------------------------
# Step 3: element fractions, pinned nodes, nodal loads, melting

def test_fractions_pinned_nodes_and_nodal_loads_keep_the_balance(coarse_sphere_mesh):
    """Halving phi on the windward owners halves their energy; deactivating some of them pins nothing (their nodes
    still belong to live elements); a nodal sink over six steps keeps the discrete balance to 1e-8."""
    from reentry_model import mesh as mesh_mod
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets)
    s = solver(m, material.Material.from_drama_json("AA7075"))
    s.set_temperature(300.0)
    surf = m.surface()
    q = np.where(surf.normals[:, 0] > 0.0, 2e6, 0.0)
    for _ in range(4):
        s.step(0.5, q, 0.0)
    hot = surf.owner[surf.normals[:, 0] > 0.8]
    E_hot = s.element_energies()[hot].sum()
    phi = np.ones(m.n_elements)
    phi[hot] = 0.5
    E_before = s.energy()
    s.set_fractions(phi)
    assert s.energy() == pytest.approx(E_before - 0.5 * E_hot, rel=1e-12) and not s.pinned.any()
    m.deactivate(hot[:20])
    phi[hot[:20]] = 0.0
    s.set_fractions(phi)
    assert m.n_active == m.n_elements - 20 and len(s.areas) == m.surface().n_patches
    assert not np.isin(np.flatnonzero(s.pinned), np.unique(m.tets[m.active])).any()          # pinned = no live element
    load = np.zeros(m.n_nodes)
    load[np.unique(m.tets[hot[20:40]])] = -50.0                               # a 50 W sink on those nodes
    q2 = np.where(m.surface().normals[:, 0] > 0.0, 2e6, 0.0)
    E0, absorbed = s.energy(), 0.0
    for _ in range(6):
        r = s.step(0.5, q2, 0.0, nodal_load=load)
        absorbed += (r.Q_conv - r.Q_rad + r.Q_extra) * 0.5
        assert r.Q_extra == pytest.approx(load.sum()) and r.Q_dropped == 0.0
    assert abs(s.energy() - E0 - absorbed) < 1e-8 * abs(absorbed)
    with pytest.raises(ValueError):
        s.set_fractions(np.full(m.n_elements, 1.5))


def test_pinned_nodes_hold_their_temperature_and_drop_their_loads():
    """A two-tet mesh: deactivating one tet pins its private node; the load on it is reported as dropped."""
    from reentry_model import mesh as mesh_mod
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float) * 1e-2
    m = mesh_mod.VolumeMesh(pts, np.array([[0, 1, 2, 3], [1, 2, 3, 4]]))
    s = thermal.thermal_solver("skfem", linear_solver="direct")
    s.setup(m, constant_material(), 0.0)
    s.set_temperature(np.array([300.0, 310.0, 320.0, 330.0, 900.0]))
    m.deactivate([1])
    s.set_fractions(np.array([1.0, 0.0]))
    assert s.pinned.tolist() == [False, False, False, False, True]
    load = np.zeros(5)
    load[4] = 100.0
    r = s.step(0.1, np.zeros(m.surface().n_patches), 0.0, nodal_load=load)
    assert r.Q_dropped == 100.0 and r.Q_extra == 100.0 and s.temperature()[4] == 900.0
    assert s.energy() == pytest.approx(RHO * CP * m.element_volumes()[0] * (315.0 - material.T_REF), rel=1e-9)   # nothing entered the live tet


def test_stefan_front_on_the_box_mesh():
    """Neumann's one-phase solution: liquid at T_w on x < X(t), X = 2 lambda sqrt(alpha t), lambda from
    lambda exp(lambda^2) erf(lambda) = St / sqrt(pi), St = c_p (T_w - T_m) / L. Box 30 x 1 x 1 mm of 0.5 mm cubes,
    wall Dirichlet at x = 0, front from the liquid-fraction contour at 0.5: within 1 % at t = 2..5 s (spec 13.4;
    measured 0.3 % on 2026-09-20)."""
    from scipy.optimize import brentq
    from scipy.special import erf
    from reentry_model import mesh as mesh_mod
    L, TM, TW = 4e5, 850.0, 1250.0
    mat = material.Material("stefan", RHO, 0.0, [100.0, 20000.0], [CP, CP], [100.0, 20000.0], [K0, K0], L, TM, TM,
                            material.LiquidProperties(2400.0, 1.3e-3, 0.86))
    alpha, St = K0 / (RHO * CP), CP * (TW - TM) / L
    lam = brentq(lambda x: x * np.exp(x * x) * erf(x) - St / np.sqrt(np.pi), 1e-6, 5.0)
    m = mesh_mod.box_mesh(0.03, 1e-3, 1e-3, 0.5e-3)
    s = solver(m, mat, linear_solver="direct")
    s.set_temperature(mat.T_solidus)
    wall = np.flatnonzero(m.points[:, 0] < 1e-9)
    axis = np.flatnonzero((np.abs(m.points[:, 1]) < 1e-12) & (np.abs(m.points[:, 2]) < 1e-12))
    x = m.points[axis, 0]
    order = np.argsort(x)
    t, dt = 0.0, 0.05
    while t < 5.0 - 1e-9:
        r = s.step(dt, np.zeros(m.surface().n_patches), 0.0, dirichlet=(wall, np.full(wall.size, TW)))
        t += dt
        assert r.iterations <= 8
        if abs(t - round(t)) < 1e-9 and t >= 2.0:
            f = mat.liquid_fraction(s.temperature())[axis][order]
            x_front = np.interp(0.5, f[::-1], x[order][::-1])
            assert x_front == pytest.approx(2.0 * lam * np.sqrt(alpha * t), rel=1e-2)


def test_melting_iteration_converges_across_the_ramp(coarse_sphere_mesh):
    """The isothermal body crossing the +-2 K ramp of the single-temperature material (the case that cycled with the
    element-mean secant iteration, 2026-09-20): every step converges in a few iterations and the enthalpy balance
    holds through the latent-heat plateau."""
    mat = material.Material.from_drama_json("AA7075")
    mat.k_table = mat.k_table * 1e4
    s = solver(coarse_sphere_mesh, mat)
    s.set_temperature(840.0)
    q = np.full(coarse_sphere_mesh.surface().n_patches, 8e5)
    E0, absorbed = s.energy(), 0.0
    for _ in range(40):
        r = s.step(0.5, q, 0.0)
        absorbed += (r.Q_conv - r.Q_rad) * 0.5
        assert r.iterations <= 8
    T = s.temperature()
    assert 848.0 < T.mean() < 852.0 and T.max() - T.min() < 0.5                 # on the plateau, isothermal
    assert abs(s.energy() - E0 - absorbed) < 1e-8 * absorbed


def test_film_mass_rides_the_boundary_nodes_with_the_solid_s_capacity(coarse_sphere_mesh):
    """A melt film handed to the solver as nodal mass: it joins the nodes' capacity with the material's own c_p, so
    the same heat raises the body less; its nodes stay live when their elements die; and the solver's nodal capacity
    reports both halves (Step 3, the film is thermally thin and has no temperature of its own)."""
    from reentry_model import mesh as mesh_mod
    m = mesh_mod.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets)
    mat = material.Material.from_drama_json("AA7075")
    s = solver(m, mat)
    s.set_temperature(300.0)
    surf = m.surface()
    q = np.where(surf.normals[:, 0] > 0.0, 5e5, 0.0)
    for _ in range(3):
        s.step(0.5, q, 0.0)
    T_dry = s.temperature().max()
    film = np.zeros(m.n_nodes)
    film[np.unique(surf.faces)] = 2.0e-4                                       # 0.2 g per boundary node
    C_dry = s.nodal_capacity()
    s.set_film_mass(film)
    assert s.nodal_capacity() == pytest.approx(C_dry + film * mat.cp_eff(s.temperature()), rel=1e-12)
    s.set_temperature(300.0)
    for _ in range(3):
        s.step(0.5, q, 0.0)
    assert s.temperature().max() < T_dry                                       # the film's capacity absorbs its share
    with pytest.raises(ValueError):
        s.set_film_mass(np.full(m.n_nodes, -1.0))
    hot = surf.owner[surf.normals[:, 0] > 0.8]
    m.deactivate(hot)
    phi = np.ones(m.n_elements)
    phi[hot] = 0.0
    s.set_fractions(phi)
    wet = np.setdiff1d(np.unique(m.tets[hot]), np.unique(m.tets[m.active]))     # nodes left with no live element
    assert wet.size and not s.pinned[wet].any()                                # ... but with film: still live
    s.set_film_mass(np.zeros(m.n_nodes))
    assert s.pinned[wet].all()                                                 # film gone: pinned again, in one breath
