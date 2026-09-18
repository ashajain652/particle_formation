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


def test_operators_match_scikit_fem_assembly(coarse_sphere_mesh):
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json())
    T = 300.0 + 400.0 * np.random.default_rng(1).random(coarse_sphere_mesh.n_nodes)
    K, M = s.operators(T)
    K_ref, M_ref = s.reference_operators(T)
    assert abs(K - K_ref).max() < 1e-10 * abs(K_ref).max() and abs(M - M_ref).max() < 1e-10 * abs(M_ref).max()
    assert abs(np.asarray(K.sum(axis=1))).max() < 1e-9 * abs(K).max()          # rows of K sum to zero
    assert M.sum() == pytest.approx((RHO * s.material.cp(T[s.tets].mean(axis=1)) * s.vol).sum(), rel=1e-12)
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
