"""thermal/fenicsx_backend.py: the same analytic checks as the scikit-fem backend and a cross-check against it.
Skipped unless dolfinx is importable (run with the fenicsx_env interpreter, spec section 12)."""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from reentry_model import material, thermal
from reentry_model.thermal import SIGMA_SB

pytest.importorskip("dolfinx")

R = 0.05
RHO, CP, K0, EPS = 2813.0, 900.0, 160.0, 0.4


def constant_material(k=K0, cp=CP, rho=RHO, eps=EPS):
    return material.Material("const", rho, eps, [100.0, 20000.0], [cp, cp], [100.0, 20000.0], [k, k])


def solver(name, mesh, mat, **options):
    s = thermal.thermal_solver(name, **options)
    s.setup(mesh, mat, mat.emissivity)
    return s


def test_constructor_options():
    assert thermal.thermal_solver("fenicsx").linear_solver == "amg"
    with pytest.raises(ValueError):
        thermal.thermal_solver("fenicsx", lumped_mass=True)


def test_temperature_round_trip_uses_the_mesh_node_order(coarse_sphere_mesh):
    s = solver("fenicsx", coarse_sphere_mesh, constant_material())
    T = 300.0 + 100.0 * coarse_sphere_mesh.points[:, 0] / R
    s.set_temperature(T)
    assert np.allclose(s.temperature(), T)
    Te = T[coarse_sphere_mesh.tets].mean(axis=1)                       # the cell set is the same, whatever dolfinx's cell order
    expected = RHO * CP * float((coarse_sphere_mesh.element_volumes() * (Te - material.T_REF)).sum())
    assert s.energy() == pytest.approx(expected, rel=1e-6)


def test_lumped_limit_matches_the_lumped_ode(coarse_sphere_mesh):
    s = solver("fenicsx", coarse_sphere_mesh, constant_material(k=K0 * 1e4), linear_solver="direct", newton_tol=1e-8)
    s.set_temperature(300.0)
    A, m = s.areas.sum(), RHO * s.cell_volumes.sum()
    q = np.full(len(s.faces), 1e6)
    for _ in range(100):
        res = s.step(0.5, q, 0.0)
    exact = solve_ivp(lambda t, y: [(1e6 * A - EPS * SIGMA_SB * A * y[0] ** 4) / (m * CP)], (0.0, 50.0), [300.0], rtol=1e-12, atol=1e-10).y[0, -1]
    assert s.energy() / (m * CP) + material.T_REF == pytest.approx(exact, rel=1e-3) and res.Q_conv == pytest.approx(1e6 * A, rel=1e-12)


def test_radiative_cooling_of_an_isothermal_sphere(coarse_sphere_mesh):
    s = solver("fenicsx", coarse_sphere_mesh, constant_material(k=K0 * 1e4), linear_solver="direct", newton_tol=1e-8)
    s.set_temperature(1500.0)
    A, m = s.areas.sum(), RHO * s.cell_volumes.sum()
    for _ in range(200):
        s.step(0.5, np.zeros(len(s.faces)), 0.0)
    exact = (1500.0 ** -3 + 3.0 * EPS * SIGMA_SB * A * 100.0 / (m * CP)) ** (-1.0 / 3.0)
    assert s.energy() / (m * CP) + material.T_REF == pytest.approx(exact, rel=2e-3)


def test_energy_balance_every_step(coarse_sphere_mesh):
    s = solver("fenicsx", coarse_sphere_mesh, material.Material.from_drama_json(), newton_tol=1e-12)
    s.set_temperature(300.0)
    theta = np.arccos(np.clip(coarse_sphere_mesh.surface().normals[:, 0], -1.0, 1.0))
    q = 2e6 * np.where(theta < math.pi / 2, np.cos(theta), 0.0)
    for _ in range(20):
        E0 = s.energy()
        res = s.step(0.5, q, 0.0)
        assert abs((s.energy() - E0) - (res.Q_conv - res.Q_rad) * 0.5) < 1e-6 * abs((res.Q_conv - res.Q_rad) * 0.5)


def test_carslaw_jaeger_with_dirichlet(uniform_test_mesh):
    s = solver("fenicsx", uniform_test_mesh, constant_material(eps=0.0), newton_tol=1e-10)
    s.set_temperature(300.0)
    alpha = K0 / (RHO * CP)
    boundary, centre = uniform_test_mesh.boundary_nodes(), uniform_test_mesh.centre_node()
    for i in range(1, 81):
        s.step(0.05, np.zeros(len(s.faces)), 0.0, dirichlet=(boundary, np.full(boundary.size, 800.0)))
    exact = 800.0 - 500.0 * 2.0 * sum((-1) ** (n + 1) * math.exp(-n * n * math.pi ** 2 * alpha * 4.0 / R ** 2) for n in range(1, 80))
    assert abs(s.temperature()[centre] - exact) < 0.01 * 500.0 and np.all(s.temperature()[boundary] == 800.0)


def test_cross_check_against_scikit_fem(coarse_sphere_mesh):
    """Same mesh, same AA7075 material, same cos(theta) loads for 20 steps of 0.5 s: energy-equivalent temperature,
    surface extremes and every nodal temperature agree to 0.1 %."""
    mat = material.Material.from_drama_json()
    theta = np.arccos(np.clip(coarse_sphere_mesh.surface().normals[:, 0], -1.0, 1.0))
    q = 1.5e6 * np.where(theta < math.pi / 2, np.cos(theta), 0.0)
    fields = {}
    for name in ("skfem", "fenicsx"):
        s = solver(name, coarse_sphere_mesh, mat, newton_tol=1e-10)
        s.set_temperature(300.0)
        for _ in range(20):
            s.step(0.5, q, 0.0)
        fields[name] = (s.temperature(), s.energy())
    T1, E1 = fields["skfem"]
    T2, E2 = fields["fenicsx"]
    assert abs(E2 / E1 - 1.0) < 1e-3 and np.abs(T2 - T1).max() / T1.max() < 1e-3
    assert abs(T2.max() / T1.max() - 1.0) < 1e-3 and abs(T2.min() / T1.min() - 1.0) < 1e-3
