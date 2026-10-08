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
        thermal.thermal_solver("fenicsx", lumped_mass=False)           # the nodal-enthalpy (lumped) capacity only


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


def test_melting_run_matches_the_skfem_backend(coarse_sphere_mesh):
    """Element fractions, pinned nodes, nodal loads and the enthalpy Newton in both backends: 6 s of a melting body
    under physics-mode loads give the same mass, sprayed mass and temperatures (measured 1e-10 / 0 K on 2026-09-20)."""
    pytest.importorskip("cantera")
    from reentry_model import body, heating, mesh
    from test_reentry_model_coupled import MASS_100MM, simulator
    out = {}
    for name in ("skfem", "fenicsx"):
        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM, T0=700.0)
        sim = simulator(b, t_max=60.0)
        sim.advance(43.0)
        model = heating.PhysicsHeating()
        for _ in range(12):
            sim.advance(0.5)
            a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
            b.advance(sim.t, 0.5, model.evaluate(a, b.theta, b.surface_temperature(), 0.05, T_mean=b.mean_temperature()), state=a)
        out[name] = (b.mass(0.0), b.sprayed_mass, b.solver.temperature(), b.mesh.n_active, b.energy_balance_residual(),
                     b.m_f.sum(), b.film_energy(), b.solver.film_mass.sum(), b.solver.film_weight().max())
    assert out["skfem"][0] == pytest.approx(out["fenicsx"][0], rel=1e-6) and out["skfem"][1] == pytest.approx(out["fenicsx"][1], rel=1e-4)
    assert np.abs(out["skfem"][2] - out["fenicsx"][2]).max() < 0.5 and out["skfem"][3] == out["fenicsx"][3]
    assert abs(out["fenicsx"][4]) < 1e-6 and out["fenicsx"][1] > 0.0
    # the film rides the boundary nodes in both backends, with the same mass, the same liquid enthalpy and the same
    # share of the nodes it owns (which is what the enthalpy Newton inverts against)
    assert out["skfem"][5] == pytest.approx(out["fenicsx"][5], rel=1e-4) and out["fenicsx"][5] > 0.0
    assert out["skfem"][6] == pytest.approx(out["fenicsx"][6], rel=1e-4)
    assert out["fenicsx"][7] == pytest.approx(out["skfem"][7], rel=1e-4) and 0.0 < out["fenicsx"][7] <= 1.0


def test_both_backends_accept_facet_temperature_with_no_argument(coarse_sphere_mesh):
    """The two backends stand behind one protocol, so a call that works on either must work on both:
    `facet_temperature()` with no argument means the solver's own stored field. The FEniCSx backend required the field
    explicitly, which nothing caught because all four callers pass it -- a trap with no upside (fixed 2026-09-24)."""
    mat = material.Material.from_drama_json()
    for name in ("skfem", "fenicsx"):
        s = solver(name, coarse_sphere_mesh, mat)
        s.set_temperature(450.0)
        assert s.facet_temperature() == pytest.approx(np.full(len(s.areas), 450.0)), name
        assert s.facet_temperature() == pytest.approx(s.facet_temperature(s.temperature())), name


def test_the_molten_cascade_matches_the_skfem_backend(coarse_sphere_mesh):
    """The molten cascade (amendment of 2026-10-03) in both backends: a pool 6 mm deep under the nose at 69.8 km and six
    coupled steps with the physics loads -- deaths expose fully molten elements and the cascade feeds them within the
    step, while the deep runoff and the spray act -- give the same cascade passes, the same active set, the same mass,
    sprayed mass, cascade mass and deep account, and the same temperatures, with both energy balances exact. numpy's
    generator is seeded for each backend because pyamg draws its starting vectors from it (plan fact 52). Measured on
    2026-10-03: passes 4, 3, 16, 4, 5, 4 in both; mass 1.3e-7, sprayed 6.3e-7, cascade mass 2.4e-7, deep account 3.7e-6,
    temperatures 0.11 K, balances 4.7e-10 and 2.4e-11 -- the stiff deep transport carries the backends' last-bit
    differences further (plan fact 49); with the deep runoff off they agree to 1e-10 and 4e-6 K. With the runoff flux
    of 2026-10-05 the passes are the same and the backends agree to 7.6e-10 in mass, 3.6e-9 in sprayed mass, 1.0e-10 in
    cascade mass, 6.9e-11 in the deep account and 2.5e-5 K: what carried the last-bit differences was the film's
    runaway emptying rates on thick patches, not the deep transport (plan fact 72)."""
    pytest.importorskip("cantera")
    from reentry_model import body, heating, mesh
    from test_reentry_model_coupled import MASS_100MM, simulator
    out = {}
    for name in ("skfem", "fenicsx"):
        np.random.seed(12345)
        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM)
        sim = simulator(b, t_max=90.0)
        sim.advance(50.0)
        p = b.mesh.points
        r = np.linalg.norm(p, axis=1)
        b.solver.set_temperature(np.where((r > 0.05 - 6e-3) & (p[:, 0] > 0.5 * r), 960.0, 850.0))
        b.energy0 = b.energy()
        model = heating.PhysicsHeating()
        passes = []
        for _ in range(6):
            sim.advance(0.5)
            a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
            b.advance(sim.t, 0.5, model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature()),
                      state=a)
            passes.append(b.last_melt["cascade_passes"])
        out[name] = (b, passes)
    (s, ps), (f, pf) = out["skfem"], out["fenicsx"]
    assert ps == pf and max(pf) >= 2 and f.cascade_mass > 0.0 and f.deep_runoff_mass > 0.0
    assert np.array_equal(s.mesh.active, f.mesh.active)
    assert f.cascade_mass == pytest.approx(s.cascade_mass, rel=1e-6) and f.mass(0.0) == pytest.approx(s.mass(0.0), rel=1e-6)
    assert f.sprayed_mass == pytest.approx(s.sprayed_mass, rel=1e-5) and f.m_d.sum() == pytest.approx(s.m_d.sum(), rel=1e-4)
    assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 0.5
    assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8


def test_the_thick_film_flux_matches_the_skfem_backend(coarse_sphere_mesh, monkeypatch):
    """The runoff flux on a thick patch (amendment of 2026-10-05) in both backends: the 6 mm pool at 69.8 km of the
    cascade test above, at ten steps of 0.05 s, where every film transport has thick patches carrying films thinner
    than delta_m -- the case the amendment changes. Both backends give the same active set, mass, sprayed mass, runoff,
    film, deep account and temperatures, with both energy balances exact, and the largest emptying rate of any edge
    stays bounded: before the amendment it grew from sub-step to sub-step and reached 3.9e25 per second in this setting
    (measured 2026-10-05); with it 2.7e7, set by the deepest piles' G b^3 / (3 mu) term, not by a vanishing film.
    Measured: the same 7 983 active elements, mass 5.8e-10, sprayed mass 5.9e-9, runoff 7.5e-10, film 1.7e-8, deep
    account 2.8e-10, temperatures 2.3e-4 K, balances 5.1e-9 and 1.1e-11.
    The rigid substrate (amendment of 2026-10-06) is off here: with the body's interior at 850 K, above AA7075_range's
    829 K half-liquid point, the whole body is slurry, every wet windward patch takes the thick branch, and films of
    about a micron that straddle the spray floor B_MIN are released whole by one backend and not by the other from the
    sixth step on -- 8 to 29 patches, the runoff then differing by 1.6e-4 and the temperatures by 0.29 K (measured
    2026-10-06). That is the model's threshold, not the backends'; the test after this one checks the two with it on."""
    pytest.importorskip("cantera")
    from reentry_model import body, film, heating, mesh
    from test_reentry_model_coupled import MASS_100MM, simulator
    seen = {"thin_on_thick": 0, "c_max": 0.0}
    lub, coeff = film.lubrication, film.Runoff.edge_coefficients

    def lubrication(tau, G, b, delta_m, mu_l, b_layer=None):
        out = lub(tau, G, b, delta_m, mu_l, b_layer)
        bb = np.asarray(b, dtype=float)
        seen["thin_on_thick"] += int((out[3] & (bb > 0.0) & (bb < np.asarray(delta_m, dtype=float))).sum())
        return out

    def edge_coefficients(self, q, b, t_hat, areas):
        c_ij, c_ji = coeff(self, q, b, t_hat, areas)
        if len(c_ij):
            seen["c_max"] = max(seen["c_max"], float(c_ij.max()), float(c_ji.max()))
        return c_ij, c_ji

    monkeypatch.setattr(film, "lubrication", lubrication)
    monkeypatch.setattr(film.Runoff, "edge_coefficients", edge_coefficients)
    out = {}
    for name in ("skfem", "fenicsx"):
        np.random.seed(12345)
        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM,
                             settings=body.MeltSettings(rigid_substrate=False))
        sim = simulator(b, t_max=90.0)
        sim.advance(50.0)
        p = b.mesh.points
        r = np.linalg.norm(p, axis=1)
        b.solver.set_temperature(np.where((r > 0.05 - 6e-3) & (p[:, 0] > 0.5 * r), 960.0, 850.0))
        b.energy0 = b.energy()
        model = heating.PhysicsHeating()
        for _ in range(10):
            sim.advance(0.05)
            a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
            b.advance(sim.t, 0.05, model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature()),
                      state=a)
        out[name] = b
    s, f = out["skfem"], out["fenicsx"]
    assert seen["thin_on_thick"] > 0 and seen["c_max"] < 1e9 and f.runoff_mass > 0.0 and f.sprayed_mass > 0.0
    assert np.array_equal(s.mesh.active, f.mesh.active)
    assert f.mass(0.0) == pytest.approx(s.mass(0.0), rel=1e-6) and f.sprayed_mass == pytest.approx(s.sprayed_mass, rel=1e-5)
    assert f.runoff_mass == pytest.approx(s.runoff_mass, rel=1e-5) and f.m_f.sum() == pytest.approx(s.m_f.sum(), rel=1e-4)
    assert f.m_d.sum() == pytest.approx(s.m_d.sum(), rel=1e-4)
    assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 0.5
    assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8


def test_the_rigid_substrate_matches_the_skfem_backend(coarse_sphere_mesh):
    """The rigid substrate (amendment of 2026-10-06) in both backends: the 6 mm pool at 69.8 km of the test above, with
    the rest of the body at 820 K, below AA7075_range's 829 K half-liquid point, so slurry forms only where the pool
    and the heating raise it -- ten steps of 0.05 s. The non-rigid depth reaches into the body, both backends make the
    same branch decision on every patch at every step, and the active set, mass, sprayed mass, runoff, film and
    temperatures agree with both energy balances exact. Measured: mass 6.7e-10, sprayed mass 1.0e-8, runoff 1.2e-8,
    film 6.1e-8, temperatures 1.8e-4 K."""
    pytest.importorskip("cantera")
    from reentry_model import body, heating, mesh
    from test_reentry_model_coupled import MASS_100MM, simulator
    out = {}
    for name in ("skfem", "fenicsx"):
        np.random.seed(12345)
        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver(name), MASS_100MM)
        assert b.settings.rigid_substrate
        sim = simulator(b, t_max=90.0)
        sim.advance(50.0)
        p = b.mesh.points
        r = np.linalg.norm(p, axis=1)
        b.solver.set_temperature(np.where((r > 0.05 - 6e-3) & (p[:, 0] > 0.5 * r), 960.0, 820.0))
        b.energy0 = b.energy()
        model = heating.PhysicsHeating()
        branches = []
        for _ in range(10):
            sim.advance(0.05)
            a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
            b.advance(sim.t, 0.05, model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature()),
                      state=a)
            branches.append(None if b.last_spray is None else b.last_spray.branch.copy())
        out[name] = (b, branches)
    (s, bs), (f, bf) = out["skfem"], out["fenicsx"]
    assert s.last_nonrigid is not None and s.last_nonrigid.max() > 6e-3 and f.sprayed_mass > 0.0
    assert all((x is None and y is None) or np.array_equal(x, y) for x, y in zip(bs, bf))
    assert np.array_equal(s.mesh.active, f.mesh.active)
    assert f.mass(0.0) == pytest.approx(s.mass(0.0), rel=1e-6) and f.sprayed_mass == pytest.approx(s.sprayed_mass, rel=1e-5)
    assert f.runoff_mass == pytest.approx(s.runoff_mass, rel=1e-5) and f.m_f.sum() == pytest.approx(s.m_f.sum(), rel=1e-4)
    assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 0.5
    assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8
