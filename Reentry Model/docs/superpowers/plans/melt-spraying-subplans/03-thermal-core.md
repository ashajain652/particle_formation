# Sub-plan: Task 3 — Thermal core: nodal enthalpy, element fractions, pinned nodes, nodal loads

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 1109–2042). Read `00-shared-context.md` first.

**Depends on:** Task 1 (mesh, active set) and Task 2 (material, enthalpy functions).
**Produces, for later tasks:** the rewritten thermal solver — nodal enthalpy, a lumped capacity matrix, an enthalpy-consistent Newton iteration, pinned dead nodes, and deferred/nodal loads — that Task 9's melting step is built on top of. Must behave identically on both the skfem and FEniCSx backends behind one shared protocol.
**Character:** core numerics — the single highest technical risk in the whole plan.
**Read before implementing:** Measured fact 3 in the shared context: the previous Step 2 approach (element-mean secant iteration) could cycle indefinitely across the latent-heat plateau instead of converging, which is why this task exists as a rewrite rather than an extension. The CG solver needs an AMG-then-direct fallback, because melting can empty elements down to near-zero conductance and break the preconditioner.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan. Explicitly include a cross-backend test (skfem vs. FEniCSx) in the test plan — the master plan requires both backends to "stand behind one protocol," not just one working implementation.

---

### Task 3: Thermal core — nodal enthalpy, element fractions, pinned nodes, nodal loads

**Files:**
- Modify (replace): `reentry_model/thermal/__init__.py`, `reentry_model/thermal/skfem_backend.py`, `reentry_model/thermal/fenicsx_backend.py`
- Modify: `reentry_model/cli.py` (the `--lumped-mass` flag becomes `--consistent-mass`, three lines)
- Test: `tests/test_reentry_model_thermal.py` (one test changed, five appended), `tests/test_reentry_model_fenicsx.py` (one line changed, two tests appended)

**Interfaces:**
- Consumes: `VolumeMesh.active`, `.surface()` (Task 1); `Material.enthalpy/enthalpy_liquid/enthalpy_mixed/cp/cp_eff/cp_mixed/temperature_from_enthalpy/temperature_from_enthalpy_mixed/melts/T_solidus/T_liquidus` (Task 2).
- Produces: `StepResult(T, Q_conv, Q_rad, iterations, Q_extra=0.0, Q_dropped=0.0)`; solver methods `set_fractions(phi)` (0 = dead; refreshes the boundary and the pinned nodes), `element_energies(T=None)` (φ_e ρ V_e mean_i h(T_i)), `step(dt, q_conv, T_amb, dirichlet=None, nodal_load=None)` (nodal_load in W, mesh node order); attributes `phi`, `pinned`, `faces`, `areas`, `last_damping`; `SkfemThermalSolver(lumped_mass=True)` default with `consistent_mass = not lumped_mass`; `element_matrices(points, tets) -> (vol, Ke, Mk)` with `Mk[e, k]` the nodal-coefficient mass matrices and `nodal_mass_weights()`; `mass_matrix(c_nodal, c_film=None)`; `operators(T, T_old=None) -> (K, M_tan[, E])`; and, for the melt film (Step 3, Task 9): `set_film_mass(mass)` (one non-negative value per node, kg -- it joins the nodes' capacity and keeps a node live even when its elements are gone), `nodal_capacity()` (J/K per node, material + film, which the body uses to bound the melt loads it defers) and `film_weight()` (the film's share of each node's mass, which the enthalpy-consistent Newton update inverts against). `facet_temperature(T=None)` (the mean of a facet's three nodal temperatures, defaulting to the solver's own stored field) is on both backends with the same signature -- they stand behind one protocol, so a call that works on either must work on both, and a test asserts it. The FEniCSx backend has the same public surface (`lumped_mass=False` raises `ValueError`).

The film's capacity is the **liquid** one -- tangent c_p(T) in M, the secant of `enthalpy_liquid` in E -- never the mixture's c_p,eff: the film has already paid its latent heat and must not pay it again on the ramp. The same weighting makes 1^T E the exact increment of (solid nodal enthalpy + film liquid enthalpy), which is what `MeltingBody.energy()` sums, so the coupled balance closes on it. CG gets one fresh AMG hierarchy and then a direct solve if it still stalls (melting drains interior elements to phi ~ 1e-3 with unscaled conduction; `direct_fallbacks` counts it).

- [ ] **Step 1: Update the conformance test and append the melting tests**

In `tests/test_reentry_model_thermal.py` replace the body of `test_operators_match_scikit_fem_assembly` — the lines

```python
def test_operators_match_scikit_fem_assembly(coarse_sphere_mesh):
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json())
    T = 300.0 + 400.0 * np.random.default_rng(1).random(coarse_sphere_mesh.n_nodes)
    K, M = s.operators(T)
    K_ref, M_ref = s.reference_operators(T)
    assert abs(K - K_ref).max() < 1e-10 * abs(K_ref).max() and abs(M - M_ref).max() < 1e-10 * abs(M_ref).max()
    assert abs(np.asarray(K.sum(axis=1))).max() < 1e-9 * abs(K).max()          # rows of K sum to zero
    assert M.sum() == pytest.approx((RHO * s.material.cp(T[s.tets].mean(axis=1)) * s.vol).sum(), rel=1e-12)
```

with

```python
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
```

(the remaining three lines of the test — the linear-field check on `K1` — stay). Then append to the file:

```python
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
```


In `tests/test_reentry_model_fenicsx.py` change the constructor test's line `thermal.thermal_solver("fenicsx", lumped_mass=True)` to `thermal.thermal_solver("fenicsx", lumped_mass=False)           # the nodal-enthalpy (lumped) capacity only` and append:

```python
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
```


(The appended FEniCSx test needs Tasks 5–9; it stays failing in `fenicsx_env` until Task 9 and is skipped in `drama_env`.)

- [ ] **Step 2: Run the thermal tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py -q`
Expected: the conformance test and the four new tests fail (no `set_fractions`, no `nodal_load`, `lumped_mass` default False).

- [ ] **Step 3: Replace `reentry_model/thermal/__init__.py`**

```python
"""Finite-element conduction solvers behind one protocol (spec section 8); `thermal_solver(name)` picks the backend.

Both backends solve rho c_p(T) dT/dt = div(k(T) grad T) with -k grad T . n = eps sigma (T^4 - T_amb^4) - q_conv on the
boundary, backward Euler in time, P1 tetrahedra in space, and report the same quantities: nodal temperatures,
stored enthalpy, radiated power. `dirichlet=(nodes, values)` in `step` is for the analytic tests only."""
from dataclasses import dataclass
from typing import Protocol

import numpy as np

SIGMA_SB = 5.670374419e-8         # Stefan-Boltzmann [W/(m2 K4)]
SOLVER_NAMES = ("skfem", "fenicsx")


class MissingBackend(RuntimeError):
    """The selected backend's library is not importable in this interpreter."""


@dataclass
class StepResult:
    T: np.ndarray                 # nodal temperatures after the step [K]
    Q_conv: float                 # W, convective power applied over the step
    Q_rad: float                  # W, radiated power at the end of the step
    iterations: int               # Newton/Picard iterations taken
    Q_extra: float = 0.0          # W, nodal load applied over the step (Step 3 deferred melt energy)
    Q_dropped: float = 0.0        # W, load that fell on pinned (material-free) nodes and was not applied


class ThermalSolver(Protocol):
    def setup(self, mesh, material, emissivity) -> None: ...
    def set_temperature(self, T) -> None: ...                       # uniform float or nodal array
    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None) -> StepResult: ...
    def temperature(self) -> np.ndarray: ...
    def energy(self) -> float: ...                                  # stored enthalpy above material.T_REF [J]
    def element_energies(self, T=None) -> np.ndarray: ...           # phi_e rho V_e h(T_e) per element [J]
    def radiated_power(self, T_amb) -> float: ...
    def set_fractions(self, phi) -> None: ...                       # element material fractions (0 = dead), Step 3
    def set_film_mass(self, mass) -> None: ...                      # nodal mass of the melt film [kg], Step 3
    def nodal_capacity(self): ...                                   # J/K per node, material + film (Step 3)
    def film_weight(self): ...                                      # the film's share of each node's mass, 0-1 (Step 3)


def thermal_solver(name, **options):
    if name == "skfem":
        from .skfem_backend import SkfemThermalSolver
        return SkfemThermalSolver(**options)
    if name == "fenicsx":
        from .fenicsx_backend import FenicsxThermalSolver      # raises MissingBackend without dolfinx
        return FenicsxThermalSolver(**options)
    raise ValueError("thermal solver must be one of {}, got {!r}".format(SOLVER_NAMES, name))
```


- [ ] **Step 4: Replace `reentry_model/thermal/skfem_backend.py`**

```python
"""scikit-fem backend: P1 tetrahedra, backward Euler, Newton on the nodal enthalpy with the radiation term.

Assembly. The P1 stiffness depends on T through one coefficient per element (k at the element-mean temperature)
and the mass through a nodal coefficient (rho c at the nodes, interpolated linearly), so the unit-coefficient
element matrices are computed once and rescaled into a fixed CSR pattern on every iteration (~20 ms for 50 k
tets). scikit-fem's generic `asm` on the same MeshTet/ElementTetP1 costs 0.2 s per iteration at 12 k nodes
(measured 2026-09-18) and is kept as the reference in `reference_operators` for the conformance test. Element
fractions phi_e (Step 3) scale both coefficients; nodes without material are pinned at their temperature.

Time stepping. Each iterate solves the tangent system of the residual R(T) = E(T)/dt + K T - F_conv + F_rad(T) -
F_extra, where E = M(rho c_sec)(T - T_old) is the exact nodal enthalpy increment (secant heat capacity per node, so
1^T E equals the change of rho int h dV with latent heat included); the tangent uses c_p,eff(T_i) per node and the
radiation Jacobian. The update is mapped through the true h(T) per node (enthalpy-consistent update) so that a node
cannot jump across the melting range, with damping as a fallback. With a smooth c_p the scheme is the Step 2
secant-capacity iteration unchanged; measured 2026-09-20: 2 iterations per step without melting, 3-5 with.

Radiation. The boundary functional uses the facet-mean temperature: eps sigma (T_f^4 - T_amb^4) A_f/3 to each of
the facet's three nodes, with its exact Jacobian 4 eps sigma T_f^3 A_f/9 on every node pair of the facet. Each
Newton iterate solves the symmetric positive definite system
    (M/dt + K + B_k) T = M T_old/dt + F_conv - F_rad(T_k) + B_k T_k
with M, K frozen at the previous iterate; because the rows of K sum to zero the discrete energy balance
1^T M (T - T_old) = dt (Q_conv - Q_rad) holds to the Newton tolerance.

Linear solver. `direct`: SciPy SuperLU; `amg` (default): CG preconditioned by a pyamg smoothed-aggregation
hierarchy rebuilt every `amg_rebuild_every` solves (the operator changes slowly). Measured on the 100 mm
2 mm/8 mm mesh (18.9 k nodes): AMG 0.03 s per solve vs SuperLU 0.6 s."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from . import SIGMA_SB, StepResult


class _Pattern:
    """Fixed CSR sparsity for repeated (rows, cols) entries; `assemble` sums values into it with one bincount."""

    def __init__(self, rows, cols, n):
        key = rows.astype(np.int64) * n + cols.astype(np.int64)
        unique, self._map = np.unique(key, return_inverse=True)
        self._indptr = np.searchsorted(unique // n, np.arange(n + 1))
        self._indices = (unique % n).astype(np.int32)
        self.n, self.nnz = n, unique.size

    def assemble(self, values):
        data = np.bincount(self._map, weights=np.asarray(values).ravel(), minlength=self.nnz)
        return sp.csr_matrix((data, self._indices, self._indptr), shape=(self.n, self.n))


def nodal_mass_weights():
    """W[k, a, b] = int phi_a phi_b phi_k dV / V_e on a tetrahedron: 1/20 (a = b = k), 1/60 (two equal), 1/120 (all
    distinct); sum_k W[k] is the consistent mass V_e/20 (1 + delta_ab)."""
    W = np.full((4, 4, 4), 1.0 / 120.0)
    for k in range(4):
        for a in range(4):
            W[k, a, a] = 1.0 / 60.0
            W[k, k, a] = W[k, a, k] = 1.0 / 60.0
        W[k, k, k] = 1.0 / 20.0
    return W


def element_matrices(points, tets):
    """Per-element volumes, unit stiffness (grad phi_a . grad phi_b V_e) and the mass matrices for a nodal coefficient:
    M_e(c) = sum_k c_k Mk[e, k] (Mk = V_e W), so that M_e(1) is the consistent mass V_e/20 (1 + delta_ab)."""
    x = points[tets]
    J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
    vol = np.abs(np.linalg.det(J)) / 6.0
    grad_ref = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    grads = np.einsum("eij,aj->eai", np.linalg.inv(J).transpose(0, 2, 1), grad_ref)
    Ke = np.einsum("eai,ebi->eab", grads, grads) * vol[:, None, None]
    Mk = nodal_mass_weights()[None] * vol[:, None, None, None]
    return vol, Ke, Mk


class SkfemThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=True, newton_tol=1e-6, max_iterations=30,
                 amg_rebuild_every=30, cg_tol=1e-10):
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        if max_iterations < 1:
            raise ValueError("max_iterations must be >= 1")
        self.linear_solver, self.lumped_mass = linear_solver, lumped_mass
        self.consistent_mass = not lumped_mass
        self.newton_tol, self.max_iterations = newton_tol, max_iterations
        self.amg_rebuild_every, self.cg_tol = amg_rebuild_every, cg_tol
        self._ml, self._solves, self.last_cg_iterations, self.direct_fallbacks = None, 0, 0, 0

    def setup(self, mesh, material, emissivity):
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        self.points, self.tets = mesh.points, mesh.tets
        self.vol, self.Ke, self.Mk = element_matrices(self.points, self.tets)
        n = len(self.points)
        self.pattern = _Pattern(np.repeat(self.tets, 4, axis=1).ravel(), np.tile(self.tets, (1, 4)).ravel(), n)
        self.T, self._T_prev = np.full(n, 300.0), None
        self.phi = np.where(mesh.active, 1.0, 0.0)
        self.film_mass = np.zeros(n)              # kg per node: the melt film riding on the boundary (Step 3)
        self.pinned = np.zeros(n, dtype=bool)
        self._refresh_surface()

    def _refresh_surface(self):
        """Facet data (loads, radiation pattern) of the mesh's current boundary and the pinned (material-free) nodes."""
        surface = self.mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        n = len(self.points)
        self.facet_pattern = _Pattern(np.repeat(self.faces, 3, axis=1).ravel(), np.tile(self.faces, (1, 3)).ravel(), n)
        self.Bf = np.ones((3, 3))[None] * (self.areas / 9.0)[:, None, None]
        self._update_pinned()

    def _update_pinned(self):
        """Nodes with neither material nor film: they keep their temperature (identity rows). Recomputed whenever
        either the fractions or the film capacity change -- a node that had film and lost it must be pinned again in
        the same breath, or its row empties and the system goes singular (measured 2026-09-22)."""
        pinned = np.ones(len(self.points), dtype=bool)
        pinned[np.unique(self.tets[self.phi > 0.0])] = False
        pinned &= self.film_mass <= 0.0                          # a node carrying film still has a heat capacity
        if not np.array_equal(pinned, self.pinned):
            self._ml = None                                   # the operator's structure changed: fresh AMG hierarchy
        self.pinned = pinned

    def set_film_mass(self, mass):
        """Nodal mass [kg] of the melt film riding on the boundary (Step 3). The film is thermally thin -- q b / k_l =
        0.22 K across a 10 um film at 2 MW/m2, and b^2/alpha = 0.3 ms against a 0.5 s macro step -- so it is given no
        temperature of its own: its mass joins the boundary nodes' capacity with the *same* heat capacity the solid
        uses (tangent c_p,eff in the operator, secant [h(T) - h(T_old)]/(T - T_old) in the enthalpy rate), so the film's
        sensible heat is part of 1^T M dT by construction and `MeltingBody.film_energy` closes the balance exactly."""
        mass = np.asarray(mass, dtype=float)
        if mass.shape != (len(self.points),) or (mass < 0.0).any() or not np.isfinite(mass).all():
            raise ValueError("film mass must be one finite non-negative value per node")
        self.film_mass = mass
        self._update_pinned()

    def set_fractions(self, phi):
        """Element material fractions phi_e in [0, 1] (0 = dead): they scale the heat capacity of every element (the
        conductivity stays that of the full element while phi > 0, see operators()); the boundary follows the mesh's
        active set and nodes without material are pinned at their temperature. Call after the mesh's active set or
        the fractions change."""
        self.phi = np.asarray(phi, dtype=float)
        if self.phi.shape != (len(self.tets),) or (self.phi < 0.0).any() or (self.phi > 1.0).any():
            raise ValueError("fractions must be one value in [0, 1] per element")
        self._refresh_surface()

    def set_temperature(self, T):
        self.T = np.full(len(self.points), float(T)) if np.ndim(T) == 0 else np.array(T, dtype=float)
        self._T_prev = None

    def temperature(self):
        return self.T.copy()

    def facet_temperature(self, T=None):
        return (self.T if T is None else T)[self.faces].mean(axis=1)

    def facet_load(self, q):
        """Nodal load vector of a per-facet flux q [W/m2]: A_f/3 to each of the facet's nodes."""
        return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=len(self.points))

    def mass_matrix(self, c_nodal, c_film=None):
        """Capacity matrix for the nodal coefficient c [J/(m3 K)] (times phi_e per element). Lumped by default:
        diag(sum_e phi_e V_e/4 c_i), so that 1^T M dT = sum_e phi_e V_e/4 sum_i c_i dT_i is exactly the increment of
        the nodal enthalpy integral (the consistent form sum_k c_k Mk integrates the product of the interpolants of c
        and dT, which is not the increment of any energy functional and left a 3e-4 balance error, measured
        2026-09-20). `consistent_mass=True` keeps the consistent form (analytic tests only)."""
        # the film is liquid: its capacity is the liquid c_p (tangent) or the secant of the liquid enthalpy, never
        # the mixture's c_p,eff -- it has already paid its latent heat and does not pay it again on the ramp
        film = self.film_mass * (c_nodal if c_film is None else c_film) / self.material.rho
        if self.consistent_mass:
            M = self.pattern.assemble(np.einsum("ek,ekab->eab", self.phi[:, None] * c_nodal[self.tets], self.Mk))
            return M + sp.diags(film, format="csr") if self.film_mass.any() else M
        return sp.diags(np.bincount(self.tets.ravel(), weights=np.repeat(self.phi * self.vol / 4.0, 4) * c_nodal[self.tets].ravel(),
                                    minlength=len(self.points)) + film, format="csr")

    def operators(self, T, T_old=None):
        """Stiffness K with k(T_e) at the element-mean temperature (active elements); the tangent mass M_tan with the nodal
        coefficient rho c_p,eff(T_i); and the enthalpy-rate vector E = M(rho c_sec) (T - T_old) with the nodal secant
        heat capacity c_sec,i = [h(T_i) - h(T_old,i)] / (T_i - T_old,i) (c_p,eff where a node has not moved). The
        enthalpy is nodal (spec Step 3 section 6, decided 2026-09-20: the element-mean enthalpy of Step 2 released the
        latent heat over a 4 K window of the element mean while the nodal temperatures span 30 K across a surface
        element, and its Newton iteration cycled): 1^T E is exactly the increment of rho int h dV with h interpolated
        linearly, whatever h(T) is. Without T_old, E is None and M_tan carries c_p,eff (the reference-operator convention)."""
        Te = T[self.tets].mean(axis=1)
        # conductivity is NOT scaled by phi_e: a partly consumed element is a thinner sliver of the same material,
        # which conducts better, not worse; scaling k with phi isolated the surface nodes of nearly consumed
        # elements and drove them to 5000 K (measured 2026-09-20). Dead elements (phi = 0) drop out.
        K = self.pattern.assemble(self.Ke * ((self.phi > 0.0) * self.material.k(Te))[:, None, None])
        c_tan, c_liq = self.material.cp_eff(T), self.material.cp(T)
        M = self.mass_matrix(self.material.rho * c_tan, self.material.rho * c_liq)
        if T_old is None:
            return K, M
        dT = T - T_old
        moved = np.abs(dT) > 1e-9
        c_sec = np.where(moved, (self.material.enthalpy(T) - self.material.enthalpy(T_old)) / np.where(moved, dT, 1.0), c_tan)
        c_sec_l = np.where(moved, (self.material.enthalpy_liquid(T) - self.material.enthalpy_liquid(T_old)) / np.where(moved, dT, 1.0), c_liq)
        # the film's mass enters both M and E through mass_matrix with its own liquid capacity, so 1^T E is exactly
        # the increment of (solid nodal enthalpy + film liquid enthalpy) -- which is what MeltingBody.energy sums
        return K, M, self.mass_matrix(self.material.rho * c_sec, self.material.rho * c_sec_l) @ dT

    def radiated_power(self, T_amb, T=None):
        Tf = self.facet_temperature(T)
        return float((self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas).sum())

    def film_weight(self):
        """Share of each node's mass that is melt film, 0 to 1. The film is liquid: its enthalpy has no latent
        plateau, so the enthalpy-consistent Newton update -- which inverts the *mixture* h(T) -- must not be applied
        to a node the film owns, or the node is pushed onto a plateau it is not on and the iteration cycles (Newton
        stopped converging at 30 iterations once a patch was mostly film, measured 2026-09-22). The update is blended
        with the plain linear one by this weight, which is exact in both limits; it changes only the iteration, never
        the residual it converges to."""
        solid = np.bincount(self.tets.ravel(), weights=np.repeat(self.phi * self.vol / 4.0, 4) * self.material.rho,
                            minlength=len(self.points))
        total = solid + self.film_mass
        return np.divide(self.film_mass, total, out=np.zeros_like(total), where=total > 0.0)

    def nodal_capacity(self):
        """Heat capacity carried by each node [J/K], material and film together (lumped, tangent c_p,eff). The body
        uses it to bound the melt loads it defers: energy booked onto a node that has melted away has nothing to
        heat, and pushing it in anyway moved drained surface nodes by 1e4-1e5 K in a step (measured 2026-09-22)."""
        c = self.material.cp_eff(self.T) * self.material.rho
        return np.bincount(self.tets.ravel(), weights=np.repeat(self.phi * self.vol / 4.0, 4) * c[self.tets].ravel(),
                           minlength=len(self.points)) + self.film_mass * self.material.cp(self.T)

    def energy(self):
        """Stored enthalpy above T_REF: rho int h dV with the nodal h(T_i) interpolated linearly, i.e.
        sum_e phi_e rho V_e mean_i h(T_i); equals 1^T M T for constant c_p (consistent or lumped mass)."""
        return float(self.element_energies().sum())

    def element_energies(self, T=None):
        """phi_e rho V_e mean_i h(T_i) per element [J]."""
        h = self.material.enthalpy(self.T if T is None else T)
        return self.phi * self.material.rho * self.vol * h[self.tets].mean(axis=1)

    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None):
        """One backward-Euler step with the per-facet convective flux q_conv [W/m2] and an optional nodal load vector
        [W] (Step 3: the deferred melt energy; loads on pinned nodes are dropped and reported in StepResult.Q_dropped)."""
        T_old = self.T
        # predictor: extrapolate the previous step (saves ~1 Newton iteration per step); plain T_old on the first step
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))
        F_extra = np.zeros(len(self.points)) if nodal_load is None else np.asarray(nodal_load, dtype=float)
        Q_dropped = float(F_conv[self.pinned].sum() + F_extra[self.pinned].sum())
        T_new, iteration = T_k, 0
        converged = False
        last_relative_change, previous_change, damping, decreases = None, None, 1.0, 0
        mat = self.material
        film_w = self.film_weight() if self.film_mass.any() else np.zeros(len(self.points))
        for iteration in range(1, self.max_iterations + 1):
            K, M, E = self.operators(T_k, T_old)
            Tf = self.facet_temperature(T_k)
            B = self.facet_pattern.assemble(self.Bf * (4.0 * self.emissivity * SIGMA_SB * Tf ** 3)[:, None, None])
            F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
            # Newton on R(T) = E(T)/dt + K T - F_conv + F_rad(T) - F_extra with the tangent M/dt + K + B:
            # A T_new = A T_k - R(T_k); with a smooth c_p (c_tan = c_sec) this is the Step 2 iteration unchanged
            A = (M / dt + K + B).tocsr()
            b = M @ T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra
            if self.pinned.any():                              # material-free nodes keep their temperature (identity rows)
                A = A + sp.diags(self.pinned.astype(float), format="csr")
                b[self.pinned] = T_old[self.pinned]
            T_new = self._solve_with_dirichlet(A, b, T_k, dirichlet)
            if mat.melts:
                # enthalpy-consistent update: the linearised step is an enthalpy increment c_p,eff(T_k) (T_new - T_k)
                # per node; inverting the true h(T) puts a node that would overshoot the melting range where the
                # latent heat actually leaves it (identity where h is linear). Without it nodes jump across the +-2 K
                # ramp of a single-temperature material and the iteration cycles (measured 2026-09-20). The enthalpy
                # inverted is the node's own mixture of material and film (see film_weight), not the material's.
                free, w = ~self.pinned, film_w[~self.pinned]
                T_new[free] = mat.temperature_from_enthalpy_mixed(
                    mat.enthalpy_mixed(T_k[free], w) + mat.cp_mixed(T_k[free], w) * (T_new[free] - T_k[free]), w)
            last_relative_change = np.linalg.norm(T_new - T_k) / np.linalg.norm(T_new)
            if last_relative_change <= self.newton_tol:
                converged = True
                break
            # damping when the change grows (a fallback for cycling iterates), released again after two decreases
            if previous_change is not None and last_relative_change > previous_change:
                damping, decreases = max(0.25, 0.5 * damping), 0
            elif damping < 1.0:
                decreases += 1
                if decreases >= 2:
                    damping, decreases = 1.0, 0
            previous_change = last_relative_change
            T_k = T_k + damping * (T_new - T_k) if damping < 1.0 else T_new
        if not converged:
            raise RuntimeError("Newton did not converge in {} iterations (last relative change {:.2e}, tol {:.1e})".format(
                self.max_iterations, last_relative_change, self.newton_tol))
        self._T_prev, self.T = T_old, T_new
        self.last_damping = damping
        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()), Q_dropped)

    def _solve_with_dirichlet(self, A, b, x0, dirichlet):
        if dirichlet is None:
            return self._solve(A, b, x0)
        nodes, values = dirichlet
        free = np.ones(A.shape[0], dtype=bool)
        free[nodes] = False
        x = np.zeros(A.shape[0])
        x[nodes] = values
        x[free] = self._solve(A[free][:, free], b[free] - A[free][:, ~free] @ x[~free], x0[free], fresh=True)
        return x

    def _solve(self, A, b, x0, fresh=False):
        if self.linear_solver == "direct":
            return spla.spsolve(A.tocsc(), b)
        for attempt in (0, 1):
            if fresh or attempt or self._ml is None or self._solves % self.amg_rebuild_every == 0:
                import pyamg
                self._ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
            self._solves += 1
            counter = []
            x, info = spla.cg(A, b, x0=x0, rtol=self.cg_tol, maxiter=500, M=self._ml.aspreconditioner(cycle="V"),
                              callback=lambda _: counter.append(1))
            self.last_cg_iterations = len(counter)
            if info == 0:
                return x
        # Melting drains interior elements to phi ~ 1e-3 while their conduction stays unscaled (spec section 9), which
        # raises the condition number by 1/phi; a hierarchy built before the drain can stall on the emptied region.
        # A fresh hierarchy clears it in every case measured except the step after the heating is switched off, where
        # the operator changes by orders of magnitude in one step -- that one system is solved directly.
        self.direct_fallbacks += 1
        return spla.spsolve(A.tocsc(), b)

    def reference_operators(self, T):
        """K (element-mean k) and M (nodal c_p) assembled by scikit-fem with the same coefficients (conformance test only)."""
        from skfem import Basis, BilinearForm, ElementTetP0, ElementTetP1, MeshTet, asm
        from skfem.helpers import dot, grad
        basis = Basis(MeshTet(self.points.T.copy(), self.tets.T.copy()), ElementTetP1(), intorder=3)   # exact for phi_i phi_j c (cubic)
        basis0 = basis.with_element(ElementTetP0())
        Te = T[self.tets].mean(axis=1)

        @BilinearForm
        def stiffness(u, v, w):
            return w.k * dot(grad(u), grad(v))

        @BilinearForm
        def mass(u, v, w):
            return w.c * u * v

        K = asm(stiffness, basis, k=basis0.interpolate(self.material.k(Te)))
        M = asm(mass, basis, c=basis.interpolate(self.material.rho * self.material.cp_eff(T)))    # nodal (P1) coefficient
        return K, M
```


- [ ] **Step 5: Replace `reentry_model/thermal/fenicsx_backend.py`**

```python
"""FEniCSx (dolfinx >= 0.11) backend: the scheme of skfem_backend on the same mesh (spec sections 8 and Step 3 6).

dolfinx assembles the stiffness K = int phi_e k(T_e) grad u . grad v dx (a DG0 coefficient refreshed every Newton
iterate); everything nodal -- the lumped capacity diag(sum_e phi_e V_e/4 rho c_p,eff(T_i)), the enthalpy-rate vector
E = diag(sum_e phi_e V_e/4 rho c_sec,i)(T - T_old), the convective loads A_f/3 per facet node, the radiation
eps sigma (T_f^4 - T_amb^4) A_f/3 with its Jacobian 4 eps sigma T_f^3 A_f/9 on the facet node pairs, and the Step 3
nodal loads -- is built in numpy on the VolumeMesh's numbering and added to the PETSc operator (the facet pairs are
edges of cells, so they lie inside K's pattern). Both backends therefore discretise identically; the radiation
follows the current boundary of the active set (element death), which a UFL `ds` measure could not. Dirichlet
nodes (analytic tests) and pinned material-free nodes are imposed with MatZeroRowsColumns. PETSc CG with hypre
BoomerAMG (hierarchy reused for `amg_rebuild_every` solves) or LU. dolfinx renumbers vertices: `node_of_dof` /
`dof_of_node` map between the VolumeMesh's node ids and the P1 dofs. Serial (the nodal vectors assume one
process). Verified with dolfinx 0.11.0 in fenicsx_env (2026-09-20). `dolfinx` is imported lazily: the constructor
raises MissingBackend without it."""
import numpy as np
import scipy.sparse as sp

from . import SIGMA_SB, MissingBackend, StepResult


class FenicsxThermalSolver:
    def __init__(self, linear_solver="amg", lumped_mass=True, newton_tol=1e-6, max_iterations=30, amg_rebuild_every=30,
                 cg_tol=1e-10, **_):
        try:
            import dolfinx  # noqa: F401
            import ufl  # noqa: F401
            from mpi4py import MPI  # noqa: F401
            from petsc4py import PETSc  # noqa: F401
        except ImportError as exc:
            raise MissingBackend("the fenicsx backend needs dolfinx, which is not importable here: create the separate "
                                 "conda environment fenicsx_env (spec section 12) and run with its interpreter") from exc
        if not lumped_mass:
            raise ValueError("the fenicsx backend implements the lumped (nodal-enthalpy) capacity only")
        if linear_solver not in ("direct", "amg"):
            raise ValueError("linear_solver must be direct or amg, got {!r}".format(linear_solver))
        if max_iterations < 1:
            raise ValueError("max_iterations must be >= 1")
        self.linear_solver, self.newton_tol, self.max_iterations, self.cg_tol = linear_solver, newton_tol, max_iterations, cg_tol
        self.amg_rebuild_every, self._solves = amg_rebuild_every, 0
        self.lumped_mass = True

    def setup(self, mesh, material, emissivity):
        import basix.ufl
        import ufl
        from dolfinx import fem
        from dolfinx import mesh as dmesh
        from dolfinx.fem import petsc
        from mpi4py import MPI
        from petsc4py import PETSc
        from scipy.spatial import cKDTree
        self.mesh, self.material, self.emissivity = mesh, material, float(emissivity)
        domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
        self.msh = dmesh.create_mesh(MPI.COMM_WORLD, mesh.tets.astype(np.int64), domain, mesh.points)   # dolfinx >= 0.9: (comm, cells, element, x)
        self.V = fem.functionspace(self.msh, ("Lagrange", 1))
        self.V0 = fem.functionspace(self.msh, ("DG", 0))
        n_local = self.V.dofmap.index_map.size_local
        _, self.node_of_dof = cKDTree(mesh.points).query(self.V.tabulate_dof_coordinates()[:n_local])
        self.dof_of_node = np.empty(n_local, dtype=np.int64)
        self.dof_of_node[self.node_of_dof] = np.arange(n_local)
        self.cell_dofs = np.asarray(self.V.dofmap.list)[:, :4]
        # dolfinx cell k is the k-th cell passed in (create_mesh keeps the order in serial): cell volumes and the map
        # from dolfinx cells to mesh elements
        x = self.V.tabulate_dof_coordinates()[self.cell_dofs]
        J = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        self.cell_volumes = np.abs(np.linalg.det(J)) / 6.0
        self.cell_nodes = self.node_of_dof[self.cell_dofs]                   # mesh node ids of each dolfinx cell
        centroids = mesh.points[mesh.tets].mean(axis=1)
        self.element_of_cell = cKDTree(centroids).query(mesh.points[self.cell_nodes].mean(axis=1))[1]
        self.coef = fem.Function(self.V0)                                   # phi_e k(T_e) per cell
        u, v = ufl.TrialFunction(self.V), ufl.TestFunction(self.V)
        self.a = fem.form(self.coef * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx)
        self.K = petsc.create_matrix(self.a)
        self.ksp = PETSc.KSP().create(self.msh.comm)
        if self.linear_solver == "direct":
            self.ksp.setType("preonly")
            self.ksp.getPC().setType("lu")
        else:
            self.ksp.setType("cg")
            self.ksp.getPC().setType("hypre")
            self.ksp.getPC().setHYPREType("boomeramg")
            self.ksp.setTolerances(rtol=self.cg_tol, max_it=500)
        self.T = np.full(mesh.n_nodes, 300.0)                               # nodal, mesh numbering
        self._T_prev = None
        self.phi = np.where(mesh.active, 1.0, 0.0)
        self.film_mass = np.zeros(mesh.n_nodes)                   # kg per node (see the skfem backend)
        self.pinned = np.zeros(mesh.n_nodes, dtype=bool)
        self._structure_changed = True
        self._refresh_surface()

    # -- geometry / fractions -----------------------------------------------------------------------------------
    def _refresh_surface(self):
        surface = self.mesh.surface()
        self.faces, self.areas = surface.faces, surface.areas
        n = self.mesh.n_nodes
        rows, cols = np.repeat(self.faces, 3, axis=1).ravel(), np.tile(self.faces, (1, 3)).ravel()
        self._facet_rows, self._facet_cols = rows, cols
        pinned = np.ones(n, dtype=bool)
        pinned[np.unique(self.mesh.tets[self.phi > 0.0])] = False
        pinned &= self.film_mass <= 0.0                           # a node carrying film keeps a heat capacity
        if not np.array_equal(pinned, self.pinned):
            self._structure_changed = True
        self.pinned = pinned

    def set_film_mass(self, mass):
        """Nodal mass [kg] of the melt film riding on the boundary (see the skfem backend)."""
        mass = np.asarray(mass, dtype=float)
        if mass.shape != (self.mesh.n_nodes,) or (mass < 0.0).any() or not np.isfinite(mass).all():
            raise ValueError("film mass must be one finite non-negative value per node")
        self.film_mass = mass
        self._refresh_surface()

    def set_fractions(self, phi):
        self.phi = np.asarray(phi, dtype=float)
        if self.phi.shape != (self.mesh.n_elements,) or (self.phi < 0.0).any() or (self.phi > 1.0).any():
            raise ValueError("fractions must be one value in [0, 1] per element")
        self._refresh_surface()

    def set_temperature(self, T):
        self.T = np.full(self.mesh.n_nodes, float(T)) if np.ndim(T) == 0 else np.array(T, dtype=float)
        self._T_prev = None

    def temperature(self):
        return self.T.copy()

    # -- nodal pieces (mesh numbering) -----------------------------------------------------------------------------
    def facet_load(self, q):
        return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=self.mesh.n_nodes)

    def facet_temperature(self, T=None):
        return (self.T if T is None else T)[self.faces].mean(axis=1)   # same default as the skfem backend: one protocol

    def lumped(self, c_nodal, c_film=None):
        """diag(sum_e phi_e V_e/4 c_i) on the mesh nodes, plus the film's own (liquid) capacity where it rides."""
        vol = self.mesh.element_volumes()
        return np.bincount(self.mesh.tets.ravel(), weights=np.repeat(self.phi * vol / 4.0, 4) * c_nodal[self.mesh.tets].ravel(),
                           minlength=self.mesh.n_nodes) + self.film_mass * (c_nodal if c_film is None else c_film) / self.material.rho

    def element_energies(self, T=None):
        h = self.material.enthalpy(self.T if T is None else T)
        return self.phi * self.material.rho * self.mesh.element_volumes() * h[self.mesh.tets].mean(axis=1)

    def film_weight(self):
        """Share of each node's mass that is melt film, 0 to 1 (see the skfem backend)."""
        solid = np.bincount(self.mesh.tets.ravel(),
                            weights=np.repeat(self.phi * self.mesh.element_volumes() / 4.0, 4) * self.material.rho,
                            minlength=self.mesh.n_nodes)
        total = solid + self.film_mass
        return np.divide(self.film_mass, total, out=np.zeros_like(total), where=total > 0.0)

    def nodal_capacity(self):
        """Heat capacity carried by each node [J/K], material and film together (see the skfem backend)."""
        T = self.temperature()
        c, vol = self.material.cp_eff(T) * self.material.rho, self.mesh.element_volumes()
        return np.bincount(self.mesh.tets.ravel(), weights=np.repeat(self.phi * vol / 4.0, 4) * c[self.mesh.tets].ravel(),
                           minlength=self.mesh.n_nodes) + self.film_mass * self.material.cp(T)

    def energy(self):
        return float(self.element_energies().sum())

    def radiated_power(self, T_amb, T=None):
        Tf = self.facet_temperature(self.T if T is None else T)
        return float((self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas).sum())

    def _assemble_stiffness(self, T):
        from dolfinx.fem.petsc import assemble_matrix
        Te = T[self.mesh.tets].mean(axis=1)
        coef = (self.phi > 0.0) * self.material.k(Te)                     # k unscaled by phi (see skfem_backend.operators)
        self.coef.x.array[:] = coef[self.element_of_cell]
        self.K.zeroEntries()
        assemble_matrix(self.K, self.a)
        self.K.assemble()

    def _petsc_from_csr(self, matrix):
        from petsc4py import PETSc
        m = matrix.tocsr()
        m.sum_duplicates()
        return PETSc.Mat().createAIJ(size=m.shape, csr=(m.indptr.astype(np.int32), m.indices.astype(np.int32), m.data), comm=self.msh.comm)

    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None):
        from petsc4py import PETSc
        mat, n = self.material, self.mesh.n_nodes
        T_old = self.T
        T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
        F_conv = self.facet_load(np.asarray(q_conv, dtype=float))
        F_extra = np.zeros(n) if nodal_load is None else np.asarray(nodal_load, dtype=float)
        Q_dropped = float(F_conv[self.pinned].sum() + F_extra[self.pinned].sum())
        fixed = self.pinned.copy()
        fixed_values = T_old.copy()
        if dirichlet is not None:
            fixed[np.asarray(dirichlet[0])] = True
            fixed_values[np.asarray(dirichlet[0])] = dirichlet[1]
        vol = self.mesh.element_volumes()
        film_w = self.film_weight() if self.film_mass.any() else np.zeros(n)
        converged, last_relative_change, previous_change, damping, decreases = False, None, None, 1.0, 0
        T_new, iteration = T_k, 0
        for iteration in range(1, self.max_iterations + 1):
            self._assemble_stiffness(T_k)
            c_tan = mat.cp_eff(T_k)
            dT = T_k - T_old
            moved = np.abs(dT) > 1e-9
            c_sec = np.where(moved, (mat.enthalpy(T_k) - mat.enthalpy(T_old)) / np.where(moved, dT, 1.0), c_tan)
            c_liq = mat.cp(T_k)                                   # the film is liquid: no latent plateau (skfem backend)
            c_sec_l = np.where(moved, (mat.enthalpy_liquid(T_k) - mat.enthalpy_liquid(T_old)) / np.where(moved, dT, 1.0), c_liq)
            M_tan = self.lumped(mat.rho * c_tan, mat.rho * c_liq)
            E = self.lumped(mat.rho * c_sec, mat.rho * c_sec_l) * dT
            Tf = self.facet_temperature(T_k)
            B = sp.csr_matrix((np.repeat(4.0 * self.emissivity * SIGMA_SB * Tf ** 3 * self.areas / 9.0, 9), (self._facet_rows, self._facet_cols)), shape=(n, n))
            F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
            b = M_tan * T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra
            # operator in dof numbering: K (dolfinx) + diag(M_tan/dt) + B
            A = self.K.copy()
            extra = sp.diags(M_tan / dt) + B
            extra = extra.tocsr()[self.node_of_dof][:, self.node_of_dof]
            A.axpy(1.0, self._petsc_from_csr(extra), structure=PETSc.Mat.Structure.SUBSET_NONZERO_PATTERN)
            rhs = PETSc.Vec().createWithArray(b[self.node_of_dof].copy(), comm=self.msh.comm)
            x = PETSc.Vec().createWithArray(fixed_values[self.node_of_dof].copy(), comm=self.msh.comm)
            fixed_dofs = self.dof_of_node[np.flatnonzero(fixed)].astype(np.int32)
            if fixed_dofs.size:
                A.zeroRowsColumns(fixed_dofs, diag=1.0, x=x, b=rhs)
            self.ksp.setOperators(A)
            reuse = self.linear_solver == "amg" and dirichlet is None and not self._structure_changed and self._solves % self.amg_rebuild_every != 0
            self.ksp.getPC().setReusePreconditioner(reuse)
            self._structure_changed = False
            self._solves += 1
            sol = PETSc.Vec().createWithArray(T_k[self.node_of_dof].copy(), comm=self.msh.comm)
            self.ksp.solve(rhs, sol)
            if self.ksp.getConvergedReason() <= 0:
                raise RuntimeError("PETSc KSP did not converge (reason {})".format(self.ksp.getConvergedReason()))
            T_new = sol.getArray()[self.dof_of_node].copy()
            T_new[fixed] = fixed_values[fixed]
            if mat.melts:
                free, w = ~fixed, film_w[~fixed]
                T_new[free] = mat.temperature_from_enthalpy_mixed(
                    mat.enthalpy_mixed(T_k[free], w) + mat.cp_mixed(T_k[free], w) * (T_new[free] - T_k[free]), w)
            last_relative_change = np.linalg.norm(T_new - T_k) / np.linalg.norm(T_new)
            if last_relative_change <= self.newton_tol:
                converged = True
                break
            if previous_change is not None and last_relative_change > previous_change:
                damping, decreases = max(0.25, 0.5 * damping), 0
            elif damping < 1.0:
                decreases += 1
                if decreases >= 2:
                    damping, decreases = 1.0, 0
            previous_change = last_relative_change
            T_k = T_k + damping * (T_new - T_k) if damping < 1.0 else T_new
        if not converged:
            raise RuntimeError("Newton did not converge in {} iterations (last relative change {:.2e}, tol {:.1e})".format(
                self.max_iterations, last_relative_change, self.newton_tol))
        self._T_prev, self.T = T_old, T_new
        self.last_damping = damping
        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()), Q_dropped)
```


- [ ] **Step 6: Rename the CLI flag**

In `reentry_model/cli.py`: replace `th.add_argument("--lumped-mass", action="store_true")` with

```python
    th.add_argument("--consistent-mass", action="store_true", help="consistent capacity matrix instead of the lumped nodal-enthalpy one (skfem only, analytic checks)")
```

replace `lumped_mass=args.lumped_mass` (in `build_thermal`) with `lumped_mass=not args.consistent_mass`, and `"lumped_mass": args.lumped_mass}` (the info dict) with `"lumped_mass": not args.consistent_mass}`.

- [ ] **Step 7: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py tests/test_reentry_model_coupled.py tests/test_reentry_model_cli.py tests/test_reentry_model_viz.py -q`
Expected: all pass — the Step 2 analytic cases (lumped limit, radiative cooling, Carslaw–Jaeger, balance) are unchanged by the lumping for constant c_p, and the melting cases converge in ≤ 8 iterations. Then in `fenicsx_env`:

```bash
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q --deselect tests/test_reentry_model_fenicsx.py::test_melting_run_matches_the_skfem_backend
```

Expected: 7 passed (the seven Step 2 conformance tests, including the skfem cross-check, with the nodal scheme).

- [ ] **Step 8: Commit**

```bash
git add reentry_model/thermal tests/test_reentry_model_thermal.py tests/test_reentry_model_fenicsx.py reentry_model/cli.py
git commit -m "Move the thermal core to the nodal enthalpy with element fractions, pinned nodes and nodal loads (Step 3 Task 3)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

