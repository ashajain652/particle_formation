"""Tests of the cut-cell thermal solve (design section 7, stage 2): mass rows and weights, the energy books with
radiation and latent heat, second-order convergence on an off-grid slab, the implicit skin base law, the
isothermal limit against the lumped ordinary differential equation, and the molten-rigid surface psi tracking a
one-phase Stefan melting front."""
from math import erf, exp, pi, sqrt

import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import inside_volumes
from reentry_model.levelset.cutthermal import SIGMA_SB, CutThermal
from reentry_model.levelset.surfaces import molten_rigid_surface
from reentry_model.levelset.transfer import interpolate, project


class ConstMaterial:
    def __init__(self, rho=2800.0, k=150.0, cp=900.0, T_ref=293.0):
        self.rho, self._k, self._cp, self.T_ref = rho, k, cp, T_ref

    def k(self, T):
        return np.full_like(np.asarray(T, float), self._k)

    def cp_eff(self, T):
        return np.full_like(np.asarray(T, float), self._cp)

    def enthalpy(self, T):
        return self._cp * (np.asarray(T, float) - self.T_ref)


class LinearMeltMaterial(ConstMaterial):
    """Constant properties plus a latent heat L released linearly between T_s and T_l (a mushy range)."""

    def __init__(self, L=390e3, T_s=880.0, T_l=908.0, **kw):
        super().__init__(**kw)
        self.L, self.T_s, self.T_l = L, T_s, T_l

    def liquid_fraction(self, T):
        return np.clip((np.asarray(T, float) - self.T_s) / (self.T_l - self.T_s), 0.0, 1.0)

    def cp_eff(self, T):
        T = np.asarray(T, float)
        return self._cp + np.where((T > self.T_s) & (T < self.T_l), self.L / (self.T_l - self.T_s), 0.0)

    def enthalpy(self, T):
        return super().enthalpy(T) + self.L * self.liquid_fraction(T)


def sphere(h=2e-3, R=10e-3, L=14e-3):
    m = build(BoxBase.around((-L,) * 3, (L,) * 3, h))
    return m, np.linalg.norm(m.points, axis=1) - R


def slab_column(h, off=0.37, W=4e-3, D=12e-3, top=0.0):
    lo = np.array([-W / 2 - h * off - h, -W / 2 - h * off - h, -D - h * off - h]); hi = -lo; hi[2] = 2 * h
    m = build(BoxBase.around(lo, hi, h))
    x, y, z = m.points.T
    sx = W / 2 + h * off
    return m, np.maximum.reduce([z - top, np.abs(x) - sx, np.abs(y) - sx, -z - D]), sx


def test_mass_rows_and_weights():
    m, phi = sphere()
    s = CutThermal(ConstMaterial()); s.setup(m, phi)
    assert np.allclose(s.w.sum(axis=1), 1.0, atol=1e-12)
    assert np.allclose(s.Me.sum(axis=(1, 2)), s.Vin, rtol=1e-12)
    assert s.metal_volume() == pytest.approx(inside_volumes(m, phi).sum(), rel=1e-13)
    assert s.metal_volume() == pytest.approx(4 / 3 * pi * 1e-6, rel=0.03)


@pytest.mark.parametrize("lumped", [False, True])
def test_energy_books_close_with_radiation_and_latent_heat(lumped):
    m, phi = sphere()
    s = CutThermal(LinearMeltMaterial(), emissivity=0.3, lumped=lumped); s.setup(m, phi); s.set_temperature(860.0)
    E0, Q = s.energy(), 0.0
    for _ in range(6):
        Q += s.step(0.05, 3e6 * np.maximum(-s.P.normal[:, 0], 0), rad=1.0, T_amb=250.0)
    assert (s.energy() - E0 - Q) / abs(Q) == pytest.approx(0.0, abs=1e-6)
    assert s.T[s.active].max() > 880.0                                   # some cells went mushy


def test_constant_flux_slab_converges_at_second_order():
    mat, q, t_end, dt = ConstMaterial(), 2e6, 0.5, 0.01
    alpha = 150 / (2800 * 900)
    exact = 2 * q / 150 * sqrt(alpha * t_end / pi)
    errs = []
    for h in (1e-3, 0.5e-3):
        m, phi, sx = slab_column(h)
        s = CutThermal(mat, linear_solver="direct"); s.setup(m, phi); s.set_temperature(300.0)
        for _ in range(int(round(t_end / dt))):
            s.step(dt, q * np.maximum(s.P.normal[:, 2], 0.0))
        up = s.P.normal[:, 2] > 0.99
        rise = (s.patch_T()[up] * s.P.area[up]).sum() / s.P.area[up].sum() - 300.0
        errs.append(abs(rise - exact) / exact)
    assert errs[0] < 0.016 and errs[1] < 0.005 and errs[0] / errs[1] > 3.0


def test_skin_base_law_is_implicit_and_reaches_its_steady_state():
    m, phi = sphere()
    s = CutThermal(ConstMaterial()); s.setup(m, phi); s.set_temperature(300.0)
    for _ in range(5):
        s.step(50.0, 3.2e5 * 900.0, q1=-3.2e5)                       # q_b = b (900 K - T_b) on every patch
    assert np.allclose(s.T[s.active], 900.0, atol=1e-6)


def test_isothermal_radiative_cooling_matches_the_lumped_ode():
    m, phi = sphere()
    mat = ConstMaterial(k=1e6)                                           # isothermal
    s = CutThermal(mat, emissivity=0.8); s.setup(m, phi); s.set_temperature(1200.0)
    V, A, T, dt = s.metal_volume(), s.P.area.sum(), 1200.0, 2.0
    for _ in range(10):
        s.step(dt, 0.0, rad=1.0, T_amb=0.0)
        T_prev, Tn = T, T
        for _ in range(60):                                              # the same backward-Euler step, lumped
            Tn = Tn - (mat.rho * mat._cp * V * (Tn - T_prev) / dt + 0.8 * SIGMA_SB * Tn ** 4 * A) / (
                mat.rho * mat._cp * V / dt + 4 * 0.8 * SIGMA_SB * Tn ** 3 * A)
        T = Tn
    assert s.T[s.active].mean() == pytest.approx(T, rel=1e-6)


def test_mass_rows_and_weights_on_a_grid_aligned_slab():
    m = build(BoxBase.around((-3e-3, -3e-3, -6e-3), (3e-3, 3e-3, 3e-3), 1e-3))
    z = m.points[:, 2]
    s = CutThermal(ConstMaterial()); s.setup(m, z)
    assert np.isfinite(s.w).all() and (s.Vin > 0).all()
    assert s.metal_volume() == pytest.approx(6e-3 * 6e-3 * 6e-3, rel=1e-12)


def test_ghost_penalty_keeps_slivers_conditioned_and_leaves_constants_alone():
    m = build(BoxBase.around((-14e-3,) * 3, (14e-3,) * 3, 1e-3))
    phi = np.linalg.norm(m.points, axis=1) - 10.001e-3
    ratios = {}
    for g in (0.0, 0.01):
        s = CutThermal(ConstMaterial(), ghost=g); s.setup(m, phi); s.set_temperature(300.0)
        K, M = s.operators(s.T, s.T)
        idx = np.flatnonzero(s.active)
        d = (M / 1e-3 + K).diagonal()[idx]
        ratios[g] = d.min() / d.max()
        assert np.allclose((K @ np.ones(len(s.T)))[idx], 0.0, atol=1e-12 * abs(K).max())
    assert ratios[0.0] < 1e-6 and ratios[0.01] > 3e-3


def test_move_books_the_enthalpy_it_removes_and_fills_new_nodes():
    m, phi = sphere()
    s = CutThermal(ConstMaterial()); s.setup(m, phi)
    s.set_temperature(300.0 + 1e4 * (m.points[:, 0] + 14e-3))
    E0 = s.energy()
    out = s.move(np.linalg.norm(m.points, axis=1) - 9.3e-3)            # shrink
    assert out == pytest.approx(E0 - s.energy(), rel=1e-14) and out > 0
    before = s.active.copy()
    s.move(np.linalg.norm(m.points, axis=1) - 10.4e-3)                   # grow: new nodes get neighbour values
    new = s.active & ~before
    assert new.any() and np.isfinite(s.T[new]).all()
    assert s.T[new].min() >= 300.0 and s.T[new].max() <= 300.0 + 1e4 * 28e-3


def test_receding_surface_through_slivers_gives_no_temperature_jumps():
    h = 1e-3
    m = build(BoxBase.around((-3.37e-3, -3.37e-3, -8e-3), (3.37e-3, 3.37e-3, 3e-3), h))
    z = m.points[:, 2]
    s = CutThermal(ConstMaterial(), linear_solver="direct"); s.setup(m, z); s.set_temperature(300.0)
    probe, Tp, top = np.array([[0.2e-3, 0.1e-3, -3.0e-3]]), [], 0.0
    for _ in range(60):
        s.step(0.02, 1e6 * np.maximum(s.P.normal[:, 2], 0))
        top -= 0.05 * h
        s.move(z - top)
        Tp.append(interpolate(m, s.T, probe)[0])
    d2 = np.abs(np.diff(Tp, 2))[3:]                                      # after the switch-on transient
    assert d2.max() < 3 * np.median(d2)


def test_remesh_closes_the_enthalpy_and_reports_the_residual():
    h = 1e-3
    base = BoxBase.around((-3.37e-3, -3.37e-3, -8e-3), (3.37e-3, 3.37e-3, 3e-3), h)
    m = build(base)
    z = m.points[:, 2]
    s = CutThermal(LinearMeltMaterial()); s.setup(m, z + 0.3e-3)
    s.set_temperature(894.0 + 2e4 * (z + 0.3e-3) + 3e7 * (z + 0.3e-3) ** 2)   # curved, crossing the melting range
    E0 = s.energy()
    m2 = build(base, lambda c, sz: np.where(np.abs(c[:, 2] + 0.3e-3) < 1.5e-3, 0.25e-3, np.inf))
    r0 = s.remesh(m2, interpolate(m, s.phi, m2.points))
    assert 0.0 < abs(r0) < 1e-2
    assert s.energy() == pytest.approx(E0, rel=1e-12)



def test_remesh_widens_the_offset_bracket_when_the_metal_volume_changes_a_lot():
    h = 1e-3
    base = BoxBase.around((-3.37e-3, -3.37e-3, -8e-3), (3.37e-3, 3.37e-3, 3e-3), h)
    m = build(base)
    z = m.points[:, 2]
    s = CutThermal(LinearMeltMaterial()); s.setup(m, z + 0.3e-3)
    s.set_temperature(894.0 + 2e4 * (z + 0.3e-3) + 3e7 * (z + 0.3e-3) ** 2)
    E0 = s.energy()
    m2 = build(base, lambda c, sz: np.where(np.abs(c[:, 2] + 0.3e-3) < 1.5e-3, 0.25e-3, np.inf))
    # the new surface is 0.6 mm higher, so the projected field holds far more enthalpy and the offset is large
    T_proj = project(m, m2, s._extend(s.active.copy(), ~s.active, rings=3, strict=False))
    s.remesh(m2, m2.points[:, 2] + 0.9e-3)
    assert abs((s.T - T_proj)[s.active].mean()) > 50.0                     # beyond the first bracket (about 71 K)
    assert s.energy() == pytest.approx(E0, rel=1e-12)


def test_psi_tracks_a_one_phase_stefan_front_within_one_fine_cell():
    """Design section 20.11 item 2: a 4.2 mm square column, 22 mm deep, on a uniform 0.7 mm mesh (the fine cell),
    starts solid at the bottom of an 890-900 K melting range; its top, 0.3 mm above a node layer, is held near
    T_hot = T_m + 100 K by the implicit base law with H = 1e7 W/m2/K; its sides and bottom are the insulated mesh
    boundary. The depth of psi = 0 (the 895.1 K isotherm) on the centre line must stay within one cell of the
    Neumann front s(t) = 2 lambda sqrt(alpha t)."""
    from scipy.optimize import brentq
    h, dt, z_top, H = 0.7e-3, 0.01, 0.3e-3, 1e7
    mat = LinearMeltMaterial(L=390e3, T_s=890.0, T_l=900.0)
    T_m = 0.5 * (mat.T_s + mat.T_l)
    T_hot = T_m + 100.0
    m = build(BoxBase.around((-3 * h, -3 * h, -31 * h), (3 * h, 3 * h, h), h))
    phi = m.points[:, 2] - z_top
    alpha = mat._k / (mat.rho * mat._cp)
    St = mat._cp * (T_hot - T_m) / mat.L
    lam = brentq(lambda x: x * exp(x * x) * erf(x) - St / sqrt(pi), 1e-6, 5.0)
    s = CutThermal(mat, linear_solver="direct"); s.setup(m, phi); s.set_temperature(mat.T_s)
    zs = z_top - np.linspace(0.0, 20e-3, 4001)
    line = np.column_stack([np.zeros_like(zs), np.zeros_like(zs), zs])
    errs = {}
    for n in range(1, 201):
        s.step(dt, H * T_hot, q1=-H)
        if n in (25, 50, 100, 200):
            t = n * dt
            psi, _ = molten_rigid_surface(m, phi, s.T)
            p = interpolate(m, psi, line)
            i = np.flatnonzero((p[:-1] > 0) & (p[1:] <= 0))[0]
            depth = (z_top - zs[i]) + (zs[i] - zs[i + 1]) * p[i] / (p[i] - p[i + 1])
            errs[t] = depth - 2 * lam * sqrt(alpha * t)
    msg = ", ".join("t = {:.2f} s: {:+.3f} mm".format(t, e * 1e3) for t, e in errs.items())
    assert all(abs(e) < h for e in errs.values()), msg
