"""The skin body (layered-skin design sections 6.6 and 7, sub-plan 20 Task 7): skins on every corner from time zero, the
trial and real passes around the 3D step, the base moving the 3D mesh, and the books."""
import numpy as np
import pytest

from reentry_model import heating, material, mesh, skin, skinbody
from reentry_model.thermal import thermal_solver

DEV = skin.SkinSettings(cell=20e-6)


def make(m, T0=880.0, name="skfem", settings=DEV, **kw):
    copy = mesh.VolumeMesh(m.points.copy(), m.tets.copy(), dict(m.params))
    mat = material.Material.from_drama_json("AA7075_scheil")
    return skinbody.SkinBody(copy, mat, thermal_solver(name), 1.47, settings, T0=T0, v_hat=(-1.0, 0.0, 0.0), **kw)


def front_loads(b, q0):
    return heating.HeatingResult(q0 * np.clip(np.cos(b.theta), 0.0, None), q0, 0.0, 0.0, 0.0)


def test_every_corner_starts_with_its_minimum_taken_from_the_3d_model(coarse_sphere_mesh):
    b = make(coarse_sphere_mesh)
    m_min = skin.minimum_mass(b.material, b.columns.area, DEV)
    V0 = coarse_sphere_mesh.volume()
    assert b.skins.n == b.columns.n == len(np.unique(coarse_sphere_mesh.surface().faces))
    assert b.mass() == pytest.approx(b.material.rho * V0, rel=1e-14)                  # nothing made or lost
    assert b.energy() == pytest.approx(b.material.rho * V0 * float(b.material.enthalpy(880.0)), rel=1e-12)
    # each column holds what its corner swept moving in by its minimum: within 2 % of it (prisms against dual areas
    # on 4 mm triangles), so the first sub-step draws any shortfall
    assert np.abs(b.skins.column_mass() / m_min - 1.0).max() < 0.02
    assert np.allclose(b.skins.T[b.skins.valid], 880.0)
    assert b.surface.n_patches == b.columns.n and b.theta[b.i_stag] < np.radians(5.0)


def test_a_body_with_no_heat_stays_uniform_and_its_mesh_and_skins_keep_their_books(coarse_sphere_mesh):
    b = make(coarse_sphere_mesh, T0=600.0)
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    b.emissivity = 0.0
    M0, E0 = b.mass(), b.energy()
    for k in range(3):
        b.advance(0.5 * (k + 1), 0.5, zero)
    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
    assert np.abs(b.solver.temperature() - 600.0).max() < 1e-6
    assert np.abs(b.skins.T[b.skins.valid] - 600.0).max() < 1e-6
    assert b.removed_mass == 0.0 and b.jump_faults == 0


def test_heating_with_removal_keeps_mass_and_energy_exact(coarse_sphere_mesh):
    """6 MW/m2 on the front for 1.5 s from 880 K: liquid forms and leaves, columns draw, the base recedes; the mass
    removed plus the mass left is the mass at start, and the energy books close to 1e-8 of the heat absorbed."""
    b = make(coarse_sphere_mesh)
    M0 = b.mass()
    for k in range(3):
        b.advance(0.5 * (k + 1), 0.5, front_loads(b, 6e6))
    assert b.removed_mass > 0.0 and b.drawn_mass > 0.0
    assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
    assert abs(b.energy_balance_residual()) < 1e-8                  # 1.6e-9 measured while planning
    stats = b.skin_stats()
    assert stats["base_recession_max_mm"] > 0.5                       # 0.88 mm at the nose (measured)
    assert b.surface_temperature().max() <= b.material.T_feed + 1.0   # liquid left as it formed
    assert stats["skin_cells"] == b.skins.count.sum() and b.last.Q_conv > 0.0


def test_hand_down_passes_rigid_bottom_metal_beyond_the_minimum():
    mat = material.Material.from_drama_json("AA7075_scheil")
    sk = skin.Skins(mat, np.full(3, 4e-6), DEV)
    m_min = skin.minimum_mass(mat, sk.area, DEV)
    sk.fill(1.5 * m_min, np.array([700.0, 700.0, 905.0]))            # two cold, one slurry (above T_rigid)
    sk.take_bottom(np.array([0.0, 0.6, 0.0]) * m_min)                 # the second below its minimum
    M0, E0 = sk.column_mass(), sk.H.sum(axis=1)
    handed, H = skin.hand_down(sk, m_min, mat.T_rigid)
    assert handed == pytest.approx(np.array([0.5 * m_min[0], 0.0, 0.0]), rel=1e-12)
    assert sk.column_mass() == pytest.approx(M0 - handed, rel=1e-13)
    assert sk.H.sum(axis=1) == pytest.approx(E0 - H, rel=1e-13)
    assert H[0] == pytest.approx(handed[0] * float(mat.enthalpy(700.0)), rel=1e-9)


def test_handing_down_moves_the_base_out_and_the_3d_model_takes_the_metal(coarse_sphere_mesh):
    """A cold skin above its minimum hands its bottom metal to the 3D model: the corners move out by its volume, the 3D
    model's volume grows by it, and mass and energy are exact."""
    b = make(coarse_sphere_mesh, T0=600.0)
    b.emissivity = 0.0
    extra = 0.5 * skin.minimum_mass(b.material, b.columns.area, DEV)
    b.skins.add_bottom(extra, float(b.material.enthalpy(600.0)))
    b.skins.regrid()
    M0, E0, V0 = b.mass(), b.energy(), b.mesh.volume()
    # every column's excess over its minimum, the added metal and what it held above its minimum from the start
    excess = np.maximum(b.skins.column_mass() - skin.minimum_mass(b.material, b.columns.area, DEV), 0.0).sum()
    b.advance(0.5, 0.5, heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0))
    assert b.handed_mass == pytest.approx(excess, rel=1e-9) and excess > extra.sum()
    assert b.mesh.volume() - V0 == pytest.approx(excess / b.material.rho, rel=0.02)
    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
    assert np.abs(b.solver.temperature() - 600.0).max() < 1e-6


# -- what triggers the interface repeat (Task 7, measured in Task 10) --------------------------------------------------
def test_every_repeat_rule_keeps_the_books_and_records_what_triggered_it(coarse_sphere_mesh):
    """Each repeat rule (measured in Task 10; the spec's per-column rule is the default) leaves the books exact, and the
    first pass's mismatch is recorded whatever the rule. 2 MW/m2 from 880 K on 4 mm elements: without a repeat, 3 MW/m2
    already runs away -- the 3D base overshoots T_feed and every sub-step removes and redraws a whole skin (measured
    while planning)."""
    results = {}
    for rule in ("column", "total", "never"):
        b = make(coarse_sphere_mesh, repeat_rule=rule)
        M0 = b.mass()
        for k in range(3):
            b.advance(0.5 * (k + 1), 0.5, front_loads(b, 2e6))
        assert b.removed_mass > 0.0
        assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
        assert abs(b.energy_balance_residual()) < 1e-8
        stats = b.skin_stats()
        assert stats["mismatch_total_ratio"] > 0.0 and stats["mismatch_worst_ratio"] >= stats["mismatch_total_ratio"]
        results[rule] = b.interface_repeats
    assert results["never"] == 0 and results["column"] >= results["total"] > 0
    with pytest.raises(ValueError):
        make(coarse_sphere_mesh, repeat_rule="sometimes")
