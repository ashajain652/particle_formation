"""analysis/skin_base_cells/column.py, the single-column test of a thick skin's base (open item 7, 2026-11-14): a 28 mm
pool at 906 K, one real skin column, on rigid metal at 894 K, against the one-dimensional reference on one grid."""
import pytest

from analysis.skin_base_cells import column


@pytest.fixture(scope="module")
def reference_front():
    """The reference's rigid front at 0.5 s on 10 um cells and 1 ms steps (3.010 mm; 3.005 mm on 5 um, 1 ms and 0.5 ms)."""
    rows = column.run_reference(column.options("reference", dx_um=10.0, dt_ms=1.0, t_end=0.5, every=0.5), log=lambda s: None)
    return rows[-1]["front_mm"]


def test_with_the_default_base_cell_a_thick_pool_s_rigid_front_follows_the_reference(reference_front):
    """With the default base cell (0.2 mm) the skin's rigid front at 0.5 s is within 2 % of the reference (1.2 % when
    measured, against 5 um), the step from the bottom cell to the base under 1 K (0.11 K), and the books close."""
    rows, closure = column.run_skins(column.options("skins", t_end=0.5, every=0.5), log=lambda s: None)
    assert 2.9 < reference_front < 3.1
    assert rows[-1]["front_mm"] == pytest.approx(reference_front, rel=0.02)
    assert rows[-1]["step_K"] < 1.0 and rows[-1]["dz_bottom_mm"] < 0.3
    assert abs(closure["mass_rel"]) < 1e-12 and abs(closure["energy_rel"]) < 1e-12 and closure["steps_unconverged"] == 0


def test_graded_from_the_top_alone_a_thick_pool_s_rigid_front_runs_far_ahead(reference_front):
    """The grading before the base cell (`--skin-base-cell 0`): the bottom cell is millimetres thick and the rigid front at
    0.5 s more than 30 % too far (64 % when measured), the step 2.2 K -- the jump-fault probe's thick skins."""
    rows, closure = column.run_skins(column.options("skins", t_end=0.5, every=0.5, base_cell_mm=0.0), log=lambda s: None)
    assert rows[-1]["front_mm"] > 1.3 * reference_front
    assert rows[-1]["dz_bottom_mm"] > 2.0 and rows[-1]["step_K"] > 1.5
    assert abs(closure["mass_rel"]) < 1e-12 and abs(closure["energy_rel"]) < 1e-12
