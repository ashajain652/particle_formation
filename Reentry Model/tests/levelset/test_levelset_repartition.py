"""Tests of the step end's re-partition (level-set stage 3 Task 9b, `reentry_model.levelset.repartition`; plan 03b Task 6:
across the division of the layer, by the columns' volumes): the transfer's rows and columns to round-off with no negative
entry and a diagonal transfer when nothing moves, and the column fit on its own, a column of zero volume included. The
skins following it are tested with the body's step end (test_levelset_skinbody.py)."""
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from scipy.spatial.transform import Rotation

from reentry_model.levelset import repartition
from reentry_model.levelset.boxmesh import BoxBase
from reentry_model.levelset.division import divide
from reentry_model.levelset.repartition import _fit_columns, surface_samples, transfer_matrix
from reentry_model.levelset.skin import SkinSettings

from cases import DESIGN_FINE, DESIGN_SPACING, refined, seeded

ROUND_OFF = 1e-14      # of a volume: the transfer's rows and columns; a column held to n eps of the mean volume at most,
                       # n the markers (the piece's summed rounding, which no flow moves: `repartition._fit_columns`)
H = SkinSettings().depth   # m: the skin depth the layer is divided to


def _sphere(R=10e-3):
    """Markers seeded one per 2 mm voxel on a 10 mm sphere on 0.7 mm cells, tessellated, and the layer divided among
    them at the skin depth."""
    base = BoxBase.around((-R - 5.6e-3,) * 3, (R + 5.6e-3,) * 3, 2.8e-3)
    mesh = refined(base, lambda x: np.linalg.norm(x, axis=1) - R, DESIGN_FINE)
    phi = np.linalg.norm(mesh.points, axis=1) - R
    mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
    mk.tessellate(mesh, phi)
    return mesh, phi, mk, divide(mesh, phi, H, mk.x, mk.normal)


def test_a_transfer_with_nothing_moved_is_the_old_volumes_on_its_diagonal():
    """The same markers, division and surface before and after: every sample keeps its owner, the stretch is one exactly
    and the transfer is diagonal, each entry the column's volume."""
    mesh, phi, mk, d = _sphere()
    samples = surface_samples(mesh, phi, DESIGN_SPACING / 16.0)
    T, stretch = transfer_matrix(mk.x, mk.normal, d.volume, samples, mk.x, d.volume, d, samples)
    off = T - csr_matrix(np.diag(T.diagonal()))
    print("nothing moved: {} markers, off-diagonal entries {}, stretch {:.17f} to {:.17f}".format(
        mk.n, off.count_nonzero(), stretch.min(), stretch.max()))
    assert off.count_nonzero() == 0
    np.testing.assert_allclose(T.diagonal(), d.volume, rtol=ROUND_OFF, atol=0.0)
    np.testing.assert_allclose(stretch, 1.0, rtol=ROUND_OFF, atol=0.0)


def test_a_transfer_after_the_markers_turn_holds_its_rows_and_columns_to_round_off():
    """The sphere's markers turned 0.3 rad about an axis through the centre, snapped, tessellated and the layer divided
    again on the same surface: the transfer's rows sum to the old columns' volumes times the stretch to 1e-14, its
    columns to the new columns' volumes to the piece's summed rounding (n eps of the mean volume); no entry is negative;
    and the stretch, on an unchanged surface, is one within the sampling's noise."""
    mesh, phi, mk, d0 = _sphere()
    samples = surface_samples(mesh, phi, DESIGN_SPACING / 16.0)
    x0, n0, V0 = mk.x.copy(), mk.normal.copy(), d0.volume.copy()
    rot = Rotation.from_rotvec([0.1, 0.2, 0.2]).as_matrix()
    mk.carry(lambda x: x @ rot.T - x, 1.0)
    mk.snap(mesh, phi)
    mk.tessellate(mesh, phi)
    d1 = divide(mesh, phi, H, mk.x, mk.normal)
    T, stretch = transfer_matrix(x0, n0, V0, samples, mk.x, d1.volume, d1, samples)
    rows, cols = np.asarray(T.sum(axis=1)).ravel(), np.asarray(T.sum(axis=0)).ravel()
    moved = (T.sum() - T.diagonal().sum()) / T.sum()
    print("turned: {:.1%} of the volume passed between columns; rows {:.1e}, columns {:.1e} off; stretch {:.4f} to "
          "{:.4f}".format(moved, np.abs(rows / (V0 * stretch) - 1).max(), np.abs(cols / d1.volume - 1).max(),
                          stretch.min(), stretch.max()))
    assert moved > 0.1
    np.testing.assert_allclose(rows, V0 * stretch, rtol=ROUND_OFF, atol=0.0)
    np.testing.assert_allclose(cols, d1.volume, rtol=0.0, atol=mk.n * np.finfo(float).eps * d1.volume.mean())
    assert T.data.min() >= 0.0
    assert stretch.sum() == pytest.approx(mk.n, rel=0.01) and np.abs(stretch - 1).max() < 0.1


def _chain(area, T0):
    """`_fit_columns` on markers in a row 1 m apart, each face 1 m long."""
    n = len(area)
    x = np.column_stack([np.arange(n, dtype=float), np.zeros(n), np.zeros(n)])
    fi, fj = np.arange(n - 1), np.arange(1, n)
    T0 = csr_matrix(T0)
    return _fit_columns(T0, np.asarray(T0.sum(axis=1)).ravel(), np.asarray(area, float), fi, fj, np.ones(n - 1), x,
                        np.zeros(n, np.int64))


def test_the_column_fit_meets_the_new_areas_keeps_the_rows_and_turns_nothing_negative():
    """A diagonal transfer of random rows on 12 markers in a row and random new areas of the same total: the fit's
    columns meet the new areas and its rows keep theirs, each to round-off, and no entry is negative."""
    rng = np.random.default_rng(4)
    R = rng.uniform(0.5, 2.0, 12)
    area = rng.uniform(0.2, 3.0, 12)
    area *= R.sum() / area.sum()
    T = _chain(area, np.diag(R))
    rows, cols = np.asarray(T.sum(axis=1)).ravel(), np.asarray(T.sum(axis=0)).ravel()
    np.testing.assert_allclose(rows, R, rtol=ROUND_OFF, atol=0.0)
    np.testing.assert_allclose(cols, area, rtol=0.0, atol=len(area) * np.finfo(float).eps * area.mean())
    assert T.data.min() >= 0.0


def test_the_column_fit_empties_a_column_whose_new_area_is_zero():
    """Fix round 1 (minor 1): six markers in a row with random rows, the third's new area zero (a marker the tessellation
    left no area), the area flowing through it: it ends empty, every column at its new area and every row at its old
    one, to round-off of the piece's mean area, no entry negative. (Before, the fit's tolerance was relative to each
    column's own area, zero here: rounding left 1e-17 in the column, and the fit raised, its message dividing by the
    zero area: "inf of an area left". The case is the first of 300 random ones that did.)"""
    rng = np.random.default_rng(197)
    R = rng.uniform(0.5, 2.0, 6)
    area = rng.uniform(0.2, 3.0, 6)
    area[2] = 0.0
    area *= R.sum() / area.sum()
    T = _chain(area, np.diag(R))
    cols = np.asarray(T.sum(axis=0)).ravel()
    print("a zero column: columns less the new areas, over their mean: {}".format((cols - area) / area.mean()))
    np.testing.assert_allclose(cols, area, rtol=0.0, atol=10 * repartition.FIT_RTOL * area.mean())
    np.testing.assert_allclose(np.asarray(T.sum(axis=1)).ravel(), R, rtol=ROUND_OFF, atol=0.0)
    assert T.data.min() >= 0.0
