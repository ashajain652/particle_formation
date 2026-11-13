"""Builders the level-set tests share (Task 9b, fix round 1: moved here from the test modules that imported them from
one another): the gas state of the loads and body tests, markers seeded as the body seeds them, a box mesh refined along a
level set, and each marker's tessellated area against the area nearest it, with the acceptance numbers of the ruling on
the neighbour test. Test-only; no physics of the level set lives here."""
import numpy as np
from scipy.spatial import cKDTree

from reentry_model.levelset.boxmesh import build
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.markers import Markers
from reentry_model.levelset.repartition import surface_samples
from reentry_model.levelset.skin import SKIN_PRESETS, SkinSettings
from reentry_model.shared import atmosphere, trajectory as tj

DESIGN_SPACING = 2e-3      # m: the body's marker spacing (skinbody.MARKER_SPACING)
DESIGN_FINE = 0.7e-3       # m: the body's fine cells (its default fine_size)
CAPPED_SHARE = 0.01        # of the cells: at most this many bounded by the cap on a continuous surface (Task 9b ruling)
OWNERSHIP_FACTOR = 1.5     # a static tessellation's cell's area within this factor of the area nearest the marker
DEV_SETTINGS = SkinSettings(cell=SKIN_PRESETS["dev"])   # the development preset's skins (20 um cells)


def gas_state(rho=2.727e-5, T=203.6, V=7500.0, kn=0.0298, ma=26.2):
    """An aero state at 77.5 km and 7.5 km/s (the freestream's density, temperature, speed, Knudsen and Mach numbers
    adjustable), the deceleration a thousandth of the speed per second, for the gas loads and the body's steps."""
    fs = atmosphere.Freestream(rho, T, rho / 28.9644e-3 * 8.314462618 * T, 28.9644 * 1.66053906660e-27, np.zeros(3))
    return tj.AeroState(77.5e3, 0.0, 0.0, fs, np.array([V, 0.0, 0.0]), V, kn, ma, 0.9, np.array([-V * 1e-3, 0, 0]),
                        0.5 * rho * V * V)


def seeded(mesh, phi, spacing, lo, material=None, settings=DEV_SETTINGS):
    """Markers seeded as the body seeds them (skinbody.LevelSetSkinBody): one per voxel of the target spacing [m] from
    the corner `lo`, on the voxel's largest patch; with a `material`, each with an empty skin."""
    P = patches(mesh, phi)
    key = np.floor((P.centroid - np.asarray(lo)) / spacing).astype(np.int64)
    order = np.lexsort((-P.area, key[:, 2], key[:, 1], key[:, 0]))
    k = key[order]
    first = np.r_[True, np.any(k[1:] != k[:-1], axis=1)]
    return Markers.create(mesh, phi, np.sort(order[first]), material=material, settings=settings, spacing=spacing)


def nearest_ratio(mk, mesh, phi):
    """Each marker's tangent-plane cell's area (`Markers.cell_area`, the cells its faces bound; its area itself is the
    area of phi it owns since Task 11's fix round 1, C1) over the area of phi nearer to it than to any other marker
    (surface samples a sixteenth of the spacing across), the two scaled to the same total: how well the tessellation
    found the true neighbours."""
    y, w = surface_samples(mesh, phi, mk.spacing / 16.0)
    _, j = cKDTree(mk.x).query(y)
    near = np.bincount(j, w, minlength=mk.n)
    return mk.cell_area / near * near.sum() / mk.cell_area.sum()


def refined(base, level, fine):
    """The box mesh of `base` refined to `fine` [m] where |level(x)| is within two fine cells (and the leaf's own size)
    of zero."""
    return build(base, lambda x, s: np.where(np.abs(level(x)) < 2.0 * fine + 0.87 * s, fine, np.inf))
