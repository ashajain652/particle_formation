"""Builders the level-set tests share (Task 9b, fix round 1: moved here from the test modules that imported them from
one another): the gas state of the loads and body tests, markers seeded as the body seeds them, a box mesh refined along a
level set, and each marker's tessellated area against the area nearest it, with the acceptance numbers of the ruling on
the neighbour test; and (plan 03b Task 2, moved here from test_levelset_division.py) the exact signed distances of a box
and of a polygonal prism, and the stage 3 box refined by the body's refinement rule around its exact distance, with (plan
03b Task 4) its division at the skin depth. Test-only; no physics of the level set lives here."""
import functools
from collections import namedtuple

import numpy as np
from scipy.spatial import cKDTree

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.criteria import kappa_max_nodal, make_target
from reentry_model.levelset.division import divide
from reentry_model.levelset.cutcell import inside_volumes, patches
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


def box_sdf(lo, hi):
    """The exact signed distance to the box from `lo` to `hi` [m], negative inside."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    c, half = 0.5 * (lo + hi), 0.5 * (hi - lo)

    def sdf(p):
        q = np.abs(p - c) - half
        return np.linalg.norm(np.maximum(q, 0.0), axis=1) + np.minimum(q.max(axis=1), 0.0)
    return sdf


def polygon_sdf(corners):
    """The exact signed distance in (x, y) to the polygon with these corners [m] (any order round it), negative inside,
    the same at every z: a prism through the mesh."""
    a = np.asarray(corners, float)
    b = np.roll(a, -1, axis=0)

    def sdf(p):
        q = p[:, None, :2]
        ab = b - a
        t = np.clip(np.einsum("nek,ek->ne", q - a, ab) / np.einsum("ek,ek->e", ab, ab), 0.0, 1.0)
        d = np.linalg.norm(q - (a + t[..., None] * ab), axis=2).min(axis=1)
        crosses = ((a[:, 1] > q[..., 1]) != (b[:, 1] > q[..., 1])) & \
            (q[..., 0] < a[:, 0] + (q[..., 1] - a[:, 1]) * ab[:, 0] / np.where(ab[:, 1] == 0.0, 1.0, ab[:, 1]))
        inside = crosses.sum(axis=1) % 2 == 1
        return np.where(inside, -d, d)
    return sdf


BODY_FINE, BODY_FLOOR = 0.7e-3, 0.25e-3    # m: the body's fine cells and cell floor (LevelSetSkinBody's defaults)
BODY_MARGIN, BODY_LEAD = 4e-3, 0.6e-3      # m: its fine margin and the fine region's lead (stepend.LEAD)
BODY_REFINES = 2                           # its refinements from the fields at construction (skinbody.INITIAL_REFINES)
BOX_H0 = 2.8e-3                            # m: the stage 3 box's base cells (analysis/levelset/stage3_cases.py H0)
BOX_LENGTH, BOX_HALF, BOX_OFF = 8e-3, 5e-3, 0.13e-3     # m: the stage 3 box and its offset off the grid (stage3_cases)

Body = namedtuple("Body", "mesh phi x normal base sdf")


def body_mesh(base, sdf, fine=BODY_FINE, floor=BODY_FLOOR):
    """`base` refined by LevelSetSkinBody's refinement rule at construction around a cold body (all rigid, so psi =
    phi), applied to the exact signed distance `sdf`: BODY_REFINES passes of criteria.make_target with the body's fine
    size, margin, lead and cell floor, the curvature of phi recovered within two base cells of it.

    This is not a body's own mesh (plan 03b Task 2, fix round 1, measured once). A body refines around its own phi, and
    rebuilds its mesh from its fields at construction. Today's body builds phi from beta moved out by the 0.8 mm stub,
    so its edges are rounded to the stub. On the stage 3 box, `box_body("dev")` has 38,324 nodes and 216,000 cells,
    against 65,296 nodes and 371,776 cells here, and the two share 34,407 nodes. Here the sharp edges' curvature puts
    203,600 cells at the 0.175 mm floor, against the body's 62,672; within 2 mm of the surface the body's phi differs
    from the box's exact distance by up to 0.69 mm."""
    mesh = build(base)
    for _ in range(BODY_REFINES):
        phi = sdf(mesh.points)
        kp = kappa_max_nodal(mesh, phi, band=2.0 * base.h0)
        mesh = build(base, make_target(mesh, phi, phi, kp, kp, fine=fine, margin=BODY_MARGIN, lead=BODY_LEAD,
                                       floor=floor))
    return mesh


def body_on(base, sdf):
    """A body's level set: the exact signed distance `sdf` on `base` refined by the body's refinement rule
    (`body_mesh`), with markers seeded on it as the body seeds them. Returns a Body (mesh, phi, marker positions, marker normals, base,
    sdf)."""
    mesh = body_mesh(base, sdf)
    phi = sdf(mesh.points)
    mk = seeded(mesh, phi, DESIGN_SPACING, base.lo)
    return Body(mesh, phi, mk.x, mk.normal, base, sdf)


def box_corners():
    """The stage 3 box's lower and upper corners [m]: 8 x 10 x 10 mm, 0.13 mm off the grid."""
    lo = np.array([BOX_OFF, BOX_OFF - BOX_HALF, BOX_OFF - BOX_HALF])
    return lo, lo + np.array([BOX_LENGTH, 2.0 * BOX_HALF, 2.0 * BOX_HALF])


@functools.lru_cache(maxsize=1)
def stage3_box():
    """The stage 3 box (analysis/levelset/stage3_cases.py's `box_body`: its BoxBase, its 0.7 mm fine cells and 0.25 mm
    cell floor) as a Body, phi its exact signed distance, on the mesh the body's refinement rule gives that distance
    (`body_on`; not `box_body`'s own mesh, `body_mesh`)."""
    lo, hi = box_corners()
    base = BoxBase.around((-BOX_H0, -BOX_HALF - BOX_H0, -BOX_HALF - BOX_H0),
                          (BOX_LENGTH + BOX_H0, BOX_HALF + BOX_H0, BOX_HALF + BOX_H0), BOX_H0)
    return body_on(base, box_sdf(lo, hi))


@functools.lru_cache(maxsize=1)
def stage3_division():
    """The stage 3 box's layer divided among its seeded markers at the skin depth (`division.divide` at
    SkinSettings().depth with its default levels), in the markers' order of `stage3_box`. Measured once (plan 03b Task 4):
    132 markers, 405 shared boundaries (108 between columns whose normals differ by more than 60 degrees, 20 with no
    length at phi but some at depth), one column pinching before H, divided in 8.3 s."""
    box = stage3_box()
    return divide(box.mesh, box.phi, SkinSettings().depth, box.x, box.normal)


def edge_length(size):
    """The total length [m] of a box's twelve edges, its sides `size` (3,) [m]."""
    return 4.0 * float(np.sum(size))


def rounding_volume(d, length):
    """What rounding an edge of length `length` [m] to radius |d| [m] changes [m3], (1 - pi/4) d^2 per unit length (plan
    03b, planning ruling 11)."""
    return (1.0 - np.pi / 4.0) * d * d * length


def volume_between(mesh, a, b):
    """The volume [m3] between the zero contours of the level sets a and b on the mesh, inside one and not the other:
    inside_volumes of their nodal minimum less that of their nodal maximum (the union less the intersection, to within
    the cells both contours cut). Planning ruling 11 judges sharpness by it."""
    return float(inside_volumes(mesh, np.minimum(a, b)).sum() - inside_volumes(mesh, np.maximum(a, b)).sum())
