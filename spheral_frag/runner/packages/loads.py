"""SurfaceLoads and Deceleration: the aerodynamic loads on the particles and the frame's deceleration (Spheral M2 plan,
Task 4 Steps 3-4; spec §7.1-7.2).

Conventions (spheral_frag.contract.V_HAT, spheral_frag.loads). The body frame's flight direction is +x: the nose is
the +x end, the oncoming flow runs toward -x, "upstream" is larger x. A particle's inclination theta is the angle
between its outward normal and +x (0 at the stagnation point, 180 deg at the tail), the angle the prepared load
tables are binned by. The pressure pushes along -n; the shear tau acts along t, the oncoming flow's direction -x
projected onto the surface and normalised, whose component along -x is sin(theta) (as `loads.drag`: the drag is
sum A (p cos theta + tau sin theta) along -x). In RZ, Spheral's positions are (x, r) with x the flight axis; a ring
particle's outward normal lies in its meridional plane, (n_x, n_r), and its force is the whole ring's: area
2 pi r ds, mass the ring's (plan fact 2), so the ring's radial acceleration is p A n_r / m and the net radial force
of the body is zero by symmetry (only the axial sums enter the deceleration).

**Surface data, once per interval** (at construction, from the current positions; spec §7.2's "each step, or every few
steps", and positions move far less than a spacing within an interval, which the driver records with
`max_displacement`):

1. *Surface particles and normals* by the colour function's gradient, G_i = h sum_j V_j grad_i W(x_i - x_j, h), over
   the neighbours within the kernel's support 2h (scipy's cKDTree over every particle, gathered to every rank;
   plan fact 1), with W the Wendland C4 shape of Spheral's kernel normalised to unit integral, h = nperh dx (the
   body's round smoothing length) and V_j = m_j / rho_j (unscaled: both scale with s). G vanishes in the interior
   of a uniform lattice and rises to 0.75-0.85 in a lattice sphere's outermost row, smoothly over the kernel's
   width (0.43-0.81 between 0.5 and 1 spacing deep, 0.19-0.68 between 1 and 1.5, at most 0.27 deeper). A particle
   is on the surface where |G| > `threshold`; its outward normal is -G/|G| with G evaluated again at the wider
   smoothing length `normal_smoothing` h. Calibrated on the lattice spheres of radius 50 mm (3D at 2.2 and 5 mm,
   RZ at 0.55, 1.1, 2.2 and 5 mm; 2026-10-09): `THRESHOLD` 0.5 marks every particle within 0.5 dx of the sphere
   and none deeper than 1.07 dx, in RZ exactly the particles with an empty lattice neighbour (63 at 2.2 mm, 127 at
   1.1 mm), in 3D those (5,258 at 2.2 mm) and 936 more within 1.07 dx (the staircase's inner corners); the
   bracket 0.3-0.6 keeps both properties (deepest 1.40 dx). Normals at h itself lie up to 12 deg (3D) and 15-18 deg
   (RZ) off radial, because the colour gradient follows the lattice's steps; at `NORMAL_SMOOTHING` = 2 h within
   2.4 deg (3D 2.2 mm), 3.0 (3D 5 mm), 4.8, 5.7 and 7.1 deg (RZ 2.2, 1.1, 0.55 mm). RZ works in the meridional half-plane with V_j the planar area m_j / (rho_j 2 pi r_j) and every particle
   mirrored to -r, so the axis is interior, not a surface.
2. *Exposure.* Lines parallel to the flight axis in bins one spacing wide across it, centred on the body stage's
   lattice lines (3D: round(y/dx), round(z/dx); RZ: floor(r/dx)). `exposure="first_body"` (default): along each line
   from upstream, the particles up to the first gap wider than `gap` spacings (`GAP`, 1.5) are the first body the
   line meets; its surface particles, windward and lee alike, take the tables' loads at their inclination, and the
   surface particles behind the gap (in a pit's shadow, in the wake of the body or of a separated piece) receive the
   base pressure only, no shear. `exposure="most_upstream"` (the plan's literal rule): only the most upstream
   surface particle per line is exposed, every other surface particle gets the base pressure. See "Deviation".
3. *Area shares.* `area_shares="volume"` (as built in Task 4; this package's keyword default, since the patch
   shares need the frame's patches, but no longer the runner's): a_i = V_i^(2/3) in 3D; in RZ the ring area
   2 pi r_i ds_i with ds_i = (planar area)^(1/2); scaled by one factor so that they sum, over every surface
   particle, to `surface_area`: the prepared frame's surface area (sum of its patch areas, `frame_surface_area`) at
   the interval's start frame k. `area_shares="patch"` (Task 8, check 6; the runner's default since decision 9,
   2026-10-11, `run.py --area-shares`): every patch of frame k's finite-element
   surface (`frame_patches`: centroid and area) gives its area to the surface particle whose *reference* position is
   nearest its centroid (the frame and the reference positions share the body coordinates, plan decision D1; in RZ
   the centroid's (x, sqrt(y^2 + z^2)) against the ring's (x, r), so a ring collects every patch around the axis
   nearest it), so the shares sum to the frame's area without a factor (`share_factor` 1) and follow the real
   surface's distribution; a surface particle no patch is nearest to (3D: the staircase's inner corners, a little
   deeper than their neighbours) gets no area and no load (`n_zero_area`). *Why* (check 6, 2026-10-10): the volume
   shares are equal per surface particle, but a lattice's staircase surface carries a number of surface particles
   per unit area that depends on the inclination (in RZ about max(|n_x|, |n_r|) / dx per unit length: the 45-degree
   steps have fewer), so after the one calibration factor the volume shares overweight the surface near the axis
   and the equator by 10-11 % and underweight it near 45 and 135 degrees by 11-13 %, at every spacing (it is a
   property of the lattice, not a resolution error): RZ drag -3 % against the table's integral (3D -0.2 to -0.6 %),
   and the particles' stresses in the core (|x| <= R/2) off the finite-element reference's by up to 21 % (RZ) and
   7-13 % (3D) of its peak, against 6-7 % and 4-8 % with the patch shares, whose drag is within 0.6 % (check 6,
   analysis/spheral_m2_check6.py, 2026-10-11). The patch shares leave a small net moment in 3D (4-6e-5 of
   drag x R: the unstructured patches against the lattice), which turns the body slowly against the damping.

**Loads, each step.** The tables of frames k and k + 1 (the prepared windward tables with the lee extension of
`spheral_frag.loads.with_lee` at the run's brackets) are evaluated at each exposed particle's theta once, at
construction; each step blends the two forces linearly in time, F_i(t) = (1 - w) F_i,k + w F_i,k+1 with
F = -p A n + tau A t, and adds F_i / m_i (the scaled mass, Review focus 1) to the hydro acceleration. The base
pressure is the lee extension's `base_fraction * p_stag` of each frame's table. A frame without loads (frame 0 of
a flight, written before any step: every bin NaN and p_stag NaN) takes the other frame's table of the interval.

**Deceleration** (spec §7.1): a_frame = sum F / sum m over every particle (M2: all are the main body), and every
particle receives -a_frame, so that sum m a = 0 and the body stays at rest in its own frame (drag along -x gives
a_frame_x < 0, i.e. every particle is pushed toward the nose). The sum of the forces is not gathered per step: the
loads package's forces are linear in w between two sums computed (over every rank) at construction, so
`SurfaceLoads.total_force(t)` is exactly the sum of the forces it adds at that t, whatever the order in which
Spheral calls the two packages, and costs no communication. In RZ only the axial component is a net force. Each
accepted step records (t_mid, sum F, a_frame) at the midpoint, where the step's forces are evaluated; the physical
deceleration of the flight is -s a_frame_x (the masses are scaled by s).

**Deviation from the plan's literal exposure rule** ("the most upstream particle per bin exposed, the rest receive
the base pressure"). Taken literally it gives the base pressure to the whole lee hemisphere (the tables' own lee
values, 19 Pa at frame 182 against a base of 3 % of p_stag = 231 Pa, are then never used) and to every surface
particle but one on the steep lines near the shoulder, where a line meets several surface particles at
theta near 90 deg; the spec's purpose ("on the undeformed body Spheral sees the finite-element model's loads
exactly", §7.2) and its examples of hidden particles (pits, the wake, a separated piece in the wake) are met by
"first body" instead. Both are settings; the drags they give are measured in the test's module docstring.

Numpy, scipy and Spheral only; `spheral_frag.loads` and `spheral_frag.frames` are the core's numpy-only modules."""
from __future__ import annotations

import math
import os

import numpy as np

from spheral_frag import loads as core_loads
from spheral_frag.runner.packages import (LocalCount, add_acceleration, allmax, allsum, comm, get_vector)

THRESHOLD = 0.5          # surface where |h grad c| exceeds this (calibrated on the lattice spheres; module docstring)
NORMAL_SMOOTHING = 2.0   # the normals' smoothing length in units of h (calibrated likewise)
GAP = 1.5                # spacings: a wider gap along a flight line ends the first body
EXPOSURES = ("first_body", "most_upstream")
AREA_SHARES = ("volume", "patch")
SUPPORT = 2.0            # the kernel's support in units of h (Spheral's WendlandC4 kernel extent)
TIME_TOL = 1e-9          # s

_CLASSES = {}


# ---------------------------------------------------------------------------------------------------------------
# The kernel and the colour gradient
# ---------------------------------------------------------------------------------------------------------------

def wendland(q):
    """Wendland C4 shape with support 2: (1 - q/2)^6 (1 + 3 q + 35/12 q^2), zero beyond."""
    q = np.asarray(q, dtype=float)
    a = np.clip(1.0 - 0.5*q, 0.0, None)
    return a**6*(1.0 + 3.0*q + 35.0/12.0*q*q)


def wendland_dq(q):
    """dW/dq of `wendland`: -(14/3) q (1 + 5 q / 2) (1 - q/2)^5."""
    q = np.asarray(q, dtype=float)
    a = np.clip(1.0 - 0.5*q, 0.0, None)
    return -(14.0/3.0)*q*(1.0 + 2.5*q)*a**5


def kernel_norm(dim, n=200001):
    """sigma with sigma integral of wendland(|x|) d^dim x = 1 (h = 1), by the trapezoid rule on q in [0, 2]."""
    q = np.linspace(0.0, SUPPORT, n)
    shell = {3: 4.0*math.pi*q*q, 2: 2.0*math.pi*q, 1: 2.0*np.ones_like(q)}[dim]
    f = wendland(q)*shell
    return 1.0/float(np.sum(0.5*(f[1:] + f[:-1])*np.diff(q)))


def colour_gradient(points, volumes, query, h, chunk=4000):
    """h grad c at each query point: h sum_j V_j grad W(x_q - x_j, h) over `points` within 2h (dimensionless). The
    query points are usually among the points (the self term is zero)."""
    from scipy.spatial import cKDTree
    points = np.asarray(points, dtype=float)
    query = np.asarray(query, dtype=float)
    dim = points.shape[1]
    sigma = kernel_norm(dim)
    tree = cKDTree(points)
    out = np.zeros((len(query), dim))
    for lo in range(0, len(query), chunk):
        q = query[lo:lo + chunk]
        nb = tree.query_ball_point(q, SUPPORT*h)
        lens = np.fromiter((len(a) for a in nb), dtype=np.int64, count=len(nb))
        if lens.sum() == 0:
            continue
        j = np.fromiter((k for a in nb for k in a), dtype=np.int64, count=int(lens.sum()))
        i = np.repeat(np.arange(len(q)), lens)
        d = q[i] - points[j]
        r = np.sqrt(np.einsum("ij,ij->i", d, d))
        with np.errstate(invalid="ignore", divide="ignore"):
            unit = np.where(r[:, None] > 0.0, d/np.where(r > 0.0, r, 1.0)[:, None], 0.0)
        # grad W = sigma / h^dim * dW/dq / h * unit;  times h: sigma / h^dim * dW/dq * unit
        g = (sigma/h**dim)*wendland_dq(r/h)[:, None]*unit*np.asarray(volumes, dtype=float)[j][:, None]
        for c in range(dim):
            out[lo:lo + chunk, c] = np.bincount(i, weights=g[:, c], minlength=len(q))
    return out


def detection_points(form, x, V):
    """The points and volumes the colour function is built from: 3D as given; RZ the planar areas
    V / (2 pi r) and every point mirrored to -r."""
    x = np.asarray(x, dtype=float)
    V = np.asarray(V, dtype=float)
    if form == "3d":
        return x, V
    Vp = V/(2.0*math.pi*x[:, 1])
    mirror = x*np.array([1.0, -1.0])
    return np.concatenate([x, mirror]), np.concatenate([Vp, Vp])


def detect_surface(form, x_all, V_all, query, h, threshold=THRESHOLD, normal_smoothing=NORMAL_SMOOTHING):
    """(|G|, surface flag, outward unit normal) at the query points (normal zero off the surface): the surface where
    |G| at h exceeds the threshold, the normal -G/|G| from the gradient at normal_smoothing * h."""
    pts, vol = detection_points(form, x_all, V_all)
    query = np.asarray(query, dtype=float)
    G = colour_gradient(pts, vol, query, h)
    g = np.sqrt(np.einsum("ij,ij->i", G, G))
    surf = g > threshold
    n = np.zeros_like(G)
    if surf.any():
        Gn = G[surf] if normal_smoothing == 1.0 else colour_gradient(pts, vol, query[surf], normal_smoothing*h)
        n[surf] = -Gn/np.sqrt(np.einsum("ij,ij->i", Gn, Gn))[:, None]
    return g, surf, n


# ---------------------------------------------------------------------------------------------------------------
# Exposure, inclination, area shares and forces (global arrays, identical on every rank)
# ---------------------------------------------------------------------------------------------------------------

def line_bins(form, x, dx):
    """Integer key of each particle's flight line: 3D (round(y/dx), round(z/dx)); RZ floor(r/dx)."""
    x = np.asarray(x, dtype=float)
    if form == "3d":
        iy = np.rint(x[:, 1]/dx).astype(np.int64)
        iz = np.rint(x[:, 2]/dx).astype(np.int64)
        return (iy - iy.min())*(int(iz.max() - iz.min()) + 1) + (iz - iz.min())
    return np.floor(x[:, 1]/dx).astype(np.int64)


def exposure(form, x, surface, dx, mode="first_body", gap=GAP):
    """Exposed surface particles (bool): see the module docstring. x and surface over every particle."""
    if mode not in EXPOSURES:
        raise ValueError(f"exposure must be one of {EXPOSURES}, got {mode!r}")
    x = np.asarray(x, dtype=float)
    surface = np.asarray(surface, dtype=bool)
    key = line_bins(form, x, dx)
    order = np.lexsort((-x[:, 0], key))           # by line, then from upstream (+x) down
    k = key[order]
    start = np.ones(len(order), dtype=bool)
    start[1:] = k[1:] != k[:-1]
    exposed = np.zeros(len(x), dtype=bool)
    if mode == "first_body":
        xs = x[order, 0]
        big = np.zeros(len(order), dtype=bool)
        big[1:] = (xs[:-1] - xs[1:] > gap*dx) & ~start[1:]
        # number of gaps passed so far within the line: cumulative sum restarted at each line's start
        cum = np.cumsum(big)
        first = np.maximum.accumulate(np.where(start, np.arange(len(order)), 0))
        run = cum - cum[first]
        exposed[order] = (run == 0) & surface[order]
        return exposed
    s = surface[order]
    idx = np.flatnonzero(s)
    if len(idx):
        ks = k[idx]
        lead = np.ones(len(idx), dtype=bool)
        lead[1:] = ks[1:] != ks[:-1]
        exposed[order[idx[lead]]] = True
    return exposed


def inclination_deg(normal):
    """theta (deg): the angle between the outward normal and the flight direction +x."""
    return np.degrees(np.arccos(np.clip(np.asarray(normal, dtype=float)[:, 0], -1.0, 1.0)))


def flow_tangent(normal):
    """t: the oncoming flow's direction -x projected onto the surface, normalised (zero where the normal is along
    the axis)."""
    n = np.asarray(normal, dtype=float)
    e = np.zeros_like(n)
    e[:, 0] = -1.0
    t = e - np.einsum("ij,ij->i", e, n)[:, None]*n
    norm = np.sqrt(np.einsum("ij,ij->i", t, t))
    ok = norm > 1e-12
    out = np.zeros_like(n)
    out[ok] = t[ok]/norm[ok][:, None]
    return out


def raw_shares(form, x, V):
    """a_i before calibration: V^(2/3) in 3D; 2 pi r (V / (2 pi r))^(1/2) in RZ."""
    V = np.asarray(V, dtype=float)
    if form == "3d":
        return V**(2.0/3.0)
    r = np.asarray(x, dtype=float)[:, 1]
    return 2.0*math.pi*r*np.sqrt(V/(2.0*math.pi*r))


def patch_shares(form, x_ref_surface, centroids, areas):
    """Area per surface particle (`area_shares="patch"`, module docstring): each patch's area to the surface particle
    whose reference position is nearest its centroid (RZ: the centroid's (x, sqrt(y^2 + z^2))). Ties go to the
    lower index (cKDTree). Arrays identical on every rank give identical shares."""
    from scipy.spatial import cKDTree
    xs = np.asarray(x_ref_surface, dtype=float)
    c = np.asarray(centroids, dtype=float)
    if form == "rz":
        c = np.stack([c[:, 0], np.hypot(c[:, 1], c[:, 2])], axis=1)
    if len(xs) == 0:
        raise ValueError("no surface particle to give the patches to")
    j = cKDTree(xs).query(c)[1]
    return np.bincount(j, weights=np.asarray(areas, dtype=float), minlength=len(xs))


def lee_table(table, lee):
    """The table with the lee extension at the run's brackets (`lee`: base_fraction, separation_deg,
    shear_factor; None for spheral_frag.loads' defaults)."""
    lee = dict(lee or {})
    return core_loads.with_lee(table, **lee)


def base_pressure(table):
    frac = (table.lee_params or {}).get("base_fraction", core_loads.BASE_FRACTION)
    return float(frac)*float(table.p_stag)


def particle_loads(table, theta_deg, exposed):
    """(p, tau) per surface particle: the table at theta where exposed, the base pressure and no shear elsewhere."""
    p, tau = core_loads.evaluate(table, theta_deg)
    p = np.where(exposed, p, base_pressure(table))
    tau = np.where(exposed, tau, 0.0)
    return p, tau


def surface_forces(normal, area, p, tau):
    """F = -p A n + tau A t per surface particle."""
    A = np.asarray(area, dtype=float)[:, None]
    return -np.asarray(p, dtype=float)[:, None]*A*normal + np.asarray(tau, dtype=float)[:, None]*A*flow_tangent(normal)


# ---------------------------------------------------------------------------------------------------------------
# Inputs from the prepared run
# ---------------------------------------------------------------------------------------------------------------

def read_tables(prepared_dir, ks):
    """The prepared windward load tables of frames `ks` (spheral_frag.loads.LoadTable, no lee extension) from
    `<prepared>/loads.npz`."""
    with np.load(os.path.join(prepared_dir, "loads.npz"), allow_pickle=False) as z:
        a = {k: z[k] for k in z.files}
    rows = []
    for k in ks:
        hit = np.flatnonzero(a["k"] == int(k))
        if len(hit) != 1:
            raise ValueError(f"frame {k} is not in {prepared_dir}/loads.npz")
        rows.append(int(hit[0]))
    sub = {key: (a[key][rows] if key not in ("theta_deg", "edges_deg") else a[key]) for key in
           ("theta_deg", "edges_deg", "p", "tau", "area", "lee", "theta_last_deg", "p_stag")}
    return core_loads.unstack(sub)


def frame_patches(prepared_dir, k):
    """(centroids (m, 3), areas (m,)) of the prepared frame k's surface patches (body coordinates, m and m^2)."""
    from spheral_frag import frames
    fr = frames.load_prepared_frame(os.path.join(prepared_dir, "frames", frames.prepared_frame_name(k)))
    return np.asarray(fr.patch["centroid"], dtype=float), np.asarray(fr.patch["area"], dtype=float)


def frame_surface_area(prepared_dir, k):
    """The prepared frame k's surface area (m^2): the sum of its patch areas (the faceted finite-element surface)."""
    from spheral_frag import frames
    fr = frames.load_prepared_frame(os.path.join(prepared_dir, "frames", frames.prepared_frame_name(k)))
    return float(np.sum(fr.patch["area"]))


# ---------------------------------------------------------------------------------------------------------------
# The packages
# ---------------------------------------------------------------------------------------------------------------

def _allgather_rows(*arrays):
    parts = comm().allgather(tuple(np.asarray(a) for a in arrays))
    return [np.concatenate([p[i] for p in parts], axis=0) for i in range(len(arrays))]


def loads_classes(S):
    if S in _CLASSES:
        return _CLASSES[S]

    class SurfaceLoads(S.Physics):
        """The aerodynamic loads on the surface particles (module docstring).

        SurfaceLoads(body, tables, times, surface_area, lee=None, threshold=THRESHOLD, exposure="first_body",
        gap=GAP, normal_smoothing=NORMAL_SMOOTHING, area_shares="volume", patches=None): `tables` the prepared
        windward tables of frames k and k + 1 (or one table, held), `times` their times (s), `surface_area` the
        frame's surface area (m^2), `lee` the lee extension's brackets, `area_shares` "volume" or "patch" (then
        `patches` = frame k's (centroids, areas), `frame_patches`).

        After construction: `n_surface`, `n_exposed` (global counts), `share_factor` (the area calibration),
        `force_k` (sum F at each frame, (2, dim)), `mass` (sum m), and per local particle `local_surface`,
        `local_normal`, `local_theta`, `local_area`, `local_exposed`; `surface_table()` gathers them by pid."""

        def __init__(self, body, tables, times, surface_area, lee=None, threshold=THRESHOLD,
                     exposure="first_body", gap=GAP, normal_smoothing=NORMAL_SMOOTHING, area_shares="volume",
                     patches=None):
            S.Physics.__init__(self)
            if exposure not in EXPOSURES:
                raise ValueError(f"exposure must be one of {EXPOSURES}, got {exposure!r}")
            if area_shares not in AREA_SHARES:
                raise ValueError(f"area_shares must be one of {AREA_SHARES}, got {area_shares!r}")
            if area_shares == "patch" and patches is None:
                raise ValueError("area_shares 'patch' needs the frame's patches")
            self.area_shares = area_shares
            tables = list(tables)
            times = [float(t) for t in times]
            if len(tables) == 1:
                tables, times = tables*2, (times*2 if len(times) == 1 else times)
            if len(tables) != 2 or len(times) != 2 or not times[1] >= times[0]:
                raise ValueError("need two tables and two non-decreasing times")
            tables = [lee_table(t, lee) for t in tables]
            if not tables[0].has_loads and not tables[1].has_loads:
                raise ValueError("neither frame's table carries loads")
            if not tables[0].has_loads:
                tables[0] = tables[1]
            if not tables[1].has_loads:
                tables[1] = tables[0]
            self.tables, self.times = tables, tuple(times)
            self.body, self.form, self.dim, self.dx = body, body.form, body.dim, body.dx
            self.count = LocalCount(body.nodes)
            self.threshold, self.exposure_mode, self.gap = float(threshold), exposure, float(gap)
            self.normal_smoothing = float(normal_smoothing)
            self.surface_area = float(surface_area)

            nodes = body.nodes
            x = get_vector(nodes.positions(), self.dim)
            from spheral_frag.runner.eos import get_internal
            m = np.array(get_internal(nodes.mass()))
            rho = np.array(get_internal(nodes.massDensity()))
            V = m/rho
            pid = body.local_pids()
            self.x0 = x
            self.local_mass = m
            x_ref = get_vector(body.x_ref, self.dim)
            x_all, V_all, pid_all, x_ref_all = _allgather_rows(x, V, pid, x_ref)
            h = body.nperh*body.dx
            g, surf, normal = detect_surface(self.form, x_all, V_all, x, h, self.threshold, self.normal_smoothing)
            surf_all, normal_all = _allgather_rows(surf, normal)
            # every rank holds the same global arrays in the same (rank, local) order: global decisions agree
            exposed_all = exposure_fn(self.form, x_all, surf_all, self.dx, self.exposure_mode, self.gap)
            if area_shares == "patch":
                raw_all = np.zeros(len(x_all))
                if not surf_all.any():
                    raise RuntimeError("SurfaceLoads: no surface particle found")
                raw_all[surf_all] = patch_shares(self.form, x_ref_all[surf_all], patches[0], patches[1])
                self.share_factor = 1.0
            else:
                raw_all = raw_shares(self.form, x_all, V_all)
                total_raw = float(np.sum(raw_all[surf_all]))
                if not total_raw > 0.0:
                    raise RuntimeError("SurfaceLoads: no surface particle found")
                self.share_factor = self.surface_area/total_raw
            self.n_zero_area = int(np.count_nonzero(surf_all & (raw_all == 0.0)))
            lo = int(np.sum(comm().allgather(len(x))[:comm().rank]))
            sl = slice(lo, lo + len(x))
            if not np.array_equal(pid_all[sl], pid):
                raise RuntimeError("SurfaceLoads: gathered rows out of order")
            self.local_grad = g
            self.local_surface = surf
            self.local_normal = normal
            self.local_exposed = exposed_all[sl]
            self.local_area = np.where(surf, raw_all[sl]*self.share_factor, 0.0)
            self.local_theta = np.where(surf, inclination_deg(normal), np.nan)
            self.n_surface = int(np.count_nonzero(surf_all))
            self.n_exposed = int(np.count_nonzero(exposed_all))
            self.idx = np.flatnonzero(surf)
            th = self.local_theta[self.idx]
            ex = self.local_exposed[self.idx]
            self.local_p, self.local_tau, F = [], [], []
            for tab in tables:
                p, tau = particle_loads(tab, th, ex)
                self.local_p.append(p)
                self.local_tau.append(tau)
                F.append(surface_forces(normal[self.idx], self.local_area[self.idx], p, tau))
            self.local_force = F
            self.acc = [f/m[self.idx][:, None] for f in F]
            fk = np.array([f.sum(axis=0) if len(f) else np.zeros(self.dim) for f in F])
            self.force_k = allsum(fk)
            if self.form == "rz":
                self.force_k[:, 1] = 0.0           # rings: no net radial force
            self.mass = allsum(m.sum())
            self.area_sum = allsum(self.local_area.sum())
            self.ramped = False

        def ramp_from_zero(self):
            """Make the first frame's loads zero, so that the loads rise linearly in time from zero at times[0] to
            the second frame's at times[1] (the replay's --loads-frame applies a held table this way over its first
            interval; run.py, "One frame's loads held"). total_force stays exactly the sum of what is added."""
            self.local_p[0] = np.zeros_like(self.local_p[0])
            self.local_tau[0] = np.zeros_like(self.local_tau[0])
            self.local_force[0] = np.zeros_like(self.local_force[0])
            self.acc[0] = np.zeros_like(self.acc[0])
            self.force_k[0] = 0.0
            self.ramped = True

        # -- time
        def weight(self, t):
            t0, t1 = self.times
            if t < t0 - TIME_TOL or t > t1 + TIME_TOL:
                raise ValueError(f"SurfaceLoads: t = {t!r} s lies outside the interval {t0!r} .. {t1!r} s")
            if t1 == t0 or t <= t0:
                return 0.0
            if t >= t1:
                return 1.0
            return (t - t0)/(t1 - t0)

        def total_force(self, t):
            """sum F (N) at time t over every particle: exactly what the package adds at t (no communication)."""
            w = self.weight(float(t))
            return (1.0 - w)*self.force_k[0] + w*self.force_k[1]

        def local_acceleration(self, t):
            w = self.weight(float(t))
            a = np.zeros((self.count.n, self.dim))
            a[self.idx] = (1.0 - w)*self.acc[0] + w*self.acc[1]
            return a

        def max_displacement(self):
            """max |x - x_construction| over every particle, in spacings (collective)."""
            d = get_vector(self.body.nodes.positions(), self.dim) - self.x0
            loc = float(np.sqrt(np.einsum("ij,ij->i", d, d)).max()) if len(d) else 0.0
            return allmax(loc)/self.dx

        def surface_table(self):
            """The surface data gathered to rank 0 sorted by pid (None elsewhere): pid, x, grad (|h grad c|),
            surface, exposed, normal, theta_deg, area, p_k, tau_k, p_k1, tau_k1, force_k, force_k1. Collective."""
            n = self.count.n
            dim = self.dim
            def full(vals, shape_tail=(), fill=np.nan):
                out = np.full((n,) + shape_tail, fill)
                out[self.idx] = vals
                return out
            local = dict(pid=self.body.local_pids(), x=self.x0, grad=self.local_grad, surface=self.local_surface,
                         exposed=self.local_exposed, normal=self.local_normal, theta_deg=self.local_theta,
                         area=self.local_area, p_k=full(self.local_p[0]), tau_k=full(self.local_tau[0]),
                         p_k1=full(self.local_p[1]), tau_k1=full(self.local_tau[1]),
                         force_k=full(self.local_force[0], (dim,), 0.0),
                         force_k1=full(self.local_force[1], (dim,), 0.0))
            parts = comm().gather(local, root=0)
            if comm().rank != 0:
                return None
            out = {k: np.concatenate([p[k] for p in parts], axis=0) for k in local}
            order = np.argsort(out["pid"], kind="stable")
            return {k: v[order] for k, v in out.items()}

        # -- Physics
        def evaluateDerivatives(self, t, dt, db, state, derivs):
            self.count.check("SurfaceLoads")
            add_acceleration(S, derivs, self.dim, self.local_acceleration(t))

        def dt(self, db, state, derivs, t):
            return (1.0e100, "SurfaceLoads: no limit")

        def registerState(self, db, state):
            return

        def registerDerivatives(self, db, derivs):
            return

        def label(self):
            return "SurfaceLoads"

    class Deceleration(S.Physics):
        """Every particle receives -a_frame, a_frame = sum F / sum m from `loads` (module docstring).

        Deceleration(loads, scale=1.0). `history()` -> (n_steps, 1 + 2 dim) array of accepted steps: t_mid,
        sum F (N), a_frame (m/s^2, Spheral's, i.e. with the scaled masses); `deceleration()` the last step's
        physical deceleration -scale a_frame_x (m/s^2)."""

        def __init__(self, loads, scale=1.0):
            S.Physics.__init__(self)
            self.loads = loads
            self.dim = loads.dim
            self.scale = float(scale)
            self.count = LocalCount(loads.body.nodes)
            self.last = None
            self.rows = []

        def frame_acceleration(self, t):
            F = self.loads.total_force(t)
            return F, F/self.loads.mass

        def evaluateDerivatives(self, t, dt, db, state, derivs):
            self.count.check("Deceleration")
            F, a = self.frame_acceleration(float(t))
            add_acceleration(S, derivs, self.dim, np.broadcast_to(-a, (self.count.n, self.dim)))
            self.last = (float(t), F, a)

        def finalize(self, t, dt, db, state, derivs):
            if self.last is not None:
                tm, F, a = self.last
                self.rows.append(np.concatenate([[tm], F, a]))

        def history(self):
            return np.array(self.rows).reshape(len(self.rows), 1 + 2*self.dim)

        def deceleration(self):
            if not self.rows:
                return float("nan")
            return -self.scale*float(self.rows[-1][1 + self.dim])

        def dt(self, db, state, derivs, t):
            return (1.0e100, "Deceleration: no limit")

        def registerState(self, db, state):
            return

        def registerDerivatives(self, db, derivs):
            return

        def label(self):
            return "Deceleration"

    _CLASSES[S] = (SurfaceLoads, Deceleration)
    return _CLASSES[S]


exposure_fn = exposure


def make_surface_loads(S, body, tables, times, surface_area, lee=None, threshold=THRESHOLD, exposure="first_body",
                       gap=GAP, normal_smoothing=NORMAL_SMOOTHING, area_shares="volume", patches=None):
    """A SurfaceLoads package for the Spheral module S (collective: every rank builds it)."""
    return loads_classes(S)[0](body, tables, times, surface_area, lee, threshold, exposure, gap, normal_smoothing,
                               area_shares, patches)


def make_deceleration(S, loads, scale=1.0):
    """A Deceleration package reading `loads` (a SurfaceLoads) for the Spheral module S."""
    return loads_classes(S)[1](loads, scale)
