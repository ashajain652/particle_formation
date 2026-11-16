"""Spraying per marker (design sections 9.5 and 20.6; level-set stage 3 Task 8): Girin's release from the skins that ride
on the markers, and its droplet records. The physics is the shared package's -- `spray.SprayModel` (Girin's thick and
thin branches, the Rayleigh-Taylor droplet, the gates), `spray.spray_layer` (Girin's Eq. (2) or the slip band's
momentum balance), `lubrication.lubrication` (the film's surface speed), `spray.slurry_properties` and
`spray.SLURRY_SIGMA` (the exposed slurry), `spray.region_extent` (the wave's domain), `spray.source_block` and
`spray.histogram` (the records), and the plan's architecture-neutral pieces both architectures call (`spray.slurry_modes`,
`plan_source_block`, `step_columns`, `droplet_mass`, `droplet_count`, `droplet_histogram`, `spraying`,
`plan_release_angle`, `no_slurry`, `no_step_spray`; moved there from the skyscrapers' skinbody.py in Task 8's fix round).
What this module holds is which depth of which skin goes into which argument, and the cascade that takes the release from
the skins' cells: it matches the skyscrapers' `girin_plan`, `cascade_mass` and `release_top` by amendment (decision D4:
the cascade lives per architecture, pinned by a test comparing the two on the same column); nothing here imports them.

The plan (`spray_plan`, once a macro step, at its start), per marker, in the markers' order:
- b, and the regime test's depth, the viscous film's floor and the lubrication's layer: the skin's flowing depth
  `depth_above(T_rigid)` (every skin is one skin depth deep, design 21.3, so b is always the skin's: plan 03b Task 5
  retired the top-zone skins' depth from phi to psi);
- delta_m from `spray_layer` (the slip band through it where the model's `slip_band` is on), v_s from
  `lubrication(tau, G, b, delta_m, mu_l, b_layer=b)`;
- the droplet's cap depth: the skin's liquid depth `depth_above(T_liquidus)`, or with `slurry` its flowing depth; the
  Rayleigh-Taylor pool the liquid depth; the film mass the run of releasable cells' mass; the wave's extent the
  contiguous region of markers with liquid on the graph of the division's shared boundaries at phi (those with a length
  at phi: `division.Division`'s i, j and length[:, 0]; plan 03b Task 6, the tangent-plane faces it read before are
  retired);
- with `slurry`, the exposed slurry's own mode per marker (`_slurry_modes`), a second `SprayModel` with the slurry's
  properties at the shear rate of its own sheared layer, found by fixed point.

No renewal limit (ruling of 2026-11-11): design 9.5's cap on liquid deeper than delta_m, rho_l delta_m A per renewal time
delta_m^2/nu_l, belongs to the two-zone skin of 2026-10-07 and is not built; the layered skins resolve the velocity
through the depth, so the release is capped only as the skyscrapers cap it, by the cap depth and the cascade.

The release (`SprayPlan.release(skins, dt)`, every flow sub-step): down each skin's top cells in time order, a liquid
cell (at or above the liquidus) at the liquid's rate and, with `slurry`, a flowing slurry cell at the slurry's own rate
while it is at least as warm as the cell that gave that rate, until dt is spent; never rigid metal (`cascade_mass`);
then exactly those cells' metal leaves at its own enthalpy and velocity (`release_top`, `Skins.take_top_mass`). The
plan books what left by the cells it came from (the liquid's mode's and the slurry's), and at the step's end gives the
records' rows (`records`) and the history's droplet columns (`step_summary`, every key the shared writer's
SKIN_SUMMARY_COLUMNS reads).

The profile path (plan 03b Task 10; design 21.7; the kinked-profile dispersion design,
docs/superpowers/specs/2026-11-11-kinked-profile-dispersion-design.md; a spray model whose `dispersion` is "profile"):
the skyscrapers' profile branch of `girin_plan`, matched by amendment. `layer_columns` reads each skin's cells as the
shared `melt_profile.LayerColumns` -- its own cells' speeds, under the gas shear, with its rigid front at
`depth_above(T_rigid)`, which is the column's own depth where no rigid metal lies in it, so that a skin flowing to its
base with no kink in it takes its base, at the skin's depth H, as its wall (the rigid front's kind: design 21.7's
"floor-like wall at depth H", the volume's melt held until the volume flow of 21.9 exists). The shared `melt_profile.profile_modes` fits the windward skins' profiles, places their walls and solves
their gradient mode for the shared spray model's profile path; each skin's release run is then fixed by the mode that
sprays it (`melt_profile.select_run`, the threshold `T_run`; NaN, today's liquid run, where no profile mode sprays), its
droplets carry the run's cells' mass-weighted density and viscosity (`melt_profile.run_properties`), and the slurry's own
mode reads the first flowing cell below that run. The cascade and the release take the mode's run (`mode_run`,
`cascade_mass(..., T_run=)`, `release_top(..., T_run=)`; a skin whose T_run is NaN is released exactly as on the table
path); the records carry the shared spray.PROFILE_COLUMNS after the 27, 36 wide, and the step summary the shared
`melt_profile.DISPERSION_COLUMNS`. On the table path, the default, nothing here changes, byte for byte."""
import types
from dataclasses import InitVar, dataclass

import numpy as np

from ..shared import melt_profile
from ..shared import spray as spray_mod
from ..shared.lubrication import lubrication
from ..shared.material import LiquidProperties
from .skin import CRUST_VISCOSITY, DZ_FLOOR, SHEAR_PASSES, SHEAR_TOL


# -- the cascade and the release (the skyscrapers' cascade_mass and release_top, matched by amendment) ---------------
@dataclass
class ReleaseBook:
    """What a release took from the skins' top cells, per skin row: `mass` [kg] (n,), `enthalpy` [J] (n,) at the cells'
    own specific enthalpies, `momentum` [kg m/s] (n, 3) at their own velocities (the droplets' momentum), `solid` [kg]
    (n,) by the cells' liquid fractions; `by_phase` the same as a ReleaseBook of (2, ...) arrays, [0] from the fully
    liquid cells and [1] from the others; `released_slurry` [kg] (n,) what the slurry's rate took (the cascade's)."""
    mass: np.ndarray
    enthalpy: np.ndarray
    momentum: np.ndarray
    solid: np.ndarray
    by_phase: object = None
    released_slurry: np.ndarray = None


def _phases(skins):
    T = skins.T()
    flowing = skins.flowing(T)
    return T, flowing, flowing & (T >= skins.material.T_liquidus)


def mode_run(skins, T_run):
    """The mode's release run in the skins' present cells (the profile path; the skyscrapers' `skin.mode_run`, matched by
    amendment), (n, K) bool: the shared `melt_profile.threshold_run` on the cells as the cascade reads them -- the
    flowing ones, and the liquid ones among them at or above the liquidus -- with each skin row's threshold T_run [K,
    (n,)]: the flowing cells from the top at least T_run warm, or where T_run is NaN the liquid run, today's."""
    T, flowing, liquid = _phases(skins)
    return melt_profile.threshold_run(T, flowing, liquid, T_run)


def _mode_cells(T, flowing, liquid, T_run):
    """(run, fast), (n, K) bool, for `cascade_mass` and `release_top` given thresholds T_run [K, (n,)] (the skyscrapers'
    `skin._mode_cells`): the mode's run (`mode_run`) and the cells at the mode's rate, booked as the liquid's mode's --
    the run where T_run is a number, and where it is NaN today's, every liquid cell, liquid beneath a slurry cell too, so
    that such a skin is released and booked as without thresholds, to the last bit."""
    T_run = np.broadcast_to(np.asarray(T_run, dtype=float), T.shape[:1])
    run = melt_profile.threshold_run(T, flowing, liquid, T_run)
    return run, np.where(np.isnan(T_run)[:, None], liquid, run)


def cascade_mass(skins, rate, rate_slurry, area, dt, slurry=False, T_slurry=None, T_run=None):
    """The mass [kg] (n,) the release takes from each skin row's top cells in dt [s] when each cell goes at its own
    phase's rate, and the part of it taken from slurry cells: down the run of cells the release may take (the liquid
    cells from the top, with `slurry` every flowing cell down to the first rigid one), in time order, a liquid cell at
    `rate` and a slurry cell at `rate_slurry` [kg/(m2 s) of top area] times the top `area` [m2], until dt is spent. A
    cell whose rate is zero ends the column's cascade, as rigid metal does, and so does a slurry cell colder than
    `T_slurry` [K] (n,), if given.

    With T_run [K, (n,)] (the profile path; the skyscrapers' `skin.cascade_mass(..., T_run=)`, matched by amendment) the
    release is heat-limited: the cells at `rate` are the mode's run (`mode_run`), the flowing cells from the top at
    least T_run warm, liquid or not. Beneath them, with `slurry`, the flowing cells go at `rate_slurry`, those at least
    T_slurry warm where it is given, the first colder one ending the skin's cascade; without it the cascade ends at the
    mode's last cell. The part returned second is then what went at `rate_slurry`, and `release_top(skins, mass, slurry,
    T_run)` takes exactly these cells. A skin whose T_run is NaN keeps today's cascade, to the last bit, whatever the
    others' thresholds."""
    T, flowing, liquid = _phases(skins)
    if T_run is None:
        fast = liquid
        run = np.cumprod(flowing if slurry else liquid, axis=1).astype(bool)
    else:
        mode, fast = _mode_cells(T, flowing, liquid, T_run)
        run = np.cumprod(flowing, axis=1).astype(bool) if slurry else mode
    rate_s = np.asarray(rate_slurry, dtype=float)[:, None]
    if T_slurry is not None:
        with np.errstate(invalid="ignore"):
            rate_s = np.where(T >= np.asarray(T_slurry, dtype=float)[:, None], rate_s, 0.0)
    r = np.where(fast, np.asarray(rate, dtype=float)[:, None], rate_s)
    r = r * np.asarray(area, dtype=float)[:, None]                         # kg/s from each cell while it is exposed
    m = skins.m
    with np.errstate(divide="ignore", invalid="ignore"):
        lasts = np.where(run & (r > 0.0), m / np.where(r > 0.0, r, 1.0), np.inf)      # s to strip each cell
    starts = np.concatenate([np.zeros((skins.n, 1)), np.cumsum(lasts, axis=1)[:, :-1]], axis=1)
    spare = np.where(run & (starts < dt), dt - starts, 0.0)
    taken = np.where(run, np.minimum(spare * r, m), 0.0)
    return taken.sum(axis=1), (taken * ~fast).sum(axis=1)


def release_top(skins, dm, slurry=False, T_run=None):
    """Release dm [kg] (n,) per skin row from its top cells: the run of liquid cells from the top, with `slurry` every
    flowing cell down to the first rigid one, never rigid metal; at most what those cells hold leaves, at its cells' own
    enthalpy and velocity (`Skins.take_top_mass`). Returns a ReleaseBook.

    With T_run [K, (n,)] (the profile path; the skyscrapers' `skin.release_top(..., T_run=)`, matched by amendment) the
    release takes the mode's run (`mode_run`), liquid or not, and with `slurry` every flowing cell, and `by_phase` splits
    it by the cells at the mode's rate (`cascade_mass`): [0] the mode's run, [1] the others. Without the slurry it takes
    nothing beneath the run. A skin whose T_run is NaN is released and split as without T_run, to the last bit."""
    T, flowing, liquid = _phases(skins)
    if T_run is None:
        admit, phase = (flowing if slurry else liquid), liquid
    else:
        run, phase = _mode_cells(T, flowing, liquid, T_run)
        admit = flowing if slurry else run
    mass, H, solid, momentum, ph = skins.take_top_mass(dm, admit, phase=phase)
    return ReleaseBook(mass, H, momentum, solid, ReleaseBook(ph[0], ph[1], ph[3], ph[2]))


# -- the skins' cells for the profile path (the skyscrapers' layer_columns, matched by amendment) -----------------------
def _series_drag(dz, valid, eta):
    """The viscous drag per unit area between neighbouring cells' centres [Pa s/m] (n, K-1) from the cells'
    viscosities eta [Pa s] (n, K): the two half-cells in series, 1 / (dz_i / (2 eta_i) + dz_j / (2 eta_j)); zero between
    two cells that do not flow and across a missing one (the skyscrapers' `skin.series_drag`, as the skins' own momentum
    step's `drags` reads it between cells)."""
    if dz.shape[1] < 2:
        return np.zeros((dz.shape[0], 0))
    resist = 0.5 * dz[:, :-1] / eta[:, :-1] + 0.5 * dz[:, 1:] / eta[:, 1:]
    link = valid[:, 1:] & (np.isfinite(eta[:, :-1]) | np.isfinite(eta[:, 1:]))
    return np.where(link, 1.0 / np.where(link, resist, 1.0), 0.0)


def _shear_rates(skins, tau, T, flowing):
    """Every cell's shear rate [1/s] (n, K) from the cells' present velocities, by the rule the skins' momentum step
    starts from (`Skins.momentum_step`, its `rates`; the skyscrapers' `skin.cell_shear_rates`, matched by amendment): a
    flowing cell's rate is the mean over its two faces of the stress there over its viscosity. The top face carries the
    gas shear |tau| [Pa, (n,)] where the top cell moves; a face between two cells their series drag times the size of
    the difference of their velocities; a flowing bottom cell's bottom face the drag of the rigid metal beneath the base
    over its half-thickness (no slip), times its speed. The velocities are those of the moving cells (the anchored rigid
    metal is at rest), and a floating crust's half-cells take CRUST_VISCOSITY. The viscosities (`Skins.viscosity`) are a
    fixed point from rest: at most SHEAR_PASSES passes after the first, ending when no flowing cell's viscosity moves by
    more than SHEAR_TOL, relatively; the rates returned are the last pass's. Zero in the cells that do not flow and
    beyond each column's cells."""
    n, K = skins.m.shape
    rows = np.arange(n)
    last = np.maximum(skins.count - 1, 0)
    valid = skins.valid
    dz = skins.dz(T)
    moving = valid & ~skins.anchored()
    crust = moving & ~flowing
    U = np.where(moving[..., None], skins.U, 0.0)
    jump = np.linalg.norm(U[:, :-1] - U[:, 1:], axis=2)
    top = np.where(moving[:, 0], np.abs(np.broadcast_to(np.asarray(tau, dtype=float), (n,))), 0.0)
    bottom_flows = flowing[rows, last] & (skins.count > 0)
    bottom_speed = np.linalg.norm(U[rows, last], axis=1)

    def viscosities(gamma):
        return np.where(crust, CRUST_VISCOSITY, skins.viscosity(gamma, T))

    def rates(eta):
        Db = np.where(bottom_flows, eta[rows, last] / (0.5 * np.maximum(dz[rows, last], DZ_FLOOR)), 0.0)
        stress = np.zeros((n, K + 1))                                     # at every face, the top face first
        stress[:, 1:K] = _series_drag(dz, valid, eta) * jump
        stress[:, 0] = top
        stress[rows, last + 1] = np.where(np.isfinite(Db), Db, 0.0) * bottom_speed
        face = 0.5 * (stress[:, :-1] + stress[:, 1:])
        return np.where(flowing, np.divide(face, eta, out=np.zeros_like(face), where=np.isfinite(eta)), 0.0)

    eta = viscosities(0.0)
    gamma = rates(eta)
    eta = viscosities(gamma)
    for _ in range(SHEAR_PASSES):
        gamma = rates(eta)
        eta_new = viscosities(gamma)
        change = np.abs(eta_new[flowing] - eta[flowing]) / eta[flowing]
        if change.size == 0 or change.max() <= SHEAR_TOL:
            break
        eta = eta_new
    return gamma


def layer_columns(skins, tau):
    """Every skin's cells as plain arrays for the shared straight-segment fit (`melt_profile.LayerColumns`; the
    skyscrapers' `skinbody.layer_columns`, matched by amendment), in the skins' row order: the cells' centre depths and
    thicknesses (`Skins.dz`, under the column's cross-section), each cell's velocity projected on the top cell's direction
    (zero where the top cell is at rest), temperatures, masses and densities (`Skins.density`), the viscosities
    `Skins.viscosity` at the cells' shear rates under the gas shear |tau| [Pa, (n,) by skin row] (`_shear_rates`, the
    rates the momentum step starts from), the contiguous runs from the top of flowing cells and of fully liquid ones
    (flowing and at or above the liquidus, as `cascade_mass` reads the liquid), and the rigid front `depth_above(T_rigid)`:
    the column's own depth where no rigid metal lies in it, so that a skin flowing to its base has its base as its front
    (design 21.7). Beyond each column's cells z, u, T and rho are NaN, dz and m 0 and mu inf."""
    mat = skins.material
    T, flowing, _ = _phases(skins)
    valid = skins.valid
    mu = skins.viscosity(_shear_rates(skins, tau, T, flowing), T)
    dz = skins.dz(T)
    top = skins.U[:, 0]
    speed = np.linalg.norm(top, axis=1)
    direction = np.divide(top, speed[:, None], out=np.zeros_like(top), where=speed[:, None] > 0.0)
    u = np.einsum("ikj,ij->ik", skins.U, direction)
    runs = np.cumprod(valid & np.isfinite(mu), axis=1).astype(bool)
    liquid = np.cumprod(runs & (T >= mat.T_liquidus), axis=1).astype(bool)
    rho = np.full(T.shape, np.nan)
    rho[valid] = skins.density(T[valid])
    cells = lambda a: np.where(valid, a, np.nan)
    return melt_profile.LayerColumns(z=cells(np.cumsum(dz, axis=1) - 0.5 * dz), dz=dz, u=cells(u), T=cells(T),
                                     m=np.where(valid, skins.m, 0.0), rho=rho, mu=mu, flowing=runs, liquid=liquid,
                                     front=skins.depth_above(mat.T_rigid))


# -- the plan ----------------------------------------------------------------------------------------------------------
@dataclass
class SprayPlan:
    """Girin's spraying from the markers' skins over one macro step, evaluated at its start (`spray_plan`), per marker in
    the markers' order: `rate` [kg/(m2 s) of top area] the liquid's mode's and `rate_slurry` the exposed slurry's;
    `result` the liquid's `spray.SprayResult`; `slurry` the exposed slurry's per marker (branch, r, wavelength, mu,
    gamma, rho, delta_m, depth, gamma_tau, T, we_s, h_floor; nan, -1 where there is none); the depths read at the start,
    `flowing` (b), `liquid` and `cap_depth` [m], and `on_slurry`; `skin_index` each marker's skin row; `slurry_spray`
    whether the release may take slurry. For the records: `theta` [rad], `closure`, `points` [m] (n, 3) the markers'
    positions, `state` the aero state, `props` the liquid's properties, `area` [m2] the skins' top areas, `centre` [m]
    and `v_hat` for Register et al.'s release angle (nan without them). The books, per marker and mode, filled by
    `release`: `mass`, `enthalpy`, `solid`, `momentum` and their `_slurry` twins; `number` and `number_slurry` the
    droplets (`count`); `hist_n`, `hist_m` the size histograms of every release at the plan's radii.

    On the profile path (module docstring; the skyscrapers' SprayPlan's fields of the same names): `T_run` [K] (n,) each
    marker's threshold for its mode's release run, which the cascade and the release take (NaN, today's liquid run, where
    no profile mode sprays); `resolved` (bool, (n,)) the markers whose profile was resolved; `profile` (n, 9) the
    records' appended columns (`melt_profile.profile_rows`, spray.PROFILE_COLUMNS); `modes` and `record` the plan's
    `spray.ProfileModes` and `melt_profile.ProfileRecord`, for the step summary; and `props` the droplets' properties per
    marker, the run's cells' density and viscosity where the marker is resolved and its run holds a cell. They are
    stored only when given, so that a plan on the table path holds, and pickles, what it did before them; read there
    they are None."""
    rate: np.ndarray
    rate_slurry: np.ndarray
    result: object
    slurry: dict
    flowing: np.ndarray
    liquid: np.ndarray
    cap_depth: np.ndarray
    on_slurry: np.ndarray
    skin_index: np.ndarray
    slurry_spray: bool = False
    theta: np.ndarray = None
    closure: np.ndarray = None
    points: np.ndarray = None
    state: object = None
    props: object = None
    area: np.ndarray = None
    centre: np.ndarray = None
    v_hat: np.ndarray = None
    T_run: InitVar[np.ndarray] = None
    resolved: InitVar[np.ndarray] = None
    profile: InitVar[np.ndarray] = None
    modes: InitVar[object] = None
    record: InitVar[object] = None

    def __post_init__(self, T_run, resolved, profile, modes, record):
        for name, value in (("T_run", T_run), ("resolved", resolved), ("profile", profile), ("modes", modes),
                            ("record", record)):
            if value is not None:
                setattr(self, name, value)
        n = len(self.rate)
        for key in ("mass", "mass_slurry", "number", "number_slurry", "enthalpy", "enthalpy_slurry", "solid",
                    "solid_slurry"):
            setattr(self, key, np.zeros(n))
        self.momentum, self.momentum_slurry = np.zeros((n, 3)), np.zeros((n, 3))
        self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)
        self.time, self.block = None, None

    def _to_rows(self, per_marker, n_rows, fill=0.0):
        out = np.full(n_rows, fill, dtype=float)
        out[self.skin_index] = per_marker
        return out

    def release(self, skins, dt):
        """The release of one flow sub-step dt [s] from `skins` (module docstring): the cascade at the plan's two rates
        (the slurry's for slurry cells at least as warm as the cell that gave it), then exactly those cells; on the
        profile path each skin's mode's run at its threshold `T_run`, found in the cells of this sub-step (NaN for a skin
        no marker of the plan owns: today's liquid run). Books what left by marker and mode, bins it by size, and
        returns the ReleaseBook in the skins' row order."""
        n = skins.n
        T_s = self.slurry["T"]
        rate = self._to_rows(self.rate, n)
        rate_s = self._to_rows(self.rate_slurry, n)
        floor = self._to_rows(np.where(np.isfinite(T_s), T_s, np.inf), n, np.inf)
        T_run = None if self.T_run is None else self._to_rows(self.T_run, n, np.nan)
        if not np.all(np.isfinite(rate) & (rate >= 0.0) & np.isfinite(rate_s) & (rate_s >= 0.0)):
            raise ValueError("a release rate must be finite and non-negative")
        dm, from_slurry = cascade_mass(skins, rate, rate_s, skins.area, dt, self.slurry_spray, floor, T_run)
        book = release_top(skins, dm, self.slurry_spray, T_run)
        book.released_slurry = from_slurry
        ph, r = book.by_phase, self.skin_index
        self.mass += ph.mass[0][r]
        self.mass_slurry += ph.mass[1][r]
        self.enthalpy += ph.enthalpy[0][r]
        self.enthalpy_slurry += ph.enthalpy[1][r]
        self.solid += ph.solid[0][r]
        self.solid_slurry += ph.solid[1][r]
        self.momentum += ph.momentum[0][r]
        self.momentum_slurry += ph.momentum[1][r]
        hn, hm = self.histogram(np.vstack([ph.mass[0][r], ph.mass[1][r]]))
        self.hist_n += hn
        self.hist_m += hm
        return book

    def _radii(self):
        return np.vstack([self.result.r, self.slurry["r"]])

    def droplet_mass(self):
        """One droplet's mass [kg] (2, n) at each mode's radius, the liquid's mode's at the plan's `props` density (the
        liquid's on the table path, each marker's run's cells' on the profile path, as the skyscrapers count with
        `p.props.rho`), the slurry's at its cell's (`spray.droplet_mass`)."""
        r = self._radii()
        rho = np.broadcast_to(np.asarray(self.props.rho, dtype=float), r.shape[1:])
        return spray_mod.droplet_mass(r, np.vstack([rho, self.slurry["rho"]]))

    def count(self):
        """The droplets of the record: each mode's mass over one droplet's at its radius (`spray.droplet_count`)."""
        one = self.droplet_mass()
        self.number = spray_mod.droplet_count(self.mass, one[0])
        self.number_slurry = spray_mod.droplet_count(self.mass_slurry, one[1])

    def histogram(self, mass):
        """(droplets, mass [kg]) per size bin of a release by mode [kg] (2, n) at the plan's radii
        (`spray.droplet_histogram`)."""
        return spray_mod.droplet_histogram(np.asarray(mass, dtype=float), self.droplet_mass(), self._radii())

    def spraying(self):
        """The markers whose plan sprays, bool (n,): a mode with a release rate."""
        return spray_mod.spraying(self.rate, self.rate_slurry)

    def release_angle(self):
        """Register et al.'s release angle [deg] (n,) of each marker from `centre`; nan without one."""
        centre = None if self.v_hat is None else self.centre
        return spray_mod.plan_release_angle(self.points, centre, self.v_hat, len(self.rate))

    def records(self, t, state=None):
        """The macro step's records as source-table rows (`spray.source_block`, SOURCE_COLUMNS), one per marker and mode
        that released, the liquid's mode's first, at the step's end `t` [s] and the aero `state` (default the plan's):
        the marker's position and angle where the plan was evaluated, its closure; the mode's branch, radius, delta_m,
        We_s, the viscous film's floor and gamma tau*; b (the liquid's: the plan's b; the slurry's: its own layer); the
        droplets (`count`) and their mass; their speed and direction the released metal's momentum over its mass; We_d,
        Oh and breakup by their phase's properties; their solid fraction and enthalpy per kilogram; the release angle.
        On the profile path each row ends with its marker's `profile` columns (spray.PROFILE_COLUMNS, the rows
        `spray.record_columns(36)` wide), the slurry's with no release depth and no growth rate, which are the liquid's
        mode's. Keeps the rows and `t` for `step_summary`."""
        self.count()
        st = self.state if state is None else state
        s = self.slurry
        extra = extra_slurry = self.profile
        if extra is not None:
            # Ruling T10-a of plan 03b Task 10: the skyscrapers' rule, written again here because the level set may not
            # import them and moving it into shared/ needs a skyscrapers edit this plan forbids (a follow-up on main):
            # their SprayPlan.source_block (skyscrapers/skinbody.py, the lines that build `extra_slurry` from `extra`)
            # blanks the slurry's rows' release depth and growth rate. The parity test's equal records pin the two.
            extra_slurry = np.array(extra, dtype=float)
            extra_slurry[:, [spray_mod.PROFILE_COLUMNS.index("release_depth_m"),
                             spray_mod.PROFILE_COLUMNS.index("growth_rate_per_s")]] = np.nan
        modes = ((self.mass, self.enthalpy, self.solid, self.momentum, self.number, self.result, self.flowing, self.props),
                 (self.mass_slurry, self.enthalpy_slurry, self.solid_slurry, self.momentum_slurry, self.number_slurry,
                  types.SimpleNamespace(**{k: s[k] for k in ("branch", "r", "delta_m", "we_s", "h_floor", "gamma_tau")}),
                  s["depth"], LiquidProperties(s["rho"], s["mu"], spray_mod.SLURRY_SIGMA)))
        block = spray_mod.plan_source_block(t, st, self.theta, self.points, self.closure, modes, self.release_angle(),
                                            extra=extra, extra_slurry=extra_slurry)
        self.time, self.block, self._state_out = t, block, st
        return self.block

    def step_summary(self, skins=None, n_before=0.0):
        """The history's droplet columns of the macro step (every key of the shared writer's SKIN_SUMMARY_COLUMNS, and
        the release angle's), after `records`: time_s, altitude_km, velocity_kms; released_step_mass_kg the step's
        release; n_released the droplets since time zero, `n_before` plus the step's (the skyscrapers' history counts
        from time zero); the radius median and largest over the rows [um]; theta_cr_deg the smallest angle and
        spraying_area_m2 the top area of the markers whose plan sprays; the release angle's mass-weighted mean, least
        and largest; liquid_depth_mean_mm and flowing_depth_mean_mm over every skin of `skins` at the step's end (the
        plan's depths at its start without them). On the profile path, last, the shared `melt_profile.summary` of the
        plan's modes and record (melt_profile.DISPERSION_COLUMNS, the skyscrapers' `dispersion_stats`), the liquid's
        mode's mass by the branch each marker took."""
        if self.block is None:
            raise ValueError("step_summary reads the step's records: call records(t) first")
        st = self._state_out
        out = spray_mod.step_columns(np.asarray(self.theta), self.area, self.spraying(), self.block)
        if skins is not None:
            mat = skins.material
            flowing, liquid = skins.depth_above(mat.T_rigid), skins.depth_above(mat.T_liquidus)
        else:
            flowing, liquid = self.flowing, self.liquid
        out.update({"time_s": float(self.time), "altitude_km": float(st.h) / 1e3, "velocity_kms": float(st.V) / 1e3,
                    "released_step_mass_kg": float(self.mass.sum() + self.mass_slurry.sum()),
                    "n_released": float(n_before) + float(self.number.sum() + self.number_slurry.sum()),
                    "liquid_depth_mean_mm": float(np.mean(liquid) * 1e3) if len(liquid) else 0.0,
                    "flowing_depth_mean_mm": float(np.mean(flowing) * 1e3) if len(flowing) else 0.0})
        if self.modes is not None:
            out.update(melt_profile.summary(self.modes, self.record, self.result.branch, self.mass))
        return out


def _region_edges(division):
    """The wave's domain graph (module docstring): the division's shared boundaries with a length at phi, (E, 2) markers'
    indices i < j; ValueError without a division."""
    if division is None:
        raise ValueError("the wave's domain graph is the division's boundaries at phi: give the division of the layer "
                         "among the markers (division=...)")
    at_phi = np.asarray(division.length, dtype=float)[:, 0] > 0.0
    return np.column_stack([np.asarray(division.i)[at_phi], np.asarray(division.j)[at_phi]]).astype(np.int64)


def spray_plan(skins, markers, flow, state, model, windward, w_n, dt, radius, nose_radius, slurry=False, theta=None,
               centre=None, v_hat=None, division=None):
    """Girin's spraying from each marker's skin for the macro step about to run (module docstring): a SprayPlan, per
    marker in the markers' order. `skins` the markers' skins (marker i's is row `markers.skin_index[i]`); `markers`
    their markers; `division` the division of the layer among them (`division.Division`, its boundaries at phi the
    wave's domain graph; needed); `flow` the shared SurfaceFlowResult per marker; `state` the
    aero state; `model` the shared SprayModel; `windward` bool (n,) and `w_n` [m/s2] (n,) the normal deceleration
    W cos(theta) (the loads' leeward test and angles); `dt` [s] the macro step; `radius` [m] the body radius of the
    quarter cap; `nose_radius` [m] the radius the slip band's integrated shear was taken at. `slurry` the exposed slurry's
    mode too (`--slurry-spray on`); `theta` [rad] (n,) the markers' angles to the flow (the slip band and the records; nan if not given); `centre` [m]
    and `v_hat` for the records' release angle.

    With `model.dispersion` "profile" (module docstring; the skyscrapers' profile branch of `girin_plan`, matched by
    amendment) the windward markers whose skin's top cell flows and moves are resolved: their skins' cells are read as
    `layer_columns` under the gas shear |tau|, and the shared `melt_profile.profile_modes` fits their profiles, places
    their walls and solves their gradient mode, for the spray model's profile path (`model.evaluate(..., profile=)`).
    There the wave must fit in the contiguous region of windward markers whose layer above the wall -- the flowing
    depth on the markers not resolved -- is at least spray.B_MIN deep, and the Rayleigh-Taylor mode's in the liquid's
    region, as on the table path; a marker not resolved keeps the table path, its wave in the liquid's region too. Each
    marker's release run is fixed by the mode that sprays it (`melt_profile.select_run`, a mode with no rate fixing
    none: today's liquid run), its droplets carry the run's cells' mass-weighted density and viscosity
    (`melt_profile.run_properties`; the liquid's on the markers not resolved, and where the run is empty), and the
    slurry's own mode reads the first flowing cell below that run."""
    edges = _region_edges(division)
    if model.slip_band and theta is None:
        raise ValueError("the slip band's depth reads the markers' angles: give theta with the model's slip band on")
    mat, liq = skins.material, model.liquid
    rows = np.asarray(markers.skin_index, dtype=np.int64)
    n = len(rows)
    T, flowing, liquid = _phases(skins)
    F_skin = skins.depth_above(mat.T_rigid)[rows]
    L = skins.depth_above(mat.T_liquidus)[rows]
    b = F_skin
    n_flow = np.cumprod(flowing, axis=1).sum(axis=1)[rows]
    n_liq = np.cumprod(liquid, axis=1).sum(axis=1)[rows]
    on_slurry = n_flow > n_liq
    windward = np.asarray(windward, dtype=bool).reshape(n)
    w_n = np.asarray(w_n, dtype=float).reshape(n)
    theta_arr = np.full(n, np.nan) if theta is None else np.asarray(theta, dtype=float).reshape(n)
    delta_m = spray_mod.spray_layer(flow, liq, model.slip_band, theta, nose_radius)
    v_s = lubrication(flow.tau, flow.G, b, delta_m, liq.mu, b_layer=b)[0]
    area = skins.area[rows]
    extent = spray_mod.region_extent(edges, L > 0.0, area)
    run = np.cumprod(flowing if slurry else liquid, axis=1).astype(bool)
    m_f = (skins.m * run).sum(axis=1)[rows]
    cap = F_skin if slurry else L
    if model.dispersion != "profile":
        res = model.evaluate(flow, state, b, delta_m, v_s, windward, dt, area, m_f, radius, b_layer=b, extent=extent,
                             deceleration_n=w_n, regime_layer=b, on_slurry=on_slurry, rt_layer=L, cap_depth=cap)
        n_above, props, on_profile = n_liq, liq, {}
    else:                                    # the profile path (the skyscrapers' girin_plan, matched by amendment)
        tau = np.zeros(skins.n)
        tau[rows] = np.abs(np.asarray(flow.tau, dtype=float))
        cols = layer_columns(skins, tau).rows(rows)
        pm, rec = melt_profile.profile_modes(cols, windward, liq.sigma)
        D = np.where(pm.resolved, pm.wall_depth, F_skin)
        waves = np.where(pm.resolved, spray_mod.region_extent(edges, windward & (D >= spray_mod.B_MIN), area), extent)
        res = model.evaluate(flow, state, b, delta_m, v_s, windward, dt, area, m_f, radius, b_layer=b, extent=waves,
                             deceleration_n=w_n, regime_layer=b, on_slurry=on_slurry, rt_layer=L, cap_depth=cap,
                             profile=pm, rt_extent=extent)
        T_run = melt_profile.select_run(np.where(res.mdot > 0.0, res.branch, -1), pm, rec)
        T_rows = np.full(skins.n, np.nan)
        T_rows[rows] = T_run
        mode = mode_run(skins, T_rows)[rows]
        n_above = mode.sum(axis=1)
        rho_run, mu_run = melt_profile.run_properties(cols, mode)
        own = pm.resolved & (n_above > 0)
        props = LiquidProperties(np.where(own, rho_run, liq.rho), np.where(own, mu_run, liq.mu), liq.sigma)
        on_profile = dict(T_run=T_run, resolved=pm.resolved, profile=melt_profile.profile_rows(rec, delta_m, v_s),
                          modes=pm, record=rec)
    if slurry:
        rate_s, modes = _slurry_modes(skins, rows, T, flow, state, model, windward, w_n, area, edges, radius, dt, b,
                                      F_skin, n_above, n_flow, theta, nose_radius)
    else:
        rate_s, modes = np.zeros(n), spray_mod.no_slurry(n)
    return SprayPlan(res.mdot, rate_s, res, modes, b, L, cap, on_slurry, rows, bool(slurry), theta=theta_arr,
                     closure=np.asarray(flow.closure), points=np.asarray(markers.x, dtype=float).copy(), state=state,
                     props=props, area=np.asarray(area, dtype=float), centre=None if centre is None else np.asarray(centre, float),
                     v_hat=None if v_hat is None else np.asarray(v_hat, dtype=float), **on_profile)


def _slurry_modes(skins, rows, T, flow, state, model, windward, w_n, area, edges, radius, dt, b, F_skin, n_liq, n_flow,
                  theta, nose_radius):
    """The exposed slurry's own mode per marker, the shared `spray.slurry_modes` (moved there from the skyscrapers for
    stage 3 Task 8): this reads the skins for it -- where a flowing cell that is not fully liquid lies beneath the liquid
    run on the windward face, that cell's temperature, solid fraction and mixture density, its layer from the cell's top
    to T_rigid, the slurry run's mass, the flowing region (b at least B_MIN) the wave must fit in and the cap (the
    flowing depth) -- with the skins' SHEAR_PASSES and SHEAR_TOL for the fixed point on the shear rate. On the profile
    path `n_liq` is the count of the mode's run's cells (`spray_plan`), so the slurry's cell is the first flowing one
    below that run, as the skyscrapers pass it. Returns (rate [kg/(m2 s)], the per-marker dict)."""
    sel = np.flatnonzero((n_flow > n_liq) & windward)
    if not sel.size:
        return np.zeros(len(rows)), spray_mod.no_slurry(len(rows))
    rr = rows[sel]
    k = n_liq[sel]
    T_cell = T[rr, k]
    dz = skins.dz(T)
    depth = F_skin[sel] - (np.cumsum(dz, axis=1) - dz)[rr, k]
    f_s = 1.0 - skins.props.liquid_fraction(T_cell)
    rho = skins.density(T_cell)
    extent = spray_mod.region_extent(edges, b >= spray_mod.B_MIN, area)
    cells = np.arange(skins.m.shape[1])[None, :]
    n_liq_r, n_flow_r = np.zeros(skins.n, dtype=np.int64), np.zeros(skins.n, dtype=np.int64)
    n_liq_r[rows], n_flow_r[rows] = n_liq, n_flow
    m_s = (skins.m * ((cells >= n_liq_r[:, None]) & (cells < n_flow_r[:, None]))).sum(axis=1)[rows]   # the slurry run
    return spray_mod.slurry_modes(sel, T_cell, f_s, rho, depth, flow, state, model, skins.material, windward, w_n, area,
                                  m_s, extent, F_skin, radius, dt, SHEAR_PASSES, SHEAR_TOL, theta=theta,
                                  nose_radius=nose_radius, law=skins.settings.rheology)
