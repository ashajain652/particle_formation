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
  `depth_above(T_rigid)` (b_mode "skin"), or on a top-zone skin (b_mode "phi_psi") the depth from phi to psi given per
  marker (nan where the marker's skin is full-depth, which then reads its own flowing depth; design 20.6);
- delta_m from `spray_layer` (the slip band through it where the model's `slip_band` is on), v_s from
  `lubrication(tau, G, b, delta_m, mu_l, b_layer=b)`;
- the droplet's cap depth: the skin's liquid depth `depth_above(T_liquidus)`, or with `slurry` its flowing depth; the
  Rayleigh-Taylor pool the liquid depth; the film mass the run of releasable cells' mass; the wave's extent the
  contiguous region of markers with liquid on the markers' face graph (`Markers.faces`);
- with `slurry`, the exposed slurry's own mode per marker (`_slurry_modes`), a second `SprayModel` with the slurry's
  properties at the shear rate of its own sheared layer, found by fixed point.

On a top-zone skin (stage 3 Task 10, decision D2 and the controller's ruling (b) of 2026-11-12) the melt beneath is the
volume's, held, and its depth is phi to psi: b, the cap depth, the Rayleigh-Taylor pool, and the slurry's layer (from
its cell's top down to psi), its cap and the region its wave must fit in (b at least B_MIN) all read it. The film mass
and the release stay the skin's own (`Skins.take_top_mass`): what leaves comes from the skin's top cells, and the skin's
base draw then refills it from the melt beneath (`modes.draw_cap`).

No renewal limit (ruling of 2026-11-11): design 9.5's cap on liquid deeper than delta_m, rho_l delta_m A per renewal time
delta_m^2/nu_l, belongs to the two-zone skin of 2026-10-07 and is not built; the layered skins resolve the velocity
through the depth, so the release is capped only as the skyscrapers cap it, by the cap depth and the cascade.

The release (`SprayPlan.release(skins, dt)`, every flow sub-step): down each skin's top cells in time order, a liquid
cell (at or above the liquidus) at the liquid's rate and, with `slurry`, a flowing slurry cell at the slurry's own rate
while it is at least as warm as the cell that gave that rate, until dt is spent; never rigid metal (`cascade_mass`);
then exactly those cells' metal leaves at its own enthalpy and velocity (`release_top`, `Skins.take_top_mass`). The
plan books what left by the cells it came from (the liquid's mode's and the slurry's), and at the step's end gives the
records' rows (`records`) and the history's droplet columns (`step_summary`, every key the shared writer's
SKIN_SUMMARY_COLUMNS reads)."""
import types
from dataclasses import dataclass

import numpy as np

from ..shared import spray as spray_mod
from ..shared.lubrication import lubrication
from ..shared.material import LiquidProperties
from .skin import SHEAR_PASSES, SHEAR_TOL

B_MODES = ("skin", "phi_psi")


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


def cascade_mass(skins, rate, rate_slurry, area, dt, slurry=False, T_slurry=None):
    """The mass [kg] (n,) the release takes from each skin row's top cells in dt [s] when each cell goes at its own
    phase's rate, and the part of it taken from slurry cells: down the run of cells the release may take (the liquid
    cells from the top, with `slurry` every flowing cell down to the first rigid one), in time order, a liquid cell at
    `rate` and a slurry cell at `rate_slurry` [kg/(m2 s) of top area] times the top `area` [m2], until dt is spent. A
    cell whose rate is zero ends the column's cascade, as rigid metal does, and so does a slurry cell colder than
    `T_slurry` [K] (n,), if given."""
    T, flowing, liquid = _phases(skins)
    run = np.cumprod(flowing if slurry else liquid, axis=1).astype(bool)
    rate_s = np.asarray(rate_slurry, dtype=float)[:, None]
    if T_slurry is not None:
        with np.errstate(invalid="ignore"):
            rate_s = np.where(T >= np.asarray(T_slurry, dtype=float)[:, None], rate_s, 0.0)
    r = np.where(liquid, np.asarray(rate, dtype=float)[:, None], rate_s)
    r = r * np.asarray(area, dtype=float)[:, None]                         # kg/s from each cell while it is exposed
    m = skins.m
    with np.errstate(divide="ignore", invalid="ignore"):
        lasts = np.where(run & (r > 0.0), m / np.where(r > 0.0, r, 1.0), np.inf)      # s to strip each cell
    starts = np.concatenate([np.zeros((skins.n, 1)), np.cumsum(lasts, axis=1)[:, :-1]], axis=1)
    spare = np.where(run & (starts < dt), dt - starts, 0.0)
    taken = np.where(run, np.minimum(spare * r, m), 0.0)
    return taken.sum(axis=1), (taken * ~liquid).sum(axis=1)


def release_top(skins, dm, slurry=False):
    """Release dm [kg] (n,) per skin row from its top cells: the run of liquid cells from the top, with `slurry` every
    flowing cell down to the first rigid one, never rigid metal; at most what those cells hold leaves, at its cells' own
    enthalpy and velocity (`Skins.take_top_mass`). Returns a ReleaseBook."""
    _, flowing, liquid = _phases(skins)
    mass, H, solid, momentum, ph = skins.take_top_mass(dm, flowing if slurry else liquid, phase=liquid)
    return ReleaseBook(mass, H, momentum, solid, ReleaseBook(ph[0], ph[1], ph[3], ph[2]))


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
    droplets (`count`); `hist_n`, `hist_m` the size histograms of every release at the plan's radii."""
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

    def __post_init__(self):
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
        (the slurry's for slurry cells at least as warm as the cell that gave it), then exactly those cells. Books what
        left by marker and mode, bins it by size, and returns the ReleaseBook in the skins' row order."""
        n = skins.n
        T_s = self.slurry["T"]
        rate = self._to_rows(self.rate, n)
        rate_s = self._to_rows(self.rate_slurry, n)
        floor = self._to_rows(np.where(np.isfinite(T_s), T_s, np.inf), n, np.inf)
        if not np.all(np.isfinite(rate) & (rate >= 0.0) & np.isfinite(rate_s) & (rate_s >= 0.0)):
            raise ValueError("a release rate must be finite and non-negative")
        dm, from_slurry = cascade_mass(skins, rate, rate_s, skins.area, dt, self.slurry_spray, floor)
        book = release_top(skins, dm, self.slurry_spray)
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
        """One droplet's mass [kg] (2, n) at each mode's radius, the liquid's at its density, the slurry's at its cell's
        (`spray.droplet_mass`)."""
        r = self._radii()
        return spray_mod.droplet_mass(r, np.vstack([np.full(r.shape[1], float(self.props.rho)), self.slurry["rho"]]))

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
        Keeps the rows and `t` for `step_summary`."""
        self.count()
        st = self.state if state is None else state
        s = self.slurry
        modes = ((self.mass, self.enthalpy, self.solid, self.momentum, self.number, self.result, self.flowing, self.props),
                 (self.mass_slurry, self.enthalpy_slurry, self.solid_slurry, self.momentum_slurry, self.number_slurry,
                  types.SimpleNamespace(**{k: s[k] for k in ("branch", "r", "delta_m", "we_s", "h_floor", "gamma_tau")}),
                  s["depth"], LiquidProperties(s["rho"], s["mu"], spray_mod.SLURRY_SIGMA)))
        block = spray_mod.plan_source_block(t, st, self.theta, self.points, self.closure, modes, self.release_angle())
        self.time, self.block, self._state_out = t, block, st
        return self.block

    def step_summary(self, skins=None, n_before=0.0):
        """The history's droplet columns of the macro step (every key of the shared writer's SKIN_SUMMARY_COLUMNS, and
        the release angle's), after `records`: time_s, altitude_km, velocity_kms; released_step_mass_kg the step's
        release; n_released the droplets since time zero, `n_before` plus the step's (the skyscrapers' history counts
        from time zero); the radius median and largest over the rows [um]; theta_cr_deg the smallest angle and
        spraying_area_m2 the top area of the markers whose plan sprays; the release angle's mass-weighted mean, least
        and largest; liquid_depth_mean_mm and flowing_depth_mean_mm over every skin of `skins` at the step's end (the
        plan's depths at its start without them)."""
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
        return out


def spray_plan(skins, markers, flow, state, model, windward, w_n, dt, radius, nose_radius, slurry=False, b_mode="skin",
               depth_phi_psi=None, top=None, theta=None, centre=None, v_hat=None):
    """Girin's spraying from each marker's skin for the macro step about to run (module docstring): a SprayPlan, per
    marker in the markers' order. `skins` the markers' skins (marker i's is row `markers.skin_index[i]`); `markers`
    tessellated (its faces are the wave's domain graph); `flow` the shared SurfaceFlowResult per marker; `state` the
    aero state; `model` the shared SprayModel; `windward` bool (n,) and `w_n` [m/s2] (n,) the normal deceleration
    W cos(theta) (the loads' leeward test and angles); `dt` [s] the macro step; `radius` [m] the body radius of the
    quarter cap; `nose_radius` [m] the radius the slip band's integrated shear was taken at. `slurry` the exposed slurry's
    mode too (`--slurry-spray on`); `b_mode` "skin" or "phi_psi" with `top` (n,) bool, the markers whose skins are
    top-zone (`Markers.top`), and `depth_phi_psi` [m] (n,), the depth from phi to psi, read where `top` (module
    docstring);
    `theta` [rad] (n,) the markers' angles to the flow (the slip band and the records; nan if not given); `centre` [m]
    and `v_hat` for the records' release angle."""
    if b_mode not in B_MODES:
        raise ValueError("b_mode must be one of {}, got {!r}".format(B_MODES, b_mode))
    if markers.faces is None:
        raise ValueError("the markers have no faces: tessellate them first")
    if model.slip_band and theta is None:
        raise ValueError("the slip band's depth reads the markers' angles: give theta with the model's slip band on")
    mat, liq = skins.material, model.liquid
    rows = np.asarray(markers.skin_index, dtype=np.int64)
    n = len(rows)
    T, flowing, liquid = _phases(skins)
    F_skin = skins.depth_above(mat.T_rigid)[rows]
    L = skins.depth_above(mat.T_liquidus)[rows]
    if b_mode == "phi_psi":
        if top is None or depth_phi_psi is None:
            raise ValueError("b_mode 'phi_psi' needs the top-zone markers (top) and the depth from phi to psi per marker")
        top = np.asarray(top, dtype=bool).reshape(n)
        d = np.asarray(depth_phi_psi, dtype=float).reshape(n)
        if not np.all(np.isfinite(d[top]) & (d[top] >= 0.0)):
            raise ValueError("b_mode 'phi_psi' needs a finite, non-negative depth from phi to psi at every top-zone marker")
        b = np.where(top, d, F_skin)
    else:
        top = np.zeros(n, dtype=bool)
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
    edges = np.column_stack([markers.faces.i, markers.faces.j]).astype(np.int64)
    extent = spray_mod.region_extent(edges, L > 0.0, area)
    run = np.cumprod(flowing if slurry else liquid, axis=1).astype(bool)
    m_f = (skins.m * run).sum(axis=1)[rows]
    cap = np.where(top, b, F_skin if slurry else L)
    res = model.evaluate(flow, state, b, delta_m, v_s, windward, dt, area, m_f, radius, b_layer=b, extent=extent,
                         deceleration_n=w_n, regime_layer=b, on_slurry=on_slurry, rt_layer=np.where(top, b, L),
                         cap_depth=cap)
    if slurry:
        rate_s, modes = _slurry_modes(skins, rows, T, flow, state, model, windward, w_n, area, edges, radius, dt, b,
                                      np.where(top, b, F_skin), n_liq, n_flow, theta, nose_radius)
    else:
        rate_s, modes = np.zeros(n), spray_mod.no_slurry(n)
    return SprayPlan(res.mdot, rate_s, res, modes, b, L, cap, on_slurry, rows, bool(slurry), theta=theta_arr,
                     closure=np.asarray(flow.closure), points=np.asarray(markers.x, dtype=float).copy(), state=state,
                     props=liq, area=np.asarray(area, dtype=float), centre=None if centre is None else np.asarray(centre, float),
                     v_hat=None if v_hat is None else np.asarray(v_hat, dtype=float))


def _slurry_modes(skins, rows, T, flow, state, model, windward, w_n, area, edges, radius, dt, b, F_skin, n_liq, n_flow,
                  theta, nose_radius):
    """The exposed slurry's own mode per marker, the shared `spray.slurry_modes` (moved there from the skyscrapers for
    stage 3 Task 8): this reads the skins for it -- where a flowing cell that is not fully liquid lies beneath the liquid
    run on the windward face, that cell's temperature, solid fraction and mixture density, its layer from the cell's top
    to T_rigid, the slurry run's mass, the flowing region (b at least B_MIN) the wave must fit in and the cap (the
    flowing depth; on a top-zone skin `F_skin` is phi to psi, module docstring) -- with the skins' SHEAR_PASSES and SHEAR_TOL for the fixed point on the shear rate. Returns (rate
    [kg/(m2 s)], the per-marker dict)."""
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
                                  nose_radius=nose_radius)
