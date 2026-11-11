"""The body of the layered-skin design (docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md sections 4,
6.6 and 7; stage 2, sub-plan 20: heat only, liquid removed as it forms).

A moving 3D mesh holds the rigid metal; every corner of its surface carries a column with a layered skin holding the
metal outward of the rigid front (`skin.Skins`). At time zero each column takes its minimum mass from the 3D model at
the initial temperature: the corners move inward by it and the skins receive exactly what they swept. Per macro step,
after the trajectory and the heating:

1. the trial pass gives every column its base law, E_base / dt = a + b T_b, with the base held at the 3D model's
   start-of-step temperature (`skin.advance`);
2. the 3D step takes that law at the corner nodes, implicitly, with no convective load and no radiation on its
   boundary (the skins' tops radiate), and last step's deferred mismatch as a nodal load;
3. the real pass, base held at the new 3D temperature, books the heat, removes the liquid and draws to the minimum;
   the difference between what it sent through the base and what the 3D step received is deferred to the node. By
   default ("damped", Asha, 2026-10-12) the step repeats, each pass re-linearised at a temperature moved only part of
   the way to the last 3D answer, until the total mismatch is within INTERFACE_TOLERANCE of the throughput; the other
   rules repeat at most once (REPEAT_RULES). Every draw takes the next layer down the 3D model's linear profile below
   the base, at that layer's enthalpy (`skin.profile_enthalpy`), the profile's slope taken at the step's start;
4. the base moves by the metal crossing it (spec 6.6): drawn (in the sub-steps) and absorbed metal move the corner in,
   handed-down metal moves it out; with the volume-neutral smoothing this is one displacement per corner, which
   `MovingMesh.move_with_policy` applies (rebuilding when needed), and each column is then reconciled with what its
   corner actually swept: its bottom gains the swept mass less what it already exchanged, and its enthalpy changes by
   the swept enthalpy less what it already exchanged, so mass and energy are exact. A sweep that takes metal out of
   the 3D model is booked along the same profile as the draws.

A macro step is run as equal parts so that no base moves more than MOVE_MAX of the mean surface edge in one part
(Asha, 2026-10-12): a part whose draws and hand-downs would move a base further is run again as shorter parts, and the
next macro step starts with as many parts as this one's moves called for. Absorption takes what is left of a part's
allowance after its draws and hand-downs, and the rest of the non-rigid metal waits for the next part (Asha,
2026-10-14: a cap of a quarter-edge per macro step stranded up to 1.3 mm of non-rigid metal in the 3D model of a plate
started at 880 K, whose exact non-rigid layer is 1.5-1.7 mm thick, and the plate melted 4-8 % slow). While any skin's top is non-rigid (at or above T_rigid at
the macro step's start), no part is longer than PART_MAX (Asha, 2026-10-12: long parts during melting were rejected for
their moves and needed many passes), and also in the macro step in which a top is about to become non-rigid: its rise
over the last macro step, carried on across this one, reaches T_rigid (Asha, 2026-10-16: with the skins' 0.8 mm minimum
the level-2 cold plate's top stood at 884.8 K at 0.5 s, so the step in which it melted ran as one part, and the plate
was 14.2 % ahead at 1.0 s; capped one step early, 3.2 %). Within a part the base heat, the draws and the move all take the part's own step;
the heating is interpolated in time between the last macro step's and this one's, at the part's middle, so a macro
step receives the trapezoid of the two (Asha, 2026-10-12).

The body's surface for the coupled run is the top surface, one 'patch' per corner (`columns.CornerSurface`), so the
heating, the history and the trajectory's loop are unchanged.

The flow (sub-plan 21 Task 6c; `flow`, off by default): every macro step the surface flow is evaluated on the corners
as the heating is, from the step's aero state and the tops' temperatures -- its wall-pressure gradient along the
meridian and its gas shear -- with the deceleration, which the metal feels along the flight direction; each part takes
them interpolated in time at its middle, as it takes the heating, and every sub-step of both passes runs the flow after
its conduction (`skin.flow_substeps`), `flow_cycles` flow sub-steps to each. With the flow on the liquid stays and flows
(`remove` False in a run; spraying is stage 4's). The momentum books: the skins' momentum changes by the flow's impulses
less what leaves with removed or handed-down metal, and `momentum_residual` reports what they miss."""
import dataclasses
import math
import time

import numpy as np

from . import columns as columns_mod
from . import movemesh
from . import skin as skin_mod
from .body import ThermalBody
from .thermal import SIGMA_SB

FLOW_DEFAULT = False              # the flow is off until stage 4's spraying gives the liquid a sink: with it on and nothing
                                  # sprayed, a flight would only pile the liquid up (Task 6c's ruling, 2026-10-22)

GRADIENT_FLOOR = 1.0              # K/m: the smallest gradient along the normal that is divided by. A corner above T_rigid
                                  # whose gradient is below it (a flat or reversed field) is given this one, so its depth
                                  # is huge and it absorbs up to its part's move allowance: it does not absorb nothing
REPEAT_RULES = ("column", "total", "never", "damped")   # when a step's interface mismatch repeats it: "damped"
                                              # repeats until converged (below); the others at most once -- "column" (the
                                              # spec's first rule) when any column's mismatch exceeds the tolerance,
                                              # "total" when the total does -- or never (sub-plan 20 Task 10)
DEFAULT_REPEAT_RULE = "damped"    # Asha, 2026-10-12: the only rule whose coupling converges at 0.5 s macro steps
# The damped repetition (2026-10-12). The skins' base heat has a cliff a little above the start-of-step base
# temperature once the skins near melting: a base held above the melting band for a whole step melts the skin, which
# then draws metal at T_rigid and melts that too, so the heat sent down falls from +5 to -7 MW/m2 over 30 K (the
# plate's centre at 0.5 s, measured). A base law fitted at one temperature cannot follow it, and the 3D surface, whose
# admittance is about four times the law's slope, overshoots it: undamped repetition alternates between two states.
# Damping moves the linearisation temperature only a fraction omega of the way to the 3D answer, Aitken's method
# choosing omega from the last two passes.
DAMPED_MAX_PASSES = 12            # passes of a macro step at most (the last pass's mismatch is deferred as before)
DAMPED_OMEGA0 = 0.5               # the first pass's relaxation factor
DAMPED_OMEGA_RANGE = (0.02, 1.0)  # Aitken's factor kept within these
DAMPED_TOL_K = 0.05               # K: converged when no base temperature moves more than this between passes
# The move limit (2026-10-12). The metal a base sweeps in one step is valued along a straight 3D profile, and the moving
# mesh refuses moves of more than about half an element (1.3 mm on 2 mm elements, 5.2 mm on 4 mm, measured): a step
# whose draws and hand-downs would move a base further than this is split.
MOVE_MAX = 0.25                   # of the mean surface edge: the most a base moves in one part (draws and hand-downs
                                  # split a part beyond it; absorption takes what they leave)
SPLIT_DEPTH_MAX = 6               # a part is split at most this many times over
MOVE_HINT_MARGIN = 1.25           # the next macro step's parts: this times the parts this one's largest move called for
# The part length while melting (2026-10-12). On the 100 mm sphere's first 1.5 s of melting, 0.5 s parts cost 1 454 s
# (53 % on parts rejected for their moves, 4.4 passes a kept part) against 341 s for 0.0625 s macro steps (1.8 passes);
# before melting short parts cost about 18 % more, so the cap holds only while some skin's top is non-rigid. (Those
# steps also removed 21.28 g against 16.6 g, but because ABSORB_MAX is per macro step, not because the parts were short.)
PART_MAX = 0.0625                 # s: the longest part while any skin's top is at or above T_rigid, or about to be
SKIN_CG_TOL = 1e-13               # the 3D solver's linear tolerance in a skin run: at its default 1e-10 the skin body's
                                  # books drifted by 1e-8 of the heat absorbed over 17 short parts (about 5e-5 J a solve
                                  # on the 4 mm test sphere); at 1e-13, as with a direct solve, 2e-11, at no cost
                                  # (measured 2026-10-14). The element model's runs keep their solver's own tolerance.


@dataclasses.dataclass
class FlowForcing:
    """The gas side's push on the skins' metal at a macro step, per corner: `dp_ds` [Pa/m, (n,)] the wall pressure's
    gradient along the meridian and `tau` [Pa, (n,)] the gas shear (the surface flow's), and `g` [m/s2, (3,)] the
    acceleration the metal feels in the body's frame -- minus the body's own: the deceleration along the flight
    direction."""
    dp_ds: np.ndarray
    tau: np.ndarray
    g: np.ndarray

    def toward(self, other, frac):
        """This forcing moved the fraction `frac` of the way to `other`, linearly."""
        return FlowForcing(*(a + (b - a) * frac for a, b in ((self.dp_ds, other.dp_ds), (self.tau, other.tau),
                                                             (self.g, other.g))))


def jump_faults(T_bottom, T_node, T_deep, crossed, tol=None):
    """The columns whose skin's bottom cell lies more than `tol` (JUMP_TOLERANCE) outside the temperatures of the metal
    that crossed their base (spec 6.6): between the base node's and the 3D profile's at the depth swept, `T_deep` (equal
    to T_node where nothing was taken from the 3D model). Metal drawn from below the node is colder than the node by the
    profile's slope times the depth, which is not a fault (Asha, 2026-10-12: comparing with the node alone counted 812-1
    694 faults on the level-2 plate at 0.1-0.5 s steps, none at 0.05 s)."""
    tol = skin_mod.JUMP_TOLERANCE if tol is None else tol
    lo, hi = np.minimum(T_node, T_deep), np.maximum(T_node, T_deep)
    return crossed & ((T_bottom < lo - tol) | (T_bottom > hi + tol))


class SkinBody(ThermalBody):
    """A ThermalBody whose surface is resolved by corner columns with layered skins on a moving mesh (module docstring).
    `remove`: liquid leaves as it forms (stage 2); `demise_fraction`: the run ends below this share of the initial mass;
    `flow`: the metal flows (sub-plan 21), `flow_cycles` flow sub-steps in each heat sub-step (None: skin.FLOW_CYCLES),
    driven by `surface_flow` (a `surface_flow.SurfaceFlow`, its defaults if None) from each macro step's aero state."""

    def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                 v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                 repeat_rule=DEFAULT_REPEAT_RULE, part_max=PART_MAX, flow=FLOW_DEFAULT, flow_cycles=None,
                 surface_flow=None):
        if repeat_rule not in REPEAT_RULES:
            raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
        self.flow_cycles = skin_mod.FLOW_CYCLES if flow_cycles is None else flow_cycles
        if not (isinstance(self.flow_cycles, (int, np.integer)) and self.flow_cycles >= 1):
            raise ValueError("flow_cycles must be an integer of at least 1, got {!r}".format(flow_cycles))
        self.flow_on, self.surface_flow = bool(flow), surface_flow
        self.repeat_rule = repeat_rule
        self.part_max = part_max                             # s, or None: no cap on the parts while melting
        self.material, self.solver, self.mass_kg = material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient, self.remove, self.demise_fraction = float(T_ambient), bool(remove), float(demise_fraction)
        self.v_hat = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
        self.settings = settings or skin_mod.SkinSettings()
        if getattr(solver, "cg_tol", None) is not None:
            solver.cg_tol = min(solver.cg_tol, SKIN_CG_TOL)
        solver.setup(mesh, material, 0.0)                    # the 3D boundary lies under the skins: it does not radiate
        solver.set_temperature(T0)
        self.mass0 = material.rho * mesh.volume()
        self.moving = movemesh.MovingMesh(mesh, solver, material, depth=deform_depth)
        self.columns = columns_mod.Columns(mesh.points[self.moving.corners], self.moving.faces_corner)
        self.skins = skin_mod.Skins(material, self.columns.area, self.settings)
        self.pending = np.zeros(self.columns.n)              # J per corner: the deferred interface mismatch
        self.column_removed = np.zeros(self.columns.n)       # kg per column removed as liquid since time zero
        self.integrated_heat = self.absorbed_heat = self.radiated_heat = 0.0
        self.removed_mass = self.removed_enthalpy = 0.0
        self.drawn_mass = self.absorbed_mass = self.handed_mass = 0.0
        self.jump_faults = self.interface_repeats = self.absorb_capped = 0
        self.last_mismatch = self.last_shortfall = 0.0
        self.last_trigger = {"columns_over": 0, "worst_ratio": 0.0, "total_ratio": 0.0}    # the first pass's mismatch
        self.last_passes, self.last_total = 1, 0.0           # the damped rule's passes and final total mismatch ratio
        self.step_parts, self._hint, self._ratio = 1, 1, 0.0   # parts of the last macro step; of the next; its move ratio
        self._q_prev = None                                  # the last macro step's heating, for the interpolation in time
        self._top_start = None                               # the tops' temperatures [K] at the last macro step's start
        self._dt_last = None                                 # s, and its length: the look-ahead of the part cap
        self.capped_ahead = 0                                # macro steps capped because a top was about to be non-rigid
        self._heating = None                                 # this macro step's (t0, dt, q at its start, q at its end)
        self.wall = {"skins": 0.0, "3d": 0.0, "move": 0.0, "rebuild": 0.0, "flow": 0.0}
        self.last_flow = None                                # the surface flow's last evaluation (SurfaceFlowResult)
        self._forcing = self._forcing_prev = None            # this macro step's (start, end) FlowForcing; the last end
        self.impulse, self.momentum_lost = np.zeros(3), np.zeros(3)   # the momentum books since time zero [N s, kg m/s]
        self.impulse_size = 0.0                              # N s: the push's, shear's and walls' impulses, in size
        self.flow_courant = self.flow_growth = 0.0           # this macro step's worst flow sub-step
        self.flow_passes = 0
        self._volume = None                                  # m3, the 3D model's volume since the last move (mass)
        self.iterations, self.last, self.last_move = [], None, None
        self._fill(T0)
        self.energy0 = self.energy()
        self.momentum0 = self.skins.momentum()
        self._top0, self._normal0 = self.surface.points.copy(), self.surface.normals.copy()

    # -- set-up -----------------------------------------------------------------------------------------------------
    def _fill(self, T0):
        """Each column takes its minimum mass from the 3D model (spec 6.1): the corners move in by it, and the skins hold
        exactly what they swept, at the initial temperature."""
        m_min = skin_mod.minimum_mass(self.material, self.columns.area, self.settings)
        u = -self.columns.normal * (m_min / (self.material.rho * self.columns.area))[:, None]
        r = self.moving.move_with_policy(u)
        self.columns.update(self.moving.mesh.points[self.moving.corners])
        self.skins.area = self.columns.area.copy()
        self.skins.fill(self.material.rho * r.swept_volume, T0)
        self.skins.exchange_bottom(np.zeros(self.columns.n), r.swept_enthalpy - self.skins.H.sum(axis=1))
        self._refresh_geometry()                       # the top surface: the bases moved out by the filled skins

    def _refresh_geometry(self):
        """The columns, their areas and the top surface after the corners moved; the 3D model's volume is measured again
        when next asked for (`mass`)."""
        mm = self.moving
        self._volume = None
        self.columns.update(mm.mesh.points[mm.corners])
        self.skins.area = self.columns.area.copy()
        self.surface = columns_mod.CornerSurface.from_points(self.columns.top_points(self.skins.thickness()),
                                                             self.columns.faces)
        self.theta = self.surface.angles_to(self.v_hat)
        self.i_stag, self.i_back = self.surface.patch_toward(self.v_hat), self.surface.patch_toward(-self.v_hat)

    @property
    def mesh(self):
        return self.moving.mesh

    # -- the macro step ---------------------------------------------------------------------------------------------
    def advance(self, t, dt, loads, state=None, forcing=None):
        """The macro step from t - dt to t (module docstring) as equal parts, as many as the last macro step's largest
        move called for. With the flow on, the gas side's push is the surface flow's on the corners from the aero state
        `state` (`flow_forcing`), or `forcing`, a FlowForcing given instead (a prescribed push, in the tests)."""
        n, self._ratio, self.step_parts = self._hint, 0.0, 0
        if self.flow_on:
            end = forcing if forcing is not None else self.flow_forcing(state)
            last = self._forcing_prev
            start = last if last is not None and last.dp_ds.shape == end.dp_ds.shape else end
            self._forcing, self._forcing_prev = (start, end), end
            self.flow_courant = self.flow_growth = 0.0
            self.flow_passes = 0
        sk, mat = self.skins, self.material
        top = np.where(sk.valid[:, 0], sk.T[:, 0], np.nan)
        if self.part_max is not None and mat.melts:
            if np.any(sk.valid[:, 0] & (sk.T[:, 0] >= mat.T_rigid)):
                n = max(n, int(math.ceil(dt / self.part_max - 1e-9)))
            elif self._top_start is not None and self._top_start.shape == top.shape:
                with np.errstate(invalid="ignore"):          # a top that rose over the last step, carried on over this one
                    ahead = top + (top - self._top_start) * (dt / self._dt_last)
                if np.any(np.nan_to_num(ahead, nan=-np.inf) >= mat.T_rigid):
                    n = max(n, int(math.ceil(dt / self.part_max - 1e-9)))
                    self.capped_ahead += 1
        self._top_start, self._dt_last = top, dt
        q_end = np.asarray(loads.q_conv, dtype=float)
        q_start = self._q_prev if self._q_prev is not None and self._q_prev.shape == q_end.shape else q_end
        self._heating = (t - dt, dt, q_start, q_end)
        for k in range(n):
            self._advance_part(t if k == n - 1 else t - dt + (k + 1) * (dt / n), dt / n, dt, loads, 0 if n == 1 else 1)
        self._hint = max(1, int(math.ceil(MOVE_HINT_MARGIN * self._ratio)))
        self._q_prev = q_end.copy()

    def flow_forcing(self, state):
        """The gas side's push on the metal at this macro step (FlowForcing): the surface flow evaluated on the corners
        from the aero state, as the heating is -- the top surface's angles to the flight direction, the nose radius and
        the tops' temperatures at the step's start -- and the deceleration along the flight direction."""
        if state is None:
            raise ValueError("a skin body with the flow on needs the macro step's aero state, or a forcing given")
        if self.surface_flow is None:
            from . import surface_flow
            self.surface_flow = surface_flow.SurfaceFlow()
        mat = self.material
        rho_l = mat.liquid.rho if mat.liquid is not None else mat.rho
        res = self.surface_flow.evaluate(state, self.theta, self.nose_radius(), rho_l, self.surface_temperature())
        self.last_flow = res
        return FlowForcing(np.asarray(res.dp_ds, dtype=float), np.asarray(res.tau, dtype=float),
                           res.deceleration * self.v_hat)

    def _advance_part(self, t, dt, dt_macro, loads, depth):
        """One part of a macro step, ending at t: the passes, then -- if its draws and hand-downs would move a base more
        than MOVE_MAX of the mean surface edge -- the part again as shorter parts, else the books and the move."""
        mat, sk, mm = self.material, self.skins, self.moving
        t0, dt_m, q_start, q_end = self._heating
        frac = ((t - 0.5 * dt) - t0) / dt_m
        q = q_start + (q_end - q_start) * frac                              # the heating at the part's middle
        A_top = self.surface.areas
        flow = None
        if self.flow_on:                                                    # the gas side's push at the part's middle
            f = self._forcing[0].toward(self._forcing[1], frac)
            sigma = mat.liquid.sigma if mat.liquid is not None else 0.0
            flow = skin_mod.Flow(self.columns, f.dp_ds, skin_mod.along_meridian(self.columns.normal, self.v_hat, f.tau),
                                 f.g, self.v_hat, sk.rho_l, sigma, A_top, self.flow_cycles)
        eps = np.full(sk.n, self.emissivity)
        m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
        side = lambda s, previous: skin_mod.Sideways(s, self.columns, previous)   # implicit (sub-plan 21)
        nodes = mm.corners
        T_b0 = self.solver.temperature()[nodes]
        g0 = self._normal_gradient()                         # the profile's slope below every base, for this part
        load = np.zeros(mm.mesh.n_nodes)
        load[nodes] = self.pending / dt
        solver_state, skin_state = self.solver.state(), sk.state()
        T_lin = T_b0
        damped = self.repeat_rule == "damped"
        omega, r_prev = DAMPED_OMEGA0, None
        flow_seconds = sum(skin_mod.FLOW_SECONDS.values())
        for attempt in range(DAMPED_MAX_PASSES if damped else 2):
            started = time.perf_counter()
            a, b = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_lin, m_min, side=side, remove=self.remove,
                                    g_base=g0, flow=flow)
            self.wall["skins"] += time.perf_counter() - started
            started = time.perf_counter()
            res = self.solver.step(dt, np.zeros(len(self.solver.areas)), self.T_ambient, nodal_load=load,
                                   nodal_law=(nodes, a, b))
            self.wall["3d"] += time.perf_counter() - started
            T_b1 = self.solver.temperature()[nodes]
            started = time.perf_counter()
            book = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_b0, m_min, T_b1=T_b1, side=side,
                                    remove=self.remove, g_base=g0, flow=flow)
            self.wall["skins"] += time.perf_counter() - started
            mismatch = book.E_base - (a + b * T_b1) * dt
            scale = np.abs(book.E_top) + np.abs(book.E_base) + 1e-30
            ratio = np.abs(mismatch) / scale
            total = float(np.abs(mismatch).sum() / scale.sum())
            if attempt == 0:
                self.last_trigger = {"columns_over": int((ratio > skin_mod.INTERFACE_TOLERANCE).sum()),
                                     "worst_ratio": float(ratio.max()), "total_ratio": total}
            if damped:
                r = T_b1 - T_lin
                self.last_passes, self.last_total = attempt + 1, total
                if (total <= skin_mod.INTERFACE_TOLERANCE or np.abs(r).max() <= DAMPED_TOL_K
                        or attempt == DAMPED_MAX_PASSES - 1):
                    break
                if r_prev is not None:                       # Aitken: omega from the change of the residual
                    dr = r - r_prev
                    denom = float(dr @ dr)
                    if denom > 0.0:
                        omega = float(np.clip(-omega * float(r_prev @ dr) / denom, *DAMPED_OMEGA_RANGE))
                self.solver.set_state(solver_state)
                sk.set_state(skin_state)
                T_lin, r_prev = T_lin + omega * r, r
                self.interface_repeats += 1
                continue
            repeat = {"column": bool((ratio > skin_mod.INTERFACE_TOLERANCE).any()),
                      "total": total > skin_mod.INTERFACE_TOLERANCE, "never": False}[self.repeat_rule]
            if attempt == 0 and repeat:
                self.solver.set_state(solver_state)
                sk.set_state(skin_state)
                T_lin = T_b1
                self.interface_repeats += 1
                continue
            break
        self.wall["flow"] += sum(skin_mod.FLOW_SECONDS.values()) - flow_seconds
        ratio = self._move_ratio(book, m_min)
        if ratio > 1.0 and depth < SPLIT_DEPTH_MAX:          # too far for one part: run it again as shorter parts
            self.solver.set_state(solver_state)
            sk.set_state(skin_state)
            n = int(math.ceil(ratio))
            for k in range(n):
                self._advance_part(t if k == n - 1 else t - dt + (k + 1) * (dt / n), dt / n, dt_macro, loads, depth + 1)
            return
        self._ratio = max(self._ratio, ratio * dt_macro / dt)
        self.step_parts += 1
        self.pending = mismatch
        self.last_mismatch = float(np.abs(mismatch).sum())
        E_top, E_rad = float(book.E_top.sum()), float(book.E_rad.sum())
        self.integrated_heat += E_top
        self.absorbed_heat += E_top - E_rad
        self.radiated_heat += E_rad
        self.removed_mass += float(book.removed.sum())
        self.column_removed += book.removed
        self.removed_enthalpy += float(book.removed_H.sum())
        self.drawn_mass += float(book.drawn.sum())
        self.momentum_lost += book.momentum_out.sum(axis=0)
        if book.flow is not None:                            # the flow's impulses (the transport's carried sum to zero)
            fb = book.flow
            self.impulse += (fb.push + fb.shear + fb.wall + fb.carried).sum(axis=0)
            self.impulse_size += float(sum(np.linalg.norm(x, axis=1).sum() for x in (fb.push, fb.shear, fb.wall)))
            self.flow_courant, self.flow_growth = max(self.flow_courant, fb.courant), max(self.flow_growth, fb.growth)
            self.flow_passes = max(self.flow_passes, fb.passes)
        self.iterations.append(res.iterations)
        self.last = dataclasses.replace(res, Q_conv=E_top / dt, Q_rad=E_rad / dt)    # the gas-side powers, for the history
        self._move_base(book, m_min, g0)
        self.last_shortfall = float(np.maximum(skin_mod.minimum_mass(mat, sk.area, self.settings) - sk.column_mass(),
                                               0.0).sum())

    def _absorb_depth(self, allowance, count=True):
        """Per corner, the depth [m] along the column's direction from the base node to where the 3D temperature falls
        to T_rigid (spec 6.6), from the node's temperature and the 3D gradient along the normal, at most `allowance` [m,
        per corner or one for all]: what the part's draws and hand-downs left of its move. `count`: add the corners the
        allowance held back to absorb_capped. (Against the 3D field traced along the normal, the line method, the
        extrapolation agrees within 0.07 mm on 2 mm elements, measured 2026-10-14.)"""
        mat, mm = self.material, self.moving
        if not mat.melts:
            return np.zeros(self.columns.n)
        T = self.solver.temperature()
        excess = T[mm.corners] - mat.T_rigid
        if not np.any(excess > 0.0):
            return np.zeros(self.columns.n)
        g = self._normal_gradient(T)
        depth = np.where(excess > 0.0, excess / np.maximum(g, GRADIENT_FLOOR), 0.0)
        if count:
            self.absorb_capped += int((depth > allowance).sum())
        return np.minimum(depth, allowance)

    def _normal_gradient(self, T=None):
        """dT/dn [K/m] at every corner along its column's outward direction, from the 3D field's nodal gradient (exact
        for a linear field; positive when the surface is hotter than beneath it)."""
        mm = self.moving
        T = self.solver.temperature() if T is None else T
        mesh = mm.mesh
        touching = np.zeros(mesh.n_nodes, dtype=bool)
        touching[mm.corners] = True
        elements = np.flatnonzero(touching[mesh.tets].any(axis=1))
        return np.einsum("ij,ij->i", movemesh.nodal_gradient(mesh.points, mesh.tets, T, elements)[mm.corners],
                         self.columns.normal)

    def _move_ratio(self, book, m_min):
        """The largest distance a base would move in this part by drawing and handing down, over MOVE_MAX of the mean
        surface edge (absorption takes what is left, so it never splits a part). The hand-down is found on a copy of
        the skins."""
        mat, sk = self.material, self.skins
        handed = np.zeros(sk.n)
        if mat.melts:
            saved = sk.state()
            handed, _ = skin_mod.hand_down(sk, m_min, mat.T_rigid)
            sk.set_state(saved)
        d = np.abs(book.drawn - handed) / (mat.rho * self.columns.area)
        return float(d.max() / (MOVE_MAX * self.columns.edge_length.mean()))

    def _move_base(self, book, m_min, g0):
        """Draw (already in the skins), absorb and hand down become one corner displacement; the move's actual sweep
        is then reconciled on every column's bottom (module docstring, step 4). Absorption takes what the draws and
        hand-downs leave of MOVE_MAX of the mean edge. Metal swept out of the 3D model is booked along the profile
        T(s) = T_node - g0 s that the draws used."""
        mat, sk, mm, col = self.material, self.skins, self.moving, self.columns
        started = time.perf_counter()
        T_node = self.solver.temperature()[mm.corners]
        handed, handed_H = skin_mod.hand_down(sk, m_min, mat.T_rigid) if mat.melts else (np.zeros(sk.n), np.zeros(sk.n))
        if mat.melts:
            self.momentum_lost += sk.momentum_out.sum(axis=0)                     # what the handed-down cells held
        moved = (book.drawn - handed) / (mat.rho * col.area)                      # m inward by drawing and handing down
        depth = self._absorb_depth(np.maximum(MOVE_MAX * col.edge_length.mean() - moved, 0.0))
        inward = book.drawn / mat.rho + depth * col.area - handed / mat.rho          # m3 per corner
        u = -col.normal * (inward / col.area)[:, None] + mm.neutral_smoothing(col.normal, col.area)
        rebuilds = mm.rebuilds
        h_inward = lambda d: skin_mod.profile_enthalpy(mat, T_node, g0, 0.0, d)
        r = mm.move_with_policy(u, h_inward)
        if mm.rebuilds > rebuilds:
            self.wall["rebuild"] += mm.last_rebuild.seconds
        crossed = (book.drawn > 0.0) | (depth > 0.0) | (handed > 0.0)
        T_bottom = sk.T[np.arange(sk.n), np.maximum(sk.count - 1, 0)]
        T_deep = T_node - g0 * (book.drawn / (mat.rho * col.area) + depth)   # the profile at the depth swept
        self.jump_faults += int(jump_faults(T_bottom, T_node, T_deep, crossed).sum())
        sk.exchange_bottom(mat.rho * r.swept_volume - (book.drawn - handed), r.swept_enthalpy - (book.drawn_H - handed_H))
        self.momentum_lost += sk.momentum_out.sum(axis=0)                         # what a reconciliation took away
        sk.regrid()
        self.absorbed_mass += float((mat.rho * depth * col.area).sum())
        self.handed_mass += float(handed.sum())
        self.last_move = r
        self._refresh_geometry()
        self.wall["move"] += time.perf_counter() - started

    # -- what the coupled run reads ---------------------------------------------------------------------------------
    def mass(self, t=None):
        """The 3D model's metal and the skins'. The 3D volume changes only when the corners move, so it is measured once
        per move (`_refresh_geometry`), not at every call: the trajectory asks at every evaluation of its right-hand
        side, and on the 1.2 million elements of the 100 mm sphere's fine band those measurements took 12 % of the time
        of a flight's first six steps (Asha, 2026-10-17)."""
        if self._volume is None:
            self._volume = self.mesh.volume()
        return self.material.rho * self._volume + self.skins.mass()

    def energy(self):
        """Enthalpy held by the 3D model and the skins, and the deferred mismatch in transit to the 3D nodes [J]."""
        return self.solver.energy() + self.skins.energy() + float(self.pending.sum())

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat + removed enthalpy) / absorbed heat: zero for exact books."""
        if not self.absorbed_heat:
            return 0.0
        return (self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy) / self.absorbed_heat

    def momentum_residual(self):
        """|P - P0 - impulses + momentum lost| over the impulses' total size: zero for exact books (spec 10.2)."""
        if self.impulse_size <= 0.0:
            return 0.0
        miss = self.skins.momentum() - self.momentum0 - self.impulse + self.momentum_lost
        return float(np.linalg.norm(miss) / self.impulse_size)

    def mean_temperature(self):
        return float(self.material.temperature_from_enthalpy(self.energy() / self.mass()))

    def surface_temperature(self):
        """The top cell of every column: what the gas sees (spec 6.3)."""
        return self.skins.T[:, 0].copy()

    def radiated_power(self):
        return float((self.emissivity * SIGMA_SB * (self.skins.T[:, 0] ** 4 - self.T_ambient ** 4) * self.surface.areas).sum())

    def field(self):
        return self.solver.temperature()

    def nose_radius(self):
        return float(self.mesh.params.get("radius_m", 0.05))

    def demised(self):
        return self.mass() < self.demise_fraction * self.mass0

    def surface_stats(self):
        T, Tt = self.solver.temperature(), self.surface_temperature()
        centre = self.mesh.points.mean(axis=0)
        i_centre = int(np.argmin(np.linalg.norm(self.mesh.points - centre, axis=1)))
        A = self.surface.areas
        return {"surface_T_max_K": float(Tt.max()), "surface_T_min_K": float(Tt.min()),
                "surface_T_mean_K": float((Tt * A).sum() / A.sum()), "T_stagnation_K": float(Tt[self.i_stag]),
                "T_back_K": float(Tt[self.i_back]), "T_centre_K": float(T[i_centre])}

    def skin_stats(self):
        """The history columns of spec 10.2 that exist at stage 2."""
        sk, mm, mat = self.skins, self.moving, self.material
        th = sk.thickness()
        m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
        flowing = sk.depth_above(mat.T_rigid) if mat.melts else np.zeros(sk.n)
        liquid = sk.depth_above(mat.T_liquidus) if mat.melts else np.zeros(sk.n)
        r = self.last_move
        i, j = self.columns.edges[:, 0], self.columns.edges[:, 1]
        return {"skin_cells": int(sk.count.sum()), "skin_thickness_max_mm": float(th.max() * 1e3),
                "skin_above_min": int((sk.column_mass() > m_min * (1.0 + 1e-9)).sum()),
                "flowing_depth_mean_mm": float(flowing.mean() * 1e3), "flowing_depth_max_mm": float(flowing.max() * 1e3),
                "liquid_depth_mean_mm": float(liquid.mean() * 1e3), "liquid_depth_max_mm": float(liquid.max() * 1e3),
                "surface_quality_min": mm.surface_quality(), "surface_fold_deg": mm.surface_fold(),
                "zone_ratio_min": float(r.worst_ratio) if r is not None else 1.0, "rebuilds": mm.rebuilds,
                "rebuild_residual_J": mm.rebuild_residual_total,
                "rebuild_residual_last_J": mm.last_rebuild.residual if mm.last_rebuild is not None else 0.0,
                "drawn_mass_kg": self.drawn_mass, "absorbed_mass_kg": self.absorbed_mass, "handed_mass_kg": self.handed_mass,
                "shortfall_kg": self.last_shortfall, "jump_faults": self.jump_faults, "absorb_capped": self.absorb_capped,
                "interface_mismatch_J": self.last_mismatch, "interface_repeats": self.interface_repeats,
                "repeat_columns_over": self.last_trigger["columns_over"],
                "mismatch_worst_ratio": self.last_trigger["worst_ratio"],
                "mismatch_total_ratio": self.last_trigger["total_ratio"],
                "removed_mass_kg": self.removed_mass, "removed_enthalpy_J": self.removed_enthalpy,
                "base_recession_max_mm": float(self._recession().max() * 1e3),
                "wall_skins_s": self.wall["skins"], "wall_3d_s": self.wall["3d"], "wall_move_s": self.wall["move"],
                "wall_rebuild_s": self.wall["rebuild"], "step_parts": self.step_parts,
                # the flow (sub-plan 21 Task 6c): its momentum books, the fastest cell, the top's roughness -- the largest
                # height difference between neighbours over their spacing -- and the flow sub-steps' worst Courant
                # number, Rayleigh-Taylor growth and linearisations in this macro step
                "momentum_residual": self.momentum_residual(),
                "speed_max_ms": float(np.linalg.norm(sk.U, axis=2).max()) if sk.U.size else 0.0,
                "roughness_max": float((np.abs(th[i] - th[j]) / self.columns.edge_length).max()),
                "flow_courant_max": self.flow_courant, "flow_growth_max": self.flow_growth,
                "flow_passes_max": self.flow_passes, "wall_flow_s": self.wall["flow"]}

    def _recession(self):
        """How far each corner's top has moved inward along its starting direction since time zero [m]."""
        return np.einsum("ij,ij->i", self._top0 - self.surface.points, self._normal0)
