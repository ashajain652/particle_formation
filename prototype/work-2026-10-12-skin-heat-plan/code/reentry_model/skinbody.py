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

A macro step is run as equal parts so that no base moves more than MOVE_MAX of the mean surface edge in one part by
drawing, absorbing and handing down (Asha, 2026-10-12): a part that would is run again as shorter parts (absorption's
cap, ABSORB_MAX per macro step, is shared among the parts by their length), and the next macro step
starts with as many parts as this one's moves called for. Within a part the base heat, the draws and the move all
take the part's own step; the heating is the macro step's.
The body's surface for the coupled run is the top surface, one 'patch' per corner (`columns.CornerSurface`), so the
heating, the history and the trajectory's loop are unchanged."""
import dataclasses
import math
import time

import numpy as np

from . import columns as columns_mod
from . import movemesh
from . import skin as skin_mod
from .body import ThermalBody
from .thermal import SIGMA_SB

ABSORB_MAX = 0.25                 # of the mean surface edge: the most a base absorbs in one macro step (the rest waits)
GRADIENT_FLOOR = 1.0              # K/m: the smallest gradient along the normal that is divided by. A corner above T_rigid
                                  # whose gradient is below it (a flat or reversed field) is given this one, so its depth
                                  # is huge and it absorbs up to the cap, ABSORB_MAX: it does not absorb nothing
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
MOVE_MAX = 0.25                   # of the mean surface edge: the most a base moves in one part (draw, absorb, hand down)
SPLIT_DEPTH_MAX = 6               # a part is split at most this many times over
MOVE_HINT_MARGIN = 1.25           # the next macro step's parts: this times the parts this one's largest move called for


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
    `remove`: liquid leaves as it forms (stage 2); `demise_fraction`: the run ends below this share of the initial mass."""

    def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                 v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                 repeat_rule=DEFAULT_REPEAT_RULE):
        if repeat_rule not in REPEAT_RULES:
            raise ValueError("repeat_rule must be one of {}, got {!r}".format(REPEAT_RULES, repeat_rule))
        self.repeat_rule = repeat_rule
        self.material, self.solver, self.mass_kg = material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient, self.remove, self.demise_fraction = float(T_ambient), bool(remove), float(demise_fraction)
        self.v_hat = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
        self.settings = settings or skin_mod.SkinSettings()
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
        self.wall = {"skins": 0.0, "3d": 0.0, "move": 0.0, "rebuild": 0.0}
        self.iterations, self.last, self.last_move = [], None, None
        self._fill(T0)
        self.energy0 = self.energy()
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
        """The columns, their areas and the top surface after the corners moved."""
        mm = self.moving
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
    def advance(self, t, dt, loads, state=None):
        """The macro step from t - dt to t (module docstring) as equal parts, as many as the last macro step's largest
        move called for."""
        n, self._ratio, self.step_parts = self._hint, 0.0, 0
        for k in range(n):
            self._advance_part(t if k == n - 1 else t - dt + (k + 1) * (dt / n), dt / n, dt, loads, 0 if n == 1 else 1)
        self._hint = max(1, int(math.ceil(MOVE_HINT_MARGIN * self._ratio)))

    def _advance_part(self, t, dt, dt_macro, loads, depth):
        """One part of a macro step, ending at t: the passes, then -- if its draws and hand-downs would move a base more
        than MOVE_MAX of the mean surface edge -- the part again as shorter parts, else the books and the move."""
        mat, sk, mm = self.material, self.skins, self.moving
        q = np.asarray(loads.q_conv, dtype=float)
        A_top = self.surface.areas
        eps = np.full(sk.n, self.emissivity)
        m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
        side = lambda s, h: skin_mod.sideways_heat(s, self.columns, h)
        nodes = mm.corners
        T_b0 = self.solver.temperature()[nodes]
        g0 = self._normal_gradient()                         # the profile's slope below every base, for this part
        load = np.zeros(mm.mesh.n_nodes)
        load[nodes] = self.pending / dt
        solver_state, skin_state = self.solver.state(), sk.state()
        T_lin = T_b0
        damped = self.repeat_rule == "damped"
        omega, r_prev = DAMPED_OMEGA0, None
        for attempt in range(DAMPED_MAX_PASSES if damped else 2):
            started = time.perf_counter()
            a, b = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_lin, m_min, side=side, remove=self.remove,
                                    g_base=g0)
            self.wall["skins"] += time.perf_counter() - started
            started = time.perf_counter()
            res = self.solver.step(dt, np.zeros(len(self.solver.areas)), self.T_ambient, nodal_load=load,
                                   nodal_law=(nodes, a, b))
            self.wall["3d"] += time.perf_counter() - started
            T_b1 = self.solver.temperature()[nodes]
            started = time.perf_counter()
            book = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_b0, m_min, T_b1=T_b1, side=side,
                                    remove=self.remove, g_base=g0)
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
        share = dt / dt_macro                                # of the macro step: absorption's cap is shared by length
        ratio = self._move_ratio(book, m_min, share)
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
        self.iterations.append(res.iterations)
        self.last = dataclasses.replace(res, Q_conv=E_top / dt, Q_rad=E_rad / dt)    # the gas-side powers, for the history
        self._move_base(book, m_min, g0, share)
        self.last_shortfall = float(np.maximum(skin_mod.minimum_mass(mat, sk.area, self.settings) - sk.column_mass(),
                                               0.0).sum())

    def _absorb_depth(self, share=1.0, count=True):
        """Per corner, the depth [m] along the column's direction from the base node to where the 3D temperature falls
        to T_rigid (spec 6.6), from the node's temperature and the 3D gradient along the normal, capped at ABSORB_MAX
        of the mean surface edge per macro step: `share` of it in a part of that length. `count`: add the capped
        corners to absorb_capped."""
        mat, mm = self.material, self.moving
        if not mat.melts:
            return np.zeros(self.columns.n)
        T = self.solver.temperature()
        excess = T[mm.corners] - mat.T_rigid
        if not np.any(excess > 0.0):
            return np.zeros(self.columns.n)
        g = self._normal_gradient(T)
        depth = np.where(excess > 0.0, excess / np.maximum(g, GRADIENT_FLOOR), 0.0)
        cap = ABSORB_MAX * self.columns.edge_length.mean() * share
        if count:
            self.absorb_capped += int((depth > cap).sum())
        return np.minimum(depth, cap)

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

    def _move_ratio(self, book, m_min, share):
        """The largest distance a base would move in this part by drawing, absorbing and handing down, over MOVE_MAX of
        the mean surface edge. The hand-down is found on a copy of the skins."""
        mat, sk = self.material, self.skins
        handed = np.zeros(sk.n)
        if mat.melts:
            saved = sk.state()
            handed, _ = skin_mod.hand_down(sk, m_min, mat.T_rigid)
            sk.set_state(saved)
        d = np.abs((book.drawn - handed) / (mat.rho * self.columns.area) + self._absorb_depth(share, count=False))
        return float(d.max() / (MOVE_MAX * self.columns.edge_length.mean()))

    def _move_base(self, book, m_min, g0, share=1.0):
        """Draw (already in the skins), absorb and hand down become one corner displacement; the move's actual sweep
        is then reconciled on every column's bottom (module docstring, step 4). Metal swept out of the 3D model is
        booked along the profile T(s) = T_node - g0 s that the draws used."""
        mat, sk, mm, col = self.material, self.skins, self.moving, self.columns
        started = time.perf_counter()
        T_node = self.solver.temperature()[mm.corners]
        depth = self._absorb_depth(share)
        handed, handed_H = skin_mod.hand_down(sk, m_min, mat.T_rigid) if mat.melts else (np.zeros(sk.n), np.zeros(sk.n))
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
        sk.regrid()
        self.absorbed_mass += float((mat.rho * depth * col.area).sum())
        self.handed_mass += float(handed.sum())
        self.last_move = r
        self._refresh_geometry()
        self.wall["move"] += time.perf_counter() - started

    # -- what the coupled run reads ---------------------------------------------------------------------------------
    def mass(self, t=None):
        return self.material.rho * self.mesh.volume() + self.skins.mass()

    def energy(self):
        """Enthalpy held by the 3D model and the skins, and the deferred mismatch in transit to the 3D nodes [J]."""
        return self.solver.energy() + self.skins.energy() + float(self.pending.sum())

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat + removed enthalpy) / absorbed heat: zero for exact books."""
        if not self.absorbed_heat:
            return 0.0
        return (self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy) / self.absorbed_heat

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
                "wall_rebuild_s": self.wall["rebuild"], "step_parts": self.step_parts}

    def _recession(self):
        """How far each corner's top has moved inward along its starting direction since time zero [m]."""
        return np.einsum("ij,ij->i", self._top0 - self.surface.points, self._normal0)
