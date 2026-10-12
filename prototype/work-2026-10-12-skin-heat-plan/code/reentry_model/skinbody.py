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
   the difference between what it sent through the base and what the 3D step received is deferred to the node, and a
   step whose mismatch exceeds INTERFACE_TOLERANCE of the column's throughput repeats once, re-linearised;
4. the base moves by the metal crossing it (spec 6.6): drawn (in the sub-steps) and absorbed metal move the corner in,
   handed-down metal moves it out; with the volume-neutral smoothing this is one displacement per corner, which
   `MovingMesh.move_with_policy` applies (rebuilding when needed), and each column is then reconciled with what its
   corner actually swept: its bottom gains the swept mass less what it already exchanged, and its enthalpy changes by
   the swept enthalpy less what it already exchanged, so mass and energy are exact.

The body's surface for the coupled run is the top surface, one 'patch' per corner (`columns.CornerSurface`), so the
heating, the history and the trajectory's loop are unchanged."""
import dataclasses
import time

import numpy as np

from . import columns as columns_mod
from . import movemesh
from . import skin as skin_mod
from .body import ThermalBody
from .thermal import SIGMA_SB

ABSORB_MAX = 0.25                 # of the mean surface edge: the most a base absorbs in one macro step (the rest waits)
GRADIENT_FLOOR = 1.0              # K/m: below this the 3D field gives no direction to absorb along
REPEAT_RULES = ("column", "total", "never")   # when a step's interface mismatch repeats it (spec 7.2: "column"; the others
                                              # are measured in sub-plan 20 Task 10 before any change of default)


class SkinBody(ThermalBody):
    """A ThermalBody whose surface is resolved by corner columns with layered skins on a moving mesh (module docstring).
    `remove`: liquid leaves as it forms (stage 2); `demise_fraction`: the run ends below this share of the initial mass."""

    def __init__(self, mesh, material, solver, mass_kg, settings=None, T0=300.0, emissivity=None, T_ambient=0.0,
                 v_hat=(1.0, 0.0, 0.0), deform_depth=movemesh.DEFORM_DEPTH, remove=True, demise_fraction=0.01,
                 repeat_rule="column"):
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
        mat, sk, mm = self.material, self.skins, self.moving
        q = np.asarray(loads.q_conv, dtype=float)
        A_top = self.surface.areas
        eps = np.full(sk.n, self.emissivity)
        m_min = skin_mod.minimum_mass(mat, sk.area, self.settings)
        side = lambda s, h: skin_mod.sideways_heat(s, self.columns, h)
        nodes = mm.corners
        T_b0 = self.solver.temperature()[nodes]
        load = np.zeros(mm.mesh.n_nodes)
        load[nodes] = self.pending / dt
        solver_state, skin_state = self.solver.state(), sk.state()
        T_lin = T_b0
        for attempt in (0, 1):
            started = time.perf_counter()
            a, b = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_lin, m_min, side=side, remove=self.remove)
            self.wall["skins"] += time.perf_counter() - started
            started = time.perf_counter()
            res = self.solver.step(dt, np.zeros(len(self.solver.areas)), self.T_ambient, nodal_load=load,
                                   nodal_law=(nodes, a, b))
            self.wall["3d"] += time.perf_counter() - started
            T_b1 = self.solver.temperature()[nodes]
            started = time.perf_counter()
            book = skin_mod.advance(sk, dt, q, A_top, self.T_ambient, eps, T_b0, m_min, T_b1=T_b1, side=side,
                                    remove=self.remove)
            self.wall["skins"] += time.perf_counter() - started
            mismatch = book.E_base - (a + b * T_b1) * dt
            scale = np.abs(book.E_top) + np.abs(book.E_base) + 1e-30
            ratio = np.abs(mismatch) / scale
            total = float(np.abs(mismatch).sum() / scale.sum())
            if attempt == 0:
                self.last_trigger = {"columns_over": int((ratio > skin_mod.INTERFACE_TOLERANCE).sum()),
                                     "worst_ratio": float(ratio.max()), "total_ratio": total}
            repeat = {"column": bool((ratio > skin_mod.INTERFACE_TOLERANCE).any()),
                      "total": total > skin_mod.INTERFACE_TOLERANCE, "never": False}[self.repeat_rule]
            if attempt == 0 and repeat:
                self.solver.set_state(solver_state)
                sk.set_state(skin_state)
                T_lin = T_b1
                self.interface_repeats += 1
                continue
            break
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
        self._move_base(book, m_min)
        self.last_shortfall = float(np.maximum(skin_mod.minimum_mass(mat, sk.area, self.settings) - sk.column_mass(),
                                               0.0).sum())

    def _absorb_depth(self):
        """Per corner, the depth [m] along the column's direction from the base node to where the 3D temperature falls
        to T_rigid (spec 6.6), from the node's temperature and the 3D gradient along the normal, capped at ABSORB_MAX
        of the mean surface edge per macro step."""
        mat, mm = self.material, self.moving
        if not mat.melts:
            return np.zeros(self.columns.n)
        T = self.solver.temperature()
        excess = T[mm.corners] - mat.T_rigid
        if not np.any(excess > 0.0):
            return np.zeros(self.columns.n)
        mesh = mm.mesh
        touching = np.zeros(mesh.n_nodes, dtype=bool)
        touching[mm.corners] = True
        elements = np.flatnonzero(touching[mesh.tets].any(axis=1))
        g = np.einsum("ij,ij->i", movemesh.nodal_gradient(mesh.points, mesh.tets, T, elements)[mm.corners],
                      self.columns.normal)
        depth = np.where(excess > 0.0, excess / np.maximum(g, GRADIENT_FLOOR), 0.0)
        cap = ABSORB_MAX * self.columns.edge_length.mean()
        self.absorb_capped += int((depth > cap).sum())
        return np.minimum(depth, cap)

    def _move_base(self, book, m_min):
        """Draw (already in the skins), absorb and hand down become one corner displacement; the move's actual sweep
        is then reconciled on every column's bottom (module docstring, step 4)."""
        mat, sk, mm, col = self.material, self.skins, self.moving, self.columns
        started = time.perf_counter()
        T_node = self.solver.temperature()[mm.corners]
        depth = self._absorb_depth()
        handed, handed_H = skin_mod.hand_down(sk, m_min, mat.T_rigid) if mat.melts else (np.zeros(sk.n), np.zeros(sk.n))
        inward = book.drawn / mat.rho + depth * col.area - handed / mat.rho          # m3 per corner
        u = -col.normal * (inward / col.area)[:, None] + mm.neutral_smoothing(col.normal, col.area)
        rebuilds = mm.rebuilds
        r = mm.move_with_policy(u)
        if mm.rebuilds > rebuilds:
            self.wall["rebuild"] += mm.last_rebuild.seconds
        crossed = (book.drawn > 0.0) | (depth > 0.0) | (handed > 0.0)
        T_bottom = sk.T[np.arange(sk.n), np.maximum(sk.count - 1, 0)]
        self.jump_faults += int((crossed & (np.abs(T_node - T_bottom) > skin_mod.JUMP_TOLERANCE)).sum())
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
                "wall_rebuild_s": self.wall["rebuild"]}

    def _recession(self):
        """How far each corner's top has moved inward along its starting direction since time zero [m]."""
        return np.einsum("ij,ij->i", self._top0 - self.surface.points, self._normal0)
