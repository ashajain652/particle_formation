"""Lockstep coupling of the trajectory stepper with the thermal body (spec section 7).

Per macro step of dt: Simulator.advance(dt) (DOP853, events truncate the last step) -> aero state at the end of the
step -> heating with that state and the wall temperatures at the start of the step -> ThermalBody.advance (radiation
implicit in the solver) -> history row, and every `frames_every` steps a VTK frame (nodal T on the volume mesh,
q_conv / q_rad / T per patch on the surface). First-order operator splitting; the dt-halving test bounds its error.
Mass is constant in Step 2; the loop already carries the body's mass and the mesh so Step 3 can change both.
`body` must be a ThermalBody: CoupledRun uses `theta`, `surface`, `radiated_power()`, `surface_stats()`,
`nose_radius()` and `integrated_heat`, beyond what the `Body` protocol declares. `ConstantBody` is for `Simulator.run()` only."""
import os
import time
from dataclasses import dataclass, field

import numpy as np

from . import trajectory as tj

THERMAL_COLUMNS = ["convective_heat_W", "rad_cooling_W", "integrated_heat_J", "absorbed_heat_J",
                   "surface_T_max_K", "surface_T_min_K", "surface_T_mean_K", "T_stagnation_K", "T_back_K",
                   "T_centre_K", "q_stag_Wm2", "heating_blend_f"]
MELT_COLUMNS = ["film_mass_kg", "sprayed_mass_kg", "runoff_mass_kg", "removed_mass_kg", "melt_front_depth_max_mm",
                "equivalent_radius_mm", "n_active_elements", "spraying_area_m2", "theta_cr_deg", "n_released",
                "released_mass_kg", "r_median_um", "r_max_um", "closure_fraction_girin", "closure_fraction_couette_slip",
                "closure_fraction_couette_fm", "rt_active", "removed_enthalpy_J", "film_thickness_max_mm",
                "film_thickness_mean_mm", "nose_radius_mm", "transverse_radius_mm", "fitted_nose_radius_mm",
                "kn_body", "kn_local_stag", "re_shock", "flow_branch", "p_w_stag_Pa", "phi_sonic_deg",
                "drag_shape_factor", "frozen_mass_kg", "film_T_max_K", "film_T_mean_K", "film_frozen_fraction",
                "unapplied_load_J", "film_blob_fraction", "rt_mass_fraction", "rt_wavelength_over_nose", "rt_growth_ms", "spray_growth_ms", "rt_region_mm", "rt_bounded_fraction", "rt_bounded_growth_ms", "molten_depth_max_mm", "molten_depth_mean_mm",
                "delta_m_mean_um", "thick_branch_fraction", "n_dead_elements", "deep_liquid_kg", "deep_mass_kg",
                "deep_runoff_mass_kg", "deep_surfaced_mass_kg", "deep_blob_fraction", "cascade_passes", "cascade_mass_kg",
                "p_w_stag_step_Pa"]      # post-reconstruction addition (2026-10-07): every step's own stagnation wall pressure
PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'


@dataclass
class CoupledSettings:
    dt: float = 0.5                  # s, macro step
    frames_every: int = 0            # VTK frame every n macro steps (0: none)
    output_dir: str = None           # directory of the VTK series (required when frames_every > 0)
    frames: list = field(default_factory=list)


class CoupledRun:
    def __init__(self, sim, body, heating_model, settings=None):
        """`body` must be a ThermalBody (uses `theta`, `surface`, `radiated_power()`, `surface_stats()` and
        `integrated_heat`, beyond the `Body` protocol); `ConstantBody` is for `Simulator.run()` only."""
        self.sim, self.body, self.heating = sim, body, heating_model
        self.settings = settings or CoupledSettings()
        if self.settings.frames_every and not self.settings.output_dir:
            raise ValueError("frames_every > 0 needs an output_dir")

    def loads_at(self, t, y):
        a = self.sim.aero_state(t, y[:3], y[3:])
        return a, self.heating.evaluate(a, self.body.theta, self.body.surface_temperature(), self.body.nose_radius(),
                                        T_mean=self.body.mean_temperature())

    def row(self, t, y, a, loads):
        body = self.body
        r = self.sim.sample_row(t, y, a)
        Q = body.last.Q_conv if body.last is not None else loads.total(body.surface.areas)     # the step's applied power (the surface may have changed since)
        r.update({"convective_heat_W": Q, "rad_cooling_W": -body.radiated_power(),
                  "integrated_heat_J": body.integrated_heat, "absorbed_heat_J": body.absorbed_heat,
                  "q_stag_Wm2": loads.q_stag, "heating_blend_f": loads.blend})
        r.update(body.surface_stats())
        if self.melting:
            r.update(body.melt_stats())
        return r

    @property
    def melting(self):
        return hasattr(self.body, "melt_stats")

    def run(self):
        started = time.perf_counter()
        sim, body, s = self.sim, self.body, self.settings
        a, loads = self.loads_at(sim.t, sim.y)
        rows, states, step = [self.row(sim.t, sim.y, a, loads)], [sim.y.copy()], 0
        self._frame(step, sim.t, loads)
        while sim.end_reason is None:
            dt = sim.advance(s.dt)
            a, loads = self.loads_at(sim.t, sim.y)
            body.advance(sim.t, dt, loads, state=a)
            step += 1
            if len(loads.q_conv) != body.surface.n_patches:       # elements died this step: loads on the new surface for the record
                a, loads = self.loads_at(sim.t, sim.y)
            rows.append(self.row(sim.t, sim.y, a, loads))
            states.append(sim.y.copy())
            self._frame(step, sim.t, loads)
            if self.melting and body.demised():
                sim.end_reason = "demise"
        self._collection()
        columns = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
        history = tj.History(columns, np.array(states), sim.end_reason)
        history.results = sim.results(history, sim.nfev, time.perf_counter() - started)
        history.results.update(self.thermal_results(history))
        if self.melting:
            history.results.update(self.melt_results(history))
        return history

    def melt_results(self, history):
        c, body = history.columns, self.body
        first = lambda mask: (float(c["time_s"][mask][0]), float(c["altitude_km"][mask][0])) if mask.any() else (None, None)
        t_melt, h_melt = first(c["removed_mass_kg"] + c["film_mass_kg"] > 0.0)
        t_spray, h_spray = first(c["sprayed_mass_kg"] > 0.0)
        i_peak = int(np.argmax(c["released_mass_kg"]))
        return {"melt_onset_time_s": t_melt, "melt_onset_altitude_km": h_melt, "spraying_onset_time_s": t_spray,
                "spraying_onset_altitude_km": h_spray, "demise_time_s": float(c["time_s"][-1]) if history.end_reason == "demise" else None,
                "demise_altitude_km": float(c["altitude_km"][-1]) if history.end_reason == "demise" else None,
                "initial_mass_kg": body.mass0, "final_mass_kg": float(c["mass_kg"][-1]), "sprayed_mass_kg": body.sprayed_mass,
                "runoff_mass_kg": body.runoff_mass, "removed_mass_kg": body.removed_mass, "film_mass_kg": float(body.m_f.sum()),
                "deep_runoff_mass_kg": body.deep_runoff_mass, "deep_surfaced_mass_kg": body.deep_surfaced_mass,
                "deep_mass_kg": float(body.m_d.sum()),
                "cascade_mass_kg": body.cascade_mass, "cascade_passes_max": int(c["cascade_passes"].max()),
                "cascade_capped_steps": body.cascade_capped_steps,
                "n_released": body.n_released, "time_of_peak_release_s": float(c["time_s"][i_peak]),
                "r_median_um": float(np.nanmedian(c["r_median_um"])) if np.isfinite(c["r_median_um"]).any() else None,
                "n_dead_elements": int(body.mesh.n_elements - body.mesh.n_active), "n_source_rows": len(body.source_rows),
                "removed_enthalpy_J": body.removed_enthalpy, "frozen_mass_kg": body.frozen_mass,
                "melt_energy_balance_residual": body.energy_balance_residual()}

    def thermal_results(self, history):
        c, body = history.columns, self.body
        i_surf, i_mean, i_q = int(np.argmax(c["surface_T_max_K"])), int(np.argmax(c["temperature_K"])), int(np.argmax(c["convective_heat_W"]))
        return {
            "peak_surface_T_K": float(c["surface_T_max_K"][i_surf]), "time_of_peak_surface_T_s": float(c["time_s"][i_surf]),
            "altitude_of_peak_surface_T_km": float(c["altitude_km"][i_surf]),
            "peak_mean_T_K": float(c["temperature_K"][i_mean]), "time_of_peak_mean_T_s": float(c["time_s"][i_mean]),
            "peak_convective_heat_W": float(c["convective_heat_W"][i_q]), "time_of_peak_heating_s": float(c["time_s"][i_q]),
            "altitude_of_peak_heating_km": float(c["altitude_km"][i_q]),
            "integrated_heat_J": body.integrated_heat, "absorbed_heat_J": body.absorbed_heat, "radiated_heat_J": body.radiated_heat,
            "energy_balance_residual": body.energy_balance_residual(),
            "n_macro_steps": len(body.iterations), "mean_newton_iterations": float(np.mean(body.iterations)) if body.iterations else 0.0,
            "n_frames": len(self.settings.frames),
        }

    def _frame(self, step, t, loads):
        s = self.settings
        if not s.frames_every or step % s.frames_every:
            return
        k = len(s.frames)
        write_vtk_frame(s.output_dir, k, self.body, loads)
        s.frames.append((t, k))

    def _collection(self):
        s = self.settings
        if not s.frames:
            return
        for name in ("field", "surface"):
            entries = "".join('<DataSet timestep="{:.6g}" file="{}_{}.{}"/>\n'.format(t, name, k, "vtu" if name == "field" else "vtp")
                              for t, k in s.frames)
            with open(os.path.join(s.output_dir, name + ".pvd"), "w") as fh:
                fh.write(PVD_TEMPLATE.format(entries))


def write_vtk_frame(output_dir, k, body, loads):
    """field_<k>.vtu: nodal T on the active volume mesh (plus liquid fraction and element fractions when melting);
    surface_<k>.vtp: q_conv, q_rad, T per patch, the derived surface's outward unit normal n_derived, and when melting the film thickness and temperature, We_s, closure,
    local Knudsen number, wall pressure, shear, droplet radius, release rate and Girin's conjugate depth delta_m -- the
    flow and spray fields of the step's own evaluation, carried across that step's element deaths by face id
    (`MeltingBody.on_current_surface`), with the defaults on faces the deaths exposed (nan for delta_m, which is also
    nan wherever the closure is not Girin's) -- and the thickness of the deep liquid beneath the film (PyVista/VTK XML).
    Since 2026-10-08 the gas-side fields (closure, kn_local, p_w, tau, delta_m) of the faces the deaths exposed come from
    the same flow model evaluated for the record on the current surface (`flow_eval` 2); the spray fields keep their
    defaults there (nothing was sprayed from a face that did not yet exist)."""
    import pyvista as pv
    from .thermal import SIGMA_SB
    os.makedirs(output_dir, exist_ok=True)
    mesh, surface, T = body.mesh, body.surface, body.field()
    melting = hasattr(body, "melt_stats")
    tets = mesh.tets[mesh.active]
    cells = np.hstack([np.full((len(tets), 1), 4), tets]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(len(tets), pv.CellType.TETRA), mesh.points)
    grid.point_data["T"] = T
    if melting:
        grid.point_data["liquid_fraction"] = body.material.liquid_fraction(T)
        grid.cell_data["phi"] = body.phi[mesh.active]
    grid.save(os.path.join(output_dir, "field_{}.vtu".format(k)))
    faces = np.hstack([np.full((surface.n_patches, 1), 3), surface.faces]).ravel()
    poly = pv.PolyData(mesh.points, faces)
    Tf = surface.facet_mean(T)
    poly.point_data["T"] = T
    poly.cell_data["q_conv"] = loads.q_conv
    poly.cell_data["q_rad"] = body.emissivity * SIGMA_SB * (Tf ** 4 - body.T_ambient ** 4)
    poly.cell_data["T_patch"] = Tf
    # post-reconstruction addition (2026-10-07, write-only; Spheral M1 decision 9): each patch's outward unit normal on
    # the derived surface of sub-plan 01 -- the Taubin-smoothed normal of `mesh.surface(derived=True)` (spec 2026-09-27
    # sec. 8 step 4) -- which the layer depths march along. The derived surface's patches are these patches one to one
    # (the same boundary faces, in the same order: checked by face id). The triangles are written with the face table's
    # node order, whose winding is inward on some patches; this field is outward by construction (oriented by the
    # patch's opposite vertex, like `surface.normals`), so it, not the winding, is the orientation to trust. Nothing in
    # the prototype displaces `derived_points` yet (sub-plan 16, plans only), so the derived surface's facets are the
    # staircase's and the field is the smoothing alone.
    derived = mesh.surface(derived=True)
    if not np.array_equal(derived.face_ids, surface.face_ids):
        raise RuntimeError("the derived surface's patches are not the body's surface patches")
    poly.cell_data["n_derived"] = derived.smoothed_normals() if surface.n_patches else np.zeros((0, 3))
    if melting:
        poly.cell_data["film_thickness"] = body.m_f / (body.liquid.rho * surface.areas)
        poly.cell_data["film_T"] = body.film_temperature()
        res, flow, carry = body.last_spray, body.last_flow, body.on_current_surface
        poly.cell_data["we_s"] = carry(res.we_s if res is not None else None, 0.0)
        # post-reconstruction addition (2026-10-07, write-only): closure, p_w and tau from the step's own flow evaluation,
        # which every melt step with an aero state makes; the spray step's `last_flow` is that same object whenever there
        # is film, and None otherwise, so these fields were the defaults before melt onset and after spraying ended. The
        # defaults (closure 1, p_w and tau 0) now remain only where the flow was not evaluated -- frame 0, before the
        # first melt step; --removal instant, which never evaluates it -- and on faces the step's deaths exposed.
        step_flow = getattr(body, "last_flow_step", None) or flow
        # post-reconstruction addition (2026-10-08, write-only): the faces this step's deaths exposed take the gas-side
        # fields of the same surface flow evaluated at the step's aero state on the current surface (`flow_for_the_record`);
        # `flow_eval` says which: 1 the step's own evaluation, 2 the record's, 0 none (frame 0, --removal instant)
        covered, record = body.flow_for_the_record() if hasattr(body, "flow_for_the_record") else (np.zeros(surface.n_patches, bool), None)

        def gas_side(values, fill, record_values):
            v = carry(values, fill)
            return np.where(covered, v, np.asarray(record_values, dtype=float)) if record is not None else v
        poly.cell_data["flow_eval"] = np.where(covered, 1.0, 2.0 if record is not None else 0.0) if step_flow is not None \
            else np.zeros(surface.n_patches)
        poly.cell_data["closure"] = gas_side(step_flow.closure.astype(float) if step_flow is not None else None, 1.0,
                                             record.closure.astype(float) if record is not None else None)
        poly.cell_data["kn_local"] = gas_side(flow.kn_local if flow is not None else None, np.nan,
                                              record.kn_local if record is not None else None) if flow is not None else \
            carry(None, np.nan)
        poly.cell_data["p_w"] = gas_side(step_flow.p_w if step_flow is not None else None, 0.0, record.p_w if record is not None else None)
        poly.cell_data["tau"] = gas_side(step_flow.tau if step_flow is not None else None, 0.0, record.tau if record is not None else None)
        poly.cell_data["r_droplet"] = carry(np.where(res.dm > 0.0, res.r, np.nan) if res is not None else None, np.nan)
        poly.cell_data["release_rate"] = carry(res.dm if res is not None else None, 0.0) / surface.areas
        # Girin's conjugate depth per patch, nan where the patch's closure is not his (no conjugate depth exists there),
        # and the deep liquid lying below it (amendment of 2026-10-02): the skin and the runoff layer beneath it
        delta_record = None
        if record is not None:
            from . import spray as spray_mod
            delta_record = spray_mod.melt_layer(record, body.liquid)[0]
        poly.cell_data["delta_m"] = gas_side(body.last_delta_m, np.nan, delta_record) if body.last_delta_m is not None else carry(None, np.nan)
        poly.cell_data["deep_thickness"] = body.m_d / (body.liquid.rho * surface.areas)
    poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))


def write_particles(run_dir, body, history, window=10.0):
    """particles.npz (the source table, SOURCE_COLUMNS as arrays), particles_summary.csv (per macro step) and
    size_distribution.csv (dn, dM per log bin per `window`-second window and in total). Returns the paths."""
    import csv
    from . import spray
    os.makedirs(run_dir, exist_ok=True)
    rows = np.array(body.source_rows, dtype=float) if body.source_rows else np.zeros((0, len(spray.SOURCE_COLUMNS)))
    npz = os.path.join(run_dir, "particles.npz")
    np.savez_compressed(npz, **{name: rows[:, i] for i, name in enumerate(spray.SOURCE_COLUMNS)})
    c = history.columns
    summary = os.path.join(run_dir, "particles_summary.csv")
    with open(summary, "w", newline="") as fh:
        w = csv.writer(fh)
        keys = ["time_s", "altitude_km", "velocity_kms", "released_mass_kg", "n_released", "r_median_um", "r_max_um", "theta_cr_deg",
                "spraying_area_m2", "film_mass_kg", "film_thickness_mean_mm"]
        w.writerow(keys)
        for i in range(len(history)):
            w.writerow(["{:.9g}".format(c[k][i]) for k in keys])
    dist = os.path.join(run_dir, "size_distribution.csv")
    edges = spray.BIN_EDGES
    t = rows[:, 0] if len(rows) else np.zeros(0)
    t_end = float(c["time_s"][-1])
    starts = np.arange(0.0, t_end + 1e-9, window)
    with open(dist, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["window_start_s", "window_end_s", "r_lo_m", "r_hi_m", "dn", "dm_kg"])
        for t0 in list(starts) + [None]:
            if t0 is None:
                sel = np.ones(len(rows), dtype=bool)
                lo, hi = 0.0, t_end
            else:
                lo, hi = t0, t0 + window
                sel = (t >= lo) & (t < hi)
            n_hist, m_hist = spray.histogram(rows[sel, 12], rows[sel, 13], rows[sel, 14]) if sel.any() else (np.zeros(spray.N_BINS), np.zeros(spray.N_BINS))
            for j in range(spray.N_BINS):
                w.writerow(["{:.9g}".format(lo), "{:.9g}".format(hi), "{:.9g}".format(edges[j]), "{:.9g}".format(edges[j + 1]),
                            "{:.9g}".format(n_hist[j]), "{:.9g}".format(m_hist[j])])
    return {"particles": npz, "particles_summary": summary, "size_distribution": dist}
