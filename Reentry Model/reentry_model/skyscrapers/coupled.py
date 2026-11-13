"""Lockstep coupling of the trajectory stepper with the thermal body (spec section 7).

Per macro step of dt: Simulator.advance(dt) (DOP853, events truncate the last step) -> aero state at the end of the
step -> heating with that state and the wall temperatures at the start of the step -> ThermalBody.advance (radiation
implicit in the solver) -> history row, and every `frames_every` steps a VTK frame (nodal T on the volume mesh,
q_conv / q_rad / T per patch on the surface). First-order operator splitting; the dt-halving test bounds its error.
Mass is constant in Step 2; the loop already carries the body's mass and the mesh so Step 3 can change both.
`body` must be a ThermalBody: CoupledRun uses `theta`, `surface`, `radiated_power()`, `surface_stats()`,
`nose_radius()` and `integrated_heat`, beyond what the `Body` protocol declares. `ConstantBody` is for `Simulator.run()` only."""
import os
import signal
import time
from dataclasses import dataclass, field

import numpy as np

from . import restart
from ..shared import trajectory as tj
from ..shared.droplets import SKIN_SUMMARY_COLUMNS, SUMMARY_COLUMNS, write_particles  # noqa: F401  re-exported
from ..shared.spray import RELEASE_ANGLE_COLUMNS
from ..shared.surface_flow import KN_BODY_SHOCK

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
                "nonrigid_depth_mean_mm", "slurry_thick_fraction", "rigid_thin_fraction", "slurry_held_mass_kg",
                "thin_invalid_patches", "thin_invalid_mass_kg",
                # the regime test's depth over the wet windward patches (fact 94); sub-plan 18's skin columns, which held
                # 0 or nan in every run of the element model, went with its one-dimensional skins (sub-plan 22 Task 3)
                "nonrigid_depth_wet_median_um"]
# both surface models' histories also carry Register et al.'s release angle over each step's released mass (its mean by
# mass, least and largest), spray.RELEASE_ANGLE_COLUMNS (Asha 2026-10-29, sub-plan 22 Task 4), after their own columns
# the layered skins on the corner columns (sub-plan 20): skinbody.SkinBody.skin_stats(), every macro step
SKINBODY_COLUMNS = ["skin_cells", "skin_thickness_max_mm", "skin_above_min", "flowing_depth_mean_mm", "flowing_depth_max_mm",
                    "liquid_depth_mean_mm", "liquid_depth_max_mm", "surface_quality_min", "surface_fold_deg", "zone_ratio_min",
                    "rebuilds", "rebuild_residual_J", "rebuild_residual_last_J", "drawn_mass_kg", "absorbed_mass_kg",
                    "handed_mass_kg", "shortfall_kg", "jump_faults", "jump_faults_bottom_cell", "absorb_capped", "interface_mismatch_J",
                    "interface_repeats", "repeat_columns_over", "mismatch_worst_ratio", "mismatch_total_ratio",
                    "removed_mass_kg", "removed_enthalpy_J", "base_recession_max_mm",
                    "wall_skins_s", "wall_3d_s", "wall_move_s", "wall_rebuild_s", "step_parts",
                    # the flow (sub-plan 21 Task 6c, spec 10.2)
                    "momentum_residual", "speed_max_ms", "roughness_max", "flow_courant_max", "flow_growth_max",
                    "flow_passes_max", "flow_cycles_max", "wall_flow_s",
                    # the droplet momentum sink (sub-plan 21 Task 7, spec 10.2): fact 29's Pi and its flag
                    "released_mass_kg", "droplet_momentum_ratio", "droplet_momentum_flagged",
                    # what the fault tolerance holds in the skins at the last base move (sub-plan 21 Task 9, spec 6.6)
                    "held_mass_kg", "held_depth_max_mm",
                    # Girin's spraying from the top cells (sub-plan 22 Task 4, spec 10.2): the step's released mass -- the
                    # released_mass_kg above is since time zero -- the droplets since time zero, the radius median and
                    # largest over the step's records, the smallest angle and the top area of the columns whose plan
                    # sprays, the element model's per-step columns of the same names, which particles_summary.csv reads
                    "released_step_mass_kg", "n_released", "r_median_um", "r_max_um", "theta_cr_deg", "spraying_area_m2",
                    # Register et al.'s release angle over the step's released mass (Asha 2026-10-29)
                    *RELEASE_ANGLE_COLUMNS,
                    # spec 10.2: b and the resolved sheared depth over Girin's delta_m on the wet windward columns, and
                    # each rebuild's mass residual beside its energy residual
                    "b_over_delta_m_median", "sheared_depth_over_delta_m_median",
                    "rebuild_mass_residual_kg", "rebuild_mass_residual_last_kg",
                    # the windward loads' nose radius (sub-plan 22 Task 5g, Asha 2026-11-10): the radius of the surface at the end of the step, which the next step's loads read, the
                    # sphere fitted to the tops' cap, the reported bound, the corners fitted and 1 where the cap was concave
                    "nose_radius_mm", "fitted_nose_radius_mm", "nose_bound_mm", "nose_cap_points", "nose_cap_concave",
                    # the leeward test (sub-plan 22 Task 5j, Asha 2026-11-10): the mass centroid's position along the flight
                    # direction relative to the 3D mesh's mean point [mm, positive forward], and the corners that face the flow
                    # from behind the plane through it
                    "centroid_offset_mm", "leeward_corners"]
# with --spread-limit only, skinbody.SPREAD_COLUMNS follow these (sub-plan 22 Task 5c): a run's columns are its first row's, so
# a resume that gives or drops a limit adds them to its saved row or takes them from it (cli.cmd_resume)
END_SNAP = 1e-6                   # of the macro step: a step that would end this close short of t_max ends at it (below)
PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'


class FlightStopped(BaseException):
    """A flight stopped by the signal `signum`: the command line's SIGTERM handler raises it during `CoupledRun.run`, so that
    `kill` stops a flight as Control-C does, saving its last completed macro step (Task 4 of the geometry stop, Asha
    2026-11-12). A BaseException, as KeyboardInterrupt is, so that no `except Exception` on the way (a save's, a fallback's,
    a wrapper that turns errors into its own) takes the stop for an error and flies on."""

    def __init__(self, signum):
        self.signum = int(signum)
        super().__init__("stopped by signal {} ({})".format(self.signum, signal.Signals(self.signum).name))


STOPS = (KeyboardInterrupt, FlightStopped)       # a stop by hand: Control-C, or a signal the command line turns into one


@dataclass
class CoupledSettings:
    dt: float = 0.5                  # s, macro step
    frames_every: int = 0            # VTK frame every n macro steps (0: none); by flight time after a continuum switch
    output_dir: str = None           # directory of the VTK series (required when frames_every > 0)
    frames: list = field(default_factory=list)
    dt_continuum: float = None       # s, the macro step from the first step whose body Knudsen number is below kn_switch,
                                     # latched (None: no switch; amendment of 2026-10-07, the continuum step)
    kn_switch: float = KN_BODY_SHOCK  # the body Knudsen number of the switch: the surface flow's own continuum boundary
    switch_on_melt: bool = False      # the switch also latches at the first macro step at whose end the body holds liquid
                                      # (`body.liquid_present()`), if that comes first: the layered skins' (sub-plan 22
                                      # Task 3b, Asha 2026-10-29); the element model's switch is the Knudsen number's alone
    snap_to_t_max: bool = False       # a macro step that would end within END_SNAP of a step short of t_max ends at it:
                                      # the steps' sum falls short of t_max by round-off, and the loop then took a last
                                      # step of about 1e-17 s, a whole macro step of the skins over which the droplet-
                                      # momentum ratio read 5.5 and 8.6 (fact 179); the layered skins' (sub-plan 22 Task
                                      # 4), the element model keeping its loop, and that step, to the last bit
    save_at: float = None             # s of flight: the run saves its whole state once, at the end of the first macro
                                      # step ending at or after it, to save_path (`restart.save`; sub-plan 22 Task 5a)
    save_path: str = None
    save_meta: dict = None            # kept with the state: what made the run (the command line's arguments and name)
    end_path: str = None              # where the run saves its whole state when the flight ends, whatever the end (the state's
                                      # "end" is `sim.end_reason`); None: no end state

    def __post_init__(self):
        if self.dt_continuum is not None and not 0.0 < self.dt_continuum <= self.dt:
            raise ValueError("dt_continuum must be > 0 and no longer than dt, got {!r}".format(self.dt_continuum))
        if self.save_at is not None and not (self.save_at > 0.0 and self.save_path):
            raise ValueError("save_at must be > 0 and needs a save_path, got {!r} and {!r}".format(self.save_at, self.save_path))


class CoupledRun:
    geometry_stopped = None            # (test, time s, value) where a skin flight ended at a geometry breach; old states have none

    def __init__(self, sim, body, heating_model, settings=None):
        """`body` must be a ThermalBody (uses `theta`, `surface`, `radiated_power()`, `surface_stats()` and
        `integrated_heat`, beyond the `Body` protocol); `ConstantBody` is for `Simulator.run()` only."""
        self.sim, self.body, self.heating = sim, body, heating_model
        self.settings = settings or CoupledSettings()
        if self.settings.frames_every and not self.settings.output_dir:
            raise ValueError("frames_every > 0 needs an output_dir")
        self.switched_at = self.switch_trigger = None              # the continuum switch: (t, Kn, altitude km), and why
        self.state_saved_at = self.state_save_error = None        # the save's flight time [s], or why it failed (Task 5a)
        self.end_saved_at = self.end_save_error = None            # the same for the state saved at the flight's end

    def loads_at(self, t, y):
        a = self.sim.aero_state(t, y[:3], y[3:])
        # the layered skins' effective angles (`SkinBody.load_angles`, sub-plan 22 Task 5j): a corner behind the centroid's
        # plane is leeward to the heating whichever way its normal faces; every other body's angles are `theta`
        angles = self.body.load_angles() if hasattr(self.body, "load_angles") else self.body.theta
        return a, self.heating.evaluate(a, angles, self.body.surface_temperature(), self.body.nose_radius(),
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
        if self.skins:
            r.update(body.skin_stats())
        return r

    @property
    def melting(self):
        return hasattr(self.body, "melt_stats")

    @property
    def skins(self):
        return hasattr(self.body, "skin_stats")

    def run(self, resume=None):
        """Fly to the end, from the start or, with `resume` -- a state this run saved (`settings.save_at`), read back by
        `restart.load` (sub-plan 22 Task 5a) -- from the end of the saved macro step: the history then starts with that
        step's row, and the continuum switch, the frames' clock and the step count are the saved run's.

        With `settings.end_path` (Task 4 of the geometry stop, Asha 2026-11-12) the run keeps a copy of its state at the
        end of the last completed step, and a step that raises, or a stop by hand (Control-C, or the command line's
        SIGTERM), ends the flight there instead of escaping: see `_failed`. Without it an error escapes as before."""
        started = time.perf_counter()
        sim, body, s = self.sim, self.body, self.settings
        if resume is None:
            self.switched_at = self.switch_trigger = None          # (t, Kn, altitude km) of the continuum switch, and why
            self._t0, self._frame_dt, self._next_frame = sim.t, s.frames_every * s.dt, sim.t
            a, loads = self.loads_at(sim.t, sim.y)
            rows, states, step = [self.row(sim.t, sim.y, a, loads)], [sim.y.copy()], 0
            self._frame(step, sim.t, loads)
        else:
            if resume.get("format") != restart.FORMAT or resume.get("run") is not self:
                raise ValueError("a run resumes only from a state of its own of format {}".format(restart.FORMAT))
            rows, states, step = [dict(resume["row"])], [sim.y.copy()], resume["step"]
            if self.skins and "jump_faults_bottom_cell" not in vars(body):
                # a state saved before sub-plan 22 Task 5k (fix round 1, controller's ruling 2026-11-11): the count of the old rule,
                # the bottom cell alone, is the saved `jump_faults`, which was made by that rule up to the save, so the old count
                # stays continuous; `jump_faults` itself continues under the new rule from the resume. Here, where `cmd_resume`
                # and a direct resume both pass, and the saved row gets the column too
                body.jump_faults_bottom_cell = body.jump_faults
                rows[0].setdefault("jump_faults_bottom_cell", rows[0].get("jump_faults", body.jump_faults))
        # the copy of the last completed step (Task 4): a local, never an attribute of the run, so that no state pickled
        # after it (the next copy, the state at save_at, the end state) carries it and grows with every step
        last = self._copy(step, rows) if s.end_path is not None else None
        try:
            while sim.end_reason is None:
                dt = sim.advance(s.dt, snap=END_SNAP * s.dt if s.snap_to_t_max else 0.0)
                a, loads = self.loads_at(sim.t, sim.y)
                body.advance(sim.t, dt, loads, state=a)
                step += 1
                self._switch(sim.t, a)
                if len(loads.q_conv) != body.surface.n_patches:       # elements died this step: loads on the new surface for the record
                    a, loads = self.loads_at(sim.t, sim.y)
                rows.append(self.row(sim.t, sim.y, a, loads))
                states.append(sim.y.copy())
                self._frame(step, sim.t, loads)
                if (self.melting or self.skins) and body.demised():
                    sim.end_reason = "demise"
                if self.skins and sim.end_reason is None and body.geometry_stop:     # any other end on this step wins (open item 4)
                    breach = body.geometry_breach()
                    if breach is not None:
                        sim.end_reason = "geometry"
                        self.geometry_stopped = (breach[0], float(sim.t), float(breach[1]))
                if (s.save_at is not None and sim.end_reason is None and self.state_saved_at is None
                        and self.state_save_error is None and sim.t >= s.save_at - END_SNAP * s.dt):
                    self._save(step, rows[-1])
                if s.end_path is not None and sim.end_reason is None:     # the step is complete and the flight goes on
                    last = self._copy(step, rows, last)
        except (Exception,) + STOPS as exc:
            if last is None:
                raise
            return self._failed(exc, last, rows, states, started)
        self._collection()
        if s.end_path is not None:
            self.end_saved_at = self.end_save_error = None
            self._save(step, rows[-1], path=s.end_path, end=sim.end_reason)
        return self._history(rows, states, started)

    def _history(self, rows, states, started):
        """The flight's history from its rows and trajectory states, with its results: what `run` returns (Task 4: the end
        of `run`, which a flight that failed also calls on the run of its last completed step)."""
        sim, s = self.sim, self.settings
        columns = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
        history = tj.History(columns, np.array(states), sim.end_reason)
        history.results = sim.results(history, sim.nfev, time.perf_counter() - started)
        history.results.update(self.thermal_results(history))
        if self.melting:
            history.results.update(self.melt_results(history))
        if self.skins:
            history.results.update(self.skin_results(history))
        if s.dt_continuum is not None:
            t_sw, kn_sw, h_sw = self.switched_at or (None, None, None)
            history.results.update({"dt_continuum_s": s.dt_continuum, "dt_switch_time_s": t_sw, "dt_switch_kn": kn_sw,
                                    "dt_switch_altitude_km": h_sw})
            if s.switch_on_melt:                                   # what latched it: knudsen, melting or both at once
                history.results["dt_switch_trigger"] = self.switch_trigger
        if s.save_at is not None:
            history.results.update({"state_saved_at_s": self.state_saved_at, "state_save_error": self.state_save_error})
        if s.end_path is not None:
            history.results.update({"end_state_saved_at_s": self.end_saved_at, "end_state_error": self.end_save_error})
        return history

    def _save(self, step, row, path=None, end=None):
        """Save the whole run at the end of macro step `step`, whose history row is `row` (sub-plan 22 Task 5a): the run,
        the step, its row and flight time, numpy's global random state, `settings.save_meta` and `end`, why the flight
        ended (None for a state saved on its way), by `restart.save` to `path` (default `settings.save_path`), which
        alters nothing the run holds. A save that fails is recorded (as the end state's, when `path` is given) and the
        flight goes on: its state is worth less than the flight."""
        s = self.settings
        ending = path is not None
        try:
            restart.save(path if ending else s.save_path, self._state(step, row, end))
            if ending:
                self.end_saved_at = float(self.sim.t)
            else:
                self.state_saved_at = float(self.sim.t)
        except Exception as exc:
            message = "{}: {}".format(type(exc).__name__, exc)
            if ending:
                self.end_save_error = message
            else:
                self.state_save_error = message

    def _state(self, step, row, end=None):
        """The state `restart` saves of the run at the end of macro step `step`, whose history row is `row`."""
        return {"format": restart.FORMAT, "run": self, "step": step, "row": row, "time_s": float(self.sim.t),
                "rng": np.random.get_state(), "meta": self.settings.save_meta, "end": end}

    def _copy(self, step, rows, last=None):
        """(the run's state at the end of macro step `step`, pickled in memory as `restart.save` pickles it, with "end"
        None; the number of history rows to that step), which a failure or a stop by hand writes as the end state (Task
        4). A copy that cannot be made leaves `last`, the previous step's, as a save that fails never fails the flight; a
        stop by hand during it leaves it too, since the new copy replaces the old only once it is complete."""
        try:
            return restart.dumps(self._state(step, rows[-1])), len(rows)
        except Exception:
            return last

    def _failed(self, exc, last, rows, states, started):
        """The end of a flight whose macro step raised `exc`, or that was stopped by hand (Task 4; Asha, 2026-11-12: "the
        code should always produce a saved state when it terminates"). `last`, the copy of the last completed step, is
        loaded back as a run of its own -- consistent, unlike this one, which the step left part way -- and saved to
        `settings.end_path` with "end" `failed`, or `stopped` for Control-C or FlightStopped. The run is saved as it was,
        its flight not ended (`sim.end_reason` None), so a resume flies on from it and repeats the failing step; only
        then is its end reason set. A save that fails, a second Control-C included, is recorded and changes nothing
        else. Returns that run's history of the completed steps, with `history.run` that run, `history.exception` `exc`,
        and `results["failure"]` = {"kind", "error": "<Type>: <message>", "time_s": the last completed step's time}."""
        kind = "stopped" if isinstance(exc, STOPS) else "failed"
        data, n = last
        state = restart.loads(data)
        state["end"] = kind
        run = state["run"]
        run.end_saved_at = run.end_save_error = None
        try:
            restart.save(run.settings.end_path, state)
            run.end_saved_at = state["time_s"]
        except (Exception,) + STOPS as error:
            run.end_save_error = "{}: {}".format(type(error).__name__, error)
        run.sim.end_reason = kind
        run._collection()
        history = run._history(rows[:n], states[:n], started)
        history.run, history.exception = run, exc
        history.results["failure"] = {"kind": kind, "error": "{}: {}".format(type(exc).__name__, exc), "time_s": state["time_s"]}
        return history

    def _switch(self, t, a):
        """The continuum step (amendment of 2026-10-07): from the first step whose aero state -- the one the body was
        just given, before any re-evaluation on a surface the step's deaths changed -- has a body Knudsen number below
        kn_switch, the macro step is dt_continuum, latched. Girin's closure begins there and the droplet population needs
        the fine step from there on (plan facts 92-96); before it the default step serves, as the 2026-10-05 harness
        that this replaces bit for bit had it. With `switch_on_melt` (the layered skins, sub-plan 22 Task 3b) it also
        latches at the first step at whose end the body holds liquid, if that comes first, since the skins spray from
        their first liquid: the trigger is recorded, `knudsen`, `melting` or `both`."""
        s = self.settings
        if s.dt_continuum is None or self.switched_at is not None:
            return
        knudsen = bool(np.isfinite(a.kn) and a.kn < s.kn_switch)
        melting = bool(s.switch_on_melt and self.body.liquid_present())
        if not (knudsen or melting):
            return
        s.dt = s.dt_continuum
        self.switched_at = (float(t), float(a.kn), float(a.h) / 1e3)
        self.switch_trigger = "both" if knudsen and melting else "knudsen" if knudsen else "melting"

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
                "thin_invalid_patch_steps": body.thin_invalid_patch_steps, "thin_invalid_mass_max_kg": body.thin_invalid_mass_max,
                "n_released": body.n_released, "time_of_peak_release_s": float(c["time_s"][i_peak]),
                "r_median_um": float(np.nanmedian(c["r_median_um"])) if np.isfinite(c["r_median_um"]).any() else None,
                "n_dead_elements": int(body.mesh.n_elements - body.mesh.n_active), "n_source_rows": len(body.source_rows),
                "removed_enthalpy_J": body.removed_enthalpy, "frozen_mass_kg": body.frozen_mass,
                "melt_energy_balance_residual": body.energy_balance_residual()}

    def skin_results(self, history):
        """The layered skins' totals (sub-plan 20)."""
        body, last = self.body, history.columns
        out = {"skin_columns": body.columns.n, "initial_mass_kg": body.mass0, "final_mass_kg": float(last["mass_kg"][-1]),
               "demise_time_s": float(last["time_s"][-1]) if history.end_reason == "demise" else None,
               "surface_fold_max_deg": float(np.max(last["surface_fold_deg"])),
               "surface_quality_min": float(np.min(last["surface_quality_min"])),
               "zone_ratio_min": float(np.min(last["zone_ratio_min"])),
               "momentum_residual": float(last["momentum_residual"][-1]),
               "flow_courant_max": float(np.max(last["flow_courant_max"])),
               "released_mass_kg": float(last["released_mass_kg"][-1]),
               "droplet_momentum_ratio_max": float(np.max(last["droplet_momentum_ratio"])),
               "droplet_momentum_flagged_steps": int((last["droplet_momentum_flagged"] > 0).sum()),
               "held_mass_max_kg": float(np.max(last["held_mass_kg"])),          # what the fault tolerance held (Task 9)
               "held_depth_max_mm": float(np.max(last["held_depth_max_mm"])),
               # Girin's spraying (sub-plan 22 Task 4): the droplet records' mass, the droplets and the rows written --
               # the element model's results of the same names
               "sprayed_mass_kg": body.sprayed_mass, "n_released": body.droplet_number,
               "n_source_rows": body.n_source_rows}
        out.update({k: float(last[k][-1]) for k in ("removed_mass_kg", "removed_enthalpy_J", "drawn_mass_kg", "absorbed_mass_kg",
                                                  "handed_mass_kg", "rebuild_residual_J", "wall_skins_s", "wall_3d_s",
                                                  "wall_move_s", "wall_rebuild_s", "wall_flow_s", "rebuild_mass_residual_kg")})
        out.update({k: int(last[k][-1]) for k in ("rebuilds", "jump_faults", "jump_faults_bottom_cell", "absorb_capped", "interface_repeats")})
        stop = self.geometry_stopped                                            # what ended the flight at a geometry breach, or None
        out["geometry_stop"] = None if stop is None else {"test": stop[0], "time_s": stop[1], "value": stop[2]}
        if body.spread_limit is not None:                                       # the deposit's spreading (sub-plan 22 Task 5c)
            out.update({"spread_limit_deg": body.spread_limit, "spread_mass_kg": body.spread_mass})
        return out

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
        """Every `frames_every` steps; after a continuum switch by flight time instead, every `frames_every` first
        steps, so that a fine step does not write a frame per step."""
        s = self.settings
        if not s.frames_every:
            return
        if self.switched_at is None or t <= self.switched_at[0] + 1e-9:
            if step % s.frames_every:
                return
        elif t + 1e-9 < self._next_frame:
            return
        k = len(s.frames)
        write_vtk_frame(s.output_dir, k, self.body, loads)
        s.frames.append((t, k))
        while self._next_frame <= t + 1e-9:
            self._next_frame += self._frame_dt

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
    surface_<k>.vtp: q_conv, q_rad, T per patch, and when melting the film thickness and temperature, We_s, closure,
    local Knudsen number, wall pressure, shear, droplet radius, release rate and Girin's conjugate depth delta_m -- the
    flow and spray fields of the step's own evaluation, carried across that step's element deaths by face id
    (`MeltingBody.on_current_surface`), with the defaults on faces the deaths exposed (nan for delta_m, which is also
    nan wherever the closure is not Girin's) -- the thickness of the deep liquid beneath the film, and the non-rigid
    depth the regime test read (PyVista/VTK XML)."""
    import pyvista as pv
    from ..thermal import SIGMA_SB
    os.makedirs(output_dir, exist_ok=True)
    if hasattr(body, "skin_stats"):
        return write_skin_frame(output_dir, k, body, loads)
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
    if melting:
        poly.cell_data["film_thickness"] = body.m_f / (body.liquid.rho * surface.areas)
        poly.cell_data["film_T"] = body.film_temperature()
        res, flow, carry = body.last_spray, body.last_flow, body.on_current_surface
        poly.cell_data["we_s"] = carry(res.we_s if res is not None else None, 0.0)
        poly.cell_data["closure"] = carry(flow.closure.astype(float) if flow is not None else None, 1.0)
        poly.cell_data["kn_local"] = carry(flow.kn_local if flow is not None else None, np.nan)
        poly.cell_data["p_w"] = carry(flow.p_w if flow is not None else None, 0.0)
        poly.cell_data["tau"] = carry(flow.tau if flow is not None else None, 0.0)
        poly.cell_data["r_droplet"] = carry(np.where(res.dm > 0.0, res.r, np.nan) if res is not None else None, np.nan)
        poly.cell_data["release_rate"] = carry(res.dm if res is not None else None, 0.0) / surface.areas
        # Girin's conjugate depth per patch, nan where the patch's closure is not his (no conjugate depth exists there),
        # and the deep liquid lying below it (amendment of 2026-10-02): the skin and the runoff layer beneath it
        poly.cell_data["delta_m"] = carry(body.last_delta_m, np.nan)
        poly.cell_data["deep_thickness"] = body.m_d / (body.liquid.rho * surface.areas)
        # the depth the regime test read under each patch: down to the first rigid point, slurry included (amendment of
        # 2026-10-06); the spray step's own, carried like delta_m, nan where it was not evaluated
        poly.cell_data["nonrigid_depth"] = carry(body.last_nonrigid, np.nan)
    poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))


def write_skin_frame(output_dir, k, body, loads):
    """A frame of the layered-skin body (sub-plan 20, spec 10.2): field_<k>.vtu, the 3D model's nodal T on its current
    (moving) mesh; surface_<k>.vtp, the top surface with per-corner point data -- top temperature, convective flux,
    liquid and flowing depths, skin thickness and cells, the base's position (the corner on the 3D mesh) and the top
    cell's speed (sub-plan 21 Task 6c)."""
    import pyvista as pv
    mesh, sk, mat = body.mesh, body.skins, body.material
    cells = np.hstack([np.full((mesh.n_elements, 1), 4), mesh.tets]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(mesh.n_elements, pv.CellType.TETRA), mesh.points)
    grid.point_data["T"] = body.field()
    grid.save(os.path.join(output_dir, "field_{}.vtu".format(k)))
    surface = body.surface
    poly = pv.PolyData(surface.points, np.hstack([np.full((len(surface.faces), 1), 3), surface.faces]).ravel())
    poly.point_data["T_top"] = sk.T[:, 0]
    poly.point_data["q_conv"] = loads.q_conv
    poly.point_data["skin_thickness"] = sk.thickness()
    poly.point_data["skin_cells"] = sk.count.astype(float)
    poly.point_data["liquid_depth"] = sk.depth_above(mat.T_liquidus) if mat.melts else np.zeros(sk.n)
    poly.point_data["flowing_depth"] = sk.depth_above(mat.T_rigid) if mat.melts else np.zeros(sk.n)
    poly.point_data["base_position"] = body.columns.points
    poly.point_data["top_speed"] = np.where(sk.count > 0, np.linalg.norm(sk.U[:, 0], axis=1), 0.0)
    poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))
