"""coupled.py: the lockstep loop with an inert body (Step 1 regression) and with the thermal body."""
import math
import os
from datetime import datetime

import numpy as np
import pytest

from reentry_model import aero, atmosphere, body, coupled, heating, material, mesh, thermal
from reentry_model import trajectory as tj

EPOCH = datetime(2024, 8, 1, 12, 53, 7)
R100 = tj.InitialState(7500.0, 77500.133, math.radians(-0.959331), math.radians(347.168296),
                       math.radians(29.546067), math.radians(-82.134333), EPOCH)
MASS_100MM = body.sphere_mass(0.1, 2813.0)


def simulator(the_body, **settings):
    return tj.Simulator(R100, the_body, atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(), aero.SesamTable(),
                        tj.Settings(diameter=0.1, **settings))


def thermal_body(the_mesh, T0=300.0, **solver_options):
    mat = material.Material.from_drama_json()
    return body.ThermalBody(the_mesh, mat, thermal.thermal_solver("skfem", **solver_options), MASS_100MM, T0=T0)


def test_coupled_full_flight_trajectory_matches_run(coarse_sphere_mesh):
    """The coupled loop (constant mass) over the whole 100 mm US76 flight, dt 2 s, vs the dense-output run() of Step 1:
    < 0.01 m/s, < 0.1 m (spec section 10, coupling)."""
    reference = simulator(body.ConstantBody(MASS_100MM), cadence=2.0).run()
    b = thermal_body(coarse_sphere_mesh)
    hist = coupled.CoupledRun(simulator(b), b, heating.SesamEquivalentHeating(), coupled.CoupledSettings(dt=2.0)).run()
    assert hist.end_reason == "ground" and reference.end_reason == "ground"
    assert abs(hist.columns["time_s"][-1] - reference.columns["time_s"][-1]) < 0.01
    t = hist.columns["time_s"]
    V_ref = np.interp(t, reference.columns["time_s"], reference.columns["velocity_kms"]) * 1e3
    h_ref = np.interp(t, reference.columns["time_s"], reference.columns["altitude_km"]) * 1e3
    assert np.abs(hist.columns["velocity_kms"] * 1e3 - V_ref).max() < 0.01
    assert np.abs(hist.columns["altitude_km"] * 1e3 - h_ref).max() < 0.1
    assert set(coupled.THERMAL_COLUMNS) <= set(hist.columns) and hist.columns["temperature_K"].max() > 500.0
    assert hist.results["end_reason"] == "ground" and hist.results["n_macro_steps"] == len(hist) - 1
    assert abs(hist.results["energy_balance_residual"]) < 1e-6 and hist.results["peak_surface_T_K"] > hist.results["peak_mean_T_K"]


def test_thermal_body_bookkeeping(coarse_sphere_mesh):
    b = thermal_body(coarse_sphere_mesh)
    assert b.mass(0.0) == MASS_100MM and b.mean_temperature() == pytest.approx(300.0) and b.field().shape == (coarse_sphere_mesh.n_nodes,)
    assert b.theta[b.i_stag] < 0.1 and b.theta[b.i_back] > 3.0 and b.energy() == pytest.approx(b.energy0)
    q = np.full(b.surface.n_patches, 1e5)
    loads = heating.HeatingResult(q, 1e5 / 0.27471, 0.0, 0.0, 0.0)
    E0 = b.energy()
    b.advance(0.5, 0.5, loads)
    Q = 1e5 * b.surface.area
    assert b.integrated_heat == pytest.approx(0.5 * Q) and b.last.Q_conv == pytest.approx(Q)
    assert b.absorbed_heat == pytest.approx(0.5 * (Q - b.last.Q_rad)) and b.radiated_heat == pytest.approx(0.5 * b.last.Q_rad)
    assert b.energy() - E0 == pytest.approx(b.absorbed_heat, rel=1e-5) and abs(b.energy_balance_residual()) < 1e-5
    stats = b.surface_stats()
    assert stats["surface_T_max_K"] > stats["T_centre_K"] > 299.9 and 300.0 < stats["surface_T_mean_K"] < stats["surface_T_max_K"] + 1e-9
    assert b.mean_temperature() > 300.0 and b.surface_temperature().shape == (b.surface.n_patches,)


def test_coupled_flight_with_the_thermal_body(coarse_sphere_mesh, tmp_path):
    """30 s of the 100 mm flight in SESAM-equivalent mode on the coarse mesh, VTK every 20 steps."""
    b = thermal_body(coarse_sphere_mesh)
    settings = coupled.CoupledSettings(dt=0.5, frames_every=20, output_dir=str(tmp_path / "vtk"))
    hist = coupled.CoupledRun(simulator(b, t_max=30.0), b, heating.SesamEquivalentHeating(), settings).run()
    c = hist.columns
    assert hist.end_reason == "t_max" and len(hist) == 61 and c["time_s"][-1] == 30.0
    assert c["convective_heat_W"][0] == pytest.approx(16244.0, rel=1e-2)        # SESAM's own t = 0 value is 16243.74 W
    assert np.all(np.diff(c["integrated_heat_J"]) > 0) and c["integrated_heat_J"][0] == 0.0
    assert c["integrated_heat_J"][-1] == pytest.approx(np.sum(c["convective_heat_W"][1:] * np.diff(c["time_s"])), rel=1e-9)
    assert np.all(c["rad_cooling_W"] < 0.0) and c["rad_cooling_W"][0] == pytest.approx(-0.4 * 5.670374419e-8 * b.surface.area * 300.0 ** 4, rel=1e-6)
    assert c["temperature_K"][-1] > c["temperature_K"][0] + 20.0 and c["surface_T_max_K"][-1] > c["temperature_K"][-1] > c["T_centre_K"][-1]
    assert np.all(c["heating_blend_f"] == [aero.SesamHeatTable()(kn) for kn in c["knudsen"]])
    assert hist.results["energy_balance_residual"] == pytest.approx(0.0, abs=1e-6)
    assert hist.results["integrated_heat_J"] == c["integrated_heat_J"][-1] and hist.results["n_frames"] == 4
    for k in range(4):
        assert os.path.isfile(tmp_path / "vtk" / "field_{}.vtu".format(k)) and os.path.isfile(tmp_path / "vtk" / "surface_{}.vtp".format(k))
    pvd = open(tmp_path / "vtk" / "field.pvd").read()
    assert 'timestep="10"' in pvd and 'file="field_1.vtu"' in pvd and os.path.isfile(tmp_path / "vtk" / "surface.pvd")


def test_frames_need_an_output_dir(coarse_sphere_mesh):
    with pytest.raises(ValueError):
        coupled.CoupledRun(simulator(body.ConstantBody(MASS_100MM)), body.ConstantBody(MASS_100MM),
                           heating.SesamEquivalentHeating(), coupled.CoupledSettings(frames_every=5))


# ---------------------------------------------------------------------------------------------------------------
# Step 3: the melting run, its columns, results, VTK fields and particle files

def test_coupled_melting_run_and_writers(coarse_sphere_mesh, tmp_path):
    """20 s of the 100 mm flight from 71 km with a warm body (700 K; physics heating, AA7075_range, Girin removal) on the
    coarse mesh: melt columns, results, demise bookkeeping, VTK melt fields and the three particle files."""
    pytest.importorskip("cantera")
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"), MASS_100MM, T0=700.0)
    sim = simulator(b, t_max=63.0)
    sim.advance(43.0)
    run_dir = str(tmp_path / "run")
    settings = coupled.CoupledSettings(dt=0.5, frames_every=20, output_dir=os.path.join(run_dir, "vtk"))
    run = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), settings)
    assert run.melting
    hist = run.run()
    c, r = hist.columns, hist.results
    assert set(coupled.MELT_COLUMNS) <= set(c) and hist.end_reason == "t_max"
    assert c["mass_kg"][-1] < c["mass_kg"][0] and c["sprayed_mass_kg"][-1] > 0.0 and np.all(np.diff(c["sprayed_mass_kg"]) >= 0.0)
    assert c["removed_mass_kg"][-1] == pytest.approx(c["sprayed_mass_kg"][-1]) and c["n_dead_elements"][-1] > 0
    assert c["mass_kg"][-1] == pytest.approx(c["mass_kg"][0] - c["sprayed_mass_kg"][-1], rel=1e-6)
    assert r["melt_onset_altitude_km"] is not None and r["spraying_onset_time_s"] >= r["melt_onset_time_s"] and r["demise_time_s"] is None
    assert r["sprayed_mass_kg"] == c["sprayed_mass_kg"][-1] and r["n_source_rows"] == len(b.source_rows) > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
    # the molten cascade (amendment of 2026-10-03): its cumulative mass and its passes per step, in the history and results
    assert r["cascade_mass_kg"] == c["cascade_mass_kg"][-1] and np.all(np.diff(c["cascade_mass_kg"]) >= 0.0)
    assert r["cascade_passes_max"] == int(c["cascade_passes"].max()) and r["cascade_capped_steps"] == b.cascade_capped_steps
    assert np.all(c["convective_heat_W"] > 0.0) and c["convective_heat_W"][-1] == pytest.approx(b.last.Q_conv)
    import pyvista as pv
    grid = pv.read(os.path.join(run_dir, "vtk", "field_1.vtu"))
    assert "liquid_fraction" in grid.point_data and "phi" in grid.cell_data and grid.n_cells == int(c["n_active_elements"][20])
    poly = pv.read(os.path.join(run_dir, "vtk", "surface_1.vtp"))
    for key in ("film_thickness", "we_s", "closure", "kn_local", "p_w", "tau", "r_droplet", "release_rate", "delta_m",
                "deep_thickness", "nonrigid_depth", "n_derived", "flow_eval"):
        assert key in poly.cell_data
    # the frame export (amendment of 2026-10-08): the derived surface's unit normals, and the record's flow evaluation on
    # the faces the frame's deaths exposed, so that every windward patch carries a wall pressure
    n_d = np.asarray(poly.cell_data["n_derived"])
    assert n_d.shape == (poly.n_cells, 3) and np.allclose(np.linalg.norm(n_d, axis=1), 1.0)
    ev = np.asarray(poly.cell_data["flow_eval"])
    assert set(np.unique(ev)) <= {1.0, 2.0} and (ev == 2.0).any() and (ev == 1.0).sum() > 0.5 * ev.size
    assert (np.asarray(poly.cell_data["p_w"])[n_d[:, 0] > 0.2] > 0.0).all()
    # the flow fields are the frame's own step's, carried across that step's element deaths: the wall pressure is positive
    # on the windward patches and nowhere above the stagnation value the history recorded from the same evaluation
    p_w = np.asarray(poly.cell_data["p_w"])
    assert c["n_dead_elements"][20] > c["n_dead_elements"][19]
    assert 0.0 < p_w.max() <= c["p_w_stag_Pa"][20] + 1e-9
    assert np.any(np.asarray(poly.cell_data["tau"]) != 0.0) and np.isfinite(np.asarray(poly.cell_data["kn_local"])).any()
    # ... and on the right patches: on the windward face the wall pressure falls with the angle from the flight direction
    from scipy.stats import spearmanr
    pts, tri = np.asarray(poly.points), np.asarray(poly.faces).reshape(-1, 4)[:, 1:]
    nrm = np.cross(pts[tri[:, 1]] - pts[tri[:, 0]], pts[tri[:, 2]] - pts[tri[:, 0]])
    nrm *= np.sign(np.einsum("ij,ij->i", nrm, pts[tri].mean(axis=1)))[:, None]          # outward on this near-sphere
    theta = np.arccos(np.clip(nrm[:, 0] / np.linalg.norm(nrm, axis=1), -1.0, 1.0))        # flight direction +x
    windward = (p_w > 0.0) & (theta < 0.5 * np.pi)
    assert spearmanr(theta[windward], p_w[windward])[0] < -0.9
    files = coupled.write_particles(run_dir, b, hist)
    for key in ("particles", "particles_summary", "size_distribution"):
        assert os.path.isfile(files[key])
    from reentry_model import spray
    npz = np.load(files["particles"])
    assert set(npz.files) == set(spray.SOURCE_COLUMNS) and npz["dm_kg"].sum() == pytest.approx(b.sprayed_mass, rel=1e-9)
    with open(files["size_distribution"]) as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "window_start_s,window_end_s,r_lo_m,r_hi_m,dn,dm_kg" and (len(lines) - 1) % spray.N_BINS == 0
    total = [l.split(",") for l in lines[1:]][-spray.N_BINS:]
    assert sum(float(row[5]) for row in total) == pytest.approx(b.sprayed_mass, rel=1e-6)
    with open(files["particles_summary"]) as fh:
        assert fh.readline().startswith("time_s,altitude_km,velocity_kms,released_mass_kg")


def test_surface_frames_carry_the_conjugate_depth_per_patch(coarse_sphere_mesh, tmp_path):
    """Girin's conjugate depth is the boundary between the skin the shear strips and the liquid that only runs off, and
    the large-fragment model that reads these frames needs it patch by patch (amendment of 2026-10-02). Three seconds
    from 69.8 km, where the windward face is under Girin's closure: the frame carries the step's own delta_m -- carried
    across the step's element deaths like p_w and tau, nan wherever the closure is not Girin's, where no conjugate depth
    exists -- and the thickness of the deep liquid on each patch. The body starts at 880 K so that it melts at once:
    the closure field is the spray step's, written only where there is film, while delta_m is written from every step
    that evaluates the flow."""
    pytest.importorskip("cantera")
    import pyvista as pv
    from reentry_model import spray, surface_flow as sf
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"), MASS_100MM, T0=880.0)
    sim = simulator(b, t_max=53.0)
    sim.advance(50.0)
    out = os.path.join(str(tmp_path), "vtk")
    hist = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), coupled.CoupledSettings(dt=0.5, frames_every=2, output_dir=out)).run()
    assert hist.end_reason == "t_max" and hist.results["n_frames"] == 4               # steps 0, 2, 4 and 6
    first, last = pv.read(os.path.join(out, "surface_0.vtp")), pv.read(os.path.join(out, "surface_3.vtp"))
    assert np.isnan(np.asarray(first.cell_data["delta_m"])).all()        # before any step there is no evaluation yet
    dm, closure = np.asarray(last.cell_data["delta_m"]), np.asarray(last.cell_data["closure"])
    girin = closure == sf.CLOSURE_GIRIN
    assert girin.sum() > 100 and (~girin).any()                           # windward patches with it, the lee without
    assert np.isfinite(dm[girin]).all() and (dm[girin] > 1e-4).all() and (dm[girin] < 1e-3).all()   # 0.1-1 mm
    assert np.isnan(dm[~girin]).all()
    expected = b.on_current_surface(spray.melt_layer(b.last_flow, b.liquid)[0], np.nan)   # the step's own, carried
    covered, record = b.flow_for_the_record()                             # and on the faces the deaths exposed, the record's
    if record is not None:                                                # (amendment of 2026-10-08 (frame export))
        expected = np.where(covered, expected, spray.melt_layer(record, b.liquid)[0])
    np.testing.assert_array_equal(dm, expected)
    assert (np.asarray(last.cell_data["deep_thickness"]) >= 0.0).all()
    # the non-rigid depth the regime test read (amendment of 2026-10-06): the spray step's own, carried; nan where none
    nr = np.asarray(last.cell_data["nonrigid_depth"])
    np.testing.assert_array_equal(nr, b.on_current_surface(b.last_nonrigid, np.nan))
    assert np.isfinite(nr).any() and (nr[np.isfinite(nr)] >= 0.0).all()


def test_demise_ends_the_run(coarse_sphere_mesh):
    """The lumped instant-removal device with a 60 % demise fraction: the loop stops with end_reason demise."""
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    mat = material.Material.from_drama_json("AA7075")
    mat.k_table = mat.k_table * 1e4
    b = body.MeltingBody(m, mat, thermal.thermal_solver("skfem"), MASS_100MM,
                         settings=body.MeltSettings(removal="instant", runoff=False, demise_fraction=0.6, size_feedback="initial"))
    hist = coupled.CoupledRun(simulator(b), b, heating.SesamEquivalentHeating(), coupled.CoupledSettings(dt=0.5)).run()
    assert hist.end_reason == "demise" and hist.results["end_reason"] == "demise" and hist.results["demise_altitude_km"] < 71.0
    assert 0.55 * MASS_100MM < hist.columns["mass_kg"][-1] < 0.6 * MASS_100MM and hist.results["final_mass_kg"] == hist.columns["mass_kg"][-1]


def test_the_continuum_step_reproduces_the_switched_step_harness(coarse_sphere_mesh, tmp_path):
    """`CoupledSettings.dt_continuum` (amendment of 2026-10-07, the continuum step): the macro step shortens to
    dt_continuum from the first step whose aero state, the one handed to the body, has a body Knudsen number below
    `kn_switch`, and stays there. That is what the 2026-10-05 harness (`dtswitch.py`) did from outside the package, and
    the option reproduces it bit for bit: from 48 s at 0.5 s to the switch (48.5 s for this hot body on the coarse mesh,
    whose Knudsen number reads its shrinking size), then 0.05 s to 50.5 s, against the harness's rule applied by hand. Frames are written every `frames_every` steps until the switch and by flight time
    after it, every `frames_every` first steps, so a fine step does not write one frame per step."""
    pytest.importorskip("cantera")
    runs, seen = {}, {}
    for label in ("option", "harness"):
        np.random.seed(12345)
        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_scheil"), thermal.thermal_solver("skfem"), MASS_100MM, T0=880.0)
        sim = simulator(b, t_max=50.5)
        sim.advance(48.0)
        if label == "option":
            s = coupled.CoupledSettings(dt=0.5, dt_continuum=0.05, kn_switch=0.01, frames_every=1,
                                        output_dir=str(tmp_path / "vtk"))
            run = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), s)
        else:
            run = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), coupled.CoupledSettings(dt=0.5))
            advance = b.advance

            def switched(t, dt, loads, state=None):          # the harness's rule: after the step, on the state it was given
                advance(t, dt, loads, state)
                if "t" not in seen and state is not None and np.isfinite(state.kn) and state.kn < 0.01:
                    run.settings.dt, seen["t"] = 0.05, t
            b.advance = switched
        runs[label] = (run.run(), run)
    (h, run), (g, _) = runs["option"], runs["harness"]
    t_sw = seen["t"]
    assert h.results["dt_switch_time_s"] == t_sw and h.results["dt_continuum_s"] == 0.05 and 48.0 < t_sw < 50.5
    assert 0.0 < h.results["dt_switch_kn"] < 0.01 and h.results["dt_switch_altitude_km"] > 0.0
    t = h.columns["time_s"]
    fine, coarse = np.diff(t)[t[1:] > t_sw + 1e-9], np.diff(t)[t[1:] <= t_sw + 1e-9]
    # the last step is whatever reaches t_max: here a sliver of ~1e-13 s, the fine steps summing to just short of it
    assert np.allclose(coarse, 0.5) and np.allclose(fine[:-1], 0.05) and 0.0 < fine[-1] <= 0.05 + 1e-12
    assert set(h.columns) == set(g.columns)
    for k in h.columns:
        np.testing.assert_array_equal(h.columns[k], g.columns[k], err_msg=k)
    frame_times = [ft for ft, _ in run.settings.frames]
    assert frame_times == pytest.approx([48.0, 48.5, 49.0, 49.5, 50.0, 50.5])          # by step, then by flight time


def test_the_continuum_step_is_off_by_default_and_rejects_a_coarser_step():
    s = coupled.CoupledSettings()
    assert s.dt_continuum is None
    with pytest.raises(ValueError):
        coupled.CoupledSettings(dt=0.5, dt_continuum=1.0)
    with pytest.raises(ValueError):
        coupled.CoupledSettings(dt=0.5, dt_continuum=0.0)


# ---------------------------------------------------------------------------------------------------------------
# The frame export for the Spheral replay (amendment of 2026-10-08; Spheral spec 2026-10-02 sec. 6.1 and 7.2)


def test_frames_carry_the_wall_loads_on_steps_without_film(coarse_sphere_mesh, tmp_path):
    """Two seconds from 69.8 km of a cold body (300 K): nothing melts, so the spray step never runs and `last_flow` stays
    None -- yet the surface flow is evaluated at the top of every melt step, and the frames carry its closure, p_w, tau
    and kn_local, and the history its stagnation wall pressure, on every step after the first (frame 0 precedes any
    evaluation). No element dies, so every patch is the step's own evaluation (flow_eval 1), and n_derived is the
    derived surface's smoothed normal, unit and outward whatever the triangle's winding."""
    pytest.importorskip("cantera")
    import pyvista as pv
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"), MASS_100MM)
    sim = simulator(b, t_max=52.0)
    sim.advance(50.0)
    out = os.path.join(str(tmp_path), "vtk")
    hist = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), coupled.CoupledSettings(dt=0.5, frames_every=1, output_dir=out)).run()
    c = hist.columns
    assert b.last_flow is None and c["sprayed_mass_kg"][-1] == 0.0 and np.isnan(c["p_w_stag_Pa"]).all()
    assert np.isnan(c["p_w_stag_step_Pa"][0]) and (c["p_w_stag_step_Pa"][1:] > 0.0).all()
    first, last = pv.read(os.path.join(out, "surface_0.vtp")), pv.read(os.path.join(out, "surface_4.vtp"))
    assert not np.asarray(first.cell_data["p_w"]).any() and not np.asarray(first.cell_data["flow_eval"]).any()
    p_w, tau = np.asarray(last.cell_data["p_w"]), np.asarray(last.cell_data["tau"])
    np.testing.assert_array_equal(p_w, b.last_flow_step.p_w)                       # no deaths: the surface is the evaluated one
    np.testing.assert_array_equal(tau, b.last_flow_step.tau)
    np.testing.assert_array_equal(np.asarray(last.cell_data["closure"]), b.last_flow_step.closure.astype(float))
    np.testing.assert_array_equal(np.asarray(last.cell_data["kn_local"]), b.last_flow_step.kn_local)
    assert (p_w > 0.0).all() and p_w.max() == pytest.approx(c["p_w_stag_step_Pa"][-1], rel=0.05) and tau.max() > 0.0
    assert (np.asarray(last.cell_data["flow_eval"]) == 1.0).all()
    n_d = np.asarray(last.cell_data["n_derived"])
    np.testing.assert_array_equal(n_d, m.surface(derived=True).smoothed_normals())
    tri = np.asarray(last.faces).reshape(-1, 4)[:, 1:]
    centre = np.asarray(last.points)[tri].mean(axis=1)
    assert np.allclose(np.linalg.norm(n_d, axis=1), 1.0) and (np.einsum("ij,ij->i", n_d, centre) > 0.0).all()


def test_the_record_evaluation_changes_nothing_the_physics_reads(coarse_sphere_mesh, tmp_path):
    """The frame export is write-only: the same 6 s of a hot body (880 K, deaths every step) run with frames at every step
    and without frames gives identical history columns, so writing a frame -- and evaluating the flow for the faces the
    deaths exposed -- changes nothing the model integrates. Both runs seed numpy as the CLI does, without which pyamg's
    random starts alone move the altitude by 1e-12."""
    pytest.importorskip("cantera")

    def run(frames):
        np.random.seed(12345)                    # pyamg's random starts, seeded as every CLI run is (amendment of 2026-10-05)
        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"),
                             MASS_100MM, T0=880.0)
        sim = simulator(b, t_max=56.0)
        sim.advance(50.0)
        st = coupled.CoupledSettings(dt=0.5, frames_every=1 if frames else 0,
                                     output_dir=os.path.join(str(tmp_path), "vtk") if frames else None)
        return coupled.CoupledRun(sim, b, heating.PhysicsHeating(), st).run(), b

    (with_frames, b), (without, _) = run(True), run(False)
    import pyvista as pv
    vtk = os.path.join(str(tmp_path), "vtk")
    n_record = sum(int((np.asarray(pv.read(os.path.join(vtk, f)).cell_data["flow_eval"]) == 2.0).sum())
                   for f in os.listdir(vtk) if f.endswith(".vtp"))
    assert with_frames.columns["n_dead_elements"][-1] > 0 and n_record > 0     # the record evaluation did run
    for key, v in without.columns.items():
        np.testing.assert_array_equal(with_frames.columns[key], v, err_msg=key)

