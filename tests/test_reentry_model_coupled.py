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
