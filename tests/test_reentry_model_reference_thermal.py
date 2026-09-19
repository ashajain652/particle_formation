"""Coupled model vs the two no-melt US76 SESAM references (marker: reference, ~25 min): SESAM-equivalent heating on
the default mesh against the acceptance thresholds of spec section 10, physics mode reported (not thresholded),
and the mesh / time-step refinement checks on the 100 mm case. Each run's metrics are also written to
reentry_model_output/verification_thermal/ for the README table."""
import json
import os

import numpy as np
import pytest

from reentry_model import aero, atmosphere, body, compare, coupled, heating, material, mesh, sesam_io, thermal
from reentry_model import trajectory as tj

NAMES = {
    "d100": "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_nowind",
    "d050": "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_mAA7075_nomelt_nowind",
}
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reentry_model_output", "verification_thermal")
# Spec section 10 expectations, except the radiated power: the resolved surface radiates at its own (hotter)
# temperature, so a 2 % temperature margin is a 4 x 2 % = 8 % margin on eps sigma T^4 (measured 6.7 % / 3.4 % on
# 2026-09-18 while T_eq was within 1.2 %). Task 12 may revise a value only together with the measured number and
# the reason, recorded in the README verification table.
THRESHOLDS = {"Q_conv_peak_rel_max": 0.03, "Q_conv_continuum_rel_max": 0.03, "integrated_heat_rel_hypersonic": 0.03,
              "dT_rel_max": 0.02, "radiated_peak_rel_max": 0.08}


def load(key):
    return sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, NAMES[key] + ".csv"))


def coupled_run(ref, heating_name, h_surface=mesh.DEFAULT_H_SURFACE, h_core=mesh.DEFAULT_H_CORE, dt=0.5, t_max=3600.0):
    initial = tj.InitialState(ref.initial.velocity, ref.initial.altitude, ref.initial.flight_path, ref.initial.heading,
                              ref.initial.lat, ref.initial.lon, ref.initial.epoch)
    settings = tj.Settings(diameter=ref.diameter, t_max=t_max)
    the_mesh = mesh.sphere_mesh(ref.diameter / 2.0, h_surface, h_core)
    the_body = body.ThermalBody(the_mesh, material.Material.from_drama_json(), thermal.thermal_solver("skfem"),
                                body.sphere_mass(ref.diameter, ref.material_density))
    model = heating.SesamEquivalentHeating() if heating_name == "sesam" else heating.PhysicsHeating()
    sim = tj.Simulator(initial, the_body, atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(), aero.SesamTable(), settings)
    return coupled.CoupledRun(sim, the_body, model, coupled.CoupledSettings(dt=dt)).run()


def record(label, history, ref):
    os.makedirs(OUTDIR, exist_ok=True)
    doc = {"case": ref.name, "results": history.results, "thermal_metrics": compare.thermal_metrics(history, ref),
           "trajectory_metrics": compare.metrics(history, ref)}
    with open(os.path.join(OUTDIR, label + ".json"), "w") as fh:
        json.dump(doc, fh, indent=2, default=str)
    return doc["thermal_metrics"]


@pytest.mark.reference
@pytest.mark.parametrize("key", list(NAMES))
def test_sesam_equivalent_mode_matches_sesam(key):
    ref = load(key)
    history = coupled_run(ref, "sesam")
    assert history.end_reason == "ground" and abs(history.results["energy_balance_residual"]) < 1e-6
    m = record(key + "__sesam", history, ref)
    assert m["Q_conv"]["max"] <= THRESHOLDS["Q_conv_peak_rel_max"], m["Q_conv"]
    assert m["Q_conv"]["continuum_rel_max"] <= THRESHOLDS["Q_conv_continuum_rel_max"], m["Q_conv"]
    assert abs(m["integrated_heat"]["rel_error_end_of_hypersonic"]) <= THRESHOLDS["integrated_heat_rel_hypersonic"], m["integrated_heat"]
    assert m["temperature"]["dT_rel_max"] <= THRESHOLDS["dT_rel_max"], m["temperature"]
    assert m["radiated"]["max"] <= THRESHOLDS["radiated_peak_rel_max"], m["radiated"]


@pytest.mark.reference
@pytest.mark.parametrize("key", list(NAMES))
def test_physics_mode_is_reported(key):
    """Lees x Fay-Riddell delivers less heat than SESAM's tumbling average: the ratio is recorded, not thresholded
    (expected 0.6-0.9: 0.196/0.2747 x Fay-Riddell/DKR x hot wall)."""
    ref = load(key)
    history = coupled_run(ref, "physics")
    m = record(key + "__physics", history, ref)
    assert history.end_reason == "ground" and 0.4 < m["integrated_heat"]["ratio_end"] < 1.1


@pytest.mark.reference
def test_mesh_and_time_step_refinement():
    """100 mm, SESAM-equivalent, to 200 s (past peak heating and peak surface temperature): halving h_surface
    (1 mm / 8 mm, 76 k nodes) changes the surface-temperature history by < 1 %; halving dt changes the peak surface
    temperature by < 0.5 %."""
    ref = load("d100")
    base = coupled_run(ref, "sesam", t_max=200.0)
    fine = coupled_run(ref, "sesam", h_surface=0.5 * mesh.DEFAULT_H_SURFACE, t_max=200.0)
    t = base.columns["time_s"]
    for col in ("surface_T_max_K", "surface_T_mean_K", "temperature_K"):
        d = np.abs(np.interp(t, fine.columns["time_s"], fine.columns[col]) - base.columns[col]) / base.columns[col]
        assert d.max() < 0.01, (col, d.max())
    half = coupled_run(ref, "sesam", dt=0.25, t_max=200.0)
    assert abs(half.results["peak_surface_T_K"] / base.results["peak_surface_T_K"] - 1.0) < 0.005
    with open(os.path.join(OUTDIR, "d100__refinement.json"), "w") as fh:
        json.dump({"base": base.results, "h_surface_halved": fine.results, "dt_halved": half.results}, fh, indent=2, default=str)
