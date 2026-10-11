"""Melting model vs the two melting US76 SESAM references (marker: reference, ~2 min): the bookkeeping device
(SESAM-equivalent heating, AA7075, instant removal, k x 1e4, D0/R0 kept, the Step 2 default mesh) against the acceptance
thresholds of spec section 13.1 -- mass within 2 % of the initial at every reference time, melt-onset altitude within
0.5 km, the 1 %-mass time within 2 % (measured 2026-09-21: 0.98 % / +0.10 km / -1.3 % for 100 mm, 1.29 % / +0.17 km /
-0.2 % for 50 mm). The resolved and physics-mode runs are analysis/melt_verification.py's business (reported). Metrics
are written to reentry_model_output/verification_melt/reference_tests/."""
import json
import os

import pytest

from reentry_model import aero, atmosphere, body, compare, coupled, heating, material, mesh, sesam_io, thermal
from reentry_model import trajectory as tj

NAMES = {"d100": "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind",
         "d050": "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind"}
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reentry_model_output", "verification_melt", "reference_tests")
THRESHOLDS = {"mass_rel_m0": 0.02, "onset_km": 0.5, "demise_time_rel": 0.02}


def bookkeeping_run(ref):
    initial = tj.InitialState(ref.initial.velocity, ref.initial.altitude, ref.initial.flight_path, ref.initial.heading,
                              ref.initial.lat, ref.initial.lon, ref.initial.epoch)
    the_mesh = mesh.sphere_mesh(ref.diameter / 2.0)
    mat = material.Material.from_drama_json("AA7075")
    mat.k_table = mat.k_table * 1e4
    the_body = body.MeltingBody(the_mesh, mat, thermal.thermal_solver("skfem"), body.sphere_mass(ref.diameter, ref.material_density),
                                settings=body.MeltSettings(removal="instant", runoff=False, size_feedback="initial"))
    sim = tj.Simulator(initial, the_body, atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(), aero.SesamTable(),
                       tj.Settings(diameter=ref.diameter))
    return coupled.CoupledRun(sim, the_body, heating.SesamEquivalentHeating(), coupled.CoupledSettings(dt=0.5)).run()


@pytest.mark.reference
@pytest.mark.parametrize("key", ["d100", "d050"])
def test_bookkeeping_device_follows_sesam(key):
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, NAMES[key] + ".csv"))
    hist = bookkeeping_run(ref)
    mm = compare.melt_metrics(hist, ref)
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, key + "__bookkeeping.json"), "w") as fh:
        json.dump({"case": ref.name, "results": hist.results, "melt_metrics": mm}, fh, indent=2, default=str)
    assert hist.end_reason == "demise" and abs(hist.results["melt_energy_balance_residual"]) < 1e-6
    assert mm["mass"]["max_rel_m0"] <= THRESHOLDS["mass_rel_m0"]
    assert abs(mm["onset_altitude_diff_km"]) <= THRESHOLDS["onset_km"]
    assert abs(mm["demise_time_rel"]) <= THRESHOLDS["demise_time_rel"]
