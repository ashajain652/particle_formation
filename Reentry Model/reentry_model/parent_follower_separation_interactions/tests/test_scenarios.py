"""Scenarios, run names and case output (Task 8, spec 4.3b-4.4): a sphere in the spacecraft's plume, released by the
follower's own regime, and a molten droplet in a solid sphere's plume."""
import json, math, numpy as np, pytest
from reentry_model.parent_follower_separation_interactions.scenarios import (ScenarioSettings, run_name, sphere_in_spacecraft_plume,
    droplet_in_sphere_plume, write_case, classify_release)
from reentry_model.parent_follower_separation_interactions.laws import LawSettings

def test_run_names_differ_for_every_setting():
    base = run_name("sphere", d_mm=50, h_km=80, settings=ScenarioSettings())
    others = [ScenarioSettings(bound="high"), ScenarioSettings(wake_flow="both"),
              ScenarioSettings(law=LawSettings(taylor_K=0.5)), ScenarioSettings(frontal_area=3.1),
              ScenarioSettings(beta_scale=1.1), ScenarioSettings(density=2700.0)]
    names = {run_name("sphere", d_mm=50, h_km=80, settings=s) for s in others}
    assert base not in names and len(names) == len(others)
    assert run_name("sphere", 50, 80, ScenarioSettings(), release="escape", theta_deg=160.0) != base

def test_far_wake_case_runs_and_writes(twin, tmp_path):
    # The brief released at theta = 90 degrees. There the shoulder gas (Newtonian pressure p_inf, expanded from the
    # stagnation state) is hot and thin: its mean free path at 60 km is 2.4 mm (measured), so a 10 mm follower has
    # Kn_f = 0.24 and is transitional, and spec 4.3c sends it to the surface path, not the far wake. At 120 degrees the
    # mean free path is 0.059 mm, Kn_f = 0.0059 (continuum), and AFFIx's zone is still the far wake (theta_bif 139.4).
    c = sphere_in_spacecraft_plume(10, 60, 120.0, ScenarioSettings(wake_flow="both"), twin=twin)
    csv_path, json_path = write_case(c, tmp_path)
    s = json.load(open(json_path))
    assert s["affix_zone"] == "far_wake" and s["t_fz"] > 0
    assert s["reason"] in ("ok", "censored", "reached_parent") and s["residence_time"] >= 0
    assert s["parent_sha256"] and "git_commit" in s and s["regime_at_release"] and s["follower_t_end"] > 0
    h = c.history                                      # starts on the axis, two parent diameters behind
    gap = np.linalg.norm([h[0]["follower_x"] - h[0]["parent_x"], h[0]["follower_y"] - h[0]["parent_y"],
                          h[0]["follower_z"] - h[0]["parent_z"]])
    assert gap == pytest.approx(2 * 2 * math.sqrt(6.215 / math.pi), rel=1e-6)

def test_follower_regime_thresholds():
    from reentry_model.parent_follower_separation_interactions.scenarios import follower_regime, AFFIX_STATUS
    assert [follower_regime(k, 0.01, 10.0) for k in (0.005, 0.5, 10.0)] == ["continuum", "transitional",
                                                                             "free_molecular"]
    assert [AFFIX_STATUS[r] for r in ("continuum", "transitional", "free_molecular")] == [
        "governing", "indicative", "not_applicable"]

def test_rarefied_follower_starts_on_the_surface_and_reports_the_shock_layer(twin):
    c = sphere_in_spacecraft_plume(1, 90, 150.0, twin=twin)
    s = c.summary
    assert s["follower_regime"] in ("transitional", "free_molecular") and s["release_kind"] == "surface"
    assert s["affix_status"] in ("indicative", "not_applicable")
    assert s["kn_f_release"] < s["kn_f_free_stream"]                     # the shock-layer gas is denser
    assert s["shock_layer_exit"] in ("through_shock", "past_base_plane", "hit_parent", "inside_at_end")
    h = c.history
    gap = np.linalg.norm([h[0]["follower_x"] - h[0]["parent_x"], h[0]["follower_y"] - h[0]["parent_y"],
                          h[0]["follower_z"] - h[0]["parent_z"]])
    assert gap == pytest.approx(math.sqrt(6.215 / math.pi), rel=1e-6)   # on the surface, D1/2 from the centroid

def test_escape_and_shock_band_and_windward_limit(twin):
    esc = sphere_in_spacecraft_plume(10, 60, 180.0, twin=twin)
    assert esc.summary["follower_regime"] == "continuum" and esc.summary["affix_status"] == "governing"
    assert esc.summary["affix_zone"] == "escape" and esc.summary["residence_time"] < 0.05
    tb = classify_release(10, 60, 90.0, twin, ScenarioSettings()).theta_bif
    band = sphere_in_spacecraft_plume(10, 60, tb, twin=twin)
    assert band.summary["reason"] == "shock_band" and band.residence is None
    with pytest.raises(ValueError):
        sphere_in_spacecraft_plume(10, 60, 45.0, twin=twin)

def test_droplet_in_sphere_plume_worked_examples(twin):
    # controller ruling R2: (50, 80, 100) may be rarefied (no law) -- assert completion and the recorded regime;
    # (100, 60, 50) is a continuum, laminar sphere parent and must give a finite residence with validity records
    c = droplet_in_sphere_plume(50, 80, 100, twin=twin)
    assert c.summary["regime_at_release"] and c.summary["kn_at_release"] > 0
    c = droplet_in_sphere_plume(100, 60, 50, twin=twin)
    assert c.summary["reason"] != "skipped_weber"
    assert math.isfinite(c.residence.residence_time) and c.residence.residence_time >= 0 and c.residence.validity
    h = c.history
    assert np.allclose([h[0]["follower_x"], h[0]["follower_y"], h[0]["follower_z"]],
                       [h[0]["parent_x"], h[0]["parent_y"], h[0]["parent_z"]], atol=1e-3)


# ---- beyond the brief's tests: the release-state conversion, the surface-release record and the file layout ----

def test_state_from_ecef_round_trips_through_the_simulator(twin):
    from reentry_model.atmosphere import US76TableAtmosphere
    from reentry_model.parent_follower_separation_interactions.adapters import state_from_ecef, sphere_simulator
    t = twin.t_at_altitude(70.0)
    r, v = twin.track.position(t), twin.track.velocity(t)
    y0 = sphere_simulator(state_from_ecef(r, v, twin.epoch), 0.01, US76TableAtmosphere()).initial_state_vector()
    assert np.allclose(y0[:3], r, rtol=0, atol=1e-6) and np.allclose(y0[3:], v, rtol=0, atol=1e-9)

def test_surface_release_summary_carries_the_shock_layer_record(twin):
    s = sphere_in_spacecraft_plume(1, 80, 120.0, twin=twin).summary      # Kn_f about 0.47: transitional
    for key in ("flux_mismatch", "flux_iterations", "flux_mismatch_above_tol", "mach_outside_billig",
                "shock_layer_nonequilibrium", "kn_f_min_in_layer", "follower_continuum_in_layer", "epsilon",
                "standoff_m", "standoff_billig_m", "q_change", "q_change_large", "cd_composed_release",
                "cd_sesam_release", "a_p_max", "a_drag_max_shock_layer", "t_shock_layer_exit"):
        assert key in s, key
    assert s["flux_mismatch"] < 0.02 and not s["flux_mismatch_above_tol"]
    assert s["transitional_follower"] == (s["follower_regime"] == "transitional")
    assert s["kn_f_min_in_layer"] <= s["kn_f_release"] * (1 + 1e-6)     # the release point is inside the layer
    assert s["follower_continuum_in_layer"] == (s["kn_f_min_in_layer"] < 0.01)
    assert s["cd_composed_release"] > 0 and s["cd_sesam_release"] > 0

def test_case_without_integration_writes_only_the_summary(twin, tmp_path):
    tb = classify_release(10, 60, 90.0, twin, ScenarioSettings()).theta_bif
    c = sphere_in_spacecraft_plume(10, 60, tb, twin=twin)
    csv_path, json_path = write_case(c, tmp_path)
    assert csv_path is None and json_path == str(tmp_path / "runs" / (c.name + ".json"))
    s = json.load(open(json_path))
    assert s["residence_time"] is None and s["follower_regime"] == "continuum"   # NaN is written as null

def test_a_twin_built_with_other_parent_settings_is_refused(twin):
    with pytest.raises(ValueError, match="frontal_area"):
        sphere_in_spacecraft_plume(10, 60, 180.0, ScenarioSettings(frontal_area=3.1), twin=twin)
