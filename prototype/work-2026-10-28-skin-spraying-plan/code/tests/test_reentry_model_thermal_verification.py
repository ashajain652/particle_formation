"""The Step 2 verification script's level-3 option (layered-skin design section 11.1, sub-plan 20 Task 10 Step 4): a
skin run's command line and label, the thresholds it is judged by, the verdicts, the progress watch and the summary that
puts a skin run beside its element run. No flight is run."""
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
import reentry_model_thermal_verification as tv  # noqa: E402
from reentry_model import coupled  # noqa: E402
from test_reentry_model_reference_thermal import THRESHOLDS  # noqa: E402

GROUND = {"end_reason": "ground", "energy_balance_residual": 1e-9}


def metrics(q=0.01, cont=0.01, hyp=0.01, dT=0.01, rad=0.05, ratio=0.9):
    return {"Q_conv": {"max": q, "continuum_rel_max": cont},
            "integrated_heat": {"rel_error_end_of_hypersonic": hyp, "rel_error_end": ratio - 1.0, "ratio_end": ratio},
            "temperature": {"dT_max_K": 5.0, "dT_rel_max": dT}, "radiated": {"max": rad}}


def test_the_thresholds_are_step_2s():
    assert tv.THRESHOLDS == THRESHOLDS


def test_a_skin_run_is_the_element_run_with_skins_and_its_own_label():
    elements = tv.argv_of("d100", "sesam", "out")
    skins = tv.argv_of("d100", "sesam", "out", surface="skin")
    assert elements[elements.index("--name") + 1] == "d100__sesam"
    assert skins[skins.index("--name") + 1] == "d100__sesam__skin"
    drop = lambda a: a[:a.index("--name")] + a[a.index("--name") + 2:]
    # each names its surface model, the skins being the model's default since sub-plan 22 Task 3
    assert elements[-2:] == ["--surface-model", "elements"] and drop(skins) == drop(elements)[:-1] + ["skin"]
    smoke = tv.argv_of("d050", "physics", "out", surface="skin", t_max=2.0)
    assert smoke[smoke.index("--name") + 1] == "d050__physics__skin__t2" and smoke[-2:] == ["--t-max", "2"]


def test_the_verdicts():
    assert tv.verdict("sesam", metrics(), GROUND) == "pass"
    assert tv.verdict("sesam", metrics(hyp=-0.031, dT=0.021), GROUND) == "FAIL (integrated_heat_rel_hypersonic, dT_rel_max)"
    assert tv.verdict("sesam", metrics(cont=float("nan")), GROUND) == "FAIL (Q_conv_continuum_rel_max)"
    assert tv.verdict("sesam", metrics(), {**GROUND, "energy_balance_residual": 2e-6}) == "FAIL (energy_balance_residual)"
    assert tv.verdict("sesam", metrics(), {**GROUND, "end_reason": "t_max"}) == "FAIL (ended by t_max)"
    assert tv.verdict("physics", metrics(q=0.5), GROUND) == "reported"
    assert tv.verdict("physics", metrics(ratio=1.2), GROUND) == "FAIL (heat ratio 1.200)"


def test_the_progress_watch_reads_the_rows_and_is_undone(tmp_path, monkeypatch):
    row = {"time_s": 1.0, "altitude_km": 70.0, "velocity_kms": 7.4, "temperature_K": 301.0, "surface_T_max_K": 305.0}
    monkeypatch.setattr(coupled.CoupledRun, "row", lambda self, t, y, a, loads: dict(row))
    plain = coupled.CoupledRun.row
    undo = tv.watch(str(tmp_path / "p.progress"), every=2)
    for _ in range(5):
        assert coupled.CoupledRun.row(None, 0.0, None, None, None) == row
    undo()
    assert coupled.CoupledRun.row is plain
    lines = (tmp_path / "p.progress").read_text().splitlines()
    assert len(lines) == 3 and lines[0].startswith("step      0") and lines[2].startswith("step      4")


def test_run_case_passes_the_skins_and_undoes_the_watch(tmp_path, monkeypatch):
    plain = coupled.CoupledRun.row
    seen = []
    monkeypatch.setattr(tv.cli, "main", lambda argv: seen.append(argv) or 0)
    assert tv.run_case("d050", "physics", str(tmp_path), False, "skin") == "d050__physics__skin"
    assert seen[0][-2:] == ["--surface-model", "skin"] and coupled.CoupledRun.row is plain
    monkeypatch.setattr(tv.cli, "main", lambda argv: 1)
    with pytest.raises(SystemExit):
        tv.run_case("d050", "physics", str(tmp_path), False, "skin")
    assert coupled.CoupledRun.row is plain


def test_the_skin_summary_puts_each_skin_run_beside_its_element_run(tmp_path, monkeypatch):
    results = {**GROUND, "peak_surface_T_K": 900.0, "time_of_peak_surface_T_s": 60.0, "runtime_s": 120.0,
               "n_macro_steps": 732, "mean_newton_iterations": 3.0}
    skins = {"skin_columns": 9417, "rebuilds": 0, "interface_repeats": 12, "jump_faults": 0, "wall_skins_s": 100.0,
             "wall_3d_s": 10.0}
    for label, res, tm in (("d100__sesam", results, metrics()), ("d100__sesam__skin", {**results, **skins}, metrics(dT=0.03))):
        with open(tmp_path / (label + ".json"), "w") as fh:
            json.dump({"comparison": {"thermal_metrics": tm}, "results": res}, fh)
    monkeypatch.setattr(tv.cli, "main", lambda argv: 0)
    assert tv.main(["--outdir", str(tmp_path), "--cases", "d100", "--modes", "sesam", "--surface-model", "skin"]) == 0
    rows = json.loads((tmp_path / "summary_skin.json").read_text())
    assert [(r["surface"], r["verdict"]) for r in rows] == [("elements", "pass"), ("skin", "FAIL (dT_rel_max)")]
    assert rows[1]["skins"].startswith("9417 columns") and rows[0]["skins"] == ""
    assert not (tmp_path / "summary.json").exists()
    monkeypatch.setattr(tv.cli, "main", lambda argv: pytest.fail("--summary-only ran a flight"))
    (tmp_path / "summary_skin.json").unlink()
    assert tv.main(["--outdir", str(tmp_path), "--surface-model", "skin", "--summary-only"]) == 0
    assert json.loads((tmp_path / "summary_skin.json").read_text()) == rows
