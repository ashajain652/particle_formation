

# ---------------------------------------------------------------------------------------------------------------
# Step 3 sub-plan 18: the melt-layer skins' flags, run names, columns and frame fields

def test_skin_settings_reach_the_run_name():
    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin")
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
                              surface_model="skin") == name + "_surface-skin"
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin", surface_model="skin",
                              skin_suffix="_skin-dev", seed=7) == name + "_surface-skin_skin-dev_seed-7"
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", None,
                              surface_model="skin") == cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table",
                                                                          "none", "physics", None)   # no melt, no skins
    assert cli.skin_name_suffix(0.4, None, 10.0, "production") == "" and cli.skin_name_suffix(0.4, 10.0, 10.0, "production") == ""
    assert cli.skin_name_suffix(0.4, None, 10.0, "dev") == "_skin-dev"
    assert cli.skin_name_suffix(0.8, 5.0, 20.0, "production") == "_skin0.8mm_cell5um_sub20ms"


def test_a_skin_run_writes_its_columns_frames_and_settings(tmp_path):
    """3 s from 71 km with a warm body and skins (dev preset), the 0.5 s step throughout (`--dt-continuum off`, as the
    short melting run passes its own step): exit 0, the run name, the skin columns, the frame fields."""
    pytest.importorskip("cantera")
    import glob
    import pyvista as pv
    argv = [a for a in MELT if a not in ("--heating", "sesam")] + [
        "--heating", "physics", "--altitude", "71", "--velocity", "7.24", "--temperature", "800", "--t-max", "3",
        "--dt-continuum", "off", "--surface-model", "skin", "--skin-preset", "dev", "--frames-every", "2", "--outdir", str(tmp_path)]
    assert cli.main(argv) == 0
    (path,) = glob.glob(str(tmp_path / "*.json"))
    doc = json.load(open(path))
    assert doc["run_name"].endswith("_melt-girin_dtcontinuum-off_surface-skin_skin-dev")
    s, r = doc["settings"], doc["results"]
    assert s["surface_model"] == "skin" and s["skin"] == {"thickness_mm": 0.4, "cell_um": 20.0, "substep_ms": 10.0, "preset": "dev"}
    assert r["skins_created"] > 0 and r["skin_drawn_mass_kg"] >= 0.0 and abs(r["melt_energy_balance_residual"]) < 1e-6
    rows = list(csv.DictReader(open(path[:-5] + ".csv")))
    assert all(k in rows[0] for k in coupled.SKIN_COLUMNS) and float(rows[-1]["skin_count"]) > 0.0
    surf = sorted(glob.glob(str(tmp_path / "**" / "surface_*.vtp"), recursive=True))
    assert surf and all(k in pv.read(surf[-1]).cell_data for k in ("skin_T_top", "skin_liquid_depth", "skin_nonrigid_depth",
                                                                    "skin_thickness"))


@pytest.mark.parametrize("argv", [
    MELT + ["--surface-model", "magic"],
    MELT + ["--surface-model", "skin", "--skin-thickness", "0"],
    MELT + ["--surface-model", "skin", "--skin-substep", "-1"],
    MELT + ["--surface-model", "skin", "--skin-cell", "300"],          # more than half of the 0.4 mm skin
    MELT + ["--surface-model", "skin", "--skin-preset", "fine"],
])
def test_bad_skin_arguments_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2
