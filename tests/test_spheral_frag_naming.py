"""spheral_frag.naming: names encode the whole configuration and round-trip (M1 plan, Task 1; spec §4.4)."""
import itertools

import pytest

from spheral_frag import naming

FE_RUN = "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam_none_fem-physics"
PREP = "prep_" + FE_RUN + "_k00000-01190"


def base(**over):
    args = dict(prepared=PREP, mode="replay", frames=(51, 450), form="rz", dx_mm=1.1, brackets={}, seed=1)
    args.update(over)
    return args


def test_prepare_name():
    assert naming.prepare_name(FE_RUN, 0, 1190) == PREP
    assert naming.prepare_name(FE_RUN, 0, 1190, every=2) == PREP + "_every2"
    assert naming.parse_prepare_name(PREP) == {"fe_run": FE_RUN, "k0": 0, "k1": 1190, "every": 1}
    assert naming.parse_prepare_name(PREP + "_every2")["every"] == 2
    for bad in [(FE_RUN, 5, 4, 1), (FE_RUN, -1, 4, 1), (FE_RUN, 0, 4, 0), ("a b", 0, 1, 1), ("a/b", 0, 1, 1),
                ("a__b", 0, 1, 1)]:
        with pytest.raises(ValueError):
            naming.prepare_name(*bad)


def test_run_name_default_brackets():
    name = naming.run_name(**base())
    assert name == PREP + "__replay_k00051-00450_rz_dx1.10mm_seed1"
    assert naming.run_name(**base(brackets=dict(naming.BRACKET_DEFAULTS))) == name   # defaults add nothing
    assert naming.analyse_name(name) == name
    assert " " not in name and "/" not in name


def test_non_default_bracket_adds_exactly_its_suffix():
    plain = naming.run_name(**base())
    assert naming.run_name(**base(brackets={"film_limit_mm": 3.0})) == plain + "_filmlimit-3mm"
    assert naming.run_name(**base(brackets={"film_limit_mm": 2.5})) == plain + "_filmlimit-2.5mm"
    assert naming.run_name(**base(brackets={"min_particles": 50})) == plain + "_minparticles-50"
    # fixed key order whatever the dict's order
    both = naming.run_name(**base(brackets={"link_h": 2.0, "film_limit_mm": 3.0}))
    assert both == plain + "_filmlimit-3mm_linkh-2"
    with pytest.raises(ValueError):
        naming.run_name(**base(brackets={"no_such_bracket": 1.0}))


NON_DEFAULT = {"film_limit_mm": 3.0, "base_pressure_pct": 1.5, "separation_deg": 120.0, "lee_shear": 0.25,
               "min_particles": 50, "link_h": 2.0, "crossings": 20, "damping": 0.2, "scheme": "crksph",
               "strength": "constant:100000000", "rate_min": 1e-07, "sampler": "radial", "azimuth_deg": 90.0}


def configurations():
    """Distinct configurations, each differing from the base in one setting (the last two in several)."""
    yield base()
    for mode in naming.MODES[1:]:
        yield base(mode=mode)
    yield base(form="3d")
    yield base(frames=(0, 450))
    yield base(frames=(51, 451))
    yield base(dx_mm=2.2)
    yield base(dx_mm=0.55)
    yield base(dx_mm=0.275)            # more than two decimals: kept exactly
    yield base(seed=2)
    yield base(seed=12345)
    yield base(prepared=PREP + "_every2")
    yield base(prepared="prep_other_k00000-00010")
    for key, value in NON_DEFAULT.items():
        yield base(brackets={key: value})
    yield base(brackets={"lee_shear": -0.5, "base_pressure_pct": 1e-05})
    yield base(brackets=dict(NON_DEFAULT))


def test_parse_round_trips():
    for cfg in configurations():
        name = naming.run_name(**cfg)
        p = naming.parse_run_name(name)
        assert p["prepared"] == cfg["prepared"] and p["mode"] == cfg["mode"] and p["form"] == cfg["form"]
        assert p["frames"] == cfg["frames"] and p["dx_mm"] == cfg["dx_mm"] and p["seed"] == cfg["seed"]
        assert p["brackets"] == {**naming.BRACKET_DEFAULTS, **cfg["brackets"]}
        again = naming.run_name(p["prepared"], p["mode"], p["frames"], p["form"], p["dx_mm"], p["brackets"],
                                p["seed"])
        assert again == name
    p = naming.parse_run_name(naming.run_name(**base()))
    assert p["prepare"] == {"fe_run": FE_RUN, "k0": 0, "k1": 1190, "every": 1}


def test_names_differ_when_any_setting_differs():
    names = [naming.run_name(**cfg) for cfg in configurations()]
    assert len(set(names)) == len(names)
    for a, b in itertools.combinations(list(configurations()), 2):
        if a != b:
            assert naming.run_name(**a) != naming.run_name(**b), (a, b)


@pytest.mark.parametrize("bad", [
    "no_double_underscore",
    PREP + "__replay_k00051-00450_rz_dx1.10mm",                 # no seed
    PREP + "__walk_k00051-00450_rz_dx1.10mm_seed1",             # unknown mode
    PREP + "__replay_k00051-00450_2d_dx1.10mm_seed1",           # unknown form
    PREP + "__replay_k00051-00450_rz_dx1.1mm_seed1",            # not canonical
    PREP + "__replay_k00051-00450_rz_dx1.10mm_seed1_filmlimit-2mm",   # a default written out
    PREP + "__replay_k00051-00450_rz_dx1.10mm_seed1_linkh-2_filmlimit-3mm",   # out of order
    PREP + "__replay_k00051-00450_rz_dx1.10mm_seed1_filmlimit-3",     # unit missing
    PREP + "__replay_k00051-00450_rz_dx1.10mm_seed1_colour-3",        # unknown bracket
    PREP + "__replay_k00051-00450_rz_dx1.10mm_seed1_strength-constant-1e8",   # strength not canonical
    PREP + "__replay_k00051-00450_rz_dx1.10mm_seed1_scheme-SPH",      # upper case
    PREP + "__replay_k00051-00450_rz_dx1.10mm_seed1_crossings-2.5",   # not an integer
])
def test_parse_rejects(bad):
    with pytest.raises(ValueError):
        naming.parse_run_name(bad)
    with pytest.raises(ValueError):
        naming.analyse_name(bad)


def test_bad_arguments():
    for over in [dict(mode="walk"), dict(form="2d"), dict(dx_mm=0.0), dict(dx_mm=float("nan")), dict(seed=-1),
                 dict(seed=1.0), dict(frames=(5, 4)), dict(prepared="a__b"), dict(prepared="a b"),
                 dict(brackets={"min_particles": 2.5}), dict(brackets={"lee_shear": float("inf")}),
                 dict(brackets={"strength": "constant:-5"}), dict(brackets={"scheme": "a_b"}),
                 dict(brackets={"scheme": 3})]:
        with pytest.raises(ValueError):
            naming.run_name(**base(**over))


# ------------------------------------------------------------------------------------- directory layout (Task 7b)
LONG_PREP = naming.prepare_name("model_" + "x" * 181, 0, 1254)          # 205 characters
EVERY_BRACKET = dict(film_limit_mm=3.0, base_pressure_pct=4.0, separation_deg=170.0, lee_shear=0.25, min_particles=50,
                     link_h=2.0, crossings=20, damping=0.2, scheme="crksph", strength="constant:1.5e8", rate_min=1e-5,
                     sampler="radialdt20r2", azimuth_deg=22.5, loads="off", mass_scaling="off", monitor_every=5,
                     monitor_limit=0.025, geometry="barl50w10h1", loads_frame=182,
                     area_shares="volume")


def test_long_run_name_nests_and_parses_back(tmp_path):
    import os
    assert len(LONG_PREP) == 205 and set(EVERY_BRACKET) == set(naming.BRACKET_DEFAULTS)
    name = naming.run_name(LONG_PREP, "replay", (0, 1254), "3d", 1.1, EVERY_BRACKET, 1)
    assert len(name) > 2 * naming.NAME_MAX                     # whole: far beyond one component
    d = naming.run_dir(str(tmp_path / "runs"), name)
    rel = os.path.relpath(d, str(tmp_path / "runs")).split(os.sep)
    assert rel[0] == LONG_PREP and len(rel) == 3                # the 357-character suffix cut once
    assert all(len(p) <= naming.NAME_MAX for p in rel) and "_".join(rel[1:]) == name.split("__")[1]
    os.makedirs(d)
    assert naming.run_name_from_dir(d) == name
    assert naming.parse_run_name(naming.run_name_from_dir(d))["brackets"] == naming.canonical_brackets(EVERY_BRACKET)
    a = naming.analyse_dir(str(tmp_path / "analyse"), name)
    os.makedirs(a)
    assert naming.run_name_from_dir(a) == name
    # the usual case: two levels, the suffix whole
    short = naming.run_name(**base(brackets={"crossings": 20}))
    assert naming.run_dir("out", short) == os.path.join("out", PREP, short.split("__")[1])
    assert naming.run_name_from_dir(naming.run_dir("out", short)) == short


def test_long_body_name_nests_and_parses_back(tmp_path):
    import os
    name = naming.body_name(LONG_PREP, 0, 1254, "rz", 0.55, -22.5, "radialdt20r2", "barl50w10h1")
    assert len(name) > naming.NAME_MAX
    d = naming.body_dir(str(tmp_path), name)
    assert d == os.path.join(str(tmp_path), LONG_PREP, "k00000-01254_rz_dx0.55mm_azimuth--22.5deg"
                             "_synthetic-radialdt20r2_geometry-barl50w10h1")
    os.makedirs(d)
    assert naming.body_name_from_dir(d) == name and naming.parse_body_name(name)["prepared"] == LONG_PREP


def test_directory_names_reject():
    with pytest.raises(ValueError):
        naming.run_dir("out", "not_a_run")
    with pytest.raises(ValueError):
        naming.run_dir("out", naming.run_name("p" * 256, "replay", (0, 1), "rz", 1.1, {}, 1))   # prepared too long
    with pytest.raises(ValueError):
        naming.run_name_from_dir("out/prep_x/not_a_suffix")
    with pytest.raises(ValueError):
        naming.body_name_from_dir("out/prep_x/k00000-00001_xx_dx1.10mm")
