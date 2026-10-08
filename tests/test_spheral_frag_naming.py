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
               "min_particles": 50, "link_h": 2.0}


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
])
def test_parse_rejects(bad):
    with pytest.raises(ValueError):
        naming.parse_run_name(bad)
    with pytest.raises(ValueError):
        naming.analyse_name(bad)


def test_bad_arguments():
    for over in [dict(mode="walk"), dict(form="2d"), dict(dx_mm=0.0), dict(dx_mm=float("nan")), dict(seed=-1),
                 dict(seed=1.0), dict(frames=(5, 4)), dict(prepared="a__b"), dict(prepared="a b"),
                 dict(brackets={"min_particles": 2.5}), dict(brackets={"lee_shear": float("inf")})]:
        with pytest.raises(ValueError):
            naming.run_name(**base(**over))
