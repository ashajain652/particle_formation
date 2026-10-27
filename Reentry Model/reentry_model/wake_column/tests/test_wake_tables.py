import math
import warnings
import pytest
from reentry_model.wake_column import species as sp, reactions as rx
from reentry_model.wake_column.settings import Settings, Case, run_name

CM3, CM6 = 1e-6, 1e-12

def test_species_table_has_every_network_species():
    s = sp.load_species()
    assert set(sp.AL_SPECIES) | set(sp.AMBIENT_SPECIES) <= set(s)
    assert s["AlO"].dipole_D == pytest.approx(4.45) and math.isnan(s["OAlO"].dipole_D)
    assert all(s[n].n_Al == 1 for n in sp.AL_SPECIES)

@pytest.mark.parametrize("T", [200.0, 300.0])
def test_plane_table_values_reproduced(T):            # spec test 8
    r = {x.rid: x for x in rx.load_reactions()}
    assert r["R1"].k(min(T, 295.0)) == pytest.approx(1.7e-10 * (min(T, 295.0) / 300) ** (1 / 6) * CM3, rel=1e-9)
    assert r["R11"].k(T) == pytest.approx(2e-10 * (T / 300) ** (1 / 6) * CM3, rel=1e-9)
    lg = math.log10(T)
    assert r["R6"].k(T) == pytest.approx(10 ** (-35.137 + 6.1052 * lg - 1.4089 * lg * lg) * CM6, rel=1e-9)
    assert r["R17"].k(T) == pytest.approx(3.3e-3) and r["R17"].day_only

def test_reverse_of_G42_is_not_duplicated():
    rids = [x.rid for x in rx.load_reactions()]
    assert "R3" in rids and "G42" in rids and len(rids) == len(set(rids))
    assert not any(x.reactants == ("OAlO", "O") for x in rx.load_reactions() if x.rid != "R3")

def test_hot_only_dissociation_not_in_box_network():
    assert "D1" not in [x.rid for x in rx.load_reactions()]
    assert "D1" in [x.rid for x in rx.load_reactions(stage=("hot", "both"))]

def test_clipping_warns_once_per_reaction():          # Review Focus 5
    r = {x.rid: x for x in rx.load_reactions()}["R6"]
    w = rx.RateWarnings()
    for _ in range(1000):
        assert r.k(50.0, w) == r.k(r.T_min)
    assert len(w.messages) == 1

def test_run_name_encodes_non_default_settings():
    c = Case(14.1, 76.0, "day")
    assert run_name(c, Settings()) == "wake_d014.10um_h076.0km_day_waccm"
    assert run_name(c, Settings(formation_scale=0.1, hot_start="droplet_temperature")).endswith("_form0.1_hotdroplet")
    assert run_name(c, Settings(ambient="msis_afgl")) != run_name(c, Settings())
