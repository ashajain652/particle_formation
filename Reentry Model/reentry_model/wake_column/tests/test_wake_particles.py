import dataclasses
import math

import numpy as np
import pytest

from reentry_model.wake_column import particles as pa, species as sp


def test_first_bin_is_a_half_nanometre_alumina_dimer():
    b = pa.Bins()
    assert b.n_atoms[0] == 2 and b.n_atoms[1] == 4
    assert pa.diameter(2, 0.0) == pytest.approx(0.4759e-9, rel=1e-3)
    assert pa.diameter(b.n_atoms[-1], 0.0) >= 400e-9 > pa.diameter(b.n_atoms[-2], 0.0)


@pytest.mark.parametrize("n", [2.0, 3.0, 5.0, 4095.0])
def test_split_conserves_number_and_atoms(n):
    b = pa.Bins(); k, z = b.split(n)
    assert z * b.n_atoms[k] + (1 - z) * b.n_atoms[k + 1] == pytest.approx(n)
    assert pa.Bins().atom_fraction(n, k, z) == pytest.approx(z * b.n_atoms[k] / n)


def test_hydroxide_is_bulkier_per_aluminium():
    assert pa.volume_per_al(3.0) / pa.volume_per_al(0.0) == pytest.approx(5.3524 / 2.8218, rel=1e-3)


def test_epstein_speed_hand_value():                        # r = 10 nm alumina at about 80 km
    assert pa.epstein_speed(1e-8, 3000.0, 1.846e-5, 198.6, 9.565) == pytest.approx(0.02929, rel=1e-3)


def test_dipole_capture_AlO_AlO_at_298K():                  # decision 7; Gomez Martin 2017 quote 1.1e-9 cm3/s
    s = sp.load_species()
    a = dataclasses.replace(s["AlO"], dipole_D=4.36)
    assert pa.capture_rate(a, a, 298.0) == pytest.approx(1.516e-15, rel=1e-3)


def test_unknown_dipole_falls_back_to_hard_sphere():
    s = sp.load_species(); a, b = s["OAlO"], s["AlO"]
    mu = a.mass_amu * b.mass_amu / (a.mass_amu + b.mass_amu) * 1.66053906660e-27
    hs = math.pi * (3.5e-10) ** 2 * math.sqrt(8 * 1.380649e-23 * 200 / (math.pi * mu))
    assert pa.capture_rate(a, b, 200.0) == pytest.approx(hs, rel=1e-9)


def test_spin_factors():
    s = sp.load_species()
    assert pa.spin_factor(s["AlO"], s["AlO"]) == 0.25
    assert pa.spin_factor(s["AlO"], s["AlOH"]) == 1.0
    assert pa.spin_factor(s["AlOH"], s["AlOH"]) == 1.0


def test_every_pair_is_its_own_pathway():
    p = pa.formation_pathways()
    assert len(p) == 28 and len({x.pid for x in p}) == 28 and "AlO+AlO" in {x.pid for x in p}


def test_smoke_uptake_is_positive_in_the_mesosphere():
    sm, s = pa.SmokeBackground.load(), sp.load_species()
    assert sm.uptake_coefficient(s["AlO"], 85.0, 200.0, 2000.0) > 0


# ---- additions beyond the brief -------------------------------------------------------------------------------

def test_bin_count_and_last_pivot_at_ratio_two():
    b = pa.Bins()
    assert b.K == len(b.n_atoms) == 31
    assert b.n_atoms[-1] == 2.0 ** 31
    assert all(b.n_atoms[k] == 2.0 ** (k + 1) for k in range(b.K))


def test_other_ratio_still_brackets_the_largest_size():
    b = pa.Bins(bin_ratio=1.5, d_max_nm=100.0)
    assert pa.diameter(b.n_atoms[-1], 0.0) >= 100e-9 > pa.diameter(b.n_atoms[-2], 0.0)
    assert b.n_atoms[1] == pytest.approx(3.0)


def test_bin_ratio_must_exceed_one():
    with pytest.raises(ValueError):
        pa.Bins(bin_ratio=1.0)


def test_split_at_a_pivot_puts_everything_in_that_pivot():
    b = pa.Bins()
    assert b.split(8.0) == (2, 1.0)
    assert b.split(2.0) == (0, 1.0)


def test_above_the_last_pivot_the_entity_goes_to_the_last_bin_with_number_scaled():
    b = pa.Bins()
    n = 3.0 * b.n_atoms[-1]
    assert b.split(n) == (b.K - 2, 0.0)                     # all of the number in the upper bin of the pair
    assert b.number_scale(n) == pytest.approx(3.0)
    assert b.place(n) == [(b.K - 1, pytest.approx(3.0))]    # 3 particles of the last size carry the atoms
    assert b.atom_fraction(n, *b.split(n)) == 0.0           # no atoms in the lower bin
    # exactly at the last pivot: one particle of the last bin, no scaling
    assert b.place(b.n_atoms[-1]) == [(b.K - 1, pytest.approx(1.0))]


def test_place_conserves_number_atoms_inside_the_range():
    b = pa.Bins(); n = 5.0
    (k0, f0), (k1, f1) = b.place(n)
    assert (k0, k1) == (1, 2) and f0 + f1 == pytest.approx(1.0)
    assert f0 * b.n_atoms[k0] + f1 * b.n_atoms[k1] == pytest.approx(n)


def test_below_the_smallest_bin_is_an_error():
    with pytest.raises(ValueError, match="below"):
        pa.Bins().split(1.0)


def test_diameter_matches_the_sphere_of_the_composition_volume():
    d = pa.diameter(10, 3.0)
    assert math.pi / 6 * d ** 3 == pytest.approx(10 * 5.3524e-29, rel=1e-12)
    # more hydrogen than Al(OH)3 cannot make a particle bulkier than the hydroxide
    assert pa.volume_per_al(9.0) == pa.volume_per_al(3.0)


def test_free_molecule_kernel_with_infinite_partner_mass_uses_the_molecule_mass():
    k_inf = pa.fm_kernel(1e-9, math.inf, 2e-10, 43 * 1.66053906660e-27, 200.0)
    m = 43 * 1.66053906660e-27
    assert k_inf == pytest.approx(math.pi * (1.2e-9) ** 2 * math.sqrt(8 * 1.380649e-23 * 200 / (math.pi * m)), rel=1e-12)
    # equal-mass limit: reduced mass is half
    k_eq = pa.fm_kernel(1e-9, m, 1e-9, m, 200.0)
    assert k_eq == pytest.approx(math.pi * (2e-9) ** 2 * math.sqrt(8 * 1.380649e-23 * 200 / (math.pi * m / 2)), rel=1e-12)


def test_capture_rate_is_symmetric_and_never_below_hard_sphere():
    s = sp.load_species()
    for T in (150.0, 200.0, 298.0):
        ab, ba = pa.capture_rate(s["AlO"], s["AlOH"], T), pa.capture_rate(s["AlOH"], s["AlO"], T)
        assert ab == pytest.approx(ba, rel=1e-14)
        mu = s["AlO"].mass_amu * s["AlOH"].mass_amu / (s["AlO"].mass_amu + s["AlOH"].mass_amu) * 1.66053906660e-27
        assert ab >= math.pi * (3.5e-10) ** 2 * math.sqrt(8 * 1.380649e-23 * T / (math.pi * mu))


def test_zero_dipole_counts_as_unknown_for_the_dipole_term():
    s = sp.load_species()
    a = dataclasses.replace(s["AlO"], dipole_D=0.0)
    b = s["AlO"]
    m = b.mass_amu * 1.66053906660e-27
    hs = pa.fm_kernel(1.75e-10, m, 1.75e-10, m, 200.0)           # sigma_AB = 3.5 A, so r1 + r2 = 3.5 A
    assert pa.capture_rate(a, b, 200.0) == pytest.approx(hs, rel=1e-12)


def test_pathway_order_follows_species_order():
    p = pa.formation_pathways()
    assert p[0].pid == "Al+Al" and p[1].pid == "Al+AlO" and p[-1].pid == "AlOH2+AlOH2"
    assert all(sp.AL_SPECIES.index(x.a) <= sp.AL_SPECIES.index(x.b) for x in p)
    assert all(x.pid == f"{x.a}+{x.b}" for x in p)


def test_smoke_uptake_equals_the_direct_sum_at_a_table_altitude():
    sm, s = pa.SmokeBackground.load(), sp.load_species()
    mol = s["AlO"]
    r_mol = mol.sigma_A * 1e-10 / 2
    m_mol = mol.mass_amu * 1.66053906660e-27
    rho = 2000.0
    rows = [(r, n) for z, r, n in sm.rows() if z == 75.0]

    def mu(r):         # reduced mass of the molecule and one smoke sphere of radius r nm
        m_smoke = 4.0 / 3.0 * math.pi * (r * 1e-9) ** 3 * rho
        return m_smoke * m_mol / (m_smoke + m_mol)
    direct = sum(n * (math.log(200) / 34) * math.pi * (r * 1e-9 + r_mol) ** 2
                 * math.sqrt(8 * 1.380649e-23 * 200.0 / (math.pi * mu(r))) for r, n in rows)
    assert len(rows) == 16
    assert sm.uptake_coefficient(mol, 75.0, 200.0, rho) == pytest.approx(direct, rel=1e-12)


def test_smoke_of_finite_mass_is_taken_up_faster_than_infinitely_heavy_smoke():
    """One 0.2 nm bin of 2000 kg/m3 smoke against AlOH: the reduced-mass kernel against a hand computation."""
    mol = sp.load_species()["AlOH"]
    r, T, rho, n = 0.2e-9, 200.0, 2000.0, 1.0e12
    sm = pa.SmokeBackground([50.0, 100.0], [r], [[n / pa.SMOKE_DLNR], [n / pa.SMOKE_DLNR]])
    m_mol = mol.mass_amu * 1.66053906660e-27
    m_smoke = 4.0 / 3.0 * math.pi * r ** 3 * rho
    r_mol = mol.sigma_A * 1e-10 / 2
    hand = n * math.pi * (r + r_mol) ** 2 * math.sqrt(8 * 1.380649e-23 * T * (1 / m_mol + 1 / m_smoke) / math.pi)
    finite = sm.uptake_coefficient(mol, 70.0, T, rho)
    infinite = n * pa.fm_kernel(r, math.inf, r_mol, m_mol, T)
    assert finite == pytest.approx(hand, rel=1e-12)
    assert finite > infinite
    assert 1.1 < finite / infinite < 1.6       # a 0.2 nm particle weighs about as much as the molecule
    # a heavier smoke approaches the infinite-mass value and a bad density is refused
    assert sm.uptake_coefficient(mol, 70.0, T, 1e9) == pytest.approx(infinite, rel=1e-4)
    with pytest.raises(ValueError, match="density"):
        sm.uptake_coefficient(mol, 70.0, T, 0.0)


def test_smoke_uptake_between_table_altitudes_is_between_the_neighbours():
    sm, s = pa.SmokeBackground.load(), sp.load_species()
    lo, hi = (sm.uptake_coefficient(s["AlO"], z, 200.0, 2000.0) for z in (75.0, 80.0))
    mid = sm.uptake_coefficient(s["AlO"], 77.5, 200.0, 2000.0)
    assert min(lo, hi) < mid < max(lo, hi)
    # log-linear bin by bin: at the midpoint each bin holds the geometric mean of its neighbours' numbers
    assert sm.number_per_bin(77.5) == pytest.approx(np.sqrt(sm.number_per_bin(75.0) * sm.number_per_bin(80.0)), rel=1e-12)


@pytest.mark.parametrize("z", [44.9, 100.1, float("nan")])
def test_smoke_uptake_outside_the_table_raises_naming_altitude_and_file(z):
    sm, s = pa.SmokeBackground.load(), sp.load_species()
    with pytest.raises(ValueError, match=r"megner2006_smoke\.csv") as e:
        sm.uptake_coefficient(s["AlO"], z, 200.0, 2000.0)
    assert str(z) in str(e.value)
