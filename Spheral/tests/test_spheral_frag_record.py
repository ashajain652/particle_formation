"""spheral_frag.record: fragment record, debris log and mass accounts (M1 plan, Task 8; spec §11; check "Record and
accounts"). Every expected number is analytic and stated where it is used; the material is a stub of the
`RecordMaterial` protocol until Task 3's `MaterialTable` exists."""
import csv
import json
import math

import numpy as np
import pytest

from spheral_frag import contract, naming, record as R

RHO_S, RHO_L = 2813.0, 2400.0          # kg/m^3, solid and liquid AA7075_scheil (plan fact 8)
MU_L, SIGMA = 1.3e-3, 0.80             # Pa s, N/m (plan fact 8)
H_FULL_LIQUID = 1.1e6                  # J/kg from 300 K: a round stub value, the record only reads it back


class StubMaterial:
    """RecordMaterial stub: 1.3 mPa s above 908 K, inf below 50 % liquid (895.1 K), 1 Pa s between; records the
    arguments of every viscosity call."""

    def __init__(self, h_fl=H_FULL_LIQUID):
        self.h_fl = h_fl
        self.calls = []

    def viscosity(self, T, shear_rate=0.0):
        self.calls.append((float(T), shear_rate))
        T = np.asarray(T, float)
        return np.where(T > 908.0, MU_L, np.where(T > 895.1, 1.0, np.inf))

    def h_full_liquid(self):
        return self.h_fl


PROV = {"run": "prep_fe_k00000-00002__synthetic_k00000-00002_3d_dx0.50mm_seed1", "window": "synthetic",
        "dx_mm": 0.5, "seed": 1, "brackets": {}}


def particles(x, mass, rho=RHO_S, T=None, f_l=None, h=None, v=None, damage=None, h_smooth=1e-3, group=None,
              ids=None, t=1.0, frame=2, check=0):
    n = len(x)
    full = lambda val, default: np.full(n, default if val is None else val, float) if np.ndim(val) == 0 \
        else np.asarray(val, float)
    return R.Particles(t=t, id=np.arange(n) if ids is None else ids, mass=np.broadcast_to(mass, (n,)).copy(),
                       x=np.asarray(x, float), v=np.zeros((n, 3)) if v is None else np.asarray(v, float),
                       T=full(T, 300.0), f_l=full(f_l, 0.0), h=full(h, 0.0), rho=full(rho, RHO_S),
                       damage=full(damage, 0.0), h_smooth=full(h_smooth, 1e-3),
                       group=np.zeros(n, int) if group is None else group, frame=frame, check=check)


def lattice(lo, hi, dx):
    axes = [np.arange(l, h + dx / 2, dx) for l, h in zip(lo, hi)]
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)


def flight(rho_inf=0.01, v_kms=5.0):
    return {"altitude_km": 60.0, "velocity_kms": v_kms, "flight_path_deg": -5.0, "density_kgm3": rho_inf,
            "T_air_K": 250.0, "deceleration_ms2": 40.0}


# ----------------------------------------------------------------------------------------------- size and shape

@pytest.mark.parametrize("rotate", [False, True])
def test_lattice_ellipsoid_principal_lengths_and_equivalent_diameter(rotate):
    # plan Task 8: semi-axes 10, 6, 4 mm at dx = 0.5 mm -> principal lengths 20, 12, 8 mm within dx, d_eq within 1 %
    a, b, c, dx = 10e-3, 6e-3, 4e-3, 0.5e-3
    pts = lattice((-12e-3,) * 3, (12e-3,) * 3, dx)
    pts = pts[(pts[:, 0] / a) ** 2 + (pts[:, 1] / b) ** 2 + (pts[:, 2] / c) ** 2 <= 1.0]
    if rotate:      # lengths are invariant under rotation and translation
        q, _ = np.linalg.qr(np.random.default_rng(3).normal(size=(3, 3)))
        pts = pts @ q.T + np.array([0.02, -0.01, 0.03])
    p = particles(pts, RHO_S * dx ** 3, f_l=0.0)
    row = R.fragment_row(p, np.ones(p.n, bool), np.zeros(p.n, bool), StubMaterial(), flight(), "separation",
                         "tearing", PROV)
    for got, want in zip((row["L1_m"], row["L2_m"], row["L3_m"]), (2 * a, 2 * b, 2 * c)):
        assert abs(got - want) < dx
    d_exact = 2.0 * (a * b * c) ** (1.0 / 3.0)
    assert abs(row["d_eq_m"] / d_exact - 1.0) < 0.01          # measured 0.10 %
    assert row["volume_m3"] == pytest.approx(math.fsum(p.mass / p.rho), rel=1e-15)
    # area is Thomsen's of the row's own semi-axes, by definition ("approximate")
    assert row["area_m2"] == pytest.approx(R.ellipsoid_area(row["L1_m"] / 2, row["L2_m"] / 2, row["L3_m"] / 2),
                                           rel=1e-12)
    assert row["ring"] is False and math.isnan(row["ring_radius_m"])


def test_thomsen_area_against_exact_sphere_and_bound():
    assert R.ellipsoid_area(5e-3, 5e-3, 5e-3) == pytest.approx(4 * math.pi * 25e-6, rel=1e-14)   # exact for a sphere
    # exact area of the (10, 6, 4) mm ellipsoid by quadrature, 5.391007e-4 m^2 (scipy dblquad, 2026-10-07)
    assert abs(R.ellipsoid_area(10e-3, 6e-3, 4e-3) / 5.391007068739417e-4 - 1.0) < 0.01061


def test_rz_torus_ring_flag_radius_lengths_volume_area():
    # a ring of circular cross-section a = 3 mm at R = 20 mm, meridional lattice dx = 0.25 mm, ring masses
    dx, Rr, a = 0.25e-3, 20e-3, 3e-3
    pts = lattice((dx / 2, -5e-3, 0.0), (25e-3, 5e-3, 0.0), dx)
    pts = pts[(pts[:, 0] - Rr) ** 2 + pts[:, 1] ** 2 <= a * a]
    p = particles(pts, RHO_L * 2 * math.pi * pts[:, 0] * dx * dx, rho=RHO_L, h_smooth=2 * dx, f_l=1.0, T=1000.0)
    row = R.fragment_row(p, np.ones(p.n, bool), np.zeros(p.n, bool), StubMaterial(), flight(), "separation",
                         "ring_release", PROV, form="rz")
    assert row["ring"] is True
    # mass-weighted r of a uniform torus: (R^2 + a^2/4) / R = 20.1125 mm (plan: "the mass-weighted r")
    assert abs(row["ring_radius_m"] - (Rr + a * a / (4 * Rr))) < dx / 10
    assert row["L1_m"] == pytest.approx(2 * math.pi * row["ring_radius_m"], rel=1e-15)
    assert abs(row["L2_m"] - 2 * a) < dx and abs(row["L3_m"] - 2 * a) < dx
    assert row["mass_kg"] == pytest.approx(math.fsum(p.mass), rel=0)          # already the whole ring's
    assert abs(row["volume_m3"] / (2 * math.pi ** 2 * Rr * a * a) - 1) < 0.01   # measured 0.09 %
    assert abs(row["area_m2"] / (4 * math.pi ** 2 * Rr * a) - 1) < 0.01         # Pappus; measured 0.05 %
    # placed in the meridional half-plane phi = 0: at the equator, 90 degrees from the nose
    assert row["theta_deg"] == pytest.approx(90.0, abs=1e-9) and math.isnan(row["phi_deg"])
    assert row["y_m"] == row["ring_radius_m"] and row["z_m"] == 0.0


def test_rz_cap_on_the_axis_is_not_a_ring_and_gets_the_body_of_revolution_lengths():
    # a spheroid cap on the axis: equatorial semi-axis 10 mm, axial 6 mm, centred 40 mm ahead of the centre
    dx, A, C, z0 = 0.25e-3, 10e-3, 6e-3, 40e-3
    pts = lattice((dx / 2, z0 - C, 0.0), (A, z0 + C, 0.0), dx)
    pts = pts[(pts[:, 0] / A) ** 2 + ((pts[:, 1] - z0) / C) ** 2 <= 1]
    p = particles(pts, RHO_S * 2 * math.pi * pts[:, 0] * dx * dx, h_smooth=2 * dx)
    row = R.fragment_row(p, np.ones(p.n, bool), np.zeros(p.n, bool), StubMaterial(), flight(), "detachment",
                         "neck_failure", PROV, form="rz")
    assert row["ring"] is False and math.isnan(row["ring_radius_m"])
    for got, want in zip((row["L1_m"], row["L2_m"], row["L3_m"]), (2 * A, 2 * A, 2 * C)):
        assert abs(got - want) < dx
    assert abs(row["volume_m3"] / (4 / 3 * math.pi * A * A * C) - 1) < 0.01      # measured 0.25 %
    assert (row["x_m"], row["y_m"], row["z_m"]) == (pytest.approx(z0, abs=1e-12), 0.0, 0.0)
    assert row["theta_deg"] == pytest.approx(0.0, abs=1e-9)


# ----------------------------------------------------------------------------------------------- phase and breakup

def drop(d=10e-3, dx=0.5e-3, f_l=1.0, T=1000.0):
    """A 10 mm drop whose particle masses make sum m / rho = pi d^3 / 6 (so d_eq = d to round-off)."""
    pts = lattice((-d / 2,) * 3, (d / 2,) * 3, dx)
    pts = pts[np.linalg.norm(pts, axis=1) <= d / 2]
    return particles(pts, RHO_L * math.pi * d ** 3 / 6 / len(pts), rho=RHO_L, f_l=f_l, T=T, h=1.2e6)


def test_weber_ohnesorge_threshold_and_flag_of_a_fluid_drop():
    # plan Task 8: rho_inf 0.01, V 5,000 m/s, d 10 mm, mu 1.3 mPa s, rho 2,400, sigma 0.80
    # -> We = 0.01 * 5000^2 * 0.01 / 0.8 = 3,125; Oh = 1.3e-3 / sqrt(2400 * 0.8 * 0.01) = 2.967e-4; threshold 12.000
    assert R.weber(0.01, 5000.0, 0.01) == pytest.approx(3125.0, rel=1e-12)
    assert R.ohnesorge(1.3e-3, 2400.0, 0.01) == pytest.approx(2.967e-4, rel=1e-3)
    p, mat = drop(), StubMaterial()
    row = R.fragment_row(p, np.ones(p.n, bool), np.zeros(p.n, bool), mat, flight(), "slurry", "slurry_breakup", PROV)
    assert row["phase"] == "fluid" and row["fluid_mass_fraction"] == 1.0
    assert row["d_eq_m"] == pytest.approx(0.01, rel=1e-12)
    assert row["weber"] == pytest.approx(3125.0, rel=1e-3)
    assert row["ohnesorge"] == pytest.approx(2.967e-4, rel=1e-3)
    assert row["breakup_threshold"] == pytest.approx(12.000, abs=5e-4)
    assert row["breakup"] == "yes" and row["mu_Pas"] == MU_L and row["sigma_Nm"] == SIGMA
    # decision 7: the viscosity at the mass-weighted temperature, low-shear limit
    assert mat.calls == [(pytest.approx(row["T_mean_K"], rel=1e-15), 0.0)]
    slow = R.fragment_row(p, np.ones(p.n, bool), np.zeros(p.n, bool), mat, flight(v_kms=0.3), "slurry",
                          "slurry_breakup", PROV)
    assert slow["weber"] == pytest.approx(11.25, rel=1e-3) and slow["breakup"] == "no"    # 11.25 < 12.000
    assert not R.breakup_flag(11.25, 2.967e-4) and R.breakup_flag(3125.0, 2.967e-4)
    assert R.breakup_threshold(0.0) == 12.0


@pytest.mark.parametrize("f_l, phase, frac", [(np.r_[np.ones(10), np.full(10, 0.4)], "mixed", 0.5),
                                              (np.full(20, 0.3), "solid_or_mush", 0.0)])
def test_mixed_or_solid_fragment_has_no_breakup_numbers(f_l, phase, frac):
    p = particles(lattice((0, 0, 0), (4e-3, 3e-3, 0), 1e-3), 1e-6, f_l=f_l, T=900.0)
    row = R.fragment_row(p, np.ones(p.n, bool), np.zeros(p.n, bool), StubMaterial(), flight(), "separation",
                         "tearing", PROV)
    assert row["phase"] == phase and row["fluid_mass_fraction"] == frac
    assert row["breakup"] == R.NOT_APPLICABLE
    assert all(math.isnan(row[k]) for k in ("weber", "ohnesorge", "mu_Pas", "breakup_threshold"))


def test_phase_state_threshold_is_strictly_more_than_half_liquid():
    assert R.phase_state([0.5, 0.9], [1.0, 3.0]) == ("mixed", 0.75)
    assert R.phase_state([0.51, 1.0], [1.0, 1.0]) == ("fluid", 1.0)
    assert R.phase_state([0.5, 0.0], [1.0, 1.0]) == ("solid_or_mush", 0.0)


def test_thermal_columns_and_heat_to_melt():
    h = np.array([0.5e6, 1.0e6, 1.2e6, 1.1e6])
    m = np.array([1e-3, 2e-3, 3e-3, 4e-3])
    p = particles(np.eye(4, 3) * 1e-2, m, h=h, T=[700.0, 850.0, 950.0, 908.0], f_l=[0, 0.2, 1, 1])
    row = R.fragment_row(p, np.ones(4, bool), np.zeros(4, bool), StubMaterial(), flight(), "separation", "tearing",
                         PROV)
    assert row["heat_to_melt_J"] == pytest.approx(1e-3 * 0.6e6 + 2e-3 * 0.1e6, rel=1e-15)   # sum m max(h_fl - h, 0)
    assert row["h_mean_Jkg"] == pytest.approx(np.dot(m, h) / m.sum(), rel=1e-15)
    assert row["T_mean_K"] == pytest.approx(np.dot(m, p.T) / m.sum(), rel=1e-15)
    assert (row["T_min_K"], row["T_max_K"]) == (700.0, 950.0)
    assert row["f_l_mean"] == pytest.approx((2e-3 * 0.2 + 7e-3) / 1e-2, rel=1e-15)


# ----------------------------------------------------------------------------------------------- groups

def block(n_side=10, dx=1e-3, h=1.3e-3, offset=(0, 0, 0)):
    return lattice(offset, np.add(offset, (n_side - 1) * dx), dx), h


def test_dust_within_link_distance_attaches_and_beyond_goes_to_the_debris_log():
    pts, h = block()
    face = np.array([9e-3, 4e-3, 4e-3])              # a lattice node on the +x face; dust placed outward from it
    near = face + np.array([1.5 * h * 0.99, 0, 0])   # within 1.5 h of the intact node (both h equal)
    far = face + np.array([1.5 * h * 1.01, 0, 0])    # just beyond
    x = np.vstack([pts, near, far])
    damage = np.r_[np.zeros(len(pts)), 1.0, 0.995]   # both dust (>= 0.99)
    group = np.r_[np.zeros(len(pts), int), 3, 4]     # the finder's labels for dust are irrelevant
    p = particles(x, 1e-6, h_smooth=h, damage=damage, group=group)
    g = R.classify_groups(p)
    assert g["main_label"] == 0
    assert g["main"][-2] and not g["main"][-1]
    assert np.flatnonzero(g["debris_dust"]).tolist() == [len(pts) + 1]
    assert g["fragments"] == {} and g["debris_groups"] == {}
    # a fragment's attached dust counts with it and is recorded as its dust mass
    frag_pts = pts + np.array([0.1, 0, 0])
    dust = frag_pts[0] - np.array([1.0 * h, 0, 0])
    x2 = np.vstack([pts, frag_pts, dust])
    p2 = particles(x2, 1e-6, h_smooth=h, damage=np.r_[np.zeros(2 * len(pts)), 1.0],
                   group=np.r_[np.zeros(len(pts), int), np.ones(len(pts), int), 0])
    p2.mass[: len(pts)] *= 2                       # the main body has the most mass
    g2 = R.classify_groups(p2)
    assert list(g2["fragments"]) == [1] and g2["fragment_dust"][1][-1] and g2["fragments"][1][-1]
    row = R.fragment_row(p2, g2["fragments"][1], g2["fragment_dust"][1], StubMaterial(), flight(), "separation",
                         "tearing", PROV)
    assert row["dust_mass_kg"] == 1e-6 and row["n_particles"] == len(pts) + 1


@pytest.mark.parametrize("n_small, where", [(29, "debris_groups"), (30, "fragments")])
def test_resolution_floor_of_thirty_particles(n_small, where):
    pts, h = block()
    small = lattice((0.2, 0, 0), (0.2 + 9e-3, 2e-3, 0), 1e-3)[:n_small]
    assert len(small) == n_small
    p = particles(np.vstack([pts, small]), 1e-6, h_smooth=h,
                  group=np.r_[np.full(len(pts), 5), np.full(n_small, 2)])
    g = R.classify_groups(p)
    assert g["main_label"] == 5
    assert list(g[where]) == [2] and np.count_nonzero(g[where][2]) == n_small
    other = "fragments" if where == "debris_groups" else "debris_groups"
    assert g[other] == {}
    assert list(R.classify_groups(p, min_particles=n_small + 1)["debris_groups"]) == [2]


def test_main_body_is_the_group_with_most_mass():
    pts, h = block(4)
    p = particles(np.vstack([pts, pts + 0.1]), np.r_[np.full(len(pts), 1.0), np.full(len(pts), 2.0)],
                  group=np.r_[np.zeros(len(pts), int), np.full(len(pts), 7)])
    assert R.classify_groups(p, min_particles=1)["main_label"] == 7


def test_match_groups_follows_permuted_labels_and_partial_removal():
    rng = np.random.default_rng(0)
    ids = rng.permutation(1000)
    prev = {0: ids[:700], 1: ids[700:900], 2: ids[900:950], 3: ids[950:]}
    curr = {11: rng.permutation(ids[:700])[:650],          # the main body, 50 particles removed
            4: ids[700:900][::2],                           # the fragment, half its particles removed
            8: ids[900:950]}                                # group 3 is gone entirely
    assert R.match_groups(prev, curr) == {0: 11, 1: 4, 2: 8}
    # a split goes to the piece holding the larger share; a merge maps two groups to one
    split = {20: ids[700:840], 21: ids[840:900], 30: ids[:700]}
    assert R.match_groups(prev, split)[1] == 20
    merged = {40: np.r_[ids[700:900], ids[900:950]], 30: ids[:700]}
    assert R.match_groups(prev, merged) == {0: 30, 1: 40, 2: 40}
    with pytest.raises(ValueError):
        R.match_groups(prev, {1: ids[:10], 2: ids[5:20]})
    # group_ids builds the inputs from masks
    p = particles(np.zeros((4, 3)), 1.0, ids=np.array([10, 11, 12, 13]))
    assert {k: v.tolist() for k, v in R.group_ids(p, {0: np.array([1, 0, 1, 0], bool)}).items()} == {0: [10, 12]}


def test_is_clear_distance_and_direction():
    pts, h = block(5)
    for gap, vx, want in ((4.0 * h, +1.0, True), (4.0 * h, -1.0, False), (2.5 * h, +1.0, False)):
        piece = lattice((0, 0, 0), (2e-3, 2e-3, 2e-3), 1e-3) + np.array([4e-3 + gap, 0, 0])
        v = np.vstack([np.zeros((len(pts), 3)), np.tile([vx, 0, 0], (len(piece), 1))])
        p = particles(np.vstack([pts, piece]), 1e-6, h_smooth=h, v=v)
        sel = np.r_[np.zeros(len(pts), bool), np.ones(len(piece), bool)]
        assert R.is_clear(p, sel, ~sel) is want, (gap, vx)


# ----------------------------------------------------------------------------------------------- flight state, columns

def test_flight_table_interpolates_linearly_and_refuses_outside(tmp_path):
    cols = {"time_s": np.array([0.0, 0.5, 1.0]), "altitude_km": np.array([77.5, 77.4, 77.2]),
            "velocity_kms": np.array([7.5, 7.5, 7.49]), "flight_path_deg": np.array([-1.0, -1.0, -1.1]),
            "density_kgm3": np.array([1e-5, 2e-5, 4e-5]), "T_air_K": np.array([200.0, 201.0, 203.0]),
            "deceleration_ms2": np.array([1.0, 2.0, 4.0]), "drag_N": np.zeros(3)}
    np.savez(tmp_path / "flight.npz", **cols)
    ft = R.FlightTable.load(tmp_path / "flight.npz")
    s = ft.at(0.75)
    assert s["density_kgm3"] == pytest.approx(3e-5, rel=1e-15) and s["T_air_K"] == pytest.approx(202.0)
    assert s["deceleration_ms2"] == pytest.approx(3.0) and set(s) == set(R.FLIGHT_STATE)
    with pytest.raises(ValueError):
        ft.at(1.01)
    # the history's part of the flight state is spelled as the contract's history keys
    history = {f.key for f in contract.fields("history")}
    assert set(R.FLIGHT_STATE) - history == {"T_air_K", "deceleration_ms2"} and R.FLIGHT_TIME in history


def test_columns_are_unique_grouped_and_rows_follow_them(tmp_path):
    assert len(R.FRAGMENT_COLUMNS) == len(set(R.FRAGMENT_COLUMNS))
    assert list(R.FRAGMENT_GROUPS) == ["When", "Where and how fast", "Size and shape", "Thermal state", "Origin",
                                       "Phase and breakup", "Provenance"]
    assert len(R.DEBRIS_COLUMNS) == len(set(R.DEBRIS_COLUMNS))
    p = drop(dx=1e-3)
    row = R.fragment_row(p, np.ones(p.n, bool), np.zeros(p.n, bool), StubMaterial(), flight(), "synthetic",
                         "ring_release", dict(PROV, brackets={"film_limit_mm": 3.0}), number=4)
    assert list(row) == R.FRAGMENT_COLUMNS and row["fragment"] == 4 and row["count"] == 1
    assert row["bracket_film_limit_mm"] == 3.0 and row["bracket_min_particles"] == 30
    path = R.write_rows(str(tmp_path / "fragments.csv"), [row], R.FRAGMENT_COLUMNS)
    with open(path, newline="") as fh:
        back = list(csv.DictReader(fh))[0]
    assert float(back["weber"]) == row["weber"] and back["resolved"] == "true" and back["ring"] == "false"
    with pytest.raises(ValueError):
        R.fragment_row(p, np.ones(p.n, bool), np.zeros(p.n, bool), StubMaterial(), flight(), "synthetic",
                       "explosion", PROV)
    with pytest.raises(TypeError):
        R.fragment_row(p, np.ones(p.n, bool), np.zeros(p.n, bool), object(), flight(), "synthetic",
                       "tearing", PROV)


def test_debris_row_size_upper_bound_is_extent_plus_one_spacing():
    x = np.array([[0.04, 0, 0], [0.043, 0.004, 0], [0.041, 0.001, 0]])     # extent 5 mm (3-4-5)
    p = particles(x, [1e-6, 2e-6, 3e-6], h=[0.0, 1e6, 2e6], T=[800.0, 900.0, 1000.0])
    row = R.debris_row(p, np.ones(3, bool), flight(), dx=1e-3, table=StubMaterial(), number=2, kind="dust")
    assert list(row) == R.DEBRIS_COLUMNS
    assert row["size_upper_m"] == pytest.approx(6e-3, rel=1e-12)
    assert row["mass_kg"] == 6e-6 and row["n_particles"] == 3 and row["kind"] == "dust"
    assert row["heat_to_melt_J"] == pytest.approx(1e-6 * 1.1e6 + 2e-6 * 0.1e6, rel=1e-14)
    assert math.isnan(R.debris_row(p, np.ones(3, bool), None, dx=1e-3)["heat_to_melt_J"])


# ----------------------------------------------------------------------------------------------- files

def test_check_removed_and_meta_files_round_trip_bitwise(tmp_path):
    rng = np.random.default_rng(1)
    n = 50
    p = R.Particles(t=12.5, id=rng.permutation(10 ** 6)[:n], mass=rng.random(n), x=rng.normal(size=(n, 3)),
                    v=rng.normal(size=(n, 3)), T=rng.random(n) * 900, f_l=rng.random(n), h=rng.random(n) * 1e6,
                    rho=np.full(n, RHO_S), damage=rng.random(n), h_smooth=np.full(n, 2e-3),
                    group=rng.integers(0, 4, n), frame=25, check=7)
    path = R.write_check(str(tmp_path), p)
    assert path.endswith("checks/check_000007.npz")
    q = R.read_check(path)
    assert (q.t, q.frame, q.check) == (12.5, 25, 7)
    for name in R.CHECK_ARRAYS:
        if name not in ("t", "frame", "check"):
            a, b = getattr(p, name), getattr(q, name)
            assert a.dtype == b.dtype and np.array_equal(a, b), name
    assert R.list_checks(str(tmp_path)) == [(7, path)]
    rem = R.Removal(id=[3, 4, 3], t=[1.0, 1.5, 2.0], check=[1, 2, 3], mass=[1e-4, 2e-4, 5e-5],
                    h=[1e6, 1.2e6, 1.3e6], reason=[R.REASON_SINK, R.REASON_CENTRE_CROSSING, R.REASON_FLOOR])
    back = R.read_removed(R.write_removed(str(tmp_path / "removed.npz"), rem))
    for name in R.REMOVED_ARRAYS:
        assert getattr(back, name).dtype == np.dtype(R.REMOVED_ARRAYS[name][0])
        assert np.array_equal(getattr(back, name), getattr(rem, name)), name
    assert len(rem.until(2)) == 2
    with np.load(path, allow_pickle=False) as z:          # plain numeric arrays only (no object dtype)
        assert all(z[k].dtype.kind in "if" for k in z.files)
    run = naming.run_name("prep_fe_k00000-00002", "synthetic", (0, 2), "3d", 0.5, {"link_h": 2.0}, 1)
    meta = R.make_meta(run, "prep_fe_k00000-00002", "synthetic", "3d", 0.5, {"link_h": 2.0}, 1, 1.47, 3.2e5, n,
                       spheral={"commit": "116c71f"})
    assert R.read_meta(R.write_meta(str(tmp_path / "meta.json"), meta)) == meta
    with pytest.raises(ValueError):          # the meta must agree with its own run name
        R.make_meta(run, "prep_fe_k00000-00002", "synthetic", "3d", 0.5, {}, 1, 1.47, 3.2e5, n)
    with pytest.raises(ValueError):
        R.validate_meta({k: v for k, v in meta.items() if k != "initial_mass_kg"})
    json.dumps(meta, allow_nan=False)


# ----------------------------------------------------------------------------------------------- accounts

def test_accounts_close_on_a_synthetic_history_and_every_removal_carries_its_enthalpy(tmp_path):
    """Removals by centre crossing, sink steps with floor removals, a released fragment and isolated debris over ten
    checks; at every check |residual| <= 1e-14 m0 (plan Task 8), and each account's enthalpy is exactly the fsum of
    m h of what it holds (the film's split by reason, the sink's at the sprayed liquid's enthalpy)."""
    rng = np.random.default_rng(42)
    n = 1500
    x = rng.uniform(-0.04, 0.04, size=(n, 3))
    mass = rng.lognormal(np.log(1e-3), 0.3, n)
    p = particles(x, mass, h=rng.uniform(2e5, 1e6, n), T=rng.uniform(300, 1000, n), h_smooth=3e-3,
                  ids=rng.permutation(10 ** 6)[:n], check=0, t=0.0)
    m0 = math.fsum(p.mass)
    m_start = dict(zip(p.id.tolist(), p.mass.tolist()))
    H_SPRAY = 1.25e6                                 # sprayed liquid's enthalpy (§14.1), not the particle's own
    removals, booked = [], {code: [] for code in R.REASONS}
    alive = np.ones(n, bool)
    for check in range(1, 11):
        t = 0.5 * check
        p.t, p.check = t, check
        p.h += 1e4                                    # the replay heats the particles (enthalpy not conserved)
        # centre crossing: 15 whole particles
        cross = rng.choice(np.flatnonzero(alive), 15, replace=False)
        for i in cross:
            removals.append((p.id[i], t, check, p.mass[i], p.h[i], R.REASON_CENTRE_CROSSING))
            booked[R.REASON_CENTRE_CROSSING].append(p.mass[i] * p.h[i])
        alive[cross] = False
        # sink: 40 particles lose 5-40 %; booked as m_old - m_new (exact by Sterbenz), at H_SPRAY
        sink = rng.choice(np.flatnonzero(alive), 40, replace=False)
        for i in sink:
            m_old = p.mass[i]
            p.mass[i] = m_old * (1.0 - rng.uniform(0.05, 0.4))
            dm = m_old - p.mass[i]
            removals.append((p.id[i], t, check, dm, H_SPRAY, R.REASON_SINK))
            booked[R.REASON_SINK].append(dm * H_SPRAY)
            if p.mass[i] < 0.25 * m_start[int(p.id[i])]:          # floor: the rest goes to the film account
                removals.append((p.id[i], t, check, p.mass[i], p.h[i], R.REASON_FLOOR))
                booked[R.REASON_FLOOR].append(p.mass[i] * p.h[i])
                alive[i] = False
        if check == 4:        # a fragment of 60 particles separates and moves off (it stays in the particles)
            frag = np.flatnonzero(alive)[:60]
            p.group[frag] = 1
            p.x[frag] += np.array([0.5, 0, 0])
        if check == 7:        # isolated dust
            dust = np.flatnonzero(alive & (p.group == 0))[:12]
            p.damage[dust] = 1.0
            p.x[dust] += np.array([0, 2.0, 0]) + rng.normal(scale=0.1, size=(12, 3))
        rem = R.Removal(*map(np.array, zip(*removals)))
        cur = R.Particles(**{f: (getattr(p, f)[alive] if isinstance(getattr(p, f), np.ndarray) else getattr(p, f))
                             for f in ("t", "id", "mass", "x", "v", "T", "f_l", "h", "rho", "damage", "h_smooth",
                                       "group", "frame", "check")})
        # through the runner -> analyse files
        cur = R.read_check(R.write_check(str(tmp_path), cur))
        rem = R.read_removed(R.write_removed(str(tmp_path / "removed.npz"), rem))
        g = R.classify_groups(cur, min_particles=30)
        acc = R.accounts_at(m0, cur, g, rem)
        assert abs(acc["residual_kg"]) <= 1e-14 * m0, (check, acc["residual_kg"])
        assert acc["mass_total_kg"] == pytest.approx(m0, rel=1e-14)
        for code, reason in R.REASONS.items():
            assert acc["enthalpy_film_{}_J".format(reason)] == math.fsum(booked[code])        # exactly
        assert acc["enthalpy_film_J"] == math.fsum(rem.mass * rem.h)
        assert acc["mass_film_kg"] == math.fsum([r[3] for r in removals])
        main = g["main"]
        assert acc["enthalpy_main_J"] == math.fsum(cur.mass[main] * cur.h[main])
        assert acc["n_fragments"] == (1 if check >= 4 else 0)
        assert acc["n_debris_dust"] == (12 if check >= 7 else 0)
    # a removal whose mass was dropped breaks the balance by exactly that mass
    short = rem.subset(np.arange(len(rem)) != 5)
    acc = R.accounts_at(m0, cur, g, short)
    assert acc["residual_kg"] == pytest.approx(-rem.mass[5], rel=1e-9)
    assert abs(acc["residual_kg"]) > 1e-12 * m0


def test_accounts_at_refuses_groups_that_do_not_partition():
    p = particles(np.zeros((3, 3)), 1.0)
    g = {"main": np.array([1, 1, 0], bool), "fragments": {}, "debris_groups": {},
         "debris_dust": np.zeros(3, bool), "dust": np.zeros(3, bool)}
    with pytest.raises(ValueError):
        R.accounts_at(3.0, p, g, R.Removal.empty())
