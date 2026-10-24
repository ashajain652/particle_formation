"""Load tables by inclination (spheral_frag.loads; Spheral M1 plan, Task 7; spec §7.2-7.3, check "Load tables").

On the committed sphere run (tests/fixtures/spheral_frag/sphere_run/), whose evaluated frames carry modified-Newtonian
loads (p_w = p_stag cos^2 theta and tau = 0.02 p_stag sin 2 theta for theta < 90 deg, 0.01 p_stag to 150 deg, 0 =
"not evaluated" beyond) and whose history's drag is the frame's own. Frame 0 is written before any evaluation and
carries no loads. The patch geometry here comes from the fixture's faces (wound outward, README), through the
builder's own `patch_geometry`; prepare (Task 9) takes it from Task 4's `derived_patch_arrays`."""
import csv
import math
import os

import numpy as np
import pytest

import spheral_frag_synthetic as syn
from helpers import REPO_ROOT
from spheral_frag import loads

pytest.importorskip("pyvista")

RUN = os.path.join(REPO_ROOT, "tests", "fixtures", "spheral_frag", "sphere_run", syn.SPHERE_RUN_NAME)
CSV_RTOL = 5e-9           # a value through "{:.9g}" is within half a unit of its 9th digit (fixtures README)
ROUNDOFF = 1e-12          # a priori: sums of ~1,200 double products
FACETING = 1e-2           # plan Task 7: faceted against continuum, "expected below 1 %"; README measures -0.55 %
BINNING = 1e-3            # plan Task 7: D_table within 0.1 % of D_patch


@pytest.fixture(scope="module")
def frames():
    with open(RUN + ".csv") as fh:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]
    out = []
    for k in range(len(syn.SPHERE_TIMES)):
        frame, _ = syn.read_frame_pyvista(RUN, k)
        area, normal, centroid, theta = syn.patch_geometry(frame.points, frame.faces)
        out.append({"frame": frame, "row": rows[k], "area": area, "normal": normal, "theta": theta,
                    "p_w": frame.patch["p_w"], "tau": frame.patch["tau"], "p_stag": rows[k]["p_w_stag_step_Pa"]})
    return out


def table_of(f):
    return loads.build_table(f["theta"], f["area"], f["p_w"], f["tau"], f["p_stag"])


# ------------------------------------------------------------------------------------------------ build_table
@pytest.mark.parametrize("k", [1, 2])
def test_table_recovers_newtonian_loads_within_the_bin_averaging_error(frames, k):
    """An area-weighted bin mean of p_stag g(theta_patch) lies within the range of g over the bin; g is monotonic on
    every 1 deg bin (cos^2 on 0-90 deg, sin 2theta on either side of 45 deg), so the range is its two edge values."""
    f = frames[k]
    t = table_of(f)
    lo, hi = np.radians(t.edges_deg[:-1]), np.radians(t.edges_deg[1:])
    filled = np.isfinite(t.p)
    wind = filled & (t.edges_deg[1:] <= 90.0)
    for law, values, rel in ((lambda x: np.cos(x) ** 2, t.p, 1.0), (lambda x: np.sin(2 * x), t.tau, syn.TAU_FRACTION)):
        a, b = rel * f["p_stag"] * law(lo[wind]), rel * f["p_stag"] * law(hi[wind])
        slack = ROUNDOFF * f["p_stag"]
        assert np.all(values[wind] >= np.minimum(a, b) - slack)
        assert np.all(values[wind] <= np.maximum(a, b) + slack)
    base = filled & (t.edges_deg[:-1] >= 90.0)
    np.testing.assert_allclose(t.p[base], syn.BASE_FRACTION * f["p_stag"], rtol=ROUNDOFF)
    np.testing.assert_array_equal(t.tau[base], 0.0)
    # the last bin with loads is the one below 150 deg; nothing beyond (README "Loads table")
    deg = np.degrees(f["theta"])
    assert t.theta_last_deg == math.floor(deg[deg < syn.THETA_EVAL_DEG].max()) + 0.5 == 149.5
    assert not np.any(filled[t.theta_deg > t.theta_last_deg])
    valid = f["p_w"] > 0
    assert np.array_equal(valid, deg < syn.THETA_EVAL_DEG)
    assert math.isclose(t.area.sum(), f["area"][valid].sum(), rel_tol=ROUNDOFF)
    assert t.p_stag == f["p_stag"]


def test_frame_without_loads_gives_an_empty_table(frames):
    """Frame 0 is written before any evaluation (p_w = tau = 0 everywhere): every bin NaN, theta_last NaN, the lee
    extension leaves it empty, evaluation gives NaN, and the comparison reports it as a frame without loads."""
    f = frames[0]
    assert not np.any(f["p_w"] > 0)
    t = table_of(f)
    assert not t.has_loads and math.isnan(t.theta_last_deg) and math.isnan(t.windward_edge_deg)
    assert np.all(np.isnan(t.p)) and np.all(np.isnan(t.tau)) and np.all(t.area == 0.0)
    lee = loads.with_lee(t)
    assert not lee.has_loads and np.all(np.isnan(lee.p)) and not lee.lee.any()
    p, tau = loads.evaluate(lee, np.array([0.0, 90.0, 180.0]))
    assert np.all(np.isnan(p)) and np.all(np.isnan(tau))
    c = loads.drag_comparison(f["theta"], f["area"], f["p_w"], f["tau"], f["row"])
    assert c["loads_valid"] is False and c["n_valid"] == 0
    assert c["D_patch"] == 0.0 and c["rel_patch_hist"] == -1.0
    assert math.isnan(c["D_table"]) and math.isnan(c["D_table_lee"]) and math.isnan(c["lee_share"])
    # the history's row 0 still records the drag the analytic loads would give (README: 6.215014915 N)
    assert math.isclose(c["D_hist"], 6.215014915, rel_tol=2 * CSV_RTOL)


# ------------------------------------------------------------------------------------------------ the three drags
def independent_drag(frame, p, tau):
    """The force on the faceted body along -v_hat from the triangles' vertices alone (no theta, no unit normals):
    pressure: sum p (1/2)((b - a) x (c - a)) . v_hat; shear along the oncoming flow (-v_hat) projected on the
    facet, of magnitude tau: its component along -v_hat is tau |t_raw| / 1 with t_raw = -v - (-v . n) n."""
    x = frame.points[frame.faces]
    N = 0.5 * np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0])      # area-weighted outward normal
    v = np.array(syn.V_HAT)
    A = np.linalg.norm(N, axis=1)
    n = N / A[:, None]
    t_raw = -v - (n @ -v)[:, None] * n
    t_norm = np.linalg.norm(t_raw, axis=1)
    shear_along = np.where(t_norm > 0, (t_raw @ -v) / np.where(t_norm > 0, t_norm, 1.0), 0.0)
    return math.fsum(p * (N @ v)) + math.fsum(tau * A * shear_along)


@pytest.mark.parametrize("k, readme_N", [(1, 6.866822220), (2, 7.311101107)])
def test_patch_drag_closed_form_and_history(frames, k, readme_N):
    f = frames[k]
    c = loads.drag_comparison(f["theta"], f["area"], f["p_w"], f["tau"], f["row"])
    D_indep = independent_drag(f["frame"], f["p_w"], f["tau"])
    assert abs(c["D_patch"] - D_indep) <= ROUNDOFF * abs(D_indep)
    assert c["D_patch"] == loads.drag(f["theta"], f["area"], f["p_w"], f["tau"])
    # the history was built from the frame's own drag, through 9-digit mass_kg and load_factor_g
    assert c["D_hist"] == f["row"]["mass_kg"] * f["row"]["load_factor_g"] * loads.G0
    assert abs(c["rel_patch_hist"]) <= CSV_RTOL
    assert math.isclose(c["D_patch"], readme_N, rel_tol=CSV_RTOL)
    # binning: the windward table evaluated at the valid patches
    assert abs(c["rel_table_patch"]) <= BINNING
    assert c["loads_valid"] and c["n_valid"] == int(np.sum(f["p_w"] > 0))


def test_continuum_newtonian_sphere_within_the_faceting_error(frames):
    """Frame 1 (intact): the windward pressure part against pi R^2 p_stag / 2 and the whole drag against
    0.5125 pi R^2 p_stag (README: -0.55 % for both, measured)."""
    f = frames[1]
    R = syn.SPHERE_R
    wind = f["theta"] < 0.5 * np.pi
    D_wp = loads.drag(f["theta"][wind], f["area"][wind], f["p_w"][wind], 0.0 * f["p_w"][wind])
    assert abs(D_wp / (0.5 * math.pi * R ** 2 * f["p_stag"]) - 1.0) < FACETING
    D = loads.drag(f["theta"], f["area"], f["p_w"], f["tau"])
    ans = syn.sphere_answers()
    assert abs(D / (ans["drag_over_p_stag"] * f["p_stag"]) - 1.0) < FACETING


@pytest.mark.parametrize("k", [1, 2])
def test_lee_share_and_windward_split(frames, k):
    """The split adds up; on the fixture the frame carries nothing beyond 150 deg, so the frame's own drag is all
    windward and the lee extension adds base pressure on the rear (a forward push: D_lee < 0) between the frame's
    last value 0.01 p_stag and the base 0.03 p_stag (interpolated across the 149.5-150.5 deg bin centres)."""
    f = frames[k]
    c = loads.drag_comparison(f["theta"], f["area"], f["p_w"], f["tau"], f["row"])
    assert math.isclose(c["D_windward_patch"], c["D_patch"], rel_tol=ROUNDOFF)   # summation order differs
    assert math.isclose(c["D_windward_table"] + c["D_lee"], c["D_table_lee"], rel_tol=ROUNDOFF)
    assert math.isclose(c["lee_share"], c["D_lee"] / c["D_table_lee"], rel_tol=ROUNDOFF)
    deg = np.degrees(f["theta"])
    lee = deg >= 150.0
    assert c["n_lee_patches"] == int(lee.sum()) > 0
    rear = math.fsum(f["area"][lee] * np.cos(f["theta"][lee])) * f["p_stag"]   # < 0
    assert 0.03 * rear <= c["D_lee"] <= syn.BASE_FRACTION * rear
    assert math.isclose(c["rel_table_lee_hist"], c["D_table_lee"] / c["D_hist"] - 1.0, rel_tol=ROUNDOFF)


@pytest.mark.parametrize("k", [1, 2])
def test_smooth_drag_without_derived_normals_is_the_lee_extended_table(frames, k):
    f = frames[k]
    c = loads.drag_comparison(f["theta"], f["area"], f["p_w"], f["tau"], f["row"])
    assert c["drag_normals"] == "facet"
    assert math.isclose(c["D_smooth"], c["D_table_lee"], rel_tol=ROUNDOFF)        # summation order differs
    assert math.isclose(c["D_smooth_lee"], c["D_lee"], rel_tol=ROUNDOFF)
    assert math.isclose(c["rel_smooth_hist"], c["D_smooth"] / c["D_hist"] - 1.0, rel_tol=ROUNDOFF)


def sawtooth(theta_d_deg, alpha_deg, n=40):
    """Facets about one smooth plane whose outward normal is at theta_d to v_hat = +x: per tooth one facet tilted by
    +alpha, one by -alpha (in the plane of v_hat and the smooth normal) and one untilted, areas 1, 2 and 0.5 m^2."""
    td, a = np.radians(theta_d_deg), np.radians(alpha_deg)
    unit = lambda t: np.array([np.cos(t), np.sin(t), 0.0])
    normals = np.array([unit(td + a), unit(td - a), unit(td)] * n)
    area = np.array([1.0, 2.0, 0.5] * n)
    n_d = np.tile(unit(td), (len(area), 1))
    return normals, area, n_d


def test_smooth_drag_reads_the_loads_at_the_smooth_inclination():
    """Analytic staircase (decision of 2026-10-08): the loads p_stag cos^2 theta and tau_0 sin 2 theta are what the
    surface flow gives each facet at its own theta; the untilted facets put the smooth inclination's values in the
    table exactly (theta_d is a bin centre). D_smooth is then the force of the smooth plane's loads, built here from
    vectors: pressure p(theta_d) on each facet along -n_facet, shear tau(theta_d) on the facet's area projected on the
    smooth plane along its tangent toward the tail. At the nose D_patch, the facets' own loads, falls short of it."""
    p_stag, tau0, v = 4000.0, 40.0, np.array([1.0, 0.0, 0.0])
    row = {"mass_kg": 1.0, "load_factor_g": 1.0, "p_w_stag_step_Pa": p_stag}
    for theta_d_deg in (0.5, 40.5):
        normals, area, n_d = sawtooth(theta_d_deg, 25.0)
        theta = np.arccos(np.clip(normals @ v, -1, 1))
        p_w, tau = p_stag * np.cos(theta) ** 2, tau0 * np.sin(2 * theta)
        theta_d = np.arccos(np.clip(n_d @ v, -1, 1))
        cos_to = np.einsum("ij,ij->i", normals, n_d)
        c = loads.drag_comparison(theta, area, p_w, tau, row, theta_smooth=theta_d, cos_to_smooth=cos_to)
        assert c["drag_normals"] == "derived" and c["n_lee_patches"] == 0
        td = np.radians(theta_d_deg)
        t_d = -v - (n_d @ -v)[:, None] * n_d
        t_d /= np.linalg.norm(t_d, axis=1)[:, None]            # the smooth tangent toward the tail (theta_d > 0)
        force = -(p_stag * np.cos(td) ** 2) * area[:, None] * normals \
            + (tau0 * np.sin(2 * td)) * (area * cos_to)[:, None] * t_d
        D_vec = math.fsum(force @ -v)
        assert abs(c["D_smooth"] - D_vec) <= ROUNDOFF * abs(D_vec)
        if theta_d_deg < 1.0:      # at the nose every tilt lowers cos^3; further aft the tilt toward the flow wins
            assert c["D_patch"] < c["D_smooth"]


def test_smooth_drag_arguments_go_together():
    th = np.zeros(3)
    with pytest.raises(ValueError):
        loads.drag_comparison(th, np.ones(3), np.ones(3), np.zeros(3), {"mass_kg": 1, "load_factor_g": 1,
                                                                        "p_w_stag_step_Pa": 1}, theta_smooth=th)
    with pytest.raises(ValueError):
        loads.drag_comparison(th, np.ones(3), np.ones(3), np.zeros(3), {"mass_kg": 1, "load_factor_g": 1,
                                                                        "p_w_stag_step_Pa": 1},
                              theta_smooth=np.zeros(2), cos_to_smooth=np.ones(2))


def test_uniform_pressure_on_the_closed_surface_gives_no_force(frames):
    for f in frames:
        A_total = f["area"].sum()
        resultant = (f["area"][:, None] * f["normal"]).sum(axis=0)
        assert np.all(np.abs(resultant) <= ROUNDOFF * A_total)
        p = 1234.5
        D = loads.drag(f["theta"], f["area"], np.full(len(f["area"]), p), np.zeros(len(f["area"])))
        assert abs(D) <= ROUNDOFF * p * A_total


# ------------------------------------------------------------------------------------------------ the lee extension
def test_with_lee_sets_exactly_the_bins_beyond_the_last_valid_one(frames):
    f = frames[1]
    t = table_of(f)
    lt = loads.with_lee(t)
    beyond = t.theta_deg > t.theta_last_deg
    assert np.array_equal(lt.lee, beyond)
    np.testing.assert_array_equal(lt.p[~beyond], t.p[~beyond])            # NaN interior bins stay NaN
    np.testing.assert_array_equal(lt.tau[~beyond], t.tau[~beyond])
    np.testing.assert_array_equal(lt.p[beyond], loads.BASE_FRACTION * f["p_stag"])
    np.testing.assert_array_equal(lt.tau[beyond], loads.LEE_SHEAR * t.tau[t.theta_deg == t.theta_last_deg][0])
    np.testing.assert_array_equal(lt.area, t.area)
    assert lt.lee_params == {"base_fraction": 0.03, "separation_deg": 180.0, "shear_factor": 0.5}
    assert not t.lee.any() and t.lee_params is None                        # the input is not modified


def synthetic_table(tau_last=40.0):
    """Valid bins 0-99 deg (p = 1000 cos^2, tau = tau_last at the last), nothing beyond."""
    theta = np.radians(np.arange(0.5, 100.0, 1.0))
    p = np.maximum(1000.0 * np.cos(theta) ** 2, 5.0)
    tau = np.full(len(theta), tau_last)
    return loads.build_table(theta, np.ones(len(theta)), p, tau, p_stag=1100.0)


def test_with_lee_brackets_separation_and_shear():
    t = synthetic_table()
    assert t.theta_last_deg == 99.5
    lt = loads.with_lee(t, base_fraction=0.05, separation_deg=150.0, shear_factor=1.0)
    beyond = t.theta_deg > 99.5
    np.testing.assert_array_equal(lt.p[beyond], 0.05 * 1100.0)
    sheared = beyond & (t.theta_deg <= 150.0)
    np.testing.assert_array_equal(lt.tau[sheared], 40.0)
    np.testing.assert_array_equal(lt.tau[t.theta_deg > 150.0], 0.0)
    assert sheared.sum() == 50 and (t.theta_deg > 150.0).sum() == 30
    np.testing.assert_array_equal(loads.with_lee(t, shear_factor=0.0).tau[beyond], 0.0)
    # evaluation: interpolation on the non-empty centres, held beyond the last one
    p, tau = loads.evaluate(lt, np.array([99.5, 100.0, 175.0, 180.0]))
    assert p[0] == t.p[99] and math.isclose(p[1], 0.5 * (t.p[99] + 55.0))
    assert p[2] == p[3] == 55.0 and tau[3] == 0.0


def test_bin_index_edges():
    deg = np.array([-1.0, 0.0, 0.999, 1.0, 89.9999, 90.0, 179.5, 180.0, 181.0])
    np.testing.assert_array_equal(loads.bin_index(deg), [0, 0, 0, 1, 89, 90, 179, 179, 179])


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        loads.build_table(np.zeros(3), np.ones(2), np.ones(3), np.ones(3), 1.0)
    with pytest.raises(ValueError):
        loads.build_table(np.zeros(3), np.ones(3), np.ones(3), np.ones(3), 1.0, edges=[0.0, 90.0, 45.0])


def test_valid_mask_ignores_nan_and_zero():
    p = np.array([0.0, np.nan, 10.0, 10.0, -1.0])
    tau = np.array([0.0, 0.0, 1.0, np.nan, 0.0])
    np.testing.assert_array_equal(loads.valid_patches(p, tau), [False, False, True, False, False])


# ------------------------------------------------------------------------------------------------ files
def test_stack_round_trip_and_header(frames, tmp_path):
    tables = [loads.with_lee(table_of(f)) for f in frames]
    arrays = loads.stack(tables)
    assert arrays["p"].shape == (3, 180) and arrays["lee"].dtype == np.int8
    path = tmp_path / "loads.npz"
    np.savez(path, **arrays)
    with np.load(path, allow_pickle=False) as z:
        back = loads.unstack({k: z[k] for k in z.files}, lee_params=tables[0].lee_params)
    for a, b in zip(tables, back):
        for name in ("theta_deg", "p", "tau", "area", "lee", "edges_deg"):
            np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
        assert (a.theta_last_deg == b.theta_last_deg) or (math.isnan(a.theta_last_deg) and math.isnan(b.theta_last_deg))
        assert a.p_stag == b.p_stag
    h = loads.table_header(tables[1])
    assert h["lee_label"] == "Step 4 §5's lee model, to be replaced by its module"
    assert h["lee_reference_pressure"] == "p_w_stag_step_Pa" and len(h["edges_deg"]) == 181
