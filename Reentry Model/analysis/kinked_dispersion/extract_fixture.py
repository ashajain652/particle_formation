"""The d050 fixture of the kinked-profile dispersion (Task 4; spec stage 3): `skinbody.layer_columns` of 36 columns of a
saved flight state, for the tests of the shared straight-segment fit (`tests/test_reentry_model_melt_profile.py`) and the
viscous check (Task 10). The state itself is git-ignored and on this laptop only; the fixture is committed.

    "$PY" -m analysis.kinked_dispersion.extract_fixture --state <flight>.state.pkl \\
        --out data/reference_values/d050_layer_columns.npz                          # from Reentry Model/

The columns, each set taken evenly spaced along the order named:
- 24 of the spec's thick set (`shear_check.select_columns`: angle from the stagnation point below 60 degrees, sheared
  depth 0.12 to 0.35 of the flowing depth), plug columns left out, in order of angle;
- 6 plug columns, in order of their fastest flowing cell's speed over the top cell's. A plug column is a wet windward
  column whose fastest cell of the flowing run moves more than PLUG_MARGIN times as fast as its top cell, both as speed
  magnitudes (the spec builder check's rule; 14 columns on d050). Only those whose fitted top segment slopes the other
  way -- the speed rising with depth, Review Focus 1 -- are taken (8 of the 14 on d050), and the script asserts it;
  with fewer than 6 it takes all it has and says so;
- 6 other wet windward columns (top cell flowing, angle below 90 degrees, in neither set above), in order of angle.

The npz holds the `melt_profile.LayerColumns` fields, cut to the deepest selected column's cells, and beside them
`theta` [rad] (`body.load_angles()`), `tau` [Pa], the size of the gas shear as the flow sub-step takes it
(`shear_check.forcing`), from which `layer_columns` took the shear rates, and `group`, each row's set: 0 thick, 1 plug,
2 other. Loading filters to the dataclass's fields.

It prints the selection, the file's size, and two checks: how many of the fixture's breakpoints sit over a cell whose
viscosity has not settled in the skins' SHEAR_PASSES passes (the cell below a breakpoint is the one `nu_below` reads;
"unsettled" means its viscosity at SHEAR_PASSES passes differs from that at UNSETTLED_PASSES by more than SHEAR_TOL,
relatively: Task 3 found about 3 % of flowing cells so, near the 367 1/s cap of Li et al.'s slurry law), and, for
comparison with the spec's builder check (0.94), the medians of the first kink's depth and of the wall's depth over the
906.4 K isotherm's on the whole thick set."""
import argparse
import dataclasses
import os

import numpy as np

from analysis.kinked_dispersion.shear_check import forcing, select_columns
from reentry_model.shared import dispersion, melt_profile, rheology
from reentry_model.skyscrapers import restart, skin, skinbody

THICK, PLUG, OTHER = 24, 6, 6       # columns from each set
GROUP_THICK, GROUP_PLUG, GROUP_OTHER = 0, 1, 2      # the npz's `group` of each row
PLUG_MARGIN = 1.05                  # a plug column's fastest flowing cell moves faster than this times its top cell
WINDWARD_DEG = 90.0                 # the wet windward columns: top cell flowing, angle below this
UNSETTLED_PASSES = 200              # the reference count of fixed-point passes for "unsettled" (Task 3's)


def evenly(idx, key, k):
    """k of the columns idx, evenly spaced along `key` (both (m,)), in column order."""
    if len(idx) <= k:
        return np.sort(idx)
    order = idx[np.argsort(key, kind="stable")]
    return np.sort(order[np.round(np.linspace(0, len(idx) - 1, k)).astype(int)])


def plug_columns(speed, flowing, candidates):
    """(the plug columns among candidates, their fastest flowing cell's speed over the top cell's): the columns whose
    fastest cell of the flowing run moves faster than PLUG_MARGIN times the top cell, from the cells' speed magnitudes
    `speed` [m/s, (n, K)] and the flowing runs `flowing` (n, K)."""
    peak = np.where(flowing, speed, -np.inf).max(axis=1)
    top = speed[:, 0]
    plug = candidates & (peak > PLUG_MARGIN * top)
    return np.flatnonzero(plug), (peak / np.where(top > 0.0, top, 1.0))[plug]


def top_slopes(cols):
    """Each column's fitted top segment slope dU/dy [1/s, (n,)], y up: negative where the speed rises with depth."""
    return melt_profile.fit_profiles(cols).profiles.slopes()[:, 0]


def unsettled_cells(sk, tau):
    """(n, K) bool: the flowing cells whose viscosity at SHEAR_PASSES passes of `cell_shear_rates` differs by more than
    SHEAR_TOL, relatively, from that at UNSETTLED_PASSES passes."""
    eta = skin.cell_viscosity(sk, skin.cell_shear_rates(sk, tau))
    passes = skin.SHEAR_PASSES
    try:
        skin.SHEAR_PASSES = UNSETTLED_PASSES
        eta_long = skin.cell_viscosity(sk, skin.cell_shear_rates(sk, tau))
    finally:
        skin.SHEAR_PASSES = passes
    flowing = np.isfinite(eta) & np.isfinite(eta_long)
    change = np.abs(np.where(flowing, eta_long, 1.0) / np.where(flowing, eta, 1.0) - 1.0)
    return flowing & (change > skin.SHEAR_TOL)


def medians_over_isotherm(cols, d906):
    """(first kink depth / d906 (m,), wall depth / d906 (m,), kinked share) over the reducible columns with d906 > 0."""
    red = melt_profile.reducible(cols)
    fit = melt_profile.fit_profiles(cols.rows(red))
    p = fit.profiles
    wall = dispersion.place_wall(p, melt_profile.nu_below(cols.rows(red), fit.cell))
    kinks = dispersion.kinks(p)
    has = kinks.any(axis=1)
    first = np.where(has, p.depth[np.arange(p.count.size), np.argmax(kinks, axis=1)], np.nan)
    d = d906[red]
    ok = d > 0.0
    return first[ok & has] / d[ok & has], wall.depth[ok] / d[ok], float(has.mean())


def extract(state_path, out_path):
    state = restart.load(state_path)
    body = state["run"].body
    sk = body.skins
    shear, _, _, source = forcing(body)
    tau = np.linalg.norm(shear, axis=1)
    theta = np.asarray(body.load_angles(), dtype=float)
    cols = skinbody.layer_columns(sk, tau)
    print("state: {} (t = {} s, step {}); gas shear from {}".format(state_path, state.get("time_s"), state.get("step"),
                                                                   source))

    thick = select_columns(body)
    wet = skin.flowing_tops(sk) & (np.degrees(theta) < WINDWARD_DEG)
    plug, peak = plug_columns(np.linalg.norm(sk.U, axis=2), cols.flowing, wet)
    is_plug = np.zeros(sk.n, dtype=bool)
    is_plug[plug] = True
    slope = top_slopes(cols.rows(plug))
    rising = slope < 0.0                                     # the speed rises with depth below the surface
    other = np.flatnonzero(wet & ~thick & ~is_plug)
    t_idx = np.flatnonzero(thick & ~is_plug)
    chosen_thick = evenly(t_idx, theta[t_idx], THICK)
    chosen_plug = evenly(plug[rising], peak[rising], PLUG)
    chosen_other = evenly(other, theta[other], OTHER)
    idx = np.concatenate([chosen_thick, chosen_plug, chosen_other])
    assert len(np.unique(idx)) == len(idx)
    group = np.repeat([GROUP_THICK, GROUP_PLUG, GROUP_OTHER], [len(chosen_thick), len(chosen_plug), len(chosen_other)])
    print("wet windward columns (top cell flowing, angle below {:g} deg): {}; thick set: {} ({} wet windward, {} plug); "
          "plug columns: {} (fastest flowing cell over the top's, as speeds: {:.4f} to {:.4f}), {} of them with a "
          "fitted top segment rising with depth; others: {}".format(
              WINDWARD_DEG, int(wet.sum()), int(thick.sum()), int((thick & wet).sum()), int((thick & is_plug).sum()),
              len(plug), peak.min(), peak.max(), int(rising.sum()), len(other)))
    if len(chosen_plug) < PLUG:
        print("ONLY {} plug columns have a fitted top segment rising with depth; the fixture takes all of them".format(
            len(chosen_plug)))
    print("fixture: {} thick, {} plug, {} other columns".format(len(chosen_thick), len(chosen_plug), len(chosen_other)))

    sub = cols.rows(idx)
    K = int(sk.count[idx].max())
    arrays = {f.name: getattr(sub, f.name) for f in dataclasses.fields(melt_profile.LayerColumns)}
    arrays = {k: (v[:, :K] if v.ndim == 2 else v) for k, v in arrays.items()}
    fixture = melt_profile.LayerColumns(**arrays)
    fixture_slope = top_slopes(fixture)
    assert (fixture_slope[group == GROUP_PLUG] < 0.0).all(), "a plug row's fitted top segment does not rise with depth"
    by_column = dict(zip(plug, peak))
    for i, s in zip(idx[group == GROUP_PLUG], fixture_slope[group == GROUP_PLUG]):
        print("  plug row, column {}: fastest flowing cell over the top's {:.4f}; fitted top slope {:.4g} 1/s".format(
            i, by_column[i], s))
    arrays["theta"], arrays["tau"], arrays["group"] = theta[idx], tau[idx], group
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    np.savez_compressed(out_path, **arrays)
    print("wrote {}: {} columns of {} cells, {:.3f} MB".format(out_path, len(idx), K, os.path.getsize(out_path) / 1e6))

    # breakpoints over cells whose viscosity has not settled (Task 3's caution)
    unsettled = unsettled_cells(sk, tau)
    flowing = cols.flowing
    print("flowing cells (wet windward) whose viscosity has not settled in {} passes: {} of {} ({:.2g} %)".format(
        skin.SHEAR_PASSES, int((unsettled & wet[:, None]).sum()), int((flowing & wet[:, None]).sum()),
        100.0 * (unsettled & wet[:, None]).sum() / max((flowing & wet[:, None]).sum(), 1)))
    fit = melt_profile.fit_profiles(fixture)
    nu = melt_profile.nu_below(fixture, fit.cell)
    kinks = dispersion.kinks(fit.profiles)
    wall = dispersion.place_wall(fit.profiles, nu)
    under = np.zeros(fit.cell.shape, dtype=bool)
    rr, jj = np.nonzero(fit.cell >= 0)
    below = fit.cell[rr, jj] + 1
    inside = below < K
    under[rr[inside], jj[inside]] = unsettled[idx][rr[inside], below[inside]]
    at_wall = np.zeros(fit.cell.shape, dtype=bool)
    at_wall[np.flatnonzero(wall.floor), wall.index[wall.floor]] = True
    print("fixture breakpoints: {}, of which over an unsettled cell {}; kinks {}, over an unsettled cell {}; "
          "floor-like walls {}, over an unsettled cell {}".format(
              int((fit.cell >= 0).sum()), int(under.sum()), int(kinks.sum()), int((under & kinks).sum()),
              int(wall.floor.sum()), int((under & at_wall).sum())))

    # the builder check's measurement on the whole thick set
    d906 = sk.depth_above(rheology.BRIDGE_LOW)
    kink_ratio, wall_ratio, share = medians_over_isotherm(cols.rows(t_idx), d906[t_idx])
    print("thick set: {:.3g} of the columns have a kink; over the depth of the {:g} K isotherm (columns where it lies "
          "below the top cell): first kink median {:.3g} (p10 {:.3g}, p90 {:.3g}, {} columns), wall median {:.3g} "
          "(p10 {:.3g}, p90 {:.3g}, {} columns)".format(share, rheology.BRIDGE_LOW, np.median(kink_ratio),
                                                        *np.percentile(kink_ratio, [10, 90]), kink_ratio.size,
                                                        np.median(wall_ratio), *np.percentile(wall_ratio, [10, 90]),
                                                        wall_ratio.size))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--state", required=True, help="a flight's saved state (restart.save)")
    p.add_argument("--out", required=True, help="the npz to write")
    a = p.parse_args(argv)
    extract(a.state, a.out)


if __name__ == "__main__":
    main()
