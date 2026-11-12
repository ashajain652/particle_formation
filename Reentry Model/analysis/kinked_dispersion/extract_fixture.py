"""The d050 fixture of the kinked-profile dispersion (Task 4; spec stage 3): `skinbody.layer_columns` of 36 columns of a
saved flight state, for the tests of the shared straight-segment fit (`tests/test_reentry_model_melt_profile.py`) and the
viscous check (Task 10). The state itself is git-ignored and on this laptop only; the fixture is committed.

    "$PY" -m analysis.kinked_dispersion.extract_fixture --state <flight>.state.pkl \\
        --out data/reference_values/d050_layer_columns.npz                          # from Reentry Model/

The columns, each set taken evenly spaced along the order named:
- 24 of the spec's thick set (`shear_check.select_columns`: angle from the stagnation point below 60 degrees, sheared
  depth 0.12 to 0.35 of the flowing depth), in order of angle;
- 6 plug columns, whose fastest flowing cell along the top cell's direction lies below the top cell, in order of that
  peak over the top cell's speed;
- 6 other wet windward columns (top cell flowing, angle below 90 degrees, in neither set above), in order of angle.

The npz holds the `melt_profile.LayerColumns` fields, cut to the deepest selected column's cells, and beside them
`theta` [rad] (`body.load_angles()`) and `tau` [Pa], the size of the gas shear as the flow sub-step takes it
(`shear_check.forcing`), from which `layer_columns` took the shear rates. Loading filters to the dataclass's fields.

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
WINDWARD_DEG = 90.0                 # the wet windward columns: top cell flowing, angle below this
UNSETTLED_PASSES = 200              # the reference count of fixed-point passes for "unsettled" (Task 3's)


def evenly(idx, key, k):
    """k of the columns idx, evenly spaced along `key` (both (m,)), in column order."""
    if len(idx) <= k:
        return np.sort(idx)
    order = idx[np.argsort(key, kind="stable")]
    return np.sort(order[np.round(np.linspace(0, len(idx) - 1, k)).astype(int)])


def plug_columns(cols, candidates):
    """(the plug columns among candidates, their fastest flowing cell's speed over the top cell's)."""
    peak = np.where(cols.flowing, cols.u, -np.inf).max(axis=1)
    top = cols.u[:, 0]
    plug = candidates & (peak > top)
    return np.flatnonzero(plug), (peak / np.where(top > 0.0, top, 1.0))[plug]


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
    plug, peak = plug_columns(cols, wet)
    is_plug = np.zeros(sk.n, dtype=bool)
    is_plug[plug] = True
    other = np.flatnonzero(wet & ~thick & ~is_plug)
    t_idx = np.flatnonzero(thick)
    chosen_thick = evenly(t_idx, theta[t_idx], THICK)
    chosen_plug = evenly(plug, peak, PLUG)
    chosen_other = evenly(other, theta[other], OTHER)
    idx = np.concatenate([chosen_thick, chosen_plug, chosen_other])
    print("wet windward columns (top cell flowing, angle below {:g} deg): {}; thick set: {} ({} wet windward); plug "
          "columns: {} (fastest flowing cell over the top's: {:.5f} to {:.4f}); others: {}".format(
              WINDWARD_DEG, int(wet.sum()), int(thick.sum()), int((thick & wet).sum()), len(plug), peak.min(),
              peak.max(), len(other)))
    print("fixture: {} thick, {} plug, {} other columns".format(len(chosen_thick), len(chosen_plug), len(chosen_other)))

    sub = cols.rows(idx)
    K = int(sk.count[idx].max())
    arrays = {f.name: getattr(sub, f.name) for f in dataclasses.fields(melt_profile.LayerColumns)}
    arrays = {k: (v[:, :K] if v.ndim == 2 else v) for k, v in arrays.items()}
    arrays["theta"], arrays["tau"] = theta[idx], tau[idx]
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    np.savez_compressed(out_path, **arrays)
    print("wrote {}: {} columns of {} cells, {:.3f} MB".format(out_path, len(idx), K, os.path.getsize(out_path) / 1e6))

    # breakpoints over cells whose viscosity has not settled (Task 3's caution)
    unsettled = unsettled_cells(sk, tau)
    flowing = cols.flowing
    print("flowing cells (wet windward) whose viscosity has not settled in {} passes: {} of {} ({:.2g} %)".format(
        skin.SHEAR_PASSES, int((unsettled & wet[:, None]).sum()), int((flowing & wet[:, None]).sum()),
        100.0 * (unsettled & wet[:, None]).sum() / max((flowing & wet[:, None]).sum(), 1)))
    fixture = melt_profile.LayerColumns(**{k: v for k, v in arrays.items() if k not in ("theta", "tau")})
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
