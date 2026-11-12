"""The gas-shear check (the kinked-profile dispersion's Task 3, Asha's decision 7; the carry-over item "MELT SURFACE
SPEED vs GAS WALL SHEAR"): do the skins' top-cell speeds in a saved flight state follow from the forcing and the cells'
own viscosities?

    "$PY" -m analysis.kinked_dispersion.shear_check --state <flight>.state.pkl      # from Reentry Model/

On the spec's columns (section 2: angle from the stagnation point below 60 degrees, sheared depth 0.12 to 0.35 of the
flowing depth) it predicts each top cell's steady speed from the face balance of `skin.momentum_step`. In the steady
state the stress at a face is what lies above it: the gas shear on the top, times A_top / A, plus the push per unit
volume times the face's depth. A cell's half-thickness resists it over dz / (2 eta), so the top cell's speed is

    u_0 = sum over the moving cells i of  dz_i / (2 eta_i) (sigma_(i+1/2) + sigma_(i-1/2) [i > 0]),

with sigma the stress along the top cell's direction at the face below (i + 1/2) and above (i - 1/2) each cell, down
to the rigid front (a floating crust's cells resist nothing). With the gas shear alone this is
tau (A_top / A) (dz_0 / (2 eta_0) + sum over i >= 1 of dz_i / eta_i). The gas shear and the push are the last macro
step's end forcing (`body._forcing[1]`; `body.last_flow.G` along the meridian if there is none), the push from
`skin.push` as the flow sub-step takes it, and the viscosities from `skin.cell_shear_rates` on the present velocities.
It prints the 10th, 50th and 90th percentiles of the actual over the predicted top speed, with the gas shear alone and
with the push, a decomposition by bands of the ratio with the push, and a self-check: the same prediction against
`momentum_step` stepped to its steady state on these columns under the same forcing (no transport, no thickness step),
which tests the probe's arithmetic rather than the skins. The probe reports; it does not decide."""
import argparse

import numpy as np

from reentry_model.shared import rheology
from reentry_model.skyscrapers import restart, skin, skinbody

ANGLE_MAX_DEG = 60.0             # the spec's windward columns
SHEARED_RANGE = (0.12, 0.35)     # sheared depth over flowing depth, the spec's 622 thick columns
STEADY = 1.0e9                   # s: a backward-Euler step this long gives the steady profile (the skin tests' STEADY)
STEADY_STEPS = 300               # steady steps at most in the self-check, to settle the viscosities (d050: about 200)
PERCENTILES = (10, 50, 90)


def select_columns(body):
    """The spec's columns, (n,) bool: angle from the stagnation point below ANGLE_MAX_DEG and the resolved sheared depth
    (`skinbody.sheared_depth`) within SHEARED_RANGE of the flowing depth (to T_rigid)."""
    sk = body.skins
    flowing_depth = sk.depth_above(body.material.T_rigid)
    sheared = skinbody.sheared_depth(sk)
    with np.errstate(invalid="ignore", divide="ignore"):
        share = sheared / flowing_depth
    angle = np.degrees(body.load_angles())
    return (angle < ANGLE_MAX_DEG) & (share >= SHEARED_RANGE[0]) & (share <= SHEARED_RANGE[1])


def forcing(body):
    """The gas shear [Pa, (n, 3)], the push [N/m3, (n, 3)], the push's terms {name: (n, 3)} and which source gave the
    push, as the flow sub-step takes them (`skinbody.SkinBody._advance_part`, `skin.flow_substeps`)."""
    sk, mat, col = body.skins, body.material, body.columns
    end = body._forcing[1] if getattr(body, "_forcing", None) is not None else None
    if end is not None:
        shear = skin.along_meridian(col.normal, body.v_hat, end.tau)
        sigma = mat.liquid.sigma if mat.liquid is not None else 0.0
        G, terms = skin.push(col, sk.thickness(), end.dp_ds, end.g, body.v_hat, sk.rho_l, sigma,
                             wet=skin.flowing_tops(sk))
        return shear, G, terms, "body._forcing[1] (skin.push)"
    flow = body.last_flow
    shear = skin.along_meridian(col.normal, body.v_hat, flow.tau)
    G = skin.along_meridian(col.normal, body.v_hat, flow.G)
    return shear, G, {"last_flow.G": G}, "body.last_flow.G (along the meridian)"


def predict(sk, shear, G, A_top, eta, direction=None):
    """Each column's top-cell speed along its direction [m/s] (actual), and the steady face balance's prediction of it
    from the gas shear alone and from the shear and the push together, (n,) each, at the viscosities eta [Pa s, (n, K)].
    `direction` (n, 3), default the top cell's own, is the unit vector the stresses are projected on."""
    n, K = sk.m.shape
    U0 = sk.U[:, 0]
    speed = np.linalg.norm(U0, axis=1)
    if direction is None:
        direction = np.divide(U0, speed[:, None], out=np.zeros_like(U0), where=speed[:, None] > 0.0)
    flowing = sk.valid & np.isfinite(skin.cell_viscosity(sk, 0.0))
    moving = sk.valid & ~skin.anchored_cells(sk, flowing)
    dz = sk.dz()
    half = np.where(moving & np.isfinite(eta), 0.5 * dz / np.where(np.isfinite(eta), eta, 1.0), 0.0)   # m / (Pa s)
    tau_along = np.einsum("ij,ij->i", shear, direction) * A_top / sk.area          # Pa, on every face
    G_along = np.einsum("ij,ij->i", G, direction)                                  # N/m3
    below = np.cumsum(dz, axis=1)                                                  # each cell's bottom face's depth
    above = below - dz
    first = np.arange(K)[None, :] > 0
    weight = (half * (1.0 + first)).sum(axis=1)                                    # dz_0/2eta_0 + sum dz_i/eta_i
    depth_weight = (half * (below + np.where(first, above, 0.0))).sum(axis=1)
    u_shear = tau_along * weight
    actual = np.einsum("ij,ij->i", U0, direction)
    return actual, u_shear, u_shear + G_along * depth_weight, direction


def steady_check(sk, idx, shear, G, A_top):
    """The prediction against `momentum_step` stepped to its steady state on the columns idx (copies) under the same
    gas shear and push: (n_idx,) the steady top speed over the predicted, both along the steady top cell's direction,
    at the steady state's own `cell_shear_rates`; the steps taken; and the columns in which some speed still changed by
    more than 1e-10 of itself at the last step."""
    sub = sk.rows(idx)
    for step in range(1, STEADY_STEPS + 1):
        before = sub.U.copy()
        skin.momentum_step(sub, STEADY, G[idx], shear[idx], A_top[idx])
        unsettled = int((np.abs(sub.U - before) > 1e-10 * np.abs(sub.U)).any(axis=(1, 2)).sum())
        if unsettled == 0:
            break
    eta = skin.cell_viscosity(sub, skin.cell_shear_rates(sub, np.linalg.norm(shear[idx], axis=1)))
    actual, _, predicted, _ = predict(sub, shear[idx], G[idx], A_top[idx], eta)
    return actual / predicted, step, unsettled


def _band_rows(r_push, edges):
    return [(lo, hi, (r_push >= np.percentile(r_push, lo)) & (r_push <= np.percentile(r_push, hi)))
            for lo, hi in zip(edges[:-1], edges[1:])]


def report(path):
    state = restart.load(path)
    body = state["run"].body
    sk, mat = body.skins, body.material
    print("state: {}  (t = {} s, step {})".format(path, state.get("time_s"), state.get("step")))
    shear, G, terms, source = forcing(body)
    print("push from: {}".format(source))
    A_top = np.asarray(body.surface.areas, dtype=float)
    tau = np.linalg.norm(shear, axis=1)
    gamma = skin.cell_shear_rates(sk, tau)
    eta = skin.cell_viscosity(sk, gamma)
    sel = select_columns(body)
    idx = np.flatnonzero(sel)
    print("columns selected: {} of {} (angle < {:g} deg, sheared depth / flowing depth in [{:g}, {:g}])".format(
        len(idx), sk.n, ANGLE_MAX_DEG, *SHEARED_RANGE))
    if len(idx) == 0:
        return
    actual, u_shear, u_push, direction = predict(sk, shear, G, A_top, eta)
    with np.errstate(divide="ignore", invalid="ignore"):          # the columns not selected may not move
        ratio_shear, ratio_push, shear_share = actual / u_shear, actual / u_push, u_shear / u_push
    r_shear, r_push = ratio_shear[idx], ratio_push[idx]
    pct = lambda x: "  ".join("p{} {:.4g}".format(p, v) for p, v in zip(PERCENTILES, np.percentile(x, PERCENTILES)))
    print()
    print("actual / predicted top speed")
    print("  gas shear alone:     " + pct(r_shear))
    print("  gas shear and push:  " + pct(r_push))
    print("  (the spec's gate: median with the push within [0.8, 1.25])")

    # what the prediction is made of, per column
    flowing_depth = sk.depth_above(mat.T_rigid)
    tau_along = np.einsum("ij,ij->i", shear, direction)
    G_along = np.einsum("ij,ij->i", G, direction)
    term_along = {k: np.einsum("ij,ij->i", v, direction) for k, v in terms.items()}
    cos = np.einsum("ij,ij->i", shear, direction) / np.where(tau > 0.0, tau, 1.0)
    angle_off = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    eta_rest = skin.cell_viscosity(sk, 0.0)
    _, u_shear_rest, u_push_rest, _ = predict(sk, shear, G, A_top, eta_rest)
    flowing = sk.valid & np.isfinite(eta)
    bridge = flowing & (sk.T > rheology.BRIDGE_LOW)
    dz = sk.dz()
    resist = np.where(flowing, dz / np.where(flowing, eta, 1.0), 0.0)
    bridge_share = (resist * bridge).sum(axis=1) / np.where(resist.sum(axis=1) > 0, resist.sum(axis=1), 1.0)
    eta_flow_max = np.where(flowing, eta, -np.inf).max(axis=1)
    med = lambda x, rows: float(np.median(x[idx][rows]))
    rng = lambda x, rows: "{:.3g}-{:.3g}".format(float(np.min(x[idx][rows])), float(np.max(x[idx][rows])))

    print()
    print("with the viscosities at rest instead (the spec's rough read, a lower bound on the predicted speed)")
    print("  gas shear alone:     " + pct(actual[idx] / u_shear_rest[idx]))
    print("  gas shear and push:  " + pct(actual[idx] / u_push_rest[idx]))

    print()
    print("decomposition by bands of the ratio with the push (medians unless a range is given)")
    edges = (0, 10, 50, 90, 100)
    for lo, hi, rows in _band_rows(r_push, edges):
        print("  p{}-p{}: {} columns".format(lo, hi, int(rows.sum())))
        print("    ratio, gas shear alone {:.3g}; with the push {:.3g}".format(med(ratio_shear, rows),
                                                                             med(ratio_push, rows)))
        print("    top speed actual {:.3g} m/s; predicted from the gas shear {:.3g}, with the push {:.3g}".format(
            med(actual, rows), med(u_shear, rows), med(u_push, rows)))
        print("    gas shear's share of the predicted speed {:.3g} (range {})".format(
            med(shear_share, rows), rng(shear_share, rows)))
        print("    gas shear {:.3g} Pa (range {}), along the top cell's direction {:.3g} Pa; top velocity {:.3g} deg "
              "off the shear's direction; A_top/A {:.3g}".format(med(tau, rows), rng(tau, rows), med(tau_along, rows),
                                                               med(angle_off, rows), med(A_top / sk.area, rows)))
        parts = ", ".join("{} {:.3g}".format(k, med(v, rows) / 1e3) for k, v in term_along.items())
        print("    push along the top cell's direction {:.3g} kPa/m: {}".format(med(G_along, rows) / 1e3, parts))
        print("    flowing depth {:.3g} um; top cell {} K; bridge cells' share of the resistance {:.3g}".format(
            med(flowing_depth, rows) * 1e6, rng(sk.T[:, 0], rows), med(bridge_share, rows)))
        print("    viscosity: top cell {} mPa s (rest {}), largest flowing cell's {} Pa s; top cell's shear rate {} "
              "1/s".format(rng(eta[:, 0] * 1e3, rows), rng(eta_rest[:, 0] * 1e3, rows), rng(eta_flow_max, rows),
                           rng(gamma[:, 0], rows)))

    print()
    check, steps, unsettled = steady_check(sk, idx, shear, G, A_top)
    off = np.abs(check - 1.0)
    print("self-check: momentum_step stepped to steady on these columns under the same forcing ({} steps, {} columns "
          "still changing by more than 1e-10), its top speed over this probe's prediction, less 1, in size: "
          "p50 {:.2g}, p90 {:.2g}, largest {:.2g} (SHEAR_TOL {:g})".format(
              steps, unsettled, *np.percentile(off, [50, 90]), off.max(), skin.SHEAR_TOL))
    end = body._forcing
    if end is not None:
        dtau = end[1].tau[idx] / np.where(end[0].tau[idx] != 0.0, end[0].tau[idx], np.nan) - 1.0
        ddp = end[1].dp_ds[idx] / np.where(end[0].dp_ds[idx] != 0.0, end[0].dp_ds[idx], np.nan) - 1.0
        print("forcing over the last macro step (end over start, less 1, median): gas shear {:.3g}, pressure gradient "
              "{:.3g}".format(float(np.nanmedian(dtau)), float(np.nanmedian(ddp))))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--state", required=True, help="a flight's saved state (restart.save)")
    report(p.parse_args(argv).state)


if __name__ == "__main__":
    main()
