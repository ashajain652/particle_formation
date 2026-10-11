"""Experiment X2 diagnostic (reads the model, changes nothing): rerun one X2 configuration and report where elements
die and where the front-facing surface is, by region of the box (lateral column, outer ring, along the side walls), at
chosen times. Usage: python diag_x2.py RUN_NAME [t1 t2 ...]"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import x2_box as X  # noqa: E402  (imports the experiment copy first)

bx = X.bx


def snapshot(b, spec, t):
    lx, ly, lz, h = spec["box"]
    m, s = b.mesh, b.surface
    ce = m.points[m.tets].mean(axis=1)
    lat_el = (ce[:, 1] > h) & (ce[:, 1] < ly - h) & (ce[:, 2] > h) & (ce[:, 2] < lz - h)
    dead = ~m.active
    dmin = np.minimum.reduce([ce[:, 1], ly - ce[:, 1], ce[:, 2], lz - ce[:, 2]])      # distance to the nearest side plane
    near_side = dmin < 0.5e-3
    c = s.centroids
    lat = (c[:, 1] > h * (1 - 1e-9)) & (c[:, 1] < ly - h * (1 - 1e-9)) & (c[:, 2] > h * (1 - 1e-9)) & (c[:, 2] < lz - h * (1 - 1e-9))
    fwd = np.maximum(0.0, -s.normals[:, 0]) * s.areas
    walls = X.wall_patches(s, spec["box"])
    on = (b.skin_of_patch >= 0) if b.skins is not None and b.skins.n else np.zeros(s.n_patches, dtype=bool)
    ns = b.smoothed_normals_x2() if spec["heating"] == "smoothed" else s.smoothed_normals()
    w, _ = X.smoothed_weights(s, ns)
    W = w * s.areas
    q, _ = X.heating_flux(b, spec)
    P = q * s.areas / (spec["q"] * ly * lz)
    band = ~lat & ((np.minimum.reduce([c[:, 1], ly - c[:, 1], c[:, 2], lz - c[:, 2]])) < 2e-3 * (1 - 1e-9))
    front_x = np.full(s.n_patches, np.nan)
    out_power = {"lateral": float(P[lat].sum()), "outer ring, x < 1 mm": float(P[band & (c[:, 0] < 1e-3)].sum()),
                 "outer ring, 1 mm <= x < 4 mm": float(P[band & (c[:, 0] >= 1e-3) & (c[:, 0] < 4e-3)].sum()),
                 "anywhere at x >= 4 mm (back half)": float(P[c[:, 0] >= 4e-3].sum()),
                 "lateral, x >= 4 mm": float(P[lat & (c[:, 0] >= 4e-3)].sum())}
    out = {"t": t, "power_share_by_region": out_power, "dead_total": int(dead.sum()), "dead_lateral_column": int((dead & lat_el).sum()),
           "dead_outer_ring": int((dead & ~lat_el).sum()), "dead_within_0.5mm_of_a_side_plane": int((dead & near_side).sum()),
           "deepest_dead_x_mm_lateral": float(ce[dead & lat_el, 0].max() * 1e3) if (dead & lat_el).any() else 0.0,
           "deepest_dead_x_mm_outer": float(ce[dead & ~lat_el, 0].max() * 1e3) if (dead & ~lat_el).any() else 0.0,
           "front_facing_projected_area_mm2": {"lateral": float(fwd[lat].sum() * 1e6), "outer": float(fwd[~lat].sum() * 1e6),
                                               "total": float(fwd.sum() * 1e6)},
           "skins": {"total": int(on.sum()), "on_box_walls": int((on & walls).sum()), "lateral": int((on & lat).sum()),
                     "outer_not_walls": int((on & ~lat & ~walls).sum())},
           "deaths_with_skins_on_walls_note": "walls = the side planes and the back plane"}
    return out


def main():
    name = sys.argv[1]
    times = [float(x) for x in sys.argv[2:]] or [0.5, 1.0, 1.5, 2.0]
    spec = X.specs()[name]
    b, mat = X.make_body(spec)
    snaps = []
    dt = spec["dt"]
    lx, ly, lz, h = spec["box"]
    for k in range(int(round(spec["t_end"] / dt))):
        t = (k + 1) * dt
        q, _ = X.heating_flux(b, spec) if spec["heating"] != "projected" else (None, None)
        b.advance(t, dt, bx.heating.HeatingResult(q, spec["q"], 0.0, 0.0, 0.0))
        L = b.m_f + b.m_d
        b.removed_mass += float(L.sum())
        b.removed_enthalpy += float((L * b.film_enthalpy()).sum())
        b.m_f[:] = 0.0
        b.m_d[:] = 0.0
        b.solver.set_film_mass(b._film_nodal())
        if any(abs(t - x) < 1e-9 for x in times):
            snaps.append(snapshot(b, spec, t))
            print(json.dumps(snaps[-1]))
    os.makedirs(os.path.join(X.RUNS, "diag"), exist_ok=True)
    with open(os.path.join(X.RUNS, "diag", name + ".json"), "w") as f:
        json.dump({"name": name, "snapshots": snaps, "residual": b.energy_balance_residual()}, f, indent=1)


if __name__ == "__main__":
    main()
