"""Sub-plan 20 Task 9: the independent one-dimensional reference of level 2 (layered-skin design section 11.1).

The centre of the level-2 plate (and the middle of sub-plan 18's fine box) is one-dimensional: a slab of `length`
heated at its front by a constant flux, adiabatic behind, no radiation, liquid removed as it forms. This solves it on a
fixed grid of `cell` cells with its own code -- backward Euler in the enthalpy, Newton in the temperature with a banded
solve, and the contiguous front cells at or above T_feed removed after every step -- sharing nothing with
`reentry_model.skin` but the material. It also records the non-rigid depth: from the front to where the slab falls to
T_rigid, the layer the skins must hold (the hot plate's reference, 2026-10-14).

`--swell` (2026-10-14): each cell holds a fixed mass (`cell` of solid metal) and is `off` one density thick, as the
level-2 references are; `all` gives every cell the skins' mixture density (spec 5.2: specific volumes add, the liquid at
its own density), so heat crosses the swollen thickness as it does in a skin -- the level-1 reference of the skins'
swelling; `nonrigid` swells only the cells at or above T_rigid, the metal the skins hold above a 3D model of one
density. Lengths are reported as solid metal (cells times `cell`).

`--remove-at liquidus` (sub-plan 21 Task 8, 2026-10-24): the front cells leave at or above the liquidus (908 K) instead
of T_feed (910 K), as sub-plan 21 Task 7's release from the top cells takes the liquid. Stage 2's T_feed keeps a liquid
layer on the front across which the last 2 K fall, about 2 K k_l / q thick (0.03 mm at 6 MW/m2, 0.09 mm at 2 MW/m2),
and the slab's melting starts that much later: at the liquidus the cold plate's slab recedes 0.048, 0.060 and 0.079 mm
further by 1.0, 1.5 and 2.0 s (6.1, 2.2 and 1.5 %) and the hot plate's 0.139, 0.149 and 0.158 mm (17.6, 8.2 and 5.4 %),
its non-rigid depth 0.13-0.14 mm (9 %) thinner (1 um, 0.5 ms).

    python analysis/ablation_1d.py [--cell-um 1] [--dt-ms 0.5] [--q 6e6] [--t0 700] [--length-mm 8] [--t-end 2]
                                   [--swell off|all|nonrigid] [--remove-at feed|liquidus] [--out PATH]
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.linalg import solve_banded

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reentry_model import material  # noqa: E402


SWELL_MODES = ("off", "all", "nonrigid")
REMOVE_AT = ("feed", "liquidus")      # the front cells leave at or above T_feed (stage 2's) or the liquidus (Task 7's)


def thickness(mat, T, cell, swell):
    """Each cell's thickness [m]: `cell` of solid metal at one density, or at the mixture density (module docstring)."""
    if swell == "off" or mat.liquid is None:
        return np.full(T.shape, cell)
    f = mat.liquid_fraction(T)
    dz = cell * mat.rho * ((1.0 - f) / mat.rho + f / mat.liquid.rho)
    return dz if swell == "all" else np.where(T >= mat.T_rigid, dz, cell)


def run(mat, q=6.0e6, T0=700.0, length=8.0e-3, cell=1.0e-6, dt=5.0e-4, t_end=2.0, every=0.05, swell="off",
        remove_at="feed"):
    """Recession [m] of the slab's front and its non-rigid depth [m] (front to T_rigid), both as solid metal, at
    multiples of `every` [s], as (times, recession, non-rigid depth) arrays; the front cells leave at or above T_feed or,
    with `remove_at` "liquidus", the liquidus."""
    if swell not in SWELL_MODES:
        raise ValueError("swell must be one of {}, got {!r}".format(SWELL_MODES, swell))
    if remove_at not in REMOVE_AT:
        raise ValueError("remove_at must be one of {}, got {!r}".format(REMOVE_AT, remove_at))
    T_remove = mat.T_feed if remove_at == "feed" else mat.T_liquidus
    n = int(round(length / cell))
    T = np.full(n, T0)
    front = 0                                            # index of the first cell still in the slab
    m = mat.rho * cell                                   # kg per m2 per cell
    times, recession, nonrigid = [0.0], [0.0], [0.0]
    steps = int(round(t_end / dt))
    record = max(1, int(round(every / dt)))
    for step in range(1, steps + 1):
        Tc = T[front:]
        k_old = Tc.size
        if k_old == 0:
            break
        H_old = m * mat.enthalpy(Tc)
        Tk = Tc.copy()
        for _ in range(80):
            k = mat.k(Tk)
            dz = thickness(mat, Tk, cell, swell)                            # at this iterate's temperatures
            G = 1.0 / (0.5 * dz[:-1] / k[:-1] + 0.5 * dz[1:] / k[1:])      # W/m2/K between neighbours
            flux = G * (Tk[:-1] - Tk[1:])
            R = (m * mat.enthalpy(Tk) - H_old) / dt
            R[0] -= q
            R[:-1] += flux
            R[1:] -= flux
            diag = m * mat.cp_eff(Tk) / dt
            diag[:-1] += G
            diag[1:] += G
            ab = np.zeros((3, k_old))
            ab[0, 1:], ab[1], ab[2, :-1] = -G, diag, -G
            dT = solve_banded((1, 1), ab, -R)
            Tk = mat.temperature_from_enthalpy(mat.enthalpy(Tk) + mat.cp_eff(Tk) * dT)   # enthalpy-consistent update
            if np.abs(dT).max() < 1e-7:
                break
        else:
            raise RuntimeError("the 1D reference did not converge at step {}".format(step))
        T[front:] = Tk
        liquid = np.cumprod(Tk >= T_remove)
        front += int(liquid.sum())
        if step % record == 0:
            times.append(step * dt)
            recession.append(front * cell)
            cold = np.flatnonzero(T[front:] <= mat.T_rigid)
            nonrigid.append((cold[0] if cold.size else n - front) * cell)
    return np.array(times), np.array(recession), np.array(nonrigid)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cell-um", type=float, default=1.0)
    p.add_argument("--dt-ms", type=float, default=0.5)
    p.add_argument("--q", type=float, default=6.0e6)
    p.add_argument("--t0", type=float, default=700.0)
    p.add_argument("--length-mm", type=float, default=8.0)
    p.add_argument("--t-end", type=float, default=2.0)
    p.add_argument("--swell", choices=SWELL_MODES, default="off")
    p.add_argument("--remove-at", choices=REMOVE_AT, default="feed")
    p.add_argument("--out", default="ablation_1d.json")
    a = p.parse_args()
    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
    t, r, nr = run(mat, a.q, a.t0, a.length_mm * 1e-3, a.cell_um * 1e-6, a.dt_ms * 1e-3, a.t_end, swell=a.swell,
                   remove_at=a.remove_at)
    doc = {"q_Wm2": a.q, "T0_K": a.t0, "length_mm": a.length_mm, "cell_um": a.cell_um, "dt_ms": a.dt_ms, "swell": a.swell,
           "remove_at": a.remove_at, "material": mat.name, "t_s": [round(x, 6) for x in t.tolist()],
           "recession_mm": [round(x, 6) for x in (r * 1e3).tolist()], "nonrigid_mm": [round(x, 6) for x in (nr * 1e3).tolist()]}
    json.dump(doc, open(a.out, "w"), indent=0)
    for t_ in (1.0, 1.5, 2.0):
        i = int(np.argmin(np.abs(t - t_)))
        print("t = {:.2f} s: recession {:.4f} mm, non-rigid depth {:.4f} mm".format(t[i], r[i] * 1e3, nr[i] * 1e3))


if __name__ == "__main__":
    main()
