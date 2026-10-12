"""Sub-plan 20 Task 9: the independent one-dimensional reference of level 2 (layered-skin design section 11.1).

The centre of the level-2 plate (and the middle of sub-plan 18's fine box) is one-dimensional: a slab of `length`
heated at its front by a constant flux, adiabatic behind, no radiation, liquid removed as it forms. This solves it on a
fixed grid of `cell` cells with its own code -- backward Euler in the enthalpy, Newton in the temperature with a banded
solve, and the contiguous front cells at or above T_feed removed after every step -- sharing nothing with
`reentry_model.skin` but the material.

    python analysis/ablation_1d.py [--cell-um 1] [--dt-ms 0.5] [--q 6e6] [--t0 700] [--length-mm 8] [--t-end 2] [--out PATH]
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.linalg import solve_banded

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reentry_model import material  # noqa: E402


def run(mat, q=6.0e6, T0=700.0, length=8.0e-3, cell=1.0e-6, dt=5.0e-4, t_end=2.0, every=0.05):
    """Recession [m] of the slab's front at multiples of `every` [s], as (times, recession) arrays."""
    n = int(round(length / cell))
    T = np.full(n, T0)
    front = 0                                            # index of the first cell still in the slab
    m = mat.rho * cell                                   # kg per m2 per cell
    times, recession = [0.0], [0.0]
    steps = int(round(t_end / dt))
    record = max(1, int(round(every / dt)))
    for step in range(1, steps + 1):
        Tc = T[front:]
        k_old = Tc.size
        if k_old == 0:
            break
        H_old = m * mat.enthalpy(Tc)
        Tk = Tc.copy()
        for _ in range(60):
            k = mat.k(Tk)
            G = 1.0 / (0.5 * cell / k[:-1] + 0.5 * cell / k[1:])          # W/m2/K between neighbours
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
        liquid = np.cumprod(Tk >= mat.T_feed)
        front += int(liquid.sum())
        if step % record == 0:
            times.append(step * dt)
            recession.append(front * cell)
    return np.array(times), np.array(recession)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cell-um", type=float, default=1.0)
    p.add_argument("--dt-ms", type=float, default=0.5)
    p.add_argument("--q", type=float, default=6.0e6)
    p.add_argument("--t0", type=float, default=700.0)
    p.add_argument("--length-mm", type=float, default=8.0)
    p.add_argument("--t-end", type=float, default=2.0)
    p.add_argument("--out", default="ablation_1d.json")
    a = p.parse_args()
    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
    t, r = run(mat, a.q, a.t0, a.length_mm * 1e-3, a.cell_um * 1e-6, a.dt_ms * 1e-3, a.t_end)
    doc = {"q_Wm2": a.q, "T0_K": a.t0, "length_mm": a.length_mm, "cell_um": a.cell_um, "dt_ms": a.dt_ms,
           "material": mat.name, "t_s": t.tolist(), "recession_mm": (r * 1e3).tolist()}
    json.dump(doc, open(a.out, "w"), indent=1)
    for t_ in (1.0, 1.5, 2.0):
        i = int(np.argmin(np.abs(t - t_)))
        print("t = {:.2f} s: recession {:.4f} mm".format(t[i], r[i] * 1e3))


if __name__ == "__main__":
    main()
