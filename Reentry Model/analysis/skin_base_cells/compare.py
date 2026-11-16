#!/usr/bin/env python3
"""Read the single-column runs (`column.py`) of a folder against its reference: the rigid front's advance and the
bottom cell's temperature against the reference's over the same layer (its mean enthalpy's temperature).

    python -I -B analysis/skin_base_cells/compare.py [--dir output/skin_base_cells]

Each case X (A, B) needs X_ref_5um_1ms.json; every X_ref_*.json is listed as a convergence check and every X_skins_*.json
or X_*.json written by `column.py skins` is compared.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from reentry_model.shared import material  # noqa: E402

MAT = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
TIMES = (0.1, 0.2, 0.5, 1.0, 2.0, 3.0)


def at(rows, t):
    return next(r for r in rows if abs(r["t"] - t) < 1e-6)


def layer_T(ref_row, z0, z1):
    """The reference's temperature of the layer z0..z1 mm below the top: its mean enthalpy's temperature."""
    z, T = np.array(ref_row["profile"]["z_mm"]), np.array(ref_row["profile"]["T"])
    sel = (z >= z0) & (z <= z1)
    if not sel.any():
        return float("nan")
    return float(MAT.temperature_from_enthalpy(np.array([MAT.enthalpy(T[sel]).mean()]))[0])


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--dir", default=os.path.join(ROOT, "output", "skin_base_cells"))
    a = p.parse_args(argv)
    load = lambda path: json.load(open(path))
    for case in ("A", "B"):
        ref_path = os.path.join(a.dir, case + "_ref_5um_1ms.json")
        if not os.path.exists(ref_path):
            continue
        ref = load(ref_path)
        r = ref["args"]
        print("\n=== case {}: pool {} mm at {} K on rigid metal at {} K, top flux {} W/m2 ===".format(
            case, r["pool_mm"], r["T_pool"], r["T_sub"], r["q"]))
        print("reference front advance [mm] and its convergence checks:")
        for path in sorted(glob.glob(os.path.join(a.dir, case + "_ref_*.json"))):
            print("  {:22s}".format(os.path.basename(path)[:-5]) + "  ".join(
                "{:.1f}s {:.3f}".format(t, at(load(path)["rows"], t)["front_mm"]) for t in TIMES))
        for path in sorted(glob.glob(os.path.join(a.dir, case + "_*.json"))):
            s = load(path)
            if s["args"]["model"] != "skins":
                continue
            c = s["closure"]
            print("\n{} (growth {}, base cell {} mm; wall {:.0f} s; closure mass {:.1e}, energy {:.1e}, mismatch at most "
                  "{:.1e}, {:.3f} passes a step, {} of {} steps unconverged)".format(
                      os.path.basename(path)[:-5], s["args"]["growth"], s["args"]["base_cell_mm"], s["wall_s"],
                      c["mass_rel"], c["energy_rel"], c["interface_mismatch_max_rel"], c["passes_per_step"],
                      c["steps_unconverged"], c["steps"]))
            print("   t s | front mm (ref, rel) | bottom cell mm | step K | T_bottom K (ref layer) | T_above K (ref layer) "
                  "| cells")
            for t in TIMES:
                row, f = at(s["rows"], t), at(ref["rows"], t)
                z1 = row["pool_mm"]
                z0 = z1 - row["dz_bottom_mm"]
                zm = z0 - row["dz_above_mm"]
                print("  {:4.1f} | {:.3f} ({:.3f}, {:+.1%}) | {:6.3f} | {:5.2f} | {:.2f} ({:.2f}) | {:.2f} ({:.2f}) | {}".format(
                    t, row["front_mm"], f["front_mm"], row["front_mm"] / f["front_mm"] - 1.0, row["dz_bottom_mm"],
                    row["step_K"], row["T_bottom"], layer_T(f, z0, z1), row["T_above"], layer_T(f, zm, z0), row["cells"]))


if __name__ == "__main__":
    main()
