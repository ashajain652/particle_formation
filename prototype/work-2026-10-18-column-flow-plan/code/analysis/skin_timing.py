"""Sub-plan 20 Task 10: the skins' cost per macro step on the default 100 mm sphere (2 mm surface, dense
band), production cells (10 um), 0.5 s macro steps, 10 ms sub-steps: three steps below melting (300 K, 1 MW/m2 at the
nose) and three melting (880 K, 2 MW/m2 at the nose), liquid removed as it forms."""
import json, sys, time
sys.path.insert(0, ".")
import numpy as np
from reentry_model import heating, material, mesh, skin, skinbody
from reentry_model.thermal import thermal_solver
out = {}
for label, T0, q0 in (("solid", 300.0, 1.0e6), ("melting", 880.0, 2.0e6)):
    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
    m = mesh.sphere_mesh(0.05, band=mesh.DEFAULT_BAND)
    t0 = time.perf_counter()
    b = skinbody.SkinBody(m, mat, thermal_solver("skfem"), mat.rho * m.volume(), skin.SkinSettings(), T0=T0, v_hat=(-1.0, 0.0, 0.0))
    setup = time.perf_counter() - t0
    steps = []
    for k in range(3):
        t1 = time.perf_counter()
        b.advance(0.5 * (k + 1), 0.5, heating.HeatingResult(q0 * np.clip(np.cos(b.theta), 0.0, None), q0, 0.0, 0.0, 0.0))
        steps.append(time.perf_counter() - t1)
    st = b.skin_stats()
    out[label] = {"columns": b.columns.n, "cells": st["skin_cells"], "n_nodes": b.mesh.n_nodes, "setup_s": setup, "step_s": steps,
                  "repeats": b.interface_repeats, "wall": b.wall, "removed_kg": b.removed_mass, "energy_residual": b.energy_balance_residual()}
    print(label, json.dumps(out[label]), flush=True)
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "skin_timing.json", "w"), indent=1)
