"""Diagnostic (2026-10-08, not part of the prototype): the drag a surface frame's own wall loads imply, against the
trajectory's drag in the history.

    python frame_drag.py <run.csv> <vtk dir> k [k ...]      (run with cwd = a prototype tree, for reentry_model.mesh)

D_frame = sum_f A_f [(p_w - p_inf) cos(theta_f) + tau sin(theta_f)] with v_hat = +x (the flight direction in the body
frame), over every face (leeward faces carry p_w = p_inf, tau = 0). theta from (a) the staircase facet's outward
normal (what the surface flow was evaluated with), (b) the derived surface's normal (n_derived if the frame has it,
else the Taubin-smoothed normal recomputed here), with A the staircase facet area in both. Faces the frame marks as not
evaluated (closure 1 and p_w 0) are reported, and also filled by interpolating p_w(theta), tau(theta) of the evaluated
faces in the staircase theta (the flow is a function of theta given the step's state). The history's drag is
q C_D A_ref with C_D the `drag` column and A_ref the hull's projected area (pi R_t^2 with R_t transverse_radius_mm).
"""
import csv
import os
import sys

import numpy as np
import pyvista as pv

sys.path.insert(0, os.getcwd())
from reentry_model import mesh as mesh_mod  # noqa: E402

rows = list(csv.DictReader(open(sys.argv[1])))
vtk = sys.argv[2]
for k in map(int, sys.argv[3:]):
    r = rows[k]
    poly = pv.read(os.path.join(vtk, "surface_{}.vtp".format(k)))
    faces = np.asarray(poly.faces).reshape(-1, 4)[:, 1:]
    pts = np.asarray(poly.points)
    a, b, c = (pts[faces[:, i]] for i in range(3))
    n = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(n, axis=1)
    n = n / (2.0 * area)[:, None]
    cen = (a + b + c) / 3.0
    centre = pts[np.unique(faces)].mean(axis=0)
    inward = np.einsum("ij,ij->i", n, cen - centre) < 0.0           # provisional; replaced below
    p_w, tau, closure = (np.asarray(poly.cell_data[x], dtype=float) for x in ("p_w", "tau", "closure"))
    if "n_derived" in poly.cell_data:
        nd = np.asarray(poly.cell_data["n_derived"])
        inward = np.einsum("ij,ij->i", n, nd) < 0.0
        n_out = np.where(inward[:, None], -n, n)     # the facet normal, oriented by n_derived
    else:
        # orient by each face's own tetrahedron in field_<k>.vtu (its opposite vertex), as mesh.surface does, then smooth
        grid = pv.read(os.path.join(vtk, "field_{}.vtu".format(k)))
        tets = np.asarray(grid.cells_dict[pv.CellType.TETRA])
        opp = {}
        for v in range(4):
            fv = np.sort(np.delete(tets, v, axis=1), axis=1)
            for key, o in zip(map(tuple, fv), tets[:, v]):
                opp[key] = o
        o_pt = pts[np.array([opp[tuple(f)] for f in np.sort(faces, axis=1)])]
        flip = np.einsum("ij,ij->i", n, cen - o_pt) < 0.0
        inward = flip
        n_out = np.where(flip[:, None], -n, n)
        s = mesh_mod.SurfaceMesh(faces, cen, n_out, area, None, None, pts)
        nd = s.smoothed_normals()
    p_inf = float(np.min(p_w[p_w > 0])) if (p_w > 0).any() else 0.0
    cos_s, cos_d = n_out[:, 0], nd[:, 0]
    sin_s, sin_d = np.sqrt(np.clip(1 - cos_s ** 2, 0, 1)), np.sqrt(np.clip(1 - cos_d ** 2, 0, 1))
    missing = (p_w == 0.0) & (closure == 1.0)
    win_s = cos_s > 0.0
    ev = ~missing & win_s
    th_s = np.arccos(np.clip(cos_s, -1, 1))
    o = np.argsort(th_s[ev])
    fill_p = np.where(missing & win_s, np.interp(th_s, th_s[ev][o], p_w[ev][o]), np.where(missing, p_inf, p_w))
    fill_t = np.where(missing & win_s, np.interp(th_s, th_s[ev][o], tau[ev][o]), np.where(missing, 0.0, tau))

    def drag(pw, tw, cs, sn):
        lee = cs <= 0.0
        return float((area * ((pw - p_inf) * cs + np.where(lee, 0.0, tw) * sn)).sum()), \
            float((area * (pw - p_inf) * cs).sum()), float((area * np.where(lee, 0.0, tw) * sn).sum())

    q, cd, Rt = float(r["dynamic_pressure_Pa"]), float(r["drag"]), float(r["transverse_radius_mm"]) * 1e-3
    proj = float((area * np.maximum(cos_s, 0)).sum())
    print("k={} t={} s: {} faces, inward-wound {} ({:.1f} %); not evaluated {} ({:.1f} % of area)".format(
        k, r["time_s"], len(faces), int(inward.sum()), 100 * inward.mean(), int(missing.sum()), 100 * area[missing].sum() / area.sum()))
    print("  history: q {:.1f} Pa, C_D {:.4f}, shape factor {:.4f}, D = q C_D pi R_t^2 = {:.2f} N; p_w_stag_step {:.1f} Pa; "
          "Cp_stag = (p_w,stag - p_inf)/q = {:.3f}".format(q, cd, float(r["drag_shape_factor"]), q * cd * np.pi * Rt ** 2,
                                                         float(r["p_w_stag_step_Pa"]), (float(r["p_w_stag_step_Pa"]) - p_inf) / q))
    print("  staircase projected area {:.4e} m2 vs pi R_t^2 {:.4e}".format(proj, np.pi * Rt ** 2))
    for label, pw, tw in (("as written", p_w, tau), ("missing faces filled", fill_p, fill_t)):
        ds = drag(pw, tw, cos_s, sin_s)
        dd = drag(pw, tw, cos_d, sin_d)
        print("  {:22s} staircase normals: D {:.2f} N (pressure {:.2f}, shear {:.2f}); derived normals: D {:.2f} N "
              "(pressure {:.2f}, shear {:.2f})".format(label, *ds, *dd))
    # p_w, tau re-read at the derived theta (the flow as a function of theta, from the evaluated faces)
    th_d = np.arccos(np.clip(cos_d, -1, 1))
    pw_d = np.where(cos_d > 0, np.interp(th_d, th_s[ev][o], p_w[ev][o]), p_inf)
    tw_d = np.where(cos_d > 0, np.interp(th_d, th_s[ev][o], tau[ev][o]), 0.0)
    print("  {:22s} derived normals: D {:.2f} N (pressure {:.2f}, shear {:.2f})".format("p_w(theta_derived)", *drag(pw_d, tw_d, cos_d, sin_d)))
    win = cos_s > 0.0
    proj_p = float((np.where(win, (pw_d - p_inf) * area * cos_s, 0.0)).sum())
    print("  {:22s} on the facets' exact projected areas A cos(theta_stair): pressure drag {:.2f} N".format("p_w(theta_derived)", proj_p))
    print("  intact-sphere reference: q Cp_stag/2 pi R_t^2 = {:.2f} N (modified Newtonian pressure drag)".format(
        q * (float(r["p_w_stag_step_Pa"]) - p_inf) / q / 2 * np.pi * Rt ** 2))
