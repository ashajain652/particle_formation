"""Picture: each surface face's own normal against its Taubin-smoothed normal, on the receded front of the box.
Runs the q6p box case (6 MW/m2 by projection, 700 K, Scheil, skins as built, dev preset, 0.05 s) to T_STOP on the Kuhn
cube box and on prism-layer boxes, then draws a slice z in [4, 6] mm projected on the x-y plane."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, "..", "x2", "code"))
from reentry_model import body, heating, material, mesh, skin, thermal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
T_STOP, DT, Q = 1.1, 0.05, 6.0e6
def run(msh):
    mat = material.Material.from_drama_json("AA7075_scheil")
    st = body.MeltSettings(removal="girin", size_feedback="initial", surface_model="skin", skin=skin.SkinSettings(cell=20e-6, substep=0.01))
    b = body.MeltingBody(msh, mat, thermal.thermal_solver("skfem", linear_solver="direct"), float(mat.rho * msh.volume()), settings=st, T0=700.0, emissivity=0.0)
    A0 = 1e-4
    for k in range(int(round(T_STOP / DT))):
        s = b.surface; w = np.maximum(0.0, -s.normals[:, 0]); q = Q * A0 * w / float((w * s.areas).sum())
        b.advance((k + 1) * DT, DT, heating.HeatingResult(q, Q, 0.0, 0.0, 0.0))
        L = b.m_f + b.m_d; b.removed_mass += float(L.sum()); b.removed_enthalpy += float((L * b.film_enthalpy()).sum())
        b.m_f[:] = 0; b.m_d[:] = 0; b.solver.set_film_mass(b._film_nodal())
    return b
cases = [("Kuhn cubes, 2 mm", mesh.box_mesh(8e-3, 10e-3, 10e-3, 2e-3)),
         ("prism layers 1.0 mm (2 mm triangles)", mesh.layered_box_mesh(8e-3, 10e-3, 10e-3, 2e-3, 1.0e-3)),
         ("prism layers 0.5 mm (2 mm triangles)", mesh.layered_box_mesh(8e-3, 10e-3, 10e-3, 2e-3, 0.5e-3))]
fig, axes = plt.subplots(1, 3, figsize=(16, 6.2))
stats = []
for ax, (title, msh) in zip(axes, cases):
    b = run(msh); s = b.surface; ns = s.smoothed_normals()
    c, n = s.centroids * 1e3, s.normals
    front = (c[:, 0] < 7.9) & (c[:, 1] > 0.05) & (c[:, 1] < 9.95) & (c[:, 2] > 0.05) & (c[:, 2] < 9.95)   # not the sides/back
    sl = front & (c[:, 2] > 4.0) & (c[:, 2] < 6.0)
    for f in np.flatnonzero(sl):
        tri = msh.points[s.faces[f]] * 1e3
        ax.fill(tri[:, 0], tri[:, 1], facecolor=(0.85, 0.88, 0.95), edgecolor=(0.45, 0.5, 0.6), lw=0.5, alpha=0.6)
    sc = 0.6
    ax.quiver(c[sl, 0], c[sl, 1], n[sl, 0], n[sl, 1], color="0.35", scale=1 / sc, scale_units="xy", width=0.004, label="face's own normal")
    ax.quiver(c[sl, 0], c[sl, 1], ns[sl, 0], ns[sl, 1], color="crimson", scale=1 / sc, scale_units="xy", width=0.004, label="smoothed normal")
    tilt = np.degrees(np.arccos(np.clip(-n[front, 0], -1, 1))); tilt_s = np.degrees(np.arccos(np.clip(-ns[front, 0], -1, 1)))
    stats.append((title, int(front.sum()), np.median(tilt), np.percentile(tilt, 90), np.median(tilt_s), np.percentile(tilt_s, 90), float((n[front, 0] > -0.1).mean()), b.mesh.n_active, b.mesh.n_elements))
    ax.set_title("{}\n{} active of {} elements at {:.1f} s".format(title, b.mesh.n_active, b.mesh.n_elements, T_STOP), fontsize=10)
    ax.set_xlabel("x, depth from the heated face (mm); heat arrives from the left"); ax.set_ylabel("y (mm)")
    ax.set_xlim(-0.3, 5.5); ax.set_ylim(-0.3, 10.3); ax.set_aspect("equal")
axes[0].legend(loc="lower right", fontsize=8)
fig.suptitle("Surface faces of the receded front (slice 4 < z < 6 mm) with their own and their smoothed normals, box at 6 MW/m2, t = 1.1 s", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "smoothed_normals_box.png"), dpi=140)
for row in stats:
    print("{:40s} front faces {:4d}  own-normal tilt from -x: median {:5.1f} deg, 90th pct {:5.1f}; smoothed: median {:5.1f}, 90th {:5.1f}; facing away {:.0%}; active {}/{}".format(*row))
