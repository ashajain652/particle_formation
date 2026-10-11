"""Own vs Taubin-smoothed normals on the receded box front: 3D view (faces coloured by the tilt of their own normal from
the heat direction, smoothed normals as arrows) and a true cross-section of the surface at z = 5 mm."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, "..", "x2", "code"))
from reentry_model import body, heating, material, mesh, skin, thermal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
T_STOP, DT, Q = 1.5, 0.05, 6.0e6
def run(msh):
    mat = material.Material.from_drama_json("AA7075_scheil")
    st = body.MeltSettings(removal="girin", size_feedback="initial", surface_model="skin", skin=skin.SkinSettings(cell=20e-6, substep=0.01))
    b = body.MeltingBody(msh, mat, thermal.thermal_solver("skfem", linear_solver="direct"), float(mat.rho * msh.volume()), settings=st, T0=700.0, emissivity=0.0)
    for k in range(int(round(T_STOP / DT))):
        s = b.surface; w = np.maximum(0.0, -s.normals[:, 0]); q = Q * 1e-4 * w / float((w * s.areas).sum())
        b.advance((k + 1) * DT, DT, heating.HeatingResult(q, Q, 0.0, 0.0, 0.0))
        L = b.m_f + b.m_d; b.removed_mass += float(L.sum()); b.removed_enthalpy += float((L * b.film_enthalpy()).sum())
        b.m_f[:] = 0; b.m_d[:] = 0; b.solver.set_film_mass(b._film_nodal())
    return b
cases = [("Kuhn cubes, 2 mm", mesh.box_mesh(8e-3, 10e-3, 10e-3, 2e-3)),
         ("prism layers 1.0 mm", mesh.layered_box_mesh(8e-3, 10e-3, 10e-3, 2e-3, 1.0e-3)),
         ("prism layers 0.5 mm", mesh.layered_box_mesh(8e-3, 10e-3, 10e-3, 2e-3, 0.5e-3))]
fig = plt.figure(figsize=(17, 11))
for i, (title, msh) in enumerate(cases):
    b = run(msh); s = b.surface; ns = s.smoothed_normals(); P = msh.points * 1e3
    c = s.centroids * 1e3; n = s.normals
    keep = (c[:, 0] < 7.95) & (np.abs(n[:, 1]) < 0.999 + 0 * c[:, 1]) | ((c[:, 0] < 7.95) & (c[:, 1] > 0.01) & (c[:, 1] < 9.99))
    keep &= ~((c[:, 0] > 7.9) & (n[:, 0] > 0.9)) & (c[:, 0] < 8.01) & (c[:, 1] > 0.01) & (c[:, 1] < 9.99) & (c[:, 2] > 0.01) & (c[:, 2] < 9.99)   # the front only
    tilt = np.degrees(np.arccos(np.clip(-n[:, 0], -1, 1)))
    ax = fig.add_subplot(2, 3, i + 1, projection="3d")
    polys = [P[f][:, [1, 2, 0]] for f in s.faces[keep]]          # plot (y, z, depth x)
    cmap = plt.get_cmap("viridis")
    pc = Poly3DCollection(polys, facecolors=cmap(np.clip(tilt[keep] / 180.0, 0, 1)), edgecolors=(0, 0, 0, 0.25), linewidths=0.3)
    ax.add_collection3d(pc)
    cs, nss = c[keep], ns[keep]
    ax.quiver(cs[:, 1], cs[:, 2], cs[:, 0], nss[:, 1], nss[:, 2], nss[:, 0], length=0.8, color="crimson", linewidth=0.6)
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_zlim(8.2, -0.5); ax.view_init(elev=28, azim=-60)
    ax.set_xlabel("y (mm)"); ax.set_ylabel("z (mm)"); ax.set_zlabel("depth x (mm)")
    ax.set_title("{}: front at 1.5 s\nface colour = own normal: dark purple = faces the heat, teal = sideways (90 deg), yellow = away\nred arrows = smoothed normals".format(title), fontsize=8)
    # cross-section z = 5 mm: intersect each front triangle with the plane
    ax2 = fig.add_subplot(2, 3, i + 4)
    z0 = 5.0
    for f in np.flatnonzero(keep):
        tri = P[s.faces[f]]; dz = tri[:, 2] - z0; pts = []
        for a, bb in ((0, 1), (1, 2), (2, 0)):
            if dz[a] * dz[bb] < 0:
                u = dz[a] / (dz[a] - dz[bb]); pts.append(tri[a] + u * (tri[bb] - tri[a]))
        if len(pts) == 2:
            pts = np.array(pts); ax2.plot(pts[:, 0], pts[:, 1], color="0.2", lw=1.6)
            m = pts.mean(axis=0)
            ax2.annotate("", xy=(m[0] + 0.6 * n[f, 0], m[1] + 0.6 * n[f, 1]), xytext=(m[0], m[1]), arrowprops=dict(arrowstyle="->", color="0.45", lw=1))
            ax2.annotate("", xy=(m[0] + 0.6 * ns[f, 0], m[1] + 0.6 * ns[f, 1]), xytext=(m[0], m[1]), arrowprops=dict(arrowstyle="->", color="crimson", lw=1.2))
    ax2.set_xlim(-1.0, 5.0); ax2.set_ylim(0, 10); ax2.set_aspect("equal")
    ax2.set_xlabel("depth x (mm); heat arrives from the left"); ax2.set_ylabel("y (mm)")
    ax2.set_title("cross-section of the surface at z = 5 mm (black)\ngrey = own normal, red = smoothed normal (in-plane parts)", fontsize=8)
    away = float((n[keep, 0] > -0.1).sum()) / keep.sum()
    print(f"{title:22s} front faces {keep.sum():4d}; own-normal tilt median {np.median(tilt[keep]):5.1f} deg; smoothed tilt median {np.median(np.degrees(np.arccos(np.clip(-ns[keep,0],-1,1)))):5.1f} deg; faces not facing the heat (own) {away:.0%}; depth range of the front {c[keep,0].min():.2f}-{c[keep,0].max():.2f} mm")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "smoothed_normals_box_3d_v2.png"), dpi=120)
