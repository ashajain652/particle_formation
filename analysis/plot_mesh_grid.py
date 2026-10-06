"""Plots of the initial sphere grid: the surface triangulation, and a plane cross-section through it.

Two figures, written to reentry_model_output/mesh_plots/:

  mesh_surface.png       the boundary triangulation of the undisturbed 100 mm sphere -- what the heating, the melt
                         film and the spray all resolve the body with -- plus the distribution of patch edge lengths.
  mesh_cross_section.png a true z = 0 cut, each element coloured by its size, for the dense band (the melting
                         default) beside the Step 2 field, and the mean element size against depth for both.

The cross-section is the plane's actual intersection polygons (interpolated along the edges that straddle z = 0),
not a slab of centroids, so a cell's size on the page is its size in the mesh.

Style and palette from analysis/plot_style.py. Magnitude is sequential single-hue (SEQ_BLUE); the two-series
depth curve uses the palette's cool/warm pair, validated 2026-09-28 (protan dE 19.2, normal-vision dE 32.2,
contrast >= 3:1 against the chart surface).

    "$PY" analysis/plot_mesh_grid.py [--radius-mm 50] [--h-surface-mm 2] [--h-core-mm 8] [--band-mm 15]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.plot_style import INK, SECOND, MUTED, SEQ_BLUE, SURFACE, apply_rcparams
from reentry_model import mesh as meshmod

REGULAR_TET = 0.117851                      # V / a^3 for a regular tetrahedron of edge a
SERIES = {"band": "#256abf", "step2": "#e34948"}       # validated cool/warm pair, see the module docstring
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reentry_model_output", "mesh_plots")


def equivalent_edge(m):
    """Edge of the regular tetrahedron of the same volume, per element [m] -- one number for 'how big is this cell'."""
    return (m.element_volumes() / REGULAR_TET) ** (1.0 / 3.0)


def surface_edge_lengths(s):
    x = s.points[s.faces]
    return np.linalg.norm(x[:, [1, 2, 0]] - x, axis=2)


def plane_section(m, sizes):
    """Polygons where the elements cross z = 0, and each one's element size. The intersection of a tetrahedron with
    a plane is a triangle or a quadrilateral; its corners are interpolated along the edges whose ends straddle."""
    EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    z = m.points[m.tets][:, :, 2]
    crossing = np.flatnonzero((z.min(axis=1) < 0.0) & (z.max(axis=1) > 0.0))
    polys, values = [], []
    for e in crossing:
        nodes, zz = m.points[m.tets[e]], z[e]
        pts = []
        for i, j in EDGES:
            if (zz[i] < 0.0) != (zz[j] < 0.0):
                t = zz[i] / (zz[i] - zz[j])
                pts.append(nodes[i] + t * (nodes[j] - nodes[i]))
        if len(pts) < 3:
            continue
        p = np.array(pts)[:, :2]
        order = np.argsort(np.arctan2(*(p - p.mean(axis=0)).T[::-1]))    # wind the corners around the centroid
        polys.append(p[order] * 1e3)
        values.append(sizes[e] * 1e3)
    return polys, np.array(values)


def figure_surface(m, path):
    s = m.surface()
    edges = surface_edge_lengths(s)
    apply_rcparams(plt)
    fig = plt.figure(figsize=(12.5, 5.0))
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    view = np.array([0.62, -0.62, 0.48])                                 # matches elev/azim below
    facing = s.normals @ view > 0.12                                     # near-side facets only: no z-fighting
    ax.add_collection3d(Poly3DCollection(s.points[s.faces][facing] * 1e3, facecolors="#e8eef7",
                                         edgecolors="#6da7ec", linewidths=0.18))
    r = m.params["radius_m"] * 1e3
    for lim in (ax.set_xlim, ax.set_ylim, ax.set_zlim):
        lim(-0.62 * r, 0.62 * r)               # crop the empty margin mplot3d leaves: the sphere fills the panel
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=28, azim=-52)
    ax.set_axis_off()
    ax.set_title("{:,} patches on the {:.0f} mm sphere".format(s.n_patches, 2 * r), color=INK, pad=2)

    # a legible patch of it: the facets around the +x pole, projected on the plane they sit in
    cx = fig.add_subplot(1, 3, 2)
    theta = s.angles_to([1.0, 0.0, 0.0])
    near = theta < np.deg2rad(11.0)
    tri = s.points[s.faces][near][:, :, 1:] * 1e3                        # (y, z) -- the tangent plane at that pole
    cx.add_collection(PolyCollection(tri, facecolors="#e8eef7", edgecolors="#3987e5", linewidths=0.5))
    span = np.abs(tri).max()
    cx.set_xlim(-span, span), cx.set_ylim(-span, span)
    cx.set_aspect("equal")
    cx.set_xlabel("mm across the surface"), cx.set_ylabel("mm")
    cx.set_title("the stagnation cap, close up\nmean edge {:.2f} mm".format(edges.mean() * 1e3), color=INK, pad=6)

    hx = fig.add_subplot(1, 3, 3)
    hx.hist(edges.ravel() * 1e3, bins=48, color="#6da7ec", edgecolor=SURFACE, linewidth=0.4)
    hx.axvline(edges.mean() * 1e3, color=INK, lw=2.0)
    hx.annotate("mean {:.2f} mm".format(edges.mean() * 1e3), xy=(edges.mean() * 1e3, 0), xycoords=("data", "axes fraction"),
                xytext=(6, 0.90), textcoords=("offset points", "axes fraction"), color=INK, va="top")
    hx.annotate("longest {:.2f} mm".format(edges.max() * 1e3), xy=(0.97, 0.72), xycoords="axes fraction",
                ha="right", color=SECOND)
    hx.set_ylim(0, hx.get_ylim()[1] * 1.18)                                  # headroom: the mean line is tall
    hx.set_xlabel("patch edge length [mm]"), hx.set_ylabel("patch edges (3 per patch)")
    hx.set_title("how uniform the surface is\nh_surface asked for {:.2f} mm".format(m.params["h_surface_m"] * 1e3),
                 color=INK, pad=6)
    fig.suptitle("(1)  The initial sphere grid: its surface", color=INK, fontsize=13, y=1.02, x=0.02, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return s.n_patches, edges.mean() * 1e3, edges.max() * 1e3


def figure_cross_section(banded, step2, path):
    apply_rcparams(plt)
    fig = plt.figure(figsize=(12.5, 9.0))
    cuts = []
    for m in (banded, step2):
        sizes = equivalent_edge(m)
        cuts.append((m, *plane_section(m, sizes)))
    vmax = max(v.max() for _, _, v in cuts)
    vmin = min(v.min() for _, _, v in cuts)
    maps = []
    labels = ("dense band, {:.0f} mm".format(banded.params["band_m"] * 1e3), "band 0 -- the Step 2 field")
    for k, ((m, polys, values), label) in enumerate(zip(cuts, labels)):
        ax = fig.add_subplot(2, 2, k + 1)
        maps.append(ax)
        pc = PolyCollection(polys, array=values, cmap=SEQ_BLUE, edgecolors="#3332301f", linewidths=0.12)
        pc.set_clim(vmin, vmax)
        ax.add_collection(pc)
        r = m.params["radius_m"] * 1e3
        ax.set_xlim(-1.04 * r, 1.04 * r), ax.set_ylim(-1.04 * r, 1.04 * r)
        ax.set_aspect("equal")
        if k == 0:
            band = m.params["band_m"] * 1e3
            ax.add_artist(plt.Circle((0, 0), r - band, fill=False, ec=SERIES["band"], lw=2.0, ls=(0, (5, 3))))
            ax.annotate("inner edge of the band", xy=(0.0, -(r - band)), xytext=(0, -20), textcoords="offset points",
                        ha="center", va="top", color=SECOND)
        ax.set_xlabel("x [mm]"), ax.set_ylabel("y [mm]")
        ax.set_title("({})  {}\n{:,} elements, {:,} cut".format(
            "2a" if k == 0 else "2b", label, m.n_elements, len(polys)), color=INK, pad=6, fontsize=11)
        last = pc

    cb = fig.colorbar(last, ax=maps, location="right", fraction=0.028, pad=0.015, shrink=0.92)
    cb.set_label("element size [mm]", color=SECOND)                          # one bar: both maps share a clim
    cb.outline.set_visible(False)

    dx = fig.add_subplot(2, 1, 2)
    for (m, _, _), key, label in zip(cuts, ("band", "step2"), ("dense band 15 mm", "band 0 (Step 2 field)")):
        d, size = m.element_depths(), equivalent_edge(m)
        live = d >= 0.0
        bins = np.linspace(0.0, m.params["radius_m"], 26)
        idx = np.digitize(d[live], bins) - 1
        centres = 0.5 * (bins[:-1] + bins[1:]) * 1e3
        mean = np.array([size[live][idx == i].mean() * 1e3 if (idx == i).any() else np.nan for i in range(len(bins) - 1)])
        dx.plot(centres, mean, lw=2.0, color=SERIES[key], marker="o", ms=4.5, label=label)
        j = int(np.argmin(np.abs(centres - 9.0)))                            # label where the two are far apart,
        dx.annotate(label, xy=(centres[j], mean[j]),                          # not at the right edge where they meet
                    xytext=(0, -17 if key == "band" else 15), textcoords="offset points",
                    ha="center", color=SECOND, fontweight="medium")     # ink, not the series colour
    band = banded.params["band_m"] * 1e3
    dx.axvline(band, color=MUTED, lw=1.2, ls=(0, (4, 3)))
    dx.annotate("band ends, {:.0f} mm".format(band), xy=(band, 0.05), xycoords=("data", "axes fraction"),
                xytext=(6, 0), textcoords=("offset points", "axes fraction"), color=MUTED)
    dx.axhline(banded.params["h_core_m"] * 1e3, color=MUTED, lw=1.0, ls=":")
    dx.annotate("h_core {:.0f} mm".format(banded.params["h_core_m"] * 1e3), xy=(0.30, banded.params["h_core_m"] * 1e3),
                xycoords=("axes fraction", "data"), xytext=(0, 5), textcoords="offset points", ha="center", color=MUTED)
    dx.set_xlabel("depth below the surface [mm]"), dx.set_ylabel("mean element size [mm]")
    dx.set_title("(2c) what the band buys: cells stay at h_surface through its whole depth", color=INK, pad=6)
    dx.legend(frameon=False, loc="upper left")
    fig.suptitle("(2)  The initial sphere grid: a cross-section through it\n"
                 "The band is the melting default; band 0 is the Step 2 field, which every run in this package pins.",
                 color=INK, fontsize=12, y=1.02, x=0.02, ha="left", linespacing=1.6)
    fig.savefig(path, dpi=170, facecolor=SURFACE, bbox_inches="tight")   # no tight_layout: it fights the shared colorbar
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--radius-mm", type=float, default=50.0)
    ap.add_argument("--h-surface-mm", type=float, default=meshmod.DEFAULT_H_SURFACE * 1e3)
    ap.add_argument("--h-core-mm", type=float, default=meshmod.DEFAULT_H_CORE * 1e3)
    ap.add_argument("--band-mm", type=float, default=meshmod.DEFAULT_BAND * 1e3)
    a = ap.parse_args(argv)
    os.makedirs(OUT_DIR, exist_ok=True)
    R, hs, hc = a.radius_mm * 1e-3, a.h_surface_mm * 1e-3, a.h_core_mm * 1e-3
    banded = meshmod.sphere_mesh(R, hs, hc, band=a.band_mm * 1e-3)
    step2 = meshmod.sphere_mesh(R, hs, hc, band=0.0)
    print("banded : {:6d} nodes, {:7d} elements".format(banded.n_nodes, banded.n_elements))
    print("band 0 : {:6d} nodes, {:7d} elements".format(step2.n_nodes, step2.n_elements))
    p1 = os.path.join(OUT_DIR, "mesh_surface.png")
    n, mean_edge, max_edge = figure_surface(banded, p1)
    print("surface: {:d} patches, mean edge {:.2f} mm, longest {:.2f} mm -> {}".format(n, mean_edge, max_edge, p1))
    p2 = os.path.join(OUT_DIR, "mesh_cross_section.png")
    figure_cross_section(banded, step2, p2)
    print("cross-section -> {}".format(p2))


if __name__ == "__main__":
    main()
