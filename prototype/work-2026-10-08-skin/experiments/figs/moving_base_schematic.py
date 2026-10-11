"""Schematic (not to scale) of a skin column on a moving base, in a cut through three corner columns.
Row 1: the treatment decided so far (option i): the 3D mesh stays fixed; a partly absorbed element is held whole with
its remaining metal fraction; the base point is found inside an element and couples by interpolation.
Row 2: the alternative: the 3D mesh's surface nodes are the corners and move inward with the base, compressing the
elements beneath (a moving-mesh, or ALE, formulation), with a rebuild of the 3D mesh every few millimetres."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch

X0, X1, YB = -1.0, 7.0, -6.0
CORNERS = [1.0, 3.0, 5.0]
SKIN = 0.45                       # skin thickness in drawing units (elements are 2 units; the skin drawn ~1/4 of that)
LIQ, SOL, ELEM, DEAD = "#f0a35e", "#cfe0f5", "#efe9dc", "#bbbbbb"
NODE, WARN, BLUE = "#222222", "#c0392b", "#0b3d91"


def layers(y_rows):
    """Nodes of the cut: rows at depths y_rows; odd rows at x = -1, 1, 3, 5, 7, even rows at x = 0, 2, 4, 6."""
    rows = []
    for k, y in enumerate(y_rows):
        xs = np.arange(-1.0, 7.01, 2.0) if k % 2 == 0 else np.arange(0.0, 6.01, 2.0)
        rows.append(np.column_stack([xs, np.full(xs.size, y)]))
    return rows


def triangles(rows):
    tris = []
    for k in range(len(rows) - 1):
        a, b = rows[k], rows[k + 1]
        top_odd = a[:, 0][0] == -1.0
        if top_odd:          # a: -1,1,3,5,7 ; b: 0,2,4,6
            for i in range(len(b)):
                tris.append([a[i], a[i + 1], b[i]])
            for i in range(len(b) - 1):
                tris.append([b[i], b[i + 1], a[i + 1]])
        else:                # a: 0,2,4,6 ; b: -1,1,3,5,7
            for i in range(len(a)):
                tris.append([b[i], b[i + 1], a[i]])
            for i in range(len(a) - 1):
                tris.append([a[i], a[i + 1], b[i + 1]])
    return [np.array(t) for t in tris]


def inside(t, p):
    (x1, y1), (x2, y2), (x3, y3) = t
    d = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    l1 = ((y2 - y3) * (p[0] - x3) + (x3 - x2) * (p[1] - y3)) / d
    l2 = ((y3 - y1) * (p[0] - x3) + (x1 - x3) * (p[1] - y3)) / d
    lam = np.array([l1, l2, 1 - l1 - l2])
    return np.all(lam >= -1e-9), lam


def skin(ax, c, y_s, liquid_frac=0.0):
    x0, x1 = c - 0.97, c + 0.97
    y_b = y_s - SKIN
    ax.add_patch(Rectangle((x0, y_b), x1 - x0, SKIN, fc=SOL, ec=BLUE, lw=1.0, zorder=4))
    if liquid_frac > 0:
        ax.add_patch(Rectangle((x0, y_s - liquid_frac * SKIN), x1 - x0, liquid_frac * SKIN, fc=LIQ, ec="none", zorder=4))
    for yy in np.linspace(y_s, y_b, 7)[1:-1]:
        ax.plot([x0, x1], [yy, yy], color=BLUE, lw=0.4, alpha=0.7, zorder=5)
    ax.plot([x0, x1], [y_s, y_s], color=BLUE, lw=2.6, zorder=6)
    return y_b


def frame(ax, title):
    ax.plot([X0, X1], [0, 0], color="#888888", lw=1.0, ls=":", zorder=1)
    ax.set_xlim(X0 - 0.1, X1 + 0.1); ax.set_ylim(YB - 0.2, 0.9); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=12.5, loc="left")


def note(ax, x, y, s, color="#333333", **kw):
    ax.text(x, y, s, fontsize=9.3, color=color, zorder=9, **kw)


fig = plt.figure(figsize=(19, 14.5))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.08)
captions = []
rows0 = layers([0.0, -2.0, -4.0, -6.0])
tris0 = triangles(rows0)

# ===================================================================== row 1: fixed mesh, phi-scaled elements
stages1 = [(0.0, 0.0, "1. Start: the skin takes its\nminimum mass"),
           (-0.75, 0.45, "2. Spray and draws: the base moves\ninto an element"),
           (-1.95, 0.45, "3. That element fully absorbed:\nit dies")]
for k, (y_s, lf, title) in enumerate(stages1):
    ax = fig.add_subplot(gs[0, k])
    frame(ax, "(1-%d) %s" % (k + 1, title[3:]))
    y_b = y_s - SKIN
    base = np.array([3.0, y_b])
    host = None
    for t in tris0:
        ok, lam = inside(t, base)
        dead = k == 2 and t[:, 1].max() > -1.0 and t[:, 1].min() >= -2.0 - 1e-9 and 0.5 < t[:, 0].mean() < 5.5
        ax.add_patch(Polygon(t, closed=True, fc="none" if dead else ELEM, ec=DEAD if dead else "#888888",
                             lw=0.9, ls=":" if dead else "-", zorder=2))
        if ok and host is None:
            host, lam_host = t, lam
    for c in CORNERS:
        skin(ax, c, y_s, lf)
    for r in rows0:
        ax.plot(r[:, 0], r[:, 1], "o", color=NODE, ms=4, zorder=3)
    # base point and its interpolation to the host element's nodes
    for (px, py), w in zip(host, lam_host):
        ax.plot([base[0], px], [base[1], py], color=BLUE, lw=1.1, ls="--", zorder=7)
        if k > 0:
            note(ax, 0.5 * (base[0] + px) + 0.1, 0.5 * (base[1] + py), "%.2f" % w, color=BLUE)
    ax.plot(*base, "o", color=BLUE, ms=8, zorder=8)
    if k == 1:
        # the host element, its absorbed part, the swept region and the node above the surface
        ax.add_patch(Polygon(host, closed=True, fc="none", ec=WARN, lw=2.0, zorder=7))
        ax.add_patch(Rectangle((2.0, y_b), 2.0, -y_b, fc="none", ec="#7a4fb5", lw=1.6, ls="--", zorder=7))
        ax.plot(3.0, 0.0, "o", ms=13, mfc="none", mec=WARN, mew=2.2, zorder=9)
        note(ax, 3.25, 0.35, "node now outside the body:\nstill in the 3D model, gets\nheat from the base law", color=WARN)
        note(ax, 4.3, -2.75, "host element: debited\nfor all the metal drawn;\nheld whole with 60 %\nof its heat capacity", color=WARN)
        note(ax, -0.9, -2.75, "metal actually swept\n(dashed purple) spans\nthree elements", color="#7a4fb5")
        ax.add_patch(FancyArrowPatch((1.0, 0.75), (1.0, 0.2), arrowstyle="->", mutation_scale=12, color=LIQ, lw=2, zorder=9))
        note(ax, 1.1, 0.55, "spray", color="#b5651d")
        ax.add_patch(FancyArrowPatch((3.0, y_b - 0.75), (3.0, y_b - 0.12), arrowstyle="->", mutation_scale=12, color=BLUE,
                                     lw=2, zorder=9))
        note(ax, 3.1, y_b - 0.7, "draw", color=BLUE)
    if k == 0:
        note(ax, 3.6, -0.15, "base point: interpolated\nfrom 3 nodes (dashed)", color=BLUE, va="top")
        note(ax, -0.9, -3.0, "3D elements (2 mm; cut\nthrough tetrahedra)", color="#666666")
    if k == 2:
        note(ax, 3.6, y_b - 0.05, "base point now in the\nnext element down", color=BLUE, va="top")
        for px in (2.0, 4.0):
            ax.plot(px, -2.0, "o", ms=12, mfc="none", mec=WARN, mew=2.0, zorder=9)
        note(ax, 4.3, -3.4, "its top nodes again sit\ninside the skin", color=WARN)
        note(ax, -0.9, -1.0, "dead\n(dotted)", color="#888888")
captions.append((fig.axes[0], "Each corner's column (blue; liquid top orange) sits on the fixed 3D mesh. Its base point "
                 "is found inside one element and its temperature is read from that element's three nodes (four in 3D) "
                 "with the weights shown; the base law's heat goes back to those nodes with the same weights."))
captions.append((fig.axes[1], "Problems this panel shows: (a) the host element keeps its full outline, so one of its "
                 "nodes is outside the body and its upper part overlaps the skin, held only by scaling its capacity; "
                 "(b) all the drawn metal is debited from the one element holding the base point, though the column "
                 "actually swept metal from three; (c) the base law's heat partly goes to that outside node."))
captions.append((fig.axes[2], "Only when an element's last metal has been drawn does it die, and the base point moves into "
                 "the element below. Each draw also re-divides the column into equal cells (the drawn metal joins at "
                 "the bottom, the sprayed liquid leaves the top), which smears the profile slightly each time."))

# ===================================================================== row 2: moving mesh (ALE)
stages2 = [(0.0, 0.0, "1. Start: the 3D surface nodes\nsit at the base"),
           (-0.75, 0.45, "2. Spray and draws: the surface\nnodes move in"),
           (-3.1, 0.45, "3. Much later: elements thin,\nthe mesh is rebuilt")]
for k, (y_s, lf, title) in enumerate(stages2):
    ax = fig.add_subplot(gs[1, k])
    frame(ax, "(2-%d) %s" % (k + 1, title[3:]))
    y_b = y_s - SKIN
    # the node stack between the moving base and a fixed anchor at y = -6 is compressed uniformly
    ys = [y_b + y0 * (YB - y_b) / YB for y0 in (0.0, -2.0, -4.0, -6.0)]
    rows = layers(ys)
    tris = triangles(rows)
    thin = (ys[0] - ys[1]) < 1.3
    for t in tris:
        top = t[:, 1].max() > ys[1] + 1e-9
        ax.add_patch(Polygon(t, closed=True, fc="#f6d5cf" if (thin and top) else ELEM, ec="#888888", lw=0.9, zorder=2))
    for c in CORNERS:
        skin(ax, c, y_s, lf)
    for r in rows:
        ax.plot(r[:, 0], r[:, 1], "o", color=NODE, ms=4, zorder=3)
    for c in CORNERS:
        ax.plot(c, y_b, "o", color=BLUE, ms=9, zorder=8)
    if k == 0:
        note(ax, 3.6, y_b - 0.1, "base point = the 3D\nsurface node itself", color=BLUE, va="top")
    if k == 1:
        ax.add_patch(FancyArrowPatch((1.0, 0.75), (1.0, 0.2), arrowstyle="->", mutation_scale=12, color=LIQ, lw=2, zorder=9))
        note(ax, 1.1, 0.55, "spray", color="#b5651d")
        for c in CORNERS:
            ax.add_patch(FancyArrowPatch((c + 0.55, -SKIN), (c + 0.55, y_b), arrowstyle="->", mutation_scale=10,
                                         color=BLUE, lw=1.4, zorder=9))
        note(ax, 3.6, y_b - 0.15, "each surface node moves in\nby its own column's draw;\nno element holds space\nthe skin occupies", color=BLUE, va="top")
        note(ax, -0.9, -4.2, "elements below\ncompressed slightly", color="#666666")
    if k == 2:
        note(ax, 3.6, y_b + 0.15, "top elements thin (pink):\nrebuild the 3D mesh and carry\nthe temperature field over\nexactly in energy", color=WARN, va="bottom")
captions.append((fig.axes[3], "The 3D mesh's surface nodes are the corners' base points: the 3D model's outer surface is "
                 "always the skins' base. Base temperature is that node's temperature; the base law's heat is a load on "
                 "that node alone."))
captions.append((fig.axes[4], "A draw moves the column's node inward by the metal drawn; the 3D model loses exactly that "
                 "volume of metal and the enthalpy it held. No element is held whole, no node is outside the body, "
                 "and no metal is debited from the wrong element. The heat equation gains a term for nodes moving "
                 "through the metal."))
captions.append((fig.axes[5], "There are no element deaths. After a few millimetres of recession the top elements are too "
                 "thin, and the 3D mesh is rebuilt inside the current base surface, with the temperature field "
                 "transferred conservatively. Over the 100 mm flight's 77 mm of nose recession that is perhaps "
                 "ten rebuilds."))

import textwrap
fig.canvas.draw()
for row in (0, 1):
    axes = [a for a, _ in captions[3 * row:3 * row + 3]]
    y_cap = min(a.get_position().y0 for a in axes) + 0.01
    for a, txt in captions[3 * row:3 * row + 3]:
        fig.text(a.get_position().x0, y_cap, textwrap.fill(txt, 60), va="top", fontsize=10.3)
r1 = min(a.get_position().y1 for a in fig.axes[:3]); r2 = min(a.get_position().y1 for a in fig.axes[3:6])
fig.text(0.125, r1 + 0.05, "Row 1 - the treatment decided so far: fixed 3D mesh; heat couples at a base point inside a partly "
         "absorbed element (cut through three corner columns; depths exaggerated; not to scale)", fontsize=12, weight="bold")
fig.text(0.125, r2 + 0.05, "Row 2 - the alternative: the 3D mesh's surface nodes are the corners and follow the base",
         fontsize=12, weight="bold")
fig.savefig("moving_base_schematic.png", dpi=125, bbox_inches="tight")
