"""Schematic (not to scale) of the redesigned Step 3 model on a heavily molten sphere: skins on every corner holding the
whole non-rigid layer, a velocity in every skin cell, a moving 3D mesh whose surface is the rigid front.
(a) the whole fragment in a cut through the flight direction; (b) five neighbouring columns on the windward face;
(c) the proposed order of one macro step."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path

LIQ, SLU, SOLSKIN, CORE = "#f0a35e", "#c98b4f", "#cfe0f5", "#efe9dc"
BLUE, WARN, PURP, GREY = "#0b3d91", "#c0392b", "#7a4fb5", "#666666"


def arrow(ax, p, q, color, lw=1.8, ms=12, style="->", z=9):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=ms, color=color, lw=lw, zorder=z))


def note(ax, x, y, s, color="#333333", fs=9.6, **kw):
    ax.text(x, y, s, fontsize=fs, color=color, zorder=12, **kw)


fig = plt.figure(figsize=(20, 15))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.33], width_ratios=[1.0, 1.15], hspace=0.04, wspace=0.06)

# =========================================================================================== (a) whole fragment
ax = fig.add_subplot(gs[0, 0])
t = np.linspace(-np.pi / 2, np.pi / 2, 200)                    # front (windward, facing -x): flattened superellipse
p = 3.2
front = np.column_stack([-0.72 * np.abs(np.cos(t)) ** (2 / p), np.sin(t)])
tb = np.linspace(np.pi / 2, 3 * np.pi / 2, 200)                # back: a rounder half
back = np.column_stack([0.95 * np.abs(np.cos(tb)), np.sin(tb)])[::-1]
outline = np.vstack([front, back[::-1]])
# order the outline counter-clockwise from the bottom of the front
pts = np.vstack([front, np.column_stack([0.95 * np.cos(np.linspace(np.pi / 2, -np.pi / 2, 200)),
                                         np.sin(np.linspace(np.pi / 2, -np.pi / 2, 200))])])
# outward normals by finite differences
d = np.gradient(pts, axis=0)
nrm = np.column_stack([d[:, 1], -d[:, 0]]); nrm /= np.linalg.norm(nrm, axis=1)[:, None]
if (nrm[len(front) // 2] @ np.array([-1, 0])) < 0:
    nrm = -nrm
ang = np.degrees(np.arctan2(pts[:, 1], -pts[:, 0] + 0.15))      # angle from the nose, roughly
windward = np.arange(len(pts)) < len(front)
# layer thickness along the outline: thick non-rigid layer on the front, a pool at the nose, runoff thickening at the
# shoulders, a thin solid skin on the back
s = np.linspace(-1, 1, len(front))
h_nr = 0.10 + 0.07 * np.exp(-(s / 0.22) ** 2) + 0.05 * np.exp(-((np.abs(s) - 0.9) / 0.08) ** 2)
h_liq = 0.45 * h_nr
H = np.concatenate([h_nr, np.full(len(pts) - len(front), 0.035)])
Hl = np.concatenate([h_liq, np.zeros(len(pts) - len(front))])
inner = pts - nrm * H[:, None]
mid = pts - nrm * Hl[:, None]
# the rigid core: a triangulated mesh inside `inner`, compressed rings near the surface (the deforming zone)
rings = []
for k, f in enumerate([0.0, 0.035, 0.08, 0.14, 0.22]):
    ring = inner - nrm * f
    rings.append(ring[::6])
centre = np.array([[0.1, 0.0]])
rng = np.random.default_rng(3)
interior = rng.uniform(-0.5, 0.7, (260, 2))
inner_path = Path(np.vstack([inner, inner[:1]]))
ring4 = Path(np.vstack([rings[-1], rings[-1][:1]]))
interior = interior[ring4.contains_points(interior)]
P = np.vstack(rings + [interior, centre])
tri = mtri.Triangulation(P[:, 0], P[:, 1])
cen = P[tri.triangles].mean(axis=1)
tri.set_mask(~inner_path.contains_points(cen))
ax.add_patch(Polygon(inner, closed=True, fc=CORE, ec="none", zorder=1))
ax.triplot(tri, color="#a89f8c", lw=0.5, zorder=2)
zone = Polygon(np.vstack([inner[windward], (inner - nrm * 0.22)[windward][::-1]]), closed=True, fc="#f6d5cf",
               ec="none", alpha=0.55, zorder=2)
ax.add_patch(zone)
# skins: slurry band and liquid band on the front, solid skin on the back
ax.add_patch(Polygon(np.vstack([mid[windward], inner[windward][::-1]]), closed=True, fc=SLU, ec="none", zorder=3))
ax.add_patch(Polygon(np.vstack([pts[windward], mid[windward][::-1]]), closed=True, fc=LIQ, ec="none", zorder=3))
lee = ~windward
ax.add_patch(Polygon(np.vstack([pts[lee], inner[lee][::-1]]), closed=True, fc=SOLSKIN, ec="none", zorder=3))
ax.plot(pts[:, 0], pts[:, 1], color=BLUE, lw=2.2, zorder=4)
ax.plot(inner[:, 0], inner[:, 1], color="#333333", lw=1.0, zorder=4)
for i in range(0, len(pts), 14):                                  # corners: columns along the normal
    ax.plot([pts[i, 0], inner[i, 0]], [pts[i, 1], inner[i, 1]], color=BLUE, lw=0.6, alpha=0.7, zorder=5)
    ax.plot(*inner[i], "o", color=BLUE, ms=3.5, zorder=6)
# freestream, runoff, spray
for y in (-0.6, 0.0, 0.6):
    arrow(ax, (-1.85, y), (-1.35, y), "#555555", lw=2.2, ms=16)
note(ax, -1.95, 0.85, "gas flow", color="#555555", fs=11)
for sgn in (1, -1):
    arrow(ax, (-0.80, sgn * 0.62), (-0.80, sgn * 0.22), PURP, lw=2.4, ms=15)       # toward the nose (deceleration)
    arrow(ax, (-0.62, sgn * 0.86), (-0.22, sgn * 1.08), PURP, lw=2.4, ms=15)       # off the shoulder (pressure drop)
    arrow(ax, (0.2, sgn * 1.1), (0.75, sgn * 0.78), PURP, lw=1.6, ms=12)           # onto the back, where it freezes
note(ax, -1.62, 0.36, "pool at the nose:\ndeceleration pulls\nliquid inward", color=PURP)
note(ax, -0.25, 1.36, "runoff over the shoulder,\nwhere the pressure falls fast", color=PURP)
note(ax, 0.75, 0.93, "reaches the back's\nsolid skins and freezes:\n3D surface moves out", color=PURP)
rng = np.random.default_rng(7)
for sgn in (1, -1):
    for _ in range(26):
        u = rng.uniform(0.0, 1.0)
        x = -0.75 + 1.6 * u; y = sgn * (0.55 + 0.55 * u ** 0.6 + rng.uniform(0.02, 0.14))
        ax.add_patch(plt.Circle((x, y), rng.uniform(0.007, 0.02), color="#d35400", alpha=0.8, zorder=7))
    for _ in range(8):                                                       # just leaving the front
        y = sgn * rng.uniform(0.05, 0.5); x = -0.72 * np.abs(np.cos(np.arcsin(min(abs(y), 0.99)))) ** (2 / p) - rng.uniform(0.03, 0.08)
        ax.add_patch(plt.Circle((x, y), rng.uniform(0.007, 0.016), color="#d35400", alpha=0.8, zorder=7))
note(ax, 0.62, -1.5, "droplets leave the front and are swept\ndownstream around the shoulders", color="#d35400")
note(ax, -0.28, 0.08, "rigid core: moving 3D\nmesh (pink: deforming\nzone under the front)", color=GREY, bbox=dict(fc="white", ec="none", alpha=0.85))
ax.set_xlim(-2.0, 1.3); ax.set_ylim(-1.6, 1.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(a) The fragment late in flight (cut along the flight direction)", fontsize=13.5, loc="left")
# legend swatches
for k, (c, lab) in enumerate([(LIQ, "liquid (in the skin)"), (SLU, "slurry, more than half liquid (in the skin)"),
                              (SOLSKIN, "solid skin (back)"), (CORE, "rigid metal (3D model)")]):
    ax.add_patch(Rectangle((-1.95, -1.05 - 0.09 * k), 0.08, 0.06, fc=c, ec="#555555", lw=0.6, zorder=8))
    note(ax, -1.83, -1.05 - 0.09 * k, lab, fs=9.2)

# =========================================================================================== (b) five columns
ax = fig.add_subplot(gs[0, 1])
W = 2.0
cols = [  # (top of skin, liquid depth, slurry depth below liquid) ; base = top - liquid - slurry
    (2.2, 0.0, 0.0),     # 0: drained: only its minimum of solid metal (blue), base drawn inward
    (3.0, 0.9, 1.0),
    (3.4, 1.1, 1.3),
    (3.9, 1.3, 1.5),     # 3: pile-up: thicker, corner moved out
    (3.3, 0.9, 1.1),     # 4: bottom freezing: base moved outward
]
MINSKIN = 0.55
bases = []
for i, (top, hl, hs) in enumerate(cols):
    x0 = i * W
    if hl == 0.0:
        base = top - MINSKIN
        ax.add_patch(Rectangle((x0 + 0.05, base), W - 0.1, MINSKIN, fc=SOLSKIN, ec=BLUE, lw=1.0, zorder=3))
    else:
        base = top - hl - hs
        ax.add_patch(Rectangle((x0 + 0.05, base), W - 0.1, hs, fc=SLU, ec=BLUE, lw=1.0, zorder=3))
        ax.add_patch(Rectangle((x0 + 0.05, top - hl), W - 0.1, hl, fc=LIQ, ec=BLUE, lw=1.0, zorder=3))
    bases.append(base)
    # graded cells: fine at the top, coarser downward
    z, dz = top, 0.06
    while z - dz > base + 1e-9:
        z -= dz
        ax.plot([x0 + 0.05, x0 + W - 0.05], [z, z], color=BLUE, lw=0.35, alpha=0.6, zorder=4)
        if top - z > 0.45:
            dz *= 1.18
    ax.plot([x0 + 0.05, x0 + W - 0.05], [top, top], color=BLUE, lw=2.6, zorder=5)
    ax.plot(x0 + W / 2, base, "o", color=BLUE, ms=8, zorder=8)
    # velocity profile: zero at the base, slow in the slurry, fast in the liquid, fastest at the top
    if hl > 0:
        for zz in np.linspace(base + 0.08, top - 0.05, 9):
            if zz < top - hl:
                u = 0.12 * (zz - base) / max(hs, 1e-9)
            else:
                u = 0.12 + 0.75 * ((zz - (top - hl)) / hl) ** 0.8
            arrow(ax, (x0 + 0.35, zz), (x0 + 0.35 + u, zz), "#2c7a2c", lw=1.2, ms=8)
# lateral exchange between neighbouring columns at each depth
for i in (1, 2, 3):
    top, hl, hs = cols[i]
    for zz, w in [(top - 0.15, 2.2), (top - 0.55 * hl, 1.6), (top - hl - 0.5 * hs, 0.9)]:
        arrow(ax, ((i + 1) * W - 0.25, zz), ((i + 1) * W + 0.25, zz), PURP, lw=w, ms=10)
# the 3D model beneath: compressed elements of the deforming zone, surface nodes = bases
xs = np.array([W / 2 + i * W for i in range(5)])
b = np.array(bases)
lvl1 = b - 0.9; lvl2 = b - 1.9
ext = lambda arr: np.concatenate([[arr[0]], arr, [arr[-1]]])
xe = np.concatenate([[0.0], xs, [5 * W]])
for (ya, yb, off) in [(ext(b), ext(lvl1), 0), (ext(lvl1), ext(lvl2), 1)]:
    ax.fill_between(xe, yb, ya, color="#f6d5cf" if off == 0 else CORE, zorder=1)
    for k in range(len(xe) - 1):
        ax.plot([xe[k], xe[k + 1]], [ya[k], yb[k + 1]], color="#a89f8c", lw=0.7, zorder=2)
        ax.plot([xe[k], xe[k + 1]], [yb[k], ya[k + 1]], color="#a89f8c", lw=0.7, zorder=2)
    ax.plot(xe, yb, color="#a89f8c", lw=0.7, zorder=2)
ax.fill_between(xe, ext(lvl2) - 0.7, ext(lvl2), color=CORE, zorder=1)
ax.plot(xe, ext(b), color="#333333", lw=1.4, zorder=6)
# spray from the tops
rng = np.random.default_rng(11)
for i in (1, 2, 3, 4):
    top = cols[i][0]
    for _ in range(5):
        x = i * W + rng.uniform(0.4, 1.8); y = top + rng.uniform(0.25, 0.85)
        ax.add_patch(plt.Circle((x, y), rng.uniform(0.04, 0.09), color="#d35400", alpha=0.85, zorder=7))
arrow(ax, (2 * W + 0.9, cols[2][0] + 0.05), (2 * W + 1.25, cols[2][0] + 0.55), "#d35400", lw=1.6)
arrow(ax, (1.9, 3.95), (2.6, 3.65), "#d35400", lw=1.2)
note(ax, 0.1, 4.05, "spray from the top cells\n(Girin: b = the whole non-rigid depth here)", color="#d35400")
note(ax, 0.1, cols[0][0] + 0.25, "drained column: drew\nrigid metal back to its\nminimum; its node moved in", color=BLUE)
arrow(ax, (W / 2, bases[0] - 0.65), (W / 2, bases[0] - 0.1), BLUE, lw=1.8)
note(ax, 3 * W + 0.1, cols[3][0] + 0.25, "pile-up: column\nthicker, corner\nmoved outward", color=PURP, bbox=dict(fc="white", ec="none", alpha=0.85))
note(ax, 4 * W + 0.15, bases[4] - 1.15, "bottom freezing:\nnode moves out", color=BLUE)
arrow(ax, (4 * W + W / 2, bases[4] - 0.55), (4 * W + W / 2, bases[4] - 0.05), BLUE, lw=1.8)
ax.plot([4 * W + 0.1, 5 * W - 0.1], [bases[4] - 0.25, bases[4] - 0.25], color=BLUE, lw=1, ls=":", zorder=6)
note(ax, 1 * W + 0.1, bases[1] - 2.3, "3D model: surface nodes (blue dots) are the columns' bases = the rigid front;\n"
     "pink = deforming zone, compressed or stretched as the nodes move; rebuilt when elements degrade",
     color=GREY)
note(ax, 1 * W + 1.05, cols[1][0] - 0.3, "speed in\nevery cell", color="#2c7a2c", fs=9.2, bbox=dict(fc="white", ec="none", alpha=0.85))
note(ax, 1.55 * W, cols[2][0] - cols[2][1] - cols[2][2] * 0.25, "", color=GREY)
note(ax, 2 * W + 0.15, cols[2][0] - cols[2][1] - 0.18, "slurry: viscous,\nmoves slowly", color="#5b3a1a", fs=9.2, bbox=dict(fc="white", ec="none", alpha=0.8))
note(ax, 2 * W + 0.15, cols[2][0] - 0.38, "liquid: fast", color="#5b3a1a", fs=9.2, bbox=dict(fc="white", ec="none", alpha=0.8))
note(ax, 3 * W + 0.3, cols[3][0] - cols[3][1] - cols[3][2] - 0.35, "", color=GREY)
note(ax, 0.1, 4.95, "Runoff direction in this cut  →   (purple arrows: mass, momentum and heat crossing at each "
     "depth; the thicker the arrow, the faster that layer)", color=PURP, fs=10)
ax.set_xlim(-0.1, 5 * W + 0.1); ax.set_ylim(-1.2, 5.2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(b) Five neighbouring columns on the front (depths exaggerated; fine 10 um cells at the top, "
             "graded below)", fontsize=13.5, loc="left")

# =========================================================================================== (c) the macro step
ax = fig.add_subplot(gs[1, :])
steps = [("1. Gas flow at the corners",
          "smoothed normals give each corner\nits heating, pressure, shear and\nthe deceleration (both components)"),
         ("2. Skin sub-steps (10 ms)",
          "heat: through the depth (implicit) and\nsideways (explicit); speed in every cell:\nviscous coupling through the depth;\n"
          "mass, momentum, heat cross to\nneighbours by layer; spray leaves\nthe top cells"),
         ("3. The base moves",
          "a column below its minimum draws\nrigid metal (node in); metal at its\nbottom turned rigid joins the 3D model\n"
          "(node out); same rule for metal\nturning non-rigid at the base"),
         ("4. 3D heat on the moving mesh",
          "nodes moved through the metal\n(moving-mesh term); heat exchanged\nwith each column at its base node"),
         ("5. Quality check",
          "every element's volume and shape;\nrebuild inside the same surface,\ntemperature carried over exactly,\n"
          "before any element degrades")]
for k, (head, body) in enumerate(steps):
    x = 0.012 + k * 0.198
    ax.add_patch(FancyBboxPatch((x, 0.05), 0.176, 0.82, boxstyle="round,pad=0.008", fc="#f4f6fb", ec=BLUE, lw=1.3,
                                transform=ax.transAxes))
    ax.text(x + 0.008, 0.80, head, transform=ax.transAxes, fontsize=10.6, weight="bold", color=BLUE, va="top")
    ax.text(x + 0.008, 0.63, body, transform=ax.transAxes, fontsize=9.4, va="top")
    if k < 4:
        ax.annotate("", xy=(x + 0.193, 0.46), xytext=(x + 0.181, 0.46), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", lw=2, color=BLUE))
ax.axis("off")
ax.set_title("(c) One macro step, in the order proposed (to be fixed in the spec)", fontsize=13.5, loc="left")
fig.savefig("molten_sphere_schematic.png", dpi=120, bbox_inches="tight")
