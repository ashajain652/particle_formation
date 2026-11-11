"""Schematic of option 3 (corners follow the columns by averaging; each column remembers its offset) and of what the
offset does over many steps, from a toy surface: a flat periodic patch of 24 x 24 squares, each split into two triangles,
columns on the triangles, every step each column recedes by its own fixed amount. Offsets are in units of one step's
recession of the slow columns. Compared: averaging without memory (option 1), averaging with memory (option 3), the best
least-squares corners with memory (option 2 plus memory), and columns placed on the corners instead of the triangles."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch

# ---------------------------------------------------------------- toy surface
n = 24
tri = []
for i in range(n):
    for j in range(n):
        a, b, c, d = i * n + j, ((i + 1) % n) * n + j, ((i + 1) % n) * n + (j + 1) % n, i * n + (j + 1) % n
        tri += [(a, b, c), (a, c, d)]
tri = np.array(tri)
T, V = len(tri), n * n
cent = np.array([[(t // 2) // n + (2 / 3 if t % 2 == 0 else 1 / 3), (t // 2) % n + (1 / 3 if t % 2 == 0 else 2 / 3)]
                 for t in range(T)]) - n / 2
A = sp.csr_matrix((np.full(3 * T, 1 / 3), (np.repeat(np.arange(T), 3), tri.ravel())), shape=(T, V))
W = sp.csr_matrix((np.ones(3 * T), (tri.ravel(), np.repeat(np.arange(T), 3))), shape=(V, T))
W = sp.diags(1.0 / np.asarray(W.sum(axis=1)).ravel()) @ W          # corner = mean of its six triangles


def run(r, rule, steps=400):
    off, hist = np.zeros(T), []
    for _ in range(steps):
        dem = r if rule == "avg" else r + off
        c = W @ dem if rule in ("avg", "avg_mem") else spl.lsqr(A, dem, atol=1e-12, btol=1e-12)[0]
        off = off + r - A @ c
        hist.append(np.abs(off).max())
    return np.array(hist)


patterns = [("A smooth bump\n(about five triangles wide, up to twice as fast)", 1 + np.exp(-(cent ** 2).sum(1) / (2 * 2.5 ** 2))),
            ("One triangle receding twice as fast\nas its neighbours", np.where(np.arange(T) == T // 2 + n, 2.0, 1.0)),
            ("A strip one triangle wide receding\ntwice as fast", np.where(np.abs(cent[:, 0]) < 0.5, 2.0, 1.0))]
curves = [[run(r, rule) for rule in ("avg", "avg_mem", "ls_mem")] for _, r in patterns]

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(18, 13))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.75], hspace=0.62, wspace=0.25)
captions = []
FAST, SLOW, EDGE = "#d9603b", "#f3d9b1", "#555555"

# (a) top view: columns on triangles share corners
ax = fig.add_subplot(gs[0, 0])
m = 4
for i in range(m):
    for j in range(m):
        p = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
        for t, (u, v, w) in enumerate([(0, 1, 2), (0, 2, 3)]):
            fast = (i, j, t) == (1, 1, 0)
            ax.add_patch(Polygon([p[u], p[v], p[w]], closed=True, fc=FAST if fast else SLOW, ec=EDGE, lw=1.0))
star = [(1, 1), (2, 1), (2, 2), (1, 2)]
for (x, y) in [(i, j) for i in range(m + 1) for j in range(m + 1)]:
    ax.plot(x, y, "o", ms=6, color="k", zorder=3)
for (x, y) in star[:3]:
    ax.plot(x, y, "o", ms=11, mfc="none", mec="#1f5fbf", mew=2.2, zorder=4)
ax.text(1.62, 1.33, "fast\ncolumn", ha="center", va="center", fontsize=10, color="white", weight="bold")
ax.set_xlim(-0.2, m + 0.2); ax.set_ylim(-0.2, m + 0.2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(a) Columns on the triangles (option A as chosen)", fontsize=13, loc="left")
captions.append((ax, "Seen from outside. Each triangle carries a column; each corner (dot) is shared by six triangles.\n"
        "The fast column's three corners (blue rings) are each shared with five slow columns, so\n"
        "they cannot sink as fast as it does. A closed surface has about two triangles per corner:\n"
        "the columns hold twice as many independent depths as the corners can express."))

# (b) cross-section: what the offset is
ax = fig.add_subplot(gs[0, 1])
S = 1.3                                            # skin thickness (drawing units; not to scale)
true = np.array([1.0, 2.6, 1.0])                   # true recession of three columns (from their absorbed mass)
xs = np.array([0.0, 2.0, 4.0, 6.0])
corner = np.array([true[0], true[:2].mean(), true[1:].mean(), true[2]])   # averaging rule, one corner between two columns
drawn = 0.5 * (corner[:-1] + corner[1:])
for gx in np.arange(-0.5, 6.6, 1.0):
    ax.plot([gx, gx], [-0.3, -5.6], color="#bbbbbb", lw=0.8, zorder=0)
for gy in np.arange(-0.3, -5.7, -1.0):
    ax.plot([-0.5, 6.5], [gy, gy], color="#bbbbbb", lw=0.8, zorder=0)
ax.text(6.55, -4.9, "3D elements\n(the reservoir)", fontsize=9.5, color="#777777", va="center")
ax.plot([-0.5, 6.5], [0, 0], color="#999999", lw=1.2, ls=":")
ax.text(-0.45, 0.12, "surface at the start", fontsize=9.5, color="#777777")
for k in range(3):
    x0, x1 = xs[k] + 0.08, xs[k + 1] - 0.08
    ax.add_patch(Rectangle((x0, -true[k] - S), x1 - x0, S, fc="#cfe0f5", ec="#1f5fbf", lw=1.0, alpha=0.9, zorder=1))
    for cy in np.linspace(-true[k], -true[k] - S, 8)[1:-1]:
        ax.plot([x0, x1], [cy, cy], color="#1f5fbf", lw=0.4, alpha=0.6, zorder=1)
    ax.plot([x0, x1], [-true[k], -true[k]], color="#0b3d91", lw=4, solid_capstyle="butt", zorder=3)
    xc = 0.5 * (xs[k] + xs[k + 1])
    ax.plot(xc, -true[k] - S, "o", ms=9, color="#0b3d91", zorder=4)
    if abs(drawn[k] - true[k]) > 0.05:
        ax.add_patch(FancyArrowPatch((xc + 0.55, -drawn[k]), (xc + 0.55, -true[k]), arrowstyle="<->", mutation_scale=12,
                                     color="#b8860b", lw=1.8, zorder=5))
        if k == 1:
            ax.text(xc + 0.68, -0.5 * (drawn[k] + true[k]), "offset", fontsize=9.5, color="#b8860b", va="center")
ax.plot(xs, -corner, "-", color=FAST, lw=2.4, zorder=4)
ax.plot(xs, -corner, "o", color=FAST, ms=7, zorder=5)
ax.text(3.0, -1.3, "drawn surface\n(corners averaged)", fontsize=9.5, color=FAST, ha="center", va="center")
ax.text(3.05, -true[1] - 0.18, "true surface\nof the column", fontsize=9.5, color="#0b3d91", va="top", ha="center")
ax.text(3.0, -true[1] - S - 0.25, "base point\n= true surface - skin", fontsize=9.5, color="#0b3d91", va="top", ha="center")
ax.set_xlim(-0.6, 7.6); ax.set_ylim(-5.7, 0.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(b) What the offset is (cut through three columns)", fontsize=13, loc="left")
captions.append((ax, "Each column's true surface comes from its absorbed mass (dark blue bar); its skin (light blue,\n"
        "40 cells of 10 um) hangs below it and its heat couples at its base point (dot), always at the\n"
        "right depth. The drawn surface, through shared corners, gives only the normals and the runoff\n"
        "graph. The offset is the gap between the two. Depths exaggerated; not to scale."))

# (c) top view: columns on the corners
ax = fig.add_subplot(gs[0, 2])
for i in range(m):
    for j in range(m):
        p = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
        for (u, v, w) in [(0, 1, 2), (0, 2, 3)]:
            ax.add_patch(Polygon([p[u], p[v], p[w]], closed=True, fc="none", ec="#aaaaaa", lw=0.8))


def dual_cell(cx, cy):
    # corners of the dual cell: centroids and edge midpoints of the six triangles around (cx, cy), in angular order
    nb = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]
    pts = []
    for k in range(6):
        a, b = nb[k], nb[(k + 1) % 6]
        pts.append((cx + a[0] / 2, cy + a[1] / 2))
        pts.append((cx + (a[0] + b[0]) / 3, cy + (a[1] + b[1]) / 3))
    return pts


for (cx, cy) in [(i, j) for i in range(1, m) for j in range(1, m)]:
    fast = (cx, cy) == (2, 2)
    ax.add_patch(Polygon(dual_cell(cx, cy), closed=True, fc=FAST if fast else SLOW, ec=EDGE, lw=1.0, alpha=0.95))
for (x, y) in [(i, j) for i in range(m + 1) for j in range(m + 1)]:
    ax.plot(x, y, "o", ms=6, color="k", zorder=3)
ax.text(2.0, 2.28, "fast\ncolumn", ha="center", va="bottom", fontsize=10, color="white", weight="bold")
ax.set_xlim(-0.2, m + 0.2); ax.set_ylim(-0.2, m + 0.2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(c) The alternative: one column per corner", fontsize=13, loc="left")
captions.append((ax, "Same triangles; each column owns the patch around one corner (a third of each touching\n"
        "triangle). One column, one corner: the corner moves by exactly its column's recession, so\n"
        "the drawn surface is the true surface and there is no offset. The enclosed volume equals the\n"
        "metal absorbed. Half as many columns (one per corner), at the same 2 mm spacing."))

# (d-f) offsets over 400 steps on the toy surface
labels = ["option 1: averaging, no memory", "option 3: averaging + memory", "best-fit corners + memory",
          "columns on corners (panel c)"]
colours = ["#999999", "#d9603b", "#7a4fb5", "#1f8f4e"]
steps = np.arange(1, 401)
for k, ((title, _), cs) in enumerate(zip(patterns, curves)):
    ax = fig.add_subplot(gs[1, k])
    for c, lab, col in zip(cs, labels[:3], colours[:3]):
        ax.plot(steps, c, color=col, lw=2.2, label=lab)
    ax.plot(steps, np.zeros_like(steps), color=colours[3], lw=2.6, label=labels[3])
    ax.set_title(title, fontsize=11.5, loc="left")
    ax.set_xlabel("steps")
    if k == 0:
        ax.set_ylabel("largest offset\n(in steps of the slow columns' recession)")
    ax.grid(alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    if k == 0:
        ax.legend(fontsize=9.5, frameon=False, loc="upper left")
fig.text(0.125, 0.035, "(d-f) Toy surface (flat, 24 x 24 squares split into triangles); each column recedes by its own fixed "
         "amount every step.\nOne step at 0.0125 s and 0.1-1 mm/s recession is about 1-12 um, so an offset of 100 steps "
         "is about 0.1-1.2 mm.", fontsize=10.5)
import textwrap
fig.canvas.draw()
y_cap = min(a.get_position().y0 for a, _ in captions) - 0.012
for a, txt in captions:
    body = " ".join(txt.split())
    fig.text(a.get_position().x0, y_cap, textwrap.fill(body, 64), va="top", fontsize=10.5)
fig.savefig("option3_offset_schematic.png", dpi=130, bbox_inches="tight")
