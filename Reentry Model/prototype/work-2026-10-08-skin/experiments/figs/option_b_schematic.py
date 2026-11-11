"""Schematic (not to scale) of option (b): the skin as its own surface, moved by its own absorbed mass, over coarse
reservoir elements that die when the skin has absorbed all their metal."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.path import Path

H, S = 2.0, 0.6                     # element size and skin thickness in drawing units (the skin drawn thicker than 400 um)
X0, X1, Y1 = -8.0, 8.0, 9.0
def recession(x, R):                 # depth of the surface (skin top) below the original surface
    return R * (0.35 + 0.65 * np.exp(-(x / 5.5) ** 2))
# coarse reservoir elements: squares of side H split into two triangles (alternating diagonals)
tris = []
for i, x in enumerate(np.arange(X0, X1, H)):
    for j, y in enumerate(np.arange(0.0, Y1, H)):
        a, b, c, d = (x, y), (x + H, y), (x + H, y + H), (x, y + H)
        tris += [(a, b, c), (a, c, d)] if (i + j) % 2 == 0 else [(a, b, d), (b, c, d)]
tris = np.array(tris)
uv = np.array([(u, v) for u in np.linspace(0, 1, 21) for v in np.linspace(0, 1, 21) if u + v <= 1])
def absorbed_fraction(tri, base):
    p = tri[0] + uv[:, :1] * (tri[1] - tri[0]) + uv[:, 1:] * (tri[2] - tri[0])
    return float(np.mean(p[:, 1] < base(p[:, 0])))

panels = [(0.0, 0.0, "1. Skin created at the solidus",
           "It takes the outer 400 um of metal\nfrom the elements beneath it.\nNothing has melted; the surface\nhas not moved."),
          (1.6, 0.35, "2. Melting and spraying",
           "The top of the skin is liquid: that is\nthe film. Spray carries liquid away;\nthe skin draws the same mass at its\nbase, so the whole skin moves inward."),
          (3.6, 0.5, "3. Elements absorbed",
           "Elements whose metal the skin has\nfully absorbed die (dotted ghosts).\nBookkeeping only: the skin surface\ndoes not jump, no new skins, no tilt."),
          (5.8, 0.5, "4. Continued recession",
           "The skin surface stays smooth and\nfollows the mass. Live elements below\n(cut by the skin base here; each is\nheld whole in the 3D model with the\nmetal it has left) are the reservoir.")]
fig = plt.figure(figsize=(18, 10))
gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 0.34, 1.05], hspace=0.05, wspace=0.08)
xs = np.linspace(X0, X1, 400)
for k, (R, liq, title, text) in enumerate(panels):
    ax = fig.add_subplot(gs[0, k])
    top = lambda x, R=R: recession(x, R)
    base = lambda x, R=R: recession(x, R) + S
    from matplotlib.patches import Rectangle
    clip_poly = np.vstack([np.column_stack([xs, base(xs)]), [[X1, Y1 + 1], [X0, Y1 + 1]]])
    clip = Polygon(clip_poly, closed=True, transform=ax.transData, fc="none", ec="none")
    ax.add_patch(clip)
    for tri in tris:
        f = absorbed_fraction(tri, base)
        ax.add_patch(Polygon(tri, closed=True, fill=False, ls=":", lw=0.6, ec="0.65"))      # the original mesh, ghosted
        if f < 0.995:
            p = Polygon(tri, closed=True, fc=(0.80, 0.84, 0.90), ec="0.3", lw=0.8)
            ax.add_patch(p); p.set_clip_path(clip)
    # the skin band, coloured: liquid top (film) then solid
    yt, yb = top(xs), base(xs)
    yl = yt + liq * S
    ax.fill_between(xs, yt, yl, color="#f08c2e", alpha=0.95, lw=0)
    ax.fill_between(xs, yl, yb, color="#6c8ebf", alpha=0.9, lw=0)
    ax.plot(xs, yt, color="#b5651d" if liq else "#3d5f8f", lw=1.6)
    ax.plot(xs, yb, color="#2b4570", lw=1.0)
    # skin columns along the local normal
    for xc in np.arange(X0 + 0.5, X1, 1.0):
        dydx = (top(xc + 1e-3) - top(xc - 1e-3)) / 2e-3
        n = np.array([-dydx, 1.0]) / np.hypot(dydx, 1.0)
        p0 = np.array([xc, top(xc)])
        ax.plot([p0[0], p0[0] + S * n[0]], [p0[1], p0[1] + S * n[1]], color="white", lw=0.5, alpha=0.8)
    ax.plot(xs, 0 * xs, ls="--", color="0.5", lw=1.0)
    # heat in at the top
    for xh in (-6, -3, 0, 3, 6):
        q = 0.35 + 0.65 * np.exp(-(xh / 5.5) ** 2)
        ax.add_patch(FancyArrowPatch((xh, top(xh) - 1.4 - 0.2), (xh, top(xh) - 0.05), arrowstyle="-|>", mutation_scale=8 + 6 * q,
                                     color="crimson", lw=0.8 + 1.4 * q))
    if R > 0:
        rng = np.random.default_rng(3)
        for xd in rng.uniform(-4.5, 4.5, 14):     # spray
            ax.plot(xd, top(xd) - rng.uniform(0.3, 1.3), "o", ms=rng.uniform(2, 4.5), color="#1f77b4")
        for xg in (-5, -2, 1, 4):                 # draw at the base
            ax.add_patch(FancyArrowPatch((xg + 0.4, base(xg) + 0.75), (xg + 0.4, base(xg) + 0.05), arrowstyle="-|>", mutation_scale=9,
                                         color="#2ca02c", lw=1.4))
    ax.set_xlim(X0, X1); ax.set_ylim(Y1 - 1.0, -2.0); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=12, loc="left")
    tx = fig.add_subplot(gs[1, k]); tx.axis("off")
    tx.text(0.02, 1.0, text, va="top", ha="left", fontsize=10.5)

# legend row inside the text axes of panel 1
lax = fig.add_subplot(gs[2, :2]); lax.axis("off")
items = [("#f08c2e", "liquid top of the skin = the film"), ("#6c8ebf", "solid part of the skin (10 um cells; drawn thicker than 400 um)"),
         ((0.80, 0.84, 0.90), "coarse element: the metal not yet absorbed (the reservoir)"), ("white", "dotted: the original element mesh; empty = absorbed")]
for i, (c, lab) in enumerate(items):
    lax.add_patch(Polygon([(0.02, 0.95 - 0.105 * i), (0.07, 0.95 - 0.105 * i), (0.07, 0.89 - 0.105 * i), (0.02, 0.89 - 0.105 * i)],
                          closed=True, fc=c, ec="0.4", ls=":" if c == "white" else "-", transform=lax.transAxes))
    lax.text(0.09, 0.92 - 0.105 * i, lab, transform=lax.transAxes, va="center", fontsize=10.5)
others = [("0.5", "--", 1.0, "original surface"),
          ("crimson", "-", 2.0, "heating (from the smoothed skin surface's normals)"), ("#1f77b4", "o", 0, "spray: liquid leaving the skin top"),
          ("#2ca02c", "-", 2.0, "draw: metal entering the skin base, equal to the mass that left")]
for i, (c, ls, lw, lab) in enumerate(others):
    y = 0.95 - 0.105 * (i + 4)
    if ls == "o":
        lax.plot([0.045], [y - 0.03], "o", color=c, ms=6, transform=lax.transAxes)
    else:
        lax.plot([0.02, 0.07], [y - 0.03, y - 0.03], ls=ls, color=c, lw=lw, transform=lax.transAxes)
    lax.text(0.09, y - 0.03, lab, transform=lax.transAxes, va="center", fontsize=10.5)

# zoom: one skin column through one step
zax = fig.add_subplot(gs[2, 2:]); zax.set_xlim(0, 10); zax.set_ylim(10.5, -1.2); zax.axis("off")
zax.set_title("One skin column over one step (zoomed, not to scale)", fontsize=12, loc="left")
ncell, ch = 16, 0.42
base0 = 1.2 + ncell * ch
for j, (x0, lab, shift) in enumerate(((1.0, "start of step", 0.0), (5.6, "end of step", 1.0))):
    y0 = 1.2 + shift
    nliq = 4 if j == 0 else 5
    zax.add_patch(Polygon([(x0 - 0.3, base0), (x0 + 1.9, base0), (x0 + 1.9, 10.3), (x0 - 0.3, 10.3)], closed=True,
                          fc=(0.80, 0.84, 0.90), ec="0.35", lw=0.8))
    if shift:
        zax.add_patch(Polygon([(x0 - 0.3, base0), (x0 + 1.9, base0), (x0 + 1.9, base0 + shift), (x0 - 0.3, base0 + shift)], closed=True,
                              fc="white", ec="#2ca02c", hatch="///", lw=1.0))
    for c in range(ncell):
        col = "#f08c2e" if c < nliq else plt.get_cmap("Blues")(0.75 - 0.4 * c / ncell)
        zax.add_patch(Polygon([(x0, y0 + c * ch), (x0 + 1.6, y0 + c * ch), (x0 + 1.6, y0 + (c + 1) * ch), (x0, y0 + (c + 1) * ch)],
                              closed=True, fc=col, ec="white", lw=0.6))
    zax.text(x0 + 0.8, 10.0, "element (reservoir,\nfixed in space)", ha="center", va="bottom", fontsize=9)
    zax.text(x0 + 0.8, -0.9, lab, ha="center", fontsize=10.5, weight="bold")
    zax.plot([x0 - 0.5, x0 + 2.1], [1.2, 1.2], ls="--", color="0.5", lw=1)
zax.text(7.9, 1.2, "surface at the start of the step", fontsize=8.5, va="center", color="0.4")
zax.add_patch(FancyArrowPatch((1.8, -0.5), (1.8, 1.15), arrowstyle="-|>", mutation_scale=14, color="crimson", lw=2))
zax.text(2.1, 0.2, "heat in", color="crimson", fontsize=9.5)
zax.add_patch(FancyArrowPatch((6.4, 2.15), (6.4, 0.2), arrowstyle="-|>", mutation_scale=14, color="#1f77b4", lw=2))
zax.text(6.7, 0.75, "spray (and runoff) out:\nmass m leaves the top", color="#1f77b4", fontsize=9.5)
zax.add_patch(FancyArrowPatch((8.05, 8.4), (7.55, 8.4), arrowstyle="-|>", mutation_scale=14, color="#2ca02c", lw=2.2))
zax.text(8.15, 8.4, "draw: the same mass m\nleaves the element (hatched)\nand enters the skin base", color="#2ca02c", fontsize=9.5, va="center")
zax.annotate("", xy=(7.25, 2.2), xytext=(7.25, 1.2), arrowprops=dict(arrowstyle="<->", color="k"))
zax.text(7.4, 2.25, "recession = m / (rho A)", fontsize=9.5)
zax.text(2.9, 3.4, "melting inside the\nskin turns solid cells\nliquid: the film grows,\nthe skin's mass does\nnot change, nothing\nis drawn", fontsize=9.5)
zax.text(2.9, 6.6, "when the element\nhas nothing left it\ndies; the base moves\non into the element\nbelow", fontsize=9.5)
fig.suptitle("Option (b): the skin surface moves with its own mass over coarse reservoir elements (schematic, not to scale)", fontsize=14, y=0.93)
fig.savefig("option_b_schematic.png", dpi=110, bbox_inches="tight")
