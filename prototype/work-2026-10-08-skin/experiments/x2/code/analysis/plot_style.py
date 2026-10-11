"""Shared style for the sweep-analysis plotting scripts in this directory.

Sequential "blue" ramp and light-mode chart palette from the dataviz skill's
reference palette (references/palette.md: validated single-hue sequential
ramp, steps 100->700), reused verbatim across every plot in analysis/ so
they read as one family.
"""
import os

from matplotlib.colors import LinearSegmentedColormap

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_SUMMARY = os.path.join(REPO_ROOT, "sphere_sweep_output", "sweep_summary_AA7075.csv")
DEFAULT_PLOT_DIR = os.path.join(REPO_ROOT, "sphere_sweep_output", "plots")

# --- validated sequential "blue" ramp, dataviz skill references/palette.md -------
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQ_BLUE = LinearSegmentedColormap.from_list("seq_blue", BLUE_RAMP, N=256)

# --- diverging blue <-> red through the palette's neutral midpoint (#f0efec): a cool->warm
#     thermal map -- blue = cooler, red = hotter -- with equal steps per arm ----------------
COOL_WARM_STEPS = ["#0d366b", "#256abf", "#6da7ec", "#cde2fb",   # blue arm (cool)
                   "#f0efec",                                     # neutral midpoint
                   "#edc6c3", "#e99492", "#e34948", "#882c2b"]    # red arm (warm)
COOL_WARM = LinearSegmentedColormap.from_list("cool_warm", COOL_WARM_STEPS, N=256)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECOND = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"


def apply_rcparams(plt, font_size=10.5):
    """Apply the shared light-mode chart style to matplotlib's rcParams."""
    plt.rcParams.update({
        "font.size": font_size,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "text.color": INK,
        "axes.labelcolor": SECOND,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.edgecolor": AXIS,
    })


def strip_top_right_spines(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def style_colorbar(cbar):
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(color=MUTED, labelcolor=MUTED)


def material_from_run_name(run_name, default="drama-AA7075"):
    """Material name a sweep summary was produced with, taken from the run-name suffix.

    sphere_reentry.run_name appends ``_m<material>`` only for non-default materials
    (``..._h077.500km_mdrama-TiAl6v4``); a run name with no suffix is the default AA7075.
    """
    tail = run_name.rsplit("km", 1)[-1]
    return tail[2:] if tail.startswith("_m") else default
