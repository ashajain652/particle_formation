"""What the plot scripts share: reading a folder of finished runs, choosing the baseline ones, and the plot style.

A sweep folder holds one sub-folder per run (summary.json, history.csv, size_distribution.csv). The same folder can
hold the baseline sweep, the NRLMSIS + AFGL cross-check, the sensitivity runs and shortened test runs, so the plot
scripts pick runs by their settings, not by their names: a run is a baseline run when every field of its Settings
that changes the physics (everything except the sampling, time-grid and solver fields in NUMERIC_FIELDS) equals the
default, apart from the ambient source, which each plot chooses. When several runs share a physics signature and a
case (for example a 12-hour run and a shortened one), the run with the longest end time wins.
"""
import dataclasses
import importlib.util
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402
import numpy as np                                                   # noqa: E402

from reentry_model import DRAMA_DIR                                  # noqa: E402
from reentry_model.wake_column.settings import Settings             # noqa: E402

# Settings fields that set how the run is sampled or integrated, not which physics it holds
NUMERIC_FIELDS = ("n_test_particles", "n_emission_samples", "seed", "t_end_s", "t_first_s", "points_per_decade",
                  "rtol")


# ---- reading runs ---------------------------------------------------------------------------------------------------

def physics_signature(settings):
    """A hashable key of the physics fields of a settings dictionary (as in summary.json, or dataclasses.asdict)."""
    settings = json.loads(json.dumps(settings))                      # tuples and lists compare equal after this
    return tuple((k, json.dumps(v, sort_keys=True)) for k, v in sorted(settings.items()) if k not in NUMERIC_FIELDS)


def is_baseline(summary, ambient=None):
    """True when the run's physics settings are the defaults (its ambient source being `ambient` when given)."""
    default = dataclasses.asdict(Settings())
    if ambient is not None:
        default["ambient"] = ambient
    elif "ambient" in summary["settings"]:
        default["ambient"] = summary["settings"]["ambient"]
    return physics_signature(summary["settings"]) == physics_signature(default)


def load_summaries(out):
    """The summary of every finished run under `out` that has an outcome: a list of dictionaries, each with its folder
    under "_dir". Runs that vaporised no aluminium ("skipped" summaries, which carry no outcome) and unreadable or
    truncated summaries are left out."""
    found = []
    for name in sorted(os.listdir(out)):
        path = os.path.join(out, name, "summary.json")
        if not os.path.isfile(path):
            continue
        try:
            with open(path) as f:
                s = json.load(f)
        except (ValueError, OSError):
            continue
        if not isinstance(s, dict) or s.get("skipped") or "condensed_fraction_end" not in s:
            continue
        s["_dir"] = os.path.join(out, name)
        found.append(s)
    return found


def pick(summaries, signature, ambient=None):
    """{(d0_um, h0_km, time_of_day): summary} over the summaries whose physics signature is `signature` (and whose
    ambient is `ambient` when given); the longest run wins when several share a case."""
    chosen = {}
    for s in summaries:
        if ambient is not None and s["settings"]["ambient"] != ambient:
            continue
        if physics_signature(s["settings"]) != signature:
            continue
        c = s["case"]
        key = (float(c["d0_um"]), float(c["h0_km"]), c["time_of_day"])
        if key not in chosen or s["settings"]["t_end_s"] > chosen[key]["settings"]["t_end_s"]:
            chosen[key] = s
    return chosen


def baseline_runs(summaries, ambient):
    """The baseline runs of one ambient source, as {(d0_um, h0_km, time_of_day): summary}."""
    default = dataclasses.asdict(Settings())
    default["ambient"] = ambient
    return pick(summaries, physics_signature(default), ambient)


def read_history(summary):
    """history.csv of a run as {column: numpy array}."""
    return read_table(os.path.join(summary["_dir"], "history.csv"))


def read_size_distribution(summary):
    """size_distribution.csv of a run as {column: numpy array}; empty cells (masked composition) become NaN."""
    return read_table(os.path.join(summary["_dir"], "size_distribution.csv"))


def read_table(path):
    with open(path) as f:
        header = f.readline().strip().split(",")
        rows = [line.rstrip("\n").split(",") for line in f if line.strip()]
    cols = {}
    for j, name in enumerate(header):
        cols[name] = np.array([float(r[j]) if r[j] != "" else np.nan for r in rows])
    return cols


# ---- style ----------------------------------------------------------------------------------------------------------

def _load_plot_style():
    """DRAMA Sphere Experiments/analysis/plot_style.py as a module, or None when it is not there."""
    path = os.path.join(DRAMA_DIR, "analysis", "plot_style.py")
    if not os.path.isfile(path):
        return None
    try:
        spec = importlib.util.spec_from_file_location("drama_plot_style", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    except Exception:                                                # a broken style file must not stop the plots
        return None


_STYLE = _load_plot_style()
INK = getattr(_STYLE, "INK", "#0b0b0b")
SECOND = getattr(_STYLE, "SECOND", "#52514e")
MUTED = getattr(_STYLE, "MUTED", "#898781")
GRID = getattr(_STYLE, "GRID", "#e1e0d9")
RAMP = getattr(_STYLE, "SEQ_BLUE", None) or matplotlib.colormaps["Blues"]
# categorical colours for up to eight series (Okabe-Ito order, colour-blind safe)
CAT = ["#0072b2", "#d55e00", "#009e73", "#e69f00", "#cc79a7", "#56b4e9", "#999999", "#000000"]
DAY_COLOR, NIGHT_COLOR = "#e69f00", "#0072b2"


def apply_style():
    """The shared chart style (plot_style when importable, plain matplotlib otherwise)."""
    if _STYLE is not None:
        _STYLE.apply_rcparams(plt)
    plt.rcParams.update({"axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                         "lines.linewidth": 1.8, "axes.titlesize": 10.5, "axes.titleweight": "bold",
                         "savefig.dpi": 140})


def ramp_colours(n):
    """n colours along the sequential ramp, light to dark but never too pale for a line on a light ground."""
    return [RAMP(x) for x in np.linspace(0.3, 1.0, n)]


def diameter_tag(d):
    """'01.0um' style text used in file names (two integer digits, as run names do)."""
    return f"d{d:04.1f}um"


def save(fig, out, name):
    plots = os.path.join(out, "plots")
    os.makedirs(plots, exist_ok=True)
    path = os.path.join(plots, name)
    fig.savefig(path)
    plt.close(fig)
    return path
