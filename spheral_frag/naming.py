"""Names of prepared inputs, Spheral runs and their analyses (spec §4.4; M1 plan, Task 1).

Names encode the whole configuration, so that runs share one output directory without overwriting each other and
studies resume by checking for a run's summary file (CLAUDE.md). The pattern is the finite-element model's (Step 3's
suffixes): the always-present parts in a fixed order, then `_<label>-<value>` for every setting that differs from its
default, in a fixed key order.

  prepared   prep_<fe run>_k<K0>-<K1>[_every<N>]                       spheral_output/prepare/<prepared>/
  run        <prepared>__<mode>_k<k0>-<k1>_<form>_dx<dx>mm_seed<n>[_<label>-<value>...]
  analysis   the run's own name                                          spheral_output/analyse/<run>/

The double underscore separates the prepared name from the run's own parts; neither a finite-element run name nor
any part written here contains one. Frame indices are zero-padded to five digits. Names never contain spaces or `/`.

Numpy-free: the standard library only."""
from __future__ import annotations

import math
import re

MODES = ("replay", "separation", "slurry", "detachment", "synthetic")
FORMS = ("rz", "3d")

# Bracket settings and their defaults. Later milestones append keys; never reorder or rename one, because the order
# is the order of the suffixes and a renamed key would rename every run that set it.
BRACKET_DEFAULTS: dict[str, object] = {
    "film_limit_mm": 2.0,        # §10: zone 3 below this layer thickness
    "base_pressure_pct": 3.0,    # §7.2 lee extension: base pressure, % of p_w_stag_Pa (Question 4)
    "separation_deg": 180.0,     # §7.2 lee extension: shear set to zero beyond this inclination (Question 4)
    "lee_shear": 0.5,            # §7.2 lee extension: shear as a fraction of the last windward bin's (Question 4)
    "min_particles": 30,         # §11.1: smallest group recorded as a fragment
    "link_h": 1.5,               # §11.1: dust within link_h smoothing lengths of an intact member attaches
}
# key -> (label in the name, unit appended to the value)
BRACKET_LABELS: dict[str, tuple[str, str]] = {
    "film_limit_mm": ("filmlimit", "mm"),
    "base_pressure_pct": ("basepressure", "pct"),
    "separation_deg": ("separation", "deg"),
    "lee_shear": ("leeshear", ""),
    "min_particles": ("minparticles", ""),
    "link_h": ("linkh", ""),
}
_KEY_OF_LABEL = {label: key for key, (label, _) in BRACKET_LABELS.items()}

_PREPARE_RE = re.compile(r"^prep_(?P<fe_run>.+)_k(?P<k0>\d{5,})-(?P<k1>\d{5,})(?:_every(?P<every>\d+))?$")
_FRAMES_RE = re.compile(r"^k(?P<k0>\d{5,})-(?P<k1>\d{5,})$")
_DX_RE = re.compile(r"^dx(?P<dx>.+)mm$")
_SEED_RE = re.compile(r"^seed(?P<seed>\d+)$")


def _check_part(text, what):
    if not isinstance(text, str) or not text:
        raise ValueError("{} must be a non-empty string, got {!r}".format(what, text))
    bad = [c for c in (" ", "/", "\\", "\t", "\n") if c in text]
    if bad or "__" in text:
        raise ValueError("{} {!r} contains {}".format(what, text, "'__'" if not bad else repr(bad[0])))


def _check_index(k, what):
    if isinstance(k, bool) or not isinstance(k, int) or k < 0:
        raise ValueError("{} must be a non-negative integer, got {!r}".format(what, k))


def _frames(k0, k1):
    _check_index(k0, "first frame")
    _check_index(k1, "last frame")
    if k1 < k0:
        raise ValueError("frame range {}:{} is reversed".format(k0, k1))
    return "k{:05d}-{:05d}".format(k0, k1)


def _format_number(value, integer):
    """Shortest exact text: integers as such, an integral float without '.0', otherwise repr (round-trips)."""
    if isinstance(value, bool):
        raise ValueError("a bracket value must be a number, got {!r}".format(value))
    if integer:
        if not isinstance(value, int) and not (isinstance(value, float) and value.is_integer()):
            raise ValueError("expected an integer, got {!r}".format(value))
        return str(int(value))
    v = float(value)
    if not math.isfinite(v):
        raise ValueError("a name cannot hold {!r}".format(value))
    return str(int(v)) if v.is_integer() and abs(v) < 1e15 else repr(v)


def _format_dx(dx_mm):
    v = float(dx_mm)
    if not (math.isfinite(v) and v > 0.0):
        raise ValueError("dx must be a positive number of millimetres, got {!r}".format(dx_mm))
    text = "{:.2f}".format(v)
    return text if float(text) == v else repr(v)     # two decimals as a rule, more only when they would lose dx


def prepare_name(fe_run: str, k0: int, k1: int, every: int = 1) -> str:
    """'prep_<fe_run>_k00000-01190' (+ '_every2' when every != 1)."""
    _check_part(fe_run, "finite-element run name")
    if isinstance(every, bool) or not isinstance(every, int) or every < 1:
        raise ValueError("every must be a positive integer, got {!r}".format(every))
    return "prep_{}_{}".format(fe_run, _frames(k0, k1)) + ("_every{}".format(every) if every != 1 else "")


def parse_prepare_name(name: str) -> dict:
    """Inverse of prepare_name: {'fe_run', 'k0', 'k1', 'every'}; ValueError if `name` is not one."""
    m = _PREPARE_RE.match(name) if isinstance(name, str) else None
    if m is None or "__" in name:
        raise ValueError("not a prepared-input name: {!r}".format(name))
    return {"fe_run": m["fe_run"], "k0": int(m["k0"]), "k1": int(m["k1"]),
            "every": int(m["every"]) if m["every"] else 1}


def bracket_suffix(brackets: dict) -> str:
    """'_<label>-<value>' for every bracket that differs from its default, in BRACKET_DEFAULTS order."""
    unknown = sorted(set(brackets) - set(BRACKET_DEFAULTS))
    if unknown:
        raise ValueError("unknown bracket setting(s): {}".format(", ".join(unknown)))
    parts = []
    for key, default in BRACKET_DEFAULTS.items():
        if key not in brackets:
            continue
        integer = isinstance(default, int) and not isinstance(default, bool)
        text = _format_number(brackets[key], integer)
        if text == _format_number(default, integer):
            continue
        label, unit = BRACKET_LABELS[key]
        parts.append("_{}-{}{}".format(label, text, unit))
    return "".join(parts)


def run_name(prepared: str, mode: str, frames: tuple[int, int], form: str, dx_mm: float,
             brackets: dict, seed: int) -> str:
    """'<prepared>__replay_k00051-00450_rz_dx1.10mm_seed1' + '_filmlimit-3mm' ... (module docstring)."""
    _check_part(prepared, "prepared-input name")
    if mode not in MODES:
        raise ValueError("mode must be one of {}, got {!r}".format(MODES, mode))
    if form not in FORMS:
        raise ValueError("form must be one of {}, got {!r}".format(FORMS, form))
    k0, k1 = frames
    _check_index(seed, "seed")       # always present: three seeds per case are routine (spec §8.4)
    return "{}__{}_{}_{}_dx{}mm_seed{}{}".format(prepared, mode, _frames(k0, k1), form, _format_dx(dx_mm), seed,
                                                 bracket_suffix(brackets or {}))


def analyse_name(run: str) -> str:
    """The analysis of run `run` is spheral_output/analyse/<run>/ (M1 plan, Global constraints): the run's own name,
    checked to be one. Analysis settings that change the result (min_particles, link_h) are brackets of the run."""
    parse_run_name(run)
    return run


def parse_run_name(name) -> dict:
    """Inverse of run_name, for resumable studies and tests.

    Returns {'prepared', 'mode', 'frames', 'form', 'dx_mm', 'seed', 'brackets', 'prepare'}: 'brackets' holds every
    key of BRACKET_DEFAULTS (the default where the name has no suffix), 'prepare' the parsed prepared name (None when
    the prepared part is not a prepare_name). ValueError if `name` is not a run name."""
    if not isinstance(name, str) or name.count("__") != 1:
        raise ValueError("not a run name (one '__' expected): {!r}".format(name))
    prepared, own = name.split("__")
    tokens = own.split("_")
    if len(tokens) < 5:
        raise ValueError("not a run name (too few parts): {!r}".format(name))
    mode, frames, form, dx, seed = tokens[:5]
    mf, mdx, ms = _FRAMES_RE.match(frames), _DX_RE.match(dx), _SEED_RE.match(seed)
    if mode not in MODES or form not in FORMS or not (mf and mdx and ms):
        raise ValueError("not a run name: {!r}".format(name))
    brackets = dict(BRACKET_DEFAULTS)
    seen = []
    for token in tokens[5:]:
        label, sep, text = token.partition("-")
        key = _KEY_OF_LABEL.get(label)
        if not sep or key is None or key in seen:
            raise ValueError("bad bracket suffix {!r} in {!r}".format(token, name))
        unit = BRACKET_LABELS[key][1]
        if unit:
            if not text.endswith(unit):
                raise ValueError("bracket suffix {!r} lacks its unit {!r}".format(token, unit))
            text = text[:-len(unit)]
        default = BRACKET_DEFAULTS[key]
        brackets[key] = int(text) if isinstance(default, int) and not isinstance(default, bool) else float(text)
        seen.append(key)
    order = [k for k in BRACKET_DEFAULTS if k in seen]
    if seen != order:
        raise ValueError("bracket suffixes out of order in {!r}".format(name))
    try:
        prepare = parse_prepare_name(prepared)
    except ValueError:
        prepare = None
    out = {"prepared": prepared, "mode": mode, "frames": (int(mf["k0"]), int(mf["k1"])), "form": form,
           "dx_mm": float(mdx["dx"]), "seed": int(ms["seed"]), "brackets": brackets, "prepare": prepare}
    if run_name(prepared, mode, out["frames"], form, out["dx_mm"], brackets, out["seed"]) != name:
        raise ValueError("not a canonical run name: {!r}".format(name))
    return out
