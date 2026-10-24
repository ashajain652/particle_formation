"""Load tables by inclination (spec §7.2-7.3; M1 plan, Task 7; check "Load tables").

The windward part of Spheral's surface loads comes from the frame's own wall pressure `p_w` and shear `tau`, binned
by the patch's inclination theta (the angle between its outward normal and the flight direction v_hat), so that on
the undeformed body Spheral sees the finite-element model's loads. Behind the last bin the frame carries loads for,
a lee extension stands in for Step 4's lee-load module (spec §7.2; plan decision 5):

  pressure  steps to the base pressure `base_fraction * p_stag` and holds it out to 180 deg, where p_stag is the
            history's wall stagnation pressure `p_w_stag_step_Pa` at the frame's row (Step 4's "pitot"); the frame's own
            p_w already contains Step 3's Prandtl-Meyer expansion up to its last bin, so no ramp is added
  shear     `shear_factor` times the shear of the last valid bin up to `separation_deg`, zero beyond

Brackets (spec §7.2): base pressure 1, 3 (default) or 5 %; separation 180 deg (default) or 150 deg; lee shear 0, 0.5
(default) or 1 times the last value. The defaults are `naming.BRACKET_DEFAULTS`' `base_pressure_pct`,
`separation_deg` and `lee_shear`. The extension is labelled `LEE_LABEL` in every table header.

Conventions:

- `theta` arrays are in **radians** (as Task 4's `patch_theta`); the table's bin centres and `evaluate`'s argument
  are in **degrees** (`theta_deg`). theta must come from the frame's **outward** normals (Task 4's
  `derived_patch_arrays`); a face wound inward would be binned at 180 deg - theta.
- A patch is valid where `p_w > 0` (contract: 0 = "not evaluated", e.g. a frame written before any step, or a face
  exposed by that step's deaths) and its shear is finite. A frame with no valid patch gives an **empty table**: every
  bin NaN, `theta_last_deg` NaN, `has_loads` False; `with_lee` returns it unchanged (still empty), `evaluate` gives
  NaN, and `drag_comparison` reports `loads_valid` False with the table drags NaN.
- Drag is the force along -v_hat (toward the tail; the body decelerates): D = sum A (p cos theta + tau sin theta).
  `tau` is the magnitude along the tangent of the oncoming flow projected onto the surface (contract), whose
  component along -v_hat is sin theta.

Numpy and the standard library only (the runner imports this module under Spheral's Python)."""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from . import contract, naming

G0 = 9.80665                     # m/s^2, the finite element's G0: its history's load_factor_g = |a_drag| / G0
THETA_EDGES_DEG = np.arange(0.0, 181.0, 1.0)
LEE_LABEL = "Step 4 §5's lee model, to be replaced by its module"
TABLE_SCHEMA_VERSION = 1

BASE_FRACTION = naming.BRACKET_DEFAULTS["base_pressure_pct"] / 100.0     # 0.03
SEPARATION_DEG = float(naming.BRACKET_DEFAULTS["separation_deg"])        # 180.0
LEE_SHEAR = float(naming.BRACKET_DEFAULTS["lee_shear"])                  # 0.5


@dataclass
class LoadTable:
    theta_deg: np.ndarray      # bin centres (n_bins,), deg
    p: np.ndarray              # Pa, area-weighted mean p_w per bin (NaN where no valid patch)
    tau: np.ndarray            # Pa, area-weighted mean tau per bin (NaN where no valid patch)
    area: np.ndarray           # m^2 of valid patches per bin
    theta_last_deg: float      # centre of the last bin with valid patches ("the last patch the frame carries") or NaN
    p_stag: float              # Pa, the history's p_w_stag_step_Pa at this frame
    edges_deg: np.ndarray = field(default_factory=lambda: THETA_EDGES_DEG.copy())   # (n_bins + 1,)
    lee: np.ndarray | None = None    # bool (n_bins,): bins set by the lee extension (all False before `with_lee`)
    lee_params: dict | None = None   # base_fraction, separation_deg, shear_factor once `with_lee` applied

    def __post_init__(self):
        if self.lee is None:
            self.lee = np.zeros(len(self.theta_deg), dtype=bool)

    @property
    def has_loads(self) -> bool:
        return bool(np.isfinite(self.theta_last_deg))

    @property
    def windward_edge_deg(self) -> float:
        """Upper edge of the last valid bin: patches at or beyond it lie in the lee extension's range."""
        if not self.has_loads:
            return float("nan")
        i = int(np.flatnonzero(self.theta_deg == self.theta_last_deg)[0])
        return float(self.edges_deg[i + 1])


def valid_patches(p_w, tau) -> np.ndarray:
    """Default validity: evaluated (p_w > 0; NaN counts as not evaluated) with a finite shear."""
    p_w = np.asarray(p_w, dtype=float)
    tau = np.asarray(tau, dtype=float)
    return (p_w > 0.0) & np.isfinite(tau)


def bin_index(theta_deg, edges=THETA_EDGES_DEG) -> np.ndarray:
    """Bin of each inclination: edges[i] <= theta < edges[i + 1]; the last bin also takes theta = edges[-1], and
    values outside the edges go to the nearest end bin."""
    edges = np.asarray(edges, dtype=float)
    i = np.searchsorted(edges, np.asarray(theta_deg, dtype=float), side="right") - 1
    return np.clip(i, 0, len(edges) - 2)


def build_table(theta, area, p_w, tau, p_stag, valid=None, edges=THETA_EDGES_DEG) -> LoadTable:
    """The frame's loads binned by inclination. `theta` in rad (outward normals against v_hat), `area` m^2, `p_w`
    and `tau` Pa per patch; `p_stag` the history's p_w_stag_step_Pa at the frame (Pa). `valid` defaults to
    `valid_patches(p_w, tau)` (plan fact 2: 0 is the "not evaluated" default)."""
    theta = np.asarray(theta, dtype=float)
    area = np.asarray(area, dtype=float)
    p_w = np.asarray(p_w, dtype=float)
    tau = np.asarray(tau, dtype=float)
    edges = np.asarray(edges, dtype=float)
    if not (theta.shape == area.shape == p_w.shape == tau.shape) or theta.ndim != 1:
        raise ValueError("theta, area, p_w and tau must be 1-D arrays of one length")
    if edges.ndim != 1 or len(edges) < 2 or np.any(np.diff(edges) <= 0.0):
        raise ValueError("edges must be strictly increasing with at least two entries")
    valid = valid_patches(p_w, tau) if valid is None else np.asarray(valid, dtype=bool)
    n_bins = len(edges) - 1
    idx = bin_index(np.degrees(theta[valid]), edges)
    a = area[valid]
    A = np.bincount(idx, weights=a, minlength=n_bins)
    sum_p = np.bincount(idx, weights=a * p_w[valid], minlength=n_bins)
    sum_tau = np.bincount(idx, weights=a * tau[valid], minlength=n_bins)
    filled = A > 0.0
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(filled, sum_p / np.where(filled, A, 1.0), np.nan)
        t = np.where(filled, sum_tau / np.where(filled, A, 1.0), np.nan)
    centres = 0.5 * (edges[:-1] + edges[1:])
    last = float(centres[np.flatnonzero(filled)[-1]]) if filled.any() else float("nan")
    return LoadTable(theta_deg=centres, p=p, tau=t, area=A, theta_last_deg=last, p_stag=float(p_stag),
                     edges_deg=edges.copy())


def with_lee(table: LoadTable, base_fraction=BASE_FRACTION, separation_deg=SEPARATION_DEG,
             shear_factor=LEE_SHEAR) -> LoadTable:
    """The table with the lee extension (spec §7.2, plan decision 5) in every bin beyond `theta_last_deg`: pressure
    `base_fraction * p_stag`; shear `shear_factor` times the last valid bin's shear for bin centres up to
    `separation_deg`, zero beyond. Bins at or before the last valid one are left as built (empty interior bins stay
    NaN; `evaluate` interpolates across them). An empty table is returned unchanged (copied, still empty)."""
    params = {"base_fraction": float(base_fraction), "separation_deg": float(separation_deg),
              "shear_factor": float(shear_factor)}
    if not table.has_loads:
        return replace(table, p=table.p.copy(), tau=table.tau.copy(), area=table.area.copy(),
                       lee=np.zeros(len(table.theta_deg), dtype=bool), lee_params=params)
    lee = table.theta_deg > table.theta_last_deg
    i_last = int(np.flatnonzero(table.theta_deg == table.theta_last_deg)[0])
    tau_last = float(table.tau[i_last])
    p = table.p.copy()
    tau = table.tau.copy()
    p[lee] = base_fraction * table.p_stag
    tau[lee] = np.where(table.theta_deg[lee] <= separation_deg, shear_factor * tau_last, 0.0)
    return replace(table, p=p, tau=tau, area=table.area.copy(), lee=lee, lee_params=params)


def evaluate(table: LoadTable, theta_deg) -> tuple[np.ndarray, np.ndarray]:
    """(p, tau) in Pa at inclinations `theta_deg` (deg), by np.interp on the centres of the non-empty bins (held
    constant outside the first and last of them). NaN everywhere for an empty table."""
    theta_deg = np.asarray(theta_deg, dtype=float)
    ok = np.isfinite(table.p)
    if not ok.any():
        nan = np.full(theta_deg.shape, np.nan)
        return nan, nan.copy()
    c = table.theta_deg[ok]
    return np.interp(theta_deg, c, table.p[ok]), np.interp(theta_deg, c, table.tau[ok])


def drag(theta, area, p, tau) -> float:
    """D = sum A (p cos(theta) + tau sin(theta)), N, along -v_hat (the flight axis, toward the tail). theta in rad."""
    theta = np.asarray(theta, dtype=float)
    return float(np.sum(np.asarray(area, dtype=float) * (np.asarray(p, dtype=float) * np.cos(theta)
                                                         + np.asarray(tau, dtype=float) * np.sin(theta))))


def smooth_drag(theta, area, theta_smooth, cos_to_smooth, p, tau) -> float:
    """The drag of the faceted body under the smooth surface's loads (N, along -v_hat), with `p` and `tau` the loads
    at the smooth inclination `theta_smooth` (rad) and `cos_to_smooth` = n_facet . n_smooth per patch:

      D = sum A [p cos(theta) + tau sin(theta_smooth) max(cos_to_smooth, 0)]

    Pressure acts normal to the facet, so each facet carries p on its own projected area A cos(theta): those areas
    sum to the body's frontal area exactly and a uniform pressure gives no force on the closed surface. Shear acts
    along the smooth surface's tangent (component sin(theta_smooth) along -v_hat) on the facet's area projected onto
    the smooth surface, A (n_facet . n_smooth). Reading p at the facet's own theta instead (`drag` on the staircase)
    loses most of the drag on a staircase whose facets scatter by tens of degrees about the smooth surface, since p
    falls like cos^2 of the facet's own tilt (prototype README, findings of 2026-10-08)."""
    theta = np.asarray(theta, dtype=float)
    shear_area = np.asarray(area, dtype=float) * np.maximum(np.asarray(cos_to_smooth, dtype=float), 0.0)
    return float(np.sum(np.asarray(area, dtype=float) * np.asarray(p, dtype=float) * np.cos(theta)
                        + shear_area * np.asarray(tau, dtype=float) * np.sin(np.asarray(theta_smooth, dtype=float))))


def history_drag(history_row, g0=G0) -> float:
    """The history's drag at a row: mass_kg * load_factor_g * g0 (N; plan fact 4: there is no drag column).
    `history_row` maps contract history keys to values."""
    return float(history_row[contract.fe_name("mass_kg")]) * float(history_row[contract.fe_name("load_factor_g")]) \
        * float(g0)


def history_p_stag(history_row) -> float:
    """The lee base pressure's reference: the history's p_w_stag_step_Pa at the row (Pa; NaN on row 0)."""
    return float(history_row[contract.fe_name("p_w_stag_step_Pa")])


def drag_comparison(theta, area, p_w, tau, history_row, valid=None, edges=THETA_EDGES_DEG,
                    base_fraction=BASE_FRACTION, separation_deg=SEPARATION_DEG, shear_factor=LEE_SHEAR,
                    g0=G0, theta_smooth=None, cos_to_smooth=None) -> dict:
    """The drags of plan Task 7, Step 2, for one frame (N), and their split into windward and lee parts:

      D_smooth           the headline against the history: the lee-extended table evaluated at every patch's smooth
                         inclination `theta_smooth` (the derived surface's normal, n_derived) and summed by
                         `smooth_drag` with `cos_to_smooth` = n_facet . n_derived; without them theta_smooth = theta
                         and cos_to_smooth = 1, and D_smooth equals D_table_lee to round-off
      D_smooth_lee       D_smooth over the patches whose smooth inclination is at or beyond the windward edge
      smooth_lee_share   D_smooth_lee / D_smooth
      rel_smooth_hist    D_smooth / D_hist - 1 (spec §7.3's comparison)
      drag_normals       "derived" or "facet": which inclination D_smooth read the loads at
      D_patch            the frame's own p_w and tau summed over its patches (staircase facets at their own theta;
                         a diagnostic since 2026-10-08)
      D_table            the windward table evaluated at each valid patch's theta and summed (the binning error)
      D_hist             the history's mass_kg * load_factor_g * g0 at the frame's row
      D_table_lee        the lee-extended table evaluated at every patch (what Spheral's undeformed body would see)
      D_windward_patch   D_patch over the patches before the last valid bin's upper edge ("windward")
      D_windward_table   D_table_lee over the same patches
      D_lee              D_table_lee over the patches at or beyond that edge (the lee extension's part)
      lee_share          D_lee / D_table_lee
      rel_table_patch    D_table / D_patch - 1;  rel_patch_hist  D_patch / D_hist - 1;
      rel_table_lee_hist D_table_lee / D_hist - 1

    plus `loads_valid` (the frame carries any loads), `n_valid`, `n_lee_patches`, `theta_last_deg` and `p_stag`.
    On a frame without loads the table drags and shares are NaN, and D_patch is 0 (rel_patch_hist = -1).

    The table is binned by the facets' own theta in every case: the frame's p_w and tau were evaluated there, so the
    table is the surface flow's p(theta), and only where it is read changes."""
    theta = np.asarray(theta, dtype=float)
    area = np.asarray(area, dtype=float)
    p_w = np.asarray(p_w, dtype=float)
    tau = np.asarray(tau, dtype=float)
    if (theta_smooth is None) != (cos_to_smooth is None):
        raise ValueError("theta_smooth and cos_to_smooth go together")
    drag_normals = "facet" if theta_smooth is None else "derived"
    theta_s = theta if theta_smooth is None else np.asarray(theta_smooth, dtype=float)
    cos_s = np.ones_like(theta) if cos_to_smooth is None else np.asarray(cos_to_smooth, dtype=float)
    if theta_s.shape != theta.shape or cos_s.shape != theta.shape:
        raise ValueError("theta_smooth and cos_to_smooth must have theta's shape")
    valid = valid_patches(p_w, tau) if valid is None else np.asarray(valid, dtype=bool)
    p_stag = history_p_stag(history_row)
    table = build_table(theta, area, p_w, tau, p_stag, valid=valid, edges=edges)
    lee_table = with_lee(table, base_fraction, separation_deg, shear_factor)
    deg = np.degrees(theta)
    own_p = np.where(valid, p_w, 0.0)
    own_tau = np.where(valid, tau, 0.0)
    D_patch = drag(theta, area, own_p, own_tau)
    D_hist = history_drag(history_row, g0)
    nan = float("nan")
    out = {"loads_valid": table.has_loads, "n_valid": int(valid.sum()), "theta_last_deg": table.theta_last_deg,
           "p_stag": p_stag, "D_patch": D_patch, "D_hist": D_hist, "drag_normals": drag_normals,
           "rel_patch_hist": D_patch / D_hist - 1.0 if D_hist != 0.0 else nan}
    if not table.has_loads:
        out.update({"D_table": nan, "D_table_lee": nan, "D_windward_patch": D_patch, "D_windward_table": nan,
                    "D_lee": nan, "lee_share": nan, "n_lee_patches": 0, "rel_table_patch": nan,
                    "rel_table_lee_hist": nan, "D_smooth": nan, "D_smooth_lee": nan, "smooth_lee_share": nan,
                    "rel_smooth_hist": nan})
        return out
    pt, tt = evaluate(table, deg[valid])
    D_table = drag(theta[valid], area[valid], pt, tt)
    pl, tl = evaluate(lee_table, deg)
    lee_side = deg >= table.windward_edge_deg
    D_windward_table = drag(theta[~lee_side], area[~lee_side], pl[~lee_side], tl[~lee_side])
    D_lee = drag(theta[lee_side], area[lee_side], pl[lee_side], tl[lee_side])
    D_table_lee = D_windward_table + D_lee
    deg_s = np.degrees(theta_s)
    ps, ts = evaluate(lee_table, deg_s)
    lee_s = deg_s >= table.windward_edge_deg
    D_smooth_wind = smooth_drag(theta[~lee_s], area[~lee_s], theta_s[~lee_s], cos_s[~lee_s], ps[~lee_s], ts[~lee_s])
    D_smooth_lee = smooth_drag(theta[lee_s], area[lee_s], theta_s[lee_s], cos_s[lee_s], ps[lee_s], ts[lee_s])
    D_smooth = D_smooth_wind + D_smooth_lee
    out.update({
        "D_smooth": D_smooth, "D_smooth_lee": D_smooth_lee,
        "smooth_lee_share": D_smooth_lee / D_smooth if D_smooth != 0.0 else nan,
        "rel_smooth_hist": D_smooth / D_hist - 1.0 if D_hist != 0.0 else nan,
        "D_table": D_table, "D_table_lee": D_table_lee,
        "D_windward_patch": drag(theta[~lee_side], area[~lee_side], own_p[~lee_side], own_tau[~lee_side]),
        "D_windward_table": D_windward_table, "D_lee": D_lee,
        "lee_share": D_lee / D_table_lee if D_table_lee != 0.0 else nan,
        "n_lee_patches": int(lee_side.sum()),
        "rel_table_patch": D_table / D_patch - 1.0 if D_patch != 0.0 else nan,
        "rel_table_lee_hist": D_table_lee / D_hist - 1.0 if D_hist != 0.0 else nan,
    })
    return out


# ------------------------------------------------------------------------------------------------ files (Task 9)
def table_header(table: LoadTable) -> dict:
    """The JSON header of a table or a stack of them: the binning, the lee extension's label and parameters."""
    return {"schema": "spheral_frag.loads", "schema_version": TABLE_SCHEMA_VERSION,
            "edges_deg": [float(e) for e in table.edges_deg], "valid": "p_w > 0 and finite tau",
            "lee_label": LEE_LABEL, "lee_params": table.lee_params,
            "lee_reference_pressure": contract.fe_name("p_w_stag_step_Pa"),
            "units": {"p": "Pa", "tau": "Pa", "area": "m^2", "theta_deg": "deg", "p_stag": "Pa"}}


def stack(tables) -> dict[str, np.ndarray]:
    """Per-frame tables (one binning) as plain numeric arrays for `loads.npz`: p, tau, area, lee (n_frames, n_bins),
    theta_last_deg, p_stag (n_frames,), theta_deg, edges_deg."""
    tables = list(tables)
    if not tables:
        raise ValueError("no tables to stack")
    edges = tables[0].edges_deg
    if any(not np.array_equal(t.edges_deg, edges) for t in tables):
        raise ValueError("tables of different binnings cannot be stacked")
    return {"theta_deg": tables[0].theta_deg.copy(), "edges_deg": edges.copy(),
            "p": np.stack([t.p for t in tables]), "tau": np.stack([t.tau for t in tables]),
            "area": np.stack([t.area for t in tables]), "lee": np.stack([t.lee for t in tables]).astype(np.int8),
            "theta_last_deg": np.array([t.theta_last_deg for t in tables]),
            "p_stag": np.array([t.p_stag for t in tables])}


def unstack(arrays, lee_params=None) -> list[LoadTable]:
    """Inverse of `stack` (numpy only; the runner's reader)."""
    return [LoadTable(theta_deg=np.array(arrays["theta_deg"], dtype=float), p=np.array(arrays["p"][i], dtype=float),
                      tau=np.array(arrays["tau"][i], dtype=float), area=np.array(arrays["area"][i], dtype=float),
                      theta_last_deg=float(arrays["theta_last_deg"][i]), p_stag=float(arrays["p_stag"][i]),
                      edges_deg=np.array(arrays["edges_deg"], dtype=float),
                      lee=np.array(arrays["lee"][i], dtype=bool), lee_params=lee_params)
            for i in range(len(arrays["p_stag"]))]
