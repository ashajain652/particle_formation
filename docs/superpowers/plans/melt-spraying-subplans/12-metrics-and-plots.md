# Sub-plan: Task 12 — Comparison metrics, plots and visualisation of the melting run

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 6130–6651). Read `00-shared-context.md` first.

**Depends on:** Task 10 (melt-history columns) and Task 11 (the two committed reference runs).
**Produces, for later tasks:** comparison metrics (including an interpolated mass-crossing time), seven new plots, and visualisation overlays (emitting patches, film-thickness frame, liquidus/solidus isolines, melt stills) that Task 13 exposes through the CLI and Task 14 uses to report results.
**Character:** plumbing plus plotting.
**Read before implementing:** none beyond the shared context; just confirm Task 11's reference files exist before this task's tests run against them.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan.

---

### Task 12: Comparison metrics, plots and visualisation of the melting run

**Files:**
- Modify: `reentry_model/sesam_io.py` (two edits), `reentry_model/compare.py` (append + two constants)
- Modify (replace): `reentry_model/viz.py`
- Test: `tests/test_reentry_model_compare.py` (append), `tests/test_reentry_model_viz.py` (append)

**Interfaces:**
- Consumes: the melting reference files (Task 11), histories with `MELT_COLUMNS` (Task 10), VTK series with the melt fields (Task 10).
- Produces: `sesam_io.Reference.mass`, `.thickness`; `compare.MELT_PLOT_NAMES` (7, with `closures_time.png` carrying the closure area fractions and both Knudsen numbers on a twin log axis), `DEMISE_FRACTION = 0.01`, `has_melt(history, reference=None)`, `melt_metrics(history, reference) -> dict` (`mass.max_rel_m0`, `onset_altitude_diff_km`, `demise_time_diff_s`, `demise_time_rel`, ...), `plot_melt(history, outdir, title, reference=None, size_distribution_csv=None) -> paths`; `viz.render_frame(..., melting=False)`, `render_film_frame`, `render_section_frame(..., melting=False)`, `animate(..., melting=False)`, `animate_section(..., melting=False)`, `animate_film(run_dir, history, radius, fps, animation, stills)`, `still_marks(history, melting=False)`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_compare.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: the melting reference, the mass metrics and plots

MELT100 = "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind"


@pytest.fixture(scope="module")
def melt_ref():
    return sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, MELT100 + ".csv"))


def melting_history(ref, lag=0.0, film=0.0):
    """A synthetic model history that reproduces the reference's mass shifted by `lag` seconds, with the melt columns."""
    h = synthetic_history(ref, dv=0.0, dh=0.0)
    c = h.columns
    n = len(ref.time)
    c["mass_kg"] = np.interp(ref.time - lag, ref.time, ref.mass) + film
    for key in ("film_mass_kg", "sprayed_mass_kg", "runoff_mass_kg", "removed_mass_kg", "spraying_area_m2", "theta_cr_deg", "n_released",
                "released_mass_kg", "r_median_um", "r_max_um", "closure_fraction_girin", "closure_fraction_couette_slip",
                "closure_fraction_couette_fm", "kn_body", "kn_local_stag", "film_thickness_mean_mm", "film_thickness_max_mm",
                "convective_heat_W", "surface_T_max_K"):
        c[key] = np.zeros(n)
    c["film_mass_kg"][:] = film
    c["removed_mass_kg"] = ref.mass[0] - c["mass_kg"] + film
    c["sprayed_mass_kg"] = c["removed_mass_kg"].copy()
    c["r_median_um"][:] = 150.0
    c["closure_fraction_girin"][:] = 0.5
    c["closure_fraction_couette_slip"][:] = 0.5
    c["kn_body"][:] = 0.01
    c["kn_local_stag"][:] = 1e-4
    return h


def test_melting_reference_and_has_melt(ref, melt_ref):
    assert ref.mass is not None and float(ref.mass.min()) == float(ref.mass[0])                         # no-melt: constant mass
    assert melt_ref.mass[0] == pytest.approx(1.473) and melt_ref.mass[-1] == 0.0 and melt_ref.thickness[0] == pytest.approx(0.05)
    assert melt_ref.thickness[-1] == 0.0 and np.all(np.diff(melt_ref.mass) <= 0.0)
    h = melting_history(melt_ref)
    assert compare.has_melt(h) and compare.has_melt(h, melt_ref) and not compare.has_melt(h, ref) and not compare.has_melt(synthetic_history(ref))


def test_melt_metrics_recover_a_known_lag(melt_ref):
    m = compare.melt_metrics(melting_history(melt_ref, lag=-1.0), melt_ref)          # the model runs 1 s ahead of SESAM
    assert m["mass"]["max_rel_m0"] < 0.06 and m["mass"]["max_abs_kg"] > 0.0
    assert m["demise_time_diff_s"] == pytest.approx(-1.0, abs=0.1) and m["demise_time_rel"] == pytest.approx(-1.0 / m["demise_time_reference_s"], rel=0.1)
    assert m["onset_time_reference_s"] == pytest.approx(43.55) and m["onset_altitude_reference_km"] == pytest.approx(71.005, abs=0.01)
    assert m["onset_altitude_diff_km"] is not None and m["demise_time_reference_s"] == pytest.approx(66.8, abs=0.1)
    exact = compare.melt_metrics(melting_history(melt_ref), melt_ref)
    assert exact["mass"]["max_rel_m0"] < 1e-12 and abs(exact["demise_time_diff_s"]) < 1e-9
    with_film = compare.melt_metrics(melting_history(melt_ref, film=0.05), melt_ref)
    assert with_film["demise_time_diff_s"] == pytest.approx(exact["demise_time_diff_s"], abs=1e-9)   # the crossing is on the body's material


def test_melt_plots_are_written(melt_ref, tmp_path):
    h = melting_history(melt_ref, lag=-1.0)
    paths = compare.plot_melt(h, str(tmp_path), "melt test", melt_ref)
    assert [os.path.basename(p) for p in paths] == list(compare.MELT_PLOT_NAMES) and all(os.path.getsize(p) > 5000 for p in paths)
    alone = compare.plot_melt(h, str(tmp_path / "alone"), "no reference")
    assert len(alone) == len(compare.MELT_PLOT_NAMES) and all(os.path.isfile(p) for p in alone)
```


Append to `tests/test_reentry_model_viz.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: overlays of the emitting patches, the film-thickness frame, the liquidus iso-line, the melt stills

def test_melting_overlays_and_film_frame(coarse_sphere_mesh, tmp_path):
    import pyvista as pv
    surface = coarse_sphere_mesh.surface()
    faces = np.hstack([np.full((surface.n_patches, 1), 3), surface.faces]).ravel()
    poly = pv.PolyData(coarse_sphere_mesh.points, faces)
    poly.point_data["T"] = 300.0 + 600.0 * np.clip(coarse_sphere_mesh.points[:, 0] / 0.05, 0.0, 1.0)
    r = np.full(surface.n_patches, np.nan)
    windward = surface.normals[:, 0] > 0.5
    r[windward] = 1e-4
    poly.cell_data["r_droplet"] = r
    poly.cell_data["film_thickness"] = np.where(windward, 2e-5, 0.0)
    plain = viz.render_frame(poly, (300.0, 900.0), "t", None, 0.05)
    overlaid = viz.render_frame(poly, (300.0, 900.0), "t", str(tmp_path / "o.png"), 0.05, melting=True)
    assert overlaid.shape == plain.shape and np.abs(overlaid.astype(int) - plain.astype(int)).mean() > 1.0     # the overlay changed the picture
    film = viz.render_film_frame(poly, None, "film", str(tmp_path / "f.png"), 0.05)
    assert film.shape == (720, 960, 3) and film.std() > 10.0 and os.path.getsize(tmp_path / "f.png") > 1000
    m = coarse_sphere_mesh
    cells = np.hstack([np.full((m.n_elements, 1), 4), m.tets]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(m.n_elements, pv.CellType.TETRA), m.points)
    grid.point_data["T"] = 300.0 + 700.0 * np.clip(m.points[:, 0] / 0.05, 0.0, 1.0)
    grid.point_data["liquid_fraction"] = np.clip((grid.point_data["T"] - 750.0) / 158.0, 0.0, 1.0)
    without = viz.render_section_frame(grid, (300.0, 1000.0), "t", None, 0.05)
    with_iso = viz.render_section_frame(grid, (300.0, 1000.0), "t", str(tmp_path / "s.png"), 0.05, melting=True)
    assert with_iso.shape == without.shape and np.abs(with_iso.astype(int) - without.astype(int)).mean() > 0.1


def test_still_marks_and_film_animation(coarse_sphere_mesh, tmp_path):
    pytest.importorskip("cantera")
    from test_reentry_model_coupled import simulator, MASS_100MM
    from reentry_model import body, coupled, heating, material, mesh, thermal
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"), MASS_100MM, T0=700.0)
    sim = simulator(b, t_max=53.0)
    sim.advance(43.0)
    run_dir = str(tmp_path / "run")
    settings = coupled.CoupledSettings(dt=0.5, frames_every=10, output_dir=run_dir)
    hist = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), settings).run()
    marks = viz.still_marks(hist, melting=True)
    assert {"start", "peak_heating", "peak_surface_T", "end", "melt_onset", "spraying_onset", "peak_release"} <= set(marks)
    assert marks["melt_onset"] <= marks["spraying_onset"] <= marks["end"]
    out = viz.animate(run_dir, hist, 0.05, fps=5, animation=False, melting=True)
    assert len(out["stills"]) == 7 and all(os.path.isfile(p) for p in out["stills"])
    film = viz.animate_film(run_dir, hist, 0.05, fps=5, animation=False)
    assert len(film["stills"]) == 7 and all("film_" in os.path.basename(p) for p in film["stills"])
    sec = viz.animate_section(run_dir, hist, 0.05, fps=5, animation=False, melting=True)
    assert any(os.path.basename(p).startswith("section_melt_onset") for p in sec["stills"])
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py tests/test_reentry_model_viz.py -q`
Expected: the five new tests fail (`Reference` has no `mass`, no `has_melt`, `render_frame` has no `melting`).

- [ ] **Step 3: Edit `reentry_model/sesam_io.py`**

In the `Reference` dataclass, after `integrated_heat: np.ndarray = None      # J, integral of the convective heat (SESAM's integrated_heat_J)` add:

```python
    mass: np.ndarray = None                 # kg (SESAM's mass_kg; decreases while the sphere melts)
    thickness: np.ndarray = None            # m (SESAM's thick_mm: the shell thickness of the melting sphere, Step 3 facts)
```

and in `load_reference`, after the `convective_heat=..., integrated_heat=optional("integrated_heat_J"),` line add:

```python
        mass=optional("mass_kg"), thickness=col("thick_mm", 1e-3) if "thick_mm" in rows[0] else None,
```

- [ ] **Step 4: Extend `reentry_model/compare.py`**

Replace the two constant lines `THERMAL_PLOT_NAMES = (...)` / `CONTINUUM_KN = 0.01` with:

```python
THERMAL_PLOT_NAMES = ("heating_time.png", "temperature_time.png", "integrated_heat.png")
MELT_PLOT_NAMES = ("mass_time.png", "mass_altitude.png", "mass_budget.png", "spraying_time.png", "closures_time.png",
                   "droplet_size_time.png", "size_distribution.png")
CONTINUUM_KN = 0.01
DEMISE_FRACTION = 0.01
```

and append to the file:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: mass loss, spraying and size distributions (spec Step 3 sections 12-13)

def has_melt(history, reference=None):
    """True when the model history carries the melt columns (and, if a reference is given, its mass varies)."""
    ok = "sprayed_mass_kg" in history.columns and "film_mass_kg" in history.columns
    if reference is None:
        return ok
    return ok and reference.mass is not None and float(reference.mass.min()) < 0.999 * float(reference.mass[0])


def _first_time(t, mask):
    return float(t[mask][0]) if mask.any() else None


def _crossing_time(t, m, level):
    """Time at which the decreasing series m first falls below `level`, interpolated linearly between samples."""
    below = np.flatnonzero(m < level)
    if below.size == 0:
        return None
    i = int(below[0])
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (m[i - 1] - level) / (m[i - 1] - m[i]) * (t[i] - t[i - 1]))


def melt_metrics(history, reference):
    """Mass vs SESAM at the reference's time stamps (relative to the initial mass), melt-onset altitude and the time
    at which the mass falls below DEMISE_FRACTION of the initial (SESAM removes the last gram at its 'demise'; the model
    stops at the fraction, so the 1 % crossing is the like-for-like end time; for the model it is taken on the body's
    material, film excluded, since a static leeward film can stay attached)."""
    c = history.columns
    t_mod, m_mod = c["time_s"], c["mass_kg"]
    mask = (reference.time >= t_mod[0]) & (reference.time <= t_mod[-1])
    t = reference.time[mask]
    m_ref = reference.mass[mask]
    m_at = np.interp(t, t_mod, m_mod)
    m0 = float(reference.mass[0])
    onset_ref = _first_time(reference.time, reference.mass < 0.999999 * m0)
    onset_mod = _first_time(t_mod, (c["removed_mass_kg"] + c["film_mass_kg"]) > 0.0)
    h_ref = float(np.interp(onset_ref, reference.time, reference.altitude)) / 1e3 if onset_ref is not None else None
    h_mod = float(np.interp(onset_mod, t_mod, c["altitude_km"])) if onset_mod is not None else None
    end_ref = _crossing_time(reference.time, reference.mass, DEMISE_FRACTION * m0)
    body_mass = m_mod - c["film_mass_kg"]                                                          # the body's material (a leeward film may stay attached)
    end_mod = _crossing_time(t_mod, body_mass, DEMISE_FRACTION * float(body_mass[0]))
    return {
        "n_points": int(t.size), "initial_mass_model_kg": float(m_mod[0]), "initial_mass_reference_kg": m0,
        "mass": {"max_abs_kg": float(np.abs(m_at - m_ref).max()), "max_rel_m0": float(np.abs(m_at - m_ref).max() / m0),
                 "rms_rel_m0": float(math.sqrt(np.mean((m_at - m_ref) ** 2)) / m0)},
        "onset_time_model_s": onset_mod, "onset_time_reference_s": onset_ref,
        "onset_altitude_model_km": h_mod, "onset_altitude_reference_km": h_ref,
        "onset_altitude_diff_km": (h_mod - h_ref) if (h_mod is not None and h_ref is not None) else None,
        "demise_time_model_s": end_mod, "demise_time_reference_s": end_ref,
        "demise_time_diff_s": (end_mod - end_ref) if (end_mod is not None and end_ref is not None) else None,
        "demise_time_rel": ((end_mod - end_ref) / end_ref) if (end_mod is not None and end_ref is not None) else None,
        "final_mass_model_kg": float(m_mod[-1]), "sprayed_mass_kg": float(c["sprayed_mass_kg"][-1]), "film_mass_kg": float(c["film_mass_kg"][-1]),
    }


def plot_melt(history, outdir, title, reference=None, size_distribution_csv=None):
    """The Step 3 plots (MELT_PLOT_NAMES); with a melting reference the mass plots carry the SESAM overlay and a
    residual panel. Returns the paths."""
    import csv
    os.makedirs(outdir, exist_ok=True)
    apply_rcparams(plt)
    c = history.columns
    t, h, m = c["time_s"], c["altitude_km"], c["mass_kg"]
    paths = [os.path.join(outdir, n) for n in MELT_PLOT_NAMES]
    # 1. mass vs time (+ residual)
    if reference is not None:
        mask = (reference.time >= t[0]) & (reference.time <= t[-1])
        tr, mr = reference.time[mask], reference.mass[mask]
        fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
        _overlay(ax, tr, mr, t, m, "mass [kg]", title=title)
        rx.plot(tr, 100.0 * (np.interp(tr, t, m) - mr) / reference.mass[0], color=MODEL_COLOR, lw=1.0)
        rx.axhline(0.0, color=MUTED, lw=0.6)
        rx.set_ylabel("model - SESAM [% of m0]"); rx.set_xlabel("time [s]")
        strip_top_right_spines(rx)
    else:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(t, m, color=MODEL_COLOR, lw=1.2, label="model")
        ax.set_ylabel("mass [kg]"); ax.set_xlabel("time [s]"); ax.set_title(title, color=SECOND, fontsize=10)
        strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[0], dpi=150); plt.close(fig)
    # 2. mass vs altitude
    fig, ax = plt.subplots(figsize=(7, 5))
    if reference is not None:
        ax.plot(reference.mass, reference.altitude / 1e3, color=REF_COLOR, lw=1.6, label="SESAM")
    ax.plot(m, h, color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_xlabel("mass [kg]"); ax.set_ylabel("altitude [km]"); ax.set_title(title, color=SECOND, fontsize=10); ax.legend(frameon=False)
    strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[1], dpi=150); plt.close(fig)
    # 3. mass budget
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, c["sprayed_mass_kg"], color=MODEL_COLOR, lw=1.2, label="sprayed (cumulative)")
    ax.plot(t, c["removed_mass_kg"] - c["sprayed_mass_kg"], color=SECOND, lw=1.0, ls=":", label="removed instantly (cumulative)")
    ax.plot(t, c["film_mass_kg"], color=MUTED, lw=1.2, label="film")
    ax.plot(t, m, color=REF_COLOR, lw=1.2, ls="--", label="remaining (body + film)")
    ax.set_xlabel("time [s]"); ax.set_ylabel("mass [kg]"); ax.set_title(title, color=SECOND, fontsize=10); ax.legend(frameon=False)
    strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[2], dpi=150); plt.close(fig)
    # 4. theta_cr and spraying area
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax.plot(t, c["theta_cr_deg"], color=MODEL_COLOR, lw=1.0); ax.set_ylabel("theta_cr [deg]"); ax.set_title(title, color=SECOND, fontsize=10)
    bx.plot(t, c["spraying_area_m2"] * 1e4, color=MODEL_COLOR, lw=1.0); bx.set_ylabel("spraying area [cm2]"); bx.set_xlabel("time [s]")
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(paths[3], dpi=150); plt.close(fig)
    # 5. closure fractions and the two Knudsen numbers
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.stackplot(t, c["closure_fraction_girin"], c["closure_fraction_couette_slip"], c["closure_fraction_couette_fm"],
                 labels=["Girin closure (Kn_local < 0.01)", "Couette, slip-corrected shear", "Couette, free-molecular shear"],
                 colors=[INK, MODEL_COLOR, MUTED], alpha=0.8)
    ax.set_xlabel("time [s]"); ax.set_ylabel("windward area fraction"); ax.set_ylim(0.0, 1.0); ax.legend(frameon=False, loc="upper left", fontsize=8)
    bx2 = ax.twinx()
    bx2.semilogy(t, c["kn_body"], color=SECOND, lw=1.0, ls="--")
    bx2.semilogy(t, c["kn_local_stag"], color=SECOND, lw=1.0, ls=":")
    bx2.axhline(0.01, color=MUTED, lw=0.6)
    bx2.set_ylabel("Kn_body (dashed), Kn_local at the nose (dotted)", color=SECOND, fontsize=8)
    ax.set_title(title, color=SECOND, fontsize=10); strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[4], dpi=150); plt.close(fig)
    # 6. droplet size and film thickness
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax.semilogy(t, c["r_median_um"], color=MODEL_COLOR, lw=1.0, label="median r")
    ax.semilogy(t, c["r_max_um"], color=MUTED, lw=0.8, ls=":", label="max r")
    ax.set_ylabel("droplet radius [um]"); ax.legend(frameon=False); ax.set_title(title, color=SECOND, fontsize=10)
    bx.semilogy(t, np.maximum(c["film_thickness_mean_mm"], 1e-6), color=MODEL_COLOR, lw=1.0, label="mean film thickness")
    bx.semilogy(t, np.maximum(c["film_thickness_max_mm"], 1e-6), color=MUTED, lw=0.8, ls=":", label="max")
    bx.set_ylabel("film thickness [mm]"); bx.set_xlabel("time [s]"); bx.legend(frameon=False)
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(paths[5], dpi=150); plt.close(fig)
    # 7. size distributions (from size_distribution.csv when present next to the run, else skipped with an empty axes)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.5))
    path = size_distribution_csv or os.path.join(outdir, "size_distribution.csv")
    if os.path.isfile(path):
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        windows = sorted({(float(r["window_start_s"]), float(r["window_end_s"])) for r in rows})
        for lo, hi in windows:
            sel = [r for r in rows if float(r["window_start_s"]) == lo and float(r["window_end_s"]) == hi]
            r_mid = np.sqrt(np.array([float(r["r_lo_m"]) for r in sel]) * np.array([float(r["r_hi_m"]) for r in sel])) * 1e6
            dn, dm = np.array([float(r["dn"]) for r in sel]), np.array([float(r["dm_kg"]) for r in sel])
            total = (lo, hi) == windows[-1] and lo == 0.0
            kw = dict(color=INK, lw=1.6, label="flight") if total else dict(color=MODEL_COLOR, lw=0.7, alpha=0.5)
            if dn.sum() > 0:
                ax.loglog(r_mid, np.maximum(dn, 1e-300), **kw)
                bx.loglog(r_mid, np.maximum(dm, 1e-300), **kw)
        ax.set_ylim(bottom=max(ax.get_ylim()[0], 1e-1)); bx.set_ylim(bottom=max(bx.get_ylim()[0], 1e-12))
    ax.set_xlabel("droplet radius [um]"); ax.set_ylabel("dn per bin"); bx.set_xlabel("droplet radius [um]"); bx.set_ylabel("dM per bin [kg]")
    ax.set_title(title, color=SECOND, fontsize=10)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False)
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(paths[6], dpi=150); plt.close(fig)
    return paths
```


- [ ] **Step 5: Replace `reentry_model/viz.py`**

```python
"""Surface-temperature and cross-section animations and stills from a coupled run's VTK series (spec section 9).

Reads only exported files: `surface.pvd` + `surface_<k>.vtp` (nodal T, per-patch q_conv) and `field.pvd` +
`field_<k>.vtu` (nodal T on the volume mesh) written by coupled.py, and the run's history (time, altitude, velocity,
surface_T_max_K for the fixed colour scale). `animate` colours the sphere's surface; `animate_section` cuts the volume
field on the meridional plane z = 0 through the flight axis (stagnation point at +x, shown on the right) so the depth
of the heated layer is visible. PyVista renders off-screen; frames go to MP4 through imageio-ffmpeg, with a GIF
fallback when the MP4 writer is unavailable."""
import os
import xml.etree.ElementTree as ET

import numpy as np

CAMERA = [(0.19, -0.14, 0.11), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]     # three-quarter view of the windward (+x) face, for R = 0.05 m
WINDOW = (960, 720)
SECTION_NORMAL, SECTION_ZOOM = (0.0, 0.0, 1.0), 1.15                  # meridional plane through the flight axis; viewed down z
SCALAR_BAR = {"fmt": "%.0f", "vertical": True, "position_x": 0.86, "position_y": 0.15, "width": 0.05, "height": 0.6,
              "title_font_size": 14, "label_font_size": 12}


def read_series(run_dir, name="surface"):
    """[(time, path)] of a PVD collection written by coupled.py."""
    root = ET.parse(os.path.join(run_dir, name + ".pvd")).getroot()
    return [(float(d.get("timestep")), os.path.join(run_dir, d.get("file"))) for d in root.iter("DataSet")]


def render_frame(poly, clim, title, path=None, radius=0.05, melting=False):
    """One off-screen frame of the surface coloured by nodal T; with `melting` the patches releasing droplets are drawn
    on top, coloured by droplet radius (log scale, um). Returns the RGB array (and writes PNG if `path`)."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW)
    plotter.add_mesh(poly, scalars="T", cmap="inferno", clim=clim, smooth_shading=True,
                     scalar_bar_args={"title": "surface T [K]", **SCALAR_BAR})
    if melting and "r_droplet" in poly.cell_data:
        r = np.asarray(poly.cell_data["r_droplet"])
        emitting = np.isfinite(r) & (r > 0.0)
        if emitting.any():
            spots = poly.extract_cells(np.flatnonzero(emitting))
            spots.cell_data["log10 r [um]"] = np.log10(r[emitting] * 1e6)
            plotter.add_mesh(spots, scalars="log10 r [um]", cmap="viridis", clim=(1.0, 3.0),
                             scalar_bar_args={"title": "released droplets: log10 r [um]", "position_x": 0.05, "position_y": 0.05,
                                              "vertical": False, "width": 0.4, "height": 0.05, "fmt": "%.1f", "title_font_size": 12, "label_font_size": 11})
    plotter.add_text(title, position="upper_left", font_size=11)
    scale = radius / 0.05
    plotter.camera_position = [tuple(scale * c for c in CAMERA[0]), CAMERA[1], CAMERA[2]]
    image = plotter.screenshot(path, return_img=True)
    plotter.close()
    return image


def render_film_frame(poly, clim, title, path=None, radius=0.05):
    """One frame of the film thickness [um] per patch (log colour scale); clim is ignored (fixed 1 um - 10 mm)."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW)
    b = np.asarray(poly.cell_data["film_thickness"]) if "film_thickness" in poly.cell_data else np.zeros(poly.n_cells)
    poly.cell_data["film [um]"] = np.log10(np.maximum(b * 1e6, 1.0))
    plotter.add_mesh(poly, scalars="film [um]", cmap="Blues", clim=(0.0, 4.0), scalar_bar_args={"title": "film thickness: log10 [um]", **SCALAR_BAR})
    plotter.add_text(title, position="upper_left", font_size=11)
    scale = radius / 0.05
    plotter.camera_position = [tuple(scale * c for c in CAMERA[0]), CAMERA[1], CAMERA[2]]
    image = plotter.screenshot(path, return_img=True)
    plotter.close()
    return image


def section_of(grid):
    """The volume field cut on the plane z = 0 (through the flight axis), as a PolyData carrying nodal T."""
    return grid.slice(normal=SECTION_NORMAL, origin=(0.0, 0.0, 0.0))


def render_section_frame(grid, clim, title, path=None, radius=0.05, melting=False):
    """One off-screen frame of the meridional cross-section coloured by T (windward +x on the right); with `melting`
    the liquidus iso-line (liquid fraction 1) and the solidus iso-line (0) are drawn. Returns the RGB array."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    section = section_of(grid)
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW)
    plotter.add_mesh(section, scalars="T", cmap="inferno", clim=clim, scalar_bar_args={"title": "T [K]", **SCALAR_BAR})
    plotter.add_mesh(section.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                                   non_manifold_edges=False), color="black", line_width=1.5)
    if melting and "liquid_fraction" in section.point_data and section.n_points:
        f = np.asarray(section.point_data["liquid_fraction"])
        for level, colour in ((0.999, "white"), (0.001, "cyan")):
            if f.min() < level < f.max():
                iso = section.contour([level], scalars="liquid_fraction")
                if iso.n_points:
                    plotter.add_mesh(iso, color=colour, line_width=2.0)
    plotter.add_text(title + "   section z = 0, flow from the right" + ("   white: liquidus, cyan: solidus" if melting else ""),
                     position="upper_left", font_size=11)
    plotter.camera_position = "xy"
    plotter.camera.zoom(SECTION_ZOOM)
    image = plotter.screenshot(path, return_img=True)
    plotter.close()
    return image


def frame_title(t, history):
    c = history.columns
    return "t = {:.1f} s   h = {:.1f} km   V = {:.2f} km/s".format(
        t, np.interp(t, c["time_s"], c["altitude_km"]), np.interp(t, c["time_s"], c["velocity_kms"]))


def animate(run_dir, history, radius, fps=10, animation=True, stills=True, melting=False):
    """MP4 (GIF fallback) of the surface temperature over the run, plus stills at the start, peak heating, peak surface
    temperature and the end (with animation=False only the stills' frames are rendered); with `melting` the emitting
    patches are overlaid and the stills add melt onset, spraying onset and peak release. Returns {"animation", "stills"}."""
    render = (lambda poly, clim, title, path, radius: render_frame(poly, clim, title, path, radius, melting=True)) if melting else render_frame
    return _animate_series(run_dir, history, radius, "surface", render, "animation", "frames", "", fps, animation, stills, melting)


def animate_section(run_dir, history, radius, fps=10, animation=True, stills=True, melting=False):
    """The same for the meridional cross-section of the volume field: `section.mp4` (GIF fallback), `frames_section/`,
    and stills named `section_<label>_t<N>s.png` next to the surface ones (with `melting`: the liquidus iso-line)."""
    render = (lambda grid, clim, title, path, radius: render_section_frame(grid, clim, title, path, radius, melting=True)) if melting else render_section_frame
    return _animate_series(run_dir, history, radius, "field", render, "section", "frames_section", "section_", fps, animation, stills, melting)


def animate_film(run_dir, history, radius, fps=10, animation=True, stills=True):
    """The film-thickness video `film.mp4` (GIF fallback), `frames_film/` and stills `film_<label>_t<N>s.png`."""
    return _animate_series(run_dir, history, radius, "surface", render_film_frame, "film", "frames_film", "film_", fps, animation, stills, True)


def still_marks(history, melting=False):
    """{label: time} of the stills: start, peak heating, peak surface temperature, end, and when melting the melt
    onset, the spraying onset and the peak release rate."""
    c = history.columns
    marks = {"start": 0.0, "peak_heating": float(c["time_s"][int(np.argmax(c["convective_heat_W"]))]),
             "peak_surface_T": float(c["time_s"][int(np.argmax(c["surface_T_max_K"]))]), "end": float(c["time_s"][-1])}
    if melting and "sprayed_mass_kg" in c:
        melted = (c["removed_mass_kg"] + c["film_mass_kg"]) > 0.0
        sprayed = c["sprayed_mass_kg"] > 0.0
        if melted.any():
            marks["melt_onset"] = float(c["time_s"][melted][0])
        if sprayed.any():
            marks["spraying_onset"] = float(c["time_s"][sprayed][0])
            marks["peak_release"] = float(c["time_s"][int(np.argmax(c["released_mass_kg"]))])
    return marks


def _animate_series(run_dir, history, radius, series_name, render, movie_name, frames_subdir, stills_prefix, fps, animation, stills, melting=False):
    import imageio.v2 as imageio
    import pyvista as pv
    series = read_series(run_dir, series_name)
    c = history.columns
    clim = (float(c["temperature_K"][0]), float(c["surface_T_max_K"].max()))
    times = np.array([t for t, _ in series])
    marks = still_marks(history, melting)
    wanted = {int(np.argmin(np.abs(times - t))) for t in marks.values()} if stills else set()
    frames, out = {}, {"animation": None, "stills": []}
    frames_dir = os.path.join(run_dir, frames_subdir)
    os.makedirs(frames_dir, exist_ok=True)
    for k, (t, path) in enumerate(series):
        if animation or k in wanted:
            frames[k] = render(pv.read(path), clim, frame_title(t, history), os.path.join(frames_dir, "frame_{:04d}.png".format(k)), radius)
    if animation and frames:
        try:
            writer = imageio.get_writer(os.path.join(run_dir, movie_name + ".mp4"), fps=fps, codec="libx264", quality=7, macro_block_size=8)
            for k in sorted(frames):
                writer.append_data(frames[k])
            writer.close()
            out["animation"] = os.path.join(run_dir, movie_name + ".mp4")
        except Exception:                                   # no ffmpeg: GIF
            imageio.mimsave(os.path.join(run_dir, movie_name + ".gif"), [frames[k] for k in sorted(frames)], duration=1.0 / fps)
            out["animation"] = os.path.join(run_dir, movie_name + ".gif")
    if stills and series:
        stills_dir = os.path.join(run_dir, "stills")
        os.makedirs(stills_dir, exist_ok=True)
        for label, t in marks.items():
            k = int(np.argmin(np.abs(times - t)))
            path = os.path.join(stills_dir, "{}{}_t{:.0f}s.png".format(stills_prefix, label, times[k]))
            imageio.imwrite(path, frames[k])
            out["stills"].append(path)
    return out
```


- [ ] **Step 6: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py tests/test_reentry_model_viz.py tests/test_reentry_model_sesam_io.py -q`
Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add reentry_model/sesam_io.py reentry_model/compare.py reentry_model/viz.py tests/test_reentry_model_compare.py tests/test_reentry_model_viz.py
git commit -m "Add the mass-loss comparison, the melt plots and the melting overlays of the videos (Step 3 Task 12)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

