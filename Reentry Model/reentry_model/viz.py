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


def render_frame(poly, clim, title, path=None, radius=0.05):
    """One off-screen frame of the surface coloured by nodal T; returns the RGB array (and writes PNG if `path`)."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW)
    plotter.add_mesh(poly, scalars="T", cmap="inferno", clim=clim, smooth_shading=True,
                     scalar_bar_args={"title": "surface T [K]", **SCALAR_BAR})
    plotter.add_text(title, position="upper_left", font_size=11)
    scale = radius / 0.05
    plotter.camera_position = [tuple(scale * c for c in CAMERA[0]), CAMERA[1], CAMERA[2]]
    image = plotter.screenshot(path, return_img=True)
    plotter.close()
    return image


def section_of(grid):
    """The volume field cut on the plane z = 0 (through the flight axis), as a PolyData carrying nodal T."""
    return grid.slice(normal=SECTION_NORMAL, origin=(0.0, 0.0, 0.0))


def render_section_frame(grid, clim, title, path=None, radius=0.05):
    """One off-screen frame of the meridional cross-section coloured by T (windward +x on the right); returns the RGB array."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    section = section_of(grid)
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW)
    plotter.add_mesh(section, scalars="T", cmap="inferno", clim=clim, scalar_bar_args={"title": "T [K]", **SCALAR_BAR})
    plotter.add_mesh(section.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                                   non_manifold_edges=False), color="black", line_width=1.5)
    plotter.add_text(title + "   section z = 0, flow from the right", position="upper_left", font_size=11)
    plotter.camera_position = "xy"
    plotter.camera.zoom(SECTION_ZOOM)
    image = plotter.screenshot(path, return_img=True)
    plotter.close()
    return image


def frame_title(t, history):
    c = history.columns
    return "t = {:.1f} s   h = {:.1f} km   V = {:.2f} km/s".format(
        t, np.interp(t, c["time_s"], c["altitude_km"]), np.interp(t, c["time_s"], c["velocity_kms"]))


def animate(run_dir, history, radius, fps=10, animation=True, stills=True):
    """MP4 (GIF fallback) of the surface temperature over the run, plus stills at the start, peak heating, peak surface
    temperature and the end (with animation=False only the stills' frames are rendered). Returns {"animation", "stills"}."""
    return _animate_series(run_dir, history, radius, "surface", render_frame, "animation", "frames", "", fps, animation, stills)


def animate_section(run_dir, history, radius, fps=10, animation=True, stills=True):
    """The same for the meridional cross-section of the volume field: `section.mp4` (GIF fallback), `frames_section/`,
    and stills named `section_<label>_t<N>s.png` next to the surface ones."""
    return _animate_series(run_dir, history, radius, "field", render_section_frame, "section", "frames_section", "section_",
                           fps, animation, stills)


def _animate_series(run_dir, history, radius, series_name, render, movie_name, frames_subdir, stills_prefix, fps, animation, stills):
    import imageio.v2 as imageio
    import pyvista as pv
    series = read_series(run_dir, series_name)
    c = history.columns
    clim = (float(c["temperature_K"][0]), float(c["surface_T_max_K"].max()))
    times = np.array([t for t, _ in series])
    marks = {"start": 0.0, "peak_heating": float(c["time_s"][int(np.argmax(c["convective_heat_W"]))]),
             "peak_surface_T": float(c["time_s"][int(np.argmax(c["surface_T_max_K"]))]), "end": float(c["time_s"][-1])}
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
