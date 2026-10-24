"""viz.py: headless rendering of one frame and the animation/stills from a short coupled run's VTK series."""
import os

import numpy as np
import pytest

from reentry_model import viz

pytest.importorskip("pyvista")


def test_render_one_frame(coarse_sphere_mesh, tmp_path):
    import pyvista as pv
    surface = coarse_sphere_mesh.surface()
    faces = np.hstack([np.full((surface.n_patches, 1), 3), surface.faces]).ravel()
    poly = pv.PolyData(coarse_sphere_mesh.points, faces)
    poly.point_data["T"] = 300.0 + 500.0 * np.clip(coarse_sphere_mesh.points[:, 0] / 0.05, 0.0, 1.0)
    image = viz.render_frame(poly, (300.0, 800.0), "t = 1.0 s", str(tmp_path / "f.png"))
    assert image.shape == (720, 960, 3) and image.dtype == np.uint8 and os.path.getsize(tmp_path / "f.png") > 1000
    assert image.std() > 10.0                                            # not a blank frame


def test_render_section_frame(coarse_sphere_mesh, tmp_path):
    """The meridional cut through the flight axis of a volume field, coloured by T."""
    import pyvista as pv
    m = coarse_sphere_mesh
    cells = np.hstack([np.full((m.n_elements, 1), 4), m.tets]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(m.n_elements, pv.CellType.TETRA), m.points)
    grid.point_data["T"] = 300.0 + 500.0 * np.clip(m.points[:, 0] / 0.05, 0.0, 1.0)
    image = viz.render_section_frame(grid, (300.0, 800.0), "t = 1.0 s", str(tmp_path / "s.png"))
    assert image.shape == (720, 960, 3) and image.dtype == np.uint8 and os.path.getsize(tmp_path / "s.png") > 1000
    assert image.std() > 10.0
    section = viz.section_of(grid)
    assert section.n_points > 100 and float(np.abs(section.points[:, 2]).max()) < 1e-9        # lies in the z = 0 plane
    assert float(section["T"].max()) > 790.0 and float(section["T"].min()) < 310.0


def test_animation_from_a_short_run(coarse_sphere_mesh, tmp_path):
    from test_reentry_model_coupled import simulator, thermal_body
    from reentry_model import coupled, heating
    b = thermal_body(coarse_sphere_mesh)
    run_dir = str(tmp_path / "run")
    settings = coupled.CoupledSettings(dt=0.5, frames_every=10, output_dir=run_dir)
    hist = coupled.CoupledRun(simulator(b, t_max=15.0), b, heating.SesamEquivalentHeating(), settings).run()
    series = viz.read_series(run_dir)
    assert [t for t, _ in series] == [0.0, 5.0, 10.0, 15.0]
    out = viz.animate(run_dir, hist, 0.05, fps=5)
    assert os.path.isfile(out["animation"]) and os.path.getsize(out["animation"]) > 1000
    assert len(out["stills"]) == 4 and all(os.path.isfile(p) for p in out["stills"])
    assert sorted(os.listdir(os.path.join(run_dir, "frames"))) == ["frame_{:04d}.png".format(k) for k in range(4)]
    assert viz.frame_title(5.0, hist).startswith("t = 5.0 s   h = ")
    # the cross-section video from the same run's volume series
    assert [t for t, _ in viz.read_series(run_dir, "field")] == [0.0, 5.0, 10.0, 15.0]
    sec = viz.animate_section(run_dir, hist, 0.05, fps=5)
    assert os.path.basename(sec["animation"]).startswith("section.") and os.path.getsize(sec["animation"]) > 1000
    assert len(sec["stills"]) == 4 and all(os.path.isfile(p) and "section_" in os.path.basename(p) for p in sec["stills"])
    assert sorted(os.listdir(os.path.join(run_dir, "frames_section"))) == ["frame_{:04d}.png".format(k) for k in range(4)]
    stills_only = viz.animate_section(run_dir, hist, 0.05, animation=False)
    assert stills_only["animation"] is None and len(stills_only["stills"]) == 4
