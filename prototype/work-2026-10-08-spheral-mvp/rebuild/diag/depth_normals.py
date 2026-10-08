"""Diagnostic (2026-10-08): on a frame of the rerun, spheral_frag's layer_depths along the staircase facets' outward
normals versus along the exported n_derived (decision 9); and n_derived against the facet normal / file winding.
    cd <repo> && python prototype/rebuild/diag/depth_normals.py <vtk dir> k [k ...]"""
import os, sys
import numpy as np
import pyvista as pv
sys.path.insert(0, os.getcwd())
from spheral_frag import contract, frames, geometry

vtk = sys.argv[1]
for k in map(int, sys.argv[2:]):
    grid, poly = pv.read(os.path.join(vtk, "field_%d.vtu" % k)), pv.read(os.path.join(vtk, "surface_%d.vtp" % k))
    pts = np.asarray(grid.points, dtype=np.float64)
    tets = np.asarray(grid.cells_dict[pv.CellType.TETRA], dtype=np.int64)
    faces = np.asarray(poly.faces, dtype=np.int64).reshape(-1, 4)[:, 1:].copy()
    fr = contract.Frame(k=k, time_s=0.0, points=pts, tets=tets, faces=faces,
                        node={"T": np.asarray(grid.point_data["T"], float), "f_l": np.asarray(grid.point_data["liquid_fraction"], float)})
    nvec, _ = frames.face_vectors(pts, faces)
    n_file = nvec / np.linalg.norm(nvec, axis=1)[:, None]
    fo, flip = frames.orient_outward(fr)
    d = frames.derived_patch_arrays(fo)
    nd = np.asarray(poly.cell_data["n_derived"], float)
    ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", nd, d["normal"]), -1, 1)))
    print("k=%d: %d faces; file winding inward on %d (%.1f %%); n_derived . file-winding normal < 0 on %d; n_derived . outward facet normal < 0 on %d; "
          "angle n_derived-facet median %.1f, 90th %.1f, max %.1f deg" % (k, len(faces), len(flip), 100 * len(flip) / len(faces),
          int((np.einsum("ij,ij->i", nd, n_file) < 0).sum()), int((np.einsum("ij,ij->i", nd, d["normal"]) < 0).sum()),
          np.median(ang), np.percentile(ang, 90), ang.max()))
    loc = geometry.TetLocator(pts, tets)
    for label, nrm in (("facet normals", d["normal"]), ("n_derived", nd)):
        liq, slu, info = geometry.layer_depths(fo, loc, nrm, d["centroid"], return_info=True)
        print("   %-14s marched %d, left body %d (%.1f %%), wrong way %d, slurry depth median %.2f mm (of marched)" % (
            label, info["n_marched"], info["n_left_body"], 100 * info["n_left_body"] / max(info["n_marched"], 1), info["n_wrong_way"],
            1e3 * np.nanmedian(slu[slu > 0]) if (slu > 0).any() else 0.0))
