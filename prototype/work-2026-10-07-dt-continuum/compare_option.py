"""Throwaway check: the model's own --dt-continuum reproduces the 2026-10-07 harness runs bit for bit. Pairs: the default
100 mm run against the harness's 0.0125 s run, the --dt-continuum off run against the harness's 0.5 s run, and the 50 mm
flight on the linear range (which never switches) against the 2026-10-07 re-check of the 2026-10-06 run. Every history
column, every result field the older run has but the run time and the frame count, every source-table column; and, where
both wrote frames, their times."""
import csv, glob, json, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SERIES = os.path.join(HERE, "..", "work-2026-10-07-scheil-default")

def load_new(name):
    d = os.path.join(HERE, "runs")
    return (list(csv.reader(open(os.path.join(d, name + ".csv")))), json.load(open(os.path.join(d, name + ".json"))),
            np.load(os.path.join(d, name, "particles.npz")), os.path.join(d, name, "vtk"))

def load_old(run_dir, name=None):
    if name is None:
        name = [f for f in os.listdir(run_dir) if f.endswith(".csv")][0][:-4]
    return (list(csv.reader(open(os.path.join(run_dir, name + ".csv")))), json.load(open(os.path.join(run_dir, name + ".json"))),
            np.load(os.path.join(run_dir, name, "particles.npz")), os.path.join(run_dir, name, "vtk"))

def frame_times(vtk):
    pvd = os.path.join(vtk, "surface.pvd")
    return [float(x) for x in re.findall(r'timestep="([^"]+)"', open(pvd).read())] if os.path.exists(pvd) else None

def compare(label, a, b):
    (ra, ja, pa, va), (rb, jb, pb, vb) = a, b
    ha, hb = ra[0], rb[0]
    common = [c for c in hb if c in ha]
    ia, ib = [ha.index(c) for c in common], [hb.index(c) for c in common]
    same_len = len(ra) == len(rb)
    diff_cols = [c for c, i, k in zip(common, ia, ib) if not same_len or any(r[i] != s[k] for r, s in zip(ra[1:], rb[1:]))]
    res_a, res_b = ja["results"], jb["results"]
    diff_res = [k for k in res_b if k not in ("runtime_s", "n_frames") and res_a.get(k) != res_b[k]]
    keys = sorted(set(pa.files) & set(pb.files))
    diff_p = [k for k in keys if not np.array_equal(pa[k], pb[k], equal_nan=True)]
    fa, fb = frame_times(va), frame_times(vb)
    same = not diff_cols and not diff_res and not diff_p
    print("%-44s rows %d/%d, %d shared columns, %d source columns: %s" % (label, len(ra) - 1, len(rb) - 1, len(common), len(keys),
          "bit-identical" if same else "DIFFERENT columns %s results %s source %s" % (diff_cols[:6], diff_res[:6], diff_p[:6])))
    if fa is not None and fb is not None:
        print("    frames: %d/%d, times %s" % (len(fa), len(fb), "identical" if fa == fb else "DIFFER %s / %s" % (fa[:8], fb[:8])))
    print("    switch: %s s at %s km, Kn %s; dt_continuum %s" % (res_a.get("dt_switch_time_s"), res_a.get("dt_switch_altitude_km"),
          res_a.get("dt_switch_kn"), res_a.get("dt_continuum_s")))

pairs = {"default_d100": ("default (0.0125 s) vs harness dt0125", os.path.join(SERIES, "runs", "dt0125"), None),
         "off_d100": ("--dt-continuum off vs harness dt0500", os.path.join(SERIES, "runs", "dt0500"), None),
         "default_d050_range": ("50 mm range, default vs 2026-10-07 re-check", os.path.join(SERIES, "recheck"), "on_d050")}
for name in sys.argv[1:]:
    label, old_dir, old_name = pairs[name]
    compare(label, load_new(name), load_old(old_dir, old_name))
