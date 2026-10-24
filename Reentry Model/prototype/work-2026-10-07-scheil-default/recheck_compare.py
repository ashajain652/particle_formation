"""Throwaway check: does the freeze-back fix of 2026-10-07 change the flights of 2026-10-06? Each re-run of the fixed
copy against its run from the unfixed copy (and dt0500 of the series against scheil_on_d100), bit for bit: every
history column they share, every result field but the run time, every source-table column."""
import csv, json, os, sys
import numpy as np
OLD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "work-2026-10-06-rigid-substrate", "runs")
NEW = os.path.join(os.path.dirname(os.path.abspath(__file__)), "recheck")

def load(d, n, csv_name=None, json_name=None, npz_dir=None):
    rows = list(csv.reader(open(os.path.join(d, csv_name or n + ".csv"))))
    j = json.load(open(os.path.join(d, json_name or n + ".json")))
    p = np.load(os.path.join(d, npz_dir or n, "particles.npz"))
    return rows, j["results"], p

def compare(label, a, b):
    ra, ja, pa = a; rb, jb, pb = b
    ha, hb = ra[0], rb[0]
    common = [c for c in ha if c in hb]
    ia, ib = [ha.index(c) for c in common], [hb.index(c) for c in common]
    diff_cols = [c for c, i, k in zip(common, ia, ib) if len(ra) != len(rb) or any(r[i] != s[k] for r, s in zip(ra[1:], rb[1:]))]
    diff_res = [k for k in ja if k in jb and k not in ("runtime_s", "n_frames") and ja[k] != jb[k]]
    diff_p = [k for k in sorted(set(pa.files) & set(pb.files)) if not np.array_equal(pa[k], pb[k], equal_nan=True)]
    same = not diff_cols and not diff_res and not diff_p
    print("%-34s rows %d/%d, %d shared columns: %s%s" % (label, len(ra) - 1, len(rb) - 1, len(common),
          "bit-identical" if same else "DIFFERENT", "" if same else " columns %s results %s source %s" % (diff_cols[:8], diff_res[:8], diff_p[:8])))
    if not same and "sprayed_mass_kg" in ja:
        print("    sprayed %r vs %r; final mass %r vs %r" % (ja["sprayed_mass_kg"], jb["sprayed_mass_kg"], ja["final_mass_kg"], jb["final_mass_kg"]))

for n in sys.argv[1:]:
    if n == "dt0500":
        run = os.path.join(os.path.dirname(NEW), "runs", "dt0500")
        csvs = [f for f in os.listdir(run) if f.endswith(".csv")]
        sub = [f for f in os.listdir(run) if os.path.isdir(os.path.join(run, f))]
        name = csvs[0][:-4]
        compare("dt0500 (series) vs scheil_on_d100", load(run, name, npz_dir=name), load(OLD, "scheil_on_d100"))
    else:
        compare(n + " fixed vs 2026-10-06", load(NEW, n), load(OLD, n))
