"""Bit-identity of two runs: every shared history column, every particles.npz array, every result field but the run time."""
import csv, glob, json, os, sys
import numpy as np
a, b = sys.argv[1], sys.argv[2]
(ja,), (jb,) = glob.glob(os.path.join(a, "*.json")), glob.glob(os.path.join(b, "*.json"))
assert os.path.basename(ja) == os.path.basename(jb), (ja, jb)
ra, rb = list(csv.DictReader(open(ja[:-5] + ".csv"))), list(csv.DictReader(open(jb[:-5] + ".csv")))
shared = [k for k in ra[0] if k in rb[0]]
extra = [k for k in rb[0] if k not in ra[0]]
assert len(ra) == len(rb), (len(ra), len(rb))
bad = [k for k in shared if any(x[k] != y[k] for x, y in zip(ra, rb))]
print("rows", len(ra), "shared columns", len(shared), "new columns", len(extra), "differing columns", bad)
pa, pb = np.load(os.path.join(ja[:-5], "particles.npz")), np.load(os.path.join(jb[:-5], "particles.npz"))
pbad = [k for k in pa.files if not np.array_equal(pa[k], pb[k], equal_nan=True)]
print("particle arrays", len(pa.files), "rows", len(pa[pa.files[0]]), "differing", pbad)
da, db = json.load(open(ja))["results"], json.load(open(jb))["results"]
rbad = [k for k in da if k in db and da[k] != db[k]]
print("result fields", len(da), "new", sorted(set(db) - set(da)), "differing", rbad)
