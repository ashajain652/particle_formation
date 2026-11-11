"""Write data/reference_values/box_q6p_fine.json from sub-plan 18's fine-box run (sub-plan 20 Task 9). Run from the
work copy's root: python ../extract_fine_box.py (or with the run file's path as the first argument)."""
import json
import os
import sys

RUN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "work-2026-10-08-skin", "analysis", "box", "runs",
                   "q6p__fine.json")
src = sys.argv[1] if len(sys.argv) > 1 else RUN
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join("data", "reference_values", "box_q6p_fine.json")
d = json.load(open(src))
doc = {"description": "The stage-1 fine box of sub-plan 18 (q6p): an 8 x 2 x 2 mm block of 0.1 mm elements at 700 K, AA7075_scheil, "
       "no radiation, 6 MW/m2 on its x = 0 face distributed by projection, liquid removed as it forms (element model, "
       "--removal instant), 0.005 s steps to 2 s. r_total_mm is the total removed mass over the 4 mm2 cross-section. "
       "Measured 2026-10-09 (prototype/work-2026-10-08-skin/analysis/box/runs/q6p__fine.json); a reference for the "
       "level-2 plate of the layered-skin design beside the one-dimensional 1 um reference (sub-plan 20 Task 9).",
       "spec": d["spec"], "t_s": d["series"]["t"], "r_total_mm": d["series"]["r_total_mm"]}
json.dump(doc, open(out, "w"), indent=0)
