"""Replay sub-plan 18's code and diff blocks on a fresh copy of the prototype and compare with the tested copy."""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.environ.get("SKIN_BASE") or ("/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/"
        "Particle Wake Evolution/prototype/work-2026-10-07-scheil-default/code")
plan = open(sys.argv[1]).read()
blocks = re.findall(r"```(python|diff)\n(.*?)\n```", plan, flags=re.S)
py = [b for k, b in blocks if k == "python"]
df = [b for k, b in blocks if k == "diff"]
assert len(py) == 13 and len(df) == 7, (len(py), len(df))
TA, SA, TB, SB, TC, SC, TD, SD, TT6, TF6, TMS, TF7, TCLI = py

dst = os.environ.get("SKIN_DST") or os.path.join(HERE, "replay", "code")
shutil.rmtree(os.path.dirname(dst), ignore_errors=True)
shutil.copytree(ORIG, dst, ignore=shutil.ignore_patterns("__pycache__", "reentry_model_output"))


def write(rel, text, mode="w"):
    with open(os.path.join(dst, rel), mode) as fh:
        fh.write(text)


def apply(d):
    out = subprocess.run(["patch", "-p1", "--no-backup-if-mismatch"], input=d + "\n", text=True, cwd=dst, capture_output=True)
    assert out.returncode == 0, out.stdout + out.stderr


# Tasks 2-5: the test file grows by appending; the module by appending methods to the class
write("tests/test_reentry_model_skin.py", TA + "\n\n\n" + TB + "\n\n\n" + TC + "\n\n\n" + TD + "\n")
write("reentry_model/skin.py", SA + "\n\n" + SB + "\n\n" + SC + "\n\n" + SD + "\n")
# Task 6
write("tests/test_reentry_model_thermal.py", "\n\n" + TT6 + "\n", "a")
write("tests/test_reentry_model_fenicsx.py", "\n\n" + TF6 + "\n", "a")
for d in df[:3]:
    apply(d)
# Task 7
write("tests/test_reentry_model_melting_skin.py", TMS + "\n")
write("tests/test_reentry_model_fenicsx.py", "\n\n" + TF7 + "\n", "a")
apply(df[3])
# Task 8
write("tests/test_reentry_model_cli.py", "\n\n" + TCLI + "\n", "a")
for d in df[4:]:
    apply(d)

out = subprocess.run(["diff", "-r", "-x", "__pycache__", "-x", "reentry_model_output", "-x", ".pytest_cache",
                      (os.environ.get("SKIN_CODE") or os.path.join(HERE, "code")), dst], capture_output=True, text=True)
print(out.stdout[:4000] or "replay identical to the tested copy")
