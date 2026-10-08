"""Apply sub-plan 18's blocks to the dt-continuum copy, keeping going past failed hunks; report rejects."""
import os, re, shutil, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ["SKIN_BASE"]
dst = os.path.join(HERE, "rebase", "code")
plan = open(sys.argv[1]).read()
blocks = re.findall(r"```(python|diff)\n(.*?)\n```", plan, flags=re.S)
py = [b for k, b in blocks if k == "python"]
df = [b for k, b in blocks if k == "diff"]
TA, SA, TB, SB, TC, SC, TD, SD, TT6, TF6, TMS, TF7, TCLI = py
shutil.rmtree(os.path.dirname(dst), ignore_errors=True)
shutil.copytree(BASE, dst, ignore=shutil.ignore_patterns("__pycache__", "reentry_model_output", ".pytest_cache"))
def write(rel, text, mode="w"):
    with open(os.path.join(dst, rel), mode) as fh:
        fh.write(text)
def apply(d, name):
    out = subprocess.run(["patch", "-p1", "--no-backup-if-mismatch"], input=d + "\n", text=True, cwd=dst, capture_output=True)
    print(name, "rc", out.returncode, out.stdout.strip().replace("\n", " | ")[:600])
write("tests/test_reentry_model_skin.py", TA + "\n\n\n" + TB + "\n\n\n" + TC + "\n\n\n" + TD + "\n")
write("reentry_model/skin.py", SA + "\n\n" + SB + "\n\n" + SC + "\n\n" + SD + "\n")
write("tests/test_reentry_model_thermal.py", "\n\n" + TT6 + "\n", "a")
write("tests/test_reentry_model_fenicsx.py", "\n\n" + TF6 + "\n", "a")
for d, n in zip(df[:3], ("init", "skfem", "fenicsx")):
    apply(d, n)
write("tests/test_reentry_model_melting_skin.py", TMS + "\n")
write("tests/test_reentry_model_fenicsx.py", "\n\n" + TF7 + "\n", "a")
apply(df[3], "body7")
write("tests/test_reentry_model_cli.py", "\n\n" + TCLI + "\n", "a")
for d, n in zip(df[4:], ("cli", "coupled", "body8")):
    apply(d, n)
