"""Fill plan_template.md's code and diff blocks from the tested scratch files: python make_plan18.py OUT.md"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.environ.get("SKIN_CODE") or os.path.join(HERE, "code")
SNAP7 = os.path.join(HERE, "snap", "after7")
ORIG = os.environ.get("SKIN_BASE") or ("/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/"
        "Particle Wake Evolution/prototype/work-2026-10-07-scheil-default/code")


def read(*p):
    return open(os.path.join(*p)).read()


def between(text, start, end=None):
    i = text.index(start) if start else 0
    j = text.index(end, i) if end else len(text)
    return text[i:j].strip("\n")


def diff(old, new, rel):
    out = subprocess.run(["diff", "-u", "--label", "a/" + rel, "--label", "b/" + rel, old, new],
                         capture_output=True, text=True)
    assert out.returncode == 1, (rel, out.returncode, out.stderr)
    return out.stdout.rstrip("\n")


skin_py = read(CODE, "reentry_model", "skin.py")
t_skin = read(CODE, "tests", "test_reentry_model_skin.py")
t_thermal = read(SNAP7, "tests", "test_reentry_model_thermal.py")
t_fx = read(SNAP7, "tests", "test_reentry_model_fenicsx.py")

blocks = {
    "SKIN_A": between(skin_py, None, "    # -- the sub-step"),
    "SKIN_B": between(skin_py, "    # -- the sub-step", "    # -- feed and draw"),
    "SKIN_C": between(skin_py, "    # -- feed and draw", "    # -- the macro step"),
    "SKIN_D": between(skin_py, "    # -- the macro step"),
    "TEST_SKIN_A": between(t_skin, None, "# -- the implicit sub-step (Task 3)"),
    "TEST_SKIN_B": between(t_skin, "# -- the implicit sub-step (Task 3)", "# -- feed and draw (Task 4)"),
    "TEST_SKIN_C": between(t_skin, "# -- feed and draw (Task 4)", "# -- the macro step's passes (Task 5)"),
    "TEST_SKIN_D": between(t_skin, "# -- the macro step's passes (Task 5)"),
    "TEST_THERMAL_T6": between(t_thermal, "def _interface_solver"),
    "TEST_FENICSX_T6": between(t_fx, "def test_backends_agree_with_an_interface_flux", "def test_skins_match_the_skfem_backend"),
    "TEST_FENICSX_T7": between(t_fx, "def test_skins_match_the_skfem_backend"),
    "TEST_MELTING_SKIN": read(CODE, "tests", "test_reentry_model_melting_skin.py").strip("\n"),
    "TEST_CLI_T8": read(HERE, "task8_cli_tests.py").strip("\n"),
    "DIFF_THERMAL_INIT": diff(os.path.join(ORIG, "reentry_model/thermal/__init__.py"),
                              os.path.join(SNAP7, "reentry_model/thermal/__init__.py"), "reentry_model/thermal/__init__.py"),
    "DIFF_SKFEM": diff(os.path.join(ORIG, "reentry_model/thermal/skfem_backend.py"),
                       os.path.join(SNAP7, "reentry_model/thermal/skfem_backend.py"), "reentry_model/thermal/skfem_backend.py"),
    "DIFF_FENICSX": diff(os.path.join(ORIG, "reentry_model/thermal/fenicsx_backend.py"),
                         os.path.join(SNAP7, "reentry_model/thermal/fenicsx_backend.py"), "reentry_model/thermal/fenicsx_backend.py"),
    "DIFF_BODY_T7": diff(os.path.join(ORIG, "reentry_model/body.py"), os.path.join(SNAP7, "body.py"), "reentry_model/body.py"),
    "DIFF_BODY_T8": diff(os.path.join(SNAP7, "body.py"), os.path.join(CODE, "reentry_model/body.py"), "reentry_model/body.py"),
    "DIFF_CLI": diff(os.path.join(ORIG, "reentry_model/cli.py"), os.path.join(CODE, "reentry_model/cli.py"), "reentry_model/cli.py"),
    "DIFF_COUPLED": diff(os.path.join(ORIG, "reentry_model/coupled.py"), os.path.join(CODE, "reentry_model/coupled.py"),
                         "reentry_model/coupled.py"),
    "UNIT_BEFORE": "258 passed, 1 skipped, 2 failed, 3 errors — the five failures are the known ones that need Task 11's "
                   "melting SESAM references",
    "UNIT_AFTER": open(os.path.join(HERE, "unit_after.txt")).read().strip(),
    "UNIT_REBASE": open(os.path.join(HERE, "unit_rebase.txt")).read().strip(),
    "BITID": open(os.path.join(HERE, "bitid_result.txt")).read().strip(),
    "CLI_EXPECTED": "7 passed in about 70 s (the short skin run takes most of it)",
}

# the tested test files must be exactly the pieces the plan hands out
assert "\n\n\n".join([blocks["TEST_SKIN_A"], blocks["TEST_SKIN_B"], blocks["TEST_SKIN_C"], blocks["TEST_SKIN_D"]]) + "\n" == t_skin
assert "\n\n\n".join([blocks["SKIN_A"], blocks["SKIN_B"], blocks["SKIN_C"], blocks["SKIN_D"]]).replace("\n\n\n    #", "\n\n    #") + "\n" == skin_py

plan = read(HERE, "plan_template.md")
for key, value in blocks.items():
    tag = "{{" + key + "}}"
    assert tag in plan, key
    plan = plan.replace(tag, value)
assert "{{" not in plan
open(sys.argv[1], "w").write(plan)
print("wrote", sys.argv[1], len(plan.splitlines()), "lines")
