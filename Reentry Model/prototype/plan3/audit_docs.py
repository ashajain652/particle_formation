"""Mechanical audit of the plan's three Task 15 documentation blocks against the prototype code.

    cd prototype/plan3 && "$PY" audit_docs.py

The blocks are prose and the generator cannot check them, so every claim that names something the
code defines -- a history column, a CLI flag, a source-table field -- is checked here, together with
a list of strings that later amendments superseded. Exit 0 clean, 1 with findings.
"""
import io, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROTO = os.environ.get("PROTO", os.path.join(HERE, os.pardir, "proto3"))
PLAN = os.path.join(HERE, os.pardir, os.pardir, "docs", "superpowers", "plans", "2026-09-20-melt-spraying.md")
rd = lambda *p: io.open(os.path.join(*p), encoding="utf-8").read()

coupled, spray, cli = (rd(PROTO, "reentry_model", f) for f in ("coupled.py", "spray.py", "cli.py"))
readme = rd(HERE, "readme_step3.md")
plan = rd(PLAN)
doc = plan[plan.index("## Physics model — `reentry_model` (Step 3"):]
fail = []

# 1. history columns: the README's Outputs list against coupled.MELT_COLUMNS
MELT = set(re.findall(r'"([^"]+)"', re.search(r"MELT_COLUMNS = \[(.*?)\]\n", coupled, re.S).group(1))) | {"mass_kg"}
named = set(re.findall(r'`([a-z][A-Za-z_0-9]*)`', readme))
SUFF = ("_kg", "_mm", "_m2", "_K", "_J", "_um", "_deg", "_Pa", "_ms", "_fraction",
        "_elements", "_released", "_active", "_factor", "_body", "_stag", "_shock", "_branch")
ghost = {n for n in named if n.endswith(SUFF)} - MELT - {"thick_mm"}   # thick_mm is SESAM's, not ours
if ghost:
    fail.append("README names history columns the code does not emit: %s" % sorted(ghost))
if MELT - named:
    fail.append("README's Outputs list omits columns: %s" % sorted(MELT - named))

# 2. CLI flags: documented vs defined. A flag named only to record that it was renamed away
#    (--lumped-mass) is a historical reference, not a claim that it exists.
actual = set(re.findall(r'add_argument\("(--[a-z-]+)"', cli))
renamed = set(re.findall(r'replacing `(--[a-z-]+)`', doc)) | set(re.findall(r'`(--[a-z-]+)` flag becomes', doc))
documented = set(re.findall(r'`(--[a-z-]+)', doc)) | set(re.findall(r'/(--[a-z-]+)`', doc))
if documented - actual - renamed:
    fail.append("documentation names CLI flags that do not exist: %s" % sorted(documented - actual - renamed))
melt_grp = cli[cli.index("me = r.add_argument_group"):cli.index('c = sub.add_parser("compare"')]
undocumented = set(re.findall(r'add_argument\("(--[a-z-]+)"', melt_grp)) - documented
if undocumented:
    fail.append("melt-group flags the documentation never mentions: %s" % sorted(undocumented))

# 3. source table
if re.search(r"source table: .*?, regime, branch", readme):
    fail.append("README calls spray.SOURCE_COLUMNS' `closure` field 'regime'")
if "closure" not in re.findall(r'"([^"]+)"', re.search(r"SOURCE_COLUMNS = \[(.*?)\]\n", spray, re.S).group(1)):
    fail.append("SOURCE_COLUMNS no longer has a closure column; the README describes one")

# 4. strings later amendments superseded. A dated amendment may keep a superseded clause only if it
#    points forward to the amendment that replaced it (the convention of amendment 20).
SUPERSEDED = {
    "regime_fraction": "column name replaced by closure_fraction_* (Task 10)",
    "+4.6 %": "droplet superheat re-measured as +2.4 % (50 mm) / +0.6 % (100 mm), fact 25",
    "reported-only Rayleigh": "the front-surface mode is applied since 2026-09-24, amendment 25",
    "inactive at 10–30 m/s": "the deceleration is 6–95 m/s2, fact 30",
    "0.5–0.7 s per macro step": "about 2 s per macro step, fact 23",
    "≈ 4 min on the default mesh": "2357 s for the surviving 100 mm flight, fact 24",
    "measured (146 s)": "same",
    "--feed-depth": "feed_depth is a MeltSettings field, not a CLI flag",
    "Regimes by Kn_δ": "superseded by the body-scale gate, amendment 15",
    "thick films: Girin 2017; thin films": "three regimes chosen in two stages, fact 36",
}
for s, why in SUPERSEDED.items():
    if s in doc:
        fail.append("superseded string %r still present (%s)" % (s, why))
for line in doc.split("\n"):
    if "Kn_δ" in line and "≥ 0.1" in line and "amendment 16" not in line:
        fail.append("Kn_delta regime claim with no supersession pointer: %s" % line.strip()[:100])

print("audit_docs: %d finding(s)" % len(fail))
for f in fail:
    print("  FAIL:", f)
if not fail:
    print("  clean — %d history columns and %d CLI flags consistent with the code, no superseded strings"
          % (len(MELT), len(actual)))
sys.exit(1 if fail else 0)
