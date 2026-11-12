"""The checks run after each stage of the skyscrapers port.  Each prints a report and exits 0 only on a pass.

    checks.py merge1 --ref <tree-ish> [--base <ref>]
    checks.py rewrite <original repo> <rewritten repo>
    checks.py outside <tree-ish> [--allow PATH]... [--base <ref>]
    checks.py same-code <build-ref> <port tree-ish> <main-ref>
    checks.py history <ref>

Git runs in the current working directory's repository, except ``rewrite``, which takes two repository paths.
``--base`` is the state of main that ``merge1`` and ``outside`` compare against: ``MAIN0``, main when the port began, unless
given, so the port's recorded commands reproduce; give a later main (for example ``--base HEAD`` inside a merge) to check
a change made after it.
"""
import argparse
import difflib
import re
import subprocess
import sys

from path_map import map_path, prototype_line_path
from resolve_both_sides import three_way

MERGE_BASE = "d16018bf6bb19944a88c19cf0987c33951ea3485"
MAIN0 = "8ed75cfeeda2364f14dba2bf73cc48f58af6e72b"
PROTOTYPE_END = "skyscrapers-prototype-end"
BUILD_END = "skyscrapers-build-end"
TRAILER = "Build-repository-commit"
GITIGNORE = "Reentry Model/.gitignore"


def git(*args, repo=None, check=True):
    cmd = ["git"] + (["-C", repo] if repo else []) + list(args)
    return subprocess.run(cmd, capture_output=True, check=check).stdout


def tree_id(ref, repo=None):
    """The tree id of any tree-ish: a bare tree id, a commit, a tag or a branch."""
    return git("rev-parse", "--verify", "--quiet", f"{ref}^{{tree}}", repo=repo).decode().strip()


def entries(ref, repo=None):
    """{path: (mode, blob id)} of the tree of ``ref``."""
    out = git("ls-tree", "-r", "-z", "--full-tree", tree_id(ref, repo), repo=repo)
    result = {}
    for entry in out.split(b"\0"):
        if entry:
            meta, path = entry.split(b"\t", 1)
            mode, _kind, blob = meta.decode().split()
            result[path.decode()] = (mode, blob)
    return result


def blob(entry, repo=None):
    return git("cat-file", "blob", entry[1], repo=repo)


def report(title, problems, notes=()):
    print(f"== {title}")
    for line in notes:
        print(line)
    for line in problems:
        print("FAIL: " + line)
    print("PASS" if not problems else f"FAILED ({len(problems)} problems)")
    return 0 if not problems else 1


def under(path, prefix):
    prefix = prefix.rstrip("/")
    return path == prefix or path.startswith(prefix + "/")


def check_merge1(ref, main_ref=MAIN0):
    tree = entries(ref)
    base, main = entries(MERGE_BASE), entries(main_ref)
    out = git("diff", "--name-status", "-M", "-z", MERGE_BASE, PROTOTYPE_END).decode().split("\0")
    changes, i = [], 0                     # (status letter, old path, new path)
    while i < len(out) and out[i]:
        status = out[i]
        if status[0] in "RC":
            changes.append((status[0], out[i + 1], out[i + 2]))
            i += 3
        else:
            changes.append((status[0], out[i + 1], out[i + 1]))
            i += 2
    problems, combined, same = [], [], 0
    for status, old, new in changes:
        if status == "D" or status == "R":
            if prototype_line_path(old) in tree:
                problems.append(f"deleted by the branch but present: {prototype_line_path(old)}")
        if status == "D":
            continue
        q = prototype_line_path(new)
        theirs = git("rev-parse", "--verify", "--quiet", f"{PROTOTYPE_END}:{new}").decode().strip()
        if q not in tree:
            problems.append(f"missing: {q}")
        elif tree[q][1] == theirs:
            same += 1
        elif main.get(prototype_line_path(old)) != base.get(old):
            combined.append(q)                       # main changed it since MERGE_BASE (old = new unless renamed)
        else:
            problems.append(f"differs from {PROTOTYPE_END}:{new} though main did not change it: {q}")
    notes = [f"paths changed by the branch: {len(changes)}; same as the branch's: {same}; combined with main's: {len(combined)}"]
    notes += ["  combined: " + q for q in sorted(combined)]
    return report("merge1", problems, notes)


def _commits(repo):
    order = git("rev-list", "--topo-order", "--reverse", "HEAD", repo=repo).decode().split()
    fmt = "%H%x00%P%x00%an%x00%ae%x00%ad%x00%cn%x00%ce%x00%cd%x00%(trailers:key=" + TRAILER + ",valueonly)"
    info = {}
    for c in order:
        fields = git("show", "-s", "--date=raw", f"--format={fmt}", c, repo=repo).decode().rstrip("\n").split("\0")
        message = git("show", "-s", "--format=%B", c, repo=repo).decode()
        info[c] = dict(parents=fields[1].split(), ident=fields[2:8], trailer=fields[8].strip(), message=message)
    return order, info


def check_rewrite(orig_repo, new_repo):
    problems = []
    o_order, o = _commits(orig_repo)
    n_order, n = _commits(new_repo)
    if len(o_order) != len(n_order):
        return report("rewrite", [f"{len(o_order)} original commits but {len(n_order)} rewritten"])
    pair = dict(zip(o_order, n_order))
    for oc, nc in zip(o_order, n_order):
        a, b = o[oc], n[nc]
        short = f"{oc[:10]} -> {nc[:10]}"
        if b["parents"] != [pair.get(p, "?") for p in a["parents"]]:
            problems.append(f"{short}: parents differ")
        if a["ident"] != b["ident"]:
            problems.append(f"{short}: author or committer name, email or date differs")
        if b["trailer"] != oc:
            problems.append(f"{short}: trailer {TRAILER} is {b['trailer']!r}, not the original identifier")
        stripped = re.sub(rf"^{TRAILER}: .*\n?", "", b["message"], flags=re.M)
        if stripped.rstrip() != a["message"].rstrip():
            problems.append(f"{short}: message differs after removing the trailer")
        want = {}
        for path, ent in entries(oc, orig_repo).items():
            q = map_path(path)
            if q in want:
                problems.append(f"{short}: two original paths map to {q}")
            want[q] = ent
        got = entries(nc, new_repo)
        for q in sorted(set(want) | set(got)):
            if q not in got:
                problems.append(f"{short}: dropped file {q}")
            elif q not in want:
                problems.append(f"{short}: extra file {q}")
            elif want[q] != got[q]:
                problems.append(f"{short}: mode or content differs: {q}")
    return report("rewrite", problems, [f"commits compared: {len(o_order)}"])


def check_outside(ref, allow, main_ref=MAIN0):
    changed = git("diff-tree", "-r", "--name-only", "--no-renames", "-z", tree_id(main_ref), tree_id(ref)).decode().split("\0")
    name = "MAIN0" if main_ref == MAIN0 else main_ref
    problems = []
    checked = 0
    for p in filter(None, changed):
        parts = p.split("/")[:-1]
        protected = not p.startswith("Reentry Model/") or "droplet_reentry" in parts or "wake_column" in parts
        if not protected:
            continue
        checked += 1
        if not any(under(p, a) for a in allow):
            problems.append(f"changed from {name} and not allowed: {p}")
    return report("outside", problems, [f"paths outside Reentry Model/ or under droplet_reentry/, wake_column/ that differ: {checked}"])


def check_same_code(build_ref, port_ref, main_ref):
    build, port, main = entries(build_ref), entries(port_ref), entries(main_ref)
    base = entries(MERGE_BASE)
    counts = {"same": 0, "main's change reapplied": 0, "main's copy kept": 0, "removed": 0, "fix-up": 0}
    problems, fixups = [], []
    for p, b in sorted(build.items()):
        q = map_path(p)
        if q not in port:
            if q == GITIGNORE:
                counts["removed"] += 1
            else:
                problems.append(f"missing from the port: {q} (build: {p})")
            continue
        if port[q] == b:
            counts["same"] += 1
            continue
        content = blob(port[q])
        if q in main and p in base:
            merged, clean = three_way(blob(main[q]), blob(base[p]), blob(b))
            if clean and merged == content:
                counts["main's change reapplied"] += 1
                continue
        if not q.startswith("Reentry Model/") and q in main and port[q] == main[q]:
            counts["main's copy kept"] += 1
            continue
        counts["fix-up"] += 1
        want = blob(b)
        try:
            diff = "".join(difflib.unified_diff(want.decode().splitlines(True), content.decode().splitlines(True),
                                                f"build:{p}", f"port:{q}"))
        except UnicodeDecodeError:
            diff = f"(binary files differ: build:{p}, port:{q})\n"
        fixups.append(f"--- fix-up: {q}\n{diff}")
    notes = [f"build files: {len(build)}"] + [f"  {k}: {v}" for k, v in counts.items()] + fixups
    return report("same-code", problems, notes)


def check_history(ref):
    path = "Reentry Model/reentry_model/skinbody.py"
    top = git("rev-parse", "--show-toplevel").decode().strip()
    blame = subprocess.run(["git", "-C", top, "blame", "--porcelain", ref, "--", path], capture_output=True)
    if blame.returncode != 0:
        return report("history", [f"git blame of {path} at {ref} failed: {blame.stderr.decode().strip()}"])
    out = blame.stdout.decode()
    shas = set()
    for line in out.splitlines():
        if re.match(r"^[0-9a-f]{40} \d+ \d+( \d+)?$", line):
            shas.add(line.split()[0])
    problems, build_commits = [], set()
    for sha in sorted(shas):
        t = git("show", "-s", f"--format=%(trailers:key={TRAILER},valueonly)", sha).decode().strip()
        if t:
            build_commits.add(t)
        else:
            problems.append(f"blame reaches {sha[:10]}, which has no {TRAILER} trailer")
    if subprocess.run(["git", "merge-base", "--is-ancestor", PROTOTYPE_END, ref], capture_output=True).returncode != 0:
        problems.append(f"{PROTOTYPE_END} is not an ancestor of {ref}")
    if subprocess.run(["git", "rev-parse", "--verify", "--quiet", f"refs/tags/{BUILD_END}"],
                      capture_output=True).returncode != 0:
        problems.append(f"tag {BUILD_END} does not exist")
    notes = [f"blame of {path}: {len(shas)} commits, {len(build_commits)} distinct build commits reached"]
    return report("history", problems, notes)


def main(argv):
    ap = argparse.ArgumentParser(prog="checks.py")
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("merge1")
    s.add_argument("--ref", required=True)
    s.add_argument("--base", default=MAIN0, help="the state of main to compare against (default: MAIN0)")
    s = sub.add_parser("rewrite")
    s.add_argument("original")
    s.add_argument("rewritten")
    s = sub.add_parser("outside")
    s.add_argument("ref")
    s.add_argument("--allow", action="append", default=[])
    s.add_argument("--base", default=MAIN0, help="the state of main to compare against (default: MAIN0)")
    s = sub.add_parser("same-code")
    s.add_argument("build_ref")
    s.add_argument("port_ref")
    s.add_argument("main_ref")
    s = sub.add_parser("history")
    s.add_argument("ref")
    a = ap.parse_args(argv[1:])
    if a.command == "merge1":
        return check_merge1(a.ref, a.base)
    if a.command == "rewrite":
        return check_rewrite(a.original, a.rewritten)
    if a.command == "outside":
        return check_outside(a.ref, a.allow, a.base)
    if a.command == "same-code":
        return check_same_code(a.build_ref, a.port_ref, a.main_ref)
    return check_history(a.ref)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
