"""Three-way resolver for the add/add conflicts of merge 2 (Task 5, Step 3).

    cd <worktree>; python resolve_both_sides.py     # inside a merge in progress

For each add/add conflict at path q: ``three_way(:2:q, MERGE_BASE:unmap_path(q), :3:q)``.  A clean result is
written and staged; the rest are left with conflict markers and listed.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from path_map import unmap_path

MERGE_BASE = "d16018bf6bb19944a88c19cf0987c33951ea3485"


def three_way(ours: bytes, base: bytes, theirs: bytes):
    """``git merge-file -p`` on the three; returns (merged bytes, whether it was clean)."""
    with tempfile.TemporaryDirectory() as tmp:
        names = []
        for label, data in (("ours", ours), ("base", base), ("theirs", theirs)):
            path = os.path.join(tmp, label)
            with open(path, "wb") as f:
                f.write(data)
            names.append(path)
        run = subprocess.run(["git", "merge-file", "-p", "-L", "ours", "-L", "base", "-L", "theirs", *names],
                             capture_output=True)
    if run.returncode < 0 or run.returncode > 127:      # 255: an error, e.g. binary files
        return ours, False
    return run.stdout, run.returncode == 0


def _git(*args, check=True):
    return subprocess.run(["git", *args], capture_output=True, check=check).stdout


def unmerged():
    """{path: {stage: (mode, blob)}} of the index's unmerged entries."""
    result = {}
    for entry in _git("ls-files", "-u", "-z").split(b"\0"):
        if entry:
            meta, path = entry.split(b"\t", 1)
            mode, blob, stage = meta.decode().split()
            result.setdefault(path.decode(), {})[int(stage)] = (mode, blob)
    return result


def main():
    toplevel = _git("rev-parse", "--show-toplevel").decode().strip()
    resolved, left = [], []
    for q, stages in sorted(unmerged().items()):
        if set(stages) != {2, 3}:
            left.append((q, "not an add/add conflict (stages %s)" % sorted(stages)))
            continue
        p = unmap_path(q)
        base = b""
        if p is not None:
            run = subprocess.run(["git", "cat-file", "blob", f"{MERGE_BASE}:{p}"], capture_output=True)
            if run.returncode == 0:
                base = run.stdout
        ours = _git("cat-file", "blob", stages[2][1])
        theirs = _git("cat-file", "blob", stages[3][1])
        merged, clean = three_way(ours, base, theirs)
        Path(toplevel, q).write_bytes(merged)
        if clean:
            _git("add", "--", q)
            resolved.append(q)
        else:
            left.append((q, "conflict markers left in the file"))
    print(f"resolved and staged: {len(resolved)}")
    for q in resolved:
        print("  " + q)
    print(f"left for hand resolution: {len(left)}")
    for q, why in left:
        print(f"  {q}: {why}")


if __name__ == "__main__":
    main()
