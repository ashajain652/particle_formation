"""Run cases of the reentry model and compare two result folders: the port against the prototype.

    same_results.py run --code <folder holding reentry_model> --out <dir> [--state <file>] --cases <list>
    same_results.py compare <dirA> <dirB> [--report]

``run`` refuses (exit 2, nothing run, nothing deleted) when a case's ``<out>/<case>`` folder already holds anything.
``run`` runs each case as ``python -B -m reentry_model <args> --outdir <out>/<case> --name <case>`` with the working
directory ``--code`` and PYTHONDONTWRITEBYTECODE=1, one after the other, and writes ``<out>/<case>.log`` (the first line
of ``pmset -g batt`` at the start and the end, the command's output, the exit code and the wall time).

``compare`` compares every file of either folder by relative path, except ``*.log *.pkl *.progress *.vtu *.vtp *.pvd``:
``.csv`` the same header and the same text in every cell except columns whose name starts with ``wall_`` or is ``runtime_s``; ``.json`` equal
throughout except the keys ``provenance`` and ``files``, any key starting with ``wall_`` or named ``runtime_s`` (the run's elapsed time), and ``settings.mesh_file`` and
``settings.material_file``, which are compared by file name and by the SHA-256 of the file they name; ``.npz`` the same
arrays, equal element for element (NaN equal to NaN); any other file byte for byte.  Exit 0 when identical; exit 2, with a message, when either folder does not exist or holds no compared
file (with
``--report``: always 0, and for every differing numeric column or field the largest absolute and relative difference
is printed).  Standard library and numpy only; importable without running anything.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

IGNORED_SUFFIXES = (".log", ".pkl", ".progress", ".vtu", ".vtp", ".pvd")
IGNORED_KEYS = ("provenance", "files")
FILE_KEYS = ("mesh_file", "material_file")      # under "settings": compared by name and content


# ------------------------------------------------------------------------------------------------ cases

COMMON = ("--diameter 100 --velocity 7.24 --altitude 65 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem "
          "--heating physics --h-surface 4 --h-core 20 --melt on --temperature 893 --dt 0.025 --t-max 0.125 "
          "--frames-every 0").split()
RESUME = ("--hand-down off --part-count-from melting --size-feedback current --spray-gates 2026-11-10 "
          "--leeward-test centroid --sideways-heat height --slip-band on --t-max 22.6 --frames-every 0").split()
NOMELT = ("--diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76 "
          "--thermal fem --heating physics --h-surface 4 --h-core 20 --t-max 20").split()


def case_arguments(state=None):
    """{case name: command-line arguments after ``-m reentry_model``, without --outdir and --name}."""
    return {
        "skins": ["run"] + COMMON + ["--skin-preset", "dev", "--fine-size", "2"],
        "element": ["run"] + COMMON + ["--surface-model", "elements"],
        "resume": ["resume", str(state)] + RESUME,
        "nomelt": ["run"] + NOMELT,
        # main's command line has no --surface-model, and the port's default for --thermal fem is the skins (which fail
        # without melting): this case is the nearest to main's non-melting run
        "nomelt_elements": ["run"] + NOMELT + ["--surface-model", "elements"],
    }


def power_line():
    out = subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout
    return out.splitlines()[0] if out else "(no pmset output)"


def run_case(code, out, name, args):
    """Run one case; returns its exit code.  The log ends with the exit code and the wall time."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, "-B", "-m", "reentry_model"] + args + ["--outdir", str(out / name), "--name", name]
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    start = time.time()
    with open(out / f"{name}.log", "w") as log:
        log.write(f"start power: {power_line()}\ncwd: {code}\ncommand: {' '.join(command)}\n")
        log.flush()
        code_returned = subprocess.run(command, cwd=code, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        log.write(f"\nend power: {power_line()}\nwall time: {time.time() - start:.1f} s\nexit {code_returned}\n")
    return code_returned


def run_cases(code, out, cases, state=None):
    table = case_arguments(state)
    for name in cases:                      # every refusal comes before anything runs; nothing is ever deleted
        if name not in table:
            raise SystemExit(f"unknown case {name!r}; cases are {sorted(table)}")
        if name == "resume" and state is None:
            raise SystemExit("the resume case needs --state")
        held = Path(out) / name
        if held.exists() and (not held.is_dir() or any(held.iterdir())):
            print(f"refused: {held} already holds results; choose another --out or remove it yourself")
            return 2
    status = 0
    for name in cases:
        returned = run_case(Path(code), out, name, table[name])
        print(f"{name}: exit {returned}", flush=True)
        status = status or returned
    return status


# ------------------------------------------------------------------------------------------------ comparison

class Difference:
    """One difference: where, how many values, and (if numeric) the largest absolute and relative difference."""

    def __init__(self, where, text, count=1, absolute=None, relative=None):
        self.where, self.text, self.count, self.absolute, self.relative = where, text, count, absolute, relative

    def __str__(self):
        line = f"{self.where}: {self.text}"
        if self.absolute is not None:
            line += f" (largest absolute difference {self.absolute:.6g}, largest relative difference {self.relative:.6g})"
        return line


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _numbers(a, b):
    """Both as floats, or None."""
    try:
        return float(a), float(b)
    except (TypeError, ValueError):
        return None


def _same_number(x, y):
    return x == y or (isinstance(x, float) and isinstance(y, float) and math.isnan(x) and math.isnan(y))


def _spread(pairs):
    """Largest absolute and relative difference of (x, y) float pairs (relative to the larger magnitude)."""
    absolute = relative = 0.0
    for x, y in pairs:
        if _same_number(x, y):
            continue
        d = abs(x - y)
        scale = max(abs(x), abs(y))
        absolute = max(absolute, d)
        relative = max(relative, d / scale if scale > 0 and math.isfinite(d) else (math.inf if d != 0 else 0.0))
    return absolute, relative


def _compare_csv(pa, pb, rel):
    with open(pa, newline="") as fa, open(pb, newline="") as fb:
        rows_a, rows_b = list(csv.reader(fa)), list(csv.reader(fb))
    if not rows_a or not rows_b or rows_a[0] != rows_b[0]:
        return [Difference(rel, "the headers differ")]
    header = rows_a[0]
    if len(rows_a) != len(rows_b):
        return [Difference(rel, f"{len(rows_a) - 1} rows against {len(rows_b) - 1} rows")]
    found = []
    for j, name in enumerate(header):
        if _is_wall(name):
            continue
        differing, pairs, numeric = 0, [], True
        for ra, rb in zip(rows_a[1:], rows_b[1:]):
            ca = ra[j] if j < len(ra) else None
            cb = rb[j] if j < len(rb) else None
            if ca == cb:
                continue
            differing += 1
            both = _numbers(ca, cb)
            if both is None:
                numeric = False
            else:
                pairs.append(both)
        if differing:
            absolute, relative = _spread(pairs) if numeric else (None, None)
            found.append(Difference(f"{rel} column {name}", f"{differing} of {len(rows_a) - 1} cells differ",
                                    differing, absolute, relative))
    return found


TIMER_KEYS = ("runtime_s",)     # the run's elapsed time (time.perf_counter in trajectory.py), besides wall_*


def _is_wall(key):
    return isinstance(key, str) and (key.startswith("wall_") or key in TIMER_KEYS)


def _named_file(value):
    """(file name, SHA-256 of the file or None if it is not there) of a path stored in a json file."""
    if not isinstance(value, str):
        return (value, None)
    path = Path(value)
    return (path.name, _sha256(path) if path.is_file() else None)


def _walk_json(a, b, where, found, in_settings=False):
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b), key=str):
            if key in IGNORED_KEYS or _is_wall(key):
                continue
            if key not in a or key not in b:
                found.append(Difference(f"{where}.{key}", "present in one folder only"))
            elif key == "state" and where.endswith(".resumed_from") and isinstance(a[key], str) \
                    and isinstance(b[key], str):
                if Path(a[key]).name != Path(b[key]).name:          # by file name only: the state is 1.4 GB, no checksum
                    found.append(Difference(f"{where}.{key}", f"file name differs ({Path(a[key]).name} against "
                                                              f"{Path(b[key]).name})"))
            elif in_settings and key in FILE_KEYS:
                fa, fb = _named_file(a[key]), _named_file(b[key])
                if fa[1] is None or fb[1] is None:
                    found.append(Difference(f"{where}.{key}", f"a named file is missing ({fa[0]}, {fb[0]})"))
                elif fa != fb:
                    found.append(Difference(f"{where}.{key}", f"file name or content differs ({fa[0]} {fa[1][:12]} "
                                                              f"against {fb[0]} {fb[1][:12]})"))
            else:
                _walk_json(a[key], b[key], f"{where}.{key}", found, in_settings=(key == "settings"))
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            found.append(Difference(where, f"lists of {len(a)} and {len(b)} items"))
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                _walk_json(x, y, f"{where}[{i}]", found)
    elif type(a) is not type(b) and not (isinstance(a, (int, float)) and isinstance(b, (int, float))
                                         and not isinstance(a, bool) and not isinstance(b, bool)):
        found.append(Difference(where, "different kinds of value"))
    elif isinstance(a, (int, float)) and not isinstance(a, bool):
        if not _same_number(float(a), float(b)) or a != b and not (math.isnan(a) and math.isnan(b)):
            absolute, relative = _spread([(float(a), float(b))])
            found.append(Difference(where, f"{a!r} against {b!r}", 1, absolute, relative))
    elif a != b:
        found.append(Difference(where, f"{a!r} against {b!r}"))


def _compare_json(pa, pb, rel):
    found = []
    _walk_json(json.loads(Path(pa).read_text()), json.loads(Path(pb).read_text()), rel, found)
    return found


def _compare_npz(pa, pb, rel):
    found = []
    with np.load(pa, allow_pickle=False) as za, np.load(pb, allow_pickle=False) as zb:
        for name in sorted(set(za.files) | set(zb.files)):
            if name not in za.files or name not in zb.files:
                found.append(Difference(f"{rel} array {name}", "present in one folder only"))
                continue
            x, y = za[name], zb[name]
            if x.shape != y.shape or x.dtype != y.dtype:
                found.append(Difference(f"{rel} array {name}", f"shape/dtype {x.shape} {x.dtype} against {y.shape} {y.dtype}"))
                continue
            if x.dtype.kind in "fc":
                same = (x == y) | (np.isnan(x) & np.isnan(y))
            else:
                same = x == y
            if not np.all(same):
                differing = int(np.count_nonzero(~same))
                absolute = relative = None
                if x.dtype.kind in "fiu":
                    xs, ys = x[~same].astype(float), y[~same].astype(float)
                    absolute, relative = _spread(zip(xs.ravel().tolist(), ys.ravel().tolist()))
                found.append(Difference(f"{rel} array {name}", f"{differing} of {x.size} elements differ",
                                        differing, absolute, relative))
    return found


def _compare_bytes(pa, pb, rel):
    return [] if _sha256(pa) == _sha256(pb) else [Difference(rel, "the files differ byte for byte")]


COMPARERS = {".csv": _compare_csv, ".json": _compare_json, ".npz": _compare_npz}


def _files(folder):
    folder = Path(folder)
    return {p.relative_to(folder).as_posix(): p for p in sorted(folder.rglob("*"))
            if p.is_file() and not p.name.endswith(IGNORED_SUFFIXES)}


def differences(dir_a, dir_b):
    """Every Difference between two result folders (an empty list when they are identical)."""
    files_a, files_b = _files(dir_a), _files(dir_b)
    found = []
    for rel in sorted(set(files_a) | set(files_b)):
        if rel not in files_a or rel not in files_b:
            found.append(Difference(rel, f"only in {'the second' if rel not in files_a else 'the first'} folder"))
            continue
        comparer = COMPARERS.get(Path(rel).suffix, _compare_bytes)
        found.extend(comparer(files_a[rel], files_b[rel], rel))
    return found


def compare(dir_a, dir_b):
    """The differences as printable lines; an empty list means identical."""
    return [str(d) for d in differences(dir_a, dir_b)]


# ------------------------------------------------------------------------------------------------ command line

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="run cases")
    run.add_argument("--code", required=True, help="the folder holding the reentry_model package")
    run.add_argument("--out", required=True)
    run.add_argument("--state", help="the saved state the resume case starts from")
    run.add_argument("--cases", required=True, help="comma-separated: skins, element, resume, nomelt, nomelt_elements")
    cmp_ = sub.add_parser("compare", help="compare two result folders")
    cmp_.add_argument("dir_a")
    cmp_.add_argument("dir_b")
    cmp_.add_argument("--report", action="store_true", help="list every difference with its size and exit 0")
    args = parser.parse_args(argv)

    if args.command == "run":
        out = os.path.abspath(args.out)         # the child runs in --code: paths are the caller's, made absolute here
        state = os.path.abspath(args.state) if args.state else None
        return run_cases(args.code, out, [c for c in args.cases.split(",") if c], state)
    for folder in (args.dir_a, args.dir_b):
        if not Path(folder).is_dir():
            print(f"refused: {folder} does not exist (or is not a folder)")
            return 2
        if not _files(folder):
            print(f"refused: {folder} holds no file to compare (a run that wrote nothing?)")
            return 2
    found = differences(args.dir_a, args.dir_b)
    if not found:
        files = len(_files(args.dir_a))
        print(f"identical: {files} files compared")
        return 0
    print(f"{len(found)} difference(s):")
    for d in found:
        print(f"  {d}")
    return 0 if args.report else 1


if __name__ == "__main__":
    sys.exit(main())
