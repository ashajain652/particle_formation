"""Path map from the build repository's layout to main's (design section 4, step 5).

Standard library only: Task 4's filter-branch runs this file under ``python -I``.

    python path_map.py index-filter   # stdin: ``git ls-files -s`` lines; stdout: the same with mapped paths
"""
import sys

REENTRY = "Reentry Model/"
DRAMA = "DRAMA Sphere Experiments/"

# The sets below were checked against MAIN0's tree (8ed75cf) on 2026-11-11.
DRAMA_DATA_FILES = ("data/fap_day.dat", "data/fap_mon.dat")
DRAMA_DATA_DIRS = ("data/user_materials/",)
DRAMA_REFERENCE_DIR = "data/reference_runs/"
DRAMA_REFERENCE_MARK = "_mAA7075_nomelt"       # the non-melting references live with the DRAMA experiments
DRAMA_TEST_DIRS = ("tests/fixtures/",)
DRAMA_TEST_FILES = ("tests/test_fixtures.py", "tests/test_integration_drama.py",
                    "tests/test_sphere_reentry.py", "tests/test_sphere_sweep.py")
DRAMA_ANALYSIS = ("plot_knudsen_number.py", "plot_max_temperature_by_diameter.py", "plot_ohnesorge_number.py",
                  "plot_outcome_by_diameter.py", "plot_outcome_vs_diameter_velocity.py", "plot_style.py",
                  "plot_weber_number.py", "summary_from_runs.py")


def map_path(build_path: str) -> str:
    """A path in any commit of the build repository to its path on main (first matching rule wins)."""
    p = build_path
    if p.startswith("reentry_model/data/"):                                        # rule 1
        return REENTRY + "data/" + p[len("reentry_model/data/"):]
    if p in DRAMA_DATA_FILES or p.startswith(DRAMA_DATA_DIRS):                     # rule 2
        return DRAMA + p
    if p.startswith(DRAMA_REFERENCE_DIR):                                          # rule 3
        return (DRAMA if DRAMA_REFERENCE_MARK in p[len(DRAMA_REFERENCE_DIR):] else REENTRY) + p
    if p.startswith(DRAMA_TEST_DIRS) or p in DRAMA_TEST_FILES:                     # rule 4
        return DRAMA + p
    if p.startswith("analysis/") and p[len("analysis/"):] in DRAMA_ANALYSIS:       # rule 5
        return DRAMA + p
    return REENTRY + p                                                             # rule 6


def unmap_path(main_path: str):
    """The build path that maps to ``main_path``, or None if no build path does.

    Under ``Reentry Model/data/`` a path is taken to come from ``reentry_model/data/`` unless it is in
    ``data/reference_runs/`` or ``data/reference_values/`` (the build's own top-level data folders).
    """
    if main_path.startswith(DRAMA):
        candidate = main_path[len(DRAMA):]
    elif main_path.startswith(REENTRY):
        rest = main_path[len(REENTRY):]
        if rest.startswith("data/") and not rest.startswith((DRAMA_REFERENCE_DIR, "data/reference_values/")):
            candidate = "reentry_model/" + rest
        else:
            candidate = rest
    else:
        return None
    return candidate if map_path(candidate) == main_path else None


def prototype_line_path(step3_path: str) -> str:
    """A path of step3-rigid-substrate-scheil to its place after merge 1."""
    if step3_path == ".gitignore":
        return ".gitignore"
    return REENTRY + step3_path


def index_filter(lines):
    for line in lines:
        meta, tab, path = line.rstrip("\n").partition("\t")
        if not tab:
            raise SystemExit(f"path_map: not a git ls-files -s line: {line!r}")
        if path.startswith('"'):
            raise SystemExit(f"path_map: quoted path (unexpected characters): {path}")
        yield f"{meta}\t{map_path(path)}\n"


def main(argv):
    if argv[1:] != ["index-filter"]:
        raise SystemExit("usage: path_map.py index-filter < `git ls-files -s`")
    sys.stdout.writelines(index_filter(sys.stdin))


if __name__ == "__main__":
    main(sys.argv)
