# Skyscrapers port tools

Tools for porting the skyscrapers code (the build repository, branch `skin-spraying`, and the prototype branch
`step3-rigid-substrate-scheil`) onto `main`. They run git in the current directory's repository, so `cd` to the port
worktree first. The tests (`test_port_tools.py`) are not collected by `pytest tests`; run them explicitly from
`Reentry Model`: `python -m pytest analysis/skyscrapers_port/test_port_tools.py`.

**`path_map.py`** holds the path map. `map_path(build_path)` sends a path in any commit of the build repository to its
path on `main` (first matching rule wins: `reentry_model/data/` to `Reentry Model/data/`, the DRAMA data, tests and
plot scripts to `DRAMA Sphere Experiments/`, the non-melting reference runs likewise, everything else under
`Reentry Model/`); `unmap_path(main_path)` is the inverse (None if no build path maps there);
`prototype_line_path(step3_path)` gives a path of the prototype branch its place after merge 1. It imports only the
standard library, because the history rewrite runs it under `python -I`:
`git ls-files -s | python -I path_map.py index-filter` writes the lines with each path mapped, for
`git update-index --index-info`.

**`resolve_both_sides.py`** resolves the add/add conflicts of a merge in progress: for each such path `q` it merges
`:2:q`, `MERGE_BASE:unmap_path(q)` and `:3:q` with `git merge-file`, stages a clean result and leaves the rest with
conflict markers, listing them. Command line: `python resolve_both_sides.py` (inside the merge). `three_way(ours, base,
theirs)` is the same merge on bytes, returning the merged bytes and whether it was clean.

**`checks.py`** holds the checks of each stage, each printing a report and exiting 0 only on a pass. A `<tree-ish>` may
be a commit, a branch, a tag or a bare tree id (for example `$(git write-tree)`).
`merge1 --ref <tree-ish>` checks that every path the prototype branch added, modified or renamed since `MERGE_BASE` is
present at its place after merge 1, equal to the branch's unless `main` also changed it (then listed as combined), and
that every path it deleted is absent. `rewrite <original repo> <rewritten repo>` checks that the rewritten history has
the same commits in the same order with the same parents, authors, dates and messages (plus the
`Build-repository-commit` trailer) and a tree equal to the original's with every path mapped. `outside <tree-ish>
[--allow PATH]...` fails on any difference from `MAIN0` outside `Reentry Model/`, or under `droplet_reentry/` or
`wake_column/`, that is not allowed. `same-code <build-ref> <port tree-ish> <main-ref>` classifies each file of the build
as same, main's change reapplied, main's copy kept, removed or fix-up (the fix-ups are printed as diffs), and fails only
on a file missing from the port. `history <ref>` checks that blame of `Reentry Model/reentry_model/skinbody.py` reaches
only commits carrying the trailer, that `skyscrapers-prototype-end` is an ancestor of the ref, and that the tag
`skyscrapers-build-end` exists.

**`same_results.py`** runs short cases of the model and compares two result folders, to show that a port (or a
refactor) computes what the code before it computed. `run --code <folder holding reentry_model> --out <dir> [--state
<file>] --cases <list>` runs each named case (`skins`, `element`, `resume`, `nomelt`, `nomelt_elements`; the first three are 0.125 s smoke
flights and the 22.6 s resume check of the runbook, the fourth is a 20 s non-melting flight; `nomelt_elements` is the same with `--surface-model elements`, which exists because main's command line has no `--surface-model` and the port's default for `--thermal fem` is the skins, which fail without melting) one after the other as
`python -B -m reentry_model <args> --outdir <out>/<case> --name <case>` in `--code`, and writes `<out>/<case>.log` with
the first line of `pmset -g batt` at the start and the end, the exit code and the wall time. `resume` needs `--state`.
`compare <dirA> <dirB> [--report]` compares every file of either folder by relative path (not `*.log`, `*.pkl`,
`*.progress`, `*.vtu`, `*.vtp`, `*.pvd`): `.csv` cell by cell as text except `wall_*` and `runtime_s` columns; `.json` equal except the
keys `provenance`, `files`, `runtime_s` (the run's elapsed time) and `wall_*`, with `settings.mesh_file` and `settings.material_file` compared by file name
and SHA-256 of the file; `.npz` element for element (NaN equal NaN); any other kind byte for byte. Exit 0 when identical;
with `--report` it lists each differing column or field with its largest absolute and relative difference and exits 0.
It uses only the standard library and numpy and can be imported without running anything. Tests: the last tests of
`test_port_tools.py`.
