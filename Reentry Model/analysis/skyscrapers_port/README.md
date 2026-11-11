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
