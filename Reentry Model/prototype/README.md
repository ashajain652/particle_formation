# Step 3 prototype and plan generator

Moved here from the session scratchpad on 2026-09-26, which had pruned it twice (ten prototype files
on 2026-09-25, and `make_plan.py` itself mid-regeneration on 2026-09-26). The prototype is a throwaway copy
of the package by design, and the plan it produces is the committed artefact; its code, the generator and the
work copy of each amendment (`work-<date>-<name>/`) are version-controlled since 2026-10-08, their run outputs
are not (`.gitignore`).

- `proto3/` — the executed prototype. A complete, runnable copy of the package with Step 3 applied,
  plus its tests, `analysis/` scripts and `data/`. The prototype run outputs
  (`reentry_model_output/`, 17 MB) were not carried over.
- `work-2026-10-08-spheral-mvp/` — the prototype as reconstructed on the Spheral machine (2026-10-07/08) plus
  the export additions Spheral M1 reads (wall loads on every step, `v_hat_body`, the derived surface's normals,
  loads on newly exposed faces, the material in the run name). It wrote the 100 mm Scheil frames that Spheral M1
  is built on as an MVP input; it is not the Step 3 line of record (its README says what it lacks).
- `work-2026-10-08-frame-export/` — the frame-export amendment (facts 101–107): `code/` is the continuum-step copy
  (`work-2026-10-07-dt-continuum/code`) with the reconstruction's export additions 02, 03 (the material now named
  always), 04 and 06 ported onto it, `amendment/` the five diffs that rebuild `code/` from that copy, `meas/`
  the frozen measurement copy (`meas_manifest.md5`) and `launch.sh` the 100 mm Scheil frames flight it ran.
- `plan3/` — the generator. `make_plan.py` assembles the plan from the prose parts
  (`part_header.py`, `part_tasks_*.py`, `readme_step3.md`, `assumptions_step3.md`) and reads every
  code block **straight out of `proto3/`**, so a regeneration that matches the committed plan proves
  the plan's code blocks are the code that was tested.

## The invariant

    cd prototype/plan3
    "$PY" make_plan.py
    diff plan.md ../../docs/superpowers/plans/2026-09-20-melt-spraying.md

must report no difference, where `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`. Byte
identity is what the workflow rule means by "the plan is generated from the prototype". It proves the
code blocks are current; it does **not** prove the prose is, which is why the three mechanical audits
in the workflow notes exist.

`plan3/audit_docs.py` is those audits as a script, for the one part of the plan the invariant cannot
reach: the three Task 15 documentation blocks are prose, so nothing stops them describing a model
that has moved on. It checks the README's output-column list against `coupled.MELT_COLUMNS`, every
`--flag` the documentation names against `cli.py`, the source-table field names against
`spray.SOURCE_COLUMNS`, and a list of strings that later amendments superseded. Run it after any
amendment, and add the superseded wording to `SUPERSEDED` when an amendment retracts something:

    cd prototype/plan3
    "$PY" audit_docs.py

A dated amendment may keep a clause a later amendment overturned, but only with a forward pointer to
the amendment that overturned it — the convention amendment 20 set and amendment 5 now follows.

To change a code block, edit the file in `proto3/`, run its tests there, regenerate, then copy
`plan.md` over the committed plan. To change prose, edit the part file and regenerate. Never edit the
committed plan by hand — the next regeneration would silently revert it.

## Running the prototype

    cd prototype/proto3
    "$PY" -m pytest -m "not drama and not reference" -q \
        --ignore=tests/test_integration_drama.py \
        --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py

Those three modules need the DRAMA wrapper scripts at the real repo root and are out of scope here.
Expect 201 passed, with 2 failures and 3 errors that all come from one cause: the two melting SESAM
reference runs (plan Task 11) are not in `proto3/data/reference_runs/`. The repo's own `pytest.ini`
sets `testpaths = tests`, so a test run at the repo root does not collect anything under
`prototype/`.

`PROTO` overrides which prototype the generator reads, if you ever need to point it elsewhere.

## If files go missing again

Every prototype file the generator reads is present in the committed plan as a fenced block, so the
tree is recoverable: patch `read()`/`tail_from()` in `make_plan.py` to return a sentinel line
`@@@MISSING::<path>@@@` when the file is absent, regenerate, then walk the fenced blocks of the
sentinel output and the committed plan in parallel — the two have the same block count — and write
each sentinel block's counterpart back to disk. A file the generator reads whole comes back exactly;
one read with `tail_from` comes back as its tail, which you append to the Step 2 version from the
repo. Regenerating to byte identity then confirms the reconstruction.
