## Amendment of 2026-10-07 — `AA7075_scheil` is the melting default

> Asha's decision of 2026-10-07 (facts 88–96 in `00-shared-context.md`). The code below is the tested code, as a diff
> against the copy of the rigid-substrate amendment of 2026-10-06, which that amendment's diff blocks rebuild exactly.

**The change.** `--melt on` without `--material` now selects `AA7075_scheil` — sub-plan 02's Scheil variant, with Step
4's latent heat of 390 kJ/kg and liquid surface tension of 0.80 N/m — instead of `AA7075_range`. The decision to give
the body one Scheil material was taken on 2026-09-27 (sub-plan 02, which named this one-line change as Task 13's); it
needs sub-plan 02's material code and its two new JSON files. `AA7075_range` and the three other packaged materials stay
selectable by name, and the help text now lists all five. A non-melting run keeps `AA7075_nomelt`. `MeltingBody`'s error
message for a material without a latent heat names the Scheil variant too.

**Run names.** Unchanged, and that is an open item (fact 96 (c)). `model_run_name` has never encoded the material, so a
default run made now carries the name a linear-range run carried before, and a run with `--material AA7075_range` the
name of a default run: in one output directory the second overwrites the first. Until the name carries the material,
compare runs by the JSON's `settings.material` and give runs of different materials their own `--name` or `--outdir`
(both analysis drivers already pass `--name`).

**Test.** `test_melting_run_writes_columns_files_and_json` expects `AA7075_scheil`, 390 kJ/kg and 0.80 N/m; it failed on
the material before the change, as it should, and passes after it. Unit tier: 251 passed, with only Task 11's five known
reference failures; FEniCSx tier 12 passed.

