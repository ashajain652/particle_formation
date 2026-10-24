## Amendment of 2026-10-07 (continuum step) — `--dt-continuum off|<s>`, by default 0.0125 s with Girin removal

> Asha's decision of 2026-10-07 on fact 96 (a) (facts 97–100 in `00-shared-context.md`; sub-plan 10's amendment of this
> date holds the switch). The code below is the tested code, as a diff against the copy of this date's earlier
> amendment.

**The flag.** `--dt-continuum off|<s>` in the Step 3 group. Its default depends on the removal mode, as
`--size-feedback`'s does: `DEFAULT_DT_CONTINUUM = 0.0125` s for melting with `--removal girin`, the model proper, where
the switched-step series converges (facts 92–94); off for `--removal instant`, the bookkeeping device, which has no film
and whose SESAM thresholds were measured at the default step. `off` keeps `--dt` throughout, as every run before this
amendment did. The switch's Knudsen number is `--kn-body-shock`, the surface flow's own continuum boundary, so moving
the gate moves the switch with it. `continuum_step(args)` resolves the value. A value of zero or below, a value longer
than `--dt`, anything else that is not a number, or `--dt-continuum` without `--melt on` exits 2.

**Run names and the JSON.** Only a value other than the mode's default changes the name: `_dtcontinuum-off` or
`_dtcontinuum-<s>`, after `_rigidsubstrate-off` and before `_seed-<n>`; the default keeps the existing name, as the
earlier switches' defaults did, so a default melting run made now carries the name one made before it carried although
it now switches. The JSON settings record `dt_continuum_s` (the value, or null) and the results the switch (sub-plan
10).

**Tests.** `test_the_continuum_step_defaults_to_0p0125_s_for_girin_removal_only`; the run-name test's four new cases;
three new bad-argument cases. `test_melting_run_writes_columns_files_and_json` (15 s from 71 km) crosses the continuum
boundary on the way, at about 9.5 s, so with the default it ran 460 macro steps instead of 30; it now passes
`--dt-continuum 0.1` and checks the switch it records, the steps on either side of it and the JSON, which keeps it at
about 30 s. Unit tier 258 passed, with Task 11's five known reference failures; FEniCSx tier 12 passed.

