## Amendment of 2026-10-07 (continuum step) — the continuum step in the README, the assumptions and the spec

> Part of the continuum-step amendment (facts 97–100 in `00-shared-context.md`; sub-plans 10 and 13 hold the code).
> These changes stack on the blocks above, which still apply except where a change below replaces their wording. The
> blocks are prose, so re-run `prototype/plan3/audit_docs.py`'s checks against them after these edits: `--dt-continuum`
> must exist in `cli.py`; no history column is added (the switch is recorded in the run JSON's results); and add to its
> `SUPERSEDED` list the clause change 3 below retracts from this date's earlier change 6, "Until the fine step is a
> model option, quote the droplet population only from a fine-step run, with its step".

**README block (`## Physics model — reentry_model (Step 3 ...)`).**

1. *Options paragraph.* After the `--rigid-substrate on|off (...)` entry insert: "`--dt-continuum off|<s>` (default
   0.0125 s with `--removal girin`, off with `instant`: from the first macro step whose body Knudsen number is below
   `--kn-body-shock`, where Girin's closure begins, the step shortens to this, latched, because the droplet population
   converges from it and not at the default step; `off` keeps `--dt` throughout, as every run before 2026-10-07 did; a
   default 100 mm flight to 120 s took 2.3 hours, against 11 minutes with `off`, measured);".
2. *Outputs.* Add to the run JSON's results `dt_continuum_s`, `dt_switch_time_s`, `dt_switch_kn` and
   `dt_switch_altitude_km`, and note that history rows stay one per macro step, so a run that switches has many more of
   them after the switch, and that frames are written by flight time after it.
3. *Findings, the time-step bullet of this date's change 6.* Replace its last sentence, "Until the fine step is a model
   option, quote the droplet population only from a fine-step run, with its step (plan facts 88–96)." with: "Since
   2026-10-07 the model switches to 0.0125 s by default at the continuum boundary (`--dt-continuum`), so a default run's
   droplet population is the converged one; a run with `--dt-continuum off` gives the masses within about 5 % and the
   droplet population six times too small (plan facts 88–100)."

**Assumptions block (`## 9. Melting and the melt film (Step 3)`).**

4. In the bullet "The droplet population and the macro step" of this date's change 10, append: "Since 2026-10-07 the
   model takes the converged step by default: from the first macro step whose body Knudsen number is below the continuum
   boundary (0.01), 0.0125 s, latched (`--dt-continuum`; off for the instant-removal device). The default 0.5 s step
   serves before it, where the flight is in the merged regime and only the thin mode and the Rayleigh–Taylor mode act;
   whether that part of the flight depends on the step was not measured."

**Spec amendments block (`## 18. Amendments`).** Append:

5. "33. §7, §13 (decided 2026-10-07) — **the continuum step.** From the first macro step whose body Knudsen number is
   below the continuum boundary (`--kn-body-shock`, 0.01) the macro step is 0.0125 s, latched (`--dt-continuum`, default
   with Girin removal, off with instant removal and `off` on request); frames are written by flight time after the
   switch, and the run's results record it. It reproduces the switched-step harness of amendment 32's series bit for
   bit, so that series is now a set of ordinary runs; the 50 mm flight never switches and is unchanged; a 100 mm flight
   to 120 s cost 2.3 hours against 11 minutes (plan facts 97–100)."

