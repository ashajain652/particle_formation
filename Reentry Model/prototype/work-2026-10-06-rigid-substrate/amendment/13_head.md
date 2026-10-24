## Amendment of 2026-10-06 — `--rigid-substrate on|off`

> Part of the rigid-substrate amendment (sub-plan 09's amendment of this date holds the design; facts 78–87 in
> `00-shared-context.md` the measurements). The code below is the tested code, as a diff against the throwaway copy
> described there.

**The flag.** `--rigid-substrate on|off` in the Step 3 group, default `on`, sets `MeltSettings.rigid_substrate`. With
`off` the regime test reads the film plus the fully molten material, as in every run before this amendment, and the
amended code reproduces the copy before the change bit for bit on both flights (fact 81).

**Run names and the JSON.** A melting run with `--rigid-substrate off` ends in `_rigidsubstrate-off`, after
`_moltencascade-off` and before `_seed-<n>`; the default keeps the existing name, as the deep runoff's and the cascade's
switches did, so a default run made after this amendment carries the name of one made before it although its physics
differs — compare runs by their JSON settings, which now record `"rigid_substrate": "on"` or `"off"`. A bad value exits
2 like every other bad argument.

