#!/usr/bin/env bash
# The Python-subclass comparison of Spheral M0 (plan Task 7, spec section 12 runner check 2), run sequentially
# through the launcher.
#
#   spheral_frag/m0/run_subclass.sh [OUTDIR]        # default spheral_output/m0/subclass
#
# 1. Method by method (one process, untimed): every method of each Python class against its built-in on synthetic
#    inputs, 3D and RZ -> OUTDIR/methods_3d.json, methods_rz.json.
# 2. Runs on Task 5's body (bench.py), damage on, 100 timed steps after 10 warm-up steps, 18 processes:
#    RZ at 1.1 mm (3.2 k particles) and 3D at 3.0 mm (19 k). Variants: builtin; eos, strength, damage and all (the
#    Python classes alone and together, each against the same builtin run); linpoly_builtin and linpoly_python (the
#    linear polynomial EOS, against each other); eos, strength and all again with --fma plain (unfused arithmetic,
#    the cost of a straightforward port and how far a few-ULP difference grows); builtin_repeat last (determinism of
#    the built-in and the drift of its step time over the session).
# 3. subclass_check.py compare -> OUTDIR/summary.json.
#
# Machine idle as run_benchmark.sh: before each run the host's 1-minute load average must be under MAX_LOAD (2.0);
# the script waits for it, passes it to the run (--host-load), and retries when the run itself finds the Docker VM
# busy. Resumable: a run whose summary.json exists is skipped. Exit 0 when every run succeeded, 1 otherwise, 2 when
# the launcher cannot run.
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo="$(cd "$here/../.." && pwd -P)"
launcher="${SPHERAL:-$repo/spheral_frag/container/spheral}"
out="${1:-spheral_output/m0/subclass}"
MAX_LOAD="${MAX_LOAD:-2.0}"
NPROC=18
STEPS=100
WARMUP=10
COMMIT=116c71f
RADIUS_MM=50

cd "$repo" || exit 2
mkdir -p "$out/logs"
case "$(cd "$out" && pwd -P)/" in "$repo"/*) ;; *) echo "run_subclass: OUTDIR must be inside $repo" >&2; exit 2 ;; esac
"$launcher" -c "import Spheral" >/dev/null 2>&1 || { echo "run_subclass: the launcher cannot import Spheral" >&2; exit 2; }

runs="$out/runs.tsv"
[[ -f "$runs" ]] || printf "run\texit\twall_s\thost_load_1m\n" > "$runs"
failures=0

host_load() { sysctl -n vm.loadavg | awk '{print $2, $3, $4}'; }
host_load_1m() { sysctl -n vm.loadavg | awk '{print $2}'; }
now() { perl -MTime::HiRes=time -e 'printf "%.3f\n", time'; }

wait_for_idle() {
    local waited=0
    until awk -v l="$(host_load_1m)" -v m="$MAX_LOAD" 'BEGIN { exit !(l < m) }'; do
        (( waited % 60 == 0 )) && echo "   waiting for the host load ($(host_load_1m)) to fall below $MAX_LOAD"
        sleep 10
        waited=$((waited + 10))
    done
}

for g in 3d rz; do
    if [[ ! -f "$out/methods_$g.json" ]]; then
        "$launcher" spheral_frag/m0/subclass_check.py methods --geometry "$g" --json "$out/methods_$g.json" \
            > "$out/logs/methods_$g.log" 2>&1 || { failures=$((failures + 1)); echo "== methods $g FAILED"; }
        grep "^METHODS" "$out/logs/methods_$g.log"
    fi
done

# run GEOMETRY DX_MM VARIANT FMA
run() {
    local geometry="$1" dx="$2" variant="$3" fma="$4"
    local name status t0 t1 load
    name="$(printf "subclass_%s_dx%.3fmm_R%smm_n%s_%s_fma%s_dmgon_steps%s_wu%s_%s" \
            "$geometry" "$dx" "$RADIUS_MM" "$NPROC" "$variant" "$fma" "$STEPS" "$WARMUP" "$COMMIT")"
    if [[ -f "$out/$name/summary.json" ]]; then
        echo "== $name: exists, skipped"
        return
    fi
    echo "== $name"
    while true; do
        wait_for_idle
        load="$(host_load)"
        t0="$(now)"
        "$launcher" -n "$NPROC" spheral_frag/m0/subclass_check.py run --geometry "$geometry" --dx-mm "$dx" \
            --variant "$variant" --fma "$fma" --steps "$STEPS" --warmup "$WARMUP" --radius-mm "$RADIUS_MM" \
            --out "$out" --max-load "$MAX_LOAD" --host-load "$load" > "$out/logs/$name.log" 2>&1
        status=$?
        t1="$(now)"
        if [[ $status -eq 2 ]] && grep -q "not starting" "$out/logs/$name.log"; then
            echo "   the Docker VM's load is above $MAX_LOAD; retrying in 30 s"
            sleep 30
            continue
        fi
        break
    done
    printf "%s\t%s\t%.1f\t%s\n" "$name" "$status" "$(echo "$t1 - $t0" | bc)" "${load%% *}" >> "$runs"
    if [[ $status -eq 0 ]]; then
        grep "^subclass_check:" "$out/logs/$name.log" | sed 's/^/   /'
    else
        failures=$((failures + 1))
        echo "   FAILED, exit $status (see $out/logs/$name.log)"
    fi
}

for case in "rz 1.1" "3d 3.0"; do
    set -- $case
    for v in builtin eos strength damage all linpoly_builtin linpoly_python; do run "$1" "$2" "$v" exact; done
    for v in eos strength all; do run "$1" "$2" "$v" plain; done
    run "$1" "$2" builtin_repeat exact
done

"$launcher" spheral_frag/m0/subclass_check.py compare --out "$out" 2>&1 | grep -v '^[|/\\]'
echo "run_subclass: $failures failed; table in $runs"
[[ $failures -eq 0 ]]
