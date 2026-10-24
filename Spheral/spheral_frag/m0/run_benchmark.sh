#!/usr/bin/env bash
# The timed benchmark matrix of Spheral M0 (plan Task 6, Asha 2026-10-05), run sequentially through the launcher.
#
#   spheral_frag/m0/run_benchmark.sh [OUTDIR]        # default spheral_output/m0/bench
#
# Each case below runs with damage off and on, 100 timed steps after 10 warm-up steps:
#
#   3D            dx 3.0, 2.2, 1.5, 1.1 mm (19 k - 394 k particles)   18 processes
#   RZ            dx 2.2, 1.1, 0.55, 0.275 mm (0.8 k - 52 k)           18 processes, fewer where a process would hold
#                                                                      under ~500 particles: n = min(18, N/500), so
#                                                                      1 at 2.2 mm and 6 at 1.1 mm (summary.json
#                                                                      records the count and the per-process spread)
#   3D scaling    dx 2.2 mm                                            1, 2, 4, 8, 12, 18 processes
#   RZ scaling    dx 0.55 mm                                           1, 2, 4, 8, 12, 18 processes
#
# 36 distinct runs: the 18-process point of each scaling series is the main matrix's run and is skipped the second
# time. The 3D 1.1 mm case drops to 30 timed steps if the 1.5 mm median step, scaled by the particle ratio, exceeds
# 30 s (the plan's rule, decided before the run rather than after it).
#
# Machine idle (plan, global constraints): before every run the host's 1-minute load average (sysctl vm.loadavg) must
# be under MAX_LOAD (2.0); the script waits for it, since an 18-process run leaves the load near 18 and it decays
# with a one-minute time constant (about two minutes to fall below 2). bench.py checks the Docker VM's own load the
# same way and exits 2 with "not starting" when it is too high; the script then waits and retries. The host load
# measured before launching is passed to bench.py and recorded in summary.json.
#
# Resumable: a run whose summary.json exists is skipped without waiting. Writes OUTDIR/<run>/summary.json and
# steps.csv (bench.py), OUTDIR/logs/<run>.log and OUTDIR/runs.tsv (run, exit status, wall seconds, host load at the
# start). Exit 0 when every run succeeded, 1 otherwise, 2 when the launcher cannot run (no Docker, no image).
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo="$(cd "$here/../.." && pwd -P)"
launcher="${SPHERAL:-$repo/spheral_frag/container/spheral}"
out="${1:-spheral_output/m0/bench}"
MAX_LOAD="${MAX_LOAD:-2.0}"
WARMUP=10
COMMIT=116c71f
RADIUS_MM=50

cd "$repo" || exit 2
mkdir -p "$out/logs"
case "$(cd "$out" && pwd -P)/" in "$repo"/*) ;; *) echo "run_benchmark: OUTDIR must be inside $repo" >&2; exit 2 ;; esac
"$launcher" -c "import Spheral" >/dev/null 2>&1 || { echo "run_benchmark: the launcher cannot import Spheral" >&2; exit 2; }

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

# bench GEOMETRY DX_MM NPROC DAMAGE STEPS
bench() {
    local geometry="$1" dx="$2" nproc="$3" damage="$4" steps="$5"
    local name status t0 t1 load
    name="$(printf "bench_%s_dx%.3fmm_R%smm_n%s_dmg%s_steps%s_wu%s_%s" \
            "$geometry" "$dx" "$RADIUS_MM" "$nproc" "$damage" "$steps" "$WARMUP" "$COMMIT")"
    if [[ -f "$out/$name/summary.json" ]]; then
        echo "== $name: exists, skipped"
        return
    fi
    echo "== $name"
    while true; do
        wait_for_idle
        load="$(host_load)"
        t0="$(now)"
        "$launcher" -n "$nproc" spheral_frag/m0/bench.py --geometry "$geometry" --dx-mm "$dx" --damage "$damage" \
            --steps "$steps" --warmup "$WARMUP" --radius-mm "$RADIUS_MM" --out "$out" --max-load "$MAX_LOAD" \
            --host-load "$load" > "$out/logs/$name.log" 2>&1
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
        grep "^bench:" "$out/logs/$name.log" | sed 's/^/   /'
    else
        failures=$((failures + 1))
        echo "   FAILED, exit $status (see $out/logs/$name.log)"
    fi
}

median_step() {   # median_step RUN_NAME_GLOB -> median step wall time in s, or empty
    local f
    f="$(ls "$out"/$1/summary.json 2>/dev/null | head -1)"
    [[ -n "$f" ]] && sed -n 's/.*"median": *\([0-9.eE+-]*\).*/\1/p' "$f" | head -1
}

rz_procs() {   # rz_procs DX_MM -> min(18, max(1, floor((pi R^2/2)/dx^2 / 500)))
    awk -v r="$RADIUS_MM" -v d="$1" 'BEGIN { n = int(3.141592653589793*r*r/2/(d*d)/500); if (n < 1) n = 1; if (n > 18) n = 18; print n }'
}

for damage in off on; do
    for dx in 3.0 2.2 1.5; do bench 3d "$dx" 18 "$damage" 100; done
    steps=100
    m="$(median_step "bench_3d_dx1.500mm_R${RADIUS_MM}mm_n18_dmg${damage}_steps100_wu${WARMUP}_${COMMIT}")"
    if [[ -n "$m" ]] && awk -v m="$m" 'BEGIN { exit !(m*393719/155331 > 30) }'; then
        steps=30
        echo "   3D 1.1 mm: the 1.5 mm median step $m s extrapolates above 30 s; 30 timed steps"
    fi
    bench 3d 1.1 18 "$damage" "$steps"
    for dx in 2.2 1.1 0.55 0.275; do bench rz "$dx" "$(rz_procs "$dx")" "$damage" 100; done
done
for damage in off on; do
    for n in 1 2 4 8 12 18; do bench 3d 2.2 "$n" "$damage" 100; done
    for n in 1 2 4 8 12 18; do bench rz 0.55 "$n" "$damage" 100; done
done

echo "run_benchmark: $failures failed; table in $runs"
[[ $failures -eq 0 ]]
