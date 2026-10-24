#!/usr/bin/env bash
# LLNL's regression examples on the pinned arm64 image (Spheral M0 plan, Task 3; spec §12 check 1).
#
#   spheral_frag/m0/run_regressions.sh [OUTDIR]        # default spheral_output/regressions
#
# Runs, through the launcher, exactly the #ATS lines of the two upstream examples that need neither a GPU nor RAJA:
#
#   TensileRod-1d.py   t10-t13 (GradyKippTensorDamageOwen) and t20-t23 (ProbabilisticDamageModel): serial run with
#                      --checkRef True against LLNL's stored reference (filearraycmp, rtol = atol = 1e-4), the same rod
#                      on 4 processes, and restarts from cycle 500 on 1 and 4 processes. testif chains are honoured:
#                      t11 runs only if t10 passed, t12/t13 only if t11 passed.
#   TaylorImpact.py    the plain-SPH level-100 lines, --steps 100 --compatibleEnergy False: 2d and RZ on 1 and 8
#                      processes, 3d on 8 only (the header has no 1-process 3d line). No stored reference exists for it
#                      upstream; these lines only generate silo snapshots, so "passes" means "runs to completion".
#
# The upstream script compares its output file with the reference only on the --checkRef True runs, and asserts the
# bitwise comparison with --comparisonFile only inside that same branch, so t11-t13 assert nothing about their output.
# This script therefore reads every TensileRod output file itself (numpy, in the container) and reports, per file,
# the largest relative difference from LLNL's reference, |a-b| / max(mean(|a|,|b|), 1e-4) (the metric commented in
# filearraycmp.py), the largest ratio of |a-b| to filearraycmp's allowance 1e-4 + 1e-4*mean(|a|,|b|) (pass <= 1), and
# whether the file is byte-identical to the serial run's.
#
# Writes OUTDIR/summary.tsv (one row per run: case, procs, exit status, outcome, wall seconds), OUTDIR/compare.tsv
# (one row per TensileRod output file) and OUTDIR/logs/<case>.log. Exit 0 when every run and comparison passed, 1
# otherwise, 2 when the launcher cannot run (no Docker, no image).
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo="$(cd "$here/../.." && pwd -P)"
launcher="${SPHERAL:-$repo/spheral_frag/container/spheral}"
out="${1:-$repo/spheral_output/regressions}"
mkdir -p "$out/logs"
out="$(cd "$out" && pwd -P)"
case "$out/" in "$repo"/*) ;; *) echo "run_regressions: OUTDIR must be inside $repo" >&2; exit 2 ;; esac

rod_src=/opt/spheral/tests/functional/Damage/TensileRod
taylor_src=/opt/spheral/tests/functional/Strength/TaylorImpact
now() { perl -MTime::HiRes=time -e 'printf "%.3f\n", time'; }

summary="$out/summary.tsv"
printf "case\tprocs\texit\toutcome\twall_s\n" > "$summary"
failures=0

# run CASE NPROC DIR SCRIPT ARGS...: one launcher run in DIR, logged and timed; sets $status.
run() {
    local name="$1" nproc="$2" dir="$3"; shift 3
    local t0 t1 outcome
    echo "== $name (np=$nproc)"
    t0="$(now)"
    (cd "$dir" && "$launcher" -n "$nproc" "$@") > "$out/logs/$name.log" 2>&1
    status=$?
    t1="$(now)"
    if [[ $status -eq 0 ]]; then outcome=pass; else outcome=FAIL; failures=$((failures + 1)); fi
    printf "%s\t%s\t%s\t%s\t%.1f\n" "$name" "$nproc" "$status" "$outcome" "$(echo "$t1 - $t0" | bc)" >> "$summary"
    printf "   %s, exit %s, %.1f s\n" "$outcome" "$status" "$(echo "$t1 - $t0" | bc)"
}
skip() {
    echo "== $1 (np=$2) skipped: $3"
    printf "%s\t%s\t-\tskipped (%s)\t-\n" "$1" "$2" "$3" >> "$summary"
    failures=$((failures + 1))
}

# The launcher must work at all; otherwise exit 2 like it does.
if ! (cd "$repo" && "$launcher" -c "import Spheral" > "$out/logs/probe.log" 2>&1); then
    echo "run_regressions: the launcher cannot import Spheral, see $out/logs/probe.log" >&2
    exit 2
fi

# ---------------------------------------------------------------------------------------------------- tensile rod
# refDir is "Reference/..." relative to the working directory, so LLNL's Reference tree is copied in next to the
# dumps. Both damage models write to their own dumps-TensileRod-1d/SPH/<model>/... directory.
rod="$out/tensilerod"
mkdir -p "$rod"
(cd "$rod" && "$launcher" -c "import shutil; shutil.copytree('$rod_src/Reference', 'Reference', dirs_exist_ok=True)") \
    > "$out/logs/tensilerod_copy_reference.log" 2>&1 || { echo "run_regressions: cannot copy Reference" >&2; exit 1; }

for dm in GradyKippTensorDamageOwen:t1 ProbabilisticDamageModel:t2; do
    model="${dm%%:*}"; t="${dm##*:}"
    common=(--DamageModelConstructor "$model" --graphics False --domainIndependent True)
    run "tensilerod_${t}0_${model}_serial_checkRef" 1 "$rod" "$rod_src/TensileRod-1d.py" "${common[@]}" \
        --clearDirectories True --outputFile TensileRod-1d-1proc.gnu --checkRef True
    if [[ $status -ne 0 ]]; then
        skip "tensilerod_${t}1_${model}_4proc" 4 "testif ${t}0 failed"
        skip "tensilerod_${t}2_${model}_serial_restart500" 1 "testif ${t}1 failed"
        skip "tensilerod_${t}3_${model}_4proc_restart500" 4 "testif ${t}1 failed"
        continue
    fi
    run "tensilerod_${t}1_${model}_4proc" 4 "$rod" "$rod_src/TensileRod-1d.py" "${common[@]}" \
        --clearDirectories False --outputFile TensileRod-1d-4proc.gnu --comparisonFile TensileRod-1d-1proc.gnu
    if [[ $status -ne 0 ]]; then
        skip "tensilerod_${t}2_${model}_serial_restart500" 1 "testif ${t}1 failed"
        skip "tensilerod_${t}3_${model}_4proc_restart500" 4 "testif ${t}1 failed"
        continue
    fi
    run "tensilerod_${t}2_${model}_serial_restart500" 1 "$rod" "$rod_src/TensileRod-1d.py" "${common[@]}" \
        --clearDirectories False --outputFile TensileRod-1d-1proc-restart.gnu --comparisonFile TensileRod-1d-1proc.gnu \
        --restoreCycle 500
    run "tensilerod_${t}3_${model}_4proc_restart500" 4 "$rod" "$rod_src/TensileRod-1d.py" "${common[@]}" \
        --clearDirectories False --outputFile TensileRod-1d-4proc-restart.gnu --comparisonFile TensileRod-1d-1proc.gnu \
        --restoreCycle 500
done

# Every output file against the reference and against the serial run.
(cd "$rod" && "$launcher" -c "
import filecmp, glob, os
import numpy as np
rows = ['model\tfile\tmax_rel_diff_vs_ref\tworst_column\tmax_ratio_to_allowance\tpass_1e-4\tidentical_to_1proc']
cols = ['x', 'rho', 'P', 'v', 'eps', 'h', 'S', 'D']
bad = 0
for model in ('GradyKippTensorDamageOwen', 'ProbabilisticDamageModel'):
    sub = os.path.join('SPH', model, 'DamageCouplingAlgorithm.PairMaxDamage', 'nx=100', 'k=652000.00_m=2.63')
    ref = np.loadtxt(os.path.join('Reference', sub, 'TensileRod-1d-1proc.gnu'))
    for name in ('TensileRod-1d-1proc.gnu', 'TensileRod-1d-4proc.gnu',
                 'TensileRod-1d-1proc-restart.gnu', 'TensileRod-1d-4proc-restart.gnu'):
        path = os.path.join('dumps-TensileRod-1d', sub, name)
        if not os.path.exists(path):
            rows.append(f'{model}\t{name}\tmissing\t-\t-\tFAIL\t-'); bad += 1; continue
        a = np.loadtxt(path)
        if a.shape != ref.shape:
            rows.append(f'{model}\t{name}\tshape {a.shape} vs {ref.shape}\t-\t-\tFAIL\t-'); bad += 1; continue
        d = np.abs(a - ref); m = 0.5 * (np.abs(a) + np.abs(ref))
        rel = d / np.clip(m, 1e-4, None)
        ratio = d / (1e-4 + 1e-4 * m)
        ok = bool(ratio.max() <= 1.0)
        bad += not ok
        same = filecmp.cmp(path, os.path.join('dumps-TensileRod-1d', sub, 'TensileRod-1d-1proc.gnu'), shallow=False)
        bad += not same
        rows.append(f'{model}\t{name}\t{rel.max():.3e}\t{cols[int(np.argmax(rel.max(axis=0)))]}\t{ratio.max():.3e}\t'
                    f'{\"pass\" if ok else \"FAIL\"}\t{same}')
open('../compare.tsv', 'w').write('\n'.join(rows) + '\n')
open('../compare.status', 'w').write(str(bad))
") > "$out/logs/tensilerod_compare.log" 2>&1
if [[ "$(cat "$out/compare.status" 2>/dev/null)" != "0" ]]; then failures=$((failures + 1)); fi
rm -f "$out/compare.status"

# -------------------------------------------------------------------------------------------------- Taylor impact
# One working directory per run, so every snapshot (written under dumps-TaylorImpact/<geometry>/.../procs=N) is kept.
for gp in 2d:1 2d:8 RZ:1 RZ:8 3d:8; do
    geom="${gp%%:*}"; np="${gp##*:}"
    lower="$(echo "$geom" | tr 'A-Z' 'a-z')"
    dir="$out/taylorimpact_sph_${lower}_${np}proc"
    mkdir -p "$dir"
    run "taylorimpact_sph_${lower}_${np}proc_steps100" "$np" "$dir" "$taylor_src/TaylorImpact.py" \
        --geometry "$geom" --hydroType SPH --steps 100 --compatibleEnergy False --clearDirectories True \
        --siloSnapShotFile "Spheral_sph_${lower}_state_snapshot_${np}proc"
done

echo
column -t -s $'\t' "$summary"
echo
column -t -s $'\t' "$out/compare.tsv" 2>/dev/null || echo "(no compare.tsv, see $out/logs/tensilerod_compare.log)"
echo
echo "summary: $summary"
[[ $failures -eq 0 ]] && exit 0 || exit 1
