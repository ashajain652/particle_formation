#!/bin/zsh
# launch_run.sh <run name> <cli args...>: one plain run (run.py: the copy's own cli.main, full-precision history and
# final state in full.npz) of the frozen measurement copy fluxfix/meas into fluxfix/runs/<run name>, under caffeinate,
# with the power source logged at the start and the end.
FF="/private/tmp/claude-501/-Users-ashajain-Documents-University-Documents--MIT-Graduate-Work-Research-Space-Sustainability-Particle-Wake-Evolution/930f32b6-e08a-42f0-a4ce-595da1b5cc42/scratchpad/fluxfix"
NAME="$1"; shift
OUT="$FF/runs/$NAME"
if [ -e "$OUT" ]; then echo "refusing: $OUT exists"; exit 3; fi
cd "$FF/meas"
echo "start $(date '+%F %T') $(pmset -g batt | head -1)" > "$FF/runs/$NAME.log"
caffeinate -i /Users/ashajain/miniforge3/envs/drama_env/bin/python "$FF/harness/run.py" "$FF/meas" "$OUT" -- "$@" >> "$FF/runs/$NAME.log" 2>&1
echo "exit $?" >> "$FF/runs/$NAME.log"
echo "end $(date '+%F %T') $(pmset -g batt | head -1)" >> "$FF/runs/$NAME.log"
