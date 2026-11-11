#!/bin/bash
# Pause the Task 9 model runs (SIGSTOP; they keep their state) when on battery below 40 %, resume them (SIGCONT) on mains.
LOG="$(dirname "$0")/power_guard.log"
paused=0
while true; do
  b=$(pmset -g batt); pct=$(echo "$b" | grep -o '[0-9]*%' | head -1 | tr -d '%')
  if echo "$b" | grep -q "Battery Power" && [ "${pct:-100}" -lt 40 ] && [ $paused = 0 ]; then
    pkill -STOP -f "reentry_model run"; pkill -STOP -f "box_ablation.py"; paused=1
    echo "$(date '+%F %T') paused at ${pct}% on battery" >> "$LOG"
  elif echo "$b" | grep -q "AC Power" && [ $paused = 1 ]; then
    pkill -CONT -f "reentry_model run"; pkill -CONT -f "box_ablation.py"; paused=0
    echo "$(date '+%F %T') resumed on mains at ${pct}%" >> "$LOG"
  fi
  sleep 60
done
