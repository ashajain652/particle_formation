#!/bin/bash
# mkdiff.sh <old tree> <new tree> <out.diff>: unified diff (3 lines of context) of every .py file that differs, with
# a/ b/ labels the rebuild script's strict applier reads
old=$1; new=$2; out=$3; : > "$out"
cd "$new"
for f in $(find reentry_model tests analysis -name '*.py' | sort); do
  if ! cmp -s "$old/$f" "$f"; then diff -u --label "a/$f" --label "b/$f" "$old/$f" "$f" >> "$out"; fi
done
