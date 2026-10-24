#!/usr/bin/env bash
# Build the native arm64 Spheral image at 116c71f. Log: spheral_output/build/build-<timestamp>.log.
# Extra arguments go to docker build (e.g. --build-arg JPY=2, --no-cache).
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tag="particle-formation/spheral:116c71f-arm64"
logdir="$repo/spheral_output/build"
mkdir -p "$logdir"
log="$logdir/build-$(date +%Y%m%d-%H%M%S).log"

start=$(date +%s)
echo "building $tag, log $log"
status=0
docker build --platform linux/arm64 --progress=plain -t "$tag" "$@" \
    "$repo/spheral_frag/container" 2>&1 | tee "$log" || status=$?
end=$(date +%s)
wall=$((end - start))
printf 'build wall time: %dh %02dm %02ds (%d s), exit %d\n' \
    $((wall / 3600)) $((wall % 3600 / 60)) $((wall % 60)) "$wall" "$status" | tee -a "$log"
exit "$status"
