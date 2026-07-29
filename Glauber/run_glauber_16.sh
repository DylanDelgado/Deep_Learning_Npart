#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-/workspace/miniforge3/envs/urqmd/bin/python}"
NUM_JOBS=16
EVENTS_PER_JOB=20000
OUTPUT_DIR="$SCRIPT_DIR/glauber_output"

mkdir -p "$OUTPUT_DIR"

pids=()
for worker_id in $(seq 0 $((NUM_JOBS - 1))); do
    "$PYTHON_BIN" "$SCRIPT_DIR/init_glauber.py" "$worker_id" "$EVENTS_PER_JOB" &
    pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        status=1
    fi
done

if (( status != 0 )); then
    echo "At least one Glauber worker failed." >&2
    exit "$status"
fi

output_count=$(find "$OUTPUT_DIR" -maxdepth 1 -type f -name 'worker_*.npz' | wc -l)
if [[ "$output_count" -ne "$NUM_JOBS" ]]; then
    echo "Expected $NUM_JOBS worker files, found $output_count." >&2
    exit 1
fi

printf 'Created %s Glauber files with %s requested events each in %s\n' \
    "$output_count" "$EVENTS_PER_JOB" "$OUTPUT_DIR"