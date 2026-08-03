#!/usr/bin/env bash
set -e

cd "$(dirname "$0")"

for worker_id in {0..15}; do
    python init_glauber.py "$worker_id" &
done

wait
