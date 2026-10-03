#!/usr/bin/env bash
# Validate inputs, run inference/evaluation, export candidates and save logs.
# Setup: docs/REPRODUCTION.md. Use --fresh to rerun model inference.
set -euo pipefail
TASK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$TASK_ROOT"
exec "${PYTHON:-python3}" scripts/run_submission.py "$@"
