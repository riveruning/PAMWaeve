#!/usr/bin/env bash
#
# One-command reproduction of the submitted candidate list.
#
# Steps:
#   1. validate the frozen manifest and every registered artifact hash
#      (no model needed);
#   2. run the system-level benchmark with real model inference;
#   3. project the frozen artifacts into results/results.csv.
#
# Prerequisites (see README.md section 1):
#   - pip install -r requirements.txt
#   - Protein2PAM cas9_full weights under $HF_HOME
#   - Protein2PAM upstream code at .reference/Protein2PAM
#
# This script does not pin a developer-machine interpreter: it uses `python3`
# from the active environment.  It never writes outside this package.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

PY="${PYTHON:-python3}"
MANIFEST="benchmarks/multi_evidence_v2/manifest_redistributable21.json"
OUTDIR="data/parsed/bench_repro21"

export PYTHONPATH="src:${PYTHONPATH:-}"
export HF_HOME="${HF_HOME:-$ROOT/data/checkpoints/hf}"

# Default to the 21-system redistributable manifest: every artifact it lists is
# shipped with this package, so validation passes out of the box.  The full
# 22-system manifest additionally needs the PAMpredict SpCas9 inputs, which are
# CC-BY-NC-ND-4.0 and are deliberately not redistributed.  Switch to the full
# manifest only when those inputs are actually present.
FULL_MANIFEST="benchmarks/multi_evidence_v2/manifest_final22.json"
if [[ -f data/raw/external_datasets/PAMpredict/Example/spacers.fna \
      && -f data/raw/cas9_full.fasta ]]; then
  echo "NOTE: PAMpredict SpCas9 inputs found; running the full 22-system manifest." >&2
  RUN_MANIFEST="$FULL_MANIFEST"
  OUTDIR="data/parsed/bench_final22"
  SYSTEMS="spcas9-pampredict-example,sp7f7-published-flanks,cj4-campylobacter-jejuni-414-protein-only,asp-acidiphilium-21-60-14-protein-only,cba-caulobacterales-protein-only"
  RESULTS_OUT="results/results.csv"
else
  echo "NOTE: PAMpredict SpCas9 paired-evidence inputs not found." >&2
  echo "      Running the 21-system redistributable subset (this is the" >&2
  echo "      documented default). The submitted results/results.csv has 25" >&2
  echo "      rows from the full 22-system run; this run yields 20 rows." >&2
  echo "      See docs/THIRD_PARTY.md section 2 to obtain the missing inputs." >&2
  RUN_MANIFEST="$MANIFEST"
  SYSTEMS="sp7f7-published-flanks,cj4-campylobacter-jejuni-414-protein-only,asp-acidiphilium-21-60-14-protein-only,cba-caulobacterales-protein-only"
  RESULTS_OUT="results/results_repro21.csv"
fi

echo "==> [1/3] validating manifest artifacts (no model)"
"$PY" scripts/validate_benchmark_manifest.py --manifest "$RUN_MANIFEST" --check-files

echo "==> [2/3] running system-level benchmark (real inference, CPU)"
"$PY" -u scripts/run_multi_evidence_benchmark.py \
  --manifest "$RUN_MANIFEST" \
  --device "${DEVICE:-cpu}" \
  --fusion-policy abstain \
  --bootstrap 2000 \
  --reuse-existing \
  --outdir "$OUTDIR"

echo "==> [3/3] building the candidate list"
"$PY" scripts/build_competition_results.py \
  --benchmark-dir "$OUTDIR" \
  --manifest "$RUN_MANIFEST" \
  --display-candidates results/display_candidates.json \
  --top-n 5 \
  --systems "$SYSTEMS" \
  --out "$RESULTS_OUT"

echo
echo "DONE. Candidate list: $ROOT/$RESULTS_OUT"
echo "Report:               $ROOT/$OUTDIR/benchmark_report.md"
