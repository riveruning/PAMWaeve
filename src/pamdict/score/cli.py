"""CLI: score a raw PAM prediction matrix and print a reliability verdict.

Usage (matrix from a JSON file, 10x4 A/C/G/T):
  PYTHONPATH=src python -m pamdict.score.cli \
      --matrix data/parsed/demo_pamreadid/...  # not wired here; see --help

Alternative (produce a matrix inline via Python). This CLI focuses on the
verdict logic; a FASTA->model->matrix front-end reuses PAM_PREDICT predict.py's
interface and is intentionally left for target #3 (finetune) integration.
"""

from __future__ import annotations

import argparse
import json
import sys

from .scoring import score_matrix


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix-json", required=True,
                    help="path to a JSON file: [[A,C,G,T], ...] (10 positions)")
    ap.add_argument("--units", choices=("bits", "prob"), default="prob")
    ap.add_argument("--neighbor-sim", type=float, default=1.0,
                    help="nearest-neighbor sequence identity to training, 0..1")
    ap.add_argument("--model", default="")
    args = ap.parse_args()

    with open(args.matrix_json, encoding="utf-8") as fh:
        matrix = json.load(fh)

    score = score_matrix(matrix, units=args.units,
                         neighbor_sim=args.neighbor_sim, model_name=args.model)
    print(json.dumps(score.to_dict(), indent=2))


if __name__ == "__main__":
    main()
