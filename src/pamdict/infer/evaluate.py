"""Evaluate Protein2PAM predictions against the clean gold set.

Loads ``gold_clean.tsv`` (Cas name + consensus PAM + provenance), resolves the
correct model input (Cas9 -> PI-domain from the official training TSV, Cas12 ->
full sequence), runs the matching cas9/cas12 model, and computes accuracy of the
predicted consensus PAM against the gold PAM.

Accuracy is reported two ways:
  * exact   — predicted consensus string equals the gold string
  * iupac   — predicted motif is consistent with the gold motif under IUPAC
              degenerate-code expansion (a specific spelling of the same motif)

A ``leakage`` column flags whether the sample's full sequence is present in the
official training TSV (i.e. the model has seen it — the accuracy is then an
upper bound, not a generalization estimate).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from pamdict.collect import consensus_pam as cp
from pamdict.infer import p2pam

WORKSPACE = Path(__file__).resolve().parents[3]
GOLD = WORKSPACE / "data/parsed/cas_pam_text/gold_clean.tsv"
OFFICIAL = WORKSPACE / "data/raw/protein2pam_train_seqs.tsv"


def build_pid_index(tsv: Path) -> dict[str, str]:
    """Map full protein_sequence -> pid_sequence from the official training TSV."""
    idx: dict[str, str] = {}
    with tsv.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("pid_sequence"):
                idx[row["protein_sequence"]] = row["pid_sequence"]
    return idx


def build_train_seq_set(tsv: Path) -> set[str]:
    """Set of all full protein sequences in the official training TSV."""
    s: set[str] = set()
    with tsv.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s.add(row["protein_sequence"])
    return s


def iupac_match(pred: str, gold: str) -> bool:
    """True if `pred` is a degenerate/specific spelling of `gold` (same length)."""
    pred, gold = pred.upper(), gold.upper()
    if len(pred) != len(gold):
        return False
    for a, b in zip(pred, gold):
        allowed = cp.IUPAC.get(a)
        gold_allowed = cp.IUPAC.get(b)
        if allowed is None or gold_allowed is None:
            return False
        # pred's allowed set must be a subset of gold's allowed set
        if not allowed <= gold_allowed:
            return False
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gold", default=str(GOLD))
    ap.add_argument("--official", default=str(OFFICIAL))
    ap.add_argument("--out", default=str(WORKSPACE / "data/parsed/cas_pam_text/eval_result.json"))
    ap.add_argument("--detail", default=str(WORKSPACE / "data/parsed/cas_pam_text/eval_detail.tsv"))
    ap.add_argument("--device", default=None)
    args = ap.parse_args()

    pid_idx = build_pid_index(Path(args.official))
    train_seqs = build_train_seq_set(Path(args.official))

    gold_rows = list(csv.DictReader(open(args.gold, encoding="utf-8"), delimiter="\t"))

    # Group by model needed.
    predictors: dict[str, p2pam.P2PAMPredictor] = {}
    def get_pred(model: str):
        if model not in predictors:
            predictors[model] = p2pam.P2PAMPredictor(model, args.device)
        return predictors[model]

    # Prepare inputs + run per-model batches to avoid reloading.
    cas9_rows = [r for r in gold_rows if r.get("cas_family") == "Cas9"]
    cas12_rows = [r for r in gold_rows if r.get("cas_family") == "Cas12"]

    results: list[dict] = []

    def run_rows(rows, model):
        nonlocal results
        if not rows:
            return
        inputs = []
        for r in rows:
            full = r["protein_sequence"]
            if model == "cas9":
                pid = pid_idx.get(full)
                inputs.append(pid if pid else full)  # fallback: full seq if pid unknown
            else:
                inputs.append(full)
        pred = get_pred(model)
        pams = pred.predict(inputs)
        for r, inp, p in zip(rows, inputs, pams):
            gold = r["pam_consensus"]
            leakage = "yes" if r["protein_sequence"] in train_seqs else "no"
            results.append({
                "cas_name": r["cas_name"], "model": model, "gold_pam": gold,
                "pred_pam": p, "exact": p == gold, "iupac": iupac_match(p, gold),
                "leakage": leakage, "doi": r.get("doi", ""),
                "input_len": len(inp),
            })

    run_rows(cas9_rows, "cas9")
    run_rows(cas12_rows, "cas12")

    # Summary.
    n = len(results)
    exact = sum(1 for r in results if r["exact"])
    iupac = sum(1 for r in results if r["iupac"])
    n_leak = sum(1 for r in results if r["leakage"] == "yes")
    summary = {
        "n_samples": n,
        "exact_match": exact,
        "iupac_match": iupac,
        "exact_accuracy": round(exact / n, 4) if n else 0.0,
        "iupac_accuracy": round(iupac / n, 4) if n else 0.0,
        "n_in_training_set": n_leak,
    }
    Path(args.out).write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    fields = ["cas_name", "model", "gold_pam", "pred_pam", "exact", "iupac",
              "leakage", "input_len", "doi"]
    with open(args.detail, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(results)

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("\n逐条明细:")
    for r in results:
        mark = "✓" if r["iupac"] else "✗"
        print(f"  {mark} {r['cas_name']:12s} [{r['model']}] gold={r['gold_pam']:10s} "
              f"pred={r['pred_pam']:10s} leak={r['leakage']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
