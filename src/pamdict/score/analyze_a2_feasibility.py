"""A2 — sequence-level feasibility analysis for multi-evidence (protein + spacer/repeat) fusion.

Question: for CRISPR systems where a protein-only model would mis-predict the
PAM, can DNA-sequence evidence (spacers, direct repeats) supply information
that the protein does not already carry?

This is a DATA-SIDE analysis only (no model). It answers the upstream question:
(1) is the DNA evidence even present/available, and (2) among protein-similar
systems do PAMs diverge, i.e., is there 'residual' that DNA might resolve.

Metrics computed on the official training TSV:
  - general statistics of the training set.
  - pairwise protein k-mer similarity vs PAM agreement, binned, to show whether
    similar proteins sometimes carry different PAMs (the 'mis-predictable'
    population).
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict

from pamdict.score.similarity import kmer_set, jaccard, build_kmer_index


def load_rows(tsv: str) -> list[dict]:
    with open(tsv, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--sample", type=int, default=3000,
                    help="max proteins to index for pairwise analysis (CPU bound)")
    args = ap.parse_args()

    rows = load_rows(args.input)
    cols = list(rows[0].keys())

    stats = {
        "total_rows": len(rows),
        "columns": cols,
        "has_spacer_or_repeat_column": any(
            "spacer" in c.lower() or "repeat" in c.lower() or "tracr" in c.lower()
            for c in cols),
        "cas_family_dist": dict(Counter(r["cas_family"] for r in rows)),
        "crispr_type_dist": dict(Counter(r["crispr_type"] for r in rows)),
    }

    # PAM diversity among protein-similar systems (residual that DNA might explain)
    # Use a bounded sample and k-mer-indexed candidate search (not O(n^2)).
    sample = rows[: args.sample]
    seqs = [r["protein_sequence"] for r in sample]
    index = build_kmer_index(seqs, 4)

    # For each protein, find nearest neighbor among candidates sharing k-mers.
    pair_records = []
    n = len(sample)
    for i in range(n):
        qi = kmer_set(seqs[i], 4)
        cand = set()
        for km in qi:
            cand |= index.get(km, set())
        cand.discard(i)
        if not cand:
            continue
        best_sim, best_j = 0.0, -1
        for j in cand:
            s = jaccard(qi, kmer_set(seqs[j], 4))
            if s > best_sim:
                best_sim, best_j = s, j
        if best_j >= 0:
            pair_records.append({
                "sim": best_sim,
                "same_pam": sample[i]["pam_consensus"] == sample[best_j]["pam_consensus"],
                "pam_i": sample[i]["pam_consensus"],
                "pam_j": sample[best_j]["pam_consensus"],
                "family": sample[i]["cas_family"],
            })

    # Bin similarity -> PAM agreement rate. A low agreement rate at HIGH
    # protein similarity is exactly the 'mis-predictable' population where DNA
    # evidence could matter.
    bins = [(0.0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 1.01)]
    bin_stats = []
    for lo, hi in bins:
        grp = [p for p in pair_records if lo <= p["sim"] < hi]
        if not grp:
            continue
        agree = sum(1 for p in grp if p["same_pam"]) / len(grp)
        bin_stats.append({
            "sim_bin": f"{lo:.2f}-{hi:.2f}",
            "n_pairs": len(grp),
            "pam_agreement_rate": round(agree, 4),
            "mean_sim": round(sum(p['sim'] for p in grp) / len(grp), 4),
        })

    # Count how many pairs are 'similar protein but different PAM' (the target
    # population that a DNA signal might rescue).
    target = [p for p in pair_records if p["sim"] >= 0.5 and not p["same_pam"]]
    report = {
        "stats": stats,
        "sampled_for_pairwise": n,
        "pairs_analyzed": len(pair_records),
        "similarity_bins": bin_stats,
        "similar_but_diff_pam_pairs": len(target),
        "conclusion_notes": [
            "If official training data has no spacer/repeat columns, DNA evidence "
            "is absent from the corpus by construction -> fusion requires external "
            "CRISPR array data (spacers, direct repeats) not present in "
            "Protein2PAM's training TSV.",
            "A low PAM-agreement rate at high protein similarity identifies the "
            "population a protein-only model is structurally likely to mis-predict.",
        ],
    }

    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
