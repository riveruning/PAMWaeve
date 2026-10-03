"""Finalize the gold set: consensus-PAM filtering + a human-review queue.

Input:  gold_set.tsv          (sequence + PAM + confidence, from merge_gold.py)
        gold_unresolved.tsv   (Cas->PAM pairs that could not be sequenced yet)

Output:
  gold_clean.tsv    — rows that are BOTH sequence-verified (verified/species)
                      AND whose PAM matches the Cas's consensus PAM. These are
                      the trustworthy (sequence, PAM) samples for protein2PAM.
  gold_review.tsv   — everything else that needs a human look: sequence from an
                      untrusted species anchor, unknown-alias PAM, or a known
                      Cas with a non-consensus PAM (recorded as filter_reason).
  cas_to_verify.tsv — the cold-start Cas aliases (unresolved / review) as a
                      manual queue: alias + paper species + suggested action.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict

from pamdict.collect import consensus_pam


def load_rows(tsv: str) -> list[dict[str, str]]:
    with open(tsv, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_rows(path: str, rows: list[dict[str, str]], fields: list[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gold", default="data/parsed/cas_pam_text/gold_set.tsv")
    ap.add_argument("--unresolved", default="data/parsed/cas_pam_text/gold_unresolved.tsv")
    ap.add_argument("--out-clean", default="data/parsed/cas_pam_text/gold_clean.tsv")
    ap.add_argument("--out-review", default="data/parsed/cas_pam_text/gold_review.tsv")
    ap.add_argument("--out-verify", default="data/parsed/cas_pam_text/cas_to_verify.tsv")
    args = ap.parse_args()

    gold = load_rows(args.gold)

    TRUSTED = {"verified", "species"}

    clean: list[dict[str, str]] = []
    review: list[dict[str, str]] = []
    verify: list[dict[str, str]] = []

    # Consensus-PAM filter over ALL resolved rows (regardless of confidence).
    kept, dropped = consensus_pam.filter_pairs(gold, pam_field="pam_consensus",
                                               cas_field="cas_name")

    for r in kept:
        if r.get("sequence_confidence") in TRUSTED:
            clean.append(r)
        else:
            # Sequence not trusted (species_review/name_only) -> review.
            review.append({**r, "review_reason": f"sequence_confidence={r.get('sequence_confidence')}"})

    for d in dropped:
        review.append({**d, "review_reason": d.get("filter_reason", "")})

    # Human-review queue of cold-start Cas aliases: aggregate unresolved + review.
    alias_info: dict[str, dict[str, str]] = {}
    for r in gold:
        alias = r.get("cas_name", "").strip()
        if not alias:
            continue
        info = alias_info.setdefault(alias, {
            "cas_name": alias,
            "paper_species": "", "pams_seen": "", "has_seq": "no",
            "seq_confidence": "", "suggested_action": "",
        })
        if r.get("species") and not info["paper_species"]:
            info["paper_species"] = r.get("species", "")
        seen = info["pams_seen"].split(",") if info["pams_seen"] else []
        p = r.get("pam_consensus", "").strip()
        if p and p not in seen:
            seen.append(p)
        info["pams_seen"] = ",".join(sorted(seen)[:12])
        if r.get("sequence"):
            info["has_seq"] = "yes"
            info["seq_confidence"] = r.get("sequence_confidence", "")

    # Alias-level queue: any alias NOT in the curated alias map is a candidate
    # for manual source-organism confirmation (its sequence is unresolved or
    # only from an untrusted paper-species anchor).
    from pamdict.collect import cas_alias_map
    queue: list[dict[str, str]] = []
    for alias, info in alias_info.items():
        curated = cas_alias_map.organism_of(alias)
        if curated is None:
            queue.append({
                "cas_name": info["cas_name"],
                "paper_species": info["paper_species"],
                "pams_seen": info["pams_seen"],
                "has_seq": info["has_seq"],
                "seq_confidence": info["seq_confidence"],
                "suggested_action": "confirm source organism, then re-resolve sequence",
            })

    # Write outputs.
    gold_fields = list(gold[0].keys()) if gold else []
    write_rows(args.out_clean, clean, gold_fields)
    write_rows(args.out_review, review, gold_fields + ["review_reason"])
    verify_fields = ["cas_name", "paper_species", "pams_seen", "has_seq",
                     "seq_confidence", "suggested_action"]
    write_rows(args.out_verify, queue, verify_fields)

    print(f"gold_clean:   {len(clean)} rows (verified seq + consensus PAM)")
    print(f"gold_review:  {len(review)} rows (needs human)")
    print(f"cas_to_verify: {len(queue)} cold-start aliases")
    return 0


if __name__ == "__main__":
    exit(main())
