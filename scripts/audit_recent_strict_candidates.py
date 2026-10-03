#!/usr/bin/env python3
"""Reproduce the 2026-09-25 recent strict protein-only candidate audit."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.finetune.similarity import ProteinNeighborIndex  # noqa: E402


CANDIDATES = (
    {
        "candidate_id": "CoCas9",
        "accession": "C7M7G9",
        "taxon": "Capnocytophaga ochracea",
        "doi": "10.1080/15476286.2023.2256578",
        "gold_pam": "NRRWC",
        "source": "data/parsed/fulltext_golds_curated.tsv",
        "source_kind": "tsv",
        "exact_experimental_sequence_link": True,
        "evidence_note": "Biochemically characterized exact protein; excluded by >=90% official-training neighbor.",
    },
    {
        "candidate_id": "AalCas9",
        "accession": "A0A3N5C319",
        "taxon": "Abyssicoccus albus",
        "doi": "10.1089/crispr.2024.0013",
        "gold_pam": "NNACR",
        "source": "data/corpus/computed_training.tsv",
        "source_kind": "tsv",
        "exact_experimental_sequence_link": False,
        "evidence_note": "Paper and UniProt both report a 1059-aa A. albus Cas9, but the publisher supplement linking the experimental construct to the accession was unavailable during this audit.",
    },
    {
        "candidate_id": "Cj4Cas9",
        "accession": "EFC33367.1",
        "taxon": "Campylobacter jejuni subsp. jejuni 414",
        "doi": "10.1038/s42003-025-09430-9",
        "gold_pam": "NNNGRY",
        "source": "benchmarks/multi_evidence_v2/systems/cj4-campylobacter-jejuni-414-protein-only/protein.faa",
        "source_kind": "fasta",
        "exact_experimental_sequence_link": True,
        "evidence_note": "Article Table 1 links Cj4Cas9 to EFC33367.1; randomized 7-bp PAM-library deep sequencing reports NNNGRY.",
    },
)


def is_cas9(row: dict[str, str]) -> bool:
    family = (row.get("cas_family") or "").strip().upper()
    crispr_type = (row.get("crispr_type") or "").replace(" ", "").upper()
    return family.startswith("CAS9") or crispr_type in {"II", "TYPEII"}


def read_sequence(candidate: dict[str, object]) -> str:
    path = WORKSPACE / str(candidate["source"])
    if candidate["source_kind"] == "fasta":
        return "".join(
            line.strip() for line in path.read_text().splitlines()
            if not line.startswith(">")
        ).upper()
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("protein_id") == candidate["accession"]:
                return row["protein_sequence"].strip().upper()
    raise ValueError(f"sequence not found for {candidate['accession']} in {path}")


def project_exposures(sequence: str, accession: str) -> list[str]:
    hits: list[str] = []
    for path in sorted((WORKSPACE / "data/corpus").glob("*.tsv")):
        try:
            with path.open(encoding="utf-8") as handle:
                rows = csv.DictReader(handle, delimiter="\t")
                if any(
                    row.get("protein_id") == accession
                    or row.get("protein_sequence", "").strip().upper() == sequence
                    for row in rows
                ):
                    hits.append(str(path.relative_to(WORKSPACE)))
        except (UnicodeDecodeError, csv.Error):
            continue
    return hits


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--official", default="data/raw/protein2pam_train_seqs.tsv")
    parser.add_argument("--identity", type=float, default=0.90)
    parser.add_argument(
        "--out",
        default="benchmarks/multi_evidence_v2/recent_strict_candidate_audit.tsv",
    )
    parser.add_argument(
        "--summary-out",
        default="benchmarks/multi_evidence_v2/recent_strict_candidate_audit_summary.json",
    )
    args = parser.parse_args()

    with (WORKSPACE / args.official).open(encoding="utf-8") as handle:
        official_rows = [
            row for row in csv.DictReader(handle, delimiter="\t") if is_cas9(row)
        ]
    references = [row["protein_sequence"].strip().upper() for row in official_rows]
    index = ProteinNeighborIndex(references)
    sequences = {str(row["candidate_id"]): read_sequence(row) for row in CANDIDATES}
    matches = index.find(list(sequences.values()), identity_threshold=args.identity)
    sensitivity_threshold = min(0.85, args.identity)
    sensitivity_matches = index.find(
        list(sequences.values()), identity_threshold=sensitivity_threshold
    )

    output: list[dict[str, object]] = []
    for candidate in CANDIDATES:
        candidate_id = str(candidate["candidate_id"])
        sequence = sequences[candidate_id]
        neighbor_rows = matches.get(sequence, [])
        sensitivity_rows = sensitivity_matches.get(sequence, [])
        sensitivity_top = sensitivity_rows[0] if sensitivity_rows else None
        exact = sequence in index.exact
        exposures = project_exposures(sequence, str(candidate["accession"]))
        exact_link = bool(candidate["exact_experimental_sequence_link"])
        if exact:
            decision = "EXCLUDE_OFFICIAL_EXACT"
        elif neighbor_rows:
            decision = "EXCLUDE_OFFICIAL_NEAR90"
        elif not exact_link:
            decision = "BLOCKED_EXACT_SEQUENCE_LINK"
        elif exposures:
            decision = "RESERVE_PROJECT_TRAINING_EXPOSURE"
        else:
            decision = "REGISTER_STRICT_PROTEIN_ONLY"
        output.append({
            "candidate_id": candidate_id,
            "accession": candidate["accession"],
            "taxon": candidate["taxon"],
            "protein_length": len(sequence),
            "doi": candidate["doi"],
            "gold_pam": candidate["gold_pam"],
            "exact_experimental_sequence_link": exact_link,
            "official_exact_match": exact,
            "official_near90_match": bool(neighbor_rows),
            "top_official_identity": neighbor_rows[0]["identity"] if neighbor_rows else "",
            "top_official_identity_at_or_above_85": (
                sensitivity_top["identity"] if sensitivity_top else ""
            ),
            "margin_below_90_threshold": (
                round(args.identity - float(sensitivity_top["identity"]), 8)
                if sensitivity_top and not neighbor_rows else ""
            ),
            "project_candidate_training_exposure": ",".join(exposures),
            "decision": decision,
            "evidence_note": candidate["evidence_note"],
        })

    out_path = WORKSPACE / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    summary = {
        "format_version": 1,
        "identity_definition": "1 - global_levenshtein_distance / max(sequence_lengths)",
        "identity_threshold": args.identity,
        "sensitivity_identity_threshold": sensitivity_threshold,
        "official_cas9_rows": len(official_rows),
        "candidate_count": len(output),
        "decision_counts": dict(sorted(Counter(row["decision"] for row in output).items())),
        "registered_ids": [
            row["candidate_id"] for row in output
            if row["decision"] == "REGISTER_STRICT_PROTEIN_ONLY"
        ],
        "output_tsv": args.out,
    }
    summary_path = WORKSPACE / args.summary_out
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
