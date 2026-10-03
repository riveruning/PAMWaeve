"""Normalize extracted samples into the unified schema and dedup against baseline.

Produces ``data/corpus/new_candidates.tsv`` with the 10 core columns (empty
logo -> placeholder) plus extension provenance, and a net-new report that
counts how many candidate rows have no protein_sequence/logo collision with
the official baseline.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pamdict.schema.constants import ALL_COLUMNS  # noqa: E402


def ingest(extracted_rows: list[dict], out_tsv: Path) -> dict:
    """Write extracted rows to the schema TSV and return summary stats."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    rows_out = []
    for r in extracted_rows:
        row = {c: "" for c in ALL_COLUMNS}
        row["crispr_type"] = r.get("crispr_type", "")
        row["cas_family"] = r.get("cas_family", "")
        row["source"] = "Literature (new; unreviewed)"
        row["protein_id"] = r.get("protein_id", "")
        row["citation"] = r.get("citation", "")
        row["doi"] = r.get("doi", "")
        row["protein_sequence"] = r.get("protein_sequence", "")
        row["pid_sequence"] = ""
        row["pam_consensus"] = r.get("pam_consensus", "")
        # Empty logo is stored as empty string (schema allows; baseline also
        # has empty-logo rows). A valid 10x4 must be filled during review.
        row["pam_logo_acgt"] = r.get("pam_logo_acgt", "")
        row["raw_ref"] = r.get("raw_ref", "")
        row["created_at"] = r.get("created_at", "")
        row["sha256"] = r.get("sha256", "")
        rows_out.append(row)

    with out_tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ALL_COLUMNS), delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows_out)

    with_consensus = sum(1 for r in rows_out if r["pam_consensus"])
    with_seq = sum(1 for r in rows_out if r["protein_sequence"])
    with_id = sum(1 for r in rows_out if r["protein_id"])
    return {
        "total_rows": len(rows_out),
        "with_consensus": with_consensus,
        "with_protein_sequence": with_seq,
        "with_protein_id": with_id,
    }


def dedup_against_baseline(new_tsv: Path, baseline_tsv: Path) -> dict:
    """Report how many new rows collide with baseline by protein_sequence."""
    def load_seqs(path: Path) -> set[str]:
        seqs: set[str] = set()
        with path.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if header_seqs := row.get("protein_sequence"):
                    seqs.add(header_seqs)
        return seqs

    base = load_seqs(baseline_tsv)
    new_seqs = load_seqs(new_tsv)
    overlap = new_seqs & base
    return {
        "baseline_seqs": len(base),
        "new_seqs": len(new_seqs),
        "overlap": len(overlap),
        "net_new": len(new_seqs) - len(overlap),
    }
