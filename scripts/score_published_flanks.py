#!/usr/bin/env python3
"""Score candidate PAMs from published per-spacer consensus flanks."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(WORKSPACE / "src"))

from scripts.score_spacer_pams import (  # noqa: E402
    file_sha256,
    load_candidates,
    resolve_path,
    write_tsv,
)
from pamdict.score.spacer import (  # noqa: E402
    SpacerHit,
    background_base_frequencies,
    flank_probability_matrix,
    information_bits,
    normalize_dna,
    probability_consensus,
    read_fasta,
    score_spacer_candidates,
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "spacer_id",
        "source_row_id",
        "oriented_spacer",
        "hit",
        "oriented_upstream_flank",
        "oriented_downstream_flank",
        "orientation",
        "orientation_source",
    }
    missing = required - set(rows[0] if rows else [])
    if missing:
        raise ValueError(f"published flank TSV missing columns: {sorted(missing)}")
    return rows


def hits_from_rows(
    rows: list[dict[str, str]],
    spacers: dict[str, str],
    required_side: str,
) -> tuple[list[SpacerHit], list[str]]:
    hits: list[SpacerHit] = []
    excluded_incomplete: list[str] = []
    for row in rows:
        spacer_id = row["spacer_id"]
        sequence = normalize_dna(row["oriented_spacer"])
        if spacer_id not in spacers:
            raise ValueError(f"{spacer_id} is absent from spacer FASTA")
        if spacers[spacer_id] != sequence:
            raise ValueError(f"{spacer_id} sequence differs between TSV and FASTA")
        if row["hit"] not in {"0", "1"}:
            raise ValueError(f"{spacer_id} hit must be 0 or 1")
        if row["hit"] == "0":
            continue
        upstream = normalize_dna(row["oriented_upstream_flank"])
        downstream = normalize_dna(row["oriented_downstream_flank"])
        required_flank = upstream if required_side == "upstream" else downstream
        if not required_flank:
            # Some source-table rows are marked hit=1 even though the
            # consensus flank on the declared PAM side is empty.  They carry
            # no scoreable PAM evidence, so exclude them explicitly while
            # retaining their source-row ids in the summary audit trail.
            excluded_incomplete.append(row["source_row_id"])
            continue
        hits.append(SpacerHit(
            spacer_id=spacer_id,
            spacer_sequence=sequence,
            target_id=f"published_consensus::{row['source_row_id']}",
            start_0=0,
            end_0=len(sequence),
            strand="published_consensus",
            mismatches=0,
            protospacer=sequence,
            upstream_flank=upstream,
            downstream_flank=downstream,
        ))
    return hits, excluded_incomplete


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Rank candidate PAMs from aggregate flank consensus records. "
            "This does not reconstruct raw protospacer hits."
        )
    )
    parser.add_argument("--spacers", required=True)
    parser.add_argument("--published-flanks", required=True)
    parser.add_argument("--system-id", required=True)
    parser.add_argument(
        "--pam-side", choices=("upstream", "downstream"), required=True
    )
    parser.add_argument("--pam-length", type=int, required=True)
    parser.add_argument("--candidates")
    parser.add_argument("--candidates-file")
    parser.add_argument("--prior-strength", type=float, default=2.0)
    parser.add_argument("--min-effective-spacers", type=int, default=5)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    if args.pam_length <= 0:
        raise ValueError("--pam-length must be positive")
    if args.prior_strength < 0:
        raise ValueError("--prior-strength must be non-negative")
    if args.min_effective_spacers <= 0:
        raise ValueError("--min-effective-spacers must be positive")

    spacer_path = resolve_path(args.spacers)
    flank_path = resolve_path(args.published_flanks)
    candidate_path = (
        resolve_path(args.candidates_file) if args.candidates_file else None
    )
    outdir = resolve_path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    spacer_records = read_fasta(spacer_path)
    spacers = dict(spacer_records)
    if len(spacers) != len(spacer_records):
        raise ValueError("spacer FASTA ids must be unique")
    rows = read_rows(flank_path)
    hits, excluded_incomplete = hits_from_rows(rows, spacers, args.pam_side)
    if not hits:
        raise ValueError("no published flank records are available")

    candidates = load_candidates(args.candidates, candidate_path)
    if not candidates:
        raise ValueError("published-flank scoring requires explicit candidates")
    if any(len(candidate) != args.pam_length for candidate in candidates):
        raise ValueError("all candidates must match --pam-length")
    background_records = [
        (
            hit.target_id,
            hit.upstream_flank + hit.downstream_flank,
        )
        for hit in hits
    ]
    base_frequencies = background_base_frequencies(background_records)
    effective_spacers = len({hit.spacer_sequence for hit in hits})
    evidence_eligible = effective_spacers >= args.min_effective_spacers
    scoring_hits = hits if evidence_eligible else []
    scored = score_spacer_candidates(
        scoring_hits,
        candidates,
        side=args.pam_side,
        base_frequencies=base_frequencies,
        prior_strength=args.prior_strength,
    )
    score_rows = []
    for score in scored:
        raw_score = score.to_dict()
        score_rows.append({
            "system_id": args.system_id,
            "pam_side": args.pam_side,
            "max_mismatches": "not_applicable",
            "evidence_status": (
                "scored" if evidence_eligible
                else "abstain_low_effective_spacers"
            ),
            **raw_score,
            "aggregate_support_record_count": raw_score["support_target_count"],
            "aggregate_eligible_record_count": raw_score["eligible_target_count"],
            # A published consensus flank is one aggregate record, not an
            # independently observed protospacer target/alignment.
            "support_target_count": 0,
            "support_alignment_count": 0,
            "eligible_target_count": 0,
            "eligible_alignment_count": 0,
            "interpretation": (
                "published aggregate spacer-flank evidence; not raw target "
                "alignments, cleavage activity, or experimental probability"
            ),
        })
    write_tsv(
        outdir / "candidate_scores.tsv",
        list(score_rows[0]),
        score_rows,
    )
    used_rows = [
        {
            "system_id": args.system_id,
            **hit.to_dict(),
            "evidence_kind": "published_per_spacer_consensus_flanks",
        }
        for hit in hits
    ]
    write_tsv(
        outdir / "published_flanks_used.tsv",
        list(used_rows[0]),
        used_rows,
    )

    side_summaries: dict[str, dict[str, object]] = {}
    flank_rows: list[dict[str, object]] = []
    for side in ("upstream", "downstream"):
        matrix, effective = flank_probability_matrix(
            hits, side=side, length=args.pam_length
        )
        infos = [information_bits(row) for row in matrix]
        side_summaries[side] = {
            "effective_spacers": effective,
            "consensus": probability_consensus(matrix),
            "information_bits": infos,
            "information_sum": sum(infos),
            "probability_matrix": matrix,
        }
        for index, (probabilities, info) in enumerate(zip(matrix, infos)):
            position = (
                index + 1 if side == "downstream"
                else index - args.pam_length
            )
            flank_rows.append({
                "system_id": args.system_id,
                "side": side,
                "position": position,
                **{
                    base: round(probabilities[base], 8)
                    for base in "ACGT"
                },
                "information_bits": round(info, 8),
                "effective_spacers": effective,
            })
    write_tsv(
        outdir / "flank_probabilities.tsv",
        [
            "system_id", "side", "position", "A", "C", "G", "T",
            "information_bits", "effective_spacers",
        ],
        flank_rows,
    )

    inferred_signal_side = max(
        side_summaries,
        key=lambda side: float(side_summaries[side]["information_sum"]),
    )
    orientation_warning = (
        inferred_signal_side != args.pam_side
        and float(side_summaries[inferred_signal_side]["information_sum"])
        > float(side_summaries[args.pam_side]["information_sum"]) * 1.25
    )
    summary = {
        "format_version": 1,
        "system_id": args.system_id,
        "evidence_track": "published_flanks",
        "inputs": {
            "spacers": str(spacer_path),
            "spacers_sha256": file_sha256(spacer_path),
            "published_flanks": str(flank_path),
            "published_flanks_sha256": file_sha256(flank_path),
        },
        "parameters": {
            "pam_side": args.pam_side,
            "pam_length": args.pam_length,
            "prior_strength": args.prior_strength,
            "min_effective_spacers": args.min_effective_spacers,
            "background": "canonical bases in published consensus flanks",
        },
        "counts": {
            "input_spacers": len(spacers),
            "unique_spacers": len(set(spacers.values())),
            "effective_spacers_for_pam_side": effective_spacers,
            "evidence_eligible": evidence_eligible,
            "mapped_unique_spacers": len({hit.spacer_sequence for hit in hits}),
            "published_flank_records": len(hits),
            "incomplete_required_flank_records": len(excluded_incomplete),
            "excluded_incomplete_source_row_ids": excluded_incomplete,
            "raw_target_records": 0,
            "alignments": 0,
        },
        "target_base_frequencies": base_frequencies,
        "side_summaries": side_summaries,
        "inferred_signal_side": inferred_signal_side,
        "orientation_warning": orientation_warning,
        "candidate_scores": [score.to_dict() for score in scored],
        "limitations": [
            "Flanks are published per-spacer aggregate consensus sequences.",
            "Raw target identifiers, per-hit alignments, and virus clusters are unavailable.",
            "One consensus record is not an independent target count.",
            "Spacer evidence is evolutionary evidence, not direct cleavage activity.",
            "Source hit rows with an empty flank on the declared PAM side are excluded and listed in counts.",
        ],
    }
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "system_id": args.system_id,
        "input_spacers": len(spacers),
        "published_flank_records": len(hits),
        "incomplete_required_flank_records": len(excluded_incomplete),
        "evidence_eligible": evidence_eligible,
        "orientation_warning": orientation_warning,
        "outdir": str(outdir),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
