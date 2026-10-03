"""Map CRISPR spacers to target genomes and score candidate PAMs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.score.spacer import (  # noqa: E402
    background_base_frequencies,
    find_spacer_hits,
    flank_probability_matrix,
    information_bits,
    pam_from_hit,
    probability_consensus,
    read_fasta,
    reverse_complement,
    score_spacer_candidates,
)
from pamdict.score.spectrum import normalize_pam, valid_pam  # noqa: E402


def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else WORKSPACE / path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_candidates(inline: str | None, path: Path | None) -> list[str]:
    values: list[str] = []
    if inline:
        values.extend(value.strip() for value in inline.split(","))
    if path:
        lines = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if lines and "\t" in lines[0]:
            reader = csv.DictReader(lines, delimiter="\t")
            fields = reader.fieldnames or []
            column = next(
                (
                    name for name in ("candidate_pam", "pam", "pam_consensus")
                    if name in fields
                ),
                None,
            )
            if column is None:
                raise ValueError(
                    "candidate TSV needs candidate_pam, pam, or pam_consensus"
                )
            values.extend((row.get(column) or "").strip() for row in reader)
        else:
            values.extend(lines)
    motifs = list(dict.fromkeys(
        normalize_pam(value) for value in values if value.strip()
    ))
    invalid = [motif for motif in motifs if not valid_pam(motif)]
    if invalid:
        raise ValueError(f"invalid IUPAC PAMs: {invalid}")
    return motifs


def load_targets(paths: list[Path]) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    used: set[str] = set()
    for path in paths:
        for original_id, sequence in read_fasta(path):
            target_id = original_id
            if target_id in used:
                target_id = f"{path.stem}:{original_id}"
            suffix = 2
            base = target_id
            while target_id in used:
                target_id = f"{base}:{suffix}"
                suffix += 1
            used.add(target_id)
            records.append((target_id, sequence))
    return records


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Find ungapped spacer/protospacer matches on both strands and score "
            "candidate IUPAC PAMs. Each unique spacer has total weight one."
        )
    )
    parser.add_argument("--spacers", required=True, help="oriented spacer FASTA")
    parser.add_argument(
        "--targets",
        required=True,
        action="append",
        help="phage/plasmid/genome FASTA; may be supplied more than once",
    )
    parser.add_argument("--system-id", default="spacer-system")
    parser.add_argument(
        "--pam-side",
        choices=("downstream", "upstream"),
        required=True,
        help="PAM side relative to the oriented spacer",
    )
    parser.add_argument("--pam-length", type=int, default=3)
    parser.add_argument("--max-mismatches", type=int, default=2)
    parser.add_argument(
        "--reverse-spacers",
        action="store_true",
        help="reverse-complement every input spacer before mapping",
    )
    parser.add_argument("--candidates", help="comma-separated IUPAC PAMs")
    parser.add_argument("--candidates-file")
    parser.add_argument("--prior-strength", type=float, default=2.0)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    if args.pam_length <= 0:
        raise ValueError("--pam-length must be positive")
    if args.max_mismatches < 0:
        raise ValueError("--max-mismatches must be non-negative")
    if args.prior_strength < 0:
        raise ValueError("--prior-strength must be non-negative")

    spacer_path = resolve_path(args.spacers)
    target_paths = [resolve_path(value) for value in args.targets]
    candidate_path = (
        resolve_path(args.candidates_file) if args.candidates_file else None
    )
    outdir = resolve_path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_spacers = read_fasta(spacer_path)
    spacers = [
        (name, reverse_complement(sequence) if args.reverse_spacers else sequence)
        for name, sequence in input_spacers
    ]
    targets = load_targets(target_paths)
    unique_spacer_count = len({sequence for _name, sequence in spacers})
    print(
        f"spacers={len(spacers)} unique={unique_spacer_count} "
        f"targets={len(targets)} max_mismatches={args.max_mismatches}",
        file=sys.stderr,
        flush=True,
    )

    hits = find_spacer_hits(
        spacers,
        targets,
        flank_length=args.pam_length,
        max_mismatches=args.max_mismatches,
    )
    print(f"alignments={len(hits)}", file=sys.stderr, flush=True)

    candidates = load_candidates(args.candidates, candidate_path)
    if candidates and any(len(motif) != args.pam_length for motif in candidates):
        raise ValueError(
            "all candidate PAMs must have the --pam-length value; run each "
            "length separately"
        )
    if not candidates:
        candidates = sorted({
            pam
            for hit in hits
            if (pam := pam_from_hit(hit, args.pam_side, args.pam_length))
            is not None
        })
    if not candidates:
        raise ValueError(
            "no candidate PAMs were supplied and no full-length flanks were found"
        )

    base_frequencies = background_base_frequencies(targets)
    scored = score_spacer_candidates(
        hits,
        candidates,
        side=args.pam_side,
        base_frequencies=base_frequencies,
        prior_strength=args.prior_strength,
    )

    match_rows = []
    for hit in hits:
        match_rows.append({
            "system_id": args.system_id,
            **hit.to_dict(),
            "upstream_pam": pam_from_hit(hit, "upstream", args.pam_length) or "",
            "downstream_pam": pam_from_hit(hit, "downstream", args.pam_length) or "",
        })
    match_fields = [
        "system_id",
        "spacer_id",
        "spacer_sequence",
        "target_id",
        "start_0",
        "end_0",
        "strand",
        "mismatches",
        "protospacer",
        "upstream_flank",
        "downstream_flank",
        "orientation_weight",
        "upstream_pam",
        "downstream_pam",
    ]
    write_tsv(outdir / "matches.tsv", match_fields, match_rows)

    score_rows = [
        {
            "system_id": args.system_id,
            "pam_side": args.pam_side,
            "max_mismatches": args.max_mismatches,
            **score.to_dict(),
            "interpretation": (
                "spacer-flank evidence; not cleavage activity or experimental probability"
            ),
        }
        for score in scored
    ]
    score_fields = list(score_rows[0])
    write_tsv(outdir / "candidate_scores.tsv", score_fields, score_rows)

    matrices: dict[str, list[dict[str, float]]] = {}
    side_summaries: dict[str, dict[str, object]] = {}
    flank_rows: list[dict[str, object]] = []
    for side in ("upstream", "downstream"):
        matrix, effective = flank_probability_matrix(
            hits, side=side, length=args.pam_length
        )
        matrices[side] = matrix
        infos = [information_bits(row) for row in matrix]
        side_summaries[side] = {
            "effective_spacers": effective,
            "consensus": probability_consensus(matrix),
            "information_bits": infos,
            "information_sum": sum(infos),
            "probability_matrix": matrix,
        }
        for index, (row, info) in enumerate(zip(matrix, infos)):
            position = (
                index + 1 if side == "downstream" else index - args.pam_length
            )
            flank_rows.append({
                "system_id": args.system_id,
                "side": side,
                "position": position,
                **{base: round(row[base], 8) for base in "ACGT"},
                "information_bits": round(info, 8),
                "effective_spacers": effective,
            })
    write_tsv(
        outdir / "flank_probabilities.tsv",
        [
            "system_id",
            "side",
            "position",
            "A",
            "C",
            "G",
            "T",
            "information_bits",
            "effective_spacers",
        ],
        flank_rows,
    )

    mapped_spacers = {hit.spacer_sequence for hit in hits}
    inferred_signal_side = max(
        side_summaries,
        key=lambda side: float(side_summaries[side]["information_sum"]),
    )
    summary = {
        "format_version": 1,
        "system_id": args.system_id,
        "inputs": {
            "spacers": str(spacer_path),
            "spacers_sha256": file_sha256(spacer_path),
            "targets": [str(path) for path in target_paths],
            "targets_sha256": {
                str(path): file_sha256(path) for path in target_paths
            },
        },
        "parameters": {
            "pam_side": args.pam_side,
            "pam_length": args.pam_length,
            "max_mismatches": args.max_mismatches,
            "reverse_spacers": args.reverse_spacers,
            "prior_strength": args.prior_strength,
            "matcher": "ungapped Hamming distance, both target strands",
        },
        "counts": {
            "input_spacers": len(spacers),
            "unique_spacers": unique_spacer_count,
            "mapped_unique_spacers": len(mapped_spacers),
            "target_records": len(targets),
            "alignments": len(hits),
        },
        "target_base_frequencies": base_frequencies,
        "side_summaries": side_summaries,
        "inferred_signal_side": inferred_signal_side,
        "orientation_warning": (
            inferred_signal_side != args.pam_side
            and float(side_summaries[inferred_signal_side]["information_sum"])
            > float(side_summaries[args.pam_side]["information_sum"]) * 1.25
        ),
        "candidate_scores": [score.to_dict() for score in scored],
        "limitations": [
            "ungapped matching only; use BLAST-like alignment for indels or large databases",
            "target contig count is not automatically an independent-phage count",
            "spacer evidence is evolutionary evidence, not direct cleavage activity",
            "array orientation and Cas-array linkage must be checked from genome context",
        ],
    }
    (outdir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"outputs -> {outdir}", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
