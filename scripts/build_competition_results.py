#!/usr/bin/env python3
"""Build the competition candidate list (``results.csv``) from benchmark output.

This script only *projects* the frozen benchmark artifacts into the competition
field layout.  It never re-ranks, never fuses, and never invents a combined
score:

* ``protein_compatibility_score`` is the protein-side compatibility score,
  renamed for clarity only -- it is **not** cleavage activity, editing
  efficiency, or a success probability;
* ``spacer_evidence_score`` is the spacer-side evidence score;
* ``evidence_status`` is a categorical state label, not a combined score;
* ranks are always ``*_rank_in_length`` because candidates of different PAM
  lengths are not comparable;
* a value that is unavailable is written as an **empty cell**, never ``0``.

The display candidate set (``--display-candidates``) limits which candidates
appear in the list so that it is readable.  It is a *presentation* subset and
must not be reported as the ranking evaluation, which always runs over the
exhaustive concrete A/C/G/T space (see EVALUATION_PROTOCOL.md).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

#: Source label for evidence built from raw spacer-to-target alignments.
SOURCE_RAW_MATCH = "raw_spacer_target_match"
#: Source label for evidence taken from a paper's aggregate per-spacer flanks.
#: These are *not* re-computed alignments and must never be described as such.
SOURCE_PUBLISHED_FLANK = "published_aggregate_flank"

#: The competition track this work is submitted to.  This is the *human-facing*
#: track named in the competition rules, and is deliberately a separate column
#: from ``track``: ``track`` records which input/evidence kind a row came from
#: (protein_only / paired_evidence / published_flanks), which is an internal
#: evaluation concept and must not be mistaken for a competition track.
COMPETITION_TRACK = "赛道二：AI 基因编辑与核酸工具设计"
COMPETITION_TRACK_ID = "track-2-ai-gene-editing-nucleic-acid-tools"

COLUMNS = [
    "candidate_id",
    "competition_track",
    "competition_track_id",
    # ``track`` is retained as the *evidence-type* label it has always been.
    # It is NOT the competition track; see competition_track above.
    "track",
    "track_definition",
    "cas9_id",
    "protein_sequence",
    "protein_sequence_sha256",
    "candidate_pam",
    "pam_length",
    "protein_compatibility_score",
    "protein_rank_in_length",
    "spacer_evidence_score",
    "spacer_rank_in_length",
    "support_spacer_count",
    "support_target_count",
    "evidence_status",
    "ranking_basis",
    "model_version",
    "code_version",
    "code_version_source",
    "code_source_sha256",
    "evidence_source",
    "notes",
]

#: Human-readable meaning of each evidence-type value in ``track``.
TRACK_DEFINITIONS = {
    "protein_only": "evidence type: protein sequence only, no spacer/target channel",
    "paired_evidence": "evidence type: oriented spacers plus raw target sequences",
    "published_flanks": "evidence type: paper-published aggregate per-spacer flanks",
}

#: Notes attached to every row, keyed by evidence source.  Keeping the caveat
#: on the row itself means the limitation travels with the number.
SOURCE_NOTES = {
    SOURCE_RAW_MATCH: (
        "spacer evidence from raw spacer-to-target alignment; "
        "target contig counts are not independent virus counts"
    ),
    SOURCE_PUBLISHED_FLANK: (
        "spacer evidence from published per-spacer aggregate flanks; "
        "raw target identifiers were never published, so support_target_count "
        "is unavailable (blank), not zero"
    ),
}

NOT_APPLICABLE_NOTES = (
    "this system provides no spacer evidence track; spacer fields are "
    "unavailable (blank), not zero"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta_single(path: Path) -> tuple[str, str]:
    name = None
    parts: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                break
            name = line[1:].split()[0]
        elif name is not None:
            parts.append("".join(line.split()))
    if name is None or not parts:
        raise ValueError(f"no FASTA record in {path}")
    return name, "".join(parts).upper()


#: Files that constitute "the analysis code" for version-recording purposes.
#: Hashing them gives a version identity that survives being unpacked without
#: Git (a plain zip has no ``.git``), so a reviewer can always tell which source
#: produced a given ``results.csv``.
CODE_IDENTITY_FILES = [
    "scripts/run_multi_evidence_benchmark.py",
    "scripts/build_competition_results.py",
    "scripts/score_candidate_pams.py",
    "scripts/score_spacer_pams.py",
    "scripts/score_published_flanks.py",
    "scripts/fuse_pam_evidence.py",
    "src/pamdict/benchmark/multi_evidence.py",
    "src/pamdict/score/candidate.py",
    "src/pamdict/score/spacer.py",
    "src/pamdict/score/fusion.py",
    "src/pamdict/infer/p2pam.py",
]


def code_version() -> tuple[str, str, str | None]:
    """Return ``(version_label, source, source_sha256)`` for the analysis code.

    ``git rev-parse HEAD`` alone is not enough: it fails (or, worse, reports a
    *different* repository's commit) once the package is unpacked outside its
    Git checkout, and it does not capture uncommitted edits.  So we always
    compute a SHA-256 over the analysis source files, and record Git HEAD only
    as supporting information when it genuinely refers to this tree.
    """
    digest = hashlib.sha256()
    hashed: list[str] = []
    for relative in sorted(CODE_IDENTITY_FILES):
        path = WORKSPACE / relative
        if not path.is_file():
            continue
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
        hashed.append(relative)
    source_sha256 = digest.hexdigest() if hashed else None

    git_head: str | None = None
    try:
        probe = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            cwd=WORKSPACE,
            capture_output=True,
            text=True,
            check=True,
        )
        # Only trust HEAD when the enclosing repo really is this tree; a zip
        # unpacked inside some unrelated repository must not adopt its commit.
        if Path(probe.stdout.strip()).resolve() == WORKSPACE.resolve():
            head = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=WORKSPACE,
                capture_output=True,
                text=True,
                check=True,
            )
            git_head = head.stdout.strip()
    except Exception:  # noqa: BLE001 - a missing git must not break the report
        git_head = None

    if git_head:
        label = f"{git_head}+src.{source_sha256[:12]}" if source_sha256 else git_head
        source = "git HEAD of this tree, plus SHA-256 over analysis sources"
    elif source_sha256:
        label = f"src.{source_sha256[:12]}"
        source = (
            "SHA-256 over analysis sources (no Git checkout enclosing this "
            "package, so no commit is claimed)"
        )
    else:
        label = "unknown"
        source = "no Git commit and no analysis sources found"
    return label, source, source_sha256


def blank(value: object) -> str:
    """Return an empty cell for unavailable values, never ``0``."""
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text == "" else text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-dir",
        default="data/parsed/bench_final22",
        help="directory written by run_multi_evidence_benchmark.py",
    )
    parser.add_argument(
        "--manifest",
        default="benchmarks/multi_evidence_v2/manifest_final22.json",
    )
    parser.add_argument("--out", default="results.csv")
    parser.add_argument(
        "--display-candidates",
        default=None,
        help=(
            "JSON file with {system_id: [PAM, ...]} restricting which "
            "candidates are listed. Presentation subset only; the ranking "
            "evaluation always uses the exhaustive candidate space."
        ),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=0,
        help=(
            "also list the top N protein-ranked candidates per system. A "
            "deterministic presentation rule that keeps the list readable "
            "without hand-picking winners, and that guarantees the model's own "
            "top candidate is present. 0 disables it."
        ),
    )
    parser.add_argument(
        "--systems",
        default=None,
        help=(
            "comma-separated system_ids to include as candidate rows. Defaults "
            "to every evaluated system. Use this to restrict the list to the "
            "frozen showcase cases."
        ),
    )
    args = parser.parse_args()

    bench_dir = WORKSPACE / args.benchmark_dir
    report_path = bench_dir / "benchmark_report.json"
    if not report_path.exists():
        raise SystemExit(
            f"benchmark report not found: {report_path}\n"
            "run scripts/run_multi_evidence_benchmark.py first"
        )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    manifest = json.loads((WORKSPACE / args.manifest).read_text(encoding="utf-8"))
    manifest_by_id = {s["system_id"]: s for s in manifest["systems"]}

    display: dict[str, list[str]] = {}
    if args.display_candidates:
        raw = json.loads(
            (WORKSPACE / args.display_candidates).read_text(encoding="utf-8")
        )
        display = {
            key: [str(v).upper() for v in value]
            for key, value in raw.items()
            if not key.startswith("_")
        }

    version, version_source, source_sha256 = code_version()
    rows: list[dict[str, str]] = []
    skipped: list[dict[str, str]] = []
    selected = (
        {value.strip() for value in args.systems.split(",") if value.strip()}
        if args.systems
        else None
    )

    for system in report["systems"]:
        system_id = system["system_id"]
        if selected is not None and system_id not in selected:
            continue
        entry = manifest_by_id.get(system_id, {})
        system_dir = bench_dir / system_id
        protein_path = system_dir / "protein.fasta"
        if not protein_path.exists():
            skipped.append({
                "system_id": system_id,
                "reason": f"missing protein FASTA: {protein_path}",
            })
            continue
        _record_id, sequence = read_fasta_single(protein_path)

        track = system.get("spacer_evidence_track") or "protein_only"
        source = (
            SOURCE_RAW_MATCH
            if track == "paired_evidence"
            else SOURCE_PUBLISHED_FLANK
            if track == "published_flanks"
            else ""
        )

        protein_rows = read_tsv(system_dir / "protein_scores.tsv")
        spacer_rows = (
            read_tsv(system_dir / "spacer_scores.tsv")
            if (system_dir / "spacer_scores.tsv").exists()
            and system["methods"].get("spacer_only", {}).get("not_applicable")
            is not True
            else []
        )
        fusion_rows = (
            read_tsv(system_dir / "fusion_scores.tsv")
            if (system_dir / "fusion_scores.tsv").exists()
            else []
        )
        spacer_by_pam = {row["candidate_pam"]: row for row in spacer_rows}
        fusion_by_pam = {row["candidate_pam"]: row for row in fusion_rows}

        allowed = display.get(system_id)
        if allowed is not None:
            allowed = set(allowed)
        if args.top_n > 0:
            # Union the explicit presentation list with the model's own top-N
            # protein candidates.  This is a fixed rule applied to every
            # system, not a per-system choice of convenient winners, and it
            # guarantees the model's top-ranked candidate always appears.
            ranked = sorted(
                protein_rows,
                key=lambda row: int(row["rank_within_length"]),
            )
            top = {row["candidate_pam"] for row in ranked[: args.top_n]}
            allowed = top if allowed is None else allowed | top
        index = 0
        for protein_row in protein_rows:
            pam = protein_row["candidate_pam"]
            if allowed is not None and pam not in allowed:
                continue
            index += 1
            spacer_row = spacer_by_pam.get(pam)
            fusion_row = fusion_by_pam.get(pam)

            if spacer_row is None:
                evidence_status = ""
                ranking_basis = (
                    "protein_only: ranked by protein compatibility within "
                    "PAM length"
                    if track == "protein_only"
                    else "spacer evidence unavailable for this candidate"
                )
            else:
                evidence_status = blank(fusion_row.get("evidence_status")) if fusion_row else ""
                gate = blank(fusion_row.get("fusion_gate")) if fusion_row else ""
                ranking_basis = (
                    f"fusion gate={gate}" if gate else "fusion gate not recorded"
                )

            note = SOURCE_NOTES.get(source, "") or NOT_APPLICABLE_NOTES
            if track == "protein_only":
                note = NOT_APPLICABLE_NOTES

            rows.append({
                "candidate_id": f"{system_id}-{index:04d}",
                "competition_track": COMPETITION_TRACK,
                "competition_track_id": COMPETITION_TRACK_ID,
                # Evidence-type label, NOT the competition track.
                "track": track,
                "track_definition": TRACK_DEFINITIONS.get(track, ""),
                "cas9_id": str(entry.get("protein_record_id") or protein_row["protein_id"]),
                "protein_sequence": sequence,
                "protein_sequence_sha256": blank(protein_row.get("sequence_sha256")),
                "candidate_pam": pam,
                "pam_length": blank(protein_row.get("pam_length")),
                "protein_compatibility_score": blank(
                    protein_row.get("specificity_adjusted_score")
                ),
                "protein_rank_in_length": blank(
                    protein_row.get("rank_within_length")
                ),
                "spacer_evidence_score": (
                    blank(spacer_row.get("spacer_evidence_score"))
                    if spacer_row else ""
                ),
                "spacer_rank_in_length": (
                    blank(spacer_row.get("spacer_rank_within_length"))
                    if spacer_row else ""
                ),
                "support_spacer_count": (
                    blank(spacer_row.get("support_spacer_count"))
                    if spacer_row else ""
                ),
                # For aggregate flanks this stays blank: the runner has already
                # emptied it, because no raw target identifiers exist.
                "support_target_count": (
                    blank(spacer_row.get("support_target_count"))
                    if spacer_row else ""
                ),
                "evidence_status": evidence_status,
                "ranking_basis": ranking_basis,
                "model_version": blank(protein_row.get("model")),
                "code_version": version,
                "code_version_source": version_source,
                "code_source_sha256": source_sha256 or "",
                "evidence_source": source,
                "notes": note,
            })

    out_path = WORKSPACE / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "rows": len(rows),
        "systems_with_rows": len({row["candidate_id"].rsplit("-", 1)[0] for row in rows}),
        "systems_skipped": skipped,
        "display_subset": bool(display),
        "source_counts": {
            label: sum(1 for row in rows if row["evidence_source"] == label)
            for label in (SOURCE_RAW_MATCH, SOURCE_PUBLISHED_FLANK, "")
        },
        "benchmark_manifest_sha256": report.get("manifest_sha256"),
        "benchmark_generated_at": report.get("generated_at"),
        "decision": report.get("decision"),
        "code_version": version,
        "code_version_source": version_source,
        "code_source_sha256": source_sha256,
        "code_identity_files": sorted(CODE_IDENTITY_FILES),
        "field_notes": {
            "competition_track": (
                "competition track this work is submitted to; distinct from "
                "'track', which is the internal evidence-type label"
            ),
            "track": (
                "evidence type (protein_only / paired_evidence / "
                "published_flanks), NOT the competition track"
            ),
            "protein_compatibility_score": (
                "protein-side model compatibility; NOT cleavage activity, "
                "editing efficiency, or a success probability"
            ),
            "evidence_status": (
                "categorical evidence state; no calibrated combined score exists"
            ),
            "missing_values": "unavailable values are empty cells, never 0",
        },
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"wrote {len(rows)} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
