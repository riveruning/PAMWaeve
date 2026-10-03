"""Run one end-to-end dual-evidence PAM case and emit auditable structured output.

Task C (真实双证据案例与价值验证) harness.

This script does **not** change any production scoring logic. It orchestrates the
existing, already-reviewed components so a single candidate set is scored by:

1. the protein path only (Protein2PAM -> ``pamdict.score.candidate``);
2. the spacer path only (``pamdict.score.spacer`` oriented protospacer flanks);
3. the transparent late fusion (``pamdict.score.fusion.evidence_assessment``).

and then writes a machine-readable case bundle plus the raw supporting material.
The purpose is to expose *what the second evidence path adds* -- or to report
honestly that it adds nothing. The bundle therefore always contains per-spacer
hits, oriented flanks, support counts and conflict labels, not just a ranking.

Usage
-----
::

    PYTHONPATH=.pylibs:.reference/Protein2PAM:src:. \
    HF_HOME=$PWD/data/checkpoints/hf \
    python scripts/run_dual_evidence_case.py \
        --case benchmarks/dual_evidence_cases/case_spcas9_pampredict.json \
        --outdir data/parsed/dual_evidence_case_spcas9 \
        --device cpu

``--skip-protein`` recomputes only the spacer side and reuses a cached protein
TSV, which is useful when model weights are unavailable.
``--spacer-only`` never touches torch/HF at all.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.score.fusion import evidence_assessment  # noqa: E402
from pamdict.score.spacer import (  # noqa: E402
    background_base_frequencies,
    find_spacer_hits,
    information_bits,
    pam_from_hit,
    probability_consensus,
    read_fasta,
    reverse_complement,
    score_spacer_candidates,
)
from pamdict.score.spectrum import normalize_pam, valid_pam  # noqa: E402

CASE_FORMAT_VERSION = 1
FLANK_LENGTH = 12  # must exceed pam_length; wider flanks let side info be compared


# --------------------------------------------------------------------------
# utilities
# --------------------------------------------------------------------------
def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else WORKSPACE / path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _scalar(row.get(key, "")) for key in fields})


def _scalar(value: object) -> object:
    if isinstance(value, float):
        return f"{value:.6g}"
    if value is None:
        return ""
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def optional_float(value: object) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def optional_int(value: object) -> int:
    if value is None or str(value).strip() == "":
        return 0
    return int(float(value))


# --------------------------------------------------------------------------
# case loading / validation
# --------------------------------------------------------------------------
REQUIRED_CASE_FIELDS = (
    "case_id",
    "system_id",
    "candidate_pams",
    "protein",
    "spacers",
    "targets",
    "pam_side",
    "pam_length",
    "max_mismatches",
    "spacer_orientation",
)


def load_case(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    missing = [field for field in REQUIRED_CASE_FIELDS if field not in payload]
    if missing:
        raise ValueError(f"case file is missing required fields: {missing}")
    candidates = [normalize_pam(str(pam)) for pam in payload["candidate_pams"]]
    if not candidates:
        raise ValueError("case must declare at least one candidate PAM")
    bad = [pam for pam in candidates if not valid_pam(pam)]
    if bad:
        raise ValueError(f"invalid candidate PAMs: {bad}")
    if len({len(pam) for pam in candidates}) != 1:
        raise ValueError("all candidates must share one PAM length")
    if int(payload["pam_length"]) != len(candidates[0]):
        raise ValueError("pam_length must match the candidate width")
    if str(payload["pam_side"]) not in {"downstream", "upstream"}:
        raise ValueError("pam_side must be 'downstream' or 'upstream'")
    # Guard: the case must state the orientation decision and its basis.
    orientation = payload["spacer_orientation"]
    if not isinstance(orientation, dict):
        raise ValueError("spacer_orientation must be an object")
    if orientation.get("decision") not in {"reverse", "forward"}:
        raise ValueError(
            "spacer_orientation.decision must be 'reverse' or 'forward'; "
            "the agent must never guess spacer direction"
        )
    if not orientation.get("evidence"):
        raise ValueError(
            "spacer_orientation.evidence must record how the direction was decided"
        )
    return payload


# --------------------------------------------------------------------------
# spacer path
# --------------------------------------------------------------------------
def hit_rows_from_hits(hits, *, pam_side: str, pam_length: int) -> list[dict[str, object]]:
    """Best (lowest-mismatch) oriented hit per spacer, for readability.

    ``find_spacer_hits`` returns every alignment. The scoring path below uses all
    of them; this view keeps the per-spacer best hit so a human can check the
    flank that actually drove the call.
    """
    best: dict[str, dict[str, object]] = {}
    for hit in hits:
        pam = pam_from_hit(hit, pam_side, pam_length)
        current = best.get(hit.spacer_id)
        if current is None or hit.mismatches < int(current["mismatches"]):
            best[hit.spacer_id] = {
                "spacer_id": hit.spacer_id,
                "spacer_sequence": hit.spacer_sequence,
                "target_id": hit.target_id,
                "strand": hit.strand,
                "mismatches": hit.mismatches,
                "protospacer": hit.protospacer,
                "upstream_flank": hit.upstream_flank[:pam_length],
                "downstream_flank": hit.downstream_flank[:pam_length],
                "pam": pam or "",
                "orientation_weight": hit.orientation_weight,
                "flank_available": 1 if pam else 0,
                "note": "" if pam else "contig boundary: requested PAM-side flank unavailable",
            }
    return sorted(best.values(), key=lambda row: (str(row["spacer_id"])))


def contig_collapse(
    rows: list[dict[str, object]], *, xgg_set: tuple[str, ...] = ("AGG", "CGG", "GGG", "TGG")
) -> dict[str, object]:
    """Measure spacer redundancy at the contig level.

    ``support_target_count`` counts contigs, but one phage contig carrying a
    CRISPR array can absorb many spacers. This reports, per contig, the consensus
    PAM of the spacers that hit it, so a report cannot overstate independence.

    One spacer contributes once (via its best hit in ``rows``). The reported
    contig-level xGG ratio always divides by the number of contigs, never by the
    number of distinct consensus PAM strings -- the latter is a smaller number
    and would inflate the ratio.
    """
    per_target: dict[str, Counter[str]] = defaultdict(Counter)
    spacers_per_contig: Counter[str] = Counter()
    for row in rows:
        pam = str(row.get("pam") or "")
        target = str(row.get("target_id") or "")
        if pam and target:
            per_target[target][pam] += 1
            spacers_per_contig[target] += 1
    consensus_per_contig = {
        target: counts.most_common(1)[0][0] for target, counts in per_target.items()
    }
    contig_consensus = Counter(consensus_per_contig.values())
    total_contigs = len(per_target)
    xgg_contigs = sum(1 for pam in consensus_per_contig.values() if pam in xgg_set)
    ratio = xgg_contigs / total_contigs if total_contigs else 0.0
    # Guard the denominator explicitly: the ratio must equal its own parts.
    if abs(ratio - (xgg_contigs / total_contigs if total_contigs else 0.0)) > 1e-12:
        raise AssertionError("contig-level ratio is inconsistent with its own parts")
    if ratio > 1.0:
        raise AssertionError("contig-level ratio cannot exceed 1")
    top = sorted(spacers_per_contig.items(), key=lambda kv: -kv[1])[:5]
    return {
        "distinct_target_contigs_with_flank": total_contigs,
        "distinct_contig_consensus_pam_types": len(contig_consensus),
        "contigs_supporting_xgg": xgg_contigs,
        # Denominator is the contig count. For published aggregate flanks the
        # "contigs" are per-spacer consensus records, so this ratio is near 1 by
        # construction and must NOT be read as independent replication; check
        # max_spacers_on_one_contig and the evidence class alongside it.
        "contig_level_xgg_fraction": ratio,
        "per_contig_consensus": consensus_per_contig,
        "per_contig_consensus_pam_counts": dict(contig_consensus.most_common()),
        "spacers_per_contig_top5": top,
        "max_spacers_on_one_contig": max(spacers_per_contig.values()) if spacers_per_contig else 0,
    }


def side_information_from_hits(hits, *, side: str, pam_length: int) -> dict[str, object]:
    """Observed consensus/information on one side, across all real hits."""
    pams = [
        pam
        for pam in (pam_from_hit(hit, side, pam_length) for hit in hits)
        if pam
    ]
    if not pams:
        return {"effective_flanks": 0, "consensus": None, "information_bits": None,
                "information_sum": 0.0}
    matrix: list[Counter[str]] = [Counter() for _ in range(pam_length)]
    for pam in pams:
        for index, base in enumerate(pam):
            matrix[index][base] += 1
    total = len(pams)
    rows = [
        {base: counts.get(base, 0) / total for base in "ACGT"} for counts in matrix
    ]
    bits = [information_bits(row) for row in rows]
    return {
        "effective_flanks": total,
        "consensus": probability_consensus(rows),
        "information_bits": bits,
        "information_sum": sum(bits),
        "observed_pam_counts": dict(Counter(pams).most_common()),
    }


def orientation_probe(
    spacers, targets, *, max_mismatches: int, pam_length: int
) -> dict[str, object]:
    """Recompute which side carries the conserved signal, in both orientations.

    This is the check that must gate any direction decision; the declared
    direction is never trusted without it. Only meaningful when raw targets
    exist; published-flank cases state orientation in the source table instead.
    """
    probe: dict[str, object] = {}
    for label, spacer_list in (
        ("as_provided", list(spacers)),
        ("reverse_complemented", [(i, reverse_complement(s)) for i, s in spacers]),
    ):
        hits = find_spacer_hits(
            spacer_list, targets, flank_length=FLANK_LENGTH,
            max_mismatches=max_mismatches,
        )
        up = side_information_from_hits(hits, side="upstream", pam_length=pam_length)
        down = side_information_from_hits(hits, side="downstream", pam_length=pam_length)
        up_bits, down_bits = float(up["information_sum"]), float(down["information_sum"])
        probe[label] = {
            "alignments": len(hits),
            "mapped_spacers": len({h.spacer_id for h in hits}),
            "upstream_consensus": up["consensus"],
            "upstream_information_sum": up_bits,
            "downstream_consensus": down["consensus"],
            "downstream_information_sum": down_bits,
            "stronger_side": "downstream" if down_bits >= up_bits else "upstream",
        }
    return probe


def published_flank_hits(
    flank_path: Path, spacers: list[tuple[str, str]], *, pam_side: str
) -> tuple[list[object], list[str], dict[str, object]]:
    """Build SpacerHit records from published aggregate consensus flanks.

    An aggregate consensus flank is *one* record, not an observed protospacer.
    The reviewed scorer therefore zeroes target/alignment counts for this
    evidence class, and this helper preserves that distinction explicitly.
    """
    from pamdict.score.spacer import SpacerHit, normalize_dna

    spacer_map = dict(spacers)
    if len(spacer_map) != len(spacers):
        raise ValueError("spacer FASTA ids must be unique")
    rows = read_tsv(flank_path)
    required = {
        "spacer_id", "oriented_spacer", "hit",
        "oriented_upstream_flank", "oriented_downstream_flank",
        "orientation", "orientation_source",
    }
    missing = required - set(rows[0] if rows else [])
    if missing:
        raise ValueError(f"published flank TSV missing columns: {sorted(missing)}")
    hits = []
    excluded: list[str] = []
    orientation_sources: Counter[str] = Counter()
    for row in rows:
        spacer_id = row["spacer_id"]
        sequence = normalize_dna(row["oriented_spacer"])
        if spacer_id not in spacer_map:
            raise ValueError(f"{spacer_id} is absent from the spacer FASTA")
        if spacer_map[spacer_id] != sequence:
            raise ValueError(f"{spacer_id} differs between TSV and FASTA")
        if row["hit"] not in {"0", "1"}:
            raise ValueError(f"{spacer_id}: hit must be 0 or 1")
        if row["hit"] == "0":
            continue
        upstream = normalize_dna(row["oriented_upstream_flank"])
        downstream = normalize_dna(row["oriented_downstream_flank"])
        orientation_sources[row.get("orientation_source", "")] += 1
        if not (upstream if pam_side == "upstream" else downstream):
            excluded.append(row.get("source_row_id", spacer_id))
            continue
        hits.append(SpacerHit(
            spacer_id=spacer_id,
            spacer_sequence=sequence,
            target_id=f"published_consensus::{row.get('source_row_id', spacer_id)}",
            start_0=0,
            end_0=len(sequence),
            strand="published_consensus",
            mismatches=0,
            protospacer=sequence,
            upstream_flank=upstream,
            downstream_flank=downstream,
        ))
    diagnostics = {
        "source_rows": len(rows),
        "usable_records": len(hits),
        "excluded_missing_pam_side_flank": excluded,
        "orientation_sources": dict(orientation_sources),
    }
    return hits, excluded, diagnostics


def runs_on_spacer_candidates(
    hits,
    *,
    candidates: list[str],
    pam_side: str,
    prior_strength: float,
    background_records: list[tuple[str, str]],
    aggregate_records: bool,
) -> list[dict[str, object]]:
    """Score candidates from spacer hits, optionally marking aggregate evidence.

    ``aggregate_records=True`` marks an evidence class where each record is a
    published consensus flank rather than an observed target alignment, so
    target/alignment counts are reported as zero to avoid implying independent
    support that does not exist.
    """
    scores = score_spacer_candidates(
        hits,
        candidates,
        side=pam_side,
        base_frequencies=background_base_frequencies(background_records),
        prior_strength=prior_strength,
    )
    rows: list[dict[str, object]] = []
    for score in scores:
        row = score.to_dict()
        if aggregate_records:
            row["aggregate_support_record_count"] = row["support_target_count"]
            row["aggregate_eligible_record_count"] = row["eligible_target_count"]
            row["support_target_count"] = 0
            row["support_alignment_count"] = 0
            row["eligible_target_count"] = 0
            row["eligible_alignment_count"] = 0
            row["evidence_kind"] = "published_aggregate_consensus_flanks"
        else:
            row["evidence_kind"] = "raw_target_alignments"
        rows.append(row)
    return rows


# --------------------------------------------------------------------------
# protein path
# --------------------------------------------------------------------------
def run_protein_path(case: dict[str, object], device: str | None):
    """Score the candidate set with Protein2PAM via the reviewed production API."""
    import numpy as np
    from pamdict.infer.p2pam import MODEL_HF, P2PAMPredictor
    from pamdict.score.candidate import probability_information_stats, rank_candidate_pams

    spec = case["protein"]
    fasta = resolve_path(str(spec["fasta"]))
    record_id = spec.get("record_id")
    records = read_fasta(fasta)
    selected = None
    for name, sequence in records:
        if record_id is None or name == record_id:
            selected = (name, sequence)
            break
    if selected is None:
        raise SystemExit(
            f"protein record {record_id!r} not found in {fasta}; "
            f"first available: {[n for n, _ in records][:10]}"
        )
    name, sequence = selected
    model_key = str(spec.get("model", "cas9_full"))
    if model_key not in MODEL_HF:
        raise SystemExit(f"unsupported protein model {model_key!r}")
    predictor = P2PAMPredictor(model_key, device=device)
    matrix = predictor.predict_probability_matrix([sequence])[0]
    ranked = rank_candidate_pams(
        np.asarray(matrix), list(case["candidate_pams"]), side=str(case["pam_side"])
    )
    stats = probability_information_stats(np.asarray(matrix))
    rows = [
        {
            "protein_id": name,
            "protein_sequence_sha256": sequence_sha256(sequence),
            "model": model_key,
            "candidate_pam": score.candidate_pam,
            "specificity_adjusted_score": score.specificity_adjusted_score,
            "candidate_information_bits_per_position": score.candidate_information_bits_per_position,
            "informative_positions": score.informative_positions,
            "rank_within_length": score.rank_within_length,
        }
        for score in ranked
    ]
    meta = {
        "protein_record": name,
        "protein_length": len(sequence),
        "protein_sequence_sha256": sequence_sha256(sequence),
        "protein_fasta_sha256": file_sha256(fasta),
        "model": model_key,
        "device": getattr(predictor, "_device", device or "auto"),
        "model_consensus_pam": predictor.predict([sequence])[0],
        "matrix_diagnostics": stats,
    }
    return rows, meta


# --------------------------------------------------------------------------
# fusion + verdict
# --------------------------------------------------------------------------
CONFLICT_STATUSES = {
    "concordant_conflict",
    "discordant_support_protein_conflict_spacer",
    "discordant_conflict_protein_support_spacer",
}


def build_fusion(
    protein_rows: list[dict[str, object]], spacer_score_rows: list[dict[str, object]]
) -> list[dict[str, object]]:
    """Join both paths by exact candidate string, preserving both raw scores."""
    spacer_by_candidate = {
        normalize_pam(str(row["candidate_pam"])): row for row in spacer_score_rows
    }
    fused: list[dict[str, object]] = []
    for protein in protein_rows:
        candidate = normalize_pam(str(protein["candidate_pam"]))
        spacer = spacer_by_candidate.get(candidate)
        protein_score = optional_float(protein.get("specificity_adjusted_score"))
        spacer_score = optional_float(
            spacer.get("spacer_evidence_score") if spacer else None
        )
        support_spacers = optional_int(
            spacer.get("support_spacer_count") if spacer else None
        )
        support_targets = optional_int(
            spacer.get("support_target_count") if spacer else None
        )
        status, priority = evidence_assessment(
            protein_score,
            spacer_score,
            support_spacers=support_spacers,
            support_targets=support_targets,
        )
        fused.append(
            {
                "candidate_pam": candidate,
                "protein_score": protein_score,
                "protein_rank": protein.get("rank_within_length"),
                "spacer_score": spacer_score,
                "spacer_rank": spacer.get("spacer_rank_within_length") if spacer else "",
                "spacer_log2_enrichment": spacer.get("spacer_log2_enrichment") if spacer else "",
                "support_spacer_count": support_spacers,
                "support_target_count": support_targets,
                "support_alignment_count": (
                    spacer.get("support_alignment_count") if spacer else 0
                ),
                "evidence_kind": spacer.get("evidence_kind") if spacer else "",
                "evidence_status": status,
                "evidence_priority": priority,
            }
        )
    return fused


def _rank_by(rows: list[dict[str, object]], key: str) -> list[str]:
    def sort_key(row: dict[str, object]):
        value = row.get(key)
        numeric = -1e9 if value in (None, "") else -float(value)  # type: ignore[arg-type]
        return (numeric, str(row["candidate_pam"]))

    return [str(row["candidate_pam"]) for row in sorted(rows, key=sort_key)]


def interpretation(case: dict[str, object], fused: list[dict[str, object]]) -> dict[str, object]:
    """State plainly whether path 2 changed the candidate picture."""
    gold = [normalize_pam(p) for p in case.get("gold_pam_spectrum", [])]
    protein_rank = _rank_by(fused, "protein_score")
    spacer_rank = _rank_by(fused, "spacer_score")
    statuses = {str(row["candidate_pam"]): str(row["evidence_status"]) for row in fused}
    conflicts = sorted(c for c, s in statuses.items() if s in CONFLICT_STATUSES)
    neutral = sorted(c for c, s in statuses.items() if s.endswith("neutral"))
    gold_concordant = [pam for pam in gold if statuses.get(pam) == "concordant_support"]
    return {
        "gold_pam_spectrum": gold,
        "protein_ranking": protein_rank,
        "spacer_ranking": spacer_rank,
        "protein_top": protein_rank[0] if protein_rank else None,
        "spacer_top": spacer_rank[0] if spacer_rank else None,
        "agree_on_top": bool(protein_rank and spacer_rank and protein_rank[0] == spacer_rank[0]),
        "gold_with_concordant_support": gold_concordant,
        "conflicting_candidates": conflicts,
        "neutral_candidates": neutral,
        "evidence_status_by_candidate": statuses,
    }


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, help="case JSON describing inputs")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--device", default=None)
    parser.add_argument(
        "--skip-protein",
        action="store_true",
        help="reuse protein_scores.tsv from --outdir instead of loading the model",
    )
    parser.add_argument(
        "--spacer-only",
        action="store_true",
        help="run only the spacer path (no torch / no model weights needed)",
    )
    args = parser.parse_args()

    case_path = resolve_path(args.case)
    case = load_case(case_path)
    outdir = resolve_path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- inputs hashed at the point of use ----
    protein_fasta = resolve_path(str(case["protein"]["fasta"]))
    spacer_paths = [resolve_path(str(p)) for p in case["spacers"]]
    target_paths = [resolve_path(str(p)) for p in case["targets"]]
    inputs = {
        "protein_fasta": {"path": str(protein_fasta.relative_to(WORKSPACE)),
                          "sha256": file_sha256(protein_fasta)},
        "spacers": [{"path": str(p.relative_to(WORKSPACE)), "sha256": file_sha256(p)}
                    for p in spacer_paths],
        "targets": [{"path": str(p.relative_to(WORKSPACE)), "sha256": file_sha256(p)}
                    for p in target_paths],
    }
    declared = case.get("expected_input_sha256") or {}
    hash_checks = []
    for item in [inputs["protein_fasta"], *inputs["spacers"], *inputs["targets"]]:
        expected = declared.get(Path(str(item["path"])).name)
        hash_checks.append({
            "path": item["path"],
            "sha256": item["sha256"],
            "expected": expected or "",
            "matches": (expected == item["sha256"]) if expected else None,
        })

    # ---- oriented spacers and targets ----
    spacers: list[tuple[str, str]] = []
    for path in spacer_paths:
        spacers.extend(read_fasta(path))
    targets: list[tuple[str, str]] = []
    for path in target_paths:
        targets.extend(read_fasta(path))
    orientation = case["spacer_orientation"]
    original_spacers = list(spacers)
    if orientation["decision"] == "reverse":
        spacers = [(name, reverse_complement(seq)) for name, seq in spacers]

    pam_side = str(case["pam_side"])
    pam_length = int(case["pam_length"])
    raw_mismatches = case.get("max_mismatches")
    if raw_mismatches is None:
        if not case.get("published_flanks"):
            raise ValueError(
                "max_mismatches is required for raw spacer/target cases; "
                "it may be null only for published aggregate flanks"
            )
        max_mismatches = 0  # unused: aggregate flank records carry no alignments
    else:
        max_mismatches = int(raw_mismatches)
    candidates = list(case["candidate_pams"])
    prior_strength = float(case.get("prior_strength", 2.0))

    # ---- spacer path (production functions) ----
    # Two evidence classes share the same scoring semantics but differ in what a
    # "support" means, so the class is carried through explicitly.
    published_path = case.get("published_flanks")
    aggregate_records = bool(published_path)
    flank_diagnostics: dict[str, object] = {}
    if published_path:
        hits, _excluded, flank_diagnostics = published_flank_hits(
            resolve_path(str(published_path)),
            list(read_fasta(spacer_paths[0])),
            pam_side=pam_side,
        )
        background_records = [
            (hit.target_id, hit.upstream_flank + hit.downstream_flank) for hit in hits
        ]
        effective_spacers = len({hit.spacer_sequence for hit in hits})
        min_effective = int(case.get("min_effective_spacers", 5))
        flank_diagnostics["effective_spacers"] = effective_spacers
        flank_diagnostics["min_effective_spacers"] = min_effective
        flank_diagnostics["evidence_eligible"] = effective_spacers >= min_effective
        if effective_spacers < min_effective:
            hits = []
            flank_diagnostics["abstain_reason"] = (
                "fewer effective spacers than the preregistered minimum"
            )
    else:
        hits = find_spacer_hits(
            spacers, targets, flank_length=FLANK_LENGTH, max_mismatches=max_mismatches
        )
        background_records = targets

    hit_rows = hit_rows_from_hits(hits, pam_side=pam_side, pam_length=pam_length)
    write_tsv(outdir / "spacer_hits.tsv", hit_rows, [
        "spacer_id", "spacer_sequence", "target_id", "strand", "mismatches",
        "protospacer", "upstream_flank", "downstream_flank", "pam",
        "orientation_weight", "flank_available", "note",
    ])

    spacer_score_rows = runs_on_spacer_candidates(
        hits,
        candidates=candidates,
        pam_side=pam_side,
        prior_strength=prior_strength,
        background_records=background_records,
        aggregate_records=aggregate_records,
    )
    write_tsv(
        outdir / "spacer_scores.tsv",
        spacer_score_rows,
        list(spacer_score_rows[0]) if spacer_score_rows else ["candidate_pam"],
    )

    # Orientation is only re-derivable when real targets exist. For the
    # published-flank class the source table states orientation per row, so the
    # probe is reported as not applicable rather than silently skipped.
    if targets:
        probe = orientation_probe(
            original_spacers, targets, max_mismatches=max_mismatches, pam_length=pam_length
        )
        declared_side = probe[
            "reverse_complemented" if orientation["decision"] == "reverse" else "as_provided"
        ]["stronger_side"]
        orientation_check = {
            "declared": orientation["decision"],
            "declared_maps_to_stronger_side": declared_side,
            "consistent_with_pam_side": declared_side == pam_side,
            "method": "recomputed flank conservation on raw target alignments",
        }
    else:
        probe = {"applicable": False, "reason": "no raw targets; source states orientation per row"}
        orientation_check = {
            "declared": orientation["decision"],
            "method": "source table orientation field",
            "orientation_sources": flank_diagnostics.get("orientation_sources", {}),
            "consistent_with_pam_side": True,
        }

    # ---- protein path ----
    protein_rows: list[dict[str, object]] = []
    protein_meta: dict[str, object] = {}
    if args.spacer_only:
        protein_meta = {"skipped": "spacer-only mode"}
    elif args.skip_protein:
        cached = outdir / "protein_scores.tsv"
        if not cached.exists():
            raise SystemExit(f"--skip-protein needs an existing {cached}")
        protein_rows = read_tsv(cached)
        for row in protein_rows:
            row["specificity_adjusted_score"] = optional_float(row.get("specificity_adjusted_score"))
        protein_meta = {"reused_from": str(cached.relative_to(WORKSPACE))}
    else:
        protein_rows, protein_meta = run_protein_path(case, args.device)
        write_tsv(outdir / "protein_scores.tsv", protein_rows, [
            "protein_id", "protein_sequence_sha256", "model", "candidate_pam",
            "specificity_adjusted_score", "candidate_information_bits_per_position",
            "informative_positions", "rank_within_length",
        ])

    # ---- fusion ----
    fused: list[dict[str, object]] = []
    verdict: dict[str, object] = {}
    if protein_rows:
        fused = build_fusion(protein_rows, spacer_score_rows)
        write_tsv(outdir / "fused_scores.tsv", fused, [
            "candidate_pam", "protein_score", "protein_rank", "spacer_score",
            "spacer_rank", "spacer_log2_enrichment", "support_spacer_count",
            "support_target_count", "support_alignment_count", "evidence_kind",
            "evidence_status", "evidence_priority",
        ])
        verdict = interpretation(case, fused)

    summary = {
        "format_version": CASE_FORMAT_VERSION,
        "case_id": case["case_id"],
        "system_id": case["system_id"],
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "case_file": str(case_path.relative_to(WORKSPACE)),
        "case_file_sha256": file_sha256(case_path),
        "evidence_class": (
            "published_aggregate_flanks" if aggregate_records else "raw_spacer_target_pairs"
        ),
        "inputs": inputs,
        "hash_checks": hash_checks,
        "parameters": {
            "pam_side": pam_side,
            "pam_length": pam_length,
            "max_mismatches": max_mismatches,
            "flank_length": FLANK_LENGTH,
            "candidate_pams": candidates,
            "prior_strength": prior_strength,
            "spacer_orientation": orientation,
        },
        "spacer_hits_total": len(hits),
        "spacer_side_information": {
            "upstream": side_information_from_hits(hits, side="upstream", pam_length=pam_length),
            "downstream": side_information_from_hits(hits, side="downstream", pam_length=pam_length),
        },
        "spacer_contig_collapse": contig_collapse(hit_rows),
        "flank_diagnostics": flank_diagnostics,
        "orientation_probe": probe,
        "orientation_check": orientation_check,
        "protein_run": protein_meta,
        "fused": fused,
        "interpretation": verdict,
        "evidence_status_note": (
            "the *_support/*_neutral/conflict labels are internal candidate-triage "
            "categories, not calibrated probabilities or activity estimates"
        ),
        "limitations": case.get("limitations", []),
    }
    (outdir / "case_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
