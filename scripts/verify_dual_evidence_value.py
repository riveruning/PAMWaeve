"""Reproduce the task-C value assessment for an existing dual-evidence case.

This is the *verification* half of task C (真实双证据案例与价值验证). The main
harness ``scripts/run_dual_evidence_case.py`` shows *what* the two evidence paths
report. This script asks the harder question: does the spacer path carry real
discriminating information, or would any spacer library have produced the same
"NGG" answer?

It answers that with controls that preserve the things that could fake a signal:

* **wrong orientation** -- the same spacers, reverse-complemented. Alignments
  still exist, so this isolates *orientation* from *homology*.
* **shuffled spacers** -- destroys spacer/protospacer homology, keeps length.
* **shuffled target contigs** -- keeps each library's base composition exactly,
  destroys homology.
* **random spacers** -- matched length, random sequence.

A real signal must survive the positive case and collapse in every control. It
also reports the contig-collapse view, because ``support_target_count`` counts
contigs and one phage contig carrying a CRISPR array can absorb many spacers.

Reads only existing case artifacts and the declared inputs; writes nothing except
its own JSON report. It does not modify any production scoring logic.

Usage
-----
::

    PYTHONPATH=src:. python scripts/verify_dual_evidence_value.py \
        --case benchmarks/dual_evidence_cases/case_spcas9_pampredict_paired.json \
        --out data/parsed/dual_evidence_case_spcas9/value_check.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.score.spacer import (  # noqa: E402
    background_base_frequencies,
    find_spacer_hits,
    pam_from_hit,
    read_fasta,
    reverse_complement,
)

SEED = 20261001

#: Concrete 3-mers that satisfy the xGG / NGG PAM family.
XGG_SET = ("AGG", "CGG", "GGG", "TGG")


def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else WORKSPACE / path


def canonical_xgg(pams: list[str]) -> float:
    """Fraction of observed concrete PAMs that match xGG."""
    if not pams:
        return 0.0
    counts = Counter(pams)
    return sum(counts.get(pam, 0) for pam in XGG_SET) / len(pams)


def contig_level_xgg(best_per_spacer: dict[str, tuple[int, str, str]]) -> dict[str, object]:
    """Collapse spacers to contigs, then ask how many *contigs* support xGG.

    The denominator is the number of contigs that produced a usable best hit.
    It must never be the number of distinct consensus PAM strings: several
    contigs can share one consensus, so that denominator is smaller and inflates
    the fraction. This function exists separately so the ratio can be tested
    against its own numerator and denominator.
    """
    per_contig: dict[str, list[str]] = {}
    for _mismatches, target_id, pam in best_per_spacer.values():
        per_contig.setdefault(target_id, []).append(pam)
    consensus: dict[str, str] = {
        target_id: Counter(pams_on_contig).most_common(1)[0][0]
        for target_id, pams_on_contig in per_contig.items()
    }
    xgg_contigs = sum(1 for pam in consensus.values() if pam in XGG_SET)
    total_contigs = len(per_contig)
    return {
        "distinct_target_contigs": total_contigs,
        "distinct_contig_consensus_pam_types": len(set(consensus.values())),
        "contigs_consensus_xgg": xgg_contigs,
        "contig_level_xgg_fraction": (
            xgg_contigs / total_contigs if total_contigs else 0.0
        ),
        "contig_consensus_counts": dict(Counter(consensus.values()).most_common()),
        "max_spacers_on_one_contig": (
            max(len(v) for v in per_contig.values()) if per_contig else 0
        ),
    }


def evaluate(
    spacers: list[tuple[str, str]],
    targets: list[tuple[str, str]],
    *,
    pam_side: str,
    pam_length: int,
    max_mismatches: int,
) -> dict[str, object]:
    hits = find_spacer_hits(
        spacers, targets, flank_length=10, max_mismatches=max_mismatches
    )
    pams = [
        pam
        for pam in (pam_from_hit(hit, pam_side, pam_length) for hit in hits)
        if pam
    ]
    counts = Counter(pams)
    # Contig-level view must use ONE PAM per spacer (its best hit). Using every
    # alignment would let a single spacer contribute several times to the same
    # contig and make the "fraction" exceed 1.
    best_per_spacer: dict[str, tuple[int, str, str]] = {}
    for hit in hits:
        pam = pam_from_hit(hit, pam_side, pam_length)
        if not pam:
            continue
        current = best_per_spacer.get(hit.spacer_id)
        if current is None or hit.mismatches < current[0]:
            best_per_spacer[hit.spacer_id] = (hit.mismatches, hit.target_id, pam)
    contig_stats = contig_level_xgg(best_per_spacer)
    result = {
        "alignments": len(hits),
        "scored_flanks": len(pams),
        "spacer_best_hits_with_pam": len(best_per_spacer),
        "mapped_spacers": len({hit.spacer_id for hit in hits}),
        "xgg_fraction": canonical_xgg(pams),
        "spacer_level_xgg_fraction": canonical_xgg(
            [pam for _m, _t, pam in best_per_spacer.values()]
        ),
        "top_concrete_pams": dict(counts.most_common(5)),
        **contig_stats,
    }
    # Self-check: the reported fraction must equal its own numerator over its own
    # denominator. This is what catches a wrong denominator at run time.
    total = result["distinct_target_contigs"]
    xgg = result["contigs_consensus_xgg"]
    expected = xgg / total if total else 0.0
    if abs(float(result["contig_level_xgg_fraction"]) - expected) > 1e-12:
        raise AssertionError(
            "contig-level fraction is inconsistent with its own numerator and "
            f"denominator: {xgg}/{total} != {result['contig_level_xgg_fraction']}"
        )
    if result["contig_level_xgg_fraction"] > 1.0:
        raise AssertionError("contig-level fraction cannot exceed 1")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--summary",
        default=None,
        help="optional case_summary.json to cross-check the pipeline against",
    )
    args = parser.parse_args()

    case_path = resolve_path(args.case)
    case = json.loads(case_path.read_text(encoding="utf-8"))

    spacer_paths = [resolve_path(str(p)) for p in case["spacers"]]
    target_paths = [resolve_path(str(p)) for p in case.get("targets") or []]
    if not target_paths:
        raise SystemExit(
            "this case has no raw target FASTA, so the homology/orientation "
            "controls cannot be computed. Published aggregate consensus flanks "
            "are one record per spacer, not observed protospacer alignments, and "
            "there is no target library to shuffle. Value assessment for this "
            "evidence class is limited to the aggregate record agreement already "
            "reported in case_summary.json."
        )
    raw_spacers: list[tuple[str, str]] = []
    for path in spacer_paths:
        raw_spacers.extend(read_fasta(path))
    targets: list[tuple[str, str]] = []
    for path in target_paths:
        targets.extend(read_fasta(path))

    orientation = case["spacer_orientation"]["decision"]
    oriented = (
        [(name, reverse_complement(seq)) for name, seq in raw_spacers]
        if orientation == "reverse"
        else list(raw_spacers)
    )
    pam_side = str(case["pam_side"])
    pam_length = int(case["pam_length"])
    max_mismatches = int(case["max_mismatches"])

    rng = random.Random(SEED)
    shuffled_spacers = [
        (name, "".join(rng.sample(list(seq), len(seq)))) for name, seq in oriented
    ]
    shuffled_targets = [
        (name, "".join(rng.sample(list(seq), len(seq)))) for name, seq in targets
    ]
    random_spacers = [
        (f"random_{index}",
         "".join(rng.choice("ACGT") for _ in range(len(seq))))
        for index, (_name, seq) in enumerate(oriented)
    ]
    wrong_orientation = [
        (name, reverse_complement(seq)) for name, seq in oriented
    ]

    controls = {
        "positive_declared_orientation": evaluate(
            oriented, targets, pam_side=pam_side, pam_length=pam_length,
            max_mismatches=max_mismatches,
        ),
        "control_reverse_orientation": evaluate(
            wrong_orientation, targets, pam_side=pam_side, pam_length=pam_length,
            max_mismatches=max_mismatches,
        ),
        "control_shuffled_spacers": evaluate(
            shuffled_spacers, targets, pam_side=pam_side, pam_length=pam_length,
            max_mismatches=max_mismatches,
        ),
        "control_shuffled_targets": evaluate(
            oriented, shuffled_targets, pam_side=pam_side, pam_length=pam_length,
            max_mismatches=max_mismatches,
        ),
        "control_random_spacers": evaluate(
            random_spacers, targets, pam_side=pam_side, pam_length=pam_length,
            max_mismatches=max_mismatches,
        ),
    }

    positive = controls["positive_declared_orientation"]
    control_alignments = [
        value["alignments"] for key, value in controls.items() if key.startswith("control_")
    ]
    report = {
        "format_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "case_id": case["case_id"],
        "case_file": str(case_path.relative_to(WORKSPACE)),
        "case_file_sha256": hashlib.sha256(case_path.read_bytes()).hexdigest(),
        "seed": SEED,
        "parameters": {
            "pam_side": pam_side,
            "pam_length": pam_length,
            "max_mismatches": max_mismatches,
            "declared_orientation": orientation,
            "n_spacers": len(oriented),
            "n_target_contigs": len(targets),
        },
        "target_base_frequencies": background_base_frequencies(targets),
        "controls": controls,
        "verdict": {
            "signal_survives_positive": positive["xgg_fraction"] > 0.5,
            "collapses_under_wrong_orientation": (
                controls["control_reverse_orientation"]["xgg_fraction"] < 0.1
            ),
            "collapses_under_homology_destroying_nulls": max(
                controls["control_shuffled_spacers"]["alignments"],
                controls["control_shuffled_targets"]["alignments"],
                controls["control_random_spacers"]["alignments"],
            ) == 0,
            "contig_redundancy_max_spacers_on_one_contig": positive[
                "max_spacers_on_one_contig"
            ],
            "contig_level_xgg_fraction": positive["contig_level_xgg_fraction"],
            "note": (
                "alignments=0 in the homology-destroying nulls means no spacer "
                "matched at all, so those controls cannot form a PAM signal; "
                "this demonstrates the positive signal requires real "
                "spacer-protospacer homology, not library base composition"
            ),
        },
        "interpretation_limits": [
            "This is a specificity/control check on one engineering-regression case.",
            "It does not establish generalisation: the protein is an exact Protein2PAM training exposure.",
            "The target library is a human gut virome collection, not an S. pyogenes ecological niche.",
            "Target contigs are not clustered to virus level; contig counts overstate independence.",
            "xGG enrichment is evolutionary spacer-protospacer evidence, not cleavage activity.",
        ],
    }

    if args.summary:
        summary_path = resolve_path(args.summary)
        if summary_path.exists():
            pipeline = json.loads(summary_path.read_text(encoding="utf-8"))
            collapse = pipeline.get("spacer_contig_collapse", {})
            report["cross_check_vs_pipeline"] = {
                "pipeline_alignments": pipeline.get("spacer_hits_total"),
                "recomputed_alignments": positive["alignments"],
                "alignments_match": pipeline.get("spacer_hits_total") == positive["alignments"],
                "pipeline_downstream_consensus": (
                    pipeline.get("spacer_side_information", {})
                    .get("downstream", {})
                    .get("consensus")
                ),
                "pipeline_max_spacers_on_one_contig": collapse.get("max_spacers_on_one_contig"),
                "recomputed_max_spacers_on_one_contig": positive["max_spacers_on_one_contig"],
            }

    out_path = resolve_path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
