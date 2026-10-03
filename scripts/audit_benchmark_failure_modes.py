#!/usr/bin/env python3
"""Classify protein/spacer benchmark failures without using post-hoc fixes."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.score.spectrum import allowed_set, normalize_pam, valid_pam  # noqa: E402


def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else WORKSPACE / path


def motifs_overlap(left: str, right: str) -> bool:
    left = normalize_pam(left)
    right = normalize_pam(right)
    return (
        len(left) == len(right)
        and valid_pam(left)
        and valid_pam(right)
        and all(allowed_set(a) & allowed_set(b) for a, b in zip(left, right))
    )


def best_rank(method: dict[str, object]) -> int | None:
    ranks = [
        details.get("best_gold_rank")
        for details in method.get("by_length", {}).values()
        if details.get("best_gold_rank") is not None
    ]
    return min(map(int, ranks)) if ranks else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", default="benchmarks/multi_evidence_v2/manifest.json"
    )
    parser.add_argument(
        "--benchmark-report",
        default="data/parsed/multi_evidence_benchmark_v2/benchmark_report.json",
    )
    parser.add_argument(
        "--replay-audit",
        default="data/parsed/p2pam_training_replay_audit.json",
    )
    parser.add_argument(
        "--out", default="data/parsed/benchmark_failure_modes.json"
    )
    parser.add_argument(
        "--markdown-out", default="docs/BENCHMARK_FAILURE_MODES.md"
    )
    args = parser.parse_args()

    manifest = json.loads(resolve_path(args.manifest).read_text(encoding="utf-8"))
    benchmark = json.loads(
        resolve_path(args.benchmark_report).read_text(encoding="utf-8")
    )
    replay = json.loads(
        resolve_path(args.replay_audit).read_text(encoding="utf-8")
    )
    manifest_by_id = {row["system_id"]: row for row in manifest["systems"]}
    replay_by_id = {row["system_id"]: row for row in replay["systems"]}
    results = []
    counts: dict[str, int] = {}
    for system in benchmark["systems"]:
        system_id = system["system_id"]
        registered = manifest_by_id[system_id]
        replay_row = replay_by_id.get(system_id, {})
        training_labels = list(replay_row.get("distinct_training_consensus_fields", []))
        replay_values = [
            float(row["model_label_flat_cosine"])
            for row in replay_row.get("training_label_replay", [])
            if row.get("model_label_flat_cosine") is not None
        ]
        best_replay_cosine = max(replay_values) if replay_values else None
        golds = [
            motif
            for values in registered["gold_pam_spectrum"].values()
            for motif in values
        ]
        label_gold_overlap = any(
            motifs_overlap(label, gold)
            for label in training_labels
            for gold in golds
        )
        protein_rank = best_rank(system["methods"]["protein_only"])
        spacer_rank = best_rank(system["methods"]["spacer_only"])
        if not training_labels:
            failure_mode = "no_exact_training_label"
            next_action = "retain as exposure control; do not diagnose checkpoint replay"
        elif best_replay_cosine is not None and best_replay_cosine < 0.50:
            failure_mode = "checkpoint_replay_failure"
            next_action = "inspect this exact label, sequence preprocessing, and checkpoint output"
        elif not label_gold_overlap:
            failure_mode = "training_label_vs_experimental_gold_mismatch"
            next_action = "curate assay/strain/PAM-window provenance before any retraining"
        elif best_replay_cosine is not None and best_replay_cosine < 0.75:
            failure_mode = "weak_checkpoint_replay"
            next_action = "inspect duplicate labels and per-position errors; keep fusion abstained"
        elif protein_rank is not None and protein_rank > 5 and spacer_rank is not None and spacer_rank <= 5:
            failure_mode = "protein_ranking_failure_despite_label_alignment"
            next_action = "audit candidate-length scoring and thresholding, not the spacer path"
        else:
            failure_mode = "no_primary_failure_detected"
            next_action = "retain as regression control"
        counts[failure_mode] = counts.get(failure_mode, 0) + 1
        results.append({
            "system_id": system_id,
            "evaluation_role": system.get("evaluation_role"),
            "gold_pams": golds,
            "protein_prediction": replay_row.get("predicted_pam"),
            "exact_training_rows": replay_row.get("exact_training_rows", 0),
            "training_consensus_fields": training_labels,
            "training_label_overlaps_gold": label_gold_overlap,
            "best_model_training_logo_cosine": best_replay_cosine,
            "protein_best_gold_rank": protein_rank,
            "spacer_best_gold_rank": spacer_rank,
            "fusion_gates": system["fusion_diagnostics"]["gate_by_length"],
            "failure_mode": failure_mode,
            "next_action": next_action,
        })
    report = {
        "format_version": 1,
        "classification_order": [
            "no exact training label",
            "checkpoint replay cosine <0.50",
            "training label has no same-length IUPAC overlap with experimental gold",
            "checkpoint replay cosine <0.75",
            "protein rank >5 while spacer rank <=5",
        ],
        "counts": counts,
        "systems": results,
    }
    out_path = resolve_path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# Multi-evidence benchmark failure-mode audit",
        "",
        "This table diagnoses observed failures; it does not retune any threshold on the holdout systems.",
        "",
        "| System | Gold | Protein prediction | Training labels | Replay cosine | Protein rank | Spacer rank | Failure mode |",
        "|---|---|---|---|---:|---:|---:|---|",
    ]
    for row in results:
        lines.append(
            "| {system} | {gold} | {prediction} | {labels} | {cosine} | "
            "{protein} | {spacer} | {mode} |".format(
                system=row["system_id"],
                gold=", ".join(row["gold_pams"]),
                prediction=row["protein_prediction"] or "—",
                labels=", ".join(row["training_consensus_fields"]) or "—",
                cosine=(
                    f"{row['best_model_training_logo_cosine']:.3f}"
                    if row["best_model_training_logo_cosine"] is not None else "—"
                ),
                protein=row["protein_best_gold_rank"] or "—",
                spacer=row["spacer_best_gold_rank"] or "—",
                mode=row["failure_mode"],
            )
        )
    lines.extend([
        "",
        "## Decision",
        "",
        "Do not launch another LoRA run yet. First resolve exact-system label/gold provenance mismatches and the small number of checkpoint replay failures. Keep spacer evidence as an independent ranked result and keep severe cross-source conflicts abstained.",
        "",
    ])
    markdown_path = resolve_path(args.markdown_out)
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_path.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"counts": counts, "systems": len(results)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
