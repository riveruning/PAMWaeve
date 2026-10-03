#!/usr/bin/env python3
"""Compare cached Protein2PAM outputs with exact matching training labels."""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.infer.p2pam import consensus_from_info, prob_to_info_numpy  # noqa: E402


def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else WORKSPACE / path


def sequence_sha256(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper()
    return hashlib.sha256(normalized.encode("ascii")).hexdigest()


def cosine(left: np.ndarray, right: np.ndarray) -> float | None:
    a = np.asarray(left, dtype=float).reshape(-1)
    b = np.asarray(right, dtype=float).reshape(-1)
    if a.shape != b.shape:
        return None
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denominator) if denominator else None


def mean(values: list[float]) -> float | None:
    return float(statistics.fmean(values)) if values else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--training", default="data/raw/protein2pam_train_seqs.tsv"
    )
    parser.add_argument(
        "--predictions-root",
        default="data/parsed/multi_evidence_benchmark_v2",
    )
    parser.add_argument(
        "--out", default="data/parsed/p2pam_training_replay_audit.json"
    )
    parser.add_argument(
        "--markdown-out", default="docs/P2PAM_TRAINING_REPLAY_AUDIT.md"
    )
    args = parser.parse_args()

    labels: dict[str, list[dict[str, str]]] = defaultdict(list)
    with resolve_path(args.training).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequence = row.get("protein_sequence", "")
            if sequence:
                labels[sequence_sha256(sequence)].append(row)

    details: list[dict[str, object]] = []
    prediction_paths = sorted(
        resolve_path(args.predictions_root).glob("*/protein_scores.json")
    )
    for prediction_path in prediction_paths:
        payload = json.loads(prediction_path.read_text(encoding="utf-8"))
        for protein in payload.get("proteins", []):
            model_info = prob_to_info_numpy(
                np.asarray(protein["probability_matrix"], dtype=np.float32)
            )
            sequence_hash = str(protein["sequence_sha256"])
            exact_rows = labels.get(sequence_hash, [])
            label_results = []
            for row in exact_rows:
                label_info = np.asarray(
                    ast.literal_eval(row["pam_logo_acgt"]), dtype=np.float32
                )
                matrix_cosine = cosine(model_info, label_info)
                per_position = [
                    value
                    for value in (
                        cosine(model_row, label_row)
                        for model_row, label_row in zip(model_info, label_info)
                    )
                    if value is not None
                ]
                label_results.append({
                    "source": row.get("source", ""),
                    "protein_id": row.get("protein_id", ""),
                    "doi": row.get("doi", ""),
                    "pam_consensus_field": row.get("pam_consensus", ""),
                    "pam_consensus_from_logo": consensus_from_info(
                        label_info, side=str(payload.get("pam_side", "downstream"))
                    ),
                    "model_label_flat_cosine": matrix_cosine,
                    "model_label_mean_position_cosine": mean(per_position),
                    "mean_absolute_information_error": (
                        float(np.mean(np.abs(model_info - label_info)))
                        if model_info.shape == label_info.shape else None
                    ),
                    "model_consensus_matches_label_field": (
                        str(protein.get("predicted_pam", ""))
                        == str(row.get("pam_consensus", ""))
                    ),
                })
            details.append({
                "system_id": (
                    protein.get("protein_id") or prediction_path.parent.name
                ),
                "prediction_path": str(prediction_path.relative_to(WORKSPACE)),
                "model": payload.get("model"),
                "sequence_sha256": sequence_hash,
                "predicted_pam": protein.get("predicted_pam"),
                "exact_training_rows": len(exact_rows),
                "distinct_training_consensus_fields": sorted({
                    row.get("pam_consensus", "") for row in exact_rows
                }),
                "training_label_replay": label_results,
            })

    cosine_values = [
        float(label["model_label_flat_cosine"])
        for item in details
        for label in item["training_label_replay"]
        if label["model_label_flat_cosine"] is not None
    ]
    consensus_checks = [
        bool(label["model_consensus_matches_label_field"])
        for item in details
        for label in item["training_label_replay"]
    ]
    labels_by_source: dict[str, list[dict[str, object]]] = defaultdict(list)
    for item in details:
        for label in item["training_label_replay"]:
            labels_by_source[str(label["source"] or "unspecified")].append(label)
    source_summary = {}
    for source, source_labels in sorted(labels_by_source.items()):
        source_cosines = [
            float(label["model_label_flat_cosine"])
            for label in source_labels
            if label["model_label_flat_cosine"] is not None
        ]
        source_summary[source] = {
            "labels": len(source_labels),
            "mean_model_label_flat_cosine": mean(source_cosines),
            "median_model_label_flat_cosine": (
                float(statistics.median(source_cosines))
                if source_cosines else None
            ),
            "exact_consensus_replay_rate": (
                sum(bool(label["model_consensus_matches_label_field"])
                    for label in source_labels) / len(source_labels)
                if source_labels else None
            ),
            "logo_cosine_below_0_50": sum(
                value < 0.50 for value in source_cosines
            ),
            "logo_cosine_below_0_75": sum(
                value < 0.75 for value in source_cosines
            ),
        }
    report = {
        "format_version": 1,
        "scope": "cached benchmark proteins with exact full-sequence training matches",
        "interpretation": (
            "A mismatch is a closed-loop label-replay diagnostic; by itself it "
            "does not distinguish underfitting, conflicting labels, or data bugs."
        ),
        "prediction_files_discovered": len(prediction_paths),
        "proteins_audited": len(details),
        "proteins_with_exact_training_match": sum(
            int(item["exact_training_rows"]) > 0 for item in details
        ),
        "exact_training_labels_compared": len(cosine_values),
        "mean_model_label_flat_cosine": mean(cosine_values),
        "median_model_label_flat_cosine": (
            float(statistics.median(cosine_values)) if cosine_values else None
        ),
        "exact_consensus_replay_rate": (
            sum(consensus_checks) / len(consensus_checks)
            if consensus_checks else None
        ),
        "logo_cosine_below_0_50": sum(value < 0.50 for value in cosine_values),
        "logo_cosine_below_0_75": sum(value < 0.75 for value in cosine_values),
        "source_summary": source_summary,
        "systems": details,
    }
    out_path = resolve_path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# Protein2PAM exact-training replay audit",
        "",
        "This is a closed-loop diagnostic, not a generalization benchmark.",
        "",
        f"- Cached proteins audited: {report['proteins_audited']}",
        f"- Exact training matches: {report['proteins_with_exact_training_match']}",
        f"- Training labels compared: {report['exact_training_labels_compared']}",
        f"- Mean model/label logo cosine: {report['mean_model_label_flat_cosine']}",
        f"- Exact consensus replay rate: {report['exact_consensus_replay_rate']}",
        f"- Logo cosine below 0.50: {report['logo_cosine_below_0_50']}",
        f"- Logo cosine below 0.75: {report['logo_cosine_below_0_75']}",
        "",
        "## Source split",
        "",
        "| Source | Labels | Mean cosine | Median cosine | Exact consensus | Cosine <0.50 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for source, values in source_summary.items():
        lines.append(
            "| {source} | {labels} | {mean:.3f} | {median:.3f} | "
            "{exact:.3f} | {low} |".format(
                source=source,
                labels=values["labels"],
                mean=values["mean_model_label_flat_cosine"],
                median=values["median_model_label_flat_cosine"],
                exact=values["exact_consensus_replay_rate"],
                low=values["logo_cosine_below_0_50"],
            )
        )
    lines.extend([
        "",
        "| System | Prediction | Exact rows | Training consensus fields | Best logo cosine |",
        "|---|---|---:|---|---:|",
    ])
    for item in details:
        row_cosines = [
            float(value["model_label_flat_cosine"])
            for value in item["training_label_replay"]
            if value["model_label_flat_cosine"] is not None
        ]
        lines.append(
            "| {system} | {prediction} | {rows} | {labels} | {cosine} |".format(
                system=item["system_id"],
                prediction=item["predicted_pam"],
                rows=item["exact_training_rows"],
                labels=", ".join(item["distinct_training_consensus_fields"]) or "—",
                cosine=(f"{max(row_cosines):.3f}" if row_cosines else "—"),
            )
        )
    lines.extend([
        "",
        "A poor exact-training replay result should pause further fine-tuning until label orientation, logo conversion, duplicate-label conflicts, and checkpoint compatibility are audited.",
        "",
    ])
    markdown_path = resolve_path(args.markdown_out)
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_path.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({
        key: value for key, value in report.items() if key != "systems"
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
