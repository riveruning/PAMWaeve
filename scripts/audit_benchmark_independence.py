"""Audit whether a benchmark can support protein-to-PAM generalization claims.

The audit keeps the source TSV files unchanged.  It joins the current
fine-tuning train set, the immutable Protein2PAM training export and the
verified >=90% protein-neighbour artifact at the benchmark-system level.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.score.spectrum import sequence_sha256, system_id  # noqa: E402


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def norm_sequence(row: dict[str, str]) -> str:
    return (row.get("protein_sequence") or "").strip().upper()


def norm_doi(value: str | None) -> str:
    result = (value or "").strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if result.startswith(prefix):
            result = result[len(prefix):]
    return result.rstrip(".")


def is_cas9(row: dict[str, str]) -> bool:
    raw_type = (row.get("crispr_type") or "").replace(" ", "").upper()
    family = (row.get("cas_family") or "").strip().upper()
    return raw_type in {"II", "TYPEII"} or (not raw_type and family.startswith("CAS9"))


def unique_nonempty(values: list[str]) -> list[str]:
    return sorted({value for value in values if value})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--train", default="data/corpus/augmented_train_v2_cas9_train.tsv"
    )
    parser.add_argument(
        "--benchmark", default="data/corpus/gold_benchmark_v2.tsv"
    )
    parser.add_argument(
        "--official", default="data/raw/protein2pam_train_seqs.tsv"
    )
    parser.add_argument(
        "--neighbors",
        default="data/parsed/cas9_official_near_neighbors_v2.json",
    )
    parser.add_argument(
        "--details-out",
        default="data/parsed/benchmark_independence_audit_v1.tsv",
    )
    parser.add_argument(
        "--out",
        default="data/parsed/benchmark_independence_audit_v1.json",
    )
    args = parser.parse_args()

    train_path = WORKSPACE / args.train
    benchmark_path = WORKSPACE / args.benchmark
    official_path = WORKSPACE / args.official
    neighbor_path = WORKSPACE / args.neighbors

    train = read_tsv(train_path)
    benchmark = read_tsv(benchmark_path)
    official = read_tsv(official_path)
    neighbor_payload = json.loads(neighbor_path.read_text(encoding="utf-8"))

    expected_benchmark_sha = neighbor_payload.get("benchmark_sha256")
    expected_official_sha = neighbor_payload.get("official_sha256")
    if expected_benchmark_sha != file_sha256(benchmark_path):
        raise ValueError("near-neighbour artifact does not match the benchmark")
    if expected_official_sha != file_sha256(official_path):
        raise ValueError("near-neighbour artifact does not match the official set")
    threshold = float(neighbor_payload.get("identity_threshold", 0))
    if threshold != 0.90:
        raise ValueError("the locked near-neighbour threshold must be 0.90")

    train_sequences = {norm_sequence(row) for row in train if norm_sequence(row)}
    official_sequences = {
        norm_sequence(row) for row in official if norm_sequence(row) and is_cas9(row)
    }
    train_dois = {
        norm_doi(row.get("doi")) for row in train if norm_doi(row.get("doi"))
    }
    matches_by_query = neighbor_payload.get("benchmark_matches") or {}

    rows_by_sequence: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in benchmark:
        sequence = norm_sequence(row)
        if not sequence:
            raise ValueError("benchmark contains an empty protein sequence")
        rows_by_sequence[sequence].append(row)

    detail_rows: list[dict[str, object]] = []
    for sequence, rows in sorted(rows_by_sequence.items()):
        cas9 = any(is_cas9(row) for row in rows)
        matches = matches_by_query.get(sequence, []) if cas9 else []
        max_identity = max(
            (float(match.get("identity", 0)) for match in matches),
            default=0.0,
        )
        if not cas9:
            base_exposure = "not_applicable_non_cas9"
        elif sequence in official_sequences:
            base_exposure = "exact"
        elif max_identity >= threshold:
            base_exposure = "near_90"
        else:
            base_exposure = "below_90"

        dois = unique_nonempty([norm_doi(row.get("doi")) for row in rows])
        overlapping_dois = sorted(set(dois) & train_dois)
        logo_rows = sum(bool(
            (row.get("pam_logo_acgt") or row.get("logo_json") or "").strip()
        ) for row in rows)
        strictly_independent = (
            cas9
            and base_exposure == "below_90"
            and sequence not in train_sequences
            and not overlapping_dois
        )
        detail_rows.append({
            "system_id": system_id(
                sequence, (rows[0].get("cas_family") or "").strip()
            ),
            "sequence_sha256": sequence_sha256(sequence),
            "family_group": "Cas9" if cas9 else "non-Cas9",
            "protein_ids_json": json.dumps(unique_nonempty([
                (row.get("protein_id") or "").strip() for row in rows
            ]), ensure_ascii=False),
            "dois_json": json.dumps(dois, ensure_ascii=False),
            "bench_layers_json": json.dumps(unique_nonempty([
                (row.get("bench_layer") or "").strip() for row in rows
            ]), ensure_ascii=False),
            "leakage_labels_json": json.dumps(unique_nonempty([
                (row.get("leakage") or "").strip() for row in rows
            ]), ensure_ascii=False),
            "observation_rows": len(rows),
            "unique_pams": len({
                (row.get("pam_consensus") or "").strip().upper() for row in rows
            }),
            "logo_observation_rows": logo_rows,
            "base_pretraining_exposure": base_exposure,
            "max_official_identity": round(max_identity, 6) if cas9 else "",
            "finetune_exact_sequence_overlap": int(sequence in train_sequences),
            "finetune_doi_overlap_json": json.dumps(
                overlapping_dois, ensure_ascii=False
            ),
            "strictly_independent": int(strictly_independent),
            "strictly_independent_with_logo": int(
                strictly_independent and logo_rows > 0
            ),
        })

    cas9_details = [row for row in detail_rows if row["family_group"] == "Cas9"]
    logo_details = [
        row for row in cas9_details if int(row["logo_observation_rows"]) > 0
    ]
    exposure_counts = Counter(
        str(row["base_pretraining_exposure"]) for row in cas9_details
    )
    finetune_doi_overlap = [
        row for row in cas9_details
        if json.loads(str(row["finetune_doi_overlap_json"]))
    ]
    summary = {
        "format_version": 1,
        "purpose": "benchmark_independence_audit",
        "inputs": {
            "train": args.train,
            "train_sha256": file_sha256(train_path),
            "benchmark": args.benchmark,
            "benchmark_sha256": file_sha256(benchmark_path),
            "official": args.official,
            "official_sha256": file_sha256(official_path),
            "neighbors": args.neighbors,
            "neighbors_sha256": file_sha256(neighbor_path),
        },
        "definitions": {
            "base_exact": "benchmark protein is an exact Protein2PAM training sequence",
            "base_near_90": "benchmark protein has a >=90% global-identity Protein2PAM training neighbour",
            "strictly_independent": (
                "Cas9, below 90% identity to official Protein2PAM training, "
                "no exact sequence in fine-tuning train, and no DOI shared "
                "with fine-tuning train"
            ),
        },
        "benchmark": {
            "rows": len(benchmark),
            "systems": len(detail_rows),
            "cas9_systems": len(cas9_details),
            "cas9_logo_systems": len(logo_details),
        },
        "base_pretraining_exposure": dict(exposure_counts),
        "fine_tuning_exposure": {
            "exact_sequence_systems": sum(
                int(row["finetune_exact_sequence_overlap"]) for row in cas9_details
            ),
            "shared_doi_systems": len(finetune_doi_overlap),
            "shared_dois": sorted({
                doi
                for row in finetune_doi_overlap
                for doi in json.loads(str(row["finetune_doi_overlap_json"]))
            }),
        },
        "strict_generalization_eligibility": {
            "cas9_systems": sum(
                int(row["strictly_independent"]) for row in cas9_details
            ),
            "cas9_logo_systems": sum(
                int(row["strictly_independent_with_logo"]) for row in cas9_details
            ),
        },
        "interpretation": (
            "Use this benchmark for regression and calibration checks.  Do not "
            "claim protein-to-PAM generalization from it unless the relevant "
            "analysis is restricted to strictly independent systems."
        ),
    }

    detail_path = WORKSPACE / args.details_out
    out_path = WORKSPACE / args.out
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(detail_rows[0])
    with detail_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="	")
        writer.writeheader()
        writer.writerows(detail_rows)
    out_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"details -> {detail_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
