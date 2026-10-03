#!/usr/bin/env python3
"""Build compact same-host Cas9/spacer/published-flank benchmark artifacts."""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import io
import json
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.score.spacer import normalize_dna, reverse_complement  # noqa: E402


def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else WORKSPACE / path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta_record(path: Path, record_id: str) -> tuple[str, str]:
    header = ""
    sequence: list[str] = []
    found = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if found:
                    break
                found = line[1:].split()[0] == record_id
                if found:
                    header = line[1:]
            elif found:
                sequence.append(line)
    result = "".join(sequence).upper()
    if not result:
        raise ValueError(f"FASTA record {record_id!r} not found in {path}")
    return header, result


def read_tsv_protein(
    path: Path,
    record_id: str,
    *,
    source_name: str,
) -> tuple[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        matches = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row.get("protein_id") == record_id
            and row.get("source") == source_name
            and row.get("protein_sequence", "").strip()
        ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one {source_name} protein {record_id!r} in {path}, "
            f"found {len(matches)}"
        )
    sequence = "".join(matches[0]["protein_sequence"].split()).upper()
    return f"{record_id} source={source_name}", sequence


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def accession_tokens(value: str) -> set[str]:
    return {part.strip() for part in value.split("|") if part.strip()}


def normalized_orientation(value: str) -> str | None:
    values = {
        part.strip().lower()
        for part in value.split("|")
        if part.strip().lower() in {"forward", "reverse"}
    }
    return next(iter(values)) if len(values) == 1 else None


def choose_orientation(row: dict[str, str]) -> tuple[str, str]:
    pam = normalized_orientation(row["orientation_PAMbased"])
    if pam:
        return pam, "orientation_PAMbased"
    crispr = normalized_orientation(row["orientation_CRISPRCasdb"])
    if crispr:
        return crispr, "orientation_CRISPRCasdb_fallback"
    raise ValueError(
        f"source row {row['Unnamed: 0']} has no unambiguous orientation"
    )


def orient_row(row: dict[str, str], boundary_overlap_nt: int) -> dict[str, object]:
    orientation, orientation_source = choose_orientation(row)
    source_spacer = normalize_dna(row["spacers"])
    spacer = (
        source_spacer if orientation == "forward"
        else reverse_complement(source_spacer)
    )
    raw_left = ""
    raw_right = ""
    if row["hit"] == "1":
        parsed = ast.literal_eval(row["consensus_flanks"])
        if not isinstance(parsed, tuple) or len(parsed) != 2:
            raise ValueError(
                f"source row {row['Unnamed: 0']} has invalid consensus_flanks"
            )
        raw_left, raw_right = map(normalize_dna, parsed)
    elif row["consensus_flanks"].strip():
        raise ValueError(
            f"source row {row['Unnamed: 0']} is hit=0 but has flanks"
        )
    if orientation == "forward":
        upstream_with_overlap, downstream_with_overlap = raw_left, raw_right
    else:
        upstream_with_overlap = reverse_complement(raw_right) if raw_right else ""
        downstream_with_overlap = reverse_complement(raw_left) if raw_left else ""
    if upstream_with_overlap and len(upstream_with_overlap) <= boundary_overlap_nt:
        raise ValueError("upstream flank is no longer than its boundary overlap")
    if downstream_with_overlap and len(downstream_with_overlap) <= boundary_overlap_nt:
        raise ValueError("downstream flank is no longer than its boundary overlap")
    upstream = (
        upstream_with_overlap[:-boundary_overlap_nt]
        if upstream_with_overlap else ""
    )
    downstream = (
        downstream_with_overlap[boundary_overlap_nt:]
        if downstream_with_overlap else ""
    )
    row_id = str(row["Unnamed: 0"])
    return {
        "spacer_id": f"source_{row_id}",
        "source_row_id": row_id,
        "source_spacer": source_spacer,
        "oriented_spacer": spacer,
        "repeat_sequence": normalize_dna(row["repeats"]),
        "matched_accessions": row["accessionnrs"],
        "subtype": row["subtype"],
        "cas_genes": row["cas_genes"],
        "hit": row["hit"],
        "raw_left_flank": raw_left,
        "raw_right_flank": raw_right,
        "oriented_upstream_flank": upstream,
        "oriented_downstream_flank": downstream,
        "orientation": orientation,
        "orientation_source": orientation_source,
        "protospacer_boundary_overlap_nt": boundary_overlap_nt,
        "published_pam_field": row["PAM"],
        "repeat_cluster": row["repeat_cluster"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rebuild the retrospective published-flank cohort."
    )
    parser.add_argument(
        "--registry",
        default="benchmarks/multi_evidence_v2/published_flank_sources.json",
    )
    parser.add_argument(
        "--supplement-zip",
        default=(
            "data/raw/external_datasets/"
            "CRISPRCasDB_PAM_repeat_PMC8482600_supplementary.zip"
        ),
    )
    parser.add_argument(
        "--ncbi-root",
        default=(
            "data/raw/external_datasets/ncbi_cas9_candidate_genomes/"
            "ncbi_dataset/data"
        ),
    )
    parser.add_argument(
        "--out-root",
        default="benchmarks/multi_evidence_v2/systems",
    )
    args = parser.parse_args()

    registry_path = resolve_path(args.registry)
    supplement_path = resolve_path(args.supplement_zip)
    ncbi_root = resolve_path(args.ncbi_root)
    out_root = resolve_path(args.out_root)
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    systems = registry["systems"]
    boundary_overlap_nt = int(registry["protospacer_boundary_overlap_nt"])
    if boundary_overlap_nt <= 0:
        raise ValueError("protospacer_boundary_overlap_nt must be positive")
    by_accession = {
        str(system["chromosome_accession"]): system for system in systems
    }
    rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    with zipfile.ZipFile(supplement_path) as archive:
        names = [
            name for name in archive.namelist()
            if name.endswith("13059_2021_2495_MOESM2_ESM.csv")
        ]
        if len(names) != 1:
            raise ValueError(f"expected one Additional file 2 CSV, got {names}")
        with archive.open(names[0]) as binary:
            text = io.TextIOWrapper(binary, encoding="utf-8-sig", newline="")
            for row in csv.DictReader(text):
                tokens = accession_tokens(row["accessionnrs"])
                for accession in tokens & by_accession.keys():
                    rows[accession].append(row)

    results = []
    for system in systems:
        system_id = str(system["system_id"])
        accession = str(system["chromosome_accession"])
        source_rows = sorted(
            rows.get(accession, []), key=lambda item: int(item["Unnamed: 0"])
        )
        source_subtype = str(system.get("source_subtype", "")).strip()
        if source_subtype:
            source_rows = [
                row for row in source_rows
                if row.get("subtype") == source_subtype
            ]
        if len(source_rows) != int(system["expected_spacer_rows"]):
            raise ValueError(
                f"{system_id}: expected {system['expected_spacer_rows']} rows, "
                f"found {len(source_rows)}"
            )
        hit_count = sum(row["hit"] == "1" for row in source_rows)
        if hit_count != int(system["expected_hit_rows"]):
            raise ValueError(
                f"{system_id}: expected {system['expected_hit_rows']} hits, "
                f"found {hit_count}"
            )
        oriented = [
            orient_row(row, boundary_overlap_nt) for row in source_rows
        ]
        sequences = [str(row["oriented_spacer"]) for row in oriented]
        if len(sequences) != len(set(sequences)):
            raise ValueError(f"{system_id}: oriented spacers are not unique")

        protein_id = str(system["protein_record_id"])
        if system.get("protein_sequence_source") == "gasiunas_gold_master":
            protein_source = resolve_path(str(system["protein_source_path"]))
            source_header, protein_sequence = read_tsv_protein(
                protein_source,
                protein_id,
                source_name="Gasiunas",
            )
            protein_source_kind = "Gasiunas exact experimental protein"
        else:
            assembly = str(system["assembly_accession"])
            protein_source = ncbi_root / assembly / "protein.faa"
            source_header, protein_sequence = read_fasta_record(
                protein_source, protein_id
            )
            protein_source_kind = "NCBI assembly protein"
        system_dir = out_root / system_id
        system_dir.mkdir(parents=True, exist_ok=True)
        protein_path = system_dir / "protein.faa"
        spacer_path = system_dir / "spacers.fna"
        flanks_path = system_dir / "published_flanks.tsv"
        summary_path = system_dir / "source_summary.json"
        write_fasta(protein_path, [(source_header, protein_sequence)])
        write_fasta(
            spacer_path,
            [
                (str(row["spacer_id"]), str(row["oriented_spacer"]))
                for row in oriented
            ],
        )
        write_tsv(flanks_path, oriented)
        summary = {
            "format_version": 1,
            "system": system,
            "source": {
                "doi": registry["source_doi"],
                "table": registry["source_table"],
                "supplement_zip_sha256": sha256_file(supplement_path),
                "protein_source": str(protein_source.relative_to(WORKSPACE)),
                "protein_source_kind": protein_source_kind,
                "source_subtype_filter": source_subtype or None,
                "flank_normalization": {
                    "published_sequence_length_nt": 26,
                    "external_flank_length_nt": 23,
                    "protospacer_boundary_overlap_removed_nt": boundary_overlap_nt,
                },
            },
            "counts": {
                "spacer_rows": len(oriented),
                "published_flank_rows": hit_count,
                "no_hit_rows": len(oriented) - hit_count,
            },
            "orientation_sources": dict(sorted({
                source: sum(row["orientation_source"] == source for row in oriented)
                for source in {str(row["orientation_source"]) for row in oriented}
            }.items())),
            "artifacts": {
                "protein_faa_sha256": sha256_file(protein_path),
                "protein_sequence_sha256": hashlib.sha256(
                    protein_sequence.encode("utf-8")
                ).hexdigest(),
                "spacers_fna_sha256": sha256_file(spacer_path),
                "published_flanks_tsv_sha256": sha256_file(flanks_path),
            },
            "limitations": [
                "Published flanks are per-spacer aggregate consensus sequences.",
                "Raw target contig identifiers and per-hit alignments are unavailable.",
                "These records are retrospective evolutionary evidence, not cleavage assays.",
            ],
        }
        summary_path.write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        results.append(summary)
        print(
            f"{system_id}: spacers={len(oriented)} published_flanks={hit_count}"
        )
    print(json.dumps({
        "systems_built": len(results),
        "out_root": str(out_root),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
