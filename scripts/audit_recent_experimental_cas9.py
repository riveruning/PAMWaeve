"""Audit recent experimentally characterized Cas9s against Protein2PAM.

The source workbook is PNAS Dataset S1 for Becker et al. (2025). It contains
the full protein sequences and the spacer-derived *predicted* PAMs. The
experimentally measured PAM labels below are transcribed separately from
Figure 4D; keeping those fields distinct prevents predicted PAMs from being
mistaken for benchmark gold labels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.finetune.similarity import ProteinNeighborIndex  # noqa: E402


DOI = "10.1073/pnas.2417674122"
SOURCE_LOCATOR = "Main article Figure 4D and Methods: In Vitro PAM Determination"
ASSAY = "TXTL expression followed by randomized 7N PAM-library cleavage and deep sequencing"
EMPIRICAL_PAMS = {
    "gallolyticus": "NNGYRAH",
    "iniae": "NGG",
    "parasanguinis": "NNAARG",
    "uberis": "NNARTA",
}
CAS_NAMES = {
    "gallolyticus": "SgaCas9",
    "iniae": "SinCas9",
    "parasanguinis": "SpaCas9",
    "uberis": "SubCas9",
}


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference)
    if letters is None:
        raise ValueError(f"invalid XLSX cell reference: {reference}")
    value = 0
    for character in letters.group(0):
        value = value * 26 + ord(character) - ord("A") + 1
    return value - 1


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    return ["".join(item.itertext()) for item in root]


def read_first_xlsx_sheet(path: Path) -> list[list[str]]:
    """Read cell values from the first worksheet using only the stdlib."""
    with zipfile.ZipFile(path) as archive:
        shared = _shared_strings(archive)
        root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows: list[list[str]] = []
    for row_node in root.iter():
        if row_node.tag.split("}")[-1] != "row":
            continue
        values: dict[int, str] = {}
        for cell in row_node:
            if cell.tag.split("}")[-1] != "c":
                continue
            index = _column_index(cell.attrib["r"])
            cell_type = cell.attrib.get("t")
            value_node = next(
                (node for node in cell if node.tag.split("}")[-1] == "v"),
                None,
            )
            if cell_type == "inlineStr":
                value = "".join(cell.itertext())
            elif value_node is None or value_node.text is None:
                value = ""
            elif cell_type == "s":
                value = shared[int(value_node.text)]
            else:
                value = value_node.text
            values[index] = value
        if values:
            width = max(values) + 1
            rows.append([values.get(index, "") for index in range(width)])
    return rows


def is_cas9(row: dict[str, str]) -> bool:
    family = (row.get("cas_family") or "").strip().upper()
    crispr_type = (row.get("crispr_type") or "").replace(" ", "").upper()
    return family.startswith("CAS9") or crispr_type in {"II", "TYPEII"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workbook",
        default=(
            "data/raw/supp_char/extracted/"
            "PMC11929499__pnas.2417674122.sd01.xlsx"
        ),
    )
    parser.add_argument("--official", default="data/raw/protein2pam_train_seqs.tsv")
    parser.add_argument(
        "--out",
        default="benchmarks/multi_evidence_v2/recent_experimental_cas9_audit.tsv",
    )
    parser.add_argument(
        "--summary-out",
        default=(
            "benchmarks/multi_evidence_v2/"
            "recent_experimental_cas9_audit_summary.json"
        ),
    )
    parser.add_argument("--identity", type=float, default=0.90)
    parser.add_argument("--candidate-cosine", type=float, default=0.30)
    parser.add_argument("--max-candidates", type=int, default=512)
    args = parser.parse_args()

    workbook_path = WORKSPACE / args.workbook
    official_path = WORKSPACE / args.official
    workbook_rows = read_first_xlsx_sheet(workbook_path)
    if not workbook_rows:
        raise ValueError(f"workbook contains no rows: {workbook_path}")
    header = workbook_rows[0]
    required = {"Genus", "Species", "Cas9 (amino acids)", "Predicted PAM"}
    if not required <= set(header):
        raise ValueError(f"unexpected workbook columns: {header}")
    workbook_records = [
        dict(zip(header, row + [""] * (len(header) - len(row))))
        for row in workbook_rows[1:]
    ]

    candidates: list[dict[str, str]] = []
    for row in workbook_records:
        species = row["Species"].strip().lower()
        if species not in EMPIRICAL_PAMS:
            continue
        sequence = "".join(row["Cas9 (amino acids)"].split()).upper()
        if not sequence:
            raise ValueError(f"missing Cas9 sequence for {species}")
        candidates.append({
            "candidate_id": CAS_NAMES[species],
            "taxon": f"{row['Genus'].strip()} {row['Species'].strip()}",
            "species_key": species,
            "protein_sequence": sequence,
            "predicted_pam": row["Predicted PAM"].strip().upper(),
            "empirical_pam": EMPIRICAL_PAMS[species],
        })
    if {row["species_key"] for row in candidates} != set(EMPIRICAL_PAMS):
        raise ValueError("not all four Figure 4D Cas9 records were found")

    with official_path.open(encoding="utf-8") as handle:
        official_rows = [
            row for row in csv.DictReader(handle, delimiter="\t") if is_cas9(row)
        ]
    references = [row["protein_sequence"] for row in official_rows]
    metadata_by_sequence: dict[str, list[dict[str, str]]] = {}
    for row in official_rows:
        sequence = row["protein_sequence"].strip().upper()
        metadata_by_sequence.setdefault(sequence, []).append(row)

    index = ProteinNeighborIndex(references, k=3)
    matches = index.find(
        [row["protein_sequence"] for row in candidates],
        identity_threshold=args.identity,
        candidate_cosine=args.candidate_cosine,
        max_candidates=args.max_candidates,
        probe_kmers=None,
    )

    output_rows: list[dict[str, object]] = []
    for candidate in candidates:
        sequence = candidate["protein_sequence"]
        candidate_matches = matches.get(sequence, [])
        exact_rows = metadata_by_sequence.get(sequence, [])
        top_match = candidate_matches[0] if candidate_matches else None
        top_metadata = (
            metadata_by_sequence.get(str(top_match["reference_sequence"]), [])
            if top_match else []
        )
        output_rows.append({
            **candidate,
            "protein_length": len(sequence),
            "protein_sha256": sequence_sha256(sequence),
            "doi": DOI,
            "assay": ASSAY,
            "source_locator": SOURCE_LOCATOR,
            "workbook_path": args.workbook,
            "workbook_sha256": file_sha256(workbook_path),
            "exact_training_match": bool(exact_rows),
            "near90_training_match": bool(candidate_matches),
            "max_confirmed_identity_at_or_above_threshold": (
                top_match["identity"] if top_match else ""
            ),
            "matching_training_pams": ",".join(sorted({
                row.get("pam_consensus", "")
                for row in (exact_rows or top_metadata)
                if row.get("pam_consensus")
            })),
            "independence_status": (
                "EXPOSED_EXACT" if exact_rows else
                "EXPOSED_NEAR90" if candidate_matches else
                "STRICT_CANDIDATE_BELOW90"
            ),
        })

    fieldnames = list(output_rows[0])
    out_path = WORKSPACE / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    status_counts: dict[str, int] = {}
    for row in output_rows:
        status = str(row["independence_status"])
        status_counts[status] = status_counts.get(status, 0) + 1
    summary = {
        "format_version": 1,
        "paper_doi": DOI,
        "evidence_scope": (
            "Four Cas9s with empirical 7N PAM-library profiles in Figure 4D; "
            "workbook predicted PAMs are retained only as non-gold provenance."
        ),
        "identity_definition": "1 - global_levenshtein_distance / max(sequence_lengths)",
        "identity_threshold": args.identity,
        "candidate_cosine": args.candidate_cosine,
        "max_candidates": args.max_candidates,
        "official_cas9_rows": len(official_rows),
        "candidate_count": len(output_rows),
        "status_counts": status_counts,
        "strict_candidate_ids": [
            row["candidate_id"] for row in output_rows
            if row["independence_status"] == "STRICT_CANDIDATE_BELOW90"
        ],
        "output_tsv": args.out,
    }
    summary_path = WORKSPACE / args.summary_out
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"audit -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
