"""End-to-end supplementary-data extraction for one paper.

Given a Europe PMC PMCID, this:
  1. downloads the supplementaryFiles bundle (zip),
  2. extracts XLSX tables,
  3. parses per-position A/C/G/T frequency matrices ("Sequence logo" sheets),
  4. converts frequencies -> information content (bits),
  5. ingests into the unified schema and writes new_candidates.tsv rows.

Usage:
  PYTHONPATH=src python -m pamdict.collect.run_one \
      --pmcid PMC12627570 --doi 10.1038/s42003-025-08984-y \
      --citation "PAM-readID (2025)" \
      --out-dir data/parsed/demo_pamreadid
"""

from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pamdict.collect.fetch_supplementary import download_file, sha256_file  # noqa: E402
from pamdict.collect.xlsx_reader import read_sheet  # noqa: E402
from pamdict.collect.parse_xlsx_supp import (  # noqa: E402
    find_logo_sheets,
    parse_sequence_logo,
    freq_to_info_content,
)
from pamdict.schema.record import pam_logo_from_matrix  # noqa: E402
from pamdict.schema.constants import PAM_NUCLEOTIDES, PAM_POSITIONS  # noqa: E402

EPMC_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"


def matrix_to_10x4(by_base: dict[str, list[float]]) -> list[list[float]]:
    """Merge {base:[positions]} frequencies into a 10x4 A/C/G/T matrix (freq)."""
    npos = max((len(v) for v in by_base.values()), default=0)
    npos = min(max(npos, 1), len(PAM_POSITIONS))
    matrix = [[0.0] * 4 for _ in range(npos)]
    for j, base in enumerate(PAM_NUCLEOTIDES):
        vals = by_base.get(base, [])
        for i in range(min(len(vals), npos)):
            matrix[i][j] = vals[i]
    return matrix


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmcid", required=True)
    ap.add_argument("--doi", default="")
    ap.add_argument("--citation", default="")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    # 1. download supplementary bundle
    zip_path = out / "supplementary.zip"
    manifest = download_file(f"{EPMC_BASE}/{args.pmcid}/supplementaryFiles", zip_path)
    manifest["sha256_file"] = sha256_file(zip_path)
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    # 2. extract
    extract_dir = out / "extracted"
    extract_dir.mkdir(exist_ok=True)
    samples = []
    with zipfile.ZipFile(zip_path) as z:
        for name in z.namelist():
            if name.lower().endswith(".xlsx"):
                dest = extract_dir / Path(name).name
                z.extract(name, extract_dir)
                rows_by_sheet = {}
                from pamdict.collect.xlsx_reader import list_sheets
                for sh in list_sheets(str(dest)):
                    if "microsoft" in sh.lower():
                        continue
                    rows = read_sheet(str(dest), sh)
                    mat = parse_sequence_logo(rows)
                    if mat.get("A") and mat.get("C") and mat.get("G") and mat.get("T"):
                        freqs = matrix_to_10x4({b: mat[b] for b in PAM_NUCLEOTIDES})
                        bits = freq_to_info_content(freqs)
                        samples.append({
                            "source_sheet": sh,
                            "source_file": Path(name).name,
                            "freq_matrix": freqs,
                            "info_matrix": bits,
                        })

    # 3. write extracted matrices as schema rows
    rows_out = []
    for s in samples:
        rows_out.append({
            "crispr_type": "",
            "cas_family": "",
            "source": "Literature (new; unreviewed)",
            "protein_id": "",
            "citation": args.citation,
            "doi": args.doi,
            "protein_sequence": "",
            "pid_sequence": "",
            "pam_consensus": "",
            "pam_logo_acgt": pam_logo_from_matrix(s["info_matrix"]),
            "raw_ref": f"{args.pmcid}::{s['source_file']}::{s['source_sheet']}",
        })

    summary = {
        "pmcid": args.pmcid,
        "extracted_xlsx": [s["source_file"] for s in samples],
        "logo_sheets": [f"{s['source_file']}::{s['source_sheet']}" for s in samples],
        "rows_written": len(rows_out),
        "units": "information content (bits); source was per-position frequency",
    }
    (out / "extraction_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    if rows_out:
        import csv
        from pamdict.schema.constants import ALL_COLUMNS
        tsv = out / "new_candidates.tsv"
        with tsv.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(ALL_COLUMNS), delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in rows_out:
                w.writerow({c: r.get(c, "") for c in ALL_COLUMNS})
        print(f"wrote {tsv}")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
