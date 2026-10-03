"""Batch runner for supplementary-data extraction across many papers.

Reads candidate_papers.csv, filters to rows carrying a PMCID, and for each:
  download supplementaryFiles zip -> extract XLSX -> parse Sequence-logo tables
  -> normalize to schema -> append to a merged new_candidates.tsv.

Every row carries provenance (pmcid, doi, source file/sheet, sha256). Papers
that fail (network / no XLSX / no logo table) are logged to failures.tsv, never
silently dropped.

Usage:
  PYTHONPATH=src python -m pamdict.collect.run_batch \
      --candidates data/parsed/discovery/candidate_papers.csv \
      --out-dir data/parsed/batch_v1 \
      --limit 0
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from pamdict.collect.fetch_supplementary import download_file, sha256_file  # noqa: E402
from pamdict.collect.xlsx_reader import read_sheet, list_sheets  # noqa: E402
from pamdict.collect.parse_xlsx_supp import parse_sequence_logo, freq_to_info_content  # noqa: E402
from pamdict.collect.parse_readcount import (  # noqa: E402
    parse_readcount_weighted,
    readcount_freq_to_matrix,
)
from pamdict.collect.parse_pam_table import parse_pam_table  # noqa: E402
from pamdict.collect.parse_position_matrix import parse_position_matrix  # noqa: E402
from pamdict.schema.record import pam_logo_from_matrix  # noqa: E402
from pamdict.schema.constants import PAM_NUCLEOTIDES, PAM_POSITIONS, ALL_COLUMNS  # noqa: E402

EPMC_SUPP = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"


def matrix_to_10x4(mat: dict[str, list[float]]) -> list[list[float]]:
    npos = min(max((len(mat.get(b, [])) for b in PAM_NUCLEOTIDES), default=0), len(PAM_POSITIONS))
    m = [[0.0] * 4 for _ in range(npos)]
    for j, b in enumerate(PAM_NUCLEOTIDES):
        for i, v in enumerate(mat.get(b, [])[:npos]):
            m[i][j] = v
    return m


def _iter_xlsx(zf: zipfile.ZipFile, rel_prefix: str = "", depth: int = 0):
    """Yield (rel_path, bytes) for every .xlsx inside a zip, recursing into
    nested zips (some journals, e.g. NAR, wrap their supplements in an inner
    zip).  Depth-limited to avoid pathological/cyclic archives."""
    if depth > 2:
        return
    for name in zf.namelist():
        low = name.lower()
        if low.endswith(".xlsx"):
            yield f"{rel_prefix}{name}", zf.read(name)
        elif low.endswith(".zip"):
            try:
                inner = zipfile.ZipFile(io.BytesIO(zf.read(name)))
            except zipfile.BadZipFile:
                continue
            yield from _iter_xlsx(inner, f"{rel_prefix}{name}::", depth + 1)


def _is_valid_zip(path: Path) -> bool:
    try:
        with zipfile.ZipFile(path) as z:
            return len(z.namelist()) >= 0  # opened OK
    except (zipfile.BadZipFile, FileNotFoundError):
        return False


def process_one(pmcid: str, doi: str, citation: str, workdir: Path) -> tuple[list[dict], bool]:
    """Return (schema rows, had_xlsx) extracted for one paper (rows may be empty)."""
    workdir.mkdir(parents=True, exist_ok=True)
    zip_path = workdir / "supplementary.zip"
    # Reuse a previously successful download (valid zip) instead of re-fetching,
    # so a re-run only retries the transiently-failed papers.
    if not (zip_path.exists() and _is_valid_zip(zip_path)):
        download_file(EPMC_SUPP.format(pmcid=pmcid), zip_path)
    zip_sha = sha256_file(zip_path)

    rows_out: list[dict] = []
    had_xlsx = False
    with zipfile.ZipFile(zip_path) as z:
        for idx, (rel_path, xlsx_bytes) in enumerate(_iter_xlsx(z)):
            had_xlsx = True
            # Write the XLSX bytes to a temp path (basename disambiguated by
            # nesting order) so the path-based sheet reader can open it.
            stem = Path(rel_path.split("::")[-1]).stem or "tables"
            dest = workdir / f"{stem}_{idx}.xlsx"
            dest.write_bytes(xlsx_bytes)
            for sheet in list_sheets(str(dest)):
                if "microsoft" in sheet.lower():
                    continue
                sheet_rows = read_sheet(str(dest), sheet)

                # Path A: direct "PAM table" (explicit validated PAM strings,
                # possibly with protein_sequence / organism / target).
                for prec in parse_pam_table(sheet_rows):
                    prec["citation"] = citation
                    prec["doi"] = doi
                    prec["sha256"] = zip_sha
                    prec["raw_ref"] = f"{pmcid}::{rel_path}::{sheet}::pam_table::" + prec.get("raw_ref", "")
                    rows_out.append(prec)

                # Path B: per-position frequency matrices. Try pre-aggregated
                # "Sequence logo", then raw read-count/enrichment, then a
                # per-position count matrix; any that yields A/C/G/T positions
                # becomes a bits matrix.
                mat = parse_sequence_logo(sheet_rows)
                logo_kind = None
                freqs = None
                if mat.get("A") and mat.get("C") and mat.get("G") and mat.get("T"):
                    freqs = matrix_to_10x4({b: mat[b] for b in PAM_NUCLEOTIDES})
                    logo_kind = "sequence_logo"
                else:
                    rc, weight_mode = parse_readcount_weighted(sheet_rows)
                    if rc and all(b in rc for b in PAM_NUCLEOTIDES):
                        freqs = readcount_freq_to_matrix(rc)
                        logo_kind = f"read_count_{weight_mode}"
                    else:
                        pm = parse_position_matrix(sheet_rows)
                        if pm and all(b in pm for b in PAM_NUCLEOTIDES):
                            freqs = readcount_freq_to_matrix(pm)
                            logo_kind = "position_matrix"

                if freqs is None:
                    continue

                bits = freq_to_info_content(freqs)
                row = {c: "" for c in ALL_COLUMNS}
                row["source"] = "Literature (new; unreviewed)"
                row["citation"] = citation
                row["doi"] = doi
                row["pam_logo_acgt"] = pam_logo_from_matrix(bits)
                row["raw_ref"] = f"{pmcid}::{rel_path}::{sheet}::{logo_kind}"
                row["sha256"] = zip_sha
                rows_out.append(row)
    return rows_out, had_xlsx


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--limit", type=int, default=0, help="0 = all")
    args = ap.parse_args()

    with args.candidates.open(newline="", encoding="utf-8") as fh:
        cands = [r for r in csv.DictReader(fh) if r["pmcid"].strip()]

    if args.limit:
        cands = cands[: args.limit]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    merged = args.out_dir / "new_candidates.tsv"
    failures = args.out_dir / "failures.tsv"
    all_rows: list[dict] = []
    failure_rows: list[dict] = []

    for i, c in enumerate(cands, 1):
        pmcid = c["pmcid"]
        doi = c["doi"]
        citation = f"{c['title'][:80]} ({c['year']})"
        wdir = args.out_dir / "per_paper" / pmcid
        try:
            rows, had_xlsx = process_one(pmcid, doi, citation, wdir)
            all_rows.extend(rows)
            if not rows:
                reason = "no_xlsx_in_supp" if not had_xlsx else "xlsx_present_but_no_pam_table"
                failure_rows.append({"pmcid": pmcid, "doi": doi, "error": reason})
                print(f"[{i}/{len(cands)}] {pmcid}: 0 rows ({reason})", flush=True)
            else:
                print(f"[{i}/{len(cands)}] {pmcid}: {len(rows)} rows", flush=True)
        except Exception as e:  # noqa: BLE001
            failure_rows.append({"pmcid": pmcid, "doi": doi, "error": str(e)[:200]})
            print(f"[{i}/{len(cands)}] {pmcid}: FAIL {type(e).__name__}: {str(e)[:120]}", flush=True)

    # write merged
    with merged.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ALL_COLUMNS), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in all_rows:
            w.writerow(r)

    with failures.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["pmcid", "doi", "error"], delimiter="\t")
        w.writeheader()
        for r in failure_rows:
            w.writerow(r)

    summary = {
        "papers_attempted": len(cands),
        "papers_failed": len(failure_rows),
        "total_rows": len(all_rows),
        "merged_file": str(merged),
        "failures_file": str(failures),
    }
    (args.out_dir / "batch_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
