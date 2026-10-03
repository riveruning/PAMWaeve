"""Parse the official Protein2PAM training TSV into a uniform corpus table.

Outputs a TSV/Parquet in ``data/corpus/`` with the 10 core columns plus
provenance extension columns, and prints per-family statistics (sample counts
and field-completeness).

Usage:
  python src/pamdict/collect/parse_baseline.py \
      --input data/raw/protein2pam_train_seqs.tsv \
      --output data/corpus/official_baseline
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pamdict.schema.constants import (  # noqa: E402
    ALL_COLUMNS,
    CORE_COLUMNS,
    EXTENSION_COLUMNS,
    PAM_NUCLEOTIDES,
    PAM_POSITIONS,
)


def validate_logo(value: str) -> bool:
    """Return True if value parses to a valid 10x4 numeric matrix."""
    try:
        rows = json.loads(value)
    except Exception:
        return False
    if not isinstance(rows, list) or len(rows) != len(PAM_POSITIONS):
        return False
    for row in rows:
        if not isinstance(row, list) or len(row) != len(PAM_NUCLEOTIDES):
            return False
        if not all(isinstance(v, (int, float)) for v in row):
            return False
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)

    # Verify column order matches the contract.
    actual = list(df.columns)
    if actual != list(CORE_COLUMNS):
        raise SystemExit(
            f"Column mismatch.\n expected={list(CORE_COLUMNS)}\n actual  ={actual}"
        )

    n = len(df)
    print(f"total rows: {n}")

    # Per-family sample counts (cross-tab of cas_family by crispr_type).
    print("\n=== per cas_family x crispr_type counts ===")
    print(
        pd.crosstab(df["cas_family"], df["crispr_type"], margins=True)
        .to_string()
    )

    # Field completeness (non-empty fraction).
    print("\n=== field completeness (non-empty %) ===")
    for col in CORE_COLUMNS:
        filled = (df[col].str.len() > 0).sum()
        print(f"  {col:20s} {filled:7d}/{n}  ({100.0*filled/n:.2f}%)")

    # Logo validity.
    valid = df["pam_logo_acgt"].map(validate_logo)
    print(f"\nvalid 10x4 logo JSON: {valid.sum()}/{n}")

    # consensus length distribution (top 10).
    print("\n=== top consensus lengths ===")
    print(df["pam_consensus"].str.len().value_counts().head(10).to_string())

    # Append provenance extension columns (empty for baseline).
    for col in EXTENSION_COLUMNS:
        if col not in df.columns:
            df[col] = ""

    # Reorder to ALL_COLUMNS.
    df = df[list(ALL_COLUMNS)]

    args.output.mkdir(parents=True, exist_ok=True)
    out_tsv = args.output / "official_baseline.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nwrote {out_tsv} ({out_tsv.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
