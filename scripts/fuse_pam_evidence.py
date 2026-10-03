"""Join Protein2PAM candidate scores with spacer-flank evidence."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.score.fusion import fuse_pam_evidence  # noqa: E402


def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else WORKSPACE / path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Late-fuse protein and spacer candidate evidence without inventing "
            "a calibrated joint probability."
        )
    )
    parser.add_argument("--protein-scores", required=True)
    parser.add_argument("--spacer-scores", required=True)
    parser.add_argument("--protein-id")
    parser.add_argument("--system-id")
    parser.add_argument(
        "--policy",
        choices=("abstain", "spacer_fallback", "legacy"),
        default="abstain",
        help=(
            "abstain on severe disagreement (default), run the exploratory "
            "spacer fallback, or reproduce the legacy joint ranking"
        ),
    )
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    protein_rows = read_tsv(resolve_path(args.protein_scores))
    spacer_rows = read_tsv(resolve_path(args.spacer_scores))
    if args.protein_id:
        protein_rows = [
            row for row in protein_rows if row.get("protein_id") == args.protein_id
        ]
    if args.system_id:
        spacer_rows = [
            row for row in spacer_rows if row.get("system_id") == args.system_id
        ]
    if not protein_rows:
        raise ValueError("no protein score rows remain after filtering")
    if not spacer_rows:
        raise ValueError("no spacer score rows remain after filtering")

    rows = fuse_pam_evidence(
        protein_rows,
        spacer_rows,
        policy=args.policy,
    )
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    out_path = resolve_path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
