#!/usr/bin/env python3
"""Validate a versioned PAMPRIDICT benchmark manifest without loading a model."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.benchmark.contract import (  # noqa: E402
    audit_manifest_document,
    verify_manifest_files,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate benchmark structure, provenance, and independence."
    )
    parser.add_argument(
        "--manifest",
        default="benchmarks/multi_evidence_v2/manifest.json",
    )
    parser.add_argument("--out")
    parser.add_argument(
        "--check-files",
        action="store_true",
        help="verify local input files against provenance SHA-256 values",
    )
    parser.add_argument(
        "--require-strict-paired",
        action="store_true",
        help="fail unless the paired strict cohort reaches the declared minimum",
    )
    args = parser.parse_args()

    path = Path(args.manifest)
    if not path.is_absolute():
        path = WORKSPACE / path
    manifest = json.loads(path.read_text(encoding="utf-8"))
    audit = audit_manifest_document(manifest)
    strict_paired = audit["strict_independent_track_counts"].get(
        "paired_evidence", 0
    )
    minimum = manifest.get("minimum_strict_systems", 5)
    audit["strict_paired_requirement_met"] = strict_paired >= minimum
    if args.check_files:
        audit["file_verification"] = verify_manifest_files(manifest, WORKSPACE)

    rendered = json.dumps(audit, indent=2, ensure_ascii=False)
    print(rendered)
    if args.out:
        output = Path(args.out)
        if not output.is_absolute():
            output = WORKSPACE / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered + "\n", encoding="utf-8")
    if not audit["valid"]:
        return 1
    if args.check_files and not audit["file_verification"]["valid"]:
        return 1
    if args.require_strict_paired and not audit["strict_paired_requirement_met"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
