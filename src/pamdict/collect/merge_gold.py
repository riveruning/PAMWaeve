"""Merge resolved Cas sequences onto Cas->PAM pairs to build a clean gold set.

Inputs
------
- ``cas_pam_pairs.tsv``  : Cas name -> PAM (from full-text extraction)
- ``cas_sequences.tsv``   : Cas name -> UniProt sequence (from resolve_uniprot)

Output
------
A ``gold_set.tsv`` where each row is a (protein sequence -> PAM) sample ready to
be fed to protein2PAM, plus full provenance (doi/pmcid/raw_ref/sha256) and a
``sequence_confidence`` column that tells downstream whether the sequence is
human-verified or still needs review. Unresolved pairs are written to a separate
``gold_unresolved.tsv`` for follow-up, never silently dropped.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_sequences(tsv: Path) -> dict[str, dict[str, str]]:
    """Map cas_name -> best sequence record (verified > species > species_review > name_only)."""
    best: dict[str, dict[str, str]] = {}
    rank = {"verified": 0, "species": 1, "species_review": 2, "name_only": 3}
    with tsv.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cas = r.get("cas_name", "").strip()
            if not cas or not r.get("sequence"):
                continue
            m = r.get("match", "none")
            rk = rank.get(m, 99)
            if cas not in best or rk < rank.get(best[cas].get("match", "none"), 99):
                best[cas] = r
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/parsed/cas_pam_text/cas_pam_pairs.tsv")
    ap.add_argument("--sequences", default="data/parsed/cas_pam_text/cas_sequences.tsv")
    ap.add_argument("--out", default="data/parsed/cas_pam_text/gold_set.tsv")
    ap.add_argument("--unresolved", default="data/parsed/cas_pam_text/gold_unresolved.tsv")
    args = ap.parse_args()

    seq_map = load_sequences(Path(args.sequences))

    out_fields = [
        "crispr_type", "cas_family", "source", "protein_id", "citation", "doi",
        "protein_sequence", "pam_consensus", "pam_logo_acgt", "created_at",
        "raw_ref", "sha256", "cas_name", "species", "uni_prot_acc",
        "sequence_confidence",
    ]

    resolved_rows: list[dict[str, str]] = []
    unresolved_rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()

    with Path(args.pairs).open(encoding="utf-8", newline="") as fh:
        for p in csv.DictReader(fh, delimiter="\t"):
            cas = (p.get("cas_name") or "").strip()
            pam = (p.get("pam") or "").strip()
            seq_rec = seq_map.get(cas)
            if seq_rec:
                key = (cas, pam, seq_rec.get("accession", ""))
                if key in seen:
                    continue
                seen.add(key)
                resolved_rows.append({
                    "crispr_type": p.get("crispr_type", "") or _crispr_type(cas),
                    "cas_family": p.get("cas_family", "") or _family(cas),
                    "source": "Literature",
                    "protein_id": seq_rec.get("accession", ""),
                    "citation": p.get("citation", ""),
                    "doi": p.get("doi", ""),
                    "protein_sequence": seq_rec.get("sequence", ""),
                    "pam_consensus": pam,
                    "pam_logo_acgt": "",  # full-text extraction yields consensus only
                    "created_at": p.get("created_at", ""),
                    "raw_ref": p.get("raw_ref", ""),
                    "sha256": p.get("sha256", ""),
                    "cas_name": cas,
                    "species": p.get("species", ""),
                    "uni_prot_acc": seq_rec.get("accession", ""),
                    "sequence_confidence": seq_rec.get("match", "none"),
                })
            else:
                unresolved_rows.append({
                    "cas_name": cas, "pam": pam, "doi": p.get("doi", ""),
                    "pmcid": p.get("pmcid", ""), "species": p.get("species", ""),
                    "raw_ref": p.get("raw_ref", ""),
                })

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        w.writerows(resolved_rows)

    unr_path = Path(args.unresolved)
    with unr_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(unresolved_rows[0].keys()) if unresolved_rows else
                           ["cas_name", "pam", "doi", "pmcid", "species", "raw_ref"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(unresolved_rows)

    print(f"resolved gold rows: {len(resolved_rows)}")
    print(f"unique (cas_name, pam, accession): {len(seen)}")
    print(f"unresolved pairs: {len(unresolved_rows)}")
    print(f"wrote {out_path} and {unr_path}")
    return 0


def _family(cas: str) -> str:
    n = cas.lower()
    if "cas12" in n or "cpf1" in n:
        return "Cas12"
    if "cas13" in n:
        return "Cas13"
    if "cas9" in n:
        return "Cas9"
    return ""


def _crispr_type(cas: str) -> str:
    f = _family(cas)
    return {"Cas12": "Type V", "Cas13": "Type VI", "Cas9": "Type II"}.get(f, "")


if __name__ == "__main__":
    exit(main())
