"""Extract "Cas ortholog -> PAM" records from characterization-paper XLSX.

Large-scale Cas characterization papers (e.g. Gasiunas 2020, Karvelis 2015,
Cas12a surveys) publish a master supplementary table where each row is one
ortholog/effector with columns for: organism/strain, PAM consensus, direct
repeat, tracrRNA/sgRNA, and the full protein (AA) sequence.

This module auto-detects those columns and emits one unified record per row
that has at least (PAM) + (protein sequence OR organism). The result is a
structured TSV in the style the training pipeline needs (protein + PAM + repeat
+ organism + provenance), every row traceable to its source paper.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import zipfile
from pathlib import Path

from .xlsx_reader import read_sheet, list_sheets

# Column-name keywords, in priority order for each semantic role.
_PAM_KEYS = ["pam", "protospacer adjacent", "consensus pam", "pam consensus"]
# Blocklist: column names that contain a PAM keyword but hold a COMPUTATIONAL
# prediction, not an experimentally-determined PAM. The user's core rule is
# "only experimentally-measured PAM, never model-predicted", so these must not
# feed the gold set.
_PAM_BLOCKLIST = [
    "predicted pam", "pam prediction", "pred pam", "pam pred",
    "iupac pam seq", "inferred pam", "machine learning",
]
_ORGANISM_KEYS = ["organism", "strain", "species", "host", "source organism",
                  "source_strain", "bacter", "host strain"]
_AA_KEYS = ["aa sequence", "amino acid", "protein sequence", "aa_seq",
            "cas9 aa", "cas12 aa", "effector protein", "amino acid seq"]
_REPEAT_KEYS = ["repeat"]
_SGRNA_KEYS = ["sgrna", "guide rna", "grna", "crrna", "tracrrna", "tracr"]
_ID_KEYS = ["abbrevation", "abbreviation", "name", "ortholog",
            "effector", "uniprot", "nuclease", "id"]


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (s or "").lower()).strip()


def _match(colname: str, keys: list[str]) -> bool:
    n = _norm(colname)
    return any(k in n for k in keys)


def _match_id(colname: str) -> bool:
    """Match an identifier/naming column without false-hitting "nucleotide" /
    "amino acid" (both contain the substring "id")."""
    n = _norm(colname)
    named = any(k in n for k in ["abrev", "abbrev", "name",
                                 "ortholog", "effector", "nuclease", "uniprot"])
    if named:
        return True
    # standalone "id" token, or "id" not inside nucleotide/acid
    if "nucleotide" in n or "amino acid" in n:
        return False
    return " id " in f" {n} " or n == "id" or n.startswith("id ")


def _is_pam_col(colname: str) -> bool:
    """True if the column indicates an experimentally-determined PAM.

    Matches a PAM keyword but excludes computational-prediction columns (see
    _PAM_BLOCKLIST)."""
    n = _norm(colname)
    if any(b in n for b in _PAM_BLOCKLIST):
        return False
    return _match(colname, _PAM_KEYS)


def _looks_like_protein(s: str) -> bool:
    s = (s or "").strip()
    return len(s) >= 80 and bool(re.fullmatch(r"[A-Za-z*]+", s))


def find_table_columns(header: list[str]) -> dict[str, int]:
    """Map semantic role -> column index. Requires a protein-seq or organism
    column (PAM is optional; separated tables may lack a PAM column)."""
    has_aa = any(_match(c, _AA_KEYS) for c in header)
    has_org = any(_match(c, _ORGANISM_KEYS) for c in header)
    if not (has_aa or has_org):
        return {}
    cols: dict[str, int] = {}
    for i, c in enumerate(header):
        cnorm = _norm(c)
        if _is_pam_col(c) and "pam" not in cols:
            cols["pam"] = i
            continue
        if _match(c, _ORGANISM_KEYS) and "organism" not in cols:
            cols["organism"] = i
            continue
        if _match(c, _AA_KEYS) and "aa" not in cols:
            cols["aa"] = i
            continue
        if _match(c, _REPEAT_KEYS) and "repeat" not in cols:
            cols["repeat"] = i
            continue
        if _match(c, _SGRNA_KEYS) and "sgrna" not in cols:
            cols["sgrna"] = i
            continue
        if _match_id(c) and "id" not in cols:
            cols["id"] = i
            continue
    return cols


def _detect_seq_and_id_columns(rows: list[list[str]]) -> dict[str, int]:
    """Heuristic: find a column whose non-empty cells are mostly long alphabetic
    strings (>80 chars) = protein sequence; the leftmost nearby text column is
    the identifier (name/nuclease/uniprot)."""
    n = len(rows) - 1
    if n <= 0:
        return {}
    colwidth = max(len(r) for r in rows)
    best_col = -1
    best_frac = 0.0
    for c in range(colwidth):
        vals = []
        for r in rows[1:]:
            if c < len(r):
                vals.append(str(r[c]))
        nonempty = [v for v in vals if v.strip()]
        if not nonempty:
            continue
        prot = [v for v in nonempty if _looks_like_protein(v)]
        frac = len(prot) / max(1, len(nonempty))
        if frac > best_frac and frac >= 0.5:
            best_frac = frac
            best_col = c
    if best_col < 0:
        return {}
    # id = leftmost column before seq with short text values
    id_col = -1
    for c in range(best_col):
        vals = [str(r[c]) for r in rows[1:] if c < len(r)]
        short = [v for v in vals if v.strip() and len(v) < 60]
        if short:
            id_col = c
    return {"aa": best_col, "id": id_col}


def extract_xlsx(path: Path, pmcid: str) -> list[dict]:
    out: list[dict] = []
    try:
        sheets = list_sheets(path)
    except Exception:
        return out
    for sheet in sheets:
        try:
            rows = read_sheet(path, sheet)
        except Exception:
            continue
        if not rows:
            continue
        header = [str(c) for c in rows[0]]
        cols = find_table_columns(header)
        if not cols:
            # fall back to heuristic sequence detection (Name + Sequence tables)
            cols = _detect_seq_and_id_columns(rows)
        if not cols or "aa" not in cols:
            continue
        for r in rows[1:]:
            def get(k):
                i = cols.get(k)
                return (r[i] if i is not None and i < len(r) else "").strip()
            pam = get("pam")
            aa = get("aa")
            org = get("organism")
            # Only keep if AA looks like protein or organism present
            if aa and not _looks_like_protein(aa):
                aa = ""
            if not aa and not org:
                continue
            out.append({
                "pmcid": pmcid,
                "sheet": sheet,
                "ortholog_id": get("id"),
                "organism": org,
                "pam": pam,
                "repeat": get("repeat"),
                "sgrna": get("sgrna"),
                "protein_sequence": aa,
            })
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True, type=Path)
    ap.add_argument("--output", default="data/parsed/ortholog_tables/ortholog_pam.tsv", type=Path)
    ap.add_argument("--pmcid", type=str, default="")
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    xlsx_files = sorted(args.input_dir.rglob("*.xlsx"))
    print(f"{len(xlsx_files)} xlsx files", flush=True)

    n_records = 0
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pmcid", "sheet", "ortholog_id", "organism", "pam",
                    "repeat", "sgrna", "protein_sequence", "sha256"])
        seen = set()
        for xp in xlsx_files:
            # derive PMCID from filename prefix "PMCxxxxx__..." or the --pmcid flag
            fname = xp.stem
            m = re.match(r"(PMC\d+)", fname)
            pmcid = args.pmcid or (m.group(1) if m else xp.parent.name)
            recs = extract_xlsx(xp, pmcid)
            for rec in recs:
                seq = rec["protein_sequence"] or ""
                # drop empty shells (no PAM and no protein seq, e.g. cutsite tables)
                if not seq and not rec["pam"]:
                    continue
                key = (rec["pmcid"], rec["organism"], rec["pam"], seq[:40])
                if key in seen:
                    continue
                seen.add(key)
                sha = hashlib.sha256(seq.encode()).hexdigest()[:16] if seq else ""
                w.writerow([rec["pmcid"], rec["sheet"], rec["ortholog_id"],
                            rec["organism"], rec["pam"], rec["repeat"],
                            rec["sgrna"], seq, sha])
                n_records += 1
    print(f"total records={n_records} -> {args.output}", flush=True)


if __name__ == "__main__":
    main()
