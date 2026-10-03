"""Parse "direct PAM table" supplementary sheets.

Some PAM experiments publish, as a plain table, a column of *validated* PAM
sequences (one row per protein/target/organism) instead of a read-count
enrichment library or a pre-aggregated Sequence-logo matrix.  Example (Casδ
paper, Table S4):

    Target ID | Species    | Target gene | Target sequence | Description | PAM
    8N        | N/A        | N/A         | GGTATAACA...   | depletion  | NNNNNNNN
    ...       | H. sapiens | EMX1        | ACTAGGG...     | targeting  | ATG

This module detects such sheets and extracts one schema row per data row that
carries a valid short IUPAC PAM.  It also, when a long protein-sequence column
is present, attaches that sequence so the row is a (protein -> PAM) sample, and
records organism/target as provenance.

The resulting rows use ``pam_consensus`` (the explicit PAM string) rather than
a per-position ``pam_logo_acgt`` matrix; the caller may convert later.
"""

from __future__ import annotations

import re
from typing import Any

from .parse import sanitize_iupac
from ..schema.constants import ALL_COLUMNS

# Column-header synonyms used to locate semantic columns in a direct-PAM table.
# _PAM_HEADERS: EXACT match only (keeps "PAM variant median count" out).
_PAM_HEADERS = ("pam", "pam sequence", "pam_sequence", "pam consensus",
                "protospacer adjacent motif", "target adjacent motif",
                "target adjancent motif", "tam")
# _PAM_HEADERS_PREFIX: PREFIX match — long column labels that get truncated or
# carry annotations (e.g. "Degenerate consensus TAM (if functional)").  Never
# put short tokens ("pam"/"tam") here.
_PAM_HEADERS_PREFIX = ("degenerate consensus", "degenerate_consensus")
_PROTEIN_SEQ_HEADERS = ("protein sequence", "protein_sequence", "amino acid sequence",
                        "aa sequence", "sequence (aa)", "protein amino acid")
_ORGANISM_HEADERS = ("species", "organism", "strain", "source", "host",
                     "bacterial strain", "organism/strain")
_PROTEIN_NAME_HEADERS = ("protein", "protein name", "effector", "effector protein",
                         "nuclease", "cas protein", "cas ortholog", "ortholog",
                         "key name", "name")
_TARGET_HEADERS = ("target", "target gene", "target id", "target_id")


def _cell(x: Any) -> str:
    return str(x).strip()


def _norm_header(x: str) -> str:
    return re.sub(r"[^a-z0-9]", "", x.lower())


def _find_col(headers: list[str], synonyms: tuple[str, ...],
              prefix_synonyms: tuple[str, ...] = ()) -> int | None:
    """Locate a semantic column by (normalized) column name.

    * ``synonyms`` must match the normalized header EXACTLY.
    * ``prefix_synonyms`` may match the normalized header by PREFIX (header
      starts with the synonym, or the synonym starts with the header).  This
      only tolerates header cells that were truncated/annotated (e.g. an XLSX
      cell "Degenerate consensus TAM (if functional)" or "Amino acid sequenc"),
      and is NEVER used for short tokens like "pam"/"tam" — otherwise a stats
      label such as "PAM variant median count" or "Number of unique PAM" would
      be mistaken for a PAM column.
    """
    norms = {_norm_header(h): j for j, h in enumerate(headers)}
    for s in synonyms:
        key = _norm_header(s)
        if key and key in norms:
            return norms[key]
    for s in prefix_synonyms:
        key = _norm_header(s)
        if not key:
            continue
        for h, j in norms.items():
            if not h:
                continue
            # Two prefix directions, each guarded:
            #  (a) header starts with full synonym  -> truncated/annotated form
            #      ("Degenerate consensus TAM" vs "degenerate consensus").
            #  (b) synonym starts with header (header is a shortened form,
            #      e.g. "Amino acid sequenc") — only when the header is long
            #      enough to be a real truncated column name, so a 1-2 char
            #      data value like "DE" can't match "degenerate consensus".
            if h.startswith(key):
                return j
            if len(h) >= 8 and key.startswith(h):
                return j
    return None


def _is_protein_seq(x: str) -> bool:
    """Heuristic: a long amino-acid string (proteins are ~>50 aa)."""
    s = x.strip().upper().replace("*", "").replace("-", "")
    if len(s) < 40:
        return False
    return all(c in "ACDEFGHIKLMNPQRSTVWY" for c in s)


# Full IUPAC DNA alphabet (U normalizes to T). These appear in real, complex
# PAM motifs such as NRRWC / BRTTTTT / NAR(A>G)TC where R/Y/B/D/H/V/W/K/M/S/N
# are standard ambiguity codes.
_IUPAC_DNA = "ACGTURYSWKMBDHVN"
_IUPAC_NORMALIZED = "ACGTT" + "RYWSKMBDHVN"  # U -> T


def _clean_pam_motif(x: str) -> str:
    """Normalize a raw PAM cell into a clean IUPAC motif.

    Strips quotes/whitespace, removes ``(X>Y)`` mutation annotations (e.g.
    ``NR(A>G)TTTT`` -> ``NRTTTT``; the parenthesised part is a per-position
    base-preference note, not part of the motif), then keeps only IUPAC DNA
    letters (U -> T).
    """
    s = x.strip().strip('"').strip("'")
    s = re.sub(r"\([^)]*\)", "", s)  # drop (X>Y) / (A>B) mutation annotations
    s = s.replace("U", "T").upper()
    return "".join(c for c in s if c in _IUPAC_NORMALIZED)


def _is_iupac_pam(x: str) -> bool:
    s = _clean_pam_motif(x)
    if not (2 <= len(s) <= 12):
        return False
    # Reject "all-N" (no specificity) and values that are pure ambiguity with
    # no concrete base at all (e.g. just "N" repeats or "-").
    if not any(c in "ACGT" for c in s):
        return False
    return True


def _is_header_cell(x: str) -> bool:
    """A cell looks like a *header* label, not a data value.

    True for text like "Target", "Species", "PAM", "AA sequence"; False for
    numeric data ("42", "0.783", "1.2E-4"), row labels ("row_5"), and empties.
    """
    s = x.strip()
    if not s:
        return False
    if not any(c.isalpha() for c in s):
        return False
    if re.fullmatch(r"[-+]?[0-9.Ee+-]+", s):
        return False  # scientific / float
    return True


def _is_plain_header_word(x: str) -> bool:
    """A cell is a real header word: >=4 letters, only letters/spaces (no
    underscores, no digits) — distinguishes "Organism"/"Target" from a row
    label like "row_5" or a 2-letter token like "No"/"NA".
    """
    s = x.strip()
    if "_" in s or any(c.isdigit() for c in s):
        return False
    return sum(c.isalpha() for c in s) >= 4


def _find_pam_header_row(rows: list[list[str]]) -> tuple[int, int] | None:
    """Return (pam_col, header_row_idx) for a bona-fide header row.

    A row only counts as the header row if the PAM column sits alongside at
    least one other *plain header word* (e.g. "Target", "Organism") — not a
    numeric cell or a row label like "row_5".  This rejects a data row whose
    only "PAM" match is a gene named "Pam" (GSEA SYMBOL column).
    """
    for i, row in enumerate(rows[:20]):
        if not row:
            continue
        cells = [_cell(c) for c in row]
        j = _find_col(cells, _PAM_HEADERS, prefix_synonyms=_PAM_HEADERS_PREFIX)
        if j is None:
            continue
        if any(_is_plain_header_word(c) for k, c in enumerate(cells) if k != j):
            return j, i
    return None


def find_pam_header(rows: list[list[str]]) -> int | None:
    """Return the column index of the PAM column, or None."""
    hit = _find_pam_header_row(rows)
    return hit[0] if hit else None


def parse_pam_table(rows: list[list[str]]) -> list[dict]:
    """Extract schema rows (dicts) from a direct-PAM table sheet.

    Returns a list of dicts with keys: pam_consensus, protein_sequence,
    protein_id, citation_extra (organism/target provenance). Empty list if the
    sheet does not look like a direct-PAM table.
    """
    hit = _find_pam_header_row(rows)
    if hit is None:
        return []
    pam_col, header_idx = hit

    headers = [_cell(c) for c in rows[header_idx]]
    seq_col = _find_col(headers, _PROTEIN_SEQ_HEADERS,
                        prefix_synonyms=_PROTEIN_SEQ_HEADERS)
    org_col = _find_col(headers, _ORGANISM_HEADERS)
    name_col = _find_col(headers, _PROTEIN_NAME_HEADERS)
    tgt_col = _find_col(headers, _TARGET_HEADERS)

    out: list[dict] = []
    for row in rows[header_idx + 1:]:
        if not row or len(row) <= pam_col:
            continue
        pam_val = _clean_pam_motif(_cell(row[pam_col]))
        if not _is_iupac_pam(pam_val):
            continue

        rec = {c: "" for c in ALL_COLUMNS}
        rec["source"] = "Literature (new; unreviewed)"
        rec["pam_consensus"] = pam_val

        if seq_col is not None and seq_col < len(row) and _is_protein_seq(_cell(row[seq_col])):
            rec["protein_sequence"] = _cell(row[seq_col])

        extra = []
        if name_col is not None and name_col < len(row) and _cell(row[name_col]):
            rec["protein_id"] = _cell(row[name_col])
            extra.append(_cell(row[name_col]))
        if org_col is not None and org_col < len(row) and _cell(row[org_col]):
            extra.append(f"org={_cell(row[org_col])}")
        if tgt_col is not None and tgt_col < len(row) and _cell(row[tgt_col]):
            extra.append(f"target={_cell(row[tgt_col])}")
        if extra:
            rec["raw_ref"] = " | ".join(extra)

        out.append(rec)
    return out
