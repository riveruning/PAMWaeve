"""Parse XLSX supplementary "Sequence logo" source tables into PAM logo matrices.

Many PAM papers publish per-position A/C/G/T frequencies as "Source data"
sheets (see PAM-readID MOESM5). This extracts those matrices and normalizes
them into the schema's 10x4 A/C/G/T layout (padding missing positions with 0).

Important: these source tables are typically **frequencies (0..1)**, NOT the
information-content (bits) used by the official ``pam_logo_acgt``. The caller
must record units and convert (freq -> bits) explicitly via
``freq_to_info_content``.
"""

from __future__ import annotations

import math
from typing import Any

from .xlsx_reader import read_sheet, list_sheets
from ..schema.constants import PAM_NUCLEOTIDES

_BASE_SYNONYMS = {
    "a": "A", "c": "C", "g": "G", "t": "T", "u": "T",
    "adenine": "A", "cytosine": "C", "guanine": "G", "thymine": "T",
}


def _norm_base(x: str) -> str:
    return _BASE_SYNONYMS.get(x.strip().lower(), "")


def find_logo_sheets(path: str) -> list[str]:
    return [s for s in list_sheets(path) if "microsoft" not in s]


def parse_sequence_logo(rows: list[list[str]]) -> dict[str, list[float]]:
    """Extract a per-position A/C/G/T matrix from a 'Sequence logo' sheet.

    Returns {base: [position frequencies...]} for A/C/G/T present in the table.
    """
    # Locate the header row containing the word 'Base'.
    header_idx = None
    for i, row in enumerate(rows):
        if any(str(c).strip().lower() == "base" for c in row):
            header_idx = i
            break
    if header_idx is None:
        return {}

    # Position columns = the numeric cells after 'Base' in the header.
    header = rows[header_idx]
    pos_cols = []
    for j, c in enumerate(header):
        val = str(c).strip()
        if val.isdigit():
            pos_cols.append(j)

    matrix: dict[str, list[float]] = {b: [] for b in PAM_NUCLEOTIDES}
    for row in rows[header_idx + 1:]:
        base = _norm_base(str(row[0])) if row else ""
        if base not in matrix:
            continue
        if matrix[base]:  # already saw this base (avoid duplicate rows)
            continue
        for j in pos_cols:
            try:
                matrix[base].append(float(str(row[j])))
            except (ValueError, IndexError):
                matrix[base].append(0.0)
    return matrix


def freq_to_info_content(freqs: list[list[float]]) -> list[list[float]]:
    """Convert per-position frequency rows (A/C/G/T) to information content (bits).

    information = p * (log2(p) + log2(K)) summed nuance follows the reference
    ``prob_to_info``: per-position per-base weighted information in bits.
    """
    K = len(PAM_NUCLEOTIDES)
    out = []
    for row in freqs:
        total = sum(row)
        p = [v / total if total > 0 else 0.0 for v in row]
        bits = []
        for pk in p:
            if pk <= 0:
                bits.append(0.0)
            else:
                # p_k * (log2(p_k) + log2(K))
                bits.append(pk * (math.log2(pk) + math.log2(K)))
        out.append(bits)
    return out
