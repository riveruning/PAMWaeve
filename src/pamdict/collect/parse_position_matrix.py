"""Parse per-position base-count matrices (PAM depletion "position weight" tables).

Some papers publish the raw PAM depletion/cleavage counts as a matrix with one
row per base (A/C/G/T) and one column per PAM position, each cell = read count
(e.g. Casδ/Type II-D paper Figure 1c source data):

    b\\p | 1        2        3        4 ...
    A    | 1713621  1434019  1617580 ...
    C    | ...
    G    | ...
    T    | ...

This module detects such a matrix and converts it to per-position A/C/G/T
frequencies (normalized by column), in the same ``{base: [freq per pos]}``
shape returned by ``parse_sequence_logo`` / ``parse_readcount``, so the caller
can feed it straight to ``freq_to_info_content``.
"""

from __future__ import annotations

from typing import Any

from ..schema.constants import PAM_NUCLEOTIDES

_BASE_SYNONYMS = {"a": "A", "c": "C", "g": "G", "t": "T", "u": "T",
                  "adenine": "A", "cytosine": "C", "guanine": "G", "thymine": "T"}


def _cell(x: Any) -> str:
    return str(x).strip()


def _norm_base(x: str) -> str:
    return _BASE_SYNONYMS.get(x.lower(), "")


def _is_number(x: str) -> bool:
    try:
        float(x.replace(",", ""))
        return True
    except (ValueError, AttributeError):
        return False


def parse_position_matrix(rows: list[list[str]]) -> dict[str, list[float]]:
    """Return {base: [freq per position]} from a per-position count matrix.

    Returns {} if no A/C/G/T x positions numeric matrix is found.
    """
    # Find the header row whose cells are position numbers (1..N).
    npos = None
    header_idx = None
    for i, row in enumerate(rows):
        nums = [_cell(c) for c in row if _is_number(_cell(c))]
        # position header = a run of consecutive small integers
        if len(nums) >= 2 and all(_is_number(n) and 1 <= int(float(n)) <= 20 for n in nums):
            header_idx = i
            npos = len(nums)
            break
    if header_idx is None:
        return {}

    # Find base rows below the header.
    accum = {b: [] for b in PAM_NUCLEOTIDES}
    found_any = False
    for row in rows[header_idx + 1:]:
        if not row:
            continue
        base = _norm_base(_cell(row[0])) if row else ""
        if base not in accum:
            continue
        vals = []
        for c in row[1:npos + 1]:
            v = _cell(c).replace(",", "")
            vals.append(float(v) if _is_number(v) else 0.0)
        if len(vals) == npos:
            accum[base] = vals[:npos]
            found_any = True

    if not found_any:
        return {}
    if not all(len(accum[b]) == npos for b in PAM_NUCLEOTIDES if accum[b]):
        return {}

    # Column-normalize to frequencies.
    freq = {b: [] for b in PAM_NUCLEOTIDES}
    for j in range(npos):
        total = sum((accum[b][j] if j < len(accum[b]) else 0.0) for b in PAM_NUCLEOTIDES)
        for b in PAM_NUCLEOTIDES:
            v = accum[b][j] if j < len(accum[b]) else 0.0
            freq[b].append(v / total if total > 0 else 0.0)
    return freq
