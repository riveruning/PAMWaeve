"""Parse high-throughput PAM determination (HT-PAMDA style) read-count tables.

Many recent PAM-profiling papers publish their raw library as a supplementary
XLSX where each row is one PAM oligo (e.g. a 3-nt sequence, all 4^3 = 64) with
per-replicate read counts for an unselected ("Unbound"/"Uncut") and a selected
("Bound"/"Cut") channel, plus normalized read counts and an enrichment/depletion
ratio.  This is a *different* format from the pre-aggregated "Sequence logo"
per-position frequency matrices handled by ``parse_xlsx_supp``.

This module detects such tables and converts them into a per-position 10x4
A/C/G/T matrix in the same normalized form the rest of the pipeline expects:

  1. detect the "positive"/selected channel (Bound for binding assays, Cut for
     cleavage assays) and its replicate count columns;
  2. weight each library sequence by its summed selected-channel counts (raw
     abundance of the selected fraction - the standard PWM-from-library signal);
  3. for each position, accumulate the weight onto the base at that position,
     normalize to a frequency (A/C/G/T sum to 1);
  4. pad to 10 positions (positions beyond the library length are left as
     uniform/zero) and return a frequency matrix.

The caller is responsible for the final freq -> bits conversion via
``freq_to_info_content`` and for recording the provenance/units, exactly as for
the Sequence-logo path.
"""

from __future__ import annotations

import math
from typing import Any

from ..schema.constants import PAM_NUCLEOTIDES, PAM_POSITIONS

# Row-0 header markers that identify a read-count/enrichment table.
_READCOUNT_MARKERS = ("read counts", "normalized read counts", "enrichment",
                      "depletion")
# Second-level channel labels for the "selected" (positive) channel.
_BOUND_ALIASES = ("bound", "cut", "cleaved", "selected", "enriched")
_UNBOUND_ALIASES = ("unbound", "uncut", "uncleaved", "unselected", "depleted")


def _cell(x: Any) -> str:
    return str(x).strip()


def _is_readcount_table(rows: list[list[str]]) -> bool:
    """True if the sheet looks like an HT-PAMDA read-count/enrichment table."""
    if not rows:
        return False
    joined = " ".join(_cell(c).lower() for row in rows[:3] for c in row)
    # Need a PAM/oligo label column AND a read-count marker.
    has_count = any(m in joined for m in _READCOUNT_MARKERS)
    has_label = any(_cell(r[0]).lower() == "pam" for r in rows[:4] if r)
    return has_count and has_label


def _norm_oligo(x: str) -> str:
    """Normalize a library oligo label to uppercase IUPAC (T for U).

    Only IUPAC letters are kept; anything else makes the label non-viable so a
    header label like "PAM" (which would otherwise degrade to "A") is rejected.
    """
    s = x.strip().upper().replace("U", "T")
    if not s or any(ch not in "ACGTN" for ch in s):
        return ""
    return s


def _find_channel_columns(rows: list[list[str]]) -> tuple[list[int], list[int], str]:
    """Locate replicate count column indices; return (positive, negative, mode).

    Scans the first few header rows for the channel labels. ``positive`` are the
    selected-channel columns (Bound/Cut), ``negative`` the unselected channel.
    Returns empty lists if the layout is not recognized.
    """
    pos, neg = [], []
    mode = "binding"
    # Header rows 0..3 may hold the merged 'Read Counts' / '#1' / channel labels.
    for row in rows[:4]:
        for j, c in enumerate(row):
            cl = _cell(c).lower()
            if cl in _BOUND_ALIASES:
                if j not in pos:
                    pos.append(j)
                if cl in ("cut", "cleaved"):
                    mode = "cleavage"
            elif cl in _UNBOUND_ALIASES:
                if j not in neg:
                    neg.append(j)
    return pos, neg, mode


def _is_oligo_row(row: list[str]) -> bool:
    """True if the row is a data row (col 0 is an oligo of pure A/C/G/T)."""
    if not row:
        return False
    oligo = _norm_oligo(row[0])
    if len(oligo) < 1:
        return False
    return all(ch in "ACGT" for ch in oligo)


def _find_mean_column(rows: list[list[str]]) -> int | None:
    """Return the column index of the enrichment/depletion 'Mean' summary
    column, or None if the table only has raw counts.

    The Mean column is the paper's per-sequence enrichment (or depletion)
    ratio averaged over replicates - the direct PAM-specificity signal. When
    present it is a better PWM weight than raw selected-channel counts, which
    also encode library base-composition bias.
    """
    for row in rows[:4]:
        for j, c in enumerate(row):
            if _cell(c).lower() == "mean":
                return j
    return None


def _parse_float(x: str) -> float:
    try:
        return float(x)
    except (ValueError, TypeError):
        return 0.0


def parse_readcount(rows: list[list[str]]) -> dict[str, list[float]]:
    """Convert an HT-PAMDA read-count table into per-position frequencies.

    Returns ``{base: [freq per position]}`` mapping over A/C/G/T, positions
    clipped to ``PAM_POSITIONS``, exactly like ``parse_sequence_logo`` so the
    two parsers are interchangeable downstream.

    When an enrichment/depletion "Mean" column is present it is used as the
    per-sequence weight (the direct PAM-specificity signal); otherwise the
    summed selected-channel (Bound/Cut) counts are used.

    Returns ``{}`` if the sheet is not recognized as a read-count table or the
    layout cannot be resolved.
    """
    freq, _mode = parse_readcount_weighted(rows)
    return freq


def parse_readcount_weighted(rows: list[list[str]]) -> tuple[dict[str, list[float]], str]:
    """Like ``parse_readcount`` but also returns the weight mode used.

    Weight mode is one of ``mean_enrichment`` (paper's Mean ratio column),
    ``selected_counts`` (summed Bound/Cut counts), or ``none`` (unrecognized).
    """
    if not _is_readcount_table(rows):
        return {}, "none"

    pos_cols, _neg_cols, _mode = _find_channel_columns(rows)
    mean_col = _find_mean_column(rows)
    if not pos_cols and mean_col is None:
        return {}, "none"

    # Find the first data (oligo) row: header depth varies across papers, so do
    # not hardcode row-index 4; scan from below the header block.
    data_start = 4
    for i in range(len(rows)):
        if _is_oligo_row(rows[i]):
            data_start = i
            break

    # What is the library length? Infer from the first valid oligo row.
    lib_len = None
    weights: dict[str, float] = {}
    for row in rows[data_start:]:
        if not _is_oligo_row(row):
            continue
        oligo = _norm_oligo(row[0])
        if lib_len is None:
            lib_len = len(oligo)
        if len(oligo) != lib_len:
            # Sanity: keep only rows matching the dominant library length.
            continue
        if mean_col is not None:
            weight = _parse_float(_cell(row[mean_col]) if mean_col < len(row) else "")
        else:
            weight = sum(_parse_float(_cell(row[j])) for j in pos_cols)
        weights[oligo] = weights.get(oligo, 0.0) + weight

    if lib_len is None or not weights:
        return {}, "none"

    # Negative/zero enrichment weights would zero a column; clip at 0 so a
    # fully-depleted position is uniform rather than NaN/negative.
    lib_len = min(max(lib_len, 1), len(PAM_POSITIONS))
    accum: dict[str, list[float]] = {b: [0.0] * lib_len for b in PAM_NUCLEOTIDES}
    for oligo, w in weights.items():
        w = max(w, 0.0)
        for i, ch in enumerate(oligo[:lib_len]):
            if ch in accum:
                accum[ch][i] += w

    # Normalize each position to a frequency vector (A/C/G/T sum = 1).
    freq: dict[str, list[float]] = {b: [] for b in PAM_NUCLEOTIDES}
    for i in range(lib_len):
        total = sum(accum[b][i] for b in PAM_NUCLEOTIDES)
        for b in PAM_NUCLEOTIDES:
            freq[b].append(accum[b][i] / total if total > 0 else 0.0)

    weight_mode = "mean_enrichment" if mean_col is not None else "selected_counts"
    return freq, weight_mode


def readcount_freq_to_matrix(freq: dict[str, list[float]]) -> list[list[float]]:
    """Convert {base:[freq per pos]} (clipped to 10) into a 10x4 freq matrix."""
    npos = min(max((len(v) for v in freq.values()), default=0), len(PAM_POSITIONS))
    m = [[0.0] * 4 for _ in range(len(PAM_POSITIONS))]
    for j, b in enumerate(PAM_NUCLEOTIDES):
        for i, v in enumerate(freq.get(b, [])[:npos]):
            m[i][j] = v
    return m
