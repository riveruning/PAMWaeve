"""Convert PAM-depletion count tables into consensus + position frequency matrix.

Large-scale PAM profiling papers (e.g. Wei 2022 eLife) publish, per effector,
a "kmer + read-count" depletion table (typically 5' flank + PAM + spacer). This
module infers the PAM region (the informative positions) and emits:

  * consensus PAM string (IUPAC)
  * per-position A/C/G/T frequency matrix (normalized, the pam_logo_acgt format)

Positions are identified by information content: fixed library flanks are
near-max conservation and are excluded; the PAM window is the run of
conserved-but-varying positions between the flanks and the random spacer tail.
"""
from __future__ import annotations

import math
from collections import Counter

_BASES = "ACGT"


def _bits(freqs: Counter) -> float:
    n = sum(freqs.values())
    if n == 0:
        return 0.0
    h = 0.0
    for f in freqs.values():
        p = f / n
        if p > 0:
            h -= p * math.log2(p)
    return 2.0 - h  # max entropy for 4 bases = 2 bits


def _iupac(freqs: Counter) -> str:
    n = sum(freqs.values())
    if n == 0:
        return "N"
    # two-letter merges
    def p(b):
        return freqs.get(b, 0) / n
    a, c, g, t = p("A"), p("C"), p("G"), p("T")
    # strong/weak/purine/pyrimidine/keto/amino
    if max(a, c, g, t) >= 0.75:
        return max("ACGT", key=lambda b: freqs.get(b, 0))
    # two-base dominates
    pairs = {
        "R": a + g, "Y": c + t, "S": c + g, "W": a + t,
        "K": g + t, "M": a + c,
    }
    for code, q in pairs.items():
        if q >= 0.85:
            return code
    # three-base (V/H/D/B by which is missing)
    triples = {"V": a + c + g, "H": a + c + t, "D": a + g + t, "B": c + g + t}
    for code, q in triples.items():
        if q >= 0.92:
            return code
    return "N"


def count_table_to_pam(kmers: list[str], counts: list[int],
                       flank_bits: float = 1.0) -> dict:
    """Given kmer strings + counts, return consensus, matrix, and PAM region."""
    L = len(kmers[0]) if kmers else 0
    cols = []
    for i in range(L):
        c = Counter()
        for kmer, ct in zip(kmers, counts):
            if i < len(kmer):
                c[kmer[i]] += ct
        cols.append((i, _bits(c), c))

    # Flank positions: near-max bits (>= flank_bits). PAM = informative window
    # between the left flank and the (random) right tail.
    # Heuristic: find the leftmost run of flank positions; PAM starts after it.
    left_flank_end = 0
    for i, bits, c in cols:
        if bits >= flank_bits:
            left_flank_end = i + 1
        else:
            break

    # Right tail: positions with bits below a low threshold (nearly random).
    right_tail_start = L
    for i, bits, c in reversed(cols):
        if bits <= 0.25:
            right_tail_start = i
        else:
            break

    pam_cols = cols[left_flank_end:right_tail_start]
    if not pam_cols:
        # fall back: everything after the flank
        pam_cols = cols[left_flank_end:]

    consensus = "".join(_iupac(c) for _, _, c in pam_cols)

    # build A/C/G/T frequency matrix over PAM positions
    matrix = []
    for _, _, c in pam_cols:
        tot = sum(c.values()) or 1
        matrix.append([round(c.get(b, 0) / tot, 4) for b in _BASES])

    return {
        "consensus": consensus,
        "matrix_acgt": matrix,
        "pam_start": left_flank_end,
        "pam_end": right_tail_start,
        "flank_len": left_flank_end,
        "n_kmers": len(kmers),
        "total_counts": sum(counts),
    }
