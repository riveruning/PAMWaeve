"""Consensus PAM per Cas alias (domain knowledge, conservative).

Used to filter out mis-attributed PAM in ``cas_pam_pairs.tsv`` — the full-text
extractor often assigns a paper's many PAM motifs to a single well-known Cas
(e.g. SpCas9 ends up with 27 distinct PAM values, but its true PAM is NGG).

Only fully-agreed consensus PAMs are listed. An alias whose consensus PAM we do
not know with certainty is ABSENT here (its PAM must be human-reviewed rather
than auto-filtered).

Equivalence is by IUPAC degenerate-code expansion: a PAM matches its consensus
if every concrete base allowed by the observed string is a subset of the
consensus (i.e. the observed PAM is a degenerate/specific spelling of the
consensus motif), OR the observed string is one of the listed explicit aliases.
"""

from __future__ import annotations

IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"C", "G"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"}, "B": {"C", "G", "T"},
    "D": {"A", "G", "T"}, "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

# alias -> (consensus PAM, explicit equivalent spellings).
# CONSERVATIVE: only entries we are confident about.
CONSENSUS_PAM: dict[str, tuple[str, frozenset[str]]] = {
    "SpCas9": ("NGG", frozenset({"NGG", "GG"})),
    "SpyCas9": ("NGG", frozenset({"NGG", "GG"})),
    "dSpCas9": ("NGG", frozenset({"NGG"})),
    "SaCas9": ("NNGRRT", frozenset({"NNGRRT", "NNGRR"})),
    "dSaCas9": ("NNGRRT", frozenset({"NNGRRT", "NNGRR"})),
    "FnCas9": ("NGG", frozenset({"NGG"})),
    "StCas9": ("NNAGAAW", frozenset({"NNAGAAW", "NNAGAA"})),
    "NmCas9": ("NNNNGATT", frozenset({"NNNNGATT", "NNNNGTTT"})),
    "NmeCas9": ("NNNNGATT", frozenset({"NNNNGATT", "NNNNGTTT"})),
    "CjCas9": ("NNNNRYAC", frozenset({"NNNNRYAC", "NNNNRYYC", "NNNVRYAC"})),
    "AsCas12a": ("TTTV", frozenset({"TTTV", "TTTN", "TTT"})),
    "FnCas12a": ("TTTV", frozenset({"TTTV", "TTTN"})),
}


def _allowed(seq: str) -> list[set[str]] | None:
    """Return per-position allowed bases, or None if any char is not IUPAC."""
    res: list[set[str]] = []
    for ch in seq.upper():
        if ch not in IUPAC:
            return None
        res.append(IUPAC[ch])
    return res


def _is_subset_of(obs: str, ref: str) -> bool:
    """True if every base allowed by `obs` is also allowed by `ref` at each
    position, AND the strings are comparable in length."""
    a = _allowed(obs)
    b = _allowed(ref)
    if a is None or b is None:
        return False
    if len(a) != len(b):
        return False
    return all(x <= y for x, y in zip(a, b))


def matches_consensus(pam: str, cons: str, aliases: frozenset[str]) -> bool:
    """Return True if `pam` is the consensus PAM or a degenerate equivalent."""
    p = pam.upper().strip()
    if not p:
        return False
    if p in aliases:
        return True
    if _is_subset_of(p, cons):
        return True
    # Also allow the consensus to be a subset of the observed (e.g. "NGGG" vs "NGG"
    # shouldn't match, so this is intentionally NOT symmetric — only specific-to-
    # general again via subset). We only check the direct degenerate match above.
    return False


def consensus_for(alias: str) -> tuple[str, frozenset[str]] | None:
    """Return (consensus, aliases) for a Cas alias, or None if unknown."""
    return CONSENSUS_PAM.get(alias)


def filter_pairs(rows: list[dict[str, str]], pam_field: str = "pam",
                 cas_field: str = "cas_name") -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Split rows into (kept, dropped) by consensus-PAM matching.

    Rows whose Cas alias has no consensus entry are KEPT untouched (they must be
    human-reviewed, not silently dropped). Rows for a known Cas whose PAM does
    not match consensus are DROPPED (mis-attribution), with the reason recorded
    in ``filter_reason``.
    """
    kept: list[dict[str, str]] = []
    dropped: list[dict[str, str]] = []
    for r in rows:
        alias = r.get(cas_field, "").strip()
        pam = r.get(pam_field, "").strip()
        info = consensus_for(alias)
        if info is None:
            kept.append(r)  # unknown alias -> keep for human review
            continue
        cons, aliases = info
        if matches_consensus(pam, cons, aliases):
            kept.append(r)
        else:
            d = dict(r)
            d["filter_reason"] = f"PAM {pam!r} not consensus {cons!r} for {alias}"
            dropped.append(d)
    return kept, dropped
