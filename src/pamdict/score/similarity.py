"""Nearest-neighbor similarity against a reference sequence set.

Mirrors (not copies) the official ``percent_identities``: estimate how novel an
input protein is by its sequence identity to the training set. For tractability
on CPU we use exact-match k-mer Jaccard as a fast proxy for Levenshtein.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Iterable


def kmer_set(seq: str, k: int = 4) -> set[str]:
    seq = seq.upper()
    if len(seq) < k:
        return {seq}
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def build_kmer_index(seqs: Iterable[str], k: int = 4) -> dict[str, set[int]]:
    """Map k-mer -> set of sequence indices."""
    idx: dict[str, set[int]] = defaultdict(set)
    for i, s in enumerate(seqs):
        for km in kmer_set(s, k):
            idx[km].add(i)
    return idx


def max_neighbor_sim(query: str, seqs: list[str], index: dict[str, set[int]],
                     k: int = 4, topk: int = 10) -> float:
    """Return mean Jaccard similarity of the top-k nearest training sequences."""
    q = kmer_set(query, k)
    if not q:
        return 0.0
    # Candidate seqs sharing any k-mer.
    cand = set()
    for km in q:
        cand |= index.get(km, set())
    if not cand:
        return 0.0
    sims = sorted((jaccard(q, kmer_set(seqs[i], k)) for i in cand), reverse=True)
    top = sims[:topk]
    return sum(top) / len(top) if top else 0.0


def load_sequences(tsv: Path, seq_col: str = "protein_sequence") -> list[str]:
    import csv
    with tsv.open(newline="", encoding="utf-8") as fh:
        return [row[seq_col] for row in csv.DictReader(fh, delimiter="\t")
                if row.get(seq_col)]
