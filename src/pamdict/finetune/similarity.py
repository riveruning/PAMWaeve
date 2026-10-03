"""Protein-sequence similarity helpers used by fine-tune leak guards.

The previous benchmark builder used a DNA/IUPAC alphabet for protein 3-mers.
This module uses sparse string k-mers over the canonical protein alphabet,
then confirms candidate matches with global edit identity. The k-mer score is
only a fast candidate filter; it is never the final 90% identity definition.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Sequence

PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWY")


def normalize_protein(sequence: str) -> str:
    return sequence.strip().upper()


def protein_kmer_counts(sequence: str, k: int = 3) -> Counter[str]:
    """Sparse counts of canonical amino-acid k-mers."""
    if k <= 0:
        raise ValueError("k must be positive")
    sequence = normalize_protein(sequence)
    counts: Counter[str] = Counter()
    for i in range(max(0, len(sequence) - k + 1)):
        token = sequence[i:i + k]
        if set(token) <= PROTEIN_ALPHABET:
            counts[token] += 1
    return counts


def sparse_cosine(left: Counter[str], right: Counter[str]) -> float:
    if not left or not right:
        return 0.0
    if len(left) > len(right):
        left, right = right, left
    dot = sum(value * right.get(key, 0) for key, value in left.items())
    if dot == 0:
        return 0.0
    left_norm = math.sqrt(sum(value * value for value in left.values()))
    right_norm = math.sqrt(sum(value * value for value in right.values()))
    return dot / (left_norm * right_norm)


def banded_edit_distance(left: str, right: str, max_distance: int) -> int | None:
    """Global Levenshtein distance, or None when the cutoff is exceeded.

    Myers' bit-vector recurrence performs the exact calculation with Python
    big integers. It is substantially faster than a Python nested-loop DP for
    the 1,000--2,000 aa proteins in this project.
    """
    if max_distance < 0:
        raise ValueError("max_distance must be non-negative")
    left = normalize_protein(left)
    right = normalize_protein(right)
    if abs(len(left) - len(right)) > max_distance:
        return None
    if len(left) > len(right):
        left, right = right, left
    if not left:
        return len(right) if len(right) <= max_distance else None
    character_masks: dict[str, int] = defaultdict(int)
    for index, character in enumerate(left):
        character_masks[character] |= 1 << index
    positive = (1 << len(left)) - 1
    negative = 0
    distance = len(left)
    last_bit = 1 << (len(left) - 1)
    for character in right:
        equal = character_masks.get(character, 0)
        vertical = equal | negative
        horizontal = (((equal & positive) + positive) ^ positive) | equal
        positive_horizontal = negative | ~(horizontal | positive)
        negative_horizontal = positive & horizontal
        if positive_horizontal & last_bit:
            distance += 1
        elif negative_horizontal & last_bit:
            distance -= 1
        positive_horizontal = (positive_horizontal << 1) | 1
        negative_horizontal <<= 1
        positive = negative_horizontal | ~(vertical | positive_horizontal)
        negative = positive_horizontal & vertical
    return distance if distance <= max_distance else None


def global_edit_identity(left: str, right: str, threshold: float = 0.90) -> float | None:
    """Return 1-edit_distance/max_length if it meets the threshold."""
    if not 0 < threshold <= 1:
        raise ValueError("threshold must be in (0, 1]")
    denominator = max(len(left), len(right))
    if denominator == 0:
        return 1.0
    max_distance = math.floor((1.0 - threshold) * denominator + 1e-9)
    distance = banded_edit_distance(left, right, max_distance)
    return None if distance is None else 1.0 - distance / denominator


class ProteinNeighborIndex:
    """Reusable sparse 3-mer index with exact edit-identity confirmation."""

    def __init__(self, references: Sequence[str], k: int = 3):
        self.k = k
        self.references = list(dict.fromkeys(
            normalize_protein(sequence)
            for sequence in references
            if sequence.strip()
        ))
        self.counts = [
            protein_kmer_counts(sequence, k) for sequence in self.references
        ]
        self.norms = [
            math.sqrt(sum(value * value for value in counts.values()))
            for counts in self.counts
        ]
        self.postings: dict[str, list[tuple[int, int]]] = defaultdict(list)
        self.exact: dict[str, list[int]] = defaultdict(list)
        for index, (sequence, counts) in enumerate(
            zip(self.references, self.counts)
        ):
            self.exact[sequence].append(index)
            for token, count in counts.items():
                self.postings[token].append((index, count))

    def find(
        self,
        queries: Sequence[str],
        *,
        identity_threshold: float = 0.90,
        candidate_cosine: float = 0.30,
        max_candidates: int = 512,
        probe_kmers: int | None = None,
    ) -> dict[str, list[dict[str, float | int | str]]]:
        if max_candidates <= 0:
            raise ValueError("max_candidates must be positive")
        if probe_kmers is not None and probe_kmers <= 0:
            raise ValueError("probe_kmers must be positive when provided")
        result: dict[str, list[dict[str, float | int | str]]] = {}
        for raw_query in dict.fromkeys(queries):
            query = normalize_protein(raw_query)
            counts = protein_kmer_counts(query, self.k)
            if probe_kmers is None:
                candidate_indices: set[int] = set()
                tokens = counts
            else:
                tokens = sorted(
                    counts,
                    key=lambda token: (len(self.postings.get(token, ())), token),
                )[:probe_kmers]
                candidate_indices = {
                    index
                    for token in tokens
                    for index, _count in self.postings.get(token, ())
                }
            dots: dict[int, int] = defaultdict(int)
            for token in tokens:
                query_count = counts[token]
                for index, ref_count in self.postings.get(token, ()):
                    dots[index] += query_count * ref_count
            if probe_kmers is None:
                candidate_indices.update(dots)
            scored: list[tuple[float, int]] = []
            for index in candidate_indices:
                cosine = sparse_cosine(counts, self.counts[index])
                if cosine >= candidate_cosine:
                    scored.append((cosine, index))
            scored.sort(reverse=True)
            ranked_indices = [
                index for _score, index in scored[:max_candidates]
            ]
            for index in self.exact.get(query, ()):
                if index not in ranked_indices:
                    ranked_indices.insert(0, index)

            matches: list[dict[str, float | int | str]] = []
            for index in ranked_indices:
                reference = self.references[index]
                identity = global_edit_identity(
                    query, reference, identity_threshold
                )
                if identity is None:
                    continue
                distance = round(
                    (1.0 - identity) * max(len(query), len(reference))
                )
                matches.append({
                    "reference_sequence": reference,
                    "identity": round(identity, 8),
                    "edit_distance": distance,
                    "kmer_cosine": round(
                        sparse_cosine(counts, self.counts[index]), 8
                    ),
                })
            matches.sort(
                key=lambda item: (
                    -float(item["identity"]),
                    -float(item["kmer_cosine"]),
                    str(item["reference_sequence"]),
                )
            )
            result[query] = matches
        return result


def find_near_neighbors(
    queries: Sequence[str],
    references: Sequence[str],
    *,
    identity_threshold: float = 0.90,
    k: int = 3,
    candidate_cosine: float = 0.30,
    max_candidates: int = 512,
    probe_kmers: int | None = None,
) -> dict[str, list[dict[str, float | int | str]]]:
    """Find reference sequences within the global-edit-identity threshold."""
    return ProteinNeighborIndex(references, k=k).find(
        queries,
        identity_threshold=identity_threshold,
        candidate_cosine=candidate_cosine,
        max_candidates=max_candidates,
        probe_kmers=probe_kmers,
    )
