"""System-level PAM spectrum utilities.

The benchmark contains one row per reported PAM observation, while the
biological prediction target is a CRISPR system with a *set* of functional
PAMs. This module keeps those two grains explicit: raw rows remain intact,
and scoring aggregates them by protein sequence before calculating a result.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from typing import Iterable


IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"C", "G"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def normalize_pam(pam: str) -> str:
    return (pam or "").strip().upper()


def allowed_set(ch: str) -> set[str]:
    return IUPAC.get(ch.upper(), set())


def valid_pam(pam: str) -> bool:
    p = normalize_pam(pam)
    return bool(p) and all(ch in IUPAC for ch in p)


def exact_match(pred: str, gold: str) -> bool:
    """Return whether a prediction is an equal-length subset of the gold."""
    pred = normalize_pam(pred)
    gold = normalize_pam(gold)
    if len(pred) != len(gold) or not pred:
        return False
    for p, g in zip(pred, gold):
        ap, ag = allowed_set(p), allowed_set(g)
        if not ap or not ag or not ap <= ag:
            return False
    return True


def spectrum_of(pams: Iterable[str]) -> dict[int, list[str]]:
    """Group unique, valid PAM observations by motif length."""
    grouped: dict[int, list[str]] = defaultdict(list)
    seen: set[tuple[int, str]] = set()
    for raw in pams:
        pam = normalize_pam(raw)
        key = (len(pam), pam)
        if not valid_pam(pam) or key in seen:
            continue
        grouped[len(pam)].append(pam)
        seen.add(key)
    return {length: values for length, values in sorted(grouped.items())}


def spectrum_cover(pred: str, spectrum: dict[int, list[str]]) -> bool:
    pred = normalize_pam(pred)
    return any(exact_match(pred, gold) for gold in spectrum.get(len(pred), []))


def pam_matrix(pam: str) -> list[list[float]]:
    """Convert an IUPAC PAM to the legacy 10x5 augmented information matrix."""
    matrix = [[0.0] * 5 for _ in range(10)]
    for i, ch in enumerate(normalize_pam(pam)[:10]):
        allowed = allowed_set(ch)
        n = len(allowed)
        if n == 1:
            matrix[i]["ACGT".index(next(iter(allowed)))] = 2.0
        elif 0 < n < 4:
            for base in allowed:
                matrix[i]["ACGT".index(base)] = 2.0 / n
        matrix[i][4] = (
            0.0 if n == 1 else 1.0 if n == 4 else 2.0 * (1 - 1 / n)
        )
    return matrix


def cosine(a: list[list[float]], b: list[list[float]]) -> float:
    dot = sum(x * y for ra, rb in zip(a, b) for x, y in zip(ra, rb))
    na = math.sqrt(sum(x * x for row in a for x in row))
    nb = math.sqrt(sum(x * x for row in b for x in row))
    return dot / (na * nb) if na and nb else 0.0


def best_spectrum_cosine(
    pred: str, spectrum: dict[int, list[str]], *, same_length: bool = True
) -> float:
    """Return the best cosine to an observed spectrum member.

    Equal-length comparison is the primary definition. It prevents a short,
    broad motif from receiving credit for a longer measured motif merely
    because both matrices are padded to ten positions.
    """
    pred = normalize_pam(pred)
    if not valid_pam(pred):
        return 0.0
    golds = (
        spectrum.get(len(pred), [])
        if same_length
        else [pam for values in spectrum.values() for pam in values]
    )
    return max(
        (cosine(pam_matrix(pred), pam_matrix(gold)) for gold in golds),
        default=0.0,
    )


def pam_information(pam: str) -> float:
    """Mean IUPAC information (bits/position); N=0 and a fixed base=2."""
    pam = normalize_pam(pam)
    if not valid_pam(pam):
        return 0.0
    return sum(2.0 - math.log2(len(allowed_set(ch))) for ch in pam) / len(pam)


def pam_expansion_size(pam: str) -> int:
    """Number of concrete DNA strings represented by an IUPAC motif."""
    pam = normalize_pam(pam)
    if not valid_pam(pam):
        return 0
    return math.prod(len(allowed_set(ch)) for ch in pam)


def sequence_sha256(sequence: str) -> str:
    normalized = "".join((sequence or "").split()).upper()
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def system_id(sequence: str, family: str = "") -> str:
    prefix = "".join(ch.lower() for ch in (family or "system") if ch.isalnum())
    return f"{prefix or 'system'}-{sequence_sha256(sequence)[:16]}"


def stable_benchmark_row_ids(rows: list[dict]) -> list[str]:
    """Build deterministic row IDs without mutating the source benchmark."""
    occurrences: dict[str, int] = defaultdict(int)
    out = []
    for row in rows:
        supplied = (row.get("row_id") or "").strip()
        if supplied:
            out.append(supplied)
            continue
        payload = json.dumps(
            row, sort_keys=True, ensure_ascii=False, separators=(",", ":")
        )
        digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
        occurrences[digest] += 1
        out.append(f"bench-{digest}-{occurrences[digest]:03d}")
    if len(out) != len(set(out)):
        raise ValueError("benchmark row IDs are not unique")
    return out
