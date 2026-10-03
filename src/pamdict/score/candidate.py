"""Score candidate IUPAC PAMs against a Protein2PAM probability matrix.

The score measures model compatibility, not cleavage activity.  It compares
the probability mass assigned to each IUPAC symbol with the mass expected
under a uniform four-base background.  This prevents broad symbols such as N
from receiving free credit merely because they accept every nucleotide.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, replace
from typing import Iterable

import numpy as np

from .spectrum import allowed_set, normalize_pam, pam_information, valid_pam

NUCLEOTIDES = "ACGT"


@dataclass(frozen=True)
class CandidatePAMScore:
    candidate_pam: str
    pam_length: int
    informative_positions: int
    allowed_probability_geomean: float | None
    specificity_adjusted_score: float | None
    specificity_adjusted_evidence_bits: float | None
    candidate_information_bits_per_position: float
    rank_within_length: int | None = None

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        for key, value in list(result.items()):
            if isinstance(value, float):
                result[key] = round(value, 6)
        return result


def candidate_window(length: int, side: str, window: int = 10) -> range:
    """Return matrix positions occupied by a PAM of the requested length."""
    if not 1 <= length <= window:
        raise ValueError(f"PAM length must be in [1, {window}], got {length}")
    if side == "downstream":
        return range(0, length)
    if side == "upstream":
        return range(window - length, window)
    raise ValueError("side must be 'downstream' or 'upstream'")


def _normalize_probability_matrix(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    if values.shape != (10, 4):
        raise ValueError(
            f"probability matrix must have shape (10, 4), got {values.shape}"
        )
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("probability matrix must contain finite non-negative values")
    totals = values.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("each probability-matrix row must have positive mass")
    return values / totals


def score_candidate_pam(
    probability_matrix: np.ndarray,
    candidate_pam: str,
    *,
    side: str = "downstream",
) -> CandidatePAMScore:
    """Score one IUPAC PAM against a 10x4 A/C/G/T probability matrix.

    specificity_adjusted_score is bounded to [0, 100]:
      - 50 means the model assigns the same mass as a uniform background;
      - above 50 means support;
      - below 50 means contradiction.

    N positions carry zero information and are excluded.  An all-N candidate
    therefore has no score instead of receiving an artificial perfect score.
    """
    pam = normalize_pam(candidate_pam)
    if not valid_pam(pam):
        raise ValueError(f"invalid IUPAC PAM: {candidate_pam!r}")
    matrix = _normalize_probability_matrix(probability_matrix)
    positions = candidate_window(len(pam), side, window=matrix.shape[0])

    masses: list[float] = []
    bounded_supports: list[float] = []
    log_enrichments: list[float] = []
    weights: list[float] = []
    eps = 1e-12
    for symbol, position in zip(pam, positions):
        allowed = allowed_set(symbol)
        size = len(allowed)
        weight = math.log2(4 / size)
        if weight == 0:
            continue
        mass = sum(matrix[position, NUCLEOTIDES.index(base)] for base in allowed)
        background = size / 4
        if mass >= background:
            bounded = (mass - background) / (1 - background)
        else:
            bounded = (mass - background) / background
        masses.append(float(mass))
        bounded_supports.append(float(bounded))
        log_enrichments.append(math.log2(max(mass, eps) / background))
        weights.append(weight)

    if not weights:
        return CandidatePAMScore(
            candidate_pam=pam,
            pam_length=len(pam),
            informative_positions=0,
            allowed_probability_geomean=None,
            specificity_adjusted_score=None,
            specificity_adjusted_evidence_bits=None,
            candidate_information_bits_per_position=pam_information(pam),
        )

    weight_sum = sum(weights)
    mean_support = sum(
        value * weight for value, weight in zip(bounded_supports, weights)
    ) / weight_sum
    evidence_bits = sum(
        value * weight for value, weight in zip(log_enrichments, weights)
    ) / weight_sum
    geomean = math.exp(sum(math.log(max(value, eps)) for value in masses) / len(masses))
    return CandidatePAMScore(
        candidate_pam=pam,
        pam_length=len(pam),
        informative_positions=len(weights),
        allowed_probability_geomean=geomean,
        specificity_adjusted_score=50 * (mean_support + 1),
        specificity_adjusted_evidence_bits=evidence_bits,
        candidate_information_bits_per_position=pam_information(pam),
    )


def rank_candidate_pams(
    probability_matrix: np.ndarray,
    candidates: Iterable[str],
    *,
    side: str = "downstream",
) -> list[CandidatePAMScore]:
    """Score candidates and assign ranks only among motifs of equal length."""
    unique = list(dict.fromkeys(normalize_pam(value) for value in candidates))
    scored = [
        score_candidate_pam(probability_matrix, value, side=side)
        for value in unique
    ]
    ranks: dict[str, int | None] = {}
    lengths = sorted({score.pam_length for score in scored})
    for length in lengths:
        group = [
            score for score in scored
            if score.pam_length == length
            and score.specificity_adjusted_score is not None
        ]
        group.sort(
            key=lambda score: (
                -float(score.specificity_adjusted_score),
                -float(score.allowed_probability_geomean),
                score.candidate_pam,
            )
        )
        for rank, score in enumerate(group, start=1):
            ranks[score.candidate_pam] = rank
    return [
        replace(score, rank_within_length=ranks.get(score.candidate_pam))
        for score in scored
    ]


def probability_information_stats(
    probability_matrix: np.ndarray,
    *,
    signal_floor_bits: float = 0.15,
) -> dict[str, object]:
    """Return transparent entropy/information diagnostics for one prediction."""
    matrix = _normalize_probability_matrix(probability_matrix)
    safe = np.clip(matrix, np.finfo(np.float64).tiny, 1.0)
    entropy = -np.sum(matrix * np.log2(safe), axis=1)
    information = 2.0 - entropy
    signal = information[information >= signal_floor_bits]
    return {
        "mean_information_bits_all_positions": float(information.mean()),
        "mean_information_bits_signal_positions": (
            float(signal.mean()) if len(signal) else float(information.mean())
        ),
        "signal_positions": [
            int(index + 1)
            for index, value in enumerate(information)
            if value >= signal_floor_bits
        ],
    }
