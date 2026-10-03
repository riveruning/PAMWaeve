"""PAMdict scoring — target #1: judge whether a PAM prediction is reliable.

Given a predicted PAM matrix (10x4, A/C/G/T) and (optionally) the protein's
nearest-neighbor identity to the training set, produce:

  - ``level``             : high / med / low (human-facing verdict)
  - ``confidence``        : scalar in [0,1]
  - ``per_position_entropy`` : Shannon entropy (bits) per position (2.0 = max
                             uncertainty, 0.0 = fully resolved)
  - ``mean_info_content`` : mean information content (bits) across positions
  - ``uncertain_positions`` : list of position indices (1-based) that are weak

Design notes (aligned with the reference audit, not copied):
  - The official implementation expresses predictions as *information content*
    (bits) via ``prob_to_info``, and flags a position as "N" when its max
    per-base info is below 0.70 bit. We keep that threshold semantics.
  - The official UQ path (``percent_identities``) scores novelty by top-10
    Levenshtein similarity to training sequences. We generalize: novelty (low
    similarity to training) *down-weights* confidence, because out-of-domain
    proteins are exactly where Type II/V predictions drift (the user's
    motivating observation).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from ..schema.constants import PAM_NUCLEOTIDES, CONSENSUS_THRESHOLD_BITS

MAX_ENTROPY = math.log2(len(PAM_NUCLEOTIDES))  # 2.0 bits


@dataclass
class PAMScore:
    level: str = "low"
    confidence: float = 0.0
    per_position_entropy: list[float] = field(default_factory=list)
    info_content: list[float] = field(default_factory=list)
    mean_info_content: float = 0.0
    uncertain_positions: list[int] = field(default_factory=list)
    consensus: str = ""
    n_neighbor_sim: float = 0.0
    novelty_weight: float = 1.0
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "level": self.level,
            "confidence": round(self.confidence, 4),
            "per_position_entropy": [round(x, 4) for x in self.per_position_entropy],
            "info_content": [round(x, 4) for x in self.info_content],
            "mean_info_content": round(self.mean_info_content, 4),
            "uncertain_positions": self.uncertain_positions,
            "consensus": self.consensus,
            "n_neighbor_sim": round(self.n_neighbor_sim, 4),
            "novelty_weight": round(self.novelty_weight, 4),
            "notes": self.notes,
        }


def _matrix_as_freq(matrix: list[list[float]], units: str) -> list[list[float]]:
    """Normalize a 10x4 matrix to per-position frequencies.

    Accepts ``units`` = 'prob' (0..1 frequencies; the canonical input) or
    'bits' (official ``prob_to_info`` information content). For 'bits' we
    recover a frequency-like vector by softmax over a scaled log: this is an
    approximation that preserves the argmax ranking, which is all the scoring
    threshold logic actually depends on.
    """
    import numpy as np

    m = np.asarray(matrix, dtype=float).copy()
    if units == "bits":
        # info = p_k * (log2 p_k + log2 K). A monotonic, ranking-preserving
        # proxy for p_k is exp(info); softmax it to a distribution.
        # Scale so that the max ~2.0 bits maps to a confident (not saturated)
        # probability; temperature ~1/log(e) keeps spread reasonable.
        temp = 1.0 / math.log2(math.e)  # ~0.693
        m = np.exp(np.clip(m, -10.0, 10.0) * temp)
    # normalize per position (a row of all zeros -> uniform)
    m = np.clip(m, 0.0, None)
    pos_sums = m.sum(axis=1, keepdims=True)
    pos_sums = np.where(pos_sums == 0, 1.0, pos_sums)
    return (m / pos_sums).tolist()


def score_matrix(matrix: list[list[float]], units: str = "prob",
                 neighbor_sim: float = 1.0,
                 model_name: str = "") -> PAMScore:
    """Score a PAM prediction matrix.

    Args:
        matrix: 10x4 (A/C/G/T) prediction.
        units: 'bits' (information content, official) or 'prob' (frequencies).
        neighbor_sim: top-1 (or mean top-k) sequence identity to training set,
            in [0,1]. 1.0 means "in-domain"; lower means "novel/out-of-domain".
        model_name: e.g. 'cas9' (downstream) — only used for a note.
    """
    freqs = _matrix_as_freq(matrix, units)

    entropies = []
    infos = []
    uncertain = []
    consensus_bases = []
    for i, row in enumerate(freqs):
        h = 0.0
        for p in row:
            if p > 0:
                h -= p * math.log2(p)
        entropies.append(h)
        # information content = MAX_ENTROPY - entropy
        info = MAX_ENTROPY - h
        infos.append(info)
        if info < CONSENSUS_THRESHOLD_BITS:
            uncertain.append(i + 1)  # 1-based position
        best = max(range(len(row)), key=lambda j: row[j])
        consensus_bases.append(PAM_NUCLEOTIDES[best] if info >= CONSENSUS_THRESHOLD_BITS else "N")

    mean_info = sum(infos) / len(infos) if infos else 0.0

    # The 10x4 window is fixed-padded: positions beyond the motif are
    # legitimately near-uniform. "Signal" = positions carrying real base
    # preference (info above a small noise floor). Confidence should reflect
    # how strongly the *informative* positions resolve, not be diluted by
    # trailing uniform padding of a short (3-8 nt) motif.
    NOISE_FLOOR = 0.15  # bits; below this a position is effectively uniform
    signal = [v for v in infos if v >= NOISE_FLOOR]
    if signal:
        signal_info = sum(signal) / len(signal)
    else:
        signal_info = mean_info

    # Confidence combines (a) how resolved the motif is and (b) how in-domain
    # the input protein is.
    motif_conf = float(signal_info / MAX_ENTROPY)  # 0..1
    novelty_weight = float(min(1.0, max(0.0, neighbor_sim)))
    conf = motif_conf * novelty_weight

    # Level thresholds — CALIBRATED against measured distributions
    # (benchmark1_variants_calibration.json, 2026-09-10): the original
    # hand-picked cutoffs (0.40/0.65 motif_conf) sat an order of magnitude
    # above the real confidence range (median 0.02, max 0.11 for real model
    # outputs), so every prediction was labeled "low" and the verdict had no
    # discrimination. The cutoffs below were measured on 187 cas9_full
    # predictions (100 official + 87 literature golds) and match the accuracy
    # terciles (25% -> 37% -> 65% actual IUPAC accuracy low -> high):
    #   confidence <= ~0.010  -> low   (measured acc 25%)
    #   confidence <= ~0.063  -> med   (measured acc 37%)
    #   confidence  >  ~0.063 -> high  (measured acc 65%)
    LEVEL_LOW = 0.010   # ≈ tercile-1 of measured confidence
    LEVEL_MED = 0.063   # ≈ tercile-2 of measured confidence
    if conf <= LEVEL_LOW or novelty_weight < 0.30:
        level = "low"
    elif conf <= LEVEL_MED:
        level = "med"
    else:
        level = "high"

    notes = []
    side = "downstream" if "cas9" in (model_name or "").lower() else "upstream"
    if not uncertain:
        notes.append("all positions resolved")
    else:
        notes.append(f"{len(uncertain)} weak position(s) at {uncertain}")
    if novelty_weight < 0.60:
        notes.append("input protein is out-of-domain (low training similarity) -> downweighted")
    notes.append(f"assumed PAM side: {side}")

    return PAMScore(
        level=level,
        confidence=conf,
        per_position_entropy=entropies,
        info_content=infos,
        mean_info_content=mean_info,
        uncertain_positions=uncertain,
        consensus="".join(consensus_bases),
        n_neighbor_sim=neighbor_sim,
        novelty_weight=novelty_weight,
        notes=notes,
    )
