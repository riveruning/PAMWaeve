"""PAMdict scoring package (target #1: reliability judgment)."""

from .scoring import PAMScore, score_matrix
from .candidate import (
    CandidatePAMScore,
    probability_information_stats,
    rank_candidate_pams,
    score_candidate_pam,
)
from .fusion import evidence_assessment, fuse_pam_evidence
from .spacer import (
    SpacerHit,
    SpacerPAMScore,
    find_spacer_hits,
    score_spacer_candidates,
)

__all__ = [
    "PAMScore",
    "score_matrix",
    "CandidatePAMScore",
    "probability_information_stats",
    "rank_candidate_pams",
    "score_candidate_pam",
    "evidence_assessment",
    "fuse_pam_evidence",
    "SpacerHit",
    "SpacerPAMScore",
    "find_spacer_hits",
    "score_spacer_candidates",
]
