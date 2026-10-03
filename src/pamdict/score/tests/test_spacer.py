"""Offline tests for spacer-derived PAM evidence and late fusion."""
from __future__ import annotations

from pamdict.score.fusion import evidence_assessment, fuse_pam_evidence
from pamdict.score.spacer import (
    background_base_frequencies,
    find_hamming_matches,
    find_spacer_hits,
    pam_from_hit,
    reverse_complement,
    score_spacer_candidates,
)


def test_reverse_complement_supports_iupac():
    assert reverse_complement("ACGTRY") == "RYACGT"


def test_hamming_match_is_complete_for_allowed_mismatches():
    assert find_hamming_matches("TTTACGTTAAA", "ACGCT", 1) == [(3, 1)]
    assert find_hamming_matches("TTTACGTTAAA", "ACGCT", 0) == []


def test_oriented_flanks_are_correct_on_both_strands():
    forward_spacer = "ACGTC"
    forward_hits = find_spacer_hits(
        [("forward", forward_spacer)],
        [("target", f"TTT{forward_spacer}AGGAAA")],
        flank_length=3,
    )
    assert len(forward_hits) == 1
    assert forward_hits[0].strand == "+"
    assert pam_from_hit(forward_hits[0], "downstream", 3) == "AGG"

    reverse_spacer = "AACCG"
    reverse_target = f"AAACCT{reverse_complement(reverse_spacer)}TTA"
    reverse_hits = find_spacer_hits(
        [("reverse", reverse_spacer)],
        [("target", reverse_target)],
        flank_length=3,
    )
    assert len(reverse_hits) == 1
    assert reverse_hits[0].strand == "-"
    assert pam_from_hit(reverse_hits[0], "downstream", 3) == "AGG"


def test_spacer_candidate_score_prefers_repeated_supported_pam():
    spacer_1 = "ACGTC"
    spacer_2 = "GGTAC"
    target = f"TTT{spacer_1}AGGAAA{spacer_2}AGGTTT"
    targets = [("phage", target)]
    hits = find_spacer_hits(
        [("s1", spacer_1), ("s2", spacer_2)],
        targets,
        flank_length=3,
    )
    scores = {
        score.candidate_pam: score
        for score in score_spacer_candidates(
            hits,
            ["NGG", "NAA", "NNN"],
            side="downstream",
            base_frequencies=background_base_frequencies(targets),
        )
    }
    assert scores["NGG"].support_spacer_count == 2
    assert scores["NGG"].spacer_evidence_score > scores["NAA"].spacer_evidence_score
    assert scores["NGG"].spacer_rank_within_length == 1
    assert scores["NNN"].spacer_evidence_score is None


def test_late_fusion_preserves_sources_and_flags_agreement():
    status, priority = evidence_assessment(
        80.0,
        75.0,
        support_spacers=3,
        support_targets=2,
    )
    assert status == "concordant_support"
    assert priority == "high_priority"
    rows = fuse_pam_evidence(
        [{"protein_id": "cas9", "candidate_pam": "NGG", "specificity_adjusted_score": "80"}],
        [{
            "system_id": "system-1",
            "candidate_pam": "NGG",
            "spacer_evidence_score": "75",
            "support_spacer_count": "3",
            "support_target_count": "2",
        }],
    )
    assert rows[0]["specificity_adjusted_score"] == "80"
    assert rows[0]["spacer_evidence_score"] == "75"
    assert rows[0]["evidence_status"] == "concordant_support"
    assert rows[0]["fusion_gate"] == "joint_ranking_top_agreement"
    weak_status, weak_priority = evidence_assessment(
        52.0,
        51.0,
        support_spacers=10,
        support_targets=10,
    )
    assert weak_status == "protein_neutral_spacer_neutral"
    assert weak_priority == "uncertain"


def test_fusion_abstains_by_default_on_severe_top_disagreement():
    protein_rows = [
        {"candidate_pam": "AAAA", "specificity_adjusted_score": "90"},
        {"candidate_pam": "AAAT", "specificity_adjusted_score": "80"},
        {"candidate_pam": "TTTA", "specificity_adjusted_score": "20"},
        {"candidate_pam": "TTTT", "specificity_adjusted_score": "10"},
    ]
    spacer_rows = [
        {"candidate_pam": "AAAA", "spacer_evidence_score": "10"},
        {"candidate_pam": "AAAT", "spacer_evidence_score": "20"},
        {"candidate_pam": "TTTA", "spacer_evidence_score": "80"},
        {"candidate_pam": "TTTT", "spacer_evidence_score": "90"},
    ]
    conservative = fuse_pam_evidence(protein_rows, spacer_rows)
    assert {
        row["fusion_gate"] for row in conservative
    } == {"abstain_severe_top_disagreement"}
    assert {
        row["fusion_top_hamming_tolerance"] for row in conservative
    } == {1}
    assert {
        row["fusion_top_distance_within_tolerance"] for row in conservative
    } == {""}

    exploratory = fuse_pam_evidence(
        protein_rows,
        spacer_rows,
        policy="spacer_fallback",
    )
    assert {
        row["fusion_gate"] for row in exploratory
    } == {"spacer_fallback_severe_top_disagreement"}

    legacy = fuse_pam_evidence(
        protein_rows,
        spacer_rows,
        policy="legacy",
    )
    assert {
        row["fusion_gate"] for row in legacy
    } == {"legacy_joint_ranking"}
