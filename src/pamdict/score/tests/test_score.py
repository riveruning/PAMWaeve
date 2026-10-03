"""Offline tests for the scoring module (target #1)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from pamdict.score.scoring import score_matrix  # noqa: E402
from pamdict.score.candidate import (  # noqa: E402
    rank_candidate_pams,
    score_candidate_pam,
)
from pamdict.score.similarity import kmer_set, jaccard, max_neighbor_sim, build_kmer_index  # noqa: E402
from pamdict.score.spectrum import (  # noqa: E402
    best_spectrum_cosine,
    exact_match,
    pam_expansion_size,
    pam_information,
    spectrum_cover,
    spectrum_of,
    stable_benchmark_row_ids,
)


def test_downstream_each_position_and_acgt_column_mapping():
    import numpy as np
    for position in range(10):
        for column, base in enumerate('ACGT'):
            matrix = np.full((10, 4), .25)
            matrix[position] = .01
            matrix[position, column] = .97
            pam = 'N' * position + base
            score = score_candidate_pam(matrix, pam, side='downstream')
            assert score.specificity_adjusted_score > 95
            wrong = 'N' * position + 'ACGT'[(column + 1) % 4]
            assert score_candidate_pam(matrix, wrong, side='downstream').specificity_adjusted_score < 5
            if position:
                assert score_candidate_pam(matrix, 'N'*(position-1)+base, side='downstream').specificity_adjusted_score == 50


def test_downstream_trailing_ns_do_not_shift_informative_positions():
    import numpy as np
    matrix = np.full((10, 4), .25)
    matrix[4] = [.01, .01, .97, .01]
    a = score_candidate_pam(matrix, 'NNNNG', side='downstream')
    b = score_candidate_pam(matrix, 'NNNNGNN', side='downstream')
    assert a.specificity_adjusted_score == b.specificity_adjusted_score
    assert a.specificity_adjusted_evidence_bits == b.specificity_adjusted_evidence_bits


def _freq_matrix(consensus: str) -> list[list[float]]:
    """Build a frequency (probability) matrix for a consensus motif.

    At each position, the consensus base gets high frequency (others low).
    Returns a full 10-position matrix padded with uniform (uncertain) rows.
    """
    bases = "ACGT"
    matrix = []
    for ch in consensus.upper():
        if ch == "N":
            matrix.append([0.25, 0.25, 0.25, 0.25])
        else:
            row = [0.02, 0.02, 0.02, 0.02]
            row[bases.index(ch)] = 0.94
            matrix.append(row)
    while len(matrix) < 10:
        matrix.append([0.25, 0.25, 0.25, 0.25])
    return matrix


def test_perfect_motif_high():
    # Fully-resolved 10-position motif (all positions confident).
    m = _freq_matrix("NGGNNNNNNN")
    # replace N positions with resolved ones so the whole 10-pos motif is strong
    m = _freq_matrix("TGGCATTACG")  # 10 distinct but fully resolved
    s = score_matrix(m, units="prob", neighbor_sim=1.0)
    assert s.level == "high"
    assert s.confidence > 0.7


def test_uniform_motif_low():
    uni = [[0.25, 0.25, 0.25, 0.25]] * 10
    s = score_matrix(uni, units="prob", neighbor_sim=1.0)
    assert s.level == "low"
    assert s.mean_info_content < 0.1


def test_novelty_downweights():
    strong = _freq_matrix("TGGCATTACG")
    in_domain = score_matrix(strong, units="prob", neighbor_sim=1.0)
    out_domain = score_matrix(strong, units="prob", neighbor_sim=0.1)
    assert out_domain.confidence < in_domain.confidence
    assert out_domain.level in ("med", "low")


def test_uncertain_positions_flagged():
    m = _freq_matrix("NGG")
    # position 0 is 'N' (uniform) -> flagged uncertain (1-based pos 1)
    m[0] = [0.25, 0.25, 0.25, 0.25]
    s = score_matrix(m, units="prob", neighbor_sim=1.0)
    assert 1 in s.uncertain_positions
    assert s.consensus[0] == "N"


def test_kmer_jaccard():
    assert jaccard({"A", "B"}, {"A", "B"}) == 1.0
    assert jaccard({"A"}, {"B"}) == 0.0


def test_max_neighbor_sim_exact():
    seqs = ["ACGTACGT", "TTTTGGGG"]
    idx = build_kmer_index(seqs, k=3)
    sim = max_neighbor_sim("ACGTACGT", seqs, idx, k=3)
    assert sim > 0.9  # exact match to first seq
    sim2 = max_neighbor_sim("AAAAAAAA", seqs, idx, k=3)
    assert sim2 < 0.9


def test_spectrum_accepts_any_observed_pam_once():
    spectrum = spectrum_of(["NGG", "NAG", "NGG", "NNGRRT"])
    assert spectrum == {3: ["NGG", "NAG"], 6: ["NNGRRT"]}
    assert spectrum_cover("AGG", spectrum)
    assert spectrum_cover("AAG", spectrum)
    assert not spectrum_cover("NNN", spectrum)
    assert not spectrum_cover("NNGRRT", {3: ["NGG"]})


def test_spectrum_cosine_requires_matching_length():
    spectrum = spectrum_of(["NGG", "NNGRRT"])
    assert best_spectrum_cosine("AGG", spectrum) > 0.0
    assert best_spectrum_cosine("AGGT", spectrum) == 0.0


def test_information_penalizes_broad_predictions():
    assert pam_information("AGG") == 2.0
    assert pam_information("NGG") < pam_information("AGG")
    assert pam_information("NNN") == 0.0
    assert pam_expansion_size("NGG") == 4
    assert pam_expansion_size("NNN") == 64


def test_iupac_subset_semantics():
    assert exact_match("AGG", "NGG")
    assert exact_match("NRG", "NNG")
    assert not exact_match("NNG", "NRG")


def test_stable_benchmark_row_ids_handle_exact_duplicates():
    rows = [{"protein_sequence": "AAA", "pam_consensus": "NGG"}] * 2
    ids = stable_benchmark_row_ids(rows)
    assert len(ids) == 2
    assert len(set(ids)) == 2


def test_candidate_score_is_neutral_on_uniform_background():
    matrix = [[0.25, 0.25, 0.25, 0.25] for _ in range(10)]
    result = score_candidate_pam(matrix, "NGG")
    assert result.specificity_adjusted_score == 50.0
    assert result.specificity_adjusted_evidence_bits == 0.0


def test_candidate_score_prefers_supported_iupac_pam():
    matrix = [[0.25, 0.25, 0.25, 0.25] for _ in range(10)]
    matrix[1] = [0.02, 0.02, 0.94, 0.02]
    matrix[2] = [0.02, 0.02, 0.94, 0.02]
    scores = {
        value.candidate_pam: value
        for value in rank_candidate_pams(matrix, ["NGG", "NAA", "NNN"])
    }
    assert scores["NGG"].specificity_adjusted_score > 90
    assert scores["NAA"].specificity_adjusted_score < 10
    assert scores["NGG"].rank_within_length == 1
    assert scores["NNN"].specificity_adjusted_score is None
    assert scores["NNN"].rank_within_length is None


def test_candidate_alignment_depends_on_pam_side():
    matrix = [[0.25, 0.25, 0.25, 0.25] for _ in range(10)]
    matrix[7] = [0.02, 0.02, 0.02, 0.94]
    matrix[8] = [0.02, 0.02, 0.02, 0.94]
    upstream = score_candidate_pam(matrix, "TTN", side="upstream")
    downstream = score_candidate_pam(matrix, "TTN", side="downstream")
    assert upstream.specificity_adjusted_score > 90
    assert downstream.specificity_adjusted_score == 50.0
