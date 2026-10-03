# Multi-evidence PAM benchmark report

Generated: 2026-10-03T05:30:21.600522+00:00

Decision: **FRAMEWORK_VALIDATED_NO_STRICT_COHORT**

## Dataset composition

- Manifest systems: 22
- Runnable systems evaluated: 22
- Strict independent systems: 3 {'protein_only': 3}
- Strict paired-evidence systems: 0 (minimum for a generalization claim: 5)
- Input tracks: {'paired_evidence': 1, 'protein_only': 22, 'published_flanks': 18}

## Method summary

| Method | Coverage | MRR (all) | Recall@1 | Recall@3 | Recall@5 | n/a | Applicable |
|---|---:|---:|---:|---:|---:|---:|---:|
| family_prior | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0 | 1 |
| fusion | 0.474 | 0.447 | 0.421 | 0.474 | 0.474 | 3 | 19 |
| protein_only | 1.000 | 0.683 | 0.636 | 0.682 | 0.727 | 0 | 22 |
| spacer_only | 0.947 | 0.700 | 0.579 | 0.842 | 0.842 | 3 | 19 |

## Evaluation-role split

| Role | Systems | Protein MRR | Spacer MRR | Fusion MRR | Fusion coverage |
|---|---:|---:|---:|---:|---:|
| development_regression | 4 | 0.533 | 0.583 | 0.500 | 0.500 |
| retrospective_expansion_2026_09_24 | 6 | 0.754 | 0.833 | 0.583 | 0.667 |
| retrospective_expansion_2026_09_25 | 6 | 1.000 | 0.576 | 0.333 | 0.333 |
| retrospective_holdout_2026_09_24 | 3 | 0.750 | 0.833 | 0.333 | 0.333 |
| strict_protein_only_expansion_2026_09_25 | 1 | 0.036 | n/a | n/a | n/a |
| strict_protein_only_expansion_2026_09_26 | 2 | 0.037 | n/a | n/a | n/a |

## Per-system results

| System | Role | Tier | Spacer evidence | Exposure | Mapped spacers | Orientation | Fusion gate | Protein MRR | Spacer MRR | Fusion MRR |
|---|---|---|---|---|---:|---|---|---:|---:|---:|
| spcas9-pampredict-example | development_regression | engineering_regression | paired_evidence | exact | 37 | pass | joint_ranking_top_agreement | 1.000 | 1.000 | 1.000 |
| sp7f7-published-flanks | development_regression | retrospective_biological | published_flanks | near | 7 | pass | joint_ranking_top_agreement | 1.000 | 1.000 | 1.000 |
| nme-m22809-published-flanks | development_regression | retrospective_biological | published_flanks | exact | 11 | pass | abstain_severe_top_disagreement | 0.067 | 0.333 | 0.000 |
| cj-nctc11168-published-flanks | development_regression | retrospective_biological | published_flanks | exact | 2 | pass | abstain_missing_scored_source | 0.067 | 0.000 | 0.000 |
| tde-atcc35405-published-flanks | retrospective_holdout_2026_09_24 | retrospective_biological | published_flanks | exact | 31 | pass | abstain_severe_top_disagreement | 0.250 | 0.500 | 0.000 |
| sgo-challis-ch1-published-flanks | retrospective_holdout_2026_09_24 | retrospective_biological | published_flanks | exact | 25 | pass | joint_ranking_top_agreement | 1.000 | 1.000 | 1.000 |
| sth1a-lmg18311-published-flanks | retrospective_holdout_2026_09_24 | retrospective_biological | published_flanks | exact | 23 | pass | abstain_severe_top_disagreement | 1.000 | 1.000 | 0.000 |
| ain-ryc-mr95-published-flanks | retrospective_expansion_2026_09_24 | retrospective_biological | published_flanks | exact | 27 | pass | abstain_severe_top_disagreement | 1.000 | 0.500 | 0.000 |
| ssa-jim8777-published-flanks | retrospective_expansion_2026_09_24 | retrospective_biological | published_flanks | exact | 25 | pass | joint_ranking_top_agreement | 0.026 | 1.000 | 1.000 |
| sga-atcc43143-published-flanks | retrospective_expansion_2026_09_24 | retrospective_biological | published_flanks | exact | 23 | pass | joint_ranking_top_agreement | 1.000 | 1.000 | 1.000 |
| smu-gs5-published-flanks | retrospective_expansion_2026_09_24 | retrospective_biological | published_flanks | exact | 18 | pass | joint_ranking_top_agreement | 0.500 | 0.500 | 0.500 |
| sdy-ac2713-published-flanks | retrospective_expansion_2026_09_24 | retrospective_biological | published_flanks | exact | 18 | pass | abstain_severe_top_disagreement | 1.000 | 1.000 | 0.000 |
| sag2-nem316-published-flanks | retrospective_expansion_2026_09_24 | retrospective_biological | published_flanks | exact | 9 | pass | joint_ranking_top_agreement | 1.000 | 1.000 | 1.000 |
| ece-nctc12421-published-flanks | retrospective_expansion_2026_09_25 | retrospective_biological | published_flanks | exact | 33 | pass | joint_ranking_top_agreement | 1.000 | 1.000 | 1.000 |
| tsp-h121-published-flanks | retrospective_expansion_2026_09_25 | retrospective_biological | published_flanks | exact | 14 | pass | abstain_severe_top_disagreement | 1.000 | 0.333 | 0.000 |
| nme2-de10444-published-flanks | retrospective_expansion_2026_09_25 | retrospective_biological | published_flanks | exact | 14 | review | abstain_severe_top_disagreement | 1.000 | 0.000 | 0.000 |
| seq2-nctc11606-published-flanks | retrospective_expansion_2026_09_25 | retrospective_biological | published_flanks | exact | 14 | pass | joint_ranking_top_agreement | 1.000 | 1.000 | 1.000 |
| psp-ct06-published-flanks | retrospective_expansion_2026_09_25 | retrospective_biological | published_flanks | exact | 12 | pass | abstain_severe_top_disagreement | 1.000 | 1.000 | 0.000 |
| wvi-dsm16922-published-flanks | retrospective_expansion_2026_09_25 | retrospective_biological | published_flanks | exact | 8 | pass | abstain_severe_top_disagreement | 1.000 | 0.125 | 0.000 |
| cj4-campylobacter-jejuni-414-protein-only | strict_protein_only_expansion_2026_09_25 | strict_independent | protein_only | none | n/a | n/a | n/a | 0.036 | n/a | n/a |
| asp-acidiphilium-21-60-14-protein-only | strict_protein_only_expansion_2026_09_26 | strict_independent | protein_only | none | n/a | n/a | n/a | 0.008 | n/a | n/a |
| cba-caulobacterales-protein-only | strict_protein_only_expansion_2026_09_26 | strict_independent | protein_only | none | n/a | n/a | n/a | 0.067 | n/a | n/a |

## Method applicability

`spacer_only` and `fusion` are `n/a` for systems that supply no spacer
evidence track (`protein_only`). Those systems are excluded from the
spacer/fusion denominators and from paired fusion-minus-baseline
comparisons, instead of being averaged in as zero scores.


## Interpretation boundary

Ranks are computed over exhaustive concrete A/C/G/T candidates within each PAM length. IUPAC gold motifs accept any concrete subset. Fusion uses an explicit evidence-state rank key, not a calibrated probability.
