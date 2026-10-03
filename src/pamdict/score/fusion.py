"""Transparent late fusion of protein and spacer PAM evidence."""
from __future__ import annotations

import itertools
from typing import Iterable

from .spectrum import normalize_pam


def _optional_float(value: object) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def _optional_int(value: object) -> int:
    if value is None or str(value).strip() == "":
        return 0
    return int(float(value))


def _minimum_top_distance_within_tolerance(
    left: set[str],
    right: set[str],
    tolerance: int,
) -> int | None:
    """Return the minimum Hamming distance up to ``tolerance``.

    Candidate spaces can contain many tied maxima. Generating bounded
    neighbours of the smaller top set avoids an unbounded Cartesian product.
    ``None`` means every top-candidate pair is farther apart than the declared
    tolerance; the exact larger distance is intentionally not needed by the
    gate.
    """
    overlap = left & right
    if overlap:
        return 0
    source, target = (left, right) if len(left) <= len(right) else (right, left)
    for distance in range(1, tolerance + 1):
        for candidate in source:
            for positions in itertools.combinations(
                range(len(candidate)), distance
            ):
                originals = {
                    position: candidate[position] for position in positions
                }
                for replacements in itertools.product("ACGT", repeat=distance):
                    if any(
                        replacement == originals[position]
                        for position, replacement in zip(
                            positions, replacements
                        )
                    ):
                        continue
                    mutated = list(candidate)
                    for position, replacement in zip(positions, replacements):
                        mutated[position] = replacement
                    if "".join(mutated) in target:
                        return distance
    return None


def _annotate_fusion_gate(
    rows: list[dict[str, object]],
    *,
    policy: str,
) -> None:
    """Annotate a per-length heuristic without reading gold at runtime.

    Developed on the initial four systems; results on those systems are
    development regressions, not independent validation.
    """
    grouped: dict[int, list[dict[str, object]]] = {}
    for row in rows:
        candidate = normalize_pam(str(row.get("candidate_pam", "")))
        if candidate:
            grouped.setdefault(len(candidate), []).append(row)

    for length, length_rows in grouped.items():
        protein_values = [
            (normalize_pam(str(row["candidate_pam"])), value)
            for row in length_rows
            if (value := _optional_float(
                row.get("specificity_adjusted_score")
            )) is not None
        ]
        spacer_values = [
            (normalize_pam(str(row["candidate_pam"])), value)
            for row in length_rows
            if (value := _optional_float(
                row.get("spacer_evidence_score")
            )) is not None
        ]
        tolerance = max(1, length // 4)
        top_distance: int | None = None
        if not protein_values or not spacer_values:
            gate = "abstain_missing_scored_source"
        elif policy == "legacy":
            gate = "legacy_joint_ranking"
        else:
            protein_max = max(value for _candidate, value in protein_values)
            spacer_max = max(value for _candidate, value in spacer_values)
            protein_top = {
                candidate
                for candidate, value in protein_values
                if abs(value - protein_max) <= 1e-12
            }
            spacer_top = {
                candidate
                for candidate, value in spacer_values
                if abs(value - spacer_max) <= 1e-12
            }
            top_distance = _minimum_top_distance_within_tolerance(
                protein_top,
                spacer_top,
                tolerance,
            )
            if top_distance is not None:
                gate = "joint_ranking_top_agreement"
            elif policy == "spacer_fallback":
                gate = "spacer_fallback_severe_top_disagreement"
            else:
                gate = "abstain_severe_top_disagreement"
        for row in length_rows:
            row["fusion_gate"] = gate
            row["fusion_top_hamming_tolerance"] = tolerance
            row["fusion_top_distance_within_tolerance"] = (
                "" if top_distance is None else top_distance
            )


def evidence_assessment(
    protein_score: float | None,
    spacer_score: float | None,
    *,
    support_spacers: int = 0,
    support_targets: int = 0,
) -> tuple[str, str]:
    """Return an interpretable status and review priority, not a probability.

    The 40/60 band is deliberately treated as neutral for triage.  Scores are
    centered at 50 but are not calibrated probabilities, so a tiny excursion
    above 50 must not become a high-priority biological claim.
    """
    if protein_score is None and spacer_score is None:
        return "no_scored_evidence", "insufficient_evidence"
    if spacer_score is None:
        return "protein_only_unverified", "needs_spacer_evidence"
    if protein_score is None:
        return "spacer_only_unverified", "needs_protein_evidence"

    protein_state = "support" if protein_score > 60 else (
        "conflict" if protein_score < 40 else "neutral"
    )
    spacer_state = "support" if spacer_score > 60 else (
        "conflict" if spacer_score < 40 else "neutral"
    )
    if protein_state == "support" and spacer_state == "support":
        if support_spacers >= 2 and support_targets >= 2:
            return "concordant_support", "high_priority"
        return "concordant_support_limited_independence", "medium_priority"
    if protein_state == "conflict" and spacer_state == "conflict":
        return "concordant_conflict", "low_priority"
    if protein_state != spacer_state and "neutral" not in (
        protein_state,
        spacer_state,
    ):
        return (
            f"discordant_{protein_state}_protein_{spacer_state}_spacer",
            "manual_review",
        )
    return f"protein_{protein_state}_spacer_{spacer_state}", "uncertain"


def fuse_pam_evidence(
    protein_rows: Iterable[dict[str, object]],
    spacer_rows: Iterable[dict[str, object]],
    *,
    policy: str = "abstain",
) -> list[dict[str, object]]:
    """Join exact candidate strings while preserving both source scores."""
    if policy not in {"abstain", "spacer_fallback", "legacy"}:
        raise ValueError("policy must be abstain, spacer_fallback, or legacy")
    spacer_by_candidate: dict[str, dict[str, object]] = {}
    for row in spacer_rows:
        candidate = normalize_pam(str(row.get("candidate_pam", "")))
        if not candidate:
            continue
        if candidate in spacer_by_candidate:
            raise ValueError(
                f"duplicate spacer evidence for candidate {candidate}; "
                "filter to one system before fusion"
            )
        spacer_by_candidate[candidate] = dict(row)

    output: list[dict[str, object]] = []
    for raw_protein in protein_rows:
        protein = dict(raw_protein)
        candidate = normalize_pam(str(protein.get("candidate_pam", "")))
        if not candidate:
            continue
        spacer = spacer_by_candidate.get(candidate)
        merged = dict(protein)
        if spacer is not None:
            for key, value in spacer.items():
                if key == "candidate_pam":
                    continue
                target_key = key if key.startswith("spacer_") else f"spacer_{key}"
                if target_key in merged:
                    target_key = f"spacer_source_{key}"
                merged[target_key] = value

        protein_score = _optional_float(
            protein.get("specificity_adjusted_score")
        )
        spacer_score = _optional_float(
            spacer.get("spacer_evidence_score") if spacer else None
        )
        support_spacers = _optional_int(
            spacer.get("support_spacer_count") if spacer else None
        )
        support_targets = _optional_int(
            spacer.get("support_target_count") if spacer else None
        )
        status, priority = evidence_assessment(
            protein_score,
            spacer_score,
            support_spacers=support_spacers,
            support_targets=support_targets,
        )
        merged["evidence_status"] = status
        merged["evidence_priority"] = priority
        merged["fusion_interpretation"] = (
            "late-fusion evidence categories; no calibrated joint probability"
        )
        output.append(merged)
    _annotate_fusion_gate(output, policy=policy)
    return output
