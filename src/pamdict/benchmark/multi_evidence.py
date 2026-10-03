"""System-level evaluation for protein, spacer, and late-fusion PAM evidence."""
from __future__ import annotations

import itertools
import math
import random
from collections import defaultdict
from typing import Iterable

from pamdict.score.spectrum import exact_match, normalize_pam, valid_pam

REQUIRED_MANIFEST_FIELDS = (
    "system_id",
    "tier",
    "status",
    "cas_family",
    "protein_model",
    "protein_fasta",
    "protein_record_id",
    "spacer_fasta",
    "target_fastas",
    "spacer_orientation",
    "pam_side",
    "pam_lengths",
    "gold_pam_spectrum",
    "gold_evidence_type",
    "protein_training_exposure",
)


def parse_pam_spectrum(value: object) -> dict[int, list[str]]:
    if not isinstance(value, dict) or not value:
        raise ValueError("gold_pam_spectrum must be a non-empty object")
    spectrum: dict[int, list[str]] = {}
    for raw_length, raw_pams in value.items():
        length = int(raw_length)
        if length <= 0 or not isinstance(raw_pams, list):
            raise ValueError("invalid PAM spectrum length group")
        pams = list(dict.fromkeys(normalize_pam(str(pam)) for pam in raw_pams))
        if not pams or any(not valid_pam(pam) or len(pam) != length for pam in pams):
            raise ValueError(f"invalid PAM spectrum for length {length}: {raw_pams}")
        spectrum[length] = pams
    return dict(sorted(spectrum.items()))


def concrete_candidates(
    lengths: Iterable[int],
    *,
    max_candidates: int = 100_000,
) -> list[str]:
    lengths = sorted(set(int(length) for length in lengths))
    if not lengths or any(length <= 0 for length in lengths):
        raise ValueError("PAM lengths must be positive")
    total = sum(4 ** length for length in lengths)
    if total > max_candidates:
        raise ValueError(
            f"exhaustive candidate set would contain {total} sequences; "
            f"limit is {max_candidates}"
        )
    return [
        "".join(values)
        for length in lengths
        for values in itertools.product("ACGT", repeat=length)
    ]


def audit_manifest_row(row: dict[str, object]) -> dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []
    missing = [field for field in REQUIRED_MANIFEST_FIELDS if field not in row]
    if missing:
        errors.append(f"missing fields: {missing}")
    try:
        spectrum = parse_pam_spectrum(row.get("gold_pam_spectrum"))
    except (TypeError, ValueError) as exc:
        spectrum = {}
        errors.append(str(exc))
    lengths = row.get("pam_lengths")
    if not isinstance(lengths, list) or sorted(set(lengths)) != sorted(spectrum):
        errors.append("pam_lengths must equal gold_pam_spectrum length keys")
    if row.get("spacer_orientation") not in ("forward", "reverse"):
        errors.append("spacer_orientation must be forward or reverse")
    if row.get("pam_side") not in ("upstream", "downstream"):
        errors.append("pam_side must be upstream or downstream")
    if not isinstance(row.get("target_fastas"), list) or not row.get("target_fastas"):
        errors.append("target_fastas must be a non-empty list")

    tier = row.get("tier")
    exposure = row.get("protein_training_exposure")
    strict_eligible = (
        tier == "strict_independent"
        and row.get("status") == "ready"
        and exposure == "none"
        and bool(row.get("gold_doi"))
        and row.get("gold_evidence_type") in {
            "experimental_activity",
            "experimental_functional_spectrum",
        }
        and row.get("target_independence") == "virus_clustered"
    )
    if tier == "strict_independent" and not strict_eligible:
        errors.append("strict_independent row does not satisfy independence contract")
    if exposure in ("exact", "near"):
        warnings.append("protein is exposed to Protein2PAM training distribution")
    if row.get("target_independence") != "virus_clustered":
        warnings.append("target records are not verified independent virus clusters")
    if row.get("status") != "ready":
        warnings.append("system is not runnable")
    return {
        "system_id": row.get("system_id", ""),
        "valid": not errors,
        "runnable": not errors and row.get("status") == "ready",
        "strict_independent_eligible": strict_eligible,
        "errors": errors,
        "warnings": warnings,
    }


def _optional_score(value: object) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    score = float(value)
    return score if math.isfinite(score) else None


def _evaluate_rank_values(
    rows: Iterable[dict[str, object]],
    spectrum: dict[int, list[str]],
    value_getter,
) -> dict[str, object]:
    grouped: dict[int, list[tuple[str, float]]] = defaultdict(list)
    for row in rows:
        candidate = normalize_pam(str(row.get("candidate_pam", "")))
        if not candidate or not valid_pam(candidate):
            continue
        value = value_getter(row)
        if value is not None:
            grouped[len(candidate)].append((candidate, value))

    by_length: dict[str, dict[str, object]] = {}
    reciprocal_ranks: list[float] = []
    recalls = {1: [], 3: [], 5: []}
    for length, golds in spectrum.items():
        candidates = grouped.get(length, [])
        gold_values = [
            value
            for candidate, value in candidates
            if any(exact_match(candidate, gold) for gold in golds)
        ]
        if not candidates or not gold_values:
            optimistic_rank = None
            non_gold_ties = 0
            best_rank = None
            reciprocal = 0.0
            top_candidates: list[str] = []
        else:
            best_gold_value = max(gold_values)
            optimistic_rank = 1 + sum(
                value > best_gold_value + 1e-12
                for _candidate, value in candidates
            )
            non_gold_ties = sum(
                abs(value - best_gold_value) <= 1e-12
                and not any(exact_match(candidate, gold) for gold in golds)
                for candidate, value in candidates
            )
            best_rank = optimistic_rank + non_gold_ties
            reciprocal = 1.0 / best_rank
            maximum = max(value for _candidate, value in candidates)
            top_candidates = sorted(
                candidate
                for candidate, value in candidates
                if abs(value - maximum) <= 1e-12
            )
        reciprocal_ranks.append(reciprocal)
        length_recalls = {k: bool(best_rank and best_rank <= k) for k in recalls}
        for k, value in length_recalls.items():
            recalls[k].append(float(value))
        by_length[str(length)] = {
            "gold_pams": golds,
            "scored_candidates": len(candidates),
            "best_gold_rank": best_rank,
            "best_gold_rank_optimistic": optimistic_rank,
            "non_gold_candidates_tied_with_best_gold": non_gold_ties,
            "reciprocal_rank": reciprocal,
            "recall_at_1": length_recalls[1],
            "recall_at_3": length_recalls[3],
            "recall_at_5": length_recalls[5],
            "top_candidates": top_candidates,
        }
    n_lengths = len(spectrum)
    scored_lengths = sum(
        details["scored_candidates"] > 0 for details in by_length.values()
    )
    return {
        "n_gold_lengths": n_lengths,
        "scored_lengths": scored_lengths,
        "abstained_lengths": n_lengths - scored_lengths,
        "coverage": scored_lengths / n_lengths if n_lengths else 0.0,
        "mrr": sum(reciprocal_ranks) / n_lengths if n_lengths else 0.0,
        "recall_at_1": sum(recalls[1]) / n_lengths if n_lengths else 0.0,
        "recall_at_3": sum(recalls[3]) / n_lengths if n_lengths else 0.0,
        "recall_at_5": sum(recalls[5]) / n_lengths if n_lengths else 0.0,
        "by_length": by_length,
    }


def evaluate_scored_rows(
    rows: Iterable[dict[str, object]],
    spectrum: dict[int, list[str]],
    *,
    score_field: str,
) -> dict[str, object]:
    return _evaluate_rank_values(
        rows,
        spectrum,
        lambda row: _optional_score(row.get(score_field)),
    )


def fusion_rank_score(row: dict[str, object]) -> float | None:
    """Deterministic late-fusion rank key; explicitly not a probability."""
    protein = _optional_score(row.get("specificity_adjusted_score"))
    spacer = _optional_score(row.get("spacer_evidence_score"))
    if protein is None or spacer is None:
        return None

    gate = str(row.get("fusion_gate", ""))
    if gate == "spacer_fallback_severe_top_disagreement":
        return spacer
    if gate in {
        "abstain_missing_scored_source",
        "abstain_severe_top_disagreement",
    }:
        return None

    def state(value: float) -> str:
        if value > 60:
            return "support"
        if value < 40:
            return "conflict"
        return "neutral"

    pair = (state(protein), state(spacer))
    classes = {
        ("support", "support"): 5,
        ("support", "neutral"): 4,
        ("neutral", "support"): 4,
        ("neutral", "neutral"): 3,
        ("support", "conflict"): 2,
        ("conflict", "support"): 2,
        ("neutral", "conflict"): 1,
        ("conflict", "neutral"): 1,
        ("conflict", "conflict"): 0,
    }
    evidence_class = classes[pair]
    return evidence_class * 1_000_000 + min(protein, spacer) * 1_000 + max(
        protein, spacer
    )


def evaluate_fusion_rows(
    rows: Iterable[dict[str, object]],
    spectrum: dict[int, list[str]],
) -> dict[str, object]:
    return _evaluate_rank_values(rows, spectrum, fusion_rank_score)


def _bootstrap_mean_ci(
    values: list[float],
    *,
    seed: int,
    iterations: int,
) -> list[float] | None:
    if not values:
        return None
    rng = random.Random(seed)
    boot = sorted(
        sum(rng.choice(values) for _ in values) / len(values)
        for _ in range(iterations)
    )
    return [
        boot[max(0, int(iterations * 0.025) - 1)],
        boot[min(iterations - 1, int(iterations * 0.975))],
    ]


def aggregate_benchmark(
    systems: list[dict[str, object]],
    *,
    seed: int = 17,
    bootstrap_iterations: int = 2000,
) -> dict[str, object]:
    if bootstrap_iterations <= 0:
        raise ValueError("bootstrap_iterations must be positive")
    methods = sorted({
        method
        for system in systems
        for method in system.get("methods", {})
        if not system.get("methods", {})[method].get("not_applicable")
    })
    metrics = ("mrr", "recall_at_1", "recall_at_3", "recall_at_5")
    summaries: dict[str, dict[str, object]] = {}
    for method in methods:
        # A system counts towards ``systems_total`` only when the method is
        # applicable to it.  A channel that a system never supplied
        # (``not_applicable=True``) must be excluded from the denominator
        # rather than averaged in as a zero, otherwise "no evidence channel"
        # would masquerade as "scored badly".
        applicable = [
            system for system in systems
            if method in system.get("methods", {})
            and not system["methods"][method].get("not_applicable")
        ]
        method_rows = [system["methods"][method] for system in applicable]
        eligible = [row for row in method_rows if row["coverage"] > 0]
        not_applicable_count = sum(
            1 for system in systems
            if method in system.get("methods", {})
            and system["methods"][method].get("not_applicable")
        )
        summary: dict[str, object] = {
            "systems_total": len(applicable),
            "systems_in_manifest": len(systems),
            "systems_not_applicable": not_applicable_count,
            "systems_with_method": len(method_rows),
            "systems_scored": len(eligible),
            "system_coverage": len(eligible) / len(applicable) if applicable else 0.0,
        }
        for metric in metrics:
            all_values = [float(row[metric]) for row in method_rows]
            eligible_values = [float(row[metric]) for row in eligible]
            summary[f"{metric}_all_systems"] = (
                sum(all_values) / len(applicable) if applicable else 0.0
            )
            summary[f"{metric}_eligible_systems"] = (
                sum(eligible_values) / len(eligible_values)
                if eligible_values else None
            )
            summary[f"{metric}_bootstrap_95_ci_eligible"] = _bootstrap_mean_ci(
                eligible_values,
                seed=seed,
                iterations=bootstrap_iterations,
            )
        summaries[method] = summary

    paired: dict[str, dict[str, object]] = {}
    for baseline in ("protein_only", "spacer_only"):
        deltas = []
        for system in systems:
            methods_for_system = system.get("methods", {})
            if baseline not in methods_for_system or "fusion" not in methods_for_system:
                continue
            left = methods_for_system[baseline]
            right = methods_for_system["fusion"]
            # Paired comparison only over systems where both methods are
            # applicable and both actually produced a ranked list.
            if left.get("not_applicable") or right.get("not_applicable"):
                continue
            if left["coverage"] > 0 and right["coverage"] > 0:
                deltas.append(float(right["mrr"]) - float(left["mrr"]))
        paired[f"fusion_minus_{baseline}"] = {
            "n": len(deltas),
            "mean_mrr_delta": sum(deltas) / len(deltas) if deltas else None,
            "bootstrap_95_ci": _bootstrap_mean_ci(
                deltas,
                seed=seed,
                iterations=bootstrap_iterations,
            ),
        }
    return {
        "system_count": len(systems),
        "methods": summaries,
        "paired_mrr": paired,
    }
