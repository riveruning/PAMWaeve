"""Offline tests for the multi-evidence benchmark contract and metrics."""
from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path

from pamdict.benchmark.contract import (
    audit_manifest_document,
    verify_manifest_files,
)
from pamdict.benchmark.multi_evidence import (
    aggregate_benchmark,
    audit_manifest_row,
    concrete_candidates,
    evaluate_fusion_rows,
    evaluate_scored_rows,
    parse_pam_spectrum,
)
from scripts.run_multi_evidence_benchmark import normalize_published_flank_rows
from scripts.score_published_flanks import hits_from_rows


def test_parse_spectrum_and_exhaustive_candidates():
    spectrum = parse_pam_spectrum({"3": ["NGG", "NAG"]})
    assert spectrum == {3: ["NGG", "NAG"]}
    candidates = concrete_candidates([3])
    assert len(candidates) == 64
    assert candidates[0] == "AAA"
    assert candidates[-1] == "TTT"


def test_rank_metrics_use_iupac_gold_and_tie_safe_rank():
    rows = [
        {"candidate_pam": candidate, "score": 0.0}
        for candidate in concrete_candidates([3])
    ]
    for row in rows:
        if row["candidate_pam"] == "AAA":
            row["score"] = 20.0
        if row["candidate_pam"] == "AGG":
            row["score"] = 10.0
    result = evaluate_scored_rows(rows, {3: ["NGG"]}, score_field="score")
    assert result["by_length"]["3"]["best_gold_rank"] == 2
    assert result["recall_at_1"] == 0.0
    assert result["recall_at_3"] == 1.0
    assert result["mrr"] == 0.5


def test_uniform_tie_does_not_receive_free_recall_at_one():
    rows = [
        {"candidate_pam": candidate, "score": 0.0}
        for candidate in concrete_candidates([3])
    ]
    result = evaluate_scored_rows(rows, {3: ["NGG"]}, score_field="score")
    details = result["by_length"]["3"]
    assert details["best_gold_rank_optimistic"] == 1
    assert details["best_gold_rank"] == 61
    assert not details["recall_at_5"]


def test_published_flank_rows_do_not_claim_independent_targets():
    rows = normalize_published_flank_rows([{
        "candidate_pam": "AGG",
        "support_target_count": "4",
        "support_alignment_count": "4",
        "eligible_target_count": "9",
        "eligible_alignment_count": "9",
    }])
    # The record counts are preserved, so nothing is lost.
    assert rows[0]["aggregate_support_record_count"] == "4"
    assert rows[0]["aggregate_eligible_record_count"] == "9"
    # The target/alignment counters must be *unavailable*, not zero.  Writing
    # "0" would assert "no target supported this candidate", which is a
    # different and unsupported claim: aggregate flanks never published raw
    # target identifiers.  An empty string keeps "unknown" distinct from
    # "measured zero".
    assert rows[0]["support_target_count"] == ""
    assert rows[0]["eligible_target_count"] == ""
    assert rows[0]["support_alignment_count"] == ""
    assert rows[0]["eligible_alignment_count"] == ""


def test_published_flank_unavailable_counts_are_not_read_as_zero_targets():
    """A blank target count must not become "2 independent targets" in fusion.

    ``evidence_assessment`` only awards ``concordant_support`` (as opposed to
    the ``_limited_independence`` variant) when at least two spacers *and* two
    targets back a candidate.  With raw targets unavailable, aggregate flanks
    must never be promoted to the full-independence label.
    """
    from pamdict.score.fusion import evidence_assessment

    rows = normalize_published_flank_rows([{
        "candidate_pam": "TGG",
        "support_spacer_count": "5",
        "support_target_count": "5",
        "eligible_target_count": "5",
    }])
    assert rows[0]["support_target_count"] == ""
    status, priority = evidence_assessment(
        90.0,
        80.0,
        support_spacers=int(rows[0]["support_spacer_count"]),
        support_targets=int(rows[0]["support_target_count"] or 0),
    )
    assert status == "concordant_support_limited_independence"
    assert priority == "medium_priority"


def test_published_flank_hit_without_required_side_is_audited_and_skipped():
    rows = [{
        "spacer_id": "spacer-1",
        "source_row_id": "source-7",
        "oriented_spacer": "ACGT",
        "hit": "1",
        "oriented_upstream_flank": "AAAA",
        "oriented_downstream_flank": "",
    }]
    hits, excluded = hits_from_rows(rows, {"spacer-1": "ACGT"}, "downstream")
    assert hits == []
    assert excluded == ["source-7"]


def test_fusion_ranking_rewards_concordance_before_one_sided_score():
    rows = [
        {
            "candidate_pam": "AGG",
            "specificity_adjusted_score": 90.0,
            "spacer_evidence_score": 80.0,
        },
        {
            "candidate_pam": "AAA",
            "specificity_adjusted_score": 99.0,
            "spacer_evidence_score": 20.0,
        },
    ]
    result = evaluate_fusion_rows(rows, {3: ["NGG"]})
    assert result["by_length"]["3"]["best_gold_rank"] == 1
    assert result["recall_at_1"] == 1.0


def test_conservative_fusion_rank_uses_spacer_fallback():
    rows = [
        {
            "candidate_pam": "AAA",
            "specificity_adjusted_score": 90.0,
            "spacer_evidence_score": 10.0,
            "fusion_gate": "spacer_fallback_severe_top_disagreement",
        },
        {
            "candidate_pam": "TTT",
            "specificity_adjusted_score": 10.0,
            "spacer_evidence_score": 90.0,
            "fusion_gate": "spacer_fallback_severe_top_disagreement",
        },
    ]
    result = evaluate_fusion_rows(rows, {3: ["TTT"]})
    assert result["coverage"] == 1.0
    assert result["mrr"] == 1.0
    assert result["recall_at_1"] == 1.0
    assert result["by_length"]["3"]["top_candidates"] == ["TTT"]

    legacy = [dict(row, fusion_gate="legacy_joint_ranking") for row in rows]
    legacy_result = evaluate_fusion_rows(legacy, {3: ["TTT"]})
    assert legacy_result["mrr"] == 0.5

    abstained = [
        dict(row, fusion_gate="abstain_severe_top_disagreement")
        for row in rows
    ]
    abstained_result = evaluate_fusion_rows(abstained, {3: ["TTT"]})
    assert abstained_result["coverage"] == 0.0
    assert abstained_result["mrr"] == 0.0


def test_strict_manifest_contract_rejects_training_exposure():
    row = {
        "system_id": "system-1",
        "tier": "strict_independent",
        "status": "ready",
        "cas_family": "Cas9",
        "protein_model": "cas9_full",
        "protein_fasta": "protein.faa",
        "protein_record_id": "cas9",
        "spacer_fasta": "spacers.fna",
        "target_fastas": ["phages.fna"],
        "spacer_orientation": "forward",
        "pam_side": "downstream",
        "pam_lengths": [3],
        "gold_pam_spectrum": {"3": ["NGG"]},
        "gold_evidence_type": "experimental_functional_spectrum",
        "gold_doi": "10.example/test",
        "protein_training_exposure": "none",
        "target_independence": "virus_clustered",
    }
    accepted = audit_manifest_row(row)
    assert accepted["valid"]
    assert accepted["strict_independent_eligible"]
    row["protein_training_exposure"] = "exact"
    rejected = audit_manifest_row(row)
    assert not rejected["valid"]
    assert not rejected["strict_independent_eligible"]


def test_aggregate_reports_coverage_and_paired_mrr():
    systems = [
        {
            "methods": {
                "protein_only": {"coverage": 1.0, "mrr": 0.5, "recall_at_1": 0.0, "recall_at_3": 1.0, "recall_at_5": 1.0},
                "spacer_only": {"coverage": 1.0, "mrr": 1.0, "recall_at_1": 1.0, "recall_at_3": 1.0, "recall_at_5": 1.0},
                "fusion": {"coverage": 1.0, "mrr": 1.0, "recall_at_1": 1.0, "recall_at_3": 1.0, "recall_at_5": 1.0},
            }
        }
    ]
    result = aggregate_benchmark(systems, bootstrap_iterations=20)
    assert result["methods"]["fusion"]["mrr_all_systems"] == 1.0
    assert result["paired_mrr"]["fusion_minus_protein_only"]["mean_mrr_delta"] == 0.5

def _v2_system() -> dict[str, object]:
    return {
        "system_id": "system-v2",
        "tier": "engineering_regression",
        "status": "ready",
        "crispr_type": "Type II-A",
        "cas_family": "Cas9",
        "host_taxon": "Test bacterium",
        "host_strain": "strain-1",
        "available_tracks": ["protein_only", "paired_evidence"],
        "protein_model": "cas9_full",
        "protein_fasta": "protein.faa",
        "protein_record_id": "cas9",
        "spacer_fasta": "spacers.fna",
        "target_fastas": ["phages.fna"],
        "spacer_orientation": "forward",
        "pam_side": "downstream",
        "pam_lengths": [3],
        "max_mismatches": 2,
        "gold_pam_spectrum": {"3": ["NGG"]},
        "gold_evidence_type": "engineering_known_pam",
        "gold_doi": "10.example/test",
        "gold_assay": {
            "assay_type": "known_reference",
            "activity_values_available": False,
            "conditions": "test",
            "source_locator": "test fixture",
        },
        "gold_system_match": "exact_protein_or_strain",
        "protein_training_exposure": "exact",
        "target_independence": "contig_only",
        "target_dataset": {
            "name": "test targets",
            "release": "v1",
            "frozen_at": "2026-09-22",
            "selection_policy": "predefined fixture",
            "cluster_policy": "not_clustered",
        },
        "split_group": "test-group",
        "provenance": {
            "protein": {
                "source_name": "test",
                "source_uri": "https://example.org/protein",
                "source_version": "v1",
                "license": "test-license",
                "sha256": "a" * 64,
                "sha256_scope": "selected_fasta_record_sequence",
            },
            "spacers": {
                "source_name": "test",
                "source_uri": "https://example.org/spacers",
                "source_version": "v1",
                "license": "test-license",
                "sha256": "b" * 64,
                "sha256_scope": "file_bytes",
            },
            "targets": [{
                "source_name": "test",
                "source_uri": "https://example.org/targets",
                "source_version": "v1",
                "license": "test-license",
                "sha256": "c" * 64,
                "sha256_scope": "file_bytes",
            }],
        },
        "notes": "test fixture",
    }


def _v2_manifest(system: dict[str, object]) -> dict[str, object]:
    return {
        "format_version": 2,
        "benchmark_name": "test benchmark",
        "frozen_at": "2026-09-22",
        "candidate_policy": "exhaustive concrete candidates",
        "minimum_strict_systems": 5,
        "systems": [system],
    }


def test_v2_engineering_manifest_is_valid_but_not_strict():
    audit = audit_manifest_document(_v2_manifest(_v2_system()))
    assert audit["valid"]
    assert audit["track_counts"] == {
        "paired_evidence": 1,
        "protein_only": 1,
    }
    assert audit["strict_independent_track_counts"] == {}
    assert audit["warnings"]


def test_v2_duplicate_system_ids_are_rejected():
    system = _v2_system()
    manifest = _v2_manifest(system)
    manifest["systems"] = [system, dict(system)]
    audit = audit_manifest_document(manifest)
    assert not audit["valid"]
    assert "duplicate system_id" in audit["errors"][0]


def test_v2_strict_paired_track_requires_clustered_targets():
    system = _v2_system()
    system["tier"] = "strict_independent"
    system["protein_training_exposure"] = "none"
    system["gold_evidence_type"] = "experimental_functional_spectrum"
    system["gold_assay"]["assay_type"] = "PAM depletion assay"
    rejected = audit_manifest_document(_v2_manifest(system))
    assert not rejected["valid"]

    system["target_independence"] = "virus_clustered"
    system["target_dataset"]["cluster_policy"] = "95pct-ANI virus clusters"
    accepted = audit_manifest_document(_v2_manifest(system))
    assert accepted["valid"]
    tracks = accepted["system_audits"][0][
        "strict_independent_eligible_tracks"
    ]
    assert tracks == ["protein_only", "paired_evidence"]


def test_v2_protein_only_track_does_not_require_spacer_inputs():
    system = _v2_system()
    system["available_tracks"] = ["protein_only"]
    system["target_independence"] = "not_applicable"
    for field in (
        "spacer_fasta",
        "target_fastas",
        "spacer_orientation",
        "max_mismatches",
        "target_dataset",
    ):
        system.pop(field, None)
    system["provenance"].pop("spacers")
    system["provenance"].pop("targets")
    audit = audit_manifest_document(_v2_manifest(system))
    assert audit["valid"]
    assert audit["track_counts"] == {"protein_only": 1}


def _published_flank_system() -> dict[str, object]:
    system = _v2_system()
    system["available_tracks"] = ["protein_only", "published_flanks"]
    system["spacer_fasta"] = "spacers.fna"
    system["published_flanks_tsv"] = "published_flanks.tsv"
    system["min_effective_spacers"] = 5
    system["target_independence"] = "not_assessable_from_aggregate_flanks"
    system["published_flank_dataset"] = {
        "name": "fixture aggregate flanks",
        "doi": "10.example/flanks",
        "release": "v1",
        "frozen_at": "2026-09-23",
        "row_selection_policy": "fixture accession",
        "orientation_policy": "published PAM orientation",
        "raw_target_identifiers_available": False,
    }
    for field in (
        "target_fastas",
        "spacer_orientation",
        "max_mismatches",
        "target_dataset",
    ):
        system.pop(field, None)
    system["provenance"].pop("targets")
    system["provenance"]["published_flanks"] = {
        "source_name": "test",
        "source_uri": "https://example.org/flanks",
        "source_version": "v1",
        "license": "test-license",
        "sha256": "d" * 64,
        "sha256_scope": "file_bytes",
    }
    return system


def test_v2_published_flanks_are_valid_but_never_strict_paired():
    system = _published_flank_system()
    audit = audit_manifest_document(_v2_manifest(system))
    assert audit["valid"]
    assert audit["track_counts"] == {
        "protein_only": 1,
        "published_flanks": 1,
    }
    assert audit["strict_independent_track_counts"] == {}
    assert "aggregate per-spacer evidence" in " ".join(
        audit["system_audits"][0]["warnings"]
    )


def test_v2_published_flank_file_verification_detects_change():
    system = _published_flank_system()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        (root / "protein.faa").write_text(">cas9\nMABC\n", encoding="utf-8")
        (root / "spacers.fna").write_text(">s1\nACGT\n", encoding="utf-8")
        flank_path = root / "published_flanks.tsv"
        flank_path.write_text("spacer_id\thit\ns1\t1\n", encoding="utf-8")
        system["provenance"]["protein"]["sha256"] = hashlib.sha256(b"MABC").hexdigest()
        system["provenance"]["spacers"]["sha256"] = hashlib.sha256((root / "spacers.fna").read_bytes()).hexdigest()
        system["provenance"]["published_flanks"]["sha256"] = hashlib.sha256(flank_path.read_bytes()).hexdigest()
        accepted = verify_manifest_files(_v2_manifest(system), root)
        assert accepted["valid"]
        assert accepted["checked_artifacts"] == 3
        flank_path.write_text("spacer_id\thit\ns1\t0\n", encoding="utf-8")
        rejected = verify_manifest_files(_v2_manifest(system), root)
        assert not rejected["valid"]
        assert rejected["failed_artifacts"] == 1


def test_complete_abstention_has_defined_rank_fields():
    result = evaluate_scored_rows([], {3: ["NGG"]}, score_field="score")
    details = result["by_length"]["3"]
    assert details["best_gold_rank"] is None
    assert details["best_gold_rank_optimistic"] is None
    assert details["non_gold_candidates_tied_with_best_gold"] == 0
    assert result["coverage"] == 0.0

def test_v2_file_verification_detects_changed_input():
    system = _v2_system()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        protein = root / "protein.faa"
        spacers = root / "spacers.fna"
        targets = root / "phages.fna"
        protein.write_text(">cas9\nMABC\n", encoding="utf-8")
        spacers.write_text(">s1\nACGT\n", encoding="utf-8")
        targets.write_text(">p1\nTTTACGTAGG\n", encoding="utf-8")
        system["provenance"]["protein"]["sha256"] = hashlib.sha256(
            b"MABC"
        ).hexdigest()
        system["provenance"]["spacers"]["sha256"] = hashlib.sha256(
            spacers.read_bytes()
        ).hexdigest()
        system["provenance"]["targets"][0]["sha256"] = hashlib.sha256(
            targets.read_bytes()
        ).hexdigest()
        manifest = _v2_manifest(system)
        accepted = verify_manifest_files(manifest, root)
        assert accepted["valid"]
        assert accepted["checked_artifacts"] == 3

        targets.write_text(">p1\nTTTACGTAAG\n", encoding="utf-8")
        rejected = verify_manifest_files(manifest, root)
        assert not rejected["valid"]
        assert rejected["failed_artifacts"] == 1

def test_v2_rejects_paths_that_escape_the_workspace():
    parent_escape = _v2_system()
    parent_escape["protein_fasta"] = "../secret.faa"
    assert not audit_manifest_document(_v2_manifest(parent_escape))["valid"]

    windows_escape = _v2_system()
    windows_escape["protein_fasta"] = r"C:\\secret.faa"
    assert not audit_manifest_document(_v2_manifest(windows_escape))["valid"]


def test_protein_only_system_has_no_spacer_or_fusion_channel():
    """A protein-only system must not be scored on channels it never supplied.

    Regression for the aggregate denominator: before this was fixed, a system
    with no spacer track was counted as a *zero* for spacer_only and fusion,
    which dragged down every spacer/fusion summary and quietly turned "no
    evidence available" into "scored badly".
    """
    from scripts.run_multi_evidence_benchmark import not_applicable_method

    marker = not_applicable_method("system provides no spacer evidence track")
    assert marker["not_applicable"] is True
    # Metrics are None, never 0.0: absence is not a failed measurement.
    assert marker["mrr"] is None
    assert marker["recall_at_1"] is None
    assert marker["coverage"] == 0.0

    systems = [
        {
            "system_id": "has-spacer",
            "methods": {
                "protein_only": {"coverage": 1.0, "mrr": 0.5,
                                 "recall_at_1": 0.5, "recall_at_3": 0.5,
                                 "recall_at_5": 0.5},
                "spacer_only": {"coverage": 0.0, "mrr": 0.0, "recall_at_1": 0.0,
                                "recall_at_3": 0.0, "recall_at_5": 0.0},
                "fusion": {"coverage": 0.0, "mrr": 0.0, "recall_at_1": 0.0,
                           "recall_at_3": 0.0, "recall_at_5": 0.0},
            },
        },
        {
            "system_id": "protein-only",
            "methods": {
                "protein_only": {"coverage": 1.0, "mrr": 0.5,
                                 "recall_at_1": 0.5, "recall_at_3": 0.5,
                                 "recall_at_5": 0.5},
                "spacer_only": marker,
                "fusion": not_applicable_method("nothing to fuse"),
            },
        },
    ]
    aggregate = aggregate_benchmark(systems, bootstrap_iterations=10)
    fusion = aggregate["methods"]["fusion"]
    # The protein-only system is excluded from the denominator, not scored 0.
    assert fusion["systems_in_manifest"] == 2
    assert fusion["systems_total"] == 1
    assert fusion["systems_not_applicable"] == 1
    # And it is excluded from the paired fusion-vs-protein comparison.
    assert aggregate["paired_mrr"]["fusion_minus_protein_only"]["n"] == 0


def test_fusion_paired_comparison_excludes_not_applicable_systems():
    """Paired deltas must only use systems where both methods produced a rank."""
    from scripts.run_multi_evidence_benchmark import not_applicable_method

    systems = [
        {
            "system_id": "both",
            "methods": {
                "protein_only": {"coverage": 1.0, "mrr": 0.2, "recall_at_1": 0.0,
                                 "recall_at_3": 0.0, "recall_at_5": 0.0},
                "spacer_only": {"coverage": 1.0, "mrr": 1.0, "recall_at_1": 1.0,
                                "recall_at_3": 1.0, "recall_at_5": 1.0},
                "fusion": {"coverage": 1.0, "mrr": 0.8, "recall_at_1": 1.0,
                           "recall_at_3": 1.0, "recall_at_5": 1.0},
            },
        },
        {
            "system_id": "protein-only",
            "methods": {
                "protein_only": {"coverage": 1.0, "mrr": 0.9, "recall_at_1": 1.0,
                                 "recall_at_3": 1.0, "recall_at_5": 1.0},
                "spacer_only": not_applicable_method("no spacer track"),
                "fusion": not_applicable_method("no spacer track"),
            },
        },
    ]
    aggregate = aggregate_benchmark(systems, bootstrap_iterations=10)
    paired = aggregate["paired_mrr"]["fusion_minus_protein_only"]
    assert paired["n"] == 1
    assert abs(paired["mean_mrr_delta"] - 0.6) < 1e-9


def test_report_renders_unavailable_metrics_as_na_not_zero():
    """The markdown report must not print 0.000 for a non-applicable method."""
    from scripts.run_multi_evidence_benchmark import _fmt

    assert _fmt(None) == "n/a"
    assert _fmt(0.0) == "0.000"
    assert _fmt(0.5) == "0.500"


def test_rel_paths_keep_absolute_dev_paths_out_of_reports():
    """Reports must record repo-relative paths, not the build machine layout."""
    from scripts.run_multi_evidence_benchmark import WORKSPACE, rel

    inside = WORKSPACE / "data" / "parsed" / "x.json"
    assert rel(inside) == "data/parsed/x.json"
    assert not rel(inside).startswith("/")
    # A path outside the repo is passed through rather than mangled.
    assert rel(Path("/definitely/outside.json")) == "/definitely/outside.json"
