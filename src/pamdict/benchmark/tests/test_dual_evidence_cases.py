"""Task C regression tests for the dual-evidence case infrastructure.

Offline: needs no model weights, no GPU and no network. Covers the parts of the
case pipeline that must not silently drift:

* case files are well formed and declare an orientation decision plus its basis;
* the orientation probe genuinely distinguishes the two orientations;
* per-hit target strand is never used as orientation evidence;
* the two evidence classes (raw target alignments vs published aggregate
  consensus flanks) stay distinguishable in the reported support counts;
* the value checker refuses to pretend it can run controls without raw targets.

These tests deliberately do not assert specific scores, because scores depend on
model weights that are not present in a bare checkout. Inputs under the
gitignored ``data/`` tree are fetched separately, so tests needing them skip
cleanly rather than failing a fresh clone.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
for _path in (ROOT, ROOT / "scripts", ROOT / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import run_dual_evidence_case as case_mod  # noqa: E402
from pamdict.score.spacer import (  # noqa: E402
    find_spacer_hits,
    pam_from_hit,
    read_fasta,
    reverse_complement,
)

CASE_DIR = ROOT / "benchmarks" / "dual_evidence_cases"
PAIRED_CASE = CASE_DIR / "case_spcas9_pampredict_paired.json"
FLANK_CASE = CASE_DIR / "case_sp7f7_published_flanks.json"

PAMPREDICT_SPACERS = ROOT / "data/raw/external_datasets/PAMpredict/Example/spacers.fna"
PAMPREDICT_TARGETS = ROOT / "data/raw/external_datasets/PAMpredict/Example/Phages/phages.fna"
SP7F7_DIR = ROOT / "benchmarks/multi_evidence_v2/systems/sp7f7-published-flanks"

_HAS_PAMPREDICT = PAMPREDICT_SPACERS.exists() and PAMPREDICT_TARGETS.exists()
_HAS_SP7F7 = FLANK_CASE.exists() and (SP7F7_DIR / "published_flanks.tsv").exists()


# --------------------------------------------------------------------------
# case file contract
# --------------------------------------------------------------------------
def test_case_files_are_loadable_and_valid():
    cases = sorted(CASE_DIR.glob("case_*.json"))
    assert cases, "no case files found"
    for path in cases:
        case = case_mod.load_case(path)
        assert case["case_id"], path.name
        assert case["system_id"], path.name
        assert case["candidate_pams"], path.name
        widths = {len(pam) for pam in case["candidate_pams"]}
        assert widths == {int(case["pam_length"])}, path.name
        assert case["pam_side"] in {"downstream", "upstream"}, path.name


def test_case_loader_rejects_undocumented_orientation():
    payload = json.loads(PAIRED_CASE.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory() as tmp:
        bad = Path(tmp) / "bad.json"
        payload["spacer_orientation"] = {"decision": "maybe"}
        bad.write_text(json.dumps(payload), encoding="utf-8")
        try:
            case_mod.load_case(bad)
        except ValueError as exc:
            assert "must be 'reverse' or 'forward'" in str(exc), str(exc)
        else:
            raise AssertionError("undocumented orientation decision was accepted")

        payload["spacer_orientation"] = {"decision": "reverse"}
        bad.write_text(json.dumps(payload), encoding="utf-8")
        try:
            case_mod.load_case(bad)
        except ValueError as exc:
            assert "evidence" in str(exc), str(exc)
        else:
            raise AssertionError("orientation decision without evidence was accepted")


def test_case_loader_rejects_mixed_pam_widths():
    payload = json.loads(PAIRED_CASE.read_text(encoding="utf-8"))
    payload["candidate_pams"] = ["NGG", "NNNN"]
    payload["pam_length"] = 3
    with tempfile.TemporaryDirectory() as tmp:
        bad = Path(tmp) / "bad.json"
        bad.write_text(json.dumps(payload), encoding="utf-8")
        try:
            case_mod.load_case(bad)
        except ValueError:
            return
    raise AssertionError("mixed candidate widths were accepted")


def test_every_case_records_provenance_and_limitations():
    """A case without provenance or stated limits must not be shipped."""
    for path in sorted(CASE_DIR.glob("case_*.json")):
        case = json.loads(path.read_text(encoding="utf-8"))
        assert case.get("limitations"), f"{path.name}: no limitations recorded"
        assert case.get("claim_scope"), f"{path.name}: no claim scope recorded"
        assert case.get("evidence_tier"), f"{path.name}: no evidence tier recorded"
        assert "independent_benchmark_eligible" in case, path.name


# --------------------------------------------------------------------------
# orientation
# --------------------------------------------------------------------------
def test_orientation_probe_moves_the_signal_when_reversed():
    if not _HAS_PAMPREDICT:
        print("SKIP test_orientation_probe_moves_the_signal_when_reversed (inputs absent)")
        return
    spacers = read_fasta(PAMPREDICT_SPACERS)
    targets = read_fasta(PAMPREDICT_TARGETS)
    probe = case_mod.orientation_probe(spacers, targets, max_mismatches=4, pam_length=3)
    as_provided, reversed_ = probe["as_provided"], probe["reverse_complemented"]

    assert as_provided["stronger_side"] == "upstream"
    assert reversed_["stronger_side"] == "downstream"
    # Reversing moves the signal without creating or destroying it.
    assert abs(
        reversed_["downstream_information_sum"] - as_provided["upstream_information_sum"]
    ) < 1e-9
    assert as_provided["alignments"] == reversed_["alignments"] > 0


def test_hit_strand_is_not_orientation_evidence():
    """Both target strands are searched, so hit strand says nothing about array orientation."""
    if not _HAS_PAMPREDICT:
        print("SKIP test_hit_strand_is_not_orientation_evidence (inputs absent)")
        return
    spacers = [(n, reverse_complement(s)) for n, s in read_fasta(PAMPREDICT_SPACERS)]
    targets = read_fasta(PAMPREDICT_TARGETS)
    hits = find_spacer_hits(spacers, targets, flank_length=10, max_mismatches=4)
    assert {hit.strand for hit in hits} == {"+", "-"}
    pams = [pam for pam in (pam_from_hit(h, "downstream", 3) for h in hits) if pam]
    xgg = sum(1 for pam in pams if pam[1:] == "GG")
    assert xgg / len(pams) > 0.5, "signal should be coherent despite mixed hit strands"


def test_case_declares_orientation_matching_its_own_probe():
    if not _HAS_PAMPREDICT:
        print("SKIP test_case_declares_orientation_matching_its_own_probe (inputs absent)")
        return
    case = case_mod.load_case(PAIRED_CASE)
    spacers = read_fasta(ROOT / case["spacers"][0])
    targets = read_fasta(ROOT / case["targets"][0])
    probe = case_mod.orientation_probe(
        spacers, targets,
        max_mismatches=int(case["max_mismatches"]),
        pam_length=int(case["pam_length"]),
    )
    key = ("reverse_complemented"
           if case["spacer_orientation"]["decision"] == "reverse" else "as_provided")
    assert probe[key]["stronger_side"] == case["pam_side"]


# --------------------------------------------------------------------------
# evidence classes must stay distinguishable
# --------------------------------------------------------------------------
def test_published_flanks_are_distinguishable_from_raw_alignments():
    if not _HAS_SP7F7:
        print("SKIP test_published_flanks_are_distinguishable_from_raw_alignments")
        return
    case = case_mod.load_case(FLANK_CASE)
    spacers = read_fasta(ROOT / case["spacers"][0])
    hits, _excluded, diagnostics = case_mod.published_flank_hits(
        ROOT / case["published_flanks"], spacers, pam_side=case["pam_side"]
    )
    assert hits, "expected usable published flank records"
    assert diagnostics["orientation_sources"], "source orientation must be recorded"
    assert diagnostics["source_rows"] >= len(hits)

    rows = case_mod.runs_on_spacer_candidates(
        hits,
        candidates=list(case["candidate_pams"]),
        pam_side=case["pam_side"],
        prior_strength=2.0,
        background_records=[
            (h.target_id, h.upstream_flank + h.downstream_flank) for h in hits
        ],
        aggregate_records=True,
    )
    for row in rows:
        assert row["evidence_kind"] == "published_aggregate_consensus_flanks"
        # An aggregate consensus flank is not an independent target.
        assert row["support_target_count"] == 0
        assert row["support_alignment_count"] == 0
        assert row["eligible_target_count"] == 0


def test_raw_case_keeps_independent_support_counts():
    if not _HAS_PAMPREDICT:
        print("SKIP test_raw_case_keeps_independent_support_counts (inputs absent)")
        return
    case = case_mod.load_case(PAIRED_CASE)
    spacers = [(n, reverse_complement(s)) for n, s in read_fasta(ROOT / case["spacers"][0])]
    targets = read_fasta(ROOT / case["targets"][0])
    hits = find_spacer_hits(spacers, targets, flank_length=12, max_mismatches=4)
    rows = case_mod.runs_on_spacer_candidates(
        hits,
        candidates=["NGG", "NAG"],
        pam_side="downstream",
        prior_strength=2.0,
        background_records=targets,
        aggregate_records=False,
    )
    ngg = next(row for row in rows if row["candidate_pam"] == "NGG")
    assert ngg["evidence_kind"] == "raw_target_alignments"
    assert ngg["support_spacer_count"] > 0
    assert ngg["support_target_count"] > 0


def test_contig_collapse_exposes_redundancy():
    """One contig absorbing many spacers must be visible in the report."""
    if not _HAS_PAMPREDICT:
        print("SKIP test_contig_collapse_exposes_redundancy (inputs absent)")
        return
    case = case_mod.load_case(PAIRED_CASE)
    spacers = [(n, reverse_complement(s)) for n, s in read_fasta(ROOT / case["spacers"][0])]
    targets = read_fasta(ROOT / case["targets"][0])
    hits = find_spacer_hits(spacers, targets, flank_length=12, max_mismatches=4)
    rows = case_mod.hit_rows_from_hits(hits, pam_side="downstream", pam_length=3)
    collapse = case_mod.contig_collapse(rows)
    assert collapse["distinct_target_contigs_with_flank"] > 0
    assert collapse["max_spacers_on_one_contig"] >= 2


# --------------------------------------------------------------------------
# denominator correctness (regression for a real defect)
# --------------------------------------------------------------------------
def test_contig_fraction_divides_by_contigs_not_pam_types():
    """The contig-level ratio must divide by contigs, never by distinct PAM types.

    A shipped version divided the xGG contig count by the number of distinct
    consensus PAM strings. Several contigs can share one consensus, so that
    denominator is smaller and the ratio was inflated (7/15 reported as 7/11,
    i.e. 63.6% instead of 46.7%). This test constructs exactly that shape: more
    contigs than distinct consensus PAMs.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "verify_dual_evidence_value", ROOT / "scripts" / "verify_dual_evidence_value.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # 5 contigs, but only 2 distinct consensus PAMs: xGG shared by 4 contigs.
    best_per_spacer = {
        "s1": (0, "contigA", "TGG"),
        "s2": (0, "contigB", "TGG"),
        "s3": (0, "contigC", "TGG"),
        "s4": (0, "contigD", "TGG"),
        "s5": (0, "contigE", "AAT"),
    }
    stats = module.contig_level_xgg(best_per_spacer)

    assert stats["distinct_target_contigs"] == 5
    assert stats["distinct_contig_consensus_pam_types"] == 2
    assert stats["contigs_consensus_xgg"] == 4
    # The correct ratio is 4/5; the old bug would have produced 4/2 = 2.0.
    assert stats["contig_level_xgg_fraction"] == 4 / 5
    assert stats["contig_level_xgg_fraction"] <= 1.0


def test_contig_fraction_is_never_above_one():
    """A fraction above 1 is the signature of the wrong denominator."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "verify_dual_evidence_value", ROOT / "scripts" / "verify_dual_evidence_value.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # Every contig shares one consensus PAM: contigs=6, types=1.
    best_per_spacer = {
        f"s{i}": (0, f"contig{i}", "AGG") for i in range(6)
    }
    stats = module.contig_level_xgg(best_per_spacer)
    assert stats["distinct_target_contigs"] == 6
    assert stats["distinct_contig_consensus_pam_types"] == 1
    assert stats["contig_level_xgg_fraction"] == 1.0


def test_harness_contig_ratio_matches_its_own_parts():
    """The harness must report a ratio equal to its own numerator/denominator."""
    if not _HAS_PAMPREDICT:
        print("SKIP test_harness_contig_ratio_matches_its_own_parts (inputs absent)")
        return
    case = case_mod.load_case(PAIRED_CASE)
    spacers = [(n, reverse_complement(s)) for n, s in read_fasta(ROOT / case["spacers"][0])]
    targets = read_fasta(ROOT / case["targets"][0])
    hits = find_spacer_hits(spacers, targets, flank_length=12, max_mismatches=4)
    rows = case_mod.hit_rows_from_hits(hits, pam_side="downstream", pam_length=3)
    collapse = case_mod.contig_collapse(rows)

    total = collapse["distinct_target_contigs_with_flank"]
    xgg = collapse["contigs_supporting_xgg"]
    assert total > 0
    assert xgg <= total, "more xGG contigs than contigs is impossible"
    assert abs(collapse["contig_level_xgg_fraction"] - xgg / total) < 1e-12
    # The trap this guards: the PAM-type denominator is a different number here.
    assert collapse["distinct_contig_consensus_pam_types"] != total or total == 0


def test_recorded_run_artifact_uses_the_correct_denominator():
    """If the case was run locally, its artifact must not show the old defect."""
    path = ROOT / "data/parsed/dual_evidence_case_spcas9/case_summary.json"
    if not path.exists():
        print("SKIP test_recorded_run_artifact_uses_the_correct_denominator (no artifact)")
        return
    collapse = json.loads(path.read_text(encoding="utf-8"))["spacer_contig_collapse"]
    total = collapse["distinct_target_contigs_with_flank"]
    xgg = collapse["contigs_supporting_xgg"]
    assert abs(collapse["contig_level_xgg_fraction"] - xgg / total) < 1e-12
    assert collapse["contig_level_xgg_fraction"] <= 1.0



def test_value_checker_refuses_cases_without_raw_targets():
    """No raw targets means no homology controls: refuse rather than fabricate."""
    payload = json.loads(FLANK_CASE.read_text(encoding="utf-8"))
    assert not payload.get("targets"), "published-flank case should have no raw targets"
    source = (ROOT / "scripts" / "verify_dual_evidence_value.py").read_text(encoding="utf-8")
    assert "no raw target FASTA" in source


def test_agent_spec_covers_every_case_and_verifies_against_outputs():
    """The agent-facing spec must describe every shipped case.

    If a recorded outcome is present and the corresponding run artifact exists,
    the spec's claim is checked against it. This keeps the handoff document from
    drifting away from what the pipeline actually produces.
    """
    spec_path = CASE_DIR / "agent_case_spec.json"
    assert spec_path.exists(), "agent case spec is missing"
    spec = json.loads(spec_path.read_text(encoding="utf-8"))

    shipped = {case_mod.load_case(p)["case_id"] for p in CASE_DIR.glob("case_*.json")}
    described = {record["case_id"] for record in spec["structured_case_records"]}
    assert described == shipped, (
        f"spec covers {sorted(described)} but shipped cases are {sorted(shipped)}"
    )

    for record in spec["structured_case_records"]:
        assert record["traps_the_agent_must_avoid"], record["case_id"]
        assert record["honest_headline"], record["case_id"]

    # Cross-check recorded outcomes against real artifacts when they are present.
    outdirs = {
        "spcas9-pampredict-paired": ROOT / "data/parsed/dual_evidence_case_spcas9/case_summary.json",
        "sp7f7-published-flanks": ROOT / "data/parsed/dual_evidence_case_sp7f7/case_summary.json",
    }
    checked = 0
    for record in spec["structured_case_records"]:
        path = outdirs.get(record["case_id"])
        if path is None or not path.exists():
            continue
        summary = json.loads(path.read_text(encoding="utf-8"))
        fused = {row["candidate_pam"]: row for row in summary["fused"]}
        expected = record["expected_outcome_when_run"]

        if "support_spacer_count_ngg" in expected:
            assert fused["NGG"]["support_spacer_count"] == expected["support_spacer_count_ngg"]
            checked += 1
        if "support_target_count_ngg" in expected:
            assert fused["NGG"]["support_target_count"] == expected["support_target_count_ngg"]
            checked += 1
        if "evidence_status_ngg" in expected:
            assert fused["NGG"]["evidence_status"] == expected["evidence_status_ngg"]
            checked += 1
        if "spacer_top_candidate" in expected:
            top = max(
                fused.values(),
                key=lambda r: (r["spacer_score"] if r["spacer_score"] not in (None, "") else -1e9),
            )
            assert top["candidate_pam"] == expected["spacer_top_candidate"], record["case_id"]
            checked += 1

        # Contig-level claims must match the artifact and respect the contig
        # denominator; this is what keeps the spec's numbers honest.
        collapse = summary.get("spacer_contig_collapse", {})
        if "contigs_supporting_xgg" in expected:
            assert (
                collapse["contigs_supporting_xgg"] == expected["contigs_supporting_xgg"]
            ), record["case_id"]
            assert (
                collapse["distinct_target_contigs_with_flank"]
                == expected["distinct_target_contigs"]
            ), record["case_id"]
            checked += 2
        if "contig_level_xgg_fraction" in expected:
            # The spec records a rounded display value, so compare at the
            # precision the spec actually claims rather than bit-exactly.
            assert round(collapse["contig_level_xgg_fraction"], 3) == round(
                float(expected["contig_level_xgg_fraction"]), 3
            ), record["case_id"]
            assert collapse["contig_level_xgg_fraction"] <= 1.0, record["case_id"]
            # And the artifact itself must still equal its own numerator/denominator.
            assert abs(
                collapse["contig_level_xgg_fraction"]
                - collapse["contigs_supporting_xgg"]
                / collapse["distinct_target_contigs_with_flank"]
            ) < 1e-12, record["case_id"]
            checked += 1

    if checked == 0:
        print("SKIP agent-spec outcome cross-check (no run artifacts present)")


def test_agent_spec_declares_the_honest_value_statement():
    """Value must be stated as corroboration, not as an accuracy improvement."""
    spec = json.loads((CASE_DIR / "agent_case_spec.json").read_text(encoding="utf-8"))
    value = spec["value_statement"]
    assert value["what_it_does_not_add"], "must state what path 2 does not add"
    joined = " ".join(value["what_it_does_not_add"]).lower()
    assert "does not change the top-1 answer" in joined
    assert "strict independent paired systems remain zero" in joined

