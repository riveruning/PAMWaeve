"""Run the versioned protein + spacer + late-fusion benchmark."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.benchmark.contract import audit_manifest_document  # noqa: E402
from pamdict.benchmark.multi_evidence import (  # noqa: E402
    aggregate_benchmark,
    concrete_candidates,
    evaluate_fusion_rows,
    evaluate_scored_rows,
    parse_pam_spectrum,
)
from pamdict.score.spectrum import exact_match  # noqa: E402


def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else WORKSPACE / path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    """Return a repo-relative path so reports carry no absolute dev paths.

    Reports are artifacts that get copied into submission packages and read on
    other machines; embedding ``/home/<user>/...`` would both leak the build
    layout and make the record meaningless elsewhere.  Paths outside the repo
    are returned unchanged.
    """
    try:
        return path.resolve().relative_to(WORKSPACE.resolve()).as_posix()
    except ValueError:
        return str(path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def normalize_published_flank_rows(
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Prevent aggregate consensus records from masquerading as targets.

    A published aggregate flank is *one record per spacer*, not an observed
    target.  The target/alignment counters are therefore **unavailable**, not
    zero, and are blanked to the empty string.  Writing a literal ``0`` here
    would assert "zero targets supported this candidate", which is a different
    and false claim: the raw target identifiers were never published.  The
    underlying record counts are preserved in ``aggregate_*`` so nothing is
    lost.
    """
    normalized: list[dict[str, str]] = []
    for row in rows:
        value = dict(row)
        value.setdefault(
            "aggregate_support_record_count",
            str(value.get("support_target_count", "")),
        )
        value.setdefault(
            "aggregate_eligible_record_count",
            str(value.get("eligible_target_count", "")),
        )
        for field in (
            "support_target_count",
            "support_alignment_count",
            "eligible_target_count",
            "eligible_alignment_count",
        ):
            # Empty string = "not available from aggregate flank evidence".
            value[field] = ""
        normalized.append(value)
    return normalized


def fasta_record(path: Path, record_id: str) -> str:
    found = False
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if found:
                    break
                found = line[1:].split()[0] == record_id
            elif found:
                parts.append("".join(line.split()))
    sequence = "".join(parts).upper()
    if not sequence:
        raise ValueError(f"FASTA record {record_id!r} not found in {path}")
    return sequence


def write_single_fasta(path: Path, record_id: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f">{record_id}"] + [
        sequence[index:index + 80] for index in range(0, len(sequence), 80)
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def not_applicable_method(reason: str) -> dict[str, object]:
    """Mark a method as not applicable rather than as a zero-scoring failure.

    The metric fields are ``None``, never ``0.0``: a system that was never
    given a spacer channel has not "failed" the spacer method, and averaging
    it in as zero would silently understate every spacer/fusion summary.
    ``aggregate_benchmark`` excludes these rows from coverage denominators.
    """
    return {
        "not_applicable": True,
        "not_applicable_reason": reason,
        "n_gold_lengths": 0,
        "scored_lengths": 0,
        "abstained_lengths": 0,
        "coverage": 0.0,
        "mrr": None,
        "recall_at_1": None,
        "recall_at_3": None,
        "recall_at_5": None,
        "by_length": {},
    }


def run_command(command: list[str], environment: dict[str, str]) -> None:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(
        command,
        cwd=WORKSPACE,
        env=environment,
        check=True,
    )


def score_prior(
    prior: str,
    spectrum: dict[int, list[str]],
    candidates: list[str],
) -> dict[str, object]:
    return evaluate_scored_rows(
        [
            {
                "candidate_pam": candidate,
                "score": (
                    float(exact_match(candidate, prior))
                    if len(candidate) == len(prior) else None
                ),
            }
            for candidate in candidates
        ],
        spectrum,
        score_field="score",
    )


def _fmt(value: object, spec: str = ".3f") -> str:
    """Render a metric, or ``n/a`` when the method does not apply.

    Missing metrics are printed as ``n/a`` rather than ``0.000`` so that a
    channel which was never supplied cannot be misread as a zero score.
    """
    if value is None:
        return "n/a"
    return format(value, spec)


def render_markdown(report: dict[str, object]) -> str:
    aggregate = report["aggregate"]
    lines = [
        "# Multi-evidence PAM benchmark report",
        "",
        f"Generated: {report['generated_at']}",
        "",
        f"Decision: **{report['decision']}**",
        "",
        "## Dataset composition",
        "",
        f"- Manifest systems: {report['manifest_system_count']}",
        f"- Runnable systems evaluated: {aggregate['system_count']}",
        f"- Strict independent systems: {report['strict_independent_systems']} "
        f"{report.get('strict_independent_track_counts', {})}",
        f"- Strict paired-evidence systems: "
        f"{report.get('strict_paired_evidence_systems', 0)} "
        f"(minimum for a generalization claim: 5)",
        f"- Input tracks: {report['track_counts']}",
        "",
        "## Method summary",
        "",
        "| Method | Coverage | MRR (all) | Recall@1 | Recall@3 | Recall@5 | n/a | Applicable |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for method, values in aggregate["methods"].items():
        lines.append(
            "| {method} | {coverage} | {mrr} | {r1} | {r3} | {r5} | "
            "{na} | {total} |".format(
                method=method,
                coverage=_fmt(values["system_coverage"]),
                mrr=_fmt(values["mrr_all_systems"]),
                r1=_fmt(values["recall_at_1_all_systems"]),
                r3=_fmt(values["recall_at_3_all_systems"]),
                r5=_fmt(values["recall_at_5_all_systems"]),
                na=values.get("systems_not_applicable", 0),
                total=values.get("systems_total", 0),
            )
        )
    lines.extend([
        "",
        "## Evaluation-role split",
        "",
        "| Role | Systems | Protein MRR | Spacer MRR | Fusion MRR | Fusion coverage |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for role, values in report["by_evaluation_role"].items():
        methods = values["methods"]

        def role_metric(name: str, key: str = "mrr_all_systems") -> object:
            # A role can consist entirely of protein-only systems, in which
            # case the spacer/fusion methods are absent from the aggregate
            # rather than present as zero.  Report that as n/a.
            return methods.get(name, {}).get(key)

        lines.append(
            "| {role} | {systems} | {protein} | {spacer} | "
            "{fusion} | {coverage} |".format(
                role=role,
                systems=values["system_count"],
                protein=_fmt(role_metric("protein_only")),
                spacer=_fmt(role_metric("spacer_only")),
                fusion=_fmt(role_metric("fusion")),
                coverage=_fmt(role_metric("fusion", "system_coverage")),
            )
        )
    lines.extend([
        "",
        "## Per-system results",
        "",
        "| System | Role | Tier | Spacer evidence | Exposure | Mapped spacers | Orientation | Fusion gate | Protein MRR | Spacer MRR | Fusion MRR |",
        "|---|---|---|---|---|---:|---|---|---:|---:|---:|",
    ])
    for system in report["systems"]:
        methods = system["methods"]
        fusion_gates = sorted({
            gate
            for gates in system.get("fusion_diagnostics", {}).get(
                "gate_by_length", {}
            ).values()
            for gate in gates
        })
        lines.append(
            "| {system_id} | {role} | {tier} | {track} | {exposure} | {mapped} | {orientation} | "
            "{gate} | {protein} | {spacer} | {fusion} |".format(
                system_id=system["system_id"],
                role=system["evaluation_role"],
                tier=system["tier"],
                track=system["spacer_evidence_track"] or "protein_only",
                exposure=system["protein_training_exposure"],
                mapped=(
                    system["spacer_diagnostics"]["mapped_unique_spacers"]
                    if system["spacer_diagnostics"]["mapped_unique_spacers"]
                    is not None else "n/a"
                ),
                orientation=(
                    "n/a"
                    if system["spacer_diagnostics"]["orientation_pass"] is None
                    else "pass"
                    if system["spacer_diagnostics"]["orientation_pass"]
                    else "review"
                ),
                gate=", ".join(fusion_gates) or "n/a",
                protein=_fmt(methods.get("protein_only", {}).get("mrr")),
                spacer=_fmt(methods.get("spacer_only", {}).get("mrr")),
                fusion=_fmt(methods.get("fusion", {}).get("mrr")),
            )
        )
    lines.extend([
        "",
        "## Method applicability",
        "",
        "`spacer_only` and `fusion` are `n/a` for systems that supply no spacer",
        "evidence track (`protein_only`). Those systems are excluded from the",
        "spacer/fusion denominators and from paired fusion-minus-baseline",
        "comparisons, instead of being averaged in as zero scores.",
        "",
    ])
    lines.extend([
        "",
        "## Interpretation boundary",
        "",
        "Ranks are computed over exhaustive concrete A/C/G/T candidates within each PAM length. IUPAC gold motifs accept any concrete subset. Fusion uses an explicit evidence-state rank key, not a calibrated probability.",
        "",
    ])
    if report["strict_independent_systems"] == 0:
        lines.append(
            "No strict independent system is registered yet. These results validate the pipeline only and must not be reported as model generalization."
        )
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run a versioned system-level multi-evidence PAM benchmark."
    )
    parser.add_argument(
        "--manifest",
        default="benchmarks/multi_evidence_v2/manifest.json",
    )
    parser.add_argument(
        "--outdir",
        default="data/parsed/multi_evidence_benchmark_v2",
    )
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--system", action="append")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--max-concrete-candidates", type=int, default=100000)
    parser.add_argument(
        "--fusion-policy",
        choices=("abstain", "spacer_fallback", "legacy"),
        default="abstain",
        help="frozen disagreement policy used for fusion evaluation",
    )
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="reuse existing per-system inference files",
    )
    parser.add_argument(
        "--skip-inference",
        action="store_true",
        help="evaluate existing files only; fail if any are missing",
    )
    args = parser.parse_args()
    if args.bootstrap <= 0:
        raise ValueError("--bootstrap must be positive")

    manifest_path = resolve_path(args.manifest)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    systems = manifest.get("systems")
    if not isinstance(systems, list):
        raise ValueError("manifest systems must be a list")
    document_audit = audit_manifest_document(manifest)
    if not document_audit["valid"]:
        raise ValueError(f"manifest contract failed: {document_audit}")
    audit_by_id = {
        audit["system_id"]: audit
        for audit in document_audit["system_audits"]
    }
    if args.system:
        selected = set(args.system)
        systems = [row for row in systems if row.get("system_id") in selected]
        missing = selected - {row.get("system_id") for row in systems}
        if missing:
            raise ValueError(f"unknown --system values: {sorted(missing)}")

    outdir = resolve_path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    audits: list[dict[str, object]] = []
    runnable: list[dict[str, object]] = []
    for row in systems:
        audit = dict(audit_by_id[str(row["system_id"])])
        audit["errors"] = list(audit["errors"])
        audit["warnings"] = list(audit["warnings"])
        tracks = set(audit.get("available_tracks", []))
        path_fields = [("protein_fasta", row.get("protein_fasta"))]
        if "paired_evidence" in tracks:
            path_fields.append(("spacer_fasta", row.get("spacer_fasta")))
            path_fields.extend(
                ("target_fastas", value)
                for value in row.get("target_fastas", [])
            )
        if "published_flanks" in tracks:
            path_fields.extend([
                ("spacer_fasta", row.get("spacer_fasta")),
                ("published_flanks_tsv", row.get("published_flanks_tsv")),
            ])
        if "end_to_end_genome" in tracks:
            path_fields.append(
                ("host_genome_fasta", row.get("host_genome_fasta"))
            )
        for field, value in path_fields:
            if value and not resolve_path(str(value)).exists():
                audit["errors"].append(f"missing {field}: {value}")
        audit["valid"] = not audit["errors"]
        audit["runner_supported"] = bool(
            tracks
            & {"paired_evidence", "published_flanks", "protein_only"}
        )
        audit["runnable"] = (
            audit["valid"]
            and row.get("status") == "ready"
            and audit["runner_supported"]
        )
        if row.get("status") == "ready" and not audit["runner_supported"]:
            audit["warnings"].append(
                "ready system has no supported spacer-evidence track; skipped by this runner"
            )
        audits.append(audit)
        if audit["runnable"]:
            runnable.append(row)
    (outdir / "manifest_document_audit.json").write_text(
        json.dumps(document_audit, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (outdir / "manifest_audit.json").write_text(
        json.dumps(audits, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    invalid_ready = [
        audit for audit, row in zip(audits, systems)
        if row.get("status") == "ready" and not audit["valid"]
    ]
    if invalid_ready:
        raise ValueError(f"ready systems failed manifest audit: {invalid_ready}")
    if not runnable:
        raise ValueError("no runnable benchmark systems")

    environment = os.environ.copy()
    # Child processes must import ``pamdict`` from this tree, but optional
    # third-party trees (the Protein2PAM upstream adapter, vendored libraries)
    # are *deployment* choices, not fixed parts of this project.  So: always
    # prepend our own ``src``, only add the conventional local directories when
    # they actually exist, and otherwise preserve whatever the operator already
    # put on PYTHONPATH.  Clobbering it would make the package unusable outside
    # this developer checkout.
    inherited = [
        entry for entry in (environment.get("PYTHONPATH") or "").split(":")
        if entry
    ]
    local_candidates = [
        str(WORKSPACE / ".pylibs"),
        str(WORKSPACE / ".reference" / "Protein2PAM"),
    ]
    python_path: list[str] = []
    for entry in [
        str(WORKSPACE / "src"),
        *inherited,
        # Only add the conventional local directories when they really exist.
        # A stale entry would shadow the operator's own (working) install with
        # a non-existent path and break imports in the child process.
        *(path for path in local_candidates if Path(path).is_dir()),
    ]:
        if entry not in python_path:
            python_path.append(entry)
    environment["PYTHONPATH"] = ":".join(python_path)
    environment.setdefault(
        "HF_HOME", str(WORKSPACE / "data/checkpoints/hf")
    )
    evaluated_systems: list[dict[str, object]] = []

    def needs_run(available: bool, label: str) -> bool:
        if args.skip_inference:
            if not available:
                raise FileNotFoundError(f"missing inference artifact: {label}")
            return False
        return not (args.reuse_existing and available)

    for row in runnable:
        system_id = str(row["system_id"])
        row_tracks = set(row.get("available_tracks", []))
        # A system with neither ``paired_evidence`` nor ``published_flanks``
        # still carries a usable protein track (``protein_only``).  Such a
        # system is evaluated on the protein path alone: ``spacer_only`` and
        # ``fusion`` are reported as not applicable rather than as abstentions,
        # so a channel that was never supplied is never counted as a scoring
        # failure, and its absence is never averaged into a method comparison.
        spacer_evidence_track = (
            "paired_evidence" if "paired_evidence" in row_tracks
            else "published_flanks" if "published_flanks" in row_tracks
            else None
        )
        system_dir = outdir / system_id
        system_dir.mkdir(parents=True, exist_ok=True)
        spectrum = parse_pam_spectrum(row["gold_pam_spectrum"])
        candidates_by_length = {
            length: concrete_candidates(
                [length], max_candidates=args.max_concrete_candidates
            )
            for length in spectrum
        }
        all_candidates = [
            candidate
            for length in sorted(candidates_by_length)
            for candidate in candidates_by_length[length]
        ]
        all_candidate_path = system_dir / "candidate_pams.txt"
        all_candidate_path.write_text(
            "\n".join(all_candidates) + "\n", encoding="utf-8"
        )

        protein_source = resolve_path(str(row["protein_fasta"]))
        protein_sequence = fasta_record(
            protein_source, str(row["protein_record_id"])
        )
        protein_input = system_dir / "protein.fasta"
        write_single_fasta(protein_input, system_id, protein_sequence)
        protein_scores_path = system_dir / "protein_scores.tsv"
        protein_json_path = system_dir / "protein_scores.json"
        need_protein = needs_run(
            protein_scores_path.exists() and protein_json_path.exists(),
            f"{system_id} protein",
        )
        if need_protein and not args.skip_inference:
            command = [
                sys.executable,
                "scripts/score_candidate_pams.py",
                "--protein",
                str(protein_input),
                "--candidates-file",
                str(all_candidate_path),
                "--model",
                str(row["protein_model"]),
                "--device",
                args.device,
                "--out",
                str(protein_scores_path),
                "--json-out",
                str(protein_json_path),
            ]
            if row.get("protein_train_fasta") is not None:
                command.extend([
                    "--train-fasta", str(row["protein_train_fasta"])
                ])
            run_command(command, environment)

        combined_spacer_rows: list[dict[str, str]] = []
        spacer_summaries: dict[str, dict] = {}
        for length, candidates in candidates_by_length.items():
            length_dir = system_dir / f"spacer_len{length}"
            candidate_path = system_dir / f"candidate_pams_len{length}.txt"
            candidate_path.write_text(
                "\n".join(candidates) + "\n", encoding="utf-8"
            )
            if spacer_evidence_track is None:
                # protein_only: no spacer evidence channel exists for this
                # system, so nothing is computed and nothing is scored.
                continue
            score_path = length_dir / "candidate_scores.tsv"
            summary_path = length_dir / "summary.json"
            need_spacer = needs_run(
                score_path.exists() and summary_path.exists(),
                f"{system_id} spacer length {length}",
            )
            if need_spacer and not args.skip_inference:
                if spacer_evidence_track == "paired_evidence":
                    command = [
                        sys.executable, "scripts/score_spacer_pams.py",
                        "--spacers", str(resolve_path(str(row["spacer_fasta"]))),
                        "--system-id", system_id,
                        "--pam-side", str(row["pam_side"]),
                        "--pam-length", str(length),
                        "--max-mismatches", str(row.get("max_mismatches", 2)),
                        "--candidates-file", str(candidate_path),
                        "--outdir", str(length_dir),
                    ]
                    for target in row["target_fastas"]:
                        command.extend([
                            "--targets", str(resolve_path(str(target)))
                        ])
                    if row["spacer_orientation"] == "reverse":
                        command.append("--reverse-spacers")
                else:
                    command = [
                        sys.executable, "scripts/score_published_flanks.py",
                        "--spacers", str(resolve_path(str(row["spacer_fasta"]))),
                        "--published-flanks",
                        str(resolve_path(str(row["published_flanks_tsv"]))),
                        "--system-id", system_id,
                        "--pam-side", str(row["pam_side"]),
                        "--pam-length", str(length),
                        "--min-effective-spacers",
                        str(row["min_effective_spacers"]),
                        "--candidates-file", str(candidate_path),
                        "--outdir", str(length_dir),
                    ]
                run_command(command, environment)
            length_rows = read_tsv(score_path)
            if spacer_evidence_track == "published_flanks":
                length_rows = normalize_published_flank_rows(length_rows)
            combined_spacer_rows.extend(length_rows)
            spacer_summaries[str(length)] = json.loads(
                summary_path.read_text(encoding="utf-8")
            )
        combined_spacer_path = system_dir / "spacer_scores.tsv"
        write_tsv(combined_spacer_path, combined_spacer_rows)

        fusion_path = system_dir / "fusion_scores.tsv"
        # Fusion is a cheap deterministic CPU transform. Always rebuild it so
        # --reuse-existing cannot silently relabel a cached policy. When no
        # spacer channel exists there is nothing to fuse, so fusion is not run
        # and remains ``not_applicable`` instead of producing a ranked list.
        if spacer_evidence_track is not None:
            run_command([
                sys.executable,
                "scripts/fuse_pam_evidence.py",
                "--protein-scores",
                str(protein_scores_path),
                "--spacer-scores",
                str(combined_spacer_path),
                "--system-id",
                system_id,
                "--policy",
                args.fusion_policy,
                "--out",
                str(fusion_path),
            ], environment)

        protein_rows = read_tsv(protein_scores_path)
        fusion_rows = read_tsv(fusion_path) if fusion_path.exists() else []
        fusion_gate_by_length = {
            str(length): sorted({
                str(candidate_row.get("fusion_gate", "unannotated"))
                for candidate_row in fusion_rows
                if len(
                    str(candidate_row.get("candidate_pam", "")).strip()
                ) == length
            })
            for length in spectrum
        }
        fallback_lengths = sorted(
            length for length, gates in fusion_gate_by_length.items()
            if "spacer_fallback_severe_top_disagreement" in gates
        )
        abstain_conflict_lengths = sorted(
            length for length, gates in fusion_gate_by_length.items()
            if "abstain_severe_top_disagreement" in gates
        )
        methods = {
            "protein_only": evaluate_scored_rows(
                protein_rows,
                spectrum,
                score_field="specificity_adjusted_score",
            ),
            "spacer_only": (
                evaluate_scored_rows(
                    combined_spacer_rows,
                    spectrum,
                    score_field="spacer_evidence_score",
                )
                if spacer_evidence_track is not None
                # ``None`` (not a zero score) so that a channel which was
                # never supplied is excluded from method comparison rather
                # than averaged in as a failure.
                else not_applicable_method(
                    "system provides no spacer evidence track"
                )
            ),
            "fusion": (
                evaluate_fusion_rows(fusion_rows, spectrum)
                if spacer_evidence_track is not None
                else not_applicable_method(
                    "system provides no spacer evidence track, so there is "
                    "nothing to fuse"
                )
            ),
        }
        normalize_prior = str(row.get("family_prior_pam", "")).strip().upper()
        if normalize_prior:
            methods["family_prior"] = score_prior(
                normalize_prior,
                spectrum,
                all_candidates,
            )

        orientation_pass = (
            all(
                summary["inferred_signal_side"] == row["pam_side"]
                and not summary["orientation_warning"]
                for summary in spacer_summaries.values()
            )
            if spacer_summaries
            else None
        )
        mapped_unique = (
            min(
                int(summary["counts"]["mapped_unique_spacers"])
                for summary in spacer_summaries.values()
            )
            if spacer_summaries
            else None
        )
        training_exact = (
            str(protein_rows[0].get("training_exact_match", "")).lower() == "true"
            if protein_rows else False
        )
        warnings = []
        if row["protein_training_exposure"] == "none" and training_exact:
            warnings.append(
                "manifest says no protein exposure but runtime found an exact training match"
            )
        if orientation_pass is False:
            warnings.append("spacer-side orientation signal requires review")
        if fallback_lengths:
            warnings.append(
                "exploratory fusion fell back to spacer-only ranking for "
                f"PAM lengths {fallback_lengths} because source tops disagreed"
            )
        if abstain_conflict_lengths:
            warnings.append(
                "fusion abstained for PAM lengths "
                f"{abstain_conflict_lengths} because source tops disagreed"
            )
        evaluated = {
            "system_id": system_id,
            "tier": row["tier"],
            "evaluation_role": row.get("evaluation_role", "unspecified"),
            "available_tracks": sorted(row.get("available_tracks", [])),
            "spacer_evidence_track": spacer_evidence_track,
            "cas_family": row["cas_family"],
            "gold_pam_spectrum": {
                str(length): values for length, values in spectrum.items()
            },
            "protein_training_exposure": row["protein_training_exposure"],
            # Strict independence is audited *per track*.  A protein-only
            # system must not inherit "not eligible" from the absence of a
            # paired-evidence cohort, and a paired system must not inherit
            # eligibility from an independent protein alone.
            "strict_independent_eligible": bool(
                next(
                    (
                        audit["strict_independent_eligible_tracks"]
                        for audit in audits
                        if audit["system_id"] == system_id
                    ),
                    [],
                )
            ),
            "strict_independent_tracks": sorted(
                next(
                    (
                        audit["strict_independent_eligible_tracks"]
                        for audit in audits
                        if audit["system_id"] == system_id
                    ),
                    [],
                )
            ),
            "candidate_definition": "exhaustive concrete A/C/G/T strings per gold PAM length",
            "candidate_count": len(all_candidates),
            "protein_diagnostics": {
                "runtime_training_exact_match": training_exact,
                "nearest_train_identity": (
                    protein_rows[0].get("nearest_train_identity")
                    if protein_rows else None
                ),
            },
            "spacer_diagnostics": {
                "mapped_unique_spacers": mapped_unique,
                "orientation_pass": orientation_pass,
                "by_length": {
                    length: {
                        "counts": summary["counts"],
                        "inferred_signal_side": summary["inferred_signal_side"],
                        "orientation_warning": summary["orientation_warning"],
                    }
                    for length, summary in spacer_summaries.items()
                },
            },
            "fusion_diagnostics": {
                "policy": args.fusion_policy,
                "gate_by_length": fusion_gate_by_length,
                "fallback_lengths": fallback_lengths,
                "abstain_conflict_lengths": abstain_conflict_lengths,
            },
            "methods": methods,
            "warnings": warnings,
            "artifacts": {
                "protein_scores": rel(protein_scores_path),
                "spacer_scores": (
                    rel(combined_spacer_path) if spacer_evidence_track else None
                ),
                "fusion_scores": (
                    rel(fusion_path) if spacer_evidence_track else None
                ),
            },
        }
        (system_dir / "evaluation.json").write_text(
            json.dumps(evaluated, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        evaluated_systems.append(evaluated)

    aggregate = aggregate_benchmark(
        evaluated_systems,
        bootstrap_iterations=args.bootstrap,
    )
    by_tier = {
        tier: aggregate_benchmark(
            [system for system in evaluated_systems if system["tier"] == tier],
            bootstrap_iterations=args.bootstrap,
        )
        for tier in sorted({system["tier"] for system in evaluated_systems})
    }
    by_evaluation_role = {
        role: aggregate_benchmark(
            [
                system for system in evaluated_systems
                if system["evaluation_role"] == role
            ],
            bootstrap_iterations=args.bootstrap,
        )
        for role in sorted({
            system["evaluation_role"] for system in evaluated_systems
        })
    }
    strict_systems = sum(
        bool(system["strict_independent_eligible"])
        for system in evaluated_systems
    )
    # The "results available" bar is only met by a strict cohort on the tracks
    # this runner can actually evaluate end to end.  Strict protein-only
    # systems are counted and reported separately: they support a protein-track
    # result, but three of them are not a strict *paired-evidence* cohort, and
    # two of those three share one paper, so they are not three independent
    # publications either.
    strict_track_counts: dict[str, int] = {}
    for system in evaluated_systems:
        for track in system.get("strict_independent_tracks", []):
            strict_track_counts[track] = strict_track_counts.get(track, 0) + 1
    strict_paired = strict_track_counts.get("paired_evidence", 0)
    decision = (
        "BENCHMARK_RESULTS_AVAILABLE"
        if strict_paired >= 5
        else "FRAMEWORK_VALIDATED_NO_STRICT_COHORT"
    )
    report = {
        "format_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "manifest": rel(manifest_path),
        "manifest_sha256": file_sha256(manifest_path),
        "manifest_format_version": manifest.get("format_version"),
        "manifest_system_count": len(systems),
        "strict_independent_systems": strict_systems,
        "strict_independent_track_counts": strict_track_counts,
        "strict_paired_evidence_systems": strict_paired,
        "track_counts": document_audit["track_counts"],
        "candidate_policy": (
            "exhaustive concrete A/C/G/T candidates, ranked within PAM length"
        ),
        "fusion_policy": args.fusion_policy,
        "metrics": ["MRR", "Recall@1", "Recall@3", "Recall@5", "coverage"],
        "aggregate": aggregate,
        "by_tier": by_tier,
        "by_evaluation_role": by_evaluation_role,
        "systems": evaluated_systems,
        "decision": decision,
        "reporting_guard": (
            "Do not claim model generalization without a strict independent cohort."
        ),
    }
    report_path = outdir / "benchmark_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (outdir / "benchmark_report.md").write_text(
        render_markdown(report), encoding="utf-8"
    )
    summary_rows = []
    for system in evaluated_systems:
        base = {
            "system_id": system["system_id"],
            "tier": system["tier"],
            "evaluation_role": system["evaluation_role"],
            "spacer_evidence_track": system["spacer_evidence_track"],
            "protein_training_exposure": system["protein_training_exposure"],
            "strict_independent_eligible": system["strict_independent_eligible"],
            "candidate_count": system["candidate_count"],
            "mapped_unique_spacers": system["spacer_diagnostics"]["mapped_unique_spacers"],
            "orientation_pass": system["spacer_diagnostics"]["orientation_pass"],
        }
        for method, values in system["methods"].items():
            base[f"{method}_coverage"] = values["coverage"]
            # Empty cell, never 0, when the method does not apply: the TSV
            # must not turn "no evidence channel" into a zero score.
            if values.get("not_applicable"):
                base[f"{method}_mrr"] = ""
                base[f"{method}_recall_at_1"] = ""
                base[f"{method}_recall_at_3"] = ""
                base[f"{method}_recall_at_5"] = ""
                continue
            base[f"{method}_mrr"] = values["mrr"]
            base[f"{method}_recall_at_1"] = values["recall_at_1"]
            base[f"{method}_recall_at_3"] = values["recall_at_3"]
            base[f"{method}_recall_at_5"] = values["recall_at_5"]
        summary_rows.append(base)
    write_tsv(outdir / "benchmark_systems.tsv", summary_rows)
    print(json.dumps({
        "decision": decision,
        "systems": len(evaluated_systems),
        "strict_independent_systems": strict_systems,
        "report": rel(report_path),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
