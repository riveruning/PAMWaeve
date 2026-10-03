"""Versioned benchmark manifest contract and independence audit."""
from __future__ import annotations

import hashlib
import re
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Mapping

from pamdict.benchmark.multi_evidence import (
    audit_manifest_row as audit_v1_manifest_row,
    parse_pam_spectrum,
)

SUPPORTED_FORMAT_VERSIONS = {1, 2}
TIERS = {
    "engineering_regression",
    "retrospective_biological",
    "strict_independent",
}
STATUSES = {"ready", "curating", "blocked", "deprecated"}
TRACKS = {
    "protein_only",
    "paired_evidence",
    "published_flanks",
    "end_to_end_genome",
}
EXPOSURES = {"none", "near", "exact", "unknown"}
TARGET_INDEPENDENCE = {
    "virus_clustered",
    "contig_only",
    "not_applicable",
    "not_assessable_from_aggregate_flanks",
    "unknown",
}
GOLD_SYSTEM_MATCHES = {
    "exact_protein_or_strain",
    "species_ortholog",
    "not_resolved",
}
EXPERIMENTAL_GOLD_TYPES = {
    "experimental_activity",
    "experimental_functional_spectrum",
}
V2_REQUIRED_FIELDS = (
    "system_id",
    "tier",
    "status",
    "crispr_type",
    "cas_family",
    "available_tracks",
    "protein_model",
    "protein_fasta",
    "protein_record_id",
    "pam_side",
    "pam_lengths",
    "gold_pam_spectrum",
    "gold_evidence_type",
    "gold_doi",
    "gold_assay",
    "gold_system_match",
    "protein_training_exposure",
    "target_independence",
    "host_taxon",
    "host_strain",
    "split_group",
    "provenance",
    "notes",
)
PROVENANCE_FIELDS = (
    "source_name",
    "source_uri",
    "source_version",
    "license",
    "sha256",
    "sha256_scope",
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _nonempty(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _safe_relative_path(value: object) -> bool:
    if not _nonempty(value):
        return False
    path = PurePosixPath(str(value))
    return not path.is_absolute() and ".." not in path.parts and "\\" not in str(value)


def _audit_provenance_entry(
    entry: object,
    *,
    label: str,
    errors: list[str],
    warnings: list[str],
) -> bool:
    if not isinstance(entry, Mapping):
        errors.append(f"provenance.{label} must be an object")
        return False
    missing = [field for field in PROVENANCE_FIELDS if not _nonempty(entry.get(field))]
    if missing:
        errors.append(f"provenance.{label} missing values: {missing}")
    digest = str(entry.get("sha256", "")).lower()
    if digest and not SHA256_RE.fullmatch(digest):
        errors.append(f"provenance.{label}.sha256 must be 64 lowercase hex characters")
    if str(entry.get("license", "")).lower() in {"", "unknown", "not_recorded"}:
        warnings.append(f"provenance.{label} license is unresolved")
    return not missing and bool(SHA256_RE.fullmatch(digest))


def audit_v2_manifest_row(row: dict[str, object]) -> dict[str, object]:
    """Audit one v2 system without assuming every track has spacer evidence."""
    errors: list[str] = []
    warnings: list[str] = []
    missing = [field for field in V2_REQUIRED_FIELDS if field not in row]
    if missing:
        errors.append(f"missing v2 fields: {missing}")

    system_id = row.get("system_id", "")
    if not _nonempty(system_id):
        errors.append("system_id must be a non-empty string")
    tier = row.get("tier")
    if tier not in TIERS:
        errors.append(f"tier must be one of {sorted(TIERS)}")
    status = row.get("status")
    if status not in STATUSES:
        errors.append(f"status must be one of {sorted(STATUSES)}")
    if row.get("pam_side") not in ("upstream", "downstream"):
        errors.append("pam_side must be upstream or downstream")
    exposure = row.get("protein_training_exposure")
    if exposure not in EXPOSURES:
        errors.append(
            f"protein_training_exposure must be one of {sorted(EXPOSURES)}"
        )
    target_independence = row.get("target_independence")
    if target_independence not in TARGET_INDEPENDENCE:
        errors.append(
            f"target_independence must be one of {sorted(TARGET_INDEPENDENCE)}"
        )
    if row.get("gold_system_match") not in GOLD_SYSTEM_MATCHES:
        errors.append(
            f"gold_system_match must be one of {sorted(GOLD_SYSTEM_MATCHES)}"
        )

    raw_tracks = row.get("available_tracks")
    tracks = set(raw_tracks) if isinstance(raw_tracks, list) else set()
    if not tracks or tracks - TRACKS:
        errors.append(f"available_tracks must be a non-empty subset of {sorted(TRACKS)}")
    if len(tracks) != len(raw_tracks or []):
        errors.append("available_tracks must not contain duplicates")

    try:
        spectrum = parse_pam_spectrum(row.get("gold_pam_spectrum"))
    except (TypeError, ValueError) as exc:
        spectrum = {}
        errors.append(str(exc))
    lengths = row.get("pam_lengths")
    if not isinstance(lengths, list) or sorted(set(lengths)) != sorted(spectrum):
        errors.append("pam_lengths must equal gold_pam_spectrum length keys")

    for field in (
        "protein_model",
        "protein_fasta",
        "protein_record_id",
        "host_taxon",
        "host_strain",
        "split_group",
    ):
        if not _nonempty(row.get(field)):
            errors.append(f"{field} must be a non-empty string")
    if row.get("protein_fasta") and not _safe_relative_path(
        row.get("protein_fasta")
    ):
        errors.append("protein_fasta must be a workspace-relative safe path")

    if "paired_evidence" in tracks:
        for field in ("spacer_fasta", "spacer_orientation", "target_fastas"):
            if field not in row:
                errors.append(f"paired_evidence requires {field}")
        if row.get("spacer_orientation") not in ("forward", "reverse"):
            errors.append("spacer_orientation must be forward or reverse")
        if not isinstance(row.get("target_fastas"), list) or not row.get(
            "target_fastas"
        ):
            errors.append("paired_evidence requires non-empty target_fastas")
        if row.get("spacer_fasta") and not _safe_relative_path(
            row.get("spacer_fasta")
        ):
            errors.append("spacer_fasta must be a workspace-relative safe path")
        for target_path in row.get("target_fastas", []):
            if not _safe_relative_path(target_path):
                errors.append(
                    "target_fastas entries must be workspace-relative safe paths"
                )
        mismatches = row.get("max_mismatches")
        if not isinstance(mismatches, int) or mismatches < 0:
            errors.append("paired_evidence requires non-negative max_mismatches")
    if "published_flanks" in tracks:
        for field in ("spacer_fasta", "published_flanks_tsv"):
            if not _nonempty(row.get(field)):
                errors.append(f"published_flanks requires {field}")
            elif not _safe_relative_path(row.get(field)):
                errors.append(f"{field} must be a workspace-relative safe path")
        minimum = row.get("min_effective_spacers")
        if not isinstance(minimum, int) or minimum < 1:
            errors.append(
                "published_flanks requires positive min_effective_spacers"
            )
    if "end_to_end_genome" in tracks:
        if not _nonempty(row.get("host_genome_fasta")):
            errors.append("end_to_end_genome requires host_genome_fasta")
        elif not _safe_relative_path(row.get("host_genome_fasta")):
            errors.append(
                "host_genome_fasta must be a workspace-relative safe path"
            )
        if not isinstance(row.get("discovery_policy"), Mapping):
            errors.append("end_to_end_genome requires discovery_policy object")

    provenance = row.get("provenance")
    provenance_valid: dict[str, bool] = {}
    if not isinstance(provenance, Mapping):
        errors.append("provenance must be an object")
        provenance = {}
    provenance_valid["protein"] = _audit_provenance_entry(
        provenance.get("protein"),
        label="protein",
        errors=errors,
        warnings=warnings,
    )
    if "paired_evidence" in tracks:
        provenance_valid["spacers"] = _audit_provenance_entry(
            provenance.get("spacers"),
            label="spacers",
            errors=errors,
            warnings=warnings,
        )
        raw_targets = provenance.get("targets")
        if not isinstance(raw_targets, list) or not raw_targets:
            errors.append("provenance.targets must be a non-empty list")
            provenance_valid["targets"] = False
        else:
            target_checks = [
                _audit_provenance_entry(
                    entry,
                    label=f"targets[{index}]",
                    errors=errors,
                    warnings=warnings,
                )
                for index, entry in enumerate(raw_targets)
            ]
            provenance_valid["targets"] = all(target_checks)
            if len(raw_targets) != len(row.get("target_fastas", [])):
                errors.append(
                    "provenance.targets must align one-to-one with target_fastas"
                )
    if "published_flanks" in tracks:
        provenance_valid["spacers"] = _audit_provenance_entry(
            provenance.get("spacers"),
            label="spacers",
            errors=errors,
            warnings=warnings,
        )
        provenance_valid["published_flanks"] = _audit_provenance_entry(
            provenance.get("published_flanks"),
            label="published_flanks",
            errors=errors,
            warnings=warnings,
        )
    if "end_to_end_genome" in tracks:
        provenance_valid["host_genome"] = _audit_provenance_entry(
            provenance.get("host_genome"),
            label="host_genome",
            errors=errors,
            warnings=warnings,
        )

    gold_assay = row.get("gold_assay")
    if not isinstance(gold_assay, Mapping):
        errors.append("gold_assay must be an object")
        gold_assay = {}
    for field in (
        "assay_type",
        "activity_values_available",
        "conditions",
        "source_locator",
    ):
        if field not in gold_assay:
            errors.append(f"gold_assay missing field: {field}")
    if not isinstance(gold_assay.get("activity_values_available"), bool):
        errors.append("gold_assay.activity_values_available must be boolean")

    target_dataset = row.get("target_dataset")
    if "paired_evidence" in tracks:
        if not isinstance(target_dataset, Mapping):
            errors.append("paired_evidence requires target_dataset object")
            target_dataset = {}
        for field in (
            "name",
            "release",
            "frozen_at",
            "selection_policy",
            "cluster_policy",
        ):
            if not _nonempty(target_dataset.get(field)):
                errors.append(f"target_dataset.{field} must be non-empty")
    if "published_flanks" in tracks:
        flank_dataset = row.get("published_flank_dataset")
        if not isinstance(flank_dataset, Mapping):
            errors.append("published_flanks requires published_flank_dataset object")
            flank_dataset = {}
        for field in (
            "name",
            "doi",
            "release",
            "frozen_at",
            "row_selection_policy",
            "orientation_policy",
        ):
            if not _nonempty(flank_dataset.get(field)):
                errors.append(f"published_flank_dataset.{field} must be non-empty")
        if not isinstance(
            flank_dataset.get("raw_target_identifiers_available"), bool
        ):
            errors.append(
                "published_flank_dataset.raw_target_identifiers_available must be boolean"
            )

    base_strict = (
        tier == "strict_independent"
        and status == "ready"
        and exposure == "none"
        and _nonempty(row.get("gold_doi"))
        and row.get("gold_evidence_type") in EXPERIMENTAL_GOLD_TYPES
        and row.get("gold_system_match") == "exact_protein_or_strain"
        and provenance_valid.get("protein", False)
        and str(gold_assay.get("assay_type", "")) not in {
            "",
            "known_reference",
            "not_recorded",
        }
    )
    eligible_tracks: list[str] = []
    if base_strict and "protein_only" in tracks:
        eligible_tracks.append("protein_only")
    if (
        base_strict
        and "paired_evidence" in tracks
        and target_independence == "virus_clustered"
        and provenance_valid.get("spacers", False)
        and provenance_valid.get("targets", False)
        and str(target_dataset.get("cluster_policy", "")) not in {
            "",
            "not_clustered",
            "not_recorded",
        }
    ):
        eligible_tracks.append("paired_evidence")
    if (
        base_strict
        and "end_to_end_genome" in tracks
        and provenance_valid.get("host_genome", False)
    ):
        eligible_tracks.append("end_to_end_genome")

    if tier == "strict_independent" and set(eligible_tracks) != tracks:
        errors.append(
            "strict_independent row does not satisfy every declared track contract"
        )
    if exposure in {"exact", "near"}:
        warnings.append("protein is exposed to Protein2PAM training distribution")
    if "paired_evidence" in tracks and target_independence != "virus_clustered":
        warnings.append("target records are not verified independent virus clusters")
    if "published_flanks" in tracks:
        warnings.append(
            "published flanks are aggregate per-spacer evidence without raw target identifiers"
        )
    if status != "ready":
        warnings.append("system is not runnable")

    return {
        "system_id": system_id,
        "format_version": 2,
        "available_tracks": sorted(tracks),
        "valid": not errors,
        "runnable": not errors and status == "ready",
        "strict_independent_eligible": bool(eligible_tracks),
        "strict_independent_eligible_tracks": eligible_tracks,
        "errors": errors,
        "warnings": list(dict.fromkeys(warnings)),
    }


def audit_manifest_document(manifest: object) -> dict[str, object]:
    """Audit a complete versioned manifest, including duplicate-system guards."""
    document_errors: list[str] = []
    document_warnings: list[str] = []
    if not isinstance(manifest, Mapping):
        return {
            "valid": False,
            "format_version": None,
            "errors": ["manifest must be an object"],
            "warnings": [],
            "system_audits": [],
        }

    version = manifest.get("format_version")
    if version not in SUPPORTED_FORMAT_VERSIONS:
        document_errors.append(
            f"format_version must be one of {sorted(SUPPORTED_FORMAT_VERSIONS)}"
        )
    if not _nonempty(manifest.get("benchmark_name")):
        document_errors.append("benchmark_name must be a non-empty string")
    systems = manifest.get("systems")
    if not isinstance(systems, list):
        document_errors.append("systems must be a list")
        systems = []

    if version == 2:
        for field in ("frozen_at", "candidate_policy", "minimum_strict_systems"):
            if field not in manifest:
                document_errors.append(f"v2 manifest missing field: {field}")
        minimum = manifest.get("minimum_strict_systems")
        if not isinstance(minimum, int) or minimum < 1:
            document_errors.append("minimum_strict_systems must be a positive integer")

    ids = [row.get("system_id") for row in systems if isinstance(row, Mapping)]
    duplicates = sorted(key for key, count in Counter(ids).items() if count > 1)
    if duplicates:
        document_errors.append(f"duplicate system_id values: {duplicates}")

    audits = []
    for index, row in enumerate(systems):
        if not isinstance(row, dict):
            audits.append({
                "system_id": "",
                "valid": False,
                "runnable": False,
                "strict_independent_eligible": False,
                "strict_independent_eligible_tracks": [],
                "errors": [f"systems[{index}] must be an object"],
                "warnings": [],
            })
        elif version == 2:
            audits.append(audit_v2_manifest_row(row))
        else:
            legacy = audit_v1_manifest_row(row)
            legacy["format_version"] = 1
            legacy["available_tracks"] = ["protein_only", "paired_evidence"]
            legacy["strict_independent_eligible_tracks"] = (
                ["protein_only", "paired_evidence"]
                if legacy["strict_independent_eligible"] else []
            )
            audits.append(legacy)

    track_counts = Counter(
        track for audit in audits for track in audit.get("available_tracks", [])
    )
    strict_track_counts = Counter(
        track
        for audit in audits
        for track in audit.get("strict_independent_eligible_tracks", [])
    )
    minimum = manifest.get("minimum_strict_systems", 5)
    if version == 2 and strict_track_counts.get("paired_evidence", 0) < minimum:
        document_warnings.append(
            "strict paired-evidence cohort is below minimum_strict_systems"
        )

    return {
        "valid": not document_errors and all(audit["valid"] for audit in audits),
        "format_version": version,
        "benchmark_name": manifest.get("benchmark_name", ""),
        "system_count": len(systems),
        "track_counts": dict(sorted(track_counts.items())),
        "strict_independent_track_counts": dict(sorted(strict_track_counts.items())),
        "errors": document_errors,
        "warnings": document_warnings,
        "system_audits": audits,
    }

def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fasta_record_sequence(path: Path, record_id: str) -> str:
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
        raise ValueError(f"FASTA record {record_id!r} not found")
    return sequence


def verify_manifest_files(
    manifest: Mapping[str, object],
    workspace: Path,
) -> dict[str, object]:
    """Verify frozen local inputs against the v2 provenance checksums."""
    workspace = workspace.resolve()
    checks: list[dict[str, object]] = []

    def verify_file(
        *,
        system_id: str,
        role: str,
        relative_path: object,
        expected: object,
        record_id: str | None = None,
    ) -> None:
        result: dict[str, object] = {
            "system_id": system_id,
            "role": role,
            "path": str(relative_path or ""),
            "valid": False,
            "expected_sha256": str(expected or "").lower(),
        }
        if not _safe_relative_path(relative_path):
            result["error"] = "unsafe or missing relative path"
            checks.append(result)
            return
        path = (workspace / str(relative_path)).resolve()
        try:
            path.relative_to(workspace)
        except ValueError:
            result["error"] = "resolved path escapes workspace"
            checks.append(result)
            return
        if not path.is_file():
            result["error"] = "file not found"
            checks.append(result)
            return
        try:
            if record_id is None:
                observed = _sha256_file(path)
            else:
                sequence = _fasta_record_sequence(path, record_id)
                observed = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
        except (OSError, UnicodeError, ValueError) as exc:
            result["error"] = str(exc)
            checks.append(result)
            return
        result["observed_sha256"] = observed
        result["valid"] = observed == result["expected_sha256"]
        if not result["valid"]:
            result["error"] = "checksum mismatch"
        checks.append(result)

    systems = manifest.get("systems", [])
    if not isinstance(systems, list):
        systems = []
    for row in systems:
        if not isinstance(row, Mapping):
            continue
        system_id = str(row.get("system_id", ""))
        tracks = set(row.get("available_tracks", []))
        provenance = row.get("provenance", {})
        if not isinstance(provenance, Mapping):
            provenance = {}
        protein = provenance.get("protein", {})
        if not isinstance(protein, Mapping):
            protein = {}
        verify_file(
            system_id=system_id,
            role="protein",
            relative_path=row.get("protein_fasta"),
            expected=protein.get("sha256"),
            record_id=str(row.get("protein_record_id", "")),
        )
        if "paired_evidence" in tracks:
            spacers = provenance.get("spacers", {})
            if not isinstance(spacers, Mapping):
                spacers = {}
            verify_file(
                system_id=system_id,
                role="spacers",
                relative_path=row.get("spacer_fasta"),
                expected=spacers.get("sha256"),
            )
            targets = row.get("target_fastas", [])
            target_provenance = provenance.get("targets", [])
            for index, target_path in enumerate(
                targets if isinstance(targets, list) else []
            ):
                entry = (
                    target_provenance[index]
                    if isinstance(target_provenance, list)
                    and index < len(target_provenance)
                    and isinstance(target_provenance[index], Mapping)
                    else {}
                )
                verify_file(
                    system_id=system_id,
                    role=f"targets[{index}]",
                    relative_path=target_path,
                    expected=entry.get("sha256"),
                )
        if "published_flanks" in tracks:
            spacers = provenance.get("spacers", {})
            if not isinstance(spacers, Mapping):
                spacers = {}
            verify_file(
                system_id=system_id,
                role="spacers",
                relative_path=row.get("spacer_fasta"),
                expected=spacers.get("sha256"),
            )
            published_flanks = provenance.get("published_flanks", {})
            if not isinstance(published_flanks, Mapping):
                published_flanks = {}
            verify_file(
                system_id=system_id,
                role="published_flanks",
                relative_path=row.get("published_flanks_tsv"),
                expected=published_flanks.get("sha256"),
            )
        if "end_to_end_genome" in tracks:
            genome = provenance.get("host_genome", {})
            if not isinstance(genome, Mapping):
                genome = {}
            verify_file(
                system_id=system_id,
                role="host_genome",
                relative_path=row.get("host_genome_fasta"),
                expected=genome.get("sha256"),
            )
    return {
        "valid": all(check["valid"] for check in checks),
        "checked_artifacts": len(checks),
        "failed_artifacts": sum(not check["valid"] for check in checks),
        "checks": checks,
    }
