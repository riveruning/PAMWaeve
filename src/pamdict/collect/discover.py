"""Discovery + conservative classification of candidate PAM papers.

Mirrors the reference project's discipline: automated labels are triage
signals only, never experimental ground truth. Deduplicates by DOI (then
PMID/PMCID, then normalized title). Classifies conservatively into
candidate_experimental / review_or_secondary / computational_only / uncertain.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

REVIEW_MARKERS = [
    "review", "systematic review", "meta-analysis", "perspective", "protocol",
    "survey of", "inspection of", "primer", "mini review", "minireview",
]
COMPUTATIONAL_MARKERS = [
    "computational", "in silico", "bioinformatic", "machine learning",
    "deep learning", "prediction of", "predicted", "algorithm",
    "web server", "database", "tool for", "benchmark",
    "neural network", "model predicts", "trained model",
]
ASSAY_TERMS = [
    "pam depletion", "pam screen", "pam screening", "interference assay",
    "interference assays", "transformation assay", "transformation assays",
    "plasmid clearance", "plasmid interference", "plasmid loss",
    "high-throughput pam", "depletion screen",
    # experimental PAM determination (endogenous/native systems)
    "pam identification", "pam characterization", "pam characterisation",
    "pam determination", "protospacer adjacent motif",
    "protospacer-adjacent motif", "pam recognition", "pam library",
    "cleavage assay", "pam-scanr", "pam-detect", "ht-pamda", "pamda",
    "functional pam", "pam specificity", "pam preference",
    "endogenous crispr", "native crispr",
]


def _norm_title(title: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", title.lower())


@dataclass
class Candidate:
    record_id: str
    title: str
    authors: str
    year: str
    journal: str
    doi: str
    pmid: str
    pmcid: str
    abstract: str
    source_api: str
    matched_queries: list[str] = field(default_factory=list)
    inferred_crispr_types: list[str] = field(default_factory=list)
    assay_keyword_hits: list[str] = field(default_factory=list)
    evidence_status: str = "uncertain"
    evidence_reason: str = ""
    full_text_available: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "record_id": self.record_id,
            "title": self.title,
            "authors": self.authors,
            "year": self.year,
            "journal": self.journal,
            "doi": self.doi,
            "pmid": self.pmid,
            "pmcid": self.pmcid,
            "abstract": self.abstract,
            "source_api": self.source_api,
            "matched_queries": ";".join(self.matched_queries),
            "inferred_crispr_types": ";".join(self.inferred_crispr_types),
            "assay_keyword_hits": ";".join(self.assay_keyword_hits),
            "evidence_status": self.evidence_status,
            "evidence_reason": self.evidence_reason,
            "full_text_available": self.full_text_available,
        }


def classify(record: dict[str, Any], query_crispr_types: list[str]) -> Candidate:
    """Build a Candidate and assign a conservative evidence_status."""
    title = (record.get("title") or "").lower()
    abstract = (record.get("abstractText") or "").lower()
    text = f"{title} {abstract}"

    review_hit = any(m in text for m in REVIEW_MARKERS)
    comp_hit = any(m in text for m in COMPUTATIONAL_MARKERS)
    assay_hits = [t for t in ASSAY_TERMS if t in text]

    # Priority: review markers are decisive for review_or_secondary.
    if review_hit:
        status, reason = "review_or_secondary", "review marker present"
    elif assay_hits and not comp_hit:
        status, reason = "candidate_experimental", f"assay term(s): {', '.join(assay_hits)}"
    elif assay_hits and comp_hit:
        status, reason = "uncertain", f"assay {', '.join(assay_hits)} but also computational marker"
    elif comp_hit:
        status, reason = "computational_only", "computational marker present, no assay term"
    else:
        status, reason = "uncertain", "no decisive assay/review/computational marker"

    doi = (record.get("doi") or "").strip().lower()
    return Candidate(
        record_id=str(record.get("id") or record.get("pmid") or ""),
        title=record.get("title") or "",
        authors=record.get("authorString") or "",
        year=str(record.get("pubYear") or ""),
        journal=(record.get("journalInfo") or {}).get("journal", {}).get("title", "") if isinstance(record.get("journalInfo"), dict) else "",
        doi=doi,
        pmid=str(record.get("pmid") or ""),
        pmcid=str(record.get("pmcid") or ""),
        abstract=record.get("abstractText") or "",
        source_api="europepmc",
        matched_queries=[query_crispr_types and query_crispr_types[0] or "?"],
        inferred_crispr_types=list(query_crispr_types) if query_crispr_types else [],
        assay_keyword_hits=assay_hits,
        evidence_status=status,
        evidence_reason=reason,
        full_text_available="yes" if record.get("fullTextIdList") else "no",
    )


def deduplicate(candidates: list[Candidate]) -> list[Candidate]:
    """DOI-first, then PMID/PMCID, then normalized-title dedup."""
    seen_doi: dict[str, Candidate] = {}
    seen_id: dict[str, Candidate] = {}
    seen_title: dict[str, Candidate] = {}
    merged: list[Candidate] = []

    for c in candidates:
        key_doi = c.doi or None
        key_id = f"pmid:{c.pmid}" if c.pmid else (f"pmcid:{c.pmcid}" if c.pmcid else None)
        key_title = _norm_title(c.title)

        target = None
        if key_doi and key_doi in seen_doi:
            target = seen_doi[key_doi]
        elif key_id and key_id in seen_id:
            target = seen_id[key_id]
        elif key_title and key_title in seen_title:
            target = seen_title[key_title]

        if target is None:
            merged.append(c)
            if key_doi:
                seen_doi[key_doi] = c
            if key_id:
                seen_id[key_id] = c
            if key_title:
                seen_title[key_title] = c
        else:
            # Merge provenance: keep a superset of matched queries/types.
            target.matched_queries = list(dict.fromkeys(target.matched_queries + c.matched_queries))
            target.inferred_crispr_types = list(dict.fromkeys(target.inferred_crispr_types + c.inferred_crispr_types))

    return merged
