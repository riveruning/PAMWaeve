"""Extract the "five elements" of a PAM training sample from literature.

The five elements that let a paper contribute a training row:
  1. protein_sequence  (the CRISPR-Cas effector/protein)
  2. pam_consensus     (the reported PAM motif)
  3. pam_logo_acgt     (base counts / probability / info matrix, when present)
  4. organism          (source organism; maps to source column)
  5. method/assay      (evidence type; maps to citation + evidence metadata)

This module parses two input shapes:
  - Europe PMC fullTextXML (``<result><fullTextUrl>`` OpenAccess XML), which
    contains body text we can regex/heuristic-scan.
  - Supplementary tables (CSV/TSV/XLSX) that list PAM motifs per system.

NOTE: automatic extraction produces CANDIDATES. Each emitted row stays
``unreviewed`` until a human confirms the five elements against primary text.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

# A PAM motif is a short IUPAC DNA string. We accept a bounded set of symbols.
IUPAC_OK = set("ACGTRYSWKMBDHVNacgtryswkmbdhvn")

# Common patterns for "PAM is NGG", "PAM = TTTN", "PAM: NGG", "NGG PAM".
PAM_PATTERNS = [
    re.compile(r"PAM\s+(?:is|was|of|:|＝|=)\s*[\"']?([ACGTURYSWKMBDHVN]{2,10})", re.I),
    re.compile(r"[\"']?([ACGTURYSMWKMBDHVN]{2,10})[\"']?\s+PAM", re.I),
    re.compile(r"PAM\s*\(?([ACGTURYSWKMBDHVN]{2,10})\)?", re.I),
]

# Protein accession / UniProt-like identifiers.
ACCESSION_RE = re.compile(r"\b(?:UniProt(?:KB)?|accession|Acc\.?|ID)[: ]\s*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-Z][0-9]{5}|[A-Z]{2}_[A-Z0-9]+)", re.I)
UNIPROT_RE = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9])\b")

ORGANISM_RE = re.compile(
    r"\b([A-Z][a-z]+(?:\s+(?:[a-z]+|sp\.|subsp\.))*(?:\s+[A-Z][a-z]+)?)\b")
TAXON_HINTS = re.compile(r"(Streptococcus|Staphylococcus|Bacillus|Escherichia|Francisella|"
                        r"Lachnospiraceae|Acidaminococcus|Campylobacter|Neisseria|"
                        r"Pseudomonas|Pyrococcus|Sulfolobus|Candidatus)", re.I)


@dataclass
class ExtractedSample:
    """A candidate five-element sample, pre-review."""
    pam_consensus: str = ""
    pam_logo_acgt: str = ""          # JSON array string or empty
    protein_sequence: str = ""
    protein_id: str = ""
    organism: str = ""
    method: str = ""
    cas_family: str = ""             # e.g. Cas9, Cas12a
    crispr_type: str = ""            # Type I / II / V ...
    citation: str = ""
    doi: str = ""
    evidence_status: str = "unreviewed"
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "pam_consensus": self.pam_consensus,
            "pam_logo_acgt": self.pam_logo_acgt,
            "protein_sequence": self.protein_sequence,
            "protein_id": self.protein_id,
            "organism": self.organism,
            "method": self.method,
            "cas_family": self.cas_family,
            "crispr_type": self.crispr_type,
            "citation": self.citation,
            "doi": self.doi,
            "evidence_status": self.evidence_status,
            "notes": ";".join(self.notes),
        }


def sanitize_iupac(motif: str) -> str:
    """Uppercase and strip characters outside the IUPAC DNA alphabet."""
    cleaned = "".join(c for c in motif.upper() if c in IUPAC_OK)
    # T -> U normalization is not needed; keep T.
    return cleaned


def extract_pam_motifs(text: str) -> list[str]:
    """Return candidate PAM motifs found in text (deduped, order-preserving)."""
    found: list[str] = []
    for pat in PAM_PATTERNS:
        for m in pat.finditer(text):
            motif = sanitize_iupac(m.group(1))
            if 2 <= len(motif) <= 10 and motif not in found:
                found.append(motif)
    return found


def extract_uniprot(text: str) -> list[str]:
    return list(dict.fromkeys(UNIPROT_RE.findall(text)))


def extract_organism(text: str) -> str:
    m = TAXON_HINTS.search(text)
    return m.group(1) if m else ""


def infer_cas_family(text: str) -> str:
    fams = ["Cas12a", "Cas12b", "Cas12f", "Cas12i", "Cas12j", "Cas9", "Cas8",
            "Cas10", "Cas3", "Cas13", "Cpf1"]
    for fam in fams:
        if re.search(rf"\b{fam}\b", text):
            return fam
    return ""


def infer_crispr_type(text: str) -> str:
    for typ, rx in [("Type I", r"\bType\s*I\b"), ("Type II", r"\bType\s*II\b"),
                    ("Type III", r"\bType\s*III\b"), ("Type IV", r"\bType\s*IV\b"),
                    ("Type V", r"\bType\s*V\b"), ("Type VI", r"\bType\s*VI\b")]:
        if re.search(rx, text, re.I):
            return typ
    return ""


def parse_fulltext_xml(xml_text: str, doi: str = "", citation: str = "") -> ExtractedSample:
    """Extract from a Europe PMC fullTextXML string."""
    # Strip tags to get plain text for regex scanning.
    body = re.sub(r"<[^>]+>", " ", xml_text)
    body = re.sub(r"\s+", " ", body)

    motifs = extract_pam_motifs(body)
    sample = ExtractedSample(
        pam_consensus=motifs[0] if motifs else "",
        pam_logo_acgt="",
        protein_id=(extract_uniprot(body) or [""])[0],
        organism=extract_organism(body),
        cas_family=infer_cas_family(body),
        crispr_type=infer_crispr_type(body),
        citation=citation,
        doi=doi,
        method="fulltext_heuristic",
    )
    if len(motifs) > 1:
        sample.notes.append(f"multiple motifs: {', '.join(motifs)}")
    return sample


def parse_supp_table(rows: list[dict[str, str]], doi: str = "",
                     citation: str = "") -> list[ExtractedSample]:
    """Parse a supplementary table listing {motif, [counts], [organism]} rows.

    Accepts a list of dicts with any of: pam/motif, organism/species,
    count/n, accession, sequence. Returns one ExtractedSample per motif row.
    """
    samples: list[ExtractedSample] = []
    for row in rows:
        motif = sanitize_iupac(row.get("pam") or row.get("motif") or "")
        if not (2 <= len(motif) <= 10):
            continue
        samples.append(ExtractedSample(
            pam_consensus=motif,
            organism=row.get("organism") or row.get("species") or "",
            protein_id=row.get("accession") or row.get("uniprot") or "",
            protein_sequence=row.get("sequence") or "",
            citation=citation,
            doi=doi,
            method="supp_table",
        ))
    return samples
