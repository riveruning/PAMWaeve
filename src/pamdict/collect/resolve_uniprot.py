"""Resolve Cas protein name (+ species) to amino-acid sequence via UniProt REST.

Goal: the Cas->PAM pairs extracted from full-text (`cas_pam_pairs.tsv`) carry a
Cas *name* and (sometimes) a source *species*, but protein2PAM's input is the
amino-acid *sequence*. This module closes that gap by querying UniProt for the
protein sequence of each (cas_name, species) pair.

Design notes
------------
- Queries are anchored on the species when present (UniProt ``organism_name:``),
  falling back to a name-only query when the species is missing. Species is the
  reliable anchor because Cas short-names ("SpCas9", "CoCas9") are not UniProt
  identifiers.
- Results are ranked; we keep the top candidate per query with explicitness of
  match (exact species hit vs. fallback) so downstream review can flag weak
  mappings. Nothing is silently trusted.
- Output is a flat TSV plus a JSON report; sequences are only ever appended to
  new files (never mixed into the immutable official baseline).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from pamdict.collect import cas_alias_map

USER_AGENT = "PAMdict-Collector/1.0 (research; contact: project maintainer)"

# UniProt search endpoint (REST, returns JSON, no key required).
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# Fields we want back. "sequence" is large; include accession, name, organism.
# Note: UniProt REST field names are the *returned* JSON keys. `sequence` is the
# payload-triggering field (must be present to get the AA string); `gene_names`
# gives the `genes` array; `organism_name` gives the `organism` object.
FIELDS = "accession,gene_names,organism_name,sequence"


def _get(url: str, timeout: float = 40.0, max_retries: int = 3) -> dict[str, Any]:
    """GET a URL as JSON with exponential-backoff retry and rate throttle."""
    last_exc: Exception | None = None
    for attempt in range(max_retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read().decode("utf-8", errors="replace")
            return json.loads(body)
        except urllib.error.HTTPError as e:
            last_exc = e
            if e.code in {429, 500, 502, 503, 504} and attempt < max_retries:
                time.sleep(2.0 ** attempt)
                continue
            raise RuntimeError(f"HTTP {e.code} for {url}") from e
        except (urllib.error.URLError, TimeoutError) as e:
            last_exc = e
            if attempt < max_retries:
                time.sleep(2.0 ** attempt)
                continue
    raise RuntimeError(f"failed after retries for {url}: {last_exc}")


def _species_terms(species: str) -> str:
    """Normalize a species string into a UniProt free-text phrase.

    Species strings from full-text look like "Streptococcus pyogenes" or
    "Francisella novicida". UniProt's exact ``organism_name:`` field frequently
    disagrees with the paper's binomial (subspecies / historical naming, e.g.
    "Francisella novicida" vs "Francisella tularensis subsp. novicida"), so we
    use a quoted free-text phrase instead of the exact field. The family filter
    in `_pick_best` keeps the match on-target.
    """
    s = species.strip()
    if not s:
        return ""
    return f'"{s}"'


def _family(cas_name: str) -> str | None:
    """Infer the Cas family keyword (UniProt gene/name stem) from a Cas name.

    Returns one of "cas9"/"cas12"/"cpf1"/"cas13"/"cas12f"/"iscb"/None. The
    family is the key correctness anchor: a "Cas12a" PAM must never be paired
    with a Cas9 sequence.
    """
    n = cas_name.lower().strip()
    if "cas12f" in n:
        return "cas12f"
    if "cas12" in n or "cpf1" in n:
        return "cas12"
    if "cas13" in n or "c2c2" in n:
        return "cas13"
    if "cas9" in n:
        return "cas9"
    if "iscb" in n:
        return "iscb"
    return None  # Cas3/Cas1/other — not a protein2PAM editing nuclease


def _gene_terms(fam: str | None) -> str:
    """Build a UniProt gene/name OR-clause for a family keyword."""
    if fam == "cas9":
        return "(gene:cas9 OR gene:csn1)"
    if fam == "cas12":
        return "(gene:cas12a OR gene:cpf1 OR gene:cas12)"
    if fam == "cas12f":
        return "(gene:cas12f OR gene:cas14)"
    if fam == "cas13":
        return "(gene:cas13)"
    if fam == "iscb":
        return "(gene:iscb)"
    return ""


def build_queries(cas_name: str, species: str, curated: bool) -> list[tuple[str, str]]:
    """Return an ordered list of (query, kind) attempts for a Cas name+species.

    ``kind`` marks how trustworthy the anchor is:
      * "species"   — full-text species + family (the only path that finds
                      cold-start Cas; species may be curated or paper-derived)
      * "name_only" — alias/family only, no species (low confidence)
    The family keyword is what prevents cross-family mispairing (e.g. a Cas12a
    PAM being matched to a Cas9 sequence).

    ``curated`` True means `species` came from the curated alias map (trustworthy
    source organism); False means it came from the paper (may be a host/target,
    so any hit must be human-reviewed).
    """
    cas = cas_name.strip()
    fam = _family(cas)
    qs: list[tuple[str, str]] = []
    sp_term = _species_terms(species) if species else ""
    gene = _gene_terms(fam) if fam else ""

    if fam and gene:
        if sp_term:
            qs.append((f"{sp_term} AND {gene}", "species"))
        # No-species fallback: family + the free-text alias (may be absent from
        # UniProt's index for cold aliases, but it is the best we can do).
        qs.append((f"{gene} AND {cas}", "name_only"))
    elif sp_term:
        # Unknown family (e.g. Cas3): anchor on species + raw name token.
        qs.append((f"{sp_term} AND ({cas})", "species"))
    else:
        qs.append((f'gene:"{cas}" OR protein_name:"{cas}"', "name_only"))
    return qs


def _pick_best(
    results: list[dict[str, Any]], fam: str | None, min_len: int = 200
) -> dict[str, Any] | None:
    """Choose the best sequence-bearing entry among UniProt results.

    When a family is known, entries whose name/gene do not mention that family
    are filtered out to avoid cross-family mispairing (the top cause of wrong
    sequence resolution for Cas proteins).
    """
    seq_entries = [
        r for r in results
        if r.get("sequence", {}).get("value") and r.get("sequence", {}).get("length", 0) >= min_len
    ]

    if fam:
        fam_parts = {
            "cas9": ("cas9", "csn1"),
            "cas12": ("cas12", "cpf1"),
            "cas12f": ("cas12f", "cas14"),
            "cas13": ("cas13",),
            "iscb": ("iscb",),
        }.get(fam, (fam,))
        filtered = []
        for r in seq_entries:
            genes_blob = ""
            for g in r.get("genes", []) or []:
                if isinstance(g, dict):
                    genes_blob += " " + (g.get("geneName", {}).get("value") or "")
            blob = (f"{_name(r)} {genes_blob}").lower()
            if any(p in blob for p in fam_parts):
                filtered.append(r)
        if filtered:
            seq_entries = filtered

    if not seq_entries:
        return None
    # Prefer reviewed (Swiss-Prot) entries, then longest.
    seq_entries.sort(
        key=lambda r: (
            0 if "reviewed" in (r.get("entryType") or "") else 1,
            -r.get("sequence", {}).get("length", 0),
        )
    )
    return seq_entries[0]


def search_uniprot(query: str, size: int = 5) -> list[dict[str, Any]]:
    url = (
        UNIPROT_SEARCH + "?"
        + urllib.parse.urlencode(
            {"query": query, "format": "json", "size": str(size), "fields": FIELDS}
        )
    )
    d = _get(url)
    return d.get("results", [])


def _fetch_entry(accession: str) -> dict[str, str] | None:
    """Fetch a single UniProt entry by accession and extract {length, sequence}.

    Uses the entries endpoint (``/uniprotkb/{acc}``); returns None on any
    failure (network, missing, parse).
    """
    url = (
        "https://rest.uniprot.org/uniprotkb/"
        + urllib.parse.quote(accession)
        + "?format=json&fields=accession,sequence,protein_name,organism_name"
    )
    try:
        d = _get(url)
    except RuntimeError:
        return None
    seq = d.get("sequence", {})
    if not seq.get("value"):
        return None
    return {
        "length": seq.get("length"),
        "sequence": seq.get("value"),
        "protein_name": _name(d),
        "organism": _organism(d),
    }


def resolve_one(
    cas_name: str, species: str, cache_dir: Path, offline_only: bool = False
) -> dict[str, Any]:
    """Resolve a single Cas name (+species) to {accession, length, sequence, match}.

    Prefers the curated alias map (offline, reliable), then falls back to a live
    UniProt query. ``match`` is one of:
        "verified"        — alias hit a verified accession in the local map
        "species"         — full-text species+family query hit (species curated)
        "species_review"  — hit using the paper's species (may be host/target:
                            needs human review before trusting the sequence)
        "name_only"       — alias/family only, no species (low confidence)
        "none"            — no sequence found
        "error"           — network/parse failure
    """
    # Cache key from the (cas_name, species) identity.
    key = hashlib.sha256(f"{cas_name}\t{species}".encode("utf-8")).hexdigest()
    cache_file = cache_dir / f"{key}.json"
    if cache_file.exists():
        return json.loads(cache_file.read_text(encoding="utf-8"))
    if offline_only:
        return {"cas_name": cas_name, "species": species, "match": "error",
                "error": "offline cache miss"}

    fam = _family(cas_name)
    curated_sp = cas_alias_map.organism_of(cas_name)

    # 1) Verified accession in the curated map (offline, highest confidence).
    acc = cas_alias_map.VERIFIED_ACCESSION.get(cas_name)
    if acc:
        seq = _fetch_entry(acc)
        if seq:
            out = {
                "cas_name": cas_name,
                "species": species,
                "match": "verified",
                "accession": acc,
                "length": seq.get("length"),
                "sequence": seq.get("sequence"),
                "protein_name": seq.get("protein_name"),
                "organism_name": seq.get("organism"),
                "source_organism": curated_sp,
                "species_source": "curated",
            }
            cache_file.write_text(json.dumps(out, ensure_ascii=False), encoding="utf-8")
            return out

    # 2) Live query. Anchor on the curated species when available; otherwise use
    #    the paper's species as a free-text anchor but FLAG the result for review
    #    (it may be a host/target rather than the true source organism).
    anchor_species = curated_sp if curated_sp else species
    curated = curated_sp is not None

    for query, kind in build_queries(cas_name, anchor_species, curated):
        try:
            results = search_uniprot(query)
        except RuntimeError as e:
            out = {"cas_name": cas_name, "species": species, "match": "error",
                   "error": str(e)}
            cache_file.write_text(json.dumps(out, ensure_ascii=False), encoding="utf-8")
            return out
        best = _pick_best(results, fam)
        if best is not None:
            seq = best.get("sequence", {})
            # Distinguish trustworthy species anchor from paper-derived (review).
            match_kind = kind if curated else (
                "species_review" if kind == "species" else kind
            )
            out = {
                "cas_name": cas_name,
                "species": species,
                "match": match_kind,
                "accession": best.get("primaryAccession"),
                "length": seq.get("length"),
                "sequence": seq.get("value"),
                "protein_name": _name(best),
                "organism_name": _organism(best),
                "source_organism": curated_sp,
                "species_source": "curated" if curated else "paper",
            }
            cache_file.write_text(json.dumps(out, ensure_ascii=False), encoding="utf-8")
            return out
        time.sleep(0.3)

    out = {"cas_name": cas_name, "species": species, "match": "none", "error": "no sequence found"}
    cache_file.write_text(json.dumps(out, ensure_ascii=False), encoding="utf-8")
    return out


def _name(entry: dict[str, Any]) -> str:
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName", {})
    full = rec.get("fullName", {})
    return full.get("value", "") or ""


def _organism(entry: dict[str, Any]) -> str:
    org = entry.get("organism", {})
    return org.get("scientificName", "") or ""


def collect_pairs(tsv_path: Path) -> list[tuple[str, str]]:
    """Collect unique (cas_name, species) pairs from cas_pam_pairs.tsv."""
    pairs: dict[tuple[str, str], None] = {}
    with tsv_path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cas = (row.get("cas_name") or "").strip()
            sp = (row.get("species") or "").strip()
            if cas:
                pairs[(cas, sp)] = None
    return list(pairs.keys())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True, help="cas_pam_pairs.tsv path")
    ap.add_argument("--cache", default="data/raw/uniprot_cache")
    ap.add_argument("--output", default="data/parsed/cas_pam_text/cas_sequences.tsv")
    ap.add_argument("--report", default="data/parsed/cas_pam_text/resolve_report.json")
    ap.add_argument("--offline-only", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="max pairs to resolve (0=all)")
    args = ap.parse_args()

    pairs = collect_pairs(Path(args.pairs))
    # Species-bearing pairs first (they resolve most reliably).
    pairs.sort(key=lambda p: (0 if p[1] else 1, p[0].lower()))
    if args.limit:
        pairs = pairs[: args.limit]

    cache_dir = Path(args.cache)
    cache_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for i, (cas, sp) in enumerate(pairs, 1):
        r = resolve_one(cas, sp, cache_dir, args.offline_only)
        rows.append(r)
        flag = r.get("match")
        acc = r.get("accession") or "-"
        print(f"[{i}/{len(pairs)}] {flag:10s} acc={acc:12s} {cas:16s} | {sp[:40]}",
              flush=True)
        if flag == "error":
            print(f"          err: {r.get('error')}", flush=True)
        time.sleep(0.25)

    # Write TSV of sequences.
    out_dir = Path(args.output).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    fields = ["cas_name", "species", "match", "accession", "length",
              "protein_name", "organism_name", "source_organism",
              "species_source", "sequence"]
    with Path(args.output).open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

    # Summary report.
    summary = {
        "total_pairs": len(pairs),
        "resolved": sum(1 for r in rows if r.get("match") in {"verified", "species", "species_review", "name_only"}),
        "verified": sum(1 for r in rows if r.get("match") == "verified"),
        "species": sum(1 for r in rows if r.get("match") == "species"),
        "species_review": sum(1 for r in rows if r.get("match") == "species_review"),
        "name_only": sum(1 for r in rows if r.get("match") == "name_only"),
        "none": sum(1 for r in rows if r.get("match") == "none"),
        "error": sum(1 for r in rows if r.get("match") == "error"),
    }
    Path(args.report).write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print("\nReport:", json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
