"""Batch download JATS fullTextXML + extract Cas->PAM statements.

Pipeline (HANDOFF todo #1+#2):
  1. Read a merged candidate CSV, collect PMCIDs (with fullText).
  2. Download each paper's fullTextXML into data/raw/fulltext/ (rate-limited,
     retried, cached locally).
  3. Extract: framed 5'-PAM-3' strings, Cas->PAM pairs, and the italicized
     species/strain name (e.g. "Cas9 from *Streptococcus pyogenes*) bound to
     each statement where possible.
  4. Emit a structured TSV with provenance (pmcid, doi, citation, raw_ref).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .extract_cas_pam_text import (
    find_pam_strings,
    find_cas_pam_pairs,
    find_cas_pam_species_triples,
    text_of,
)
from .apiclient import USER_AGENT, RETRYABLE

FULLTEXT_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"

# Italicized binomial species names: <italic>Streptococcus pyogenes</italic>
_ITALIC = re.compile(r"<italic>([^<]{4,80})</italic>", re.S)
_BINOMIAL = re.compile(r"([A-Z][a-z]+)\s+([a-z]{3,})")


_throttle_lock = threading.Lock()
_throttle_last = 0.0


def _acquire_slot(rate_limit: float) -> None:
    """Global rate limiter so parallel workers stay polite to the API."""
    global _throttle_last
    with _throttle_lock:
        elapsed = time.monotonic() - _throttle_last
        if elapsed < rate_limit:
            time.sleep(rate_limit - elapsed)
        _throttle_last = time.monotonic()


def download_fulltext(pmcid: str, dest_dir: Path, rate_limit: float,
                      max_retries: int, timeout: float) -> Path | None:
    """Download one fullTextXML, cache by PMCID. Returns path or None on fail."""
    dest = dest_dir / f"{pmcid}.xml"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = f"{FULLTEXT_BASE}/{pmcid}/fullTextXML"
    last_exc: Exception | None = None
    for attempt in range(max_retries + 1):
        _acquire_slot(rate_limit)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            if not data.strip():
                return None
            dest.write_bytes(data)
            return dest
        except urllib.error.HTTPError as e:
            last_exc = e
            if e.code in RETRYABLE and attempt < max_retries:
                time.sleep(2.0 ** attempt)
                continue
            return None
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last_exc = e
            if attempt < max_retries:
                time.sleep(2.0 ** attempt)
                continue
            return None
    return None


def _species_from_italic(xml: str) -> list[str]:
    """Distinct italicized binomial species names in the article."""
    out: list[str] = []
    seen: set[str] = set()
    for m in _ITALIC.finditer(xml):
        s = m.group(1)
        mm = _BINOMIAL.search(s)
        if mm:
            sp = f"{mm.group(1)} {mm.group(2)}"
            if sp not in seen:
                seen.add(sp)
                out.append(sp)
    return out


def extract_record(pmcid: str, xml_path: Path) -> dict:
    xml = xml_path.read_text(encoding="utf-8", errors="ignore")
    text = text_of(xml)
    triples = find_cas_pam_species_triples(xml)
    return {
        "pmcid": pmcid,
        "pam_strings": find_pam_strings(text),
        "cas_pam_pairs": [(n, p) for n, p, _ in triples],
        "triples": triples,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--outdir", default="data/raw/fulltext", type=Path)
    ap.add_argument("--output", default="data/parsed/cas_pam_text/cas_pam_pairs.tsv", type=Path)
    ap.add_argument("--rate-limit", type=float, default=0.5)
    ap.add_argument("--max-retries", type=int, default=3)
    ap.add_argument("--timeout", type=float, default=60.0)
    ap.add_argument("--cap", type=int, default=0, help="0 = no cap")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(args.candidates, encoding="utf-8")))
    pmcids: list[str] = []
    seen: set[str] = set()
    for r in rows:
        pmcid = (r.get("pmcid") or "").strip()
        if pmcid.startswith("PMC") and pmcid not in seen:
            seen.add(pmcid)
            pmcids.append(pmcid)
    if args.cap:
        pmcids = pmcids[: args.cap]
    print(f"candidates: {len(pmcids)} PMCIDs to process")

    # Prefer candidate_experimental first (best signal), then the rest.
    exp = {r["pmcid"] for r in rows if r.get("evidence_status") == "candidate_experimental"}
    pmcids = sorted(pmcids, key=lambda p: (p not in exp, p))
    meta = {r["pmcid"]: r for r in rows if r.get("pmcid")}

    downloaded = 0
    failed = 0
    no_data = 0
    n_pairs = 0
    writers: list[dict] = []

    # Phase 1: parallel download (cached files are skipped instantly).
    paths: dict[str, Path | None] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {
            ex.submit(download_fulltext, p, args.outdir, args.rate_limit,
                      args.max_retries, args.timeout): p
            for p in pmcids
        }
        for i, fut in enumerate(as_completed(futs), 1):
            p = futs[fut]
            paths[p] = fut.result()
            if i % 50 == 0 or i == len(pmcids):
                print(f"  downloaded {i}/{len(pmcids)}", flush=True)

    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "pmcid", "doi", "title", "year", "journal", "citation",
            "cas_name", "pam", "species", "pam_source", "raw_ref", "sha256",
        ])
        for pmcid in pmcids:
            path = paths[pmcid]
            if path is None:
                failed += 1
                continue
            downloaded += 1
            rec = extract_record(pmcid, path)
            m = meta.get(pmcid, {})
            doi = (m.get("doi") or "").strip()
            title = (m.get("title") or "").strip()
            year = (m.get("year") or "").strip()
            journal = (m.get("journal") or "").strip()
            citation = "; ".join(x for x in [m.get("authors") or "", title, year] if x)
            raw_ref = f"{pmcid}::fullTextXML::body::prose"
            sha = hashlib.sha256(path.read_bytes()).hexdigest()[:16]

            if not rec["cas_pam_pairs"]:
                no_data += 1
            seen_rows: set[tuple] = set()
            for name, pam, species in rec["triples"]:
                key = (pmcid, name, pam, species)
                if key in seen_rows:  # drop intra-paper duplicate statements
                    continue
                seen_rows.add(key)
                n_pairs += 1
                w.writerow([pmcid, doi, title, year, journal, citation,
                            name, pam, species, "fullTextXML", raw_ref, sha])
                writers.append({
                    "pmcid": pmcid, "cas_name": name, "pam": pam,
                    "species": species, "doi": doi,
                })

    print(f"downloaded={downloaded} failed={failed} no_pairs={no_data}")
    print(f"total cas_pam pairs={n_pairs}")
    # Dedup report
    uniq = set()
    for x in writers:
        uniq.add((x["pmcid"], x["cas_name"], x["pam"]))
    print(f"unique (pmcid,cas,pam)={len(uniq)}")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
