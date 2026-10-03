"""Discovery runner: execute queries, classify, dedup, emit candidate table.

Usage:
  python src/pamdict/collect/run_discover.py \
      --queries src/pamdict/collect/queries_recent.json \
      --cache data/raw/epmc_cache \
      --output data/parsed/discovery \
      [--offline-only] [--cap 200]
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .apiclient import ApiClient
from .discover import Candidate, classify, deduplicate

FIELDS = [
    "record_id", "title", "authors", "year", "journal", "doi", "pmid", "pmcid",
    "abstract", "source_api", "matched_queries", "inferred_crispr_types",
    "assay_keyword_hits", "evidence_status", "evidence_reason", "full_text_available",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries", required=True, type=Path)
    ap.add_argument("--cache", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--offline-only", action="store_true")
    ap.add_argument("--cap", type=int, default=200)
    ap.add_argument("--page-size", type=int, default=100)
    args = ap.parse_args()

    queries = json.loads(args.queries.read_text(encoding="utf-8"))
    client = ApiClient(args.cache, offline_only=args.offline_only)

    all_candidates: list[Candidate] = []
    query_counts: dict[str, int] = {}
    for q in queries:
        qid = q["id"]
        resp = client.search_europepmc(q["query"], page_size=args.page_size)
        results = resp.get("resultList", {}).get("result", [])
        query_counts[qid] = len(results)
        for r in results:
            all_candidates.append(classify(r, q.get("crispr_types", [])))

    unique = deduplicate(all_candidates)
    emitted = unique[: args.cap]

    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "candidate_papers.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for c in emitted:
            w.writerow(c.to_dict())

    summary = {
        "queries": query_counts,
        "raw_results": sum(query_counts.values()),
        "unique_after_dedup": len(unique),
        "emitted": len(emitted),
        "evidence_status": {s: sum(1 for c in emitted if c.evidence_status == s) for s in
                            ["candidate_experimental", "review_or_secondary", "computational_only", "uncertain"]},
    }
    (args.output / "discovery_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
