"""Fetch supplementary data files (XLSX/CSV/TSV) for a candidate paper.

The core empirical finding (2026-09-05) is that PAM characterization data
almost always lives in *supplementary data files*, not in JATS XML body text.
This module locates and downloads those files.

Sources (in order of preference):
  1. Europe PMC supplementaryFiles endpoint (if listed).
  2. PMC "Supplementary Data" links from the fullTextXML supplementary-material.
  3. bioRxiv/medRxiv media links.

Each downloaded file gets a manifest entry with URL, size, sha256.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path
from urllib.parse import urljoin

from .apiclient import RETRYABLE, USER_AGENT


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, dest: Path, max_retries: int = 3,
                  timeout: float = 60.0) -> dict:
    """Download a binary file with retry+backoff for transient failures.

    Europe PMC's supplementaryFiles endpoint is intermittent (documented flaky
    proxy/network), so we retry timeouts, connection resets and 5xx a few times
    with exponential backoff before surfacing the error. Non-retryable 4xx
    (e.g. 404) fail immediately.
    """
    last_exc: Exception | None = None
    for attempt in range(max_retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
            return {
                "url": url,
                "local_path": str(dest),
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            }
        except urllib.error.HTTPError as e:
            last_exc = e
            if e.code in RETRYABLE and attempt < max_retries:
                time.sleep(2.0 ** attempt)
                continue
            raise
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last_exc = e
            if attempt < max_retries:
                time.sleep(2.0 ** attempt)
                continue
            raise
    raise RuntimeError(f"failed after retries for {url}: {last_exc}")


def extract_supplementary_links(xml_text: str, base_url: str) -> list[dict]:
    """Find supplementary-material <a>/<oasis:table> links in JATS XML."""
    links: list[dict] = []
    # <supplementary-material xlink:href="..."> with <media xlink:href="...">
    for m in re.finditer(
        r'<supplementary-material[^>]*xlink:href="([^"]+)"[^>]*>(.*?)</supplementary-material>',
        xml_text, re.S,
    ):
        href, inner = m.group(1), m.group(2)
        medias = re.findall(r'<media[^>]*xlink:href="([^"]+)"', inner)
        if medias:
            for mm in medias:
                links.append({"url": urljoin(base_url, mm), "label": href})
        else:
            links.append({"url": urljoin(base_url, href), "label": href})
    return links
