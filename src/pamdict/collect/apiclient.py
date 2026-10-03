"""Minimal HTTP client for Europe PMC / Crossref with cache, retry, integrity.

Uses stdlib ``urllib`` (via a thin curl-free path) so it is dependency-light and
works in the sandbox. Responses are cached by SHA-256 of the exact URL; a
``--offline-only`` mode refuses to hit the network on a cache miss.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

USER_AGENT = "PAMdict-Collector/1.0 (research; contact: project maintainer)"

RETRYABLE = {429, 500, 502, 503, 504}


class ApiClient:
    def __init__(self, cache_dir: Path, offline_only: bool = False,
                 timeout: float = 30.0, max_retries: int = 3,
                 rate_limit_delay: float = 0.5) -> None:
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.offline_only = offline_only
        self.timeout = timeout
        self.max_retries = max_retries
        self.rate_limit_delay = rate_limit_delay
        self._last_request = 0.0

    def _cache_path(self, url: str) -> Path:
        digest = hashlib.sha256(url.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.json"

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_request
        if elapsed < self.rate_limit_delay:
            time.sleep(self.rate_limit_delay - elapsed)

    def get(self, url: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
        if params:
            sep = "&" if "?" in url else "?"
            url = url + sep + urllib.parse.urlencode(params)

        cache_path = self._cache_path(url)
        if cache_path.exists():
            return json.loads(cache_path.read_text(encoding="utf-8"))

        if self.offline_only:
            raise RuntimeError(f"offline-only cache miss for {url}")

        last_exc: Exception | None = None
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    body = resp.read().decode("utf-8", errors="replace")
                payload = json.loads(body)
                cache_path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
                return payload
            except urllib.error.HTTPError as e:
                last_exc = e
                if e.code in RETRYABLE and attempt < self.max_retries:
                    time.sleep(2.0 ** attempt)
                    continue
                if e.code not in RETRYABLE:
                    raise RuntimeError(f"non-retryable HTTP {e.code} for {url}") from e
            except (urllib.error.URLError, TimeoutError) as e:
                last_exc = e
                if attempt < self.max_retries:
                    time.sleep(2.0 ** attempt)
                    continue
        raise RuntimeError(f"failed after retries for {url}: {last_exc}")

    def search_europepmc(self, query: str, page_size: int = 25) -> dict[str, Any]:
        return self.get(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            params={
                "query": query,
                "format": "json",
                "resultType": "core",
                "pageSize": str(page_size),
            },
        )
