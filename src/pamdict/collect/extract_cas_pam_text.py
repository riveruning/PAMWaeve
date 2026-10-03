"""Extract "Cas protein -> PAM sequence" statements from JATS full-text XML.

A paper that engineers a bacterium or uses a specific CRISPR system almost
always states, in prose, which Cas protein was used and which PAM it
recognizes (e.g. "SpCas9 requires a 5'-NGG-3' PAM").  These statements are
clean one-to-one (Cas -> PAM) triples, unlike supplementary tables which are
target-site enumerations.

This module pulls (a) the raw PAM strings with 5'/3' framing and (b) pairs of
(Cas-name, PAM) that appear close together, as candidates for downstream
curation.
"""
from __future__ import annotations

import re
from pathlib import Path

# IUPAC nucleotide codes (DNA), including ambiguity codes.
_IUPAC = "ACGTUNRYWSKMBDHV"

# Uppercase strings that are real English words (or start codons) that the
# IUPAC alphabet can accidentally spell out.  These are NOT PAM sequences.
_FALSE_PAM = {"AND", "ATG", "TAG", "TAA", "TGA", "ARE", "ANY", "AC", "CA", "NOT", "GAT"}


def text_of(xml: str) -> str:
    """Strip XML tags to plain text, collapse whitespace."""
    t = re.sub(r"<[^>]+>", " ", xml)
    return re.sub(r"\s+", " ", t)


# Frame like 5'-NGG-3' / 5′-NNGRRT-3′ (curly quotes, prime marks).  Note the
# actual text is "5 ’ - NRRWC - 3 ’" (a hyphen separates the quote from the
# bases and the "3").  Use actual Unicode quote chars, not \u escapes in a raw
# string (which would match a literal backslash).
_QUOTES = "'\u2019\u2018\u2032\u2033"
_FRAME = re.compile(
    "5\\s*[" + _QUOTES + "]\\s*-?\\s*([" + _IUPAC + "]{2,10})\\s*-?\\s*3\\s*["
    + _QUOTES + "]",
    re.IGNORECASE,
)


def _is_real_pam(s: str) -> bool:
    """Heuristic gate: drop start codons / English words the alphabet can spell."""
    return len(s) >= 3 and s not in _FALSE_PAM


def find_pam_strings(text: str) -> list[str]:
    """All framed PAM-like strings (5'-...-3'), uppercased, dedup in order."""
    out: list[str] = []
    seen: set[str] = set()
    for m in _FRAME.finditer(text):
        s = m.group(1).upper()
        if s not in seen and _is_real_pam(s):
            seen.add(s)
            out.append(s)
    return out


# Italicized binomial species names in the *raw* XML.
_ITALIC = re.compile(r"<italic>([^<]{4,80})</italic>", re.S)
_BINOMIAL = re.compile(r"\b([A-Z][a-z]+)\s+([a-z]{3,})\b")


def find_species_spans(xml: str) -> list[tuple[int, int, str]]:
    """(start, end, name) spans of italicized binomial names, on plain text.

    Positions are computed on the tag-stripped text so a caller can find the
    binomial closest to a PAM statement by character distance.
    """
    spans: list[tuple[int, int, str]] = []
    text = text_of(xml)
    for m in _ITALIC.finditer(xml):
        raw = m.group(1)
        mm = _BINOMIAL.search(raw)
        if not mm:
            continue
        name = f"{mm.group(1)} {mm.group(2)}"
        # Locate this name in the plain text, choose the occurrence closest to
        # the tag's position.  A global .find() picks the first; refine by
        # searching from an approximate mapped offset.
        approx = int(m.start() / max(1, len(xml)) * len(text))
        plain_pos = text.find(name, max(0, approx - 500))
        if plain_pos < 0:
            plain_pos = text.find(name)
        if plain_pos < 0:
            continue
        spans.append((plain_pos, plain_pos + len(name), name))
    return spans


def species_at(text_span_center: int, spans: list[tuple[int, int, str]],
               max_dist: int = 400) -> str:
    """Nearest italic binomial to a text position; empty string if none close."""
    best = ""
    best_d = max_dist + 1
    for s, e, name in spans:
        d = max(s - text_span_center, text_span_center - e, 0)
        if d < best_d:
            best_d = d
            best = name
    return best if best_d <= max_dist else ""


# A Cas-like protein token (SpCas9, Nme1Cas9, Cpf1, IscB, TnpB, Cas12a, ...).
_CAS = re.compile(
    r"\b([A-Za-z][A-Za-z0-9]*(?:Cas\d*[A-Za-z]?|Cpf1|IscB|TnpB|Cas13[abc]?|"
    r"Cas12[a-z]?|Cascade[Cc]as[A-Za-z0-9]*))\b"
)


def find_cas_pam_pairs(text: str, window: int = 80) -> list[tuple[str, str]]:
    """Pairs of (Cas-name, PAM-string) that co-occur within `window` chars.

    Heuristic: walk all framed PAMs, look backwards for the nearest Cas token,
    and forwards for a nearby Cas token, to bind a name to the motif.
    """
    pairs: list[tuple[str, str]] = []
    for m in _FRAME.finditer(text):
        pam = m.group(1).upper()
        lo = max(0, m.start() - window)
        hi = min(len(text), m.end() + window)
        before = text[lo:m.start()]
        after = text[m.end():hi]
        name = None
        mb = _CAS.findall(before)
        if mb:
            name = mb[-1]
        else:
            ma = _CAS.findall(after)
            if ma:
                name = ma[0]
        if name and _is_real_pam(pam):
            pairs.append((name, pam))
    return pairs


def find_cas_pam_species_triples(xml: str, window: int = 80) -> list[tuple[str, str, str]]:
    """(cas, pam, species) triples with the nearest italic binomial bound.

    Species binding is per-statement (not global): each framed PAM looks for the
    binomial closest to its own span, fixing the earlier bug where every row got
    the article's first italicized species.
    """
    text = text_of(xml)
    spans = find_species_spans(xml)
    triples: list[tuple[str, str, str]] = []
    for m in _FRAME.finditer(text):
        pam = m.group(1).upper()
        if not _is_real_pam(pam):
            continue
        lo = max(0, m.start() - window)
        hi = min(len(text), m.end() + window)
        before = text[lo:m.start()]
        after = text[m.end():hi]
        name = None
        mb = _CAS.findall(before)
        if mb:
            name = mb[-1]
        else:
            ma = _CAS.findall(after)
            if ma:
                name = ma[0]
        if not name:
            continue
        center = (m.start() + m.end()) // 2
        sp = species_at(center, spans)
        triples.append((name, pam, sp))
    return triples


def extract_file(path: str | Path) -> dict:
    p = Path(path)
    text = text_of(p.read_text(encoding="utf-8", errors="ignore"))
    return {
        "file": str(p),
        "pam_strings": find_pam_strings(text),
        "cas_pam_pairs": find_cas_pam_pairs(text),
    }


if __name__ == "__main__":
    import sys
    for arg in sys.argv[1:]:
        r = extract_file(arg)
        print(f"== {r['file']} ==")
        print("  framed PAM strings (%d):" % len(r["pam_strings"]))
        for s in r["pam_strings"][:30]:
            print("    -", s)
        print("  Cas->PAM pairs (%d):" % len(r["cas_pam_pairs"]))
        for name, pam in r["cas_pam_pairs"][:30]:
            print("    - %s -> %s" % (name, pam))
