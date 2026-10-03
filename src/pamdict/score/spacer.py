"""Spacer-to-protospacer evidence for candidate PAM ranking.

This module intentionally implements a small, auditable ungapped matcher.  It
is suitable for exact or low-Hamming-distance searches against a bounded
phage/plasmid FASTA.  Large databases and indel-tolerant discovery should use
BLAST or another aligner, then feed equivalent oriented flanks into the same
scoring semantics.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence

from .spectrum import allowed_set, normalize_pam, pam_information, valid_pam

DNA_BASES = frozenset("ACGT")
COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)
IUPAC_FROM_BASES = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("T"): "T",
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("CG"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
    frozenset("CGT"): "B",
    frozenset("AGT"): "D",
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}


def normalize_dna(sequence: str) -> str:
    return "".join((sequence or "").split()).upper().replace("U", "T")


def reverse_complement(sequence: str) -> str:
    sequence = normalize_dna(sequence)
    unknown = set(sequence) - set("ACGTRYSWKMBDHVN")
    if unknown:
        raise ValueError(f"unsupported DNA symbols: {sorted(unknown)}")
    return sequence.translate(COMPLEMENT)[::-1]


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    path = Path(path)
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, normalize_dna("".join(parts))))
                name = line[1:].split()[0] or f"record_{len(records) + 1}"
                parts = []
            elif name is None:
                raise ValueError(f"FASTA sequence precedes its header: {path}")
            else:
                parts.append(line)
    if name is not None:
        records.append((name, normalize_dna("".join(parts))))
    if not records or any(not sequence for _name, sequence in records):
        raise ValueError(f"FASTA has no complete records: {path}")
    return records


def find_hamming_matches(
    reference: str,
    query: str,
    max_mismatches: int = 0,
) -> list[tuple[int, int]]:
    """Return zero-based starts and Hamming distances for ungapped matches.

    Query partitioning into ``max_mismatches + 1`` exact seeds makes candidate
    discovery complete by the pigeonhole principle while keeping the code
    dependency-free.
    """
    reference = normalize_dna(reference)
    query = normalize_dna(query)
    if not query or set(query) - DNA_BASES:
        raise ValueError("query must contain only A/C/G/T")
    if max_mismatches < 0 or max_mismatches >= len(query):
        raise ValueError("max_mismatches must be in [0, len(query))")
    if len(query) > len(reference):
        return []

    candidates: set[int] = set()
    parts = max_mismatches + 1
    for part in range(parts):
        left = part * len(query) // parts
        right = (part + 1) * len(query) // parts
        seed = query[left:right]
        position = reference.find(seed)
        while position >= 0:
            start = position - left
            if 0 <= start <= len(reference) - len(query):
                candidates.add(start)
            position = reference.find(seed, position + 1)

    matches: list[tuple[int, int]] = []
    for start in sorted(candidates):
        target = reference[start:start + len(query)]
        mismatches = sum(
            left != right or right not in DNA_BASES
            for left, right in zip(query, target)
        )
        if mismatches <= max_mismatches:
            matches.append((start, mismatches))
    return matches


@dataclass(frozen=True)
class SpacerHit:
    spacer_id: str
    spacer_sequence: str
    target_id: str
    start_0: int
    end_0: int
    strand: str
    mismatches: int
    protospacer: str
    upstream_flank: str
    downstream_flank: str
    orientation_weight: float = 1.0

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def find_spacer_hits(
    spacers: Sequence[tuple[str, str]],
    targets: Sequence[tuple[str, str]],
    *,
    flank_length: int = 10,
    max_mismatches: int = 0,
) -> list[SpacerHit]:
    """Find both-strand spacer matches and return oriented target flanks."""
    if flank_length <= 0:
        raise ValueError("flank_length must be positive")

    unique_spacers: dict[str, str] = {}
    for spacer_id, raw_sequence in spacers:
        sequence = normalize_dna(raw_sequence)
        if not sequence or set(sequence) - DNA_BASES:
            raise ValueError(f"spacer {spacer_id!r} must contain only A/C/G/T")
        unique_spacers.setdefault(sequence, spacer_id)

    normalized_targets: list[tuple[str, str]] = []
    seen_target_ids: set[str] = set()
    for target_id, raw_sequence in targets:
        if target_id in seen_target_ids:
            raise ValueError(f"duplicate target FASTA id: {target_id}")
        seen_target_ids.add(target_id)
        sequence = normalize_dna(raw_sequence)
        if not sequence or set(sequence) - set("ACGTN"):
            raise ValueError(
                f"target {target_id!r} must contain only A/C/G/T/N"
            )
        normalized_targets.append((target_id, sequence))

    hits: list[SpacerHit] = []
    for spacer_sequence, spacer_id in unique_spacers.items():
        reverse = reverse_complement(spacer_sequence)
        palindromic = reverse == spacer_sequence
        orientations = [("+", spacer_sequence), ("-", reverse)]
        for target_id, target in normalized_targets:
            for strand, query in orientations:
                for start, mismatches in find_hamming_matches(
                    target, query, max_mismatches
                ):
                    end = start + len(query)
                    if strand == "+":
                        protospacer = target[start:end]
                        upstream = target[max(0, start - flank_length):start]
                        downstream = target[end:end + flank_length]
                    else:
                        protospacer = reverse_complement(target[start:end])
                        upstream = reverse_complement(
                            target[end:end + flank_length]
                        )
                        downstream = reverse_complement(
                            target[max(0, start - flank_length):start]
                        )
                    hits.append(SpacerHit(
                        spacer_id=spacer_id,
                        spacer_sequence=spacer_sequence,
                        target_id=target_id,
                        start_0=start,
                        end_0=end,
                        strand=strand,
                        mismatches=mismatches,
                        protospacer=protospacer,
                        upstream_flank=upstream,
                        downstream_flank=downstream,
                        orientation_weight=0.5 if palindromic else 1.0,
                    ))
    return sorted(
        hits,
        key=lambda hit: (
            hit.spacer_id,
            hit.target_id,
            hit.start_0,
            hit.strand,
        ),
    )


def pam_from_hit(hit: SpacerHit, side: str, length: int) -> str | None:
    if length <= 0:
        raise ValueError("PAM length must be positive")
    if side == "downstream":
        pam = hit.downstream_flank[:length]
    elif side == "upstream":
        pam = hit.upstream_flank[-length:]
    else:
        raise ValueError("side must be 'downstream' or 'upstream'")
    if len(pam) != length or set(pam) - DNA_BASES:
        return None
    return pam


def concrete_matches_iupac(sequence: str, motif: str) -> bool:
    motif = normalize_pam(motif)
    sequence = normalize_dna(sequence)
    return (
        len(sequence) == len(motif)
        and valid_pam(motif)
        and all(base in allowed_set(symbol) for base, symbol in zip(sequence, motif))
    )


def background_base_frequencies(
    targets: Sequence[tuple[str, str]],
) -> dict[str, float]:
    counts = Counter(
        base
        for _name, sequence in targets
        for base in normalize_dna(sequence)
        if base in DNA_BASES
    )
    total = sum(counts.values())
    if total == 0:
        raise ValueError("target FASTA contains no canonical DNA bases")
    return {base: counts[base] / total for base in "ACGT"}


def candidate_background_probability(
    motif: str,
    base_frequencies: dict[str, float],
) -> float:
    motif = normalize_pam(motif)
    if not valid_pam(motif):
        raise ValueError(f"invalid IUPAC PAM: {motif!r}")
    probability = 1.0
    for symbol in motif:
        probability *= sum(base_frequencies[base] for base in allowed_set(symbol))
    return probability


def _weighted_observations(
    hits: Sequence[SpacerHit],
    side: str,
    length: int,
) -> list[tuple[SpacerHit, str, float]]:
    grouped: dict[str, list[tuple[SpacerHit, str]]] = defaultdict(list)
    for hit in hits:
        pam = pam_from_hit(hit, side, length)
        if pam is not None:
            grouped[hit.spacer_sequence].append((hit, pam))

    weighted: list[tuple[SpacerHit, str, float]] = []
    for observations in grouped.values():
        denominator = sum(hit.orientation_weight for hit, _pam in observations)
        for hit, pam in observations:
            weighted.append((hit, pam, hit.orientation_weight / denominator))
    return weighted


@dataclass(frozen=True)
class SpacerPAMScore:
    candidate_pam: str
    pam_length: int
    spacer_rank_within_length: int | None
    spacer_empirical_probability: float | None
    spacer_posterior_probability: float | None
    spacer_background_probability: float
    spacer_evidence_score: float | None
    spacer_log2_enrichment: float | None
    effective_spacer_hits: int
    support_spacer_count: int
    support_target_count: int
    support_alignment_count: int
    eligible_target_count: int
    eligible_alignment_count: int
    candidate_information_bits_per_position: float

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        for key, value in list(result.items()):
            if isinstance(value, float):
                result[key] = round(value, 8)
        return result


def score_spacer_candidates(
    hits: Sequence[SpacerHit],
    candidates: Iterable[str],
    *,
    side: str,
    base_frequencies: dict[str, float],
    prior_strength: float = 2.0,
) -> list[SpacerPAMScore]:
    """Score candidates after giving each unique spacer total weight one."""
    if prior_strength < 0:
        raise ValueError("prior_strength must be non-negative")
    motifs = list(dict.fromkeys(normalize_pam(value) for value in candidates))
    if not motifs:
        raise ValueError("at least one candidate PAM is required")
    if any(not valid_pam(motif) for motif in motifs):
        invalid = [motif for motif in motifs if not valid_pam(motif)]
        raise ValueError(f"invalid IUPAC PAMs: {invalid}")

    observation_cache = {
        length: _weighted_observations(hits, side, length)
        for length in {len(motif) for motif in motifs}
    }
    scores: list[SpacerPAMScore] = []
    for motif in motifs:
        observations = observation_cache[len(motif)]
        total_weight = float(len({hit.spacer_sequence for hit, _pam, _w in observations}))
        background = candidate_background_probability(motif, base_frequencies)
        supported = [
            (hit, weight)
            for hit, pam, weight in observations
            if concrete_matches_iupac(pam, motif)
        ]
        support_weight = sum(weight for _hit, weight in supported)
        empirical = support_weight / total_weight if total_weight else None
        if total_weight and pam_information(motif) > 0:
            posterior = (
                support_weight + prior_strength * background
            ) / (total_weight + prior_strength)
            if posterior >= background:
                bounded = (
                    (posterior - background) / (1.0 - background)
                    if background < 1.0 else 0.0
                )
            else:
                bounded = (
                    (posterior - background) / background
                    if background > 0.0 else 0.0
                )
            evidence_score = 50.0 * (bounded + 1.0)
            enrichment = math.log2(
                max(posterior, 1e-12) / max(background, 1e-12)
            )
        else:
            posterior = None
            evidence_score = None
            enrichment = None
        scores.append(SpacerPAMScore(
            candidate_pam=motif,
            pam_length=len(motif),
            spacer_rank_within_length=None,
            spacer_empirical_probability=empirical,
            spacer_posterior_probability=posterior,
            spacer_background_probability=background,
            spacer_evidence_score=evidence_score,
            spacer_log2_enrichment=enrichment,
            effective_spacer_hits=int(total_weight),
            support_spacer_count=len({hit.spacer_sequence for hit, _w in supported}),
            support_target_count=len({hit.target_id for hit, _w in supported}),
            support_alignment_count=len(supported),
            eligible_target_count=len({hit.target_id for hit, _pam, _w in observations}),
            eligible_alignment_count=len(observations),
            candidate_information_bits_per_position=pam_information(motif),
        ))

    ranks: dict[str, int] = {}
    for length in sorted({score.pam_length for score in scores}):
        ranked = [
            score for score in scores
            if score.pam_length == length and score.spacer_evidence_score is not None
        ]
        ranked.sort(key=lambda score: (
            -float(score.spacer_evidence_score),
            -score.support_spacer_count,
            -score.support_target_count,
            score.candidate_pam,
        ))
        for rank, score in enumerate(ranked, start=1):
            ranks[score.candidate_pam] = rank
    return [
        replace(score, spacer_rank_within_length=ranks.get(score.candidate_pam))
        for score in scores
    ]


def flank_probability_matrix(
    hits: Sequence[SpacerHit],
    *,
    side: str,
    length: int,
) -> tuple[list[dict[str, float]], int]:
    observations = _weighted_observations(hits, side, length)
    effective_spacers = len({hit.spacer_sequence for hit, _pam, _w in observations})
    rows: list[dict[str, float]] = []
    for position in range(length):
        counts = {base: 0.0 for base in "ACGT"}
        for _hit, pam, weight in observations:
            counts[pam[position]] += weight
        total = sum(counts.values())
        rows.append({
            base: counts[base] / total if total else 0.25
            for base in "ACGT"
        })
    return rows, effective_spacers


def information_bits(row: dict[str, float]) -> float:
    return 2.0 + sum(
        probability * math.log2(probability)
        for probability in row.values()
        if probability > 0
    )


def probability_consensus(
    matrix: Sequence[dict[str, float]],
    *,
    information_floor: float = 0.15,
) -> str:
    symbols: list[str] = []
    for row in matrix:
        if information_bits(row) < information_floor:
            symbols.append("N")
            continue
        maximum = max(row.values())
        selected = frozenset(
            base for base, probability in row.items()
            if probability >= max(0.25, maximum * 0.5)
        )
        symbols.append(IUPAC_FROM_BASES.get(selected, "N"))
    return "".join(symbols)
