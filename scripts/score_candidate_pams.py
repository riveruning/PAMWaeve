"""Score candidate PAMs for one or more Cas proteins without retraining."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import torch

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE / "src"))

from pamdict.finetune.similarity import ProteinNeighborIndex  # noqa: E402
from pamdict.infer.p2pam import (  # noqa: E402
    CONSENSUS_THRESHOLD_BITS,
    MODEL_HF,
    P2PAMPredictor,
    consensus_from_info,
    prob_to_info_numpy,
)
from pamdict.score.candidate import (  # noqa: E402
    probability_information_stats,
    rank_candidate_pams,
)
from pamdict.score.spectrum import sequence_sha256  # noqa: E402


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence).upper()))
                name = line[1:].split()[0] or f"sequence_{len(records) + 1}"
                sequence = []
            elif name is None:
                raise ValueError(f"FASTA sequence appears before a header: {path}")
            else:
                sequence.append("".join(line.split()))
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    if not records or any(not sequence for _name, sequence in records):
        raise ValueError(f"FASTA contains no complete protein records: {path}")
    return records


def load_candidates(inline: str | None, path: Path | None) -> list[str]:
    candidates: list[str] = []
    if inline:
        candidates.extend(value.strip() for value in inline.split(","))
    if path:
        lines = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if lines and "	" in lines[0]:
            reader = csv.DictReader(lines, delimiter="	")
            if reader.fieldnames is None:
                raise ValueError(f"candidate file has no header: {path}")
            column = next(
                (
                    value for value in ("candidate_pam", "pam", "pam_consensus")
                    if value in reader.fieldnames
                ),
                None,
            )
            if column is None:
                raise ValueError(
                    "candidate TSV needs candidate_pam, pam, or pam_consensus"
                )
            candidates.extend((row.get(column) or "").strip() for row in reader)
        else:
            candidates.extend(lines)
    candidates = list(dict.fromkeys(value.upper() for value in candidates if value))
    if not candidates:
        raise ValueError("provide candidates with --candidates or --candidates-file")
    return candidates


def neighbor_stats(
    proteins: list[str],
    train_fasta: Path | None,
    *,
    minimum_identity: float,
) -> dict[str, dict[str, object]]:
    if train_fasta is None:
        return {
            sequence: {
                "nearest_train_identity": None,
                "mean_top10_train_identity": None,
                "training_exact_match": None,
                "neighbor_note": "training FASTA not supplied",
            }
            for sequence in proteins
        }
    references = [sequence for _name, sequence in read_fasta(train_fasta)]
    index = ProteinNeighborIndex(references, k=3)
    matches = index.find(
        proteins,
        identity_threshold=minimum_identity,
        candidate_cosine=0.05,
        max_candidates=256,
        probe_kmers=None,
    )
    results: dict[str, dict[str, object]] = {}
    for sequence in proteins:
        values = matches.get(sequence.upper(), [])
        identities = [float(value["identity"]) for value in values[:10]]
        results[sequence] = {
            "nearest_train_identity": identities[0] if identities else None,
            "mean_top10_train_identity": (
                sum(identities) / len(identities) if identities else None
            ),
            "training_exact_match": bool(identities and identities[0] == 1.0),
            "neighbor_note": (
                f"confirmed among indexed candidates at identity >= "
                f"{minimum_identity:.2f}"
                if identities else
                f"no confirmed training neighbor at identity >= "
                f"{minimum_identity:.2f}"
            ),
        }
    return results


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Rank IUPAC PAM candidates by compatibility with a Protein2PAM "
            "probability matrix. Scores are not cleavage-activity estimates."
        )
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--protein", help="protein FASTA")
    source.add_argument(
        "--sequence",
        action="append",
        help="literal protein sequence; may be supplied more than once",
    )
    parser.add_argument("--candidates", help="comma-separated IUPAC PAMs")
    parser.add_argument("--candidates-file", help="one PAM per line or TSV")
    parser.add_argument("--auto-candidates", action="store_true", help="score the model consensus only; not exhaustive PAM discovery")
    parser.add_argument("--model", default="cas9_full", choices=sorted(MODEL_HF))
    parser.add_argument("--device", default=None)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument(
        "--train-fasta",
        default="auto",
        help=(
            "training proteins for sequence-domain diagnostics; 'auto' uses "
            "data/raw/cas9_full.fasta only with cas9_full; use '' to skip"
        ),
    )
    parser.add_argument(
        "--neighbor-min-identity",
        type=float,
        default=0.30,
        help="minimum exact global identity reported by the candidate search",
    )
    parser.add_argument("--out", default="-", help="output TSV or '-'")
    parser.add_argument("--json-out", default=None)
    args = parser.parse_args()
    if args.batch_size <= 0:
        raise ValueError("--batch-size must be positive")
    if not 0 < args.neighbor_min_identity <= 1:
        raise ValueError("--neighbor-min-identity must be in (0, 1]")

    candidate_path = (
        WORKSPACE / args.candidates_file if args.candidates_file else None
    )
    if args.protein:
        records = read_fasta(WORKSPACE / args.protein)
    else:
        records = [
            (f"sequence_{index}", "".join(value.split()).upper())
            for index, value in enumerate(args.sequence, start=1)
        ]
        if any(not sequence for _name, sequence in records):
            raise ValueError("--sequence cannot be empty")
    if args.auto_candidates and (args.candidates or args.candidates_file):
        parser.error('--auto-candidates cannot be combined with a candidate list')
    candidates = [] if args.auto_candidates else load_candidates(args.candidates, candidate_path)
    names = [name for name, _sequence in records]
    proteins = [sequence for _name, sequence in records]

    device = torch.device(
        args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    )
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but torch.cuda.is_available() is false")
    predictor = P2PAMPredictor(args.model, device="cpu", local_files_only=True)
    use_bf16 = device.type == "cuda" and torch.cuda.is_bf16_supported()
    if use_bf16:
        predictor.model.to(dtype=torch.bfloat16)
    predictor.model.to(device).eval()
    predictor._device = device
    print(
        f"inference device={device} bf16={use_bf16} proteins={len(proteins)}",
        file=sys.stderr,
        flush=True,
    )

    probabilities: list[np.ndarray] = []
    for start in range(0, len(proteins), args.batch_size):
        probabilities.extend(
            predictor.predict_probability_matrix(
                proteins[start:start + args.batch_size]
            )
        )

    train_fasta = None
    train_fasta_argument = args.train_fasta
    if train_fasta_argument == "auto":
        if args.model == "cas9_full":
            train_fasta_argument = "data/raw/cas9_full.fasta"
        else:
            train_fasta_argument = ""
            print(
                "training-neighbor diagnostics skipped: no automatic "
                f"reference is defined for {args.model}",
                file=sys.stderr,
            )
    if train_fasta_argument:
        train_fasta = WORKSPACE / train_fasta_argument
        if not train_fasta.exists():
            print(
                f"warning: training FASTA not found; skipping neighbors: "
                f"{train_fasta}",
                file=sys.stderr,
            )
            train_fasta = None
    neighbors = neighbor_stats(
        proteins,
        train_fasta,
        minimum_identity=args.neighbor_min_identity,
    )

    fieldnames = [
        "protein_id",
        "sequence_sha256",
        "model",
        "pam_side",
        "predicted_pam",
        "candidate_pam",
        "pam_length",
        "rank_within_length",
        "informative_positions",
        "specificity_adjusted_score",
        "specificity_adjusted_evidence_bits",
        "allowed_probability_geomean",
        "candidate_information_bits_per_position",
        "mean_information_bits_signal_positions",
        "mean_information_bits_all_positions",
        "signal_positions",
        "nearest_train_identity",
        "mean_top10_train_identity",
        "training_exact_match",
        "neighbor_note",
        "interpretation",
    ]
    rows: list[dict[str, object]] = []
    protein_summaries: list[dict[str, object]] = []
    for name, sequence, probability in zip(names, proteins, probabilities):
        info_matrix = prob_to_info_numpy(probability)
        predicted_pam = consensus_from_info(info_matrix, side=predictor.side)
        information = probability_information_stats(probability)
        neighbor = neighbors[sequence]
        # A consensus containing 'N' is still a scorable IUPAC pattern: 'N' means
        # "any base" (see IUPAC in pamdict.score.spectrum), and a canonical 3-nt
        # Cas9 consensus such as NGG has an undetermined position 1 yet remains
        # the expected SpCas9 result.  Dropping such a consensus would break the
        # downstream spacer stage, which needs a candidate length to work with.
        # Only a consensus with no determined position at all is undecidable, and
        # that case is what must not be turned into a candidate.
        unresolved = not predicted_pam or not predicted_pam.replace("N", "")
        # A position is *called* only when the consensus rule accepted it, i.e.
        # its max base probability reached the threshold.  This is deliberately
        # the same rule consensus_from_info used, not the separate 0.15-bit
        # information diagnostic in `information["signal_positions"]`.
        call_threshold = CONSENSUS_THRESHOLD_BITS
        determined_positions = [
            i + 1 for i, row in enumerate(info_matrix[:len(predicted_pam)])
            if float(max(row)) >= call_threshold
        ]
        partial = not unresolved and len(determined_positions) < len(predicted_pam)
        if args.auto_candidates and unresolved:
            scored = []
        else:
            scored = rank_candidate_pams(
                probability, candidates or [predicted_pam], side=predictor.side
            )
        protein_summaries.append({
            "protein_id": name,
            "sequence_sha256": sequence_sha256(sequence),
            "predicted_pam": predicted_pam or None,
            "no_clear_pam": bool(unresolved),
            "partial_consensus": partial,
            "determined_positions": determined_positions,
            "no_clear_pam_note": (
                "The model determined no position for this protein; this is 'no clear PAM "
                "obtained', not a biological PAM result."
                if unresolved else
                "Some positions of the consensus are 'N' (any base), so it denotes a set of "
                "PAMs rather than one sequence. It is still scored as an IUPAC pattern; "
                "determined_positions lists the positions with a concrete base call."
                if partial else None
            ),
            "candidate_source": (
                "user_candidates" if candidates else
                "model_consensus" if not unresolved else "none"
            ),
            "probability_matrix": probability.tolist(),
            "information": information,
            "neighbor": neighbor,
        })
        for score in scored:
            item = score.to_dict()
            rows.append({
                "protein_id": name,
                "sequence_sha256": sequence_sha256(sequence),
                "model": args.model,
                "pam_side": predictor.side,
                "predicted_pam": predicted_pam,
                **item,
                "mean_information_bits_signal_positions": round(
                    float(information["mean_information_bits_signal_positions"]),
                    6,
                ),
                "mean_information_bits_all_positions": round(
                    float(information["mean_information_bits_all_positions"]),
                    6,
                ),
                "signal_positions": ",".join(
                    str(value) for value in information["signal_positions"]
                ),
                **neighbor,
                "interpretation": (
                    "no clear PAM obtained; not a biological PAM result"
                    if args.auto_candidates and unresolved else
                    "consensus scored as an IUPAC pattern; N positions mean any base "
                    "and are not a single-base call"
                    if args.auto_candidates and partial else
                    "model compatibility; not cleavage activity or experimental probability"
                ),
            })

    if args.out == "-":
        handle = sys.stdout
        close_handle = False
    else:
        out_path = WORKSPACE / args.out
        out_path.parent.mkdir(parents=True, exist_ok=True)
        handle = out_path.open("w", encoding="utf-8", newline="")
        close_handle = True
    try:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="	", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if close_handle:
            handle.close()

    if args.json_out:
        json_path = WORKSPACE / args.json_out
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(
            json.dumps({
                "format_version": 1,
                "model": args.model,
                "pam_side": predictor.side,
                "candidate_score_definition": (
                    "IUPAC probability mass relative to a uniform four-base background"
                ),
                "interpretation": (
                    "compatibility score only; not cleavage activity"
                ),
                "candidate_source": (
                    "user_candidates" if candidates else
                    "model_consensus_only" if any(
                        s.get("candidate_source") == "model_consensus"
                        for s in protein_summaries
                    ) else "none"
                ),
                "candidate_scope_note": (
                    "Only the supplied or consensus candidates were scored. This is a "
                    "presentation of model outputs, not exhaustive PAM discovery and not "
                    "independent validation."
                ),
                "no_clear_pam": all(
                    s.get("no_clear_pam") for s in protein_summaries
                ),
                "partial_consensus": [
                    s.get("predicted_pam") for s in protein_summaries
                    if s.get("partial_consensus")
                ],
                "proteins": protein_summaries,
                "scores": rows,
            }, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
