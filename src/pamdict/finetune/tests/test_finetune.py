from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import torch

WORKSPACE = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WORKSPACE))

from scripts.finetune_p2pam import (  # noqa: E402
    _epoch_batches,
    _trainable_state_dict,
    load_rows,
    logo_of,
    position_weights_of,
    training_fingerprint,
)
from pamdict.finetune.similarity import (  # noqa: E402
    PROTEIN_ALPHABET,
    ProteinNeighborIndex,
    banded_edit_distance,
    find_near_neighbors,
    global_edit_identity,
    protein_kmer_counts,
)


def test_logo_zero_information_is_uniform():
    target = logo_of({
        "logo_json": json.dumps([[0.0, 0.0, 0.0, 0.0] for _ in range(10)]),
        "pam_consensus": "",
    })
    assert tuple(target.shape) == (10, 4)
    assert torch.allclose(target, torch.full((10, 4), 0.25))


def test_load_rows_reports_corrupt_header():
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "bad.tsv"
        path.write_text("control-protein_sequence\tpam_consensus\nAAA\tNGG\n")
        try:
            load_rows(path, "cas9_full")
        except ValueError as exc:
            assert "protein_sequence" in str(exc)
            assert "parsed columns" in str(exc)
        else:
            raise AssertionError("corrupt training header was accepted")


def test_trainable_state_includes_classifier():
    model = torch.nn.Module()
    model.base = torch.nn.Linear(3, 3)
    model.classifier = torch.nn.Linear(3, 2)
    for parameter in model.base.parameters():
        parameter.requires_grad_(False)
    state = _trainable_state_dict(model)
    assert sorted(state) == ["classifier.bias", "classifier.weight"]


def test_epoch_batches_are_deterministic_and_complete():
    rows = [{"protein_sequence": "A" * length} for length in (8, 2, 6, 3, 7, 4)]
    first = _epoch_batches(rows, batch_size=2, seed=13)
    second = _epoch_batches(rows, batch_size=2, seed=13)
    assert first == second
    assert sorted(i for batch in first for i in batch) == list(range(len(rows)))
    for batch in first:
        lengths = [len(rows[i]["protein_sequence"]) for i in batch]
        assert max(lengths) - min(lengths) <= 2


def test_protein_kmers_use_canonical_amino_acid_alphabet():
    assert len(PROTEIN_ALPHABET) == 20
    counts = protein_kmer_counts("AEFILPQ")
    assert counts["AEF"] == 1
    assert counts["LPQ"] == 1


def test_banded_edit_identity_and_neighbor_confirmation():
    reference = "ACDEFGHIKLMNPQRSTVWY" * 5
    one_mutation = reference[:20] + "Y" + reference[21:]
    assert banded_edit_distance(reference, one_mutation, 1) == 1
    identity = global_edit_identity(reference, one_mutation, 0.90)
    assert identity is not None and identity >= 0.99
    matches = find_near_neighbors(
        [one_mutation], [reference, "Y" * len(reference)],
        identity_threshold=0.90, candidate_cosine=0.10,
    )
    assert [item["reference_sequence"] for item in matches[one_mutation]] == [reference]
    indexed = ProteinNeighborIndex([reference, "Y" * len(reference)]).find(
        [one_mutation],
        identity_threshold=0.90,
        candidate_cosine=0.10,
        probe_kmers=4,
    )
    assert [item["reference_sequence"] for item in indexed[one_mutation]] == [reference]


def test_load_rows_rejects_wrong_family():
    logo = json.dumps([[0.25, 0.25, 0.25, 0.25] for _ in range(10)])
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "wrong_family.tsv"
        path.write_text(
            "crispr_type\tprotein_sequence\tlogo_json\n"
            f"Type I\tAAAA\t{logo}\n",
            encoding="utf-8",
        )
        try:
            load_rows(path, "cas9_full")
        except ValueError as exc:
            assert "requires Type II" in str(exc)
        else:
            raise AssertionError("wrong-family row was accepted")


def test_position_weights_and_fingerprint_are_auditable():
    weights = position_weights_of({
        "position_weight_json": json.dumps([0.1 + i for i in range(10)])
    })
    assert tuple(weights.shape) == (10,)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "train.tsv"
        path.write_text(
            "crispr_type\tprotein_sequence\tlogo_json\n"
            "Type II\tAAAA\t[]\n",
            encoding="utf-8",
        )
        fingerprint = training_fingerprint(path, "cas9_full")
        assert fingerprint["rows"] == 1
        assert fingerprint["required_family"] == "Type II"
        assert fingerprint["sequence_column"] == "protein_sequence"
        assert len(fingerprint["sha256"]) == 64
