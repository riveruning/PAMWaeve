"""Integration tests for the genome agent that need the real model weights.

These are deliberately kept OUT of ``scripts/run_tests.py`` so that the offline
unit suite stays runnable on a clean machine with no network and no checkpoints.
They are run separately:

    HF_HOME=$PWD/data/checkpoints/hf PYTHONPATH=src:$PYTHONPATH \\
    python3 scripts/run_integration_tests.py

Every test here calls :func:`_require_model` first, so a missing checkpoint
produces one clear message instead of a transformers stack trace. Nothing here
silently skips: if the prerequisites are absent the run FAILS and says why.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
MODEL_REPO = "Profluent-Bio/protein2pam-cas9_full"
# The verified real Cas9 whose consensus is NGG. Tests must not depend on
# historical run output, so this is regenerated from the input genome instead.
GENOME = (
    ROOT / "data/raw/external_datasets/ncbi_cas9_candidate_genomes/ncbi_dataset"
    / "data/GCF_000767505.1/GCF_000767505.1_ASM76750v1_genomic.fna"
)
REAL_CAS9 = (
    ROOT / "data/parsed/genome_agent_real_20260928/final_single/selected_cas9.faa"
)


class PrerequisiteMissing(RuntimeError):
    """Raised when the integration environment is not usable."""


def _hf_home() -> Path:
    return Path(os.environ.get("HF_HOME") or (ROOT / "data/checkpoints/hf"))


def _model_snapshot() -> Path | None:
    hub = _hf_home() / "hub" / f"models--{MODEL_REPO.replace('/', '--')}" / "snapshots"
    if not hub.is_dir():
        return None
    for snap in sorted(hub.iterdir()):
        if (snap / "config.json").is_file():
            return snap
    return None


def _require_model() -> Path:
    """Fail loudly with an actionable message when the model is unavailable."""
    snap = _model_snapshot()
    if snap is None:
        raise PrerequisiteMissing(
            "integration test needs the protein2pam model weights.\n"
            f"  looked in : {_hf_home() / 'hub'}\n"
            "  fix       : export HF_HOME=$PWD/data/checkpoints/hf "
            "(and make sure the snapshot is present), then re-run "
            "scripts/run_integration_tests.py"
        )
    # transformers caches are addressed via HF_HOME, so make it explicit for the
    # child process regardless of how the caller set it.
    os.environ["HF_HOME"] = str(_hf_home())
    return snap


def _require_real_cas9() -> Path:
    if not REAL_CAS9.is_file():
        raise PrerequisiteMissing(
            f"integration test needs the verified real Cas9 fixture: {REAL_CAS9}\n"
            "  fix: run scripts/run_genome_agent.py on "
            f"{GENOME.relative_to(ROOT)} to regenerate it"
        )
    return REAL_CAS9


def _python() -> str:
    return sys.executable


def test_integration_ngg_auto_candidates_reach_spacer_stage():
    """The real CLI path: --auto-candidates on a real Cas9 must yield a PAM.

    This is the exact protein->spacer hand-off that was once broken: a consensus
    containing N was dropped, so the spacer stage never received a candidate.
    """
    _require_model()
    cas9 = _require_real_cas9()
    with tempfile.TemporaryDirectory() as d:
        result = subprocess.run(
            [_python(), str(ROOT / "scripts/score_candidate_pams.py"),
             "--protein", str(cas9), "--model", "cas9_full",
             "--auto-candidates", "--device", "cpu", "--train-fasta", "",
             "--out", str(Path(d) / "o.tsv"), "--json-out", str(Path(d) / "o.json")],
            capture_output=True, text=True, timeout=900, cwd=str(ROOT))
        assert result.returncode == 0, result.stderr[-2000:]
        payload = json.loads((Path(d) / "o.json").read_text())
        entry = payload["proteins"][0]
        assert entry["predicted_pam"] == "NGG", entry["predicted_pam"]
        assert entry["determined_positions"] == [2, 3], entry["determined_positions"]
        assert not entry["no_clear_pam"]
        motifs = [r["candidate_pam"] for r in payload["scores"]]
        assert motifs == ["NGG"], f"spacer stage would receive {motifs}"
        assert len(motifs[0]) == 3


def test_integration_determined_positions_match_consensus_rule():
    """A real divergent Cas9 must report only positions the 0.70 rule called."""
    _require_model()
    from pamdict.infer.p2pam import CONSENSUS_THRESHOLD_BITS
    divergent = (
        ROOT / "data/parsed/genome_agent_real_20260928/final_override/selected_cas9.faa"
    )
    if not divergent.is_file():
        raise PrerequisiteMissing(
            f"integration test needs the divergent Cas9 fixture: {divergent}"
        )
    with tempfile.TemporaryDirectory() as d:
        result = subprocess.run(
            [_python(), str(ROOT / "scripts/score_candidate_pams.py"),
             "--protein", str(divergent), "--model", "cas9_full",
             "--auto-candidates", "--device", "cpu", "--train-fasta", "",
             "--out", str(Path(d) / "o.tsv"), "--json-out", str(Path(d) / "o.json")],
            capture_output=True, text=True, timeout=900, cwd=str(ROOT))
        assert result.returncode == 0, result.stderr[-2000:]
        entry = json.loads((Path(d) / "o.json").read_text())["proteins"][0]
        # Recompute the answer independently from the published matrix.
        matrix = entry["probability_matrix"]
        called = [
            i + 1 for i, row in enumerate(matrix[: len(entry["predicted_pam"])])
            if max(row) >= CONSENSUS_THRESHOLD_BITS
        ]
        assert entry["determined_positions"] == called, (
            entry["predicted_pam"], entry["determined_positions"], called)
        # The historic bug reported every information-signal position, which for
        # this protein was [2,3,4,5,6,7] while only [5,7] were real calls.
        assert entry["determined_positions"] != [2, 3, 4, 5, 6, 7]
