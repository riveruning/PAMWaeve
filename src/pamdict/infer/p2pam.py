"""Run the official Protein2PAM model (HF ESM2 implementation) locally.

This module loads the official ``EsmForSequenceClassification`` model from the
Profluent-AI/Protein2PAM repo (cloned under ``.reference/Protein2PAM``) and the
pretrained weights from the HuggingFace Hub, runs forward inference on protein
sequences, and converts logits -> probability -> information content (bits) ->
consensus PAM, reproducing the official conversion semantics exactly.

Why importlib: the official repo's top-level ``protein2pam/__init__.py`` pulls
in the full Oracle stack (Bio, Levenshtein, fair-esm, logomaker) which we do not
need for bare PAM prediction. We load only the ``huggingface`` submodule via
source files, registering ``protein2pam.huggingface.configuration_esm`` in
``sys.modules`` so the intra-package import inside ``modeling_esm.py`` resolves.

Environment
-----------
Run with the d2l-zh interpreter (torch 2.12) plus the workspace .pylibs on the
path, e.g.::

    export PYTHONPATH=.pylibs:.reference/Protein2PAM:src
    export HF_HOME=$PWD/data/checkpoints/hf
    python -m pamdict.infer.p2pam --model cas9 --input ... --output ...

Input format: cas9 expects the Cas9 *PI-domain* (pid_sequence), cas12 expects
the full-length protein. cas8 expects Cas8/Cas10d.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

# --- Resolve the official upstream adapter location ------------------------
# The Protein2PAM adapter tree is a *deployment* dependency, not a file of this
# project, so it is located in this order:
#   1. ``PAMPRIDICT_PROTEIN2PAM_HOME`` -- explicit override, pointing either at
#      the cloned repo root or directly at its ``protein2pam`` package;
#   2. an importable ``protein2pam`` package already on ``sys.path``;
#   3. the conventional in-checkout location ``.reference/Protein2PAM``.
# A hardcoded relative path alone would break any deployment that keeps the
# upstream tree elsewhere (see SUBMISSION.md verification notes).
_WORKSPACE = Path(__file__).resolve().parents[3]  # src/pamdict/infer/p2pam.py -> root


def _resolve_upstream_adapter() -> Path | None:
    override = os.environ.get("PAMPRIDICT_PROTEIN2PAM_HOME")
    if override:
        root = Path(override).expanduser()
        # Accept either the repo root or the inner package directory.
        for candidate in (
            root / "protein2pam" / "huggingface",
            root / "huggingface",
            root,
        ):
            if candidate.is_dir():
                return candidate
        return root
    for entry in sys.path:
        if not entry:
            continue
        candidate = Path(entry) / "protein2pam" / "huggingface"
        if candidate.is_dir():
            return candidate
    conventional = _WORKSPACE / ".reference" / "Protein2PAM"
    if conventional.is_dir():
        return conventional / "protein2pam" / "huggingface"
    return None


_REPO = _resolve_upstream_adapter()

# Model names available on the HF hub.
#   cas9      -> Cas9 PI-domain input (webserver model)
#   cas9_full -> full-length Cas9 input
#   cas12     -> full-length Cas12 input (webserver model)
#   cas8      -> Cas8/Cas10d input
#   *_nolit   -> ablation models (no literature-augmented training data)
MODEL_HF = {
    "cas8": "Profluent-Bio/protein2pam-cas8",
    "cas9": "Profluent-Bio/protein2pam-cas9",
    "cas9_full": "Profluent-Bio/protein2pam-cas9_full",
    "cas12": "Profluent-Bio/protein2pam-cas12",
    "cas9_full_nolit": "Profluent-Bio/protein2pam-cas9_full_nolit",
    "cas9_nolit": "Profluent-Bio/protein2pam-cas9_pid_nolit",
    "cas12_nolit": "Profluent-Bio/protein2pam-cas12_nolit",
}

# Local snapshot dirs (under HF_HOME) used when ``local_files_only`` is set or
# the network is down. The models are also mirrored here to survive HF outages.
def _local_snapshot_dir(hf_name: str) -> str:
    """Best-effort local snapshot path for a hub model (HF_HOME cache layout)."""
    hf_home = os.environ.get("HF_HOME") or os.path.expanduser("~/.cache/huggingface")
    slug = hf_name.replace("/", "--")
    base = Path(hf_home) / f"models--{slug}"
    snaps = base / "snapshots"
    if snaps.exists():
        snaps = sorted(p for p in snaps.iterdir() if p.is_dir())
        if snaps:
            return str(snaps[0])
    # New hub layout fallback.
    base2 = Path(hf_home) / "hub" / f"models--{slug}"
    snaps2 = base2 / "snapshots"
    if snaps2.exists():
        s2 = sorted(p for p in snaps2.iterdir() if p.is_dir())
        if s2:
            return str(s2[0])
    return hf_name  # fall back to hub name (requires network)

# Position threshold for consensus (bits), matching official PAM.consensus_pam.
CONSENSUS_THRESHOLD_BITS = 0.70
NUCLEOTIDES = "ACGT"


def _load_official_model_cls():
    """Import the official EsmForSequenceClassification + EsmConfig via importlib."""
    if _REPO is None:
        raise RuntimeError(
            "Protein2PAM upstream adapter not found. Provide it in one of "
            "these ways:\n"
            "  1. clone the upstream repo and set "
            "PAMPRIDICT_PROTEIN2PAM_HOME=/path/to/Protein2PAM\n"
            "  2. put a directory containing 'protein2pam/' on PYTHONPATH\n"
            "  3. place it at <package-root>/.reference/Protein2PAM\n"
            "See models/MODEL_CARD.md section 1 for the pinned commit.\n"
            "Searched override, sys.path, and "
            f"{_WORKSPACE / '.reference' / 'Protein2PAM'}"
        )

    def _load(name: str, path: str):
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return mod

    cfg = _load(
        "protein2pam.huggingface.configuration_esm",
        str(_REPO / "configuration_esm.py"),
    )
    model_mod = _load(
        "protein2pam.huggingface.modeling_esm",
        str(_REPO / "modeling_esm.py"),
    )
    return model_mod.EsmForSequenceClassification, cfg.EsmConfig


def prob_to_info(p: torch.Tensor) -> torch.Tensor:
    """Convert a probability matrix (batch, pos, 4) to information content (bits).

    Reproduces ``protein2pam.utils.common.prob_to_info`` exactly:
        info = p * (sum_k p_k * (log2 p_k + log2 K)) / log2(e) ... equivalently
    the official code computes, per position, the weighted information content.
    """
    logp = torch.log(p + 1e-12)
    K = p.shape[-1]
    # official: p * nansum(p * (logp + log K), dim=-1, keepdim) / log(2)
    inner = p * (logp + np.log(K))
    info = p * torch.nansum(inner, dim=-1, keepdims=True) / np.log(2)
    return info


def prob_to_info_numpy(p: np.ndarray) -> np.ndarray:
    """NumPy equivalent of prob_to_info for cached probabilities."""
    values = np.asarray(p, dtype=np.float32)
    k = values.shape[-1]
    inner = values * (np.log(values + 1e-12) + np.log(k))
    return values * np.nansum(inner, axis=-1, keepdims=True) / np.log(2)


def consensus_from_info(info: np.ndarray, side: str = "downstream",
                        min_info: float = CONSENSUS_THRESHOLD_BITS) -> str:
    """Build the consensus PAM from a (pos, 4) information matrix.

    Matches official ``PAM.consensus_pam``: a position is 'N' if its max info is
    below the threshold; otherwise the argmax nucleotide. Downstream (cas9)
    strips trailing N; upstream (cas12/cas8) strips leading N.
    """
    nts = []
    for row in info:
        if np.max(row) < min_info:
            nts.append("N")
        else:
            nts.append(NUCLEOTIDES[int(np.argmax(row))])
    s = "".join(nts)
    if side == "downstream":
        return s.rstrip("N")
    else:
        return s.lstrip("N")


class P2PAMPredictor:
    def __init__(self, model_name: str = "cas9", device: str = None,
                 local_files_only: bool = True):
        if model_name not in MODEL_HF:
            raise ValueError(f"unknown model {model_name!r}; choose {sorted(MODEL_HF)}")
        self.model_name = model_name
        self.hf_name = MODEL_HF[model_name]
        self.side = "downstream" if "cas9" in model_name else "upstream"
        self._device = device or ("cuda" if torch.cuda.is_available() else "cpu")

        EsmCls, _ = _load_official_model_cls()
        # Load tokenizer from the official repo (tokenizer.json).
        from tokenizers import Tokenizer
        self.tokenizer = Tokenizer.from_file(str(_REPO / "tokenizer.json"))

        # Prefer the local snapshot (survives HF outages); fall back to hub name.
        src = _local_snapshot_dir(self.hf_name)
        eprint(f"loading model weights from {src} (local_files_only={local_files_only}) ...")
        self.model = EsmCls.from_pretrained(
            src, local_files_only=True
        )
        self.model.to(self._device).eval()
        eprint("model loaded onto", self._device)

    def predict_probability_matrix(self, sequences) -> np.ndarray:
        """Return A/C/G/T probability matrices with shape (batch, 10, 4)."""
        if isinstance(sequences, str):
            sequences = [sequences]
        encs = self.tokenizer.encode_batch(sequences)
        input_ids = torch.tensor([e.ids for e in encs], device=self._device)
        attention_mask = torch.tensor(
            [e.attention_mask for e in encs], device=self._device
        )
        with torch.no_grad():
            logits = self.model(
                input_ids=input_ids, attention_mask=attention_mask
            ).logits
            # Always normalize in FP32. BF16 softmax can leave each four-base
            # row a few 1e-3 away from one after serialization, which is large
            # enough to perturb KL comparisons.
            prob = F.softmax(logits.float(), dim=-1)
        return prob.detach().float().cpu().numpy()

    def predict_matrix(self, sequences) -> np.ndarray:
        """Return information matrices (batch, 10, 4) in bits for the sequences."""
        probabilities = self.predict_probability_matrix(sequences)
        return prob_to_info_numpy(probabilities)

    def predict(self, sequences):
        """Return list of consensus PAM strings."""
        mats = self.predict_matrix(sequences)
        return [consensus_from_info(m, side=self.side) for m in mats]


def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="cas9", choices=sorted(MODEL_HF))
    ap.add_argument("--device", default=None)
    ap.add_argument("--input", help="TSV with a 'sequence' column (one per row)")
    ap.add_argument("--seq-col", default="sequence")
    ap.add_argument("--output", default="-")
    ap.add_argument("--matrix", action="store_true",
                    help="also write the full 10x4 info matrix (JSON-encoded)")
    args = ap.parse_args()

    pred = P2PAMPredictor(args.model, args.device)

    if args.input and args.input != "-":
        rows = []
        with open(args.input, encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows.append(r)
        seqs = [r[args.seq_col] for r in rows]
    else:
        seqs = [line.strip() for line in sys.stdin if line.strip()]
        rows = [{"sequence": s} for s in seqs]

    mats = pred.predict_matrix(seqs)
    pams = [consensus_from_info(m, side=pred.side) for m in mats]

    out_fields = ("sequence", "consensus_pam")
    if rows and any("cas_name" in r for r in rows):
        out_fields = ("cas_name", "consensus_pam")

    if args.output == "-":
        w = csv.DictWriter(sys.stdout, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        for r, p in zip(rows, pams):
            w.writerow({k: r.get(k, "") if k != "consensus_pam" else p for k in out_fields})
    else:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t")
            w.writeheader()
            for r, p in zip(rows, pams):
                w.writerow({k: r.get(k, "") if k != "consensus_pam" else p for k in out_fields})

    if args.matrix:
        np.save("_p2pam_matrices.npy", mats)
        eprint(f"matrices saved to _p2pam_matrices.npy shape={mats.shape}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
