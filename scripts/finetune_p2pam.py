"""Fine-tune the official Protein2PAM model (ESM2-650M class) with LoRA.

Designed for the user's 8GB GPU (RTX 4060, WSL):
  - bf16 weights + LoRA adapters (attention q/k/v + classifier head trained).
    On an 8GB GPU, verify one optimizer step with micro-batch 1 before a full
    run; full-parameter training of a 650M model is not supported here.
  - 10x4 soft-target loss: KL(softmax(logits) || experimental/atlas logo) per
    position — treats PAM prediction as a 10-position x 4-base distribution
    matching problem (docs/PAM_DATABASE_LANDSCAPE.md route 3: spectrum-level
    target; fixes "variant outputs only main consensus").
  - Resume: restores trainable weights, optimizer/RNG state, and row position
    from ``last_checkpoint.pt``. A JSON state file alone is never resumed.

Usage (on the GPU machine, from the workspace root):
  export PYTHONPATH="$PWD:$PWD/.pylibs:$PWD/.reference/Protein2PAM:$PWD/src"
  export HF_HOME=$PWD/data/checkpoints/hf
  python scripts/finetune_p2pam.py --train data/corpus/augmented_train_v2.tsv \
      --model cas9_full --epochs 1 --batch-size 4 --grad-accum 8 --device cuda
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
import time
from pathlib import Path

import torch
import torch.nn.functional as F

WORKSPACE = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", default="data/corpus/augmented_train_v2_cas9.tsv")
    ap.add_argument("--model", default="cas9_full",
                    choices=["cas9_full", "cas9", "cas12"])
    ap.add_argument("--epochs", type=int, default=3)
    ap.add_argument("--batch-size", type=int, default=4,
                    help="micro-batch size; 4 is verified on this 8GB GPU")
    ap.add_argument("--grad-accum", type=int, default=8)
    ap.add_argument("--lr", type=float, default=1e-4, help="LoRA lr")
    ap.add_argument("--classifier-mode", choices=["freeze", "train"], default="freeze",
                    help="freeze is the conservative F0 control; train enables F1")
    ap.add_argument("--classifier-lr", type=float, default=1e-5,
                    help="classifier LR when --classifier-mode train")
    ap.add_argument("--loss-mode", choices=["equal", "info-weighted"],
                    default="equal")
    ap.add_argument("--label-smoothing", type=float, default=0.0,
                    help="research ablation only; keep 0 for family-fix runs")
    ap.add_argument("--lora-rank", type=int, default=8)
    ap.add_argument("--lora-alpha", type=int, default=16)
    ap.add_argument("--device", default=None, help="cuda / cpu (auto)")
    ap.add_argument("--max-steps", type=int, default=0, help="0 = full epochs")
    ap.add_argument("--out", default="data/checkpoints/finetune_v3_cas9")
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--log-every", type=int, default=1,
                    help="print an optimizer-step summary every N steps")
    ap.add_argument("--save-every", type=int, default=20,
                    help="write a resumable checkpoint every N optimizer steps")
    ap.add_argument("--no-resume", action="store_true",
                    help="ignore last_checkpoint.pt and start from the base model")
    ap.add_argument("--keep-step-checkpoints", action="store_true",
                    help="also retain immutable checkpoint_stepNNNNNN.pt files")
    return ap.parse_args()


# ---------------------------------------------------------------- LoRA layer
class LoRALinear(torch.nn.Module):
    """Wraps nn.Linear with a low-rank residual (w + B@A), A/B trainable."""

    def __init__(self, base: torch.nn.Linear, rank: int, alpha: int):
        super().__init__()
        self.base = base
        for p in self.base.parameters():
            p.requires_grad_(False)
        r = rank
        self.lora_a = torch.nn.Parameter(torch.empty(r, base.in_features))
        self.lora_b = torch.nn.Parameter(torch.zeros(base.out_features, r))
        torch.nn.init.kaiming_uniform_(self.lora_a, a=math.sqrt(5))
        self.scale = alpha / r
        self.dropout = torch.nn.Dropout(0.05)

    def forward(self, x):
        return self.base(x) + self.scale * (self.dropout(x) @ self.lora_a.T @ self.lora_b.T)


def apply_lora(model, target_names=("query", "key", "value"), rank=8, alpha=16):
    n = 0
    for name, module in model.named_modules():
        for tn in target_names:
            if name.endswith(tn) and isinstance(module, torch.nn.Linear):
                parent = model.get_submodule(name.rsplit(".", 1)[0]) if "." in name else model
                setattr(parent, tn, LoRALinear(module, rank, alpha))
                n += 1
    return n


# ---------------------------------------------------------------- data
MODEL_DATA_SPEC = {
    "cas9_full": ("Type II", "protein_sequence"),
    "cas9": ("Type II", "pid_sequence"),
    "cas12": ("Type V", "protein_sequence"),
}


def _normalize_crispr_type(value: str) -> str:
    compact = value.replace(" ", "").upper()
    return {"II": "Type II", "TYPEII": "Type II",
            "V": "Type V", "TYPEV": "Type V"}.get(compact, value.strip())


def training_fingerprint(path: Path, model_name: str) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        schema = [name.strip() for name in next(reader)]
        row_count = sum(1 for row in reader if row)
    return {
        "sha256": digest.hexdigest(),
        "rows": row_count,
        "model_name": model_name,
        "required_family": MODEL_DATA_SPEC[model_name][0],
        "sequence_column": MODEL_DATA_SPEC[model_name][1],
        "schema": schema,
    }


def load_rows(path: Path, family: str) -> list[dict]:
    required_type, sequence_column = MODEL_DATA_SPEC[family]
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"training TSV has no header: {path}")
        # Generated smoke files may contain a BOM or surrounding whitespace.
        # Normalize those harmless cases, but report the actual parsed header
        # when shell control bytes or a wrong delimiter corrupted the file.
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        missing = {sequence_column, "crispr_type"}.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                f"training TSV is missing columns {sorted(missing)}; "
                f"parsed columns are {reader.fieldnames!r}: {path}"
            )
        for line_no, r in enumerate(reader, start=2):
            actual_type = _normalize_crispr_type(r.get("crispr_type") or "")
            if actual_type != required_type:
                raise ValueError(
                    f"{family} requires {required_type}, got "
                    f"{r.get('crispr_type')!r} at {path}:{line_no}"
                )
            sequence = (r.get(sequence_column) or "").strip()
            if not sequence:
                raise ValueError(f"empty {sequence_column} at {path}:{line_no}")
            if not any((r.get(k) or "").strip()
                       for k in ("logo_json", "pam_logo_acgt", "pam_consensus")):
                raise ValueError(f"missing PAM target at {path}:{line_no}")
            logo_json = (r.get("logo_json") or "").strip()
            if logo_json:
                try:
                    raw_target = torch.tensor(json.loads(logo_json), dtype=torch.float32)
                except (json.JSONDecodeError, TypeError, ValueError) as exc:
                    raise ValueError(f"invalid logo_json at {path}:{line_no}") from exc
                if (
                    raw_target.shape != (10, 4)
                    or not torch.isfinite(raw_target).all()
                    or not (raw_target >= 0).all()
                    or not (raw_target.sum(-1) > 0).all()
                ):
                    raise ValueError(
                        f"logo_json must be a finite non-negative 10x4 matrix "
                        f"with positive row sums at {path}:{line_no}"
                    )
            r[sequence_column] = sequence
            r["_sequence"] = sequence
            target = logo_of(r)
            if (
                target.shape != (10, 4)
                or not torch.isfinite(target).all()
                or not torch.allclose(
                    target.sum(-1), torch.ones(10), rtol=0, atol=1e-5
                )
            ):
                raise ValueError(f"invalid PAM target at {path}:{line_no}")
            position_weights_of(r)
            rows.append(r)
    if not rows:
        raise ValueError(f"training TSV contains no data rows: {path}")
    return rows


def assert_no_benchmark_leakage(
    rows: list[dict], benchmark_path: Path, sequence_column: str
) -> None:
    required_type = _normalize_crispr_type(rows[0]["crispr_type"])
    with benchmark_path.open(encoding="utf-8", newline="") as handle:
        benchmark_sequences = {
            (row.get(sequence_column) or "").strip()
            for row in csv.DictReader(handle, delimiter="\t")
            if _normalize_crispr_type(row.get("crispr_type") or "") == required_type
        }
    overlap = {row["_sequence"] for row in rows} & benchmark_sequences
    if overlap:
        raise ValueError(
            f"benchmark leakage: {len(overlap)} exact protein sequences from "
            f"{benchmark_path} are present in training data"
        )


def logo_of(r: dict) -> torch.Tensor:
    """Return a valid 10x4 probability target.

    The upstream ``pam_logo_acgt`` values are information contributions rather
    than normalized probabilities.  Dividing a non-empty row by its sum
    recovers the base distribution; a zero-information row is unconstrained
    and therefore maps to the uniform distribution.
    """
    raw = (r.get("logo_json") or r.get("pam_logo_acgt") or "").strip()
    if raw.startswith("["):
        try:
            arr = json.loads(raw)
            if isinstance(arr, list) and len(arr) == 10 and isinstance(arr[0], list):
                t = torch.tensor(arr, dtype=torch.float32)
                if t.shape == (10, 4) and torch.isfinite(t).all() and (t >= 0).all():
                    s = t.sum(-1, keepdim=True)
                    uniform = torch.full_like(t, 0.25)
                    return torch.where(s > 0, t / s.clamp(min=1e-9), uniform)
        except (json.JSONDecodeError, TypeError, ValueError):
            pass
    # fallback: uniform-over-allowed-set from the consensus string
    IUPAC = {
        "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
        "R": ["A", "G"], "Y": ["C", "T"], "S": ["C", "G"], "W": ["A", "T"],
        "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"], "D": ["A", "G", "T"],
        "H": ["A", "C", "T"], "V": ["A", "C", "G"], "N": ["A", "C", "G", "T"],
    }
    pam = (r.get("pam_consensus") or "").upper()
    logo = []
    for i in range(10):
        ch = pam[i] if i < len(pam) else "N"
        allowed = IUPAC.get(ch, ["A", "C", "G", "T"])
        p = 1.0 / len(allowed)
        logo.append([p if n in allowed else 0.0 for n in "ACGT"])
    return torch.tensor(logo, dtype=torch.float32)


def position_weights_of(row: dict) -> torch.Tensor:
    raw = (row.get("position_weight_json") or "").strip()
    if not raw:
        return torch.ones(10, dtype=torch.float32)
    try:
        values = torch.tensor(json.loads(raw), dtype=torch.float32)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        raise ValueError("invalid position_weight_json") from exc
    if values.shape != (10,) or not torch.isfinite(values).all() or not (values > 0).all():
        raise ValueError("position_weight_json must contain 10 finite positive values")
    return values


def _trainable_state_dict(model: torch.nn.Module) -> dict[str, torch.Tensor]:
    """Copy every trainable parameter (LoRA *and* classifier) to CPU."""
    state = model.state_dict()
    names = {name for name, p in model.named_parameters() if p.requires_grad}
    return {name: state[name].detach().cpu() for name in sorted(names)}


def _atomic_torch_save(payload: dict, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, tmp)
    tmp.replace(path)


def _write_state_json(state: dict, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(state, indent=2), encoding="utf-8")
    tmp.replace(path)


def _move_optimizer_state(opt: torch.optim.Optimizer, device: torch.device) -> None:
    for values in opt.state.values():
        for key, value in values.items():
            if torch.is_tensor(value):
                values[key] = value.to(device)


def _epoch_batches(rows: list[dict], batch_size: int, seed: int) -> list[list[int]]:
    """Build deterministic length-bucketed batches to limit padding memory."""
    rng = random.Random(seed)
    order = list(range(len(rows)))
    rng.shuffle(order)
    order.sort(key=lambda i: len(rows[i].get("_sequence") or rows[i]["protein_sequence"]))
    batches = [order[i:i + batch_size] for i in range(0, len(order), batch_size)]
    rng.shuffle(batches)
    return batches


def main() -> int:
    args = parse_args()
    for name in ("epochs", "batch_size", "grad_accum", "lora_rank",
                 "lora_alpha", "log_every", "save_every"):
        if getattr(args, name) <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be positive")
    if args.lr <= 0 or args.classifier_lr <= 0:
        raise ValueError("learning rates must be positive")
    if not 0 <= args.label_smoothing < 1:
        raise ValueError("--label-smoothing must be in [0, 1)")
    random.seed(args.seed)
    torch.manual_seed(args.seed)
    device = torch.device(args.device or ("cuda" if torch.cuda.is_available() else "cpu"))
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but torch.cuda.is_available() is false")

    sys.path.insert(0, str(WORKSPACE))
    sys.path.insert(0, str(WORKSPACE / "src"))
    from pamdict.infer.p2pam import P2PAMPredictor  # reuse the official loader

    train_path = WORKSPACE / args.train
    fingerprint = training_fingerprint(train_path, args.model)
    rows = load_rows(train_path, args.model)
    benchmark_path = WORKSPACE / "data/corpus/gold_benchmark_v2.tsv"
    assert_no_benchmark_leakage(
        rows, benchmark_path, fingerprint["sequence_column"]
    )
    print(
        f"train rows: {len(rows)}  family={fingerprint['required_family']} "
        f"sequence_column={fingerprint['sequence_column']} "
        f"sha256={fingerprint['sha256'][:12]}  device={device}",
        flush=True,
    )

    # load base model (frozen) via the official loader
    predictor = P2PAMPredictor(args.model, device="cpu", local_files_only=True)
    model = predictor.model
    # freeze the WHOLE base first (LoRA residual + classifier head stay trainable)
    for p in model.parameters():
        p.requires_grad_(False)
    # swap classifier output shape already 10x4 (config.output_shape)
    n_lora = apply_lora(model, rank=args.lora_rank, alpha=args.lora_alpha)
    # The classifier is a 3.33M-parameter component. Keep it frozen for the F0
    # control, or train it at a separately configured LR for F1.
    for p in model.classifier.parameters():
        p.requires_grad_(args.classifier_mode == "train")
    trainable = [p for p in model.parameters() if p.requires_grad]
    classifier_parameters = sum(
        p.numel() for p in model.classifier.parameters() if p.requires_grad
    )
    print(
        f"LoRA adapters: {n_lora}  classifier={args.classifier_mode} "
        f"classifier_trainable={classifier_parameters / 1e6:.3f}M "
        f"total_trainable={sum(p.numel() for p in trainable) / 1e6:.3f}M "
        f"loss={args.loss_mode}",
        flush=True,
    )

    out_dir = WORKSPACE / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    state_path = out_dir / "trainer_state.json"
    checkpoint_path = out_dir / "last_checkpoint.pt"
    state = {
        "epoch": 0,
        "step": 0,
        "optimizer_steps": 0,
        "loss_sum": 0.0,
        "rows_seen": 0,
    }
    resume_payload = None
    if checkpoint_path.exists() and not args.no_resume:
        resume_payload = torch.load(checkpoint_path, map_location="cpu")
        if not isinstance(resume_payload, dict) or resume_payload.get("format_version") != 3:
            raise ValueError(
                f"unsupported resume checkpoint: {checkpoint_path}; "
                "use --no-resume to start cleanly"
            )
        expected = {
            "model_name": args.model,
            "lora_rank": args.lora_rank,
            "lora_alpha": args.lora_alpha,
            "batch_size": args.batch_size,
            "grad_accum": args.grad_accum,
            "train_path": args.train,
            "classifier_mode": args.classifier_mode,
            "classifier_lr": args.classifier_lr,
            "loss_mode": args.loss_mode,
            "label_smoothing": args.label_smoothing,
            "train_fingerprint": fingerprint,
        }
        mismatches = {
            key: (resume_payload.get(key), value)
            for key, value in expected.items()
            if resume_payload.get(key) != value
        }
        if mismatches:
            raise ValueError(
                f"resume arguments differ from {checkpoint_path}: {mismatches}; "
                "use the original arguments or pass --no-resume"
            )
        saved = resume_payload.get("trainable_state")
        if not isinstance(saved, dict) or not saved:
            raise ValueError(f"checkpoint has no trainable_state: {checkpoint_path}")
        current = model.state_dict()
        expected_trainable = {
            name for name, parameter in model.named_parameters()
            if parameter.requires_grad
        }
        missing_trainable = sorted(expected_trainable.difference(saved))
        unknown = sorted(set(saved).difference(current))
        wrong_shape = sorted(
            key for key, value in saved.items()
            if key in current and tuple(value.shape) != tuple(current[key].shape)
        )
        if missing_trainable or unknown or wrong_shape:
            raise ValueError(
                f"checkpoint/model mismatch; missing_trainable="
                f"{missing_trainable[:5]}, unknown={unknown[:5]}, "
                f"wrong_shape={wrong_shape[:5]}"
            )
        model.load_state_dict(saved, strict=False)
        state.update(resume_payload.get("trainer_state") or {})
        print(
            f"resuming from epoch {state['epoch']} batch {state['step']} "
            f"optimizer_step {state['optimizer_steps']}",
            flush=True,
        )
    elif state_path.exists() and not checkpoint_path.exists() and not args.no_resume:
        print(
            f"warning: ignoring legacy {state_path}; it has no model/optimizer "
            "weights and is not a valid resume checkpoint",
            flush=True,
        )
    elif checkpoint_path.exists() and args.no_resume:
        raise FileExistsError(
            f"--no-resume refuses to overwrite existing {checkpoint_path}; "
            "choose a new --out directory"
        )

    # Convert on CPU first. Moving FP32 to CUDA and converting there creates a
    # large temporary allocation and left the 8GB GPU almost completely full.
    use_bf16 = device.type == "cuda" and torch.cuda.is_bf16_supported()
    if use_bf16:
        model.to(dtype=torch.bfloat16)
        print("weights: bf16", flush=True)
    model.to(device)
    model.train()

    attention = model.esm.encoder.layer[0].attention.self
    print(
        f"flash attention: {bool(getattr(attention, 'use_flash_attn', False))}  "
        f"gradient checkpointing: "
        f"{bool(getattr(model.esm.encoder, 'gradient_checkpointing', False))}",
        flush=True,
    )
    if device.type == "cuda":
        print(
            f"cuda memory after load: allocated="
            f"{torch.cuda.memory_allocated(device) / 2**30:.2f} GiB  reserved="
            f"{torch.cuda.memory_reserved(device) / 2**30:.2f} GiB",
            flush=True,
        )
        torch.cuda.reset_peak_memory_stats(device)

    # Separate groups prevent the large classifier from moving at the LoRA LR.
    lora_parameters = [
        parameter for name, parameter in model.named_parameters()
        if parameter.requires_grad and "lora_" in name
    ]
    parameter_groups: list[dict] = [
        {"params": lora_parameters, "lr": args.lr, "group_name": "lora"}
    ]
    if args.classifier_mode == "train":
        parameter_groups.append({
            "params": [
                parameter for parameter in model.classifier.parameters()
                if parameter.requires_grad
            ],
            "lr": args.classifier_lr,
            "group_name": "classifier",
        })
    opt = torch.optim.AdamW(parameter_groups, weight_decay=0.01)
    if resume_payload is not None:
        opt.load_state_dict(resume_payload["optimizer_state"])
        _move_optimizer_state(opt, device)
        if "torch_rng_state" in resume_payload:
            torch.set_rng_state(resume_payload["torch_rng_state"])
        if device.type == "cuda" and resume_payload.get("cuda_rng_state_all"):
            torch.cuda.set_rng_state_all(resume_payload["cuda_rng_state_all"])
    total_steps = int(state["optimizer_steps"])
    if args.max_steps and total_steps >= args.max_steps:
        print(
            f"checkpoint already reached --max-steps={args.max_steps}; "
            "nothing to do",
            flush=True,
        )
        return 0

    tokenizer = predictor.tokenizer
    t0 = time.time()

    def checkpoint_payload(include_optimizer: bool) -> dict:
        payload = {
            "format_version": 3,
            "model_name": args.model,
            "lora_rank": args.lora_rank,
            "lora_alpha": args.lora_alpha,
            "batch_size": args.batch_size,
            "grad_accum": args.grad_accum,
            "train_path": args.train,
            "train_fingerprint": fingerprint,
            "classifier_mode": args.classifier_mode,
            "classifier_lr": args.classifier_lr,
            "loss_mode": args.loss_mode,
            "label_smoothing": args.label_smoothing,
            "trainable_state": _trainable_state_dict(model),
            "trainer_state": dict(state),
        }
        if include_optimizer:
            payload.update({
                "optimizer_state": opt.state_dict(),
                "torch_rng_state": torch.get_rng_state(),
                "cuda_rng_state_all": (
                    torch.cuda.get_rng_state_all() if device.type == "cuda" else []
                ),
            })
        return payload

    def save_resume_checkpoint() -> None:
        _atomic_torch_save(checkpoint_payload(include_optimizer=True), checkpoint_path)
        _write_state_json(state, state_path)
        if args.keep_step_checkpoints:
            immutable_path = out_dir / (
                f"checkpoint_step{int(state['optimizer_steps']):06d}.pt"
            )
            # The save interval and a controlled max-steps exit can request
            # the same optimizer step twice. Retaining the first artifact is
            # idempotent and still preserves the no-overwrite guarantee.
            if not immutable_path.exists():
                _atomic_torch_save(
                    checkpoint_payload(include_optimizer=True), immutable_path
                )

    trace_first_batch = True
    opt.zero_grad(set_to_none=True)
    try:
        for epoch in range(int(state["epoch"]), args.epochs):
            if epoch != int(state["epoch"]):
                state.update({"epoch": epoch, "step": 0,
                              "loss_sum": 0.0, "rows_seen": 0})
            batches = _epoch_batches(rows, args.batch_size, args.seed + epoch)
            if int(state["step"]) > len(batches):
                raise ValueError(
                    f"resume batch {state['step']} exceeds epoch size {len(batches)}"
                )
            pending_loss_sum = 0.0
            pending_rows = 0
            for step_in_epoch, batch_idx in enumerate(batches):
                if step_in_epoch < int(state["step"]):
                    continue
                batch = [rows[i] for i in batch_idx]
                seqs = [r["_sequence"] for r in batch]
                trace = trace_first_batch
                if trace:
                    lengths = [len(seq) for seq in seqs]
                    print(
                        f"[first batch] epoch={epoch} batch={step_in_epoch} "
                        f"items={len(batch)} aa_len={min(lengths)}..{max(lengths)}; "
                        "tokenizing...",
                        flush=True,
                    )

                phase = time.perf_counter()
                encs = tokenizer.encode_batch(seqs)
                if trace:
                    print(
                        f"[first batch] tokenize={time.perf_counter() - phase:.3f}s "
                        f"tokens={len(encs[0].ids)}; moving tensors...",
                        flush=True,
                    )
                phase = time.perf_counter()
                input_ids = torch.tensor([e.ids for e in encs], device=device)
                attn = torch.tensor([e.attention_mask for e in encs], device=device)
                targets = torch.stack([logo_of(r) for r in batch]).to(device)
                if args.label_smoothing:
                    targets = (
                        targets * (1.0 - args.label_smoothing)
                        + args.label_smoothing / 4.0
                    )
                position_weights = torch.stack(
                    [position_weights_of(r) for r in batch]
                ).to(device)
                if trace and device.type == "cuda":
                    torch.cuda.synchronize(device)
                if trace:
                    print(
                        f"[first batch] tensors={time.perf_counter() - phase:.3f}s "
                        f"shape={tuple(input_ids.shape)}; forward...",
                        flush=True,
                    )

                phase = time.perf_counter()
                with torch.autocast(device_type=device.type, dtype=torch.bfloat16,
                                    enabled=use_bf16):
                    logits = model(input_ids=input_ids, attention_mask=attn).logits
                if trace and device.type == "cuda":
                    torch.cuda.synchronize(device)
                if trace:
                    print(
                        f"[first batch] forward={time.perf_counter() - phase:.3f}s "
                        f"logits={tuple(logits.shape)}; backward...",
                        flush=True,
                    )

                logq = F.log_softmax(logits.float(), dim=-1)
                per_position_kl = F.kl_div(
                    logq, targets, reduction="none"
                ).sum(-1)
                if args.loss_mode == "info-weighted":
                    per_row = (
                        (per_position_kl * position_weights).sum(-1)
                        / position_weights.sum(-1)
                    )
                    raw_loss = per_row.mean()
                else:
                    raw_loss = per_position_kl.mean()
                group_start = (step_in_epoch // args.grad_accum) * args.grad_accum
                group_end = min(group_start + args.grad_accum, len(batches))
                accumulation_size = group_end - group_start
                loss = raw_loss / accumulation_size

                phase = time.perf_counter()
                loss.backward()
                if trace and device.type == "cuda":
                    torch.cuda.synchronize(device)
                if trace:
                    peak = (
                        f" peak_cuda={torch.cuda.max_memory_allocated(device) / 2**30:.2f}GiB"
                        if device.type == "cuda" else ""
                    )
                    print(
                        f"[first batch] backward={time.perf_counter() - phase:.3f}s "
                        f"loss={raw_loss.item():.5f}{peak}",
                        flush=True,
                    )
                    trace_first_batch = False

                pending_loss_sum += raw_loss.item() * len(batch)
                pending_rows += len(batch)
                if step_in_epoch + 1 == group_end:
                    grad_norm = torch.nn.utils.clip_grad_norm_(trainable, 1.0)
                    opt.step()
                    opt.zero_grad(set_to_none=True)
                    total_steps += 1
                    state["step"] = step_in_epoch + 1
                    state["optimizer_steps"] = total_steps
                    state["loss_sum"] += pending_loss_sum
                    state["rows_seen"] += pending_rows
                    pending_loss_sum = 0.0
                    pending_rows = 0

                    if total_steps == 1 or total_steps % args.log_every == 0:
                        mean_loss = state["loss_sum"] / max(1, state["rows_seen"])
                        print(
                            f"epoch {epoch} optimizer_step {total_steps} "
                            f"next_batch {state['step']}/{len(batches)} "
                            f"loss {mean_loss:.5f} grad_norm {float(grad_norm):.3f} "
                            f"({state['rows_seen']} rows, {time.time() - t0:.0f}s)",
                            flush=True,
                        )
                    if total_steps % args.save_every == 0:
                        save_resume_checkpoint()
                        print(f"resume checkpoint -> {checkpoint_path}", flush=True)
                    if args.max_steps and total_steps >= args.max_steps:
                        break

            epoch_complete = int(state["step"]) == len(batches)
            if not epoch_complete:
                save_resume_checkpoint()
                print(
                    f"stopped at optimizer_step {total_steps} before epoch {epoch} "
                    f"completed; resume checkpoint -> {checkpoint_path}",
                    flush=True,
                )
                break

            completed_epoch = epoch
            state.update({"epoch": epoch + 1, "step": 0,
                          "loss_sum": 0.0, "rows_seen": 0})
            epoch_path = out_dir / f"lora_epoch{completed_epoch}.pt"
            _atomic_torch_save(checkpoint_payload(include_optimizer=False), epoch_path)
            save_resume_checkpoint()
            print(f"epoch {completed_epoch} done -> {epoch_path}", flush=True)
            if args.max_steps and total_steps >= args.max_steps:
                break
    except KeyboardInterrupt:
        opt.zero_grad(set_to_none=True)
        save_resume_checkpoint()
        print(f"interrupted; safe resume checkpoint -> {checkpoint_path}", flush=True)
        return 130

    print("training complete", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
