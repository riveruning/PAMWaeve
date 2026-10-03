#!/usr/bin/env python3
"""Capture actual execution logs and provenance without changing scoring."""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default=os.environ.get("DEVICE", "cpu"))
    parser.add_argument("--fresh", action="store_true", help="rerun inference without cached outputs")
    parser.add_argument("--outdir", help="output directory relative to package root")
    parser.add_argument("--logdir", help="new log directory relative to package root")
    args = parser.parse_args()
    started = datetime.now(timezone.utc)
    stamp = started.strftime("%Y%m%dT%H%M%S_%fZ")
    logdir = ROOT / (args.logdir or f"logs/runs/{stamp}")
    logdir.mkdir(parents=True, exist_ok=False)
    full = all((ROOT / path).is_file() for path in ["data/raw/cas9_full.fasta", "data/raw/external_datasets/PAMpredict/Example/spacers.fna", "data/raw/external_datasets/PAMpredict/Example/Phages/phages.fna"])
    manifest = f"benchmarks/multi_evidence_v2/manifest_{'final22' if full else 'redistributable21'}.json"
    outdir = args.outdir or f"data/parsed/bench_{'final22' if full else 'repro21'}"
    results_path = f"results/results{'' if full else '_repro21'}.csv"
    systems = ["sp7f7-published-flanks", "cj4-campylobacter-jejuni-414-protein-only", "asp-acidiphilium-21-60-14-protein-only", "cba-caulobacterales-protein-only"]
    if full:
        systems.insert(0, "spcas9-pampredict-example")
    env = dict(os.environ)
    env.update(PYTHONPATH=str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", ""), HF_HOME=env.get("HF_HOME", str(ROOT / "data/checkpoints/hf")), PYTHONHASHSEED="17", PYTHONUNBUFFERED="1")
    # The benchmark uses this actual default; do not invent a training seed.
    from pamdict.benchmark.multi_evidence import aggregate_benchmark
    bootstrap_seed = aggregate_benchmark.__kwdefaults__["seed"]
    replacements = [(str(ROOT), "<PACKAGE_ROOT>"), (sys.executable, "<PYTHON>")]
    for variable, label in [("HF_HOME", "<MODEL_CACHE>"), ("PAMPRIDICT_PROTEIN2PAM_HOME", "<UPSTREAM_ADAPTER>")]:
        if env.get(variable):
            replacements.append((env[variable], label))
    for entry in env["PYTHONPATH"].split(os.pathsep):
        if entry and Path(entry).is_absolute() and not Path(entry).is_relative_to(ROOT):
            replacements.append((entry, "<EXTERNAL_PYTHONPATH>"))
    replacements.sort(key=lambda item: len(item[0]), reverse=True)
    def portable(value: str) -> str:
        for original, replacement in replacements:
            value = value.replace(original, replacement)
        return value
    packages = {}
    for package in ("numpy", "torch", "transformers", "tokenizers", "safetensors", "huggingface_hub"):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = "not installed"
    source_files = sorted([*ROOT.glob("src/**/*.py"), *ROOT.glob("scripts/*.py"), ROOT / "run_all.sh", ROOT / "requirements.txt"])
    cached = [row["system_id"] for row in json.loads((ROOT / manifest).read_text())["systems"] if (ROOT / outdir / row["system_id"] / "protein_scores.json").is_file()]
    record = {
        "format_version": 1, "started_utc": started.isoformat(), "status": "running",
        "verification_scope": "local execution; not independent-machine acceptance",
        "manifest": manifest, "manifest_sha256": sha256(ROOT / manifest),
        "device": args.device, "fresh_requested": args.fresh,
        "systems_with_existing_protein_outputs_at_start": cached,
        "bootstrap_iterations": 2000, "bootstrap_seed": bootstrap_seed, "python_hash_seed": 17,
        "inference_randomness": "pretrained model in eval mode; no training or new PyTorch inference seed is claimed",
        "model": "Profluent-Bio/protein2pam-cas9_full",
        "documented_model_revision": "407f7fc32146a4c0db13c05284f9f3a7cf0ff612",
        "environment": {"python": platform.python_version(), "os": platform.system(), "os_release": platform.release(), "architecture": platform.machine(), "cpu_logical_count": os.cpu_count(), "packages": packages},
        "source_sha256": {path.relative_to(ROOT).as_posix(): sha256(path) for path in source_files},
        "stages": [], "outputs": {},
    }
    def save_record() -> None:
        (logdir / "run.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    def run_stage(name: str, command: list[str]) -> int:
        stage = {"name": name, "command": [portable(item) for item in command], "started_utc": datetime.now(timezone.utc).isoformat(), "log": f"{name}.log"}
        tick = time.monotonic()
        print(f"==> {name}: {portable(shlex.join(command))}", flush=True)
        with (logdir / stage["log"]).open("w", encoding="utf-8") as log:
            log.write("Started UTC: " + stage["started_utc"] + "\nCommand: " + portable(shlex.join(command)) + "\n")
            log.flush()
            process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", bufsize=1)
            assert process.stdout is not None
            for line in process.stdout:
                line = portable(line)
                log.write(line); log.flush()
                print(line, end="", flush=True)
            code = process.wait()
            stage.update(finished_utc=datetime.now(timezone.utc).isoformat(), seconds=round(time.monotonic() - tick, 3), exit_code=code)
            log.write(f"\nExit code: {code}\nElapsed seconds: {stage['seconds']}\n")
        record["stages"].append(stage); save_record()
        return code
    benchmark = [sys.executable, "-u", "scripts/run_multi_evidence_benchmark.py", "--manifest", manifest, "--device", args.device, "--fusion-policy", "abstain", "--bootstrap", "2000", "--outdir", outdir]
    if not args.fresh:
        benchmark.append("--reuse-existing")
    stages = [
        ("01_validate", [sys.executable, "scripts/validate_benchmark_manifest.py", "--manifest", manifest, "--check-files"]),
        ("02_benchmark", benchmark),
        ("03_candidates", [sys.executable, "scripts/build_competition_results.py", "--benchmark-dir", outdir, "--manifest", manifest, "--display-candidates", "results/display_candidates.json", "--top-n", "5", "--systems", ",".join(systems), "--out", results_path]),
    ]
    save_record(); overall = time.monotonic(); code = 0
    try:
        for name, command in stages:
            code = run_stage(name, command)
            if code:
                break
        if code == 0:
            for relative in [results_path, f"{outdir}/benchmark_report.json", f"{outdir}/benchmark_report.md", f"{outdir}/benchmark_systems.tsv", f"{outdir}/manifest_audit.json"]:
                path = ROOT / relative
                record["outputs"][relative] = {"sha256": sha256(path), "bytes": path.stat().st_size}
            with (ROOT / results_path).open(encoding="utf-8", newline="") as handle:
                record["candidate_rows"] = sum(1 for _ in csv.DictReader(handle))
    except BaseException as error:
        code = 130 if isinstance(error, KeyboardInterrupt) else 1
        record["error"] = portable(f"{type(error).__name__}: {error}")
        raise
    finally:
        record.update(status="completed" if code == 0 else "failed", exit_code=code, finished_utc=datetime.now(timezone.utc).isoformat(), seconds=round(time.monotonic() - overall, 3))
        save_record()
    print(f"Candidate list: {results_path}\nExecution record: {logdir.relative_to(ROOT)}\nExit code: {code}", flush=True)
    return code

if __name__ == "__main__":
    raise SystemExit(main())
