"""Local bounded genome-to-PAM workflow, without an LLM or arbitrary commands.

State machine::

    validated -> annotating -> (needs_system_selection | no_supported_cas9)
             -> (needs_array_confirmation) -> predicting -> complete
    any stage may end in blocked_dependency or failed

Every stage writes ``job.json`` atomically and a per-stage log.  A failed or
blocked run never produces a "success" report.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import html
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pamdict.genome import (  # noqa: E402
    check_genome,
    collect_systems,
    fasta,
    select_system,
    sha256,
)

DEPENDENCY_HINT = (
    "CRISPRCasTyper is not on PATH. Install it into an isolated environment and "
    "pass the executable explicitly, e.g.\n"
    "  conda create -y --prefix <workspace>/.envs/cctyper -c conda-forge -c bioconda "
    "-c russel88 cctyper\n"
    "  python scripts/run_genome_agent.py --genome G.fna --outdir OUT "
    "--annotator <workspace>/.envs/cctyper/bin/cctyper\n"
    "No prediction was made."
)

STATUS_LABELS = {
    "validated": "输入校验通过",
    "annotating": "正在运行 CRISPR-Cas 注释",
    "needs_system_selection": "存在多个合格系统，等待选择",
    "needs_array_confirmation": "等待确认 array 与 spacer 方向",
    "no_supported_cas9": "本次未识别到可用完整 Cas9",
    "predicting": "正在运行 Cas9 模型推理",
    "complete": "完成",
    "blocked_dependency": "依赖缺失，已阻塞",
    "failed": "失败",
}


def now():
    return datetime.now(timezone.utc).isoformat()


def save(out, job):
    """Atomically persist job.json and rewrite the human-readable index."""
    job["updated_utc"] = now()
    tmp = out / "job.json.tmp"
    tmp.write_text(json.dumps(job, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp.replace(out / "job.json")

    status = job.get("status", "unknown")
    label = STATUS_LABELS.get(status, status)
    banner = (
        '<p class="warn">仅支持完整 Cas9；注释与预测均不等于实验验证。</p>'
        if status != "complete"
        else '<p class="ok">运行完成；结果仍不是实验验证。</p>'
    )
    links = []
    if (out / "pam" / "report.html").is_file():
        links.append('<a href="pam/report.html">打开 PAM 报告</a>')
    if (out / "pam" / "run.json").is_file():
        links.append('<a href="pam/run.json">查看推理运行记录</a>')
    for name in ("annotating.log", "predicting.log"):
        if (out / name).is_file():
            links.append(f'<a href="{name}">查看 {name}</a>')
    message = job.get("message") or job.get("error") or ""
    page = (
        '<!doctype html><html lang="zh-CN"><meta charset="utf-8">'
        "<title>基因组到 PAM 分析</title><style>"
        "body{font-family:system-ui,sans-serif;margin:32px;color:#182938;max-width:1100px}"
        "pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f4f6f8;padding:14px}"
        ".warn{background:#fff2cf;padding:12px}.ok{background:#e4f6e8;padding:12px}"
        "a{margin-right:14px}</style>"
        "<h1>基因组 → CRISPR-Cas → PAM</h1>"
        f'<h2>状态：{html.escape(label)}</h2>{banner}'
        + (f'<p class="warn">{html.escape(str(message))}</p>' if message else "")
        + "<p>" + "".join(links) + "</p>"
        "<h2>作业记录</h2><pre>"
        + html.escape(json.dumps(job, ensure_ascii=False, indent=2))
        + "</pre></html>"
    )
    (out / "index.html").write_text(page, encoding="utf-8")


def run_stage(command, out, job, stage, timeout, env=None):
    """Run one child stage in its own process group, cleaning it up on timeout.

    ``subprocess.run(timeout=...)`` only kills the direct child; a killed
    annotator can leave grandchildren (HMMER/BLAST workers) running.  Starting
    the child in a new session and killing that whole group bounds the damage
    to this job's processes only.
    """
    job["status"] = stage
    job.setdefault("commands", []).append(command)
    save(out, job)
    print(f"[{stage}] {' '.join(command)}", flush=True)
    log_path = out / f"{stage}.log"
    started = time.monotonic()
    try:
        with log_path.open("w", encoding="utf-8") as log:
            process = subprocess.Popen(
                command,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                cwd=ROOT,
                start_new_session=True,
                text=True,
            )
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                _terminate_group(process)
                raise TimeoutError(
                    f"Stage {stage} exceeded {timeout}s and its process group was terminated"
                )
    finally:
        job.setdefault("stage_seconds", {})[stage] = round(time.monotonic() - started, 3)
    if code != 0:
        raise RuntimeError(f"Stage {stage} exited with code {code}; see {log_path.name}")


def _terminate_group(process):
    """Terminate the child's whole process group, escalating if needed.

    Only the group created for this stage is signalled, so unrelated jobs on
    the machine are untouched.
    """
    if process.poll() is not None:
        return
    try:
        pgid = os.getpgid(process.pid)
    except ProcessLookupError:
        return
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(pgid, sig)
        except ProcessLookupError:
            return
        try:
            process.wait(timeout=10)
            return
        except subprocess.TimeoutExpired:
            continue


def annotation_hashes(directory):
    return {
        str(p.relative_to(directory)): sha256(p)
        for p in sorted(Path(directory).rglob("*"))
        if p.is_file()
    }


def annotation_upstream_error(directory, log_path):
    """Recognize a known crash signature, never a biological negative.

    Missing outputs and a KeyError do not establish that the search completed.
    This helper only adds an actionable diagnostic to a failed annotation.
    """
    directory = Path(directory)
    log_path = Path(log_path)
    if (directory / "hmmer.tab").exists():
        return None
    if any((directory / n).exists() for n in (
        "cas_operons.tab", "cas_operons_putative.tab",
        "cas_operons_orphan.tab", "crisprs_all.tab", "CRISPR_Cas.tab",
    )):
        return None
    # CCTyper creates an empty hmmer/ scratch directory before crashing, so its
    # mere existence means nothing; only actual hits there would.
    hmmer_dir = directory / "hmmer"
    if hmmer_dir.is_dir() and any(hmmer_dir.iterdir()):
        return None
    if not log_path.is_file():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    if "['Hmm', 'ORF'] not in index" not in text:
        return None
    if "KeyError" not in text:
        return None
    return (
        "上游 CRISPRCasTyper 注释错误（HMMER 表处理 KeyError）；搜索未完成，"
        "无法判断是否存在 Cas9 或 CRISPR-Cas 系统。不是阴性结果，未生成 PAM。"
        "请检查 annotating.log，修复注释工具后用新输出目录重跑。"
    )


def annotation_is_scientifically_empty(directory):
    """Detect "exit 0 but no usable output", which is not a reliable negative.

    Only outputs CCTyper *always* writes are required.  CRISPR array tables are
    deliberately optional: when a genome has a Cas operon but no array, CCTyper
    logs "No CRISPRs found.", deletes ``crisprs_all.tab`` and writes no
    ``CRISPR_Cas.tab``.  That is a legitimate orphan-Cas9 result which must still
    reach protein-only inference, so its absence is not treated as failure.
    Whether arrays exist is decided later, from the parsed systems.
    """
    directory = Path(directory)
    # genes.tab comes from gene prediction and proteins.faa from translation;
    # both are produced before any CRISPR/Cas search, and a run that lacks them
    # produced nothing usable.
    required = ["genes.tab", "proteins.faa"]
    missing = [n for n in required if not (directory / n).is_file()]
    if missing:
        return f"CCTyper exited successfully but did not write: {', '.join(missing)}"
    # At least one operon table must exist, otherwise no Cas call was made at all.
    if not (directory / "cas_operons.tab").exists() and not (
        directory / "cas_operons_putative.tab"
    ).exists() and not (directory / "cas_operons_orphan.tab").exists():
        return (
            "CCTyper exited successfully but wrote no cas_operons.tab, "
            "cas_operons_putative.tab or cas_operons_orphan.tab, so no Cas operon "
            "call was produced"
        )
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", type=Path, required=True, help="assembled nucleotide FASTA, <=30 MB")
    ap.add_argument("--outdir", type=Path, required=True, help="new directory; refuses to overwrite")
    ap.add_argument("--annotator", default="cctyper", help="installed cctyper executable path (not a shell command)")
    ap.add_argument("--annotation-timeout", type=int, default=1800)
    ap.add_argument("--inference-timeout", type=int, default=900)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--device", choices=["cpu", "cuda", "auto"], default="auto")
    ap.add_argument("--resume", action="store_true", help="reuse an unchanged completed annotation after selection")
    ap.add_argument("--system-id", help="Cas9 system ID from job.json")
    ap.add_argument("--array-id", help="linked trusted array ID from job.json")
    ap.add_argument("--spacer-orientation", choices=["input", "reverse"], help="explicit user confirmation; never inferred from the Cas strand")
    ap.add_argument("--targets", type=Path, action="append", help="optional local phage/plasmid FASTA; never uploaded or downloaded")
    ap.add_argument("--candidates", help="optional same-length IUPAC list; exclusive with --auto-candidates")
    ap.add_argument("--auto-candidates", action="store_true", help="report the model consensus only (default when no --candidates)")
    args = ap.parse_args()

    # --- validation ------------------------------------------------------
    try:
        if min(args.annotation_timeout, args.inference_timeout, args.threads) <= 0:
            raise ValueError("Timeouts and threads must be positive")
        if args.candidates and args.auto_candidates:
            raise ValueError(
                "--candidates and --auto-candidates are mutually exclusive; "
                "pass one, not both, so neither is silently ignored"
            )
        genome = args.genome.resolve()
        summary = check_genome(genome)
        target_paths = [p.resolve() for p in (args.targets or [])]
        for target in target_paths:
            check_genome(target)
        if args.spacer_orientation and not target_paths:
            raise ValueError("--spacer-orientation requires --targets")
        out = args.outdir.resolve()
    except (ValueError, OSError) as exc:
        ap.error(str(exc))

    job = None
    try:
        if args.resume:
            job_path = out / "job.json"
            if not job_path.is_file():
                raise ValueError(f"--resume requires an existing {job_path}")
            job = json.loads(job_path.read_text(encoding="utf-8"))
            if job.get("genome", {}).get("sha256") != summary["sha256"]:
                raise ValueError("Resume rejected: the genome file changed since annotation")
            # A pam/ directory from a run that never reached 'complete' is a
            # partial prediction.  It is preserved, never overwritten or reused.
            if (out / "pam").exists() and job.get("status") != "complete":
                raise ValueError(
                    "Resume rejected: a partial prediction directory already exists. "
                    "It was preserved, not overwritten; inspect it and start a new job directory."
                )
            # Integrity of the cached annotation is checked before the job's own
            # status, so tampering is reported as tampering rather than being
            # masked by an "already complete" message.
            expected = job.get("annotation_sha256")
            annotation_dir = out / "annotation"
            if not expected or not annotation_dir.is_dir():
                raise ValueError("Resume requires a completed annotation to reuse")
            for rel, digest in expected.items():
                path = annotation_dir / rel
                if not path.is_file():
                    raise ValueError(f"Resume rejected: annotation output is missing: {rel}")
                if sha256(path) != digest:
                    raise ValueError(f"Resume rejected: annotation output changed: {rel}")
            actual = annotation_hashes(annotation_dir)
            if set(actual) != set(expected):
                extra = sorted(set(actual) - set(expected))
                gone = sorted(set(expected) - set(actual))
                raise ValueError(
                    "Resume rejected: annotation file set changed "
                    f"(added: {extra}, removed: {gone})"
                )
            if job.get("status") == "complete":
                raise ValueError(
                    "Job is already complete and its annotation is unchanged; "
                    "use a new output directory to change inputs"
                )
            job.setdefault("resume_history", []).append({"resumed_utc": now()})
        else:
            out.mkdir(parents=True, exist_ok=False)
            job = {
                "format_version": 2,
                "status": "validated",
                "genome": summary,
                "genome_path": str(genome),
                "commands": [],
                "limitations": [
                    "Cas9 (Type II) only; Cas12/Cas13/Type I and other systems are listed but not scored",
                    "The predicted protein and PAM are not experimentally validated",
                    "No automatic external target search; no fusion superiority claim",
                    "Array linkage is a CCTyper annotation hypothesis and is not experimental evidence",
                    "Spacer orientation is never inferred from the Cas strand; the user must confirm it",
                    "Absence of a reported Cas9 does not prove the organism lacks CRISPR systems",
                ],
                "workflow_sha256": sha256(Path(__file__)),
                "adapter_sha256": sha256(ROOT / "src/pamdict/genome.py"),
            }
        job["request"] = {k: (str(v) if not isinstance(v, list) else [str(x) for x in v]) for k, v in vars(args).items()}
        save(out, job)
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        # A rejected --resume must not touch the existing job: overwriting it with
        # status='failed' would destroy the record of a completed successful run.
        # For a brand-new job the failure is recorded so the user can inspect it.
        if job is not None and not args.resume and out.is_dir():
            job.update(status="failed", error=str(exc),
                       message="Request rejected before any computation; no outputs were produced.")
            save(out, job)
        ap.error(str(exc))

    # --- annotation ------------------------------------------------------
    try:
        if not args.resume:
            executable = shutil.which(args.annotator) if os.sep not in args.annotator else args.annotator
            if not executable or not Path(executable).exists():
                job.update(status="blocked_dependency", message=DEPENDENCY_HINT)
                save(out, job)
                print(DEPENDENCY_HINT, file=sys.stderr)
                return 3
            executable = str(Path(executable).resolve())
            env = os.environ.copy()
            env["PATH"] = str(Path(executable).parent) + os.pathsep + env.get("PATH", "")
            # Keep annotation dependencies out of the model environment.
            env.pop("PYTHONPATH", None)
            env.pop("PYTHONHOME", None)
            annotation_dir = out / "annotation"
            command = [
                executable, str(genome), str(annotation_dir),
                "--keep_tmp", "-t", str(args.threads),
            ]
            try:
                run_stage(command, out, job, "annotating", args.annotation_timeout, env)
            except RuntimeError as exc:
                diagnostic = annotation_upstream_error(
                    annotation_dir, out / "annotating.log"
                )
                if diagnostic is None:
                    raise
                job.update(status="failed", error=str(exc), message=diagnostic,
                           annotation_outcome="unknown",
                           annotated_system_count=None,
                           upstream_annotator_bug="KeyError: ['Hmm', 'ORF'] not in index")
                save(out, job)
                print(diagnostic, file=sys.stderr)
                return 1
            if sha256(genome) != summary["sha256"]:
                raise ValueError("The genome file changed while annotation was running")
            empty_reason = annotation_is_scientifically_empty(annotation_dir)
            if empty_reason:
                raise RuntimeError(
                    empty_reason
                    + ". An exit code of 0 without usable output is not a reliable negative result."
                )
            job["annotation_sha256"] = annotation_hashes(annotation_dir)
            job["annotator"] = {
                "executable": executable,
                "version_note": "CRISPRCasTyper 1.9.0 (bioconda); parses pyrodigal-gv headers without partial=",
            }
            save(out, job)

        # --- system selection --------------------------------------------
        systems = collect_systems(out / "annotation", summary, genome_path=genome)
        job["systems"] = [
            {k: v for k, v in s.items() if k != "protein_sequence"} for s in systems
        ]
        job["annotated_system_count"] = len(systems)
        selected = select_system(systems, args.system_id, allow_flagged=True)

        if selected is None:
            eligible = [s for s in systems if not s["issues"]]
            if eligible:
                status, message, code = (
                    "needs_system_selection",
                    "Multiple eligible Cas9 systems. Re-run with --resume --system-id <ID>. "
                    "Candidates: " + ", ".join(s["system_id"] for s in eligible),
                    2,
                )
            elif systems:
                status, message, code = (
                    "no_supported_cas9",
                    "A Cas9-like gene was found but no system passed automatic quality checks. "
                    "CCTyper subtype/annotation flags: " + "; ".join(
                        f"{s['system_id']}[{','.join(s['issues'])}]" for s in systems
                    ) + ". You may review a specific one explicitly with "
                    "--resume --system-id <ID>; that overrides only annotation-confidence "
                    "flags, never a protein that failed genome-level integrity verification.",
                    2,
                )
            else:
                status, message, code = (
                    "no_supported_cas9",
                    "CRISPRCasTyper reported no Cas9/Csn1 gene for this genome. This means only "
                    "that this tool run found no usable complete Cas9; it does not prove the "
                    "organism has no CRISPR system.",
                    0,
                )
            job.update(status=status, message=message)
            save(out, job)
            print(message, flush=True)
            return code

        job["selected_system"] = selected["system_id"]
        job["selected_system_evidence"] = {
            "subtype": selected["subtype"],
            "operon_table": selected["operon_table"],
            "integrity_verified": selected["integrity_verified"],
            "integrity_notes": selected["integrity_notes"],
            "protein_length": selected["protein_length"],
            "arrays": selected["arrays"],
        }
        if selected["issues"]:
            # Reached only through an explicit --system-id, since automatic
            # selection never returns a flagged system.
            job["selection_override"] = {
                "requested_system_id": args.system_id,
                "overridden_issues": selected["issues"],
                "note": (
                    "The user explicitly selected this system despite the listed annotation "
                    "flags. CCTyper's subtype call is uncertain for it. The protein itself "
                    "was reproduced from the genome, but the system assignment is not "
                    "confirmed and the PAM result inherits that uncertainty."
                ),
            }
            job.setdefault("warnings", []).append(job["selection_override"]["note"])
        if not selected["arrays"]:
            job["warnings"] = job.get("warnings", []) + [
                "Selected Cas9 has no linked CRISPR array: this is an orphan Cas9 and the "
                "result is protein-only, with no array evidence."
            ]

        # --- optional spacer evidence ------------------------------------
        spacer_args = []
        if target_paths:
            trusted = [
                a for a in selected["arrays"]
                if a["trusted"] and not a["issues"] and a["spacer_file_present"]
            ]
            if args.array_id:
                trusted = [a for a in trusted if a["array_id"] == args.array_id]
            if len(trusted) != 1 or not args.spacer_orientation:
                job.update(
                    status="needs_array_confirmation",
                    message=(
                        "Spacer evidence needs exactly one linked trusted array and an explicit "
                        "--spacer-orientation input|reverse. The array direction is NOT inferred "
                        "from the Cas strand. Available trusted arrays: "
                        + (", ".join(a["array_id"] for a in trusted) or "none")
                        + ". Omit --targets for a protein-only result."
                    ),
                )
                save(out, job)
                print(job["message"], flush=True)
                return 2
            array = trusted[0]
            spacer_path = Path(array["spacer_file"])
            fasta(spacer_path, "ACGTN")
            job["selected_array"] = array
            job["confirmed_spacer_orientation"] = args.spacer_orientation
            spacer_args = ["--spacers", str(spacer_path)]
            for p in target_paths:
                spacer_args += ["--targets", str(p)]
            if args.spacer_orientation == "reverse":
                spacer_args.append("--reverse-spacers")

        # --- inference ----------------------------------------------------
        protein_path = out / "selected_cas9.faa"
        protein_path.write_text(
            ">" + selected["system_id"] + "\n" + selected["protein_sequence"] + "\n",
            encoding="utf-8",
        )
        job["selected_cas9_faa"] = str(protein_path.relative_to(out))
        candidate_args = (
            ["--candidates", args.candidates] if args.candidates else ["--auto-candidates"]
        )
        job["candidate_mode"] = "explicit_user_candidates" if args.candidates else "model_consensus_only"
        device = args.device
        if device == "auto":
            device = "cuda" if _cuda_available() else "cpu"
        job["inference_device"] = device
        # Tests substitute a stub predictor via this variable so the full state
        # machine (resume, integrity, failure handling) can be exercised without
        # the multi-GB model weights. Unset in normal use, so production always
        # runs the real inference script.
        predictor_script = os.environ.get("PAMDICT_PREDICTOR_SCRIPT")
        if predictor_script:
            predictor_script = Path(predictor_script)
            if not predictor_script.is_file():
                raise RuntimeError(
                    f"PAMDICT_PREDICTOR_SCRIPT does not exist: {predictor_script}"
                )
        else:
            predictor_script = ROOT / "scripts/run_pam_report.py"
        command = [
            sys.executable, str(predictor_script),
            "--protein", str(protein_path),
            "--outdir", str(out / "pam"),
            "--device", device,
            "--timeout", str(args.inference_timeout),
            *candidate_args, *spacer_args,
        ]
        run_stage(command, out, job, "predicting", args.inference_timeout * 2 + 120)

        _summarise_pam_result(out / "pam", job)
        job.update(
            status="complete",
            report="pam/report.html",
            evidence_mode="protein_and_spacer" if spacer_args else "protein_only",
        )
        save(out, job)
        print(f"完成：{out / 'index.html'}")
        return 0
    except TimeoutError as exc:
        job.update(status="failed", error=str(exc),
                   message="Timeout: the stage was terminated with its process group. No completed result.")
        save(out, job)
        print(job["message"], file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001 - reported to the user, never swallowed
        job.update(
            status="failed", error=f"{type(exc).__name__}: {exc}",
            message="No completed result. Inspect the stage logs; existing outputs were preserved.",
        )
        save(out, job)
        print(f"{job['message']} {exc}", file=sys.stderr)
        return 1


def _cuda_available():
    try:
        import torch
        return bool(torch.cuda.is_available())
    except Exception:  # noqa: BLE001 - torch is optional for --device cpu
        return False


def _summarise_pam_result(pam_dir, job):
    """Record the model consensus, and say exactly how much of it was called.

    ``score_candidate_pams`` builds the consensus with the official rule: a
    position becomes ``N`` when no base reaches the confidence threshold, and
    trailing ``N`` positions are stripped for a downstream (Cas9) PAM.  ``N`` is
    a valid IUPAC symbol meaning "any base", so a consensus containing ``N`` is
    still a usable pattern and stays in the pipeline -- a canonical 3-nt Cas9
    consensus such as ``NGG`` has an undetermined position 1 yet is still the
    expected result and must still reach the spacer stage.

    Three outcomes are distinguished:

    ``no_clear_pam``
        No position was called at all.  Nothing usable was obtained, so this
        must never be shown as a biological PAM or as a successful candidate
        screen.  Only this case yields no candidates.
    ``partial_consensus``
        Some positions were called and others were not, so the consensus denotes
        a *set* of PAMs rather than one sequence.  ``determined_positions`` lists
        the positions carrying a concrete base call.  This is a normal outcome
        and is reported as scope, not as failure.
    ``model_consensus_only``
        Every retained position was called.  Still only model compatibility.
    """
    report_path = Path(pam_dir) / "report.json"
    if not report_path.is_file():
        raise RuntimeError("Inference finished without writing report.json")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    proteins = report.get("proteins", [])
    consensus = [str(p.get("predicted_pam") or "") for p in proteins]
    job["model_consensus"] = consensus

    # `determined_positions` is computed by score_candidate_pams.py using the same
    # probability threshold that built the consensus.  It is consumed from there
    # rather than re-derived here, because deriving it from the separate
    # information diagnostic (a 0.15-bit floor) silently disagreed with the
    # consensus rule (a 0.70 probability threshold) and overstated how many
    # positions were actually called.
    determined: list[list[int]] = []
    for entry in proteins:
        positions = entry.get("determined_positions")
        if positions is None:
            raise RuntimeError(
                "protein.json is missing determined_positions; it was written by an older "
                "score_candidate_pams.py and cannot be graded reliably"
            )
        determined.append([int(x) for x in positions])
    job["determined_positions"] = determined

    empty, partial = [], []
    for value, positions in zip(consensus, determined):
        if not value or not positions:
            empty.append(value or "(empty)")
        elif len(positions) < len(value):
            partial.append(value)

    if empty:
        job["pam_conclusion"] = "no_clear_pam"
        job["pam_conclusion_note"] = (
            "The model called no position for this protein (consensus "
            f"{empty}). This is reported as 'no clear PAM obtained', not as a biological "
            "PAM finding, and not as a successful candidate screen. Supply an explicit "
            "--candidates list to rank a user-chosen PAM set instead."
        )
        job.setdefault("warnings", []).append(job["pam_conclusion_note"])
    elif partial:
        job["pam_conclusion"] = "partial_consensus"
        job["pam_conclusion_note"] = (
            f"The consensus {partial} contains 'N' positions, so it denotes a set of PAMs "
            "rather than one sequence. Only the positions in determined_positions carry a "
            "concrete base call; an 'N' means any base and must not be read as one specific "
            "nucleotide. It is still scored as an IUPAC pattern. This is a scope statement, "
            "not a claim of failure."
        )
    else:
        job["pam_conclusion"] = "model_consensus_only"
        job["pam_conclusion_note"] = (
            "Every retained position was determined, but this reflects the model's own "
            "probability matrix only. It is a compatibility estimate, not cleavage activity, "
            "and it is not exhaustive PAM discovery or independent validation."
        )


if __name__ == "__main__":
    raise SystemExit(main())
