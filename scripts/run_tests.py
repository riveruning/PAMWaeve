#!/usr/bin/env python3
"""Offline regression suite for the submitted analysis package.

This is the *submission-package* test entry.  It intentionally differs from the
main repository's ``scripts/run_tests.py``, which also runs the local web
interface, the optional LLM assistant and the reproducible-run substrate.  Those
modules are **not shipped here** (see docs/SUBMISSION.md "what was excluded"),
so listing them would abort the run with ``ModuleNotFoundError``.

What this entry guarantees instead:

* every module it lists is actually present in this package, verified before
  the first test is imported, and a missing module is reported as an explicit
  failure naming the module rather than a bare traceback;
* the suite needs no model weights, no network and no CRISPRCasTyper;
* the exit status is non-zero if anything fails.

Usage:
    PYTHONPATH=src:$PYTHONPATH python3 scripts/run_tests.py
"""
from __future__ import annotations

import importlib
import importlib.util
import sys
import traceback
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(WORKSPACE / "src"))
sys.path.insert(0, str(WORKSPACE / "scripts"))

#: Modules using module-level ``test_*`` functions.
#: Every entry must exist in this package; ``verify_modules`` enforces that.
MODULES = [
    "pamdict.genome_tests.test_genome",
    "pamdict.score.tests.test_delivery",
    "pamdict.benchmark.tests.test_pammla_dataset",
    "pamdict.benchmark.tests.test_pammla_pilot",
    "pamdict.benchmark.tests.test_candidate_audits",
    "pamdict.benchmark.tests.test_multi_evidence",
    "pamdict.benchmark.tests.test_dual_evidence_cases",
    "pamdict.collect.tests.test_collect",
    "pamdict.finetune.tests.test_finetune",
    "pamdict.score.tests.test_score",
    "pamdict.score.tests.test_spacer",
]

#: Modules that the main repository runs but that are deliberately excluded
#: from this package.  Recorded (not silently dropped) so a reader can see the
#: coverage difference and knows those guarantees were verified in the main
#: repository rather than here.
EXCLUDED_MODULES = [
    "pamdict.app.tests.test_app_tasks",
    "pamdict.app.tests.test_app_submission",
    "pamdict.app.tests.test_app_security",
    "pamdict.app.tests.test_app_entries",
    "pamdict.app.tests.test_assistant",
    "pamdict.agent.tests.test_agent_tools",
    "pamdict.agent.tests.test_agent_review_fixes",
    "pamdict.repro.tests.test_repro",
]


def module_path(modname: str) -> Path | None:
    try:
        spec = importlib.util.find_spec(modname)
    except (ImportError, ValueError):
        return None
    return Path(spec.origin) if spec and spec.origin else None


def verify_modules() -> list[str]:
    """Fail fast, naming the module, if a listed module is not in this package."""
    problems: list[str] = []
    for modname in MODULES:
        path = module_path(modname)
        if path is None or not path.is_file():
            problems.append(
                f"{modname}: not found in this package (expected under "
                f"{WORKSPACE / 'src'}). This entry must only list modules that "
                f"are actually shipped."
            )
    return problems


def run_function_module(modname: str, passed: int, failed: int) -> tuple[int, int]:
    mod = importlib.import_module(modname)
    names = sorted(n for n in dir(mod) if n.startswith("test_"))
    if not names:
        print(f"WARN {modname}: no test_ functions found")
    for name in names:
        fn = getattr(mod, name)
        if not callable(fn):
            continue
        try:
            fn()
            passed += 1
            print(f"PASS {modname}.{name}")
        except Exception:  # noqa: BLE001 - a failing test must not stop the run
            failed += 1
            print(f"FAIL {modname}.{name}")
            traceback.print_exc()
    return passed, failed


def main() -> int:
    problems = verify_modules()
    if problems:
        print("SUITE CANNOT RUN: listed module(s) missing from this package.")
        for problem in problems:
            print(f"  - {problem}")
        return 2

    print("Modules excluded from this package (verified in the main repository; "
          "see docs/SUBMISSION.md):")
    for modname in EXCLUDED_MODULES:
        print(f"  - {modname}")
    print()

    passed = failed = 0
    for modname in MODULES:
        try:
            passed, failed = run_function_module(modname, passed, failed)
        except Exception:  # noqa: BLE001
            failed += 1
            print(f"FAIL {modname} (module import/collection error)")
            traceback.print_exc()

    print()
    print(f"--- {passed} passed, {failed} failed, {passed + failed} total ---")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
