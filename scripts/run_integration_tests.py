#!/usr/bin/env python3
"""Integration test runner (needs the real protein2pam weights).

Kept separate from ``scripts/run_tests.py`` so the offline unit suite stays
runnable with no network and no checkpoints.

Usage:
    HF_HOME=$PWD/data/checkpoints/hf \\
    PYTHONPATH=src:$PYTHONPATH \\
    python3 scripts/run_integration_tests.py

The upstream Protein2PAM adapter is located via
``PAMPRIDICT_PROTEIN2PAM_HOME``, an importable ``protein2pam`` package on
``PYTHONPATH``, or the conventional ``.reference/Protein2PAM`` (see README).

Prerequisites are checked up front; a missing environment fails with one clear
message rather than a transformers stack trace. Nothing here skips silently.
"""

import importlib
import sys
import traceback
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE))
sys.path.insert(0, str(WORKSPACE / "src"))

MODULES = [
    "pamdict.genome_tests.test_genome_integration",
]


def main() -> int:
    # Fail fast with an actionable message if the environment is unusable.
    try:
        preflight = importlib.import_module(MODULES[0])
        preflight._require_model()
        preflight._require_real_cas9()
    except Exception as e:  # noqa: BLE001
        print(f"PREREQUISITE MISSING\n{e}\n")
        return 1

    passed = failed = 0
    for modname in MODULES:
        mod = importlib.import_module(modname)
        for name in sorted(n for n in dir(mod) if n.startswith("test_")):
            fn = getattr(mod, name)
            try:
                fn()
                passed += 1
                print(f"PASS {modname}.{name}")
            except Exception as e:  # noqa: BLE001
                failed += 1
                print(f"FAIL {modname}.{name}: {e!r}")
                traceback.print_exc()
    print(f"\n--- integration: {passed} passed, {failed} failed, {passed + failed} total ---")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
