#!/usr/bin/env python3
"""Which recorded runs read a file that has since changed on disk?

This is the cheap replacement for a dependency graph: the run ledger already recorded a
digest of every input at the moment the numbers were produced, so staleness is a
comparison, not a build system. The Stop hook runs this check automatically for runs
quoted in the current session; this script sweeps the whole ledger directory.

    python trust/check_ledgers.py [ledger_dir]
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from provenance import _digest  # noqa: E402

d = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("TRUST_LEDGER_DIR", ".trust/runs")
if not os.path.isdir(d):
    print(f"no ledger directory at {d} -- no runs have been recorded yet")
    sys.exit(0)
stale = fresh = gone = 0
for fn in sorted(os.listdir(d)):
    if not fn.endswith(".json"):
        continue
    led = json.load(open(os.path.join(d, fn)))
    for i in led.get("inputs", []):
        if not os.path.exists(i["path"]):
            print(f"GONE  {fn}: {i['path']}")
            gone += 1
        elif _digest(i["path"]) != i["digest"]:
            print(f"STALE {fn}: {i['path']}")
            print(f"        run at {led['started']} -- {led['question'][:70]}")
            stale += 1
        else:
            fresh += 1
print(f"\n{fresh} input(s) unchanged, {stale} changed since the run, {gone} missing")
sys.exit(1 if (stale or gone) else 0)
