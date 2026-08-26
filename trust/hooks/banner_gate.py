#!/usr/bin/env python3
"""PostToolUse hook (matcher: Bash): a script that printed a table must have printed
its provenance.

Closes: "labels that named nothing checkable" and "silent input substitution". It cannot
undo the run -- PostToolUse fires after the tool has executed -- but exit code 2 feeds
the message back to the model, so the correction happens in the same turn, before the
number is quoted into prose.

Fires only when ALL of:
  * the command invoked a .py under analysis/ or scripts/, and
  * the output looks like a results table (>= 3 lines carrying >= 2 numeric fields), and
  * the output carries no `== PROVENANCE v1 ==` banner.

Deliberately silent for one-liners, greps, file listings and pipeline runs.

Test: python trust/hooks/banner_gate.py --selftest
"""
from __future__ import annotations

import json
import os
import re
import sys

BANNER_START = "== PROVENANCE v1 =="
SCRIPT_RE = re.compile(r"(?:^|[\s/])(?:analysis|scripts)/\S*\.py\b")
ROW_RE = re.compile(r"(?:(?<![\w.])[+-]?\d[\d,]*(?:\.\d+)?(?![\w])[^\n]*){2,}")


def looks_like_table(out: str) -> bool:
    rows = [l for l in out.splitlines() if len(re.findall(r"(?<![\w.])[+-]?\d[\d,]*(?:\.\d+)?", l)) >= 2]
    return len(rows) >= 3


def verdict(command: str, out: str):
    if not SCRIPT_RE.search(command or ""):
        return None
    if BANNER_START in (out or ""):
        return None
    if not looks_like_table(out or ""):
        return None
    return (
        "TRUST GATE (provenance) -- this script printed a results table with no provenance "
        "banner, so its numbers carry no checkable label.\n"
        "Before quoting ANY number from it, wire it through trust/provenance.py:\n"
        "    from trust.provenance import Run\n"
        "    R = Run('<the question this run answers>')\n"
        "    R.input(path, role='<what the numbers in this file MEAN>')   # for every input\n"
        "    R.param(name, value, src='andy'|'repo'|'MINE')               # for every parameter\n"
        "    R.mask('stable', m, defn='<the cuts applied>')               # for every selection\n"
        "    R.column(name, '<definition, with units>')                   # for every column\n"
        "    R.banner(); R.table(header, rows); R.done()\n"
        "The banner is what makes the label come from the code that made the number, "
        "rather than from you afterwards. Andy has to be able to read the input paths "
        "next to the table, not take your word for which file was read."
    )


def main():
    if os.environ.get("TRUST_OFF") == "1" or os.environ.get("TRUST_LEVEL", "").lower() in ("off", "core"):
        return 0
    try:
        d = json.load(sys.stdin)
    except Exception:
        return 0
    cmd = (d.get("tool_input") or {}).get("command", "")
    res = d.get("tool_response") or d.get("toolUseResult") or {}
    out = res if isinstance(res, str) else (res.get("stdout", "") or "") + (res.get("stderr", "") or "")
    v = verdict(cmd, out)
    if v:
        sys.stderr.write(v + "\n")
        return 2
    return 0


def _selftest():
    tbl = "  bin   n    dz\n  a   12  -3.4\n  b   40  -8.1\n  c   11  -1.2\n"
    cases = [
        ("script + table + no banner -> block", "python analysis/x.py", tbl, True),
        ("script + table + banner -> pass", "python analysis/x.py", BANNER_START + "\n" + tbl, False),
        ("script + no table -> pass", "python analysis/x.py", "wrote figure\n", False),
        ("non-script one-liner -> pass", "wc -l data/*.npy", tbl, False),
        ("scripts/ dir counts", "./lidar-icp/bin/python scripts/y.py --a 1", tbl, True),
    ]
    ok = True
    for name, cmd, out, want in cases:
        got = verdict(cmd, out) is not None
        print(("OK  " if got == want else "FAIL") + f" {name}")
        ok &= got == want
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else main())
