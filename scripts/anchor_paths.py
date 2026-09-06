#!/usr/bin/env python3
"""Make __file__-relative sys.path walks depth-independent, so the file can be moved.

28 scripts reach the repo root, src/ or analysis/ by counting directories upward:

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

which silently means something different after a move. The import lint cannot follow a
COMPUTED path (it reads literal strings only), so a break would surface at runtime, later.
That is why the reorganization refused to move them.

This resolves each walk AS IT IS NOW, expresses the same absolute target relative to the
repository root, and rewrites it against an anchor found by walking up to pyproject.toml.
Semantics are preserved exactly; only the dependence on depth is removed.

    anchor_paths.py --list      show what each walk resolves to; change nothing
    anchor_paths.py --apply     rewrite, one commit per file, gated on the full suite
"""
import argparse
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANCHOR = ('_REPO = os.path.dirname(os.path.abspath(__file__))\n'
          'while _REPO != "/" and not os.path.exists(os.path.join(_REPO, "pyproject.toml")):\n'
          '    _REPO = os.path.dirname(_REPO)   # depth-independent: find the repo root\n')

#: sys.path.insert(0, <expr>) where <expr> walks up from __file__.
WALK = re.compile(
    r"sys\.path\.(?:insert\(\s*\d+\s*,|append\()\s*("
    r"os\.path\.(?:abspath|join|dirname)\((?:[^()]|\([^()]*\))*\)"
    r")\s*\)")


def _resolve(expr, path):
    """What this expression evaluates to for a file at `path`, as an absolute path."""
    env = {"os": os, "__file__": os.path.join(REPO, path)}
    try:
        return os.path.normpath(eval(expr, env))       # noqa: S307 -- our own source
    except Exception:
        return None


def _tracked():
    out = subprocess.run(["git", "-C", REPO, "ls-files", "analysis/*.py"],
                         capture_output=True, text=True).stdout.split()
    return [p for p in out if os.path.exists(os.path.join(REPO, p))]


def rewrite(path, dry=True):
    """Return (new_source, [(expr, target_relpath)]) or (None, []) if nothing to do."""
    src = open(os.path.join(REPO, path)).read()
    hits = []
    for m in WALK.finditer(src):
        expr = m.group(1)
        if "__file__" not in expr:
            continue
        target = _resolve(expr, path)
        if target is None or not target.startswith(REPO):
            continue                                    # outside the repo: leave alone
        rel = os.path.relpath(target, REPO)
        if rel == ".":
            new = "_REPO"
        else:
            new = "os.path.join(_REPO, " + ", ".join(f'"{c}"' for c in rel.split("/")) + ")"
        hits.append((expr, rel, new))
    if not hits:
        return None, []
    out = src
    for expr, _, new in hits:
        out = out.replace(expr, new)
    if "_REPO = os.path.dirname" not in out:            # inject the anchor before first use
        i = out.index("sys.path.")
        j = out.rindex("\n", 0, i) + 1
        out = out[:j] + ANCHOR + out[j:]
    return out, hits


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--list", action="store_true")
    g.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)

    todo = []
    for p in _tracked():
        new, hits = rewrite(p)
        if hits:
            todo.append((p, new, hits))
    print(f"{len(todo)} file(s) with a depth-dependent sys.path walk\n")
    for p, _, hits in todo:
        print(f"  {p}")
        for expr, rel, _ in hits:
            print(f"      -> {rel or '<repo root>'}")
    if a.list:
        return 0

    done = failed = 0
    for p, new, hits in todo:
        open(os.path.join(REPO, p), "w").write(new)
        r = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=REPO,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  RED on {p}; reverting")
            subprocess.run(["git", "-C", REPO, "checkout", "--", p], check=True)
            failed += 1
            continue
        subprocess.run(["git", "-C", REPO, "add", p], check=True)
        subprocess.run(["git", "-C", REPO, "commit", "-q", "-m",
                        f"anchor: {p} no longer depends on its depth\n\n"
                        + "\n".join(f"{e}  ->  {rel or '<repo root>'}" for e, rel, _ in hits)
                        + "\n\nSame absolute target, expressed against a repo root found by "
                          "walking up to pyproject.toml. Full suite green.\n\n"
                          "Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"],
                       check=True)
        done += 1
    print(f"\nanchored {done}; {failed} reverted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
