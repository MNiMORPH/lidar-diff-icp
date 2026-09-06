#!/usr/bin/env python3
"""Emit the repository's import graph and citation graph as JSON.

PHASE 0c OF analysis/REORGANIZATION_PLAN.md. Phase 2 moves 143 files and rewrites the
citations that name them. This is the artifact that makes "did the move change anything?"
answerable rather than argued: the graphs are recorded BEFORE, and re-emitted after, and
the only thing allowed to differ is a path.

The graph is keyed by BASENAME on purpose. A move changes a script's directory, so a
path-keyed graph would differ everywhere and prove nothing; a basename-keyed one is
invariant under exactly the operation Phase 2 performs, so any difference is a real
change in structure.

    ./lidar-icp/bin/python scripts/repo_graph.py > analysis/repo_graph.json
"""
import ast
import json
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CITE = re.compile(r"(?<![\w/.])((?:analysis|scripts|src|ground_control|trust)"
                  r"/[A-Za-z0-9_./-]+\.py)")


def tracked(pattern):
    out = subprocess.run(["git", "-C", REPO, "ls-files", pattern],
                         capture_output=True, text=True).stdout.split()
    return [p for p in out if os.path.exists(os.path.join(REPO, p))]


def imports_of(path):
    try:
        with open(os.path.join(REPO, path), errors="ignore") as fh:
            tree = ast.parse(fh.read(), filename=path)
    except SyntaxError:
        return []
    mods = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            mods |= {a.name.split(".")[0] for a in n.names}
        elif isinstance(n, ast.ImportFrom) and n.module and not n.level:
            mods.add(n.module.split(".")[0])
    return sorted(mods)


def main():
    pys = tracked("*.py")
    local = {os.path.splitext(os.path.basename(p))[0] for p in pys}

    # SIBLING IMPORTS: the edges that pin a file to its directory, and so the edges that
    # decide which files must move together in Phase 2.
    sibling = {}
    for p in pys:
        here = os.path.dirname(p)
        hits = [m for m in imports_of(p)
                if m in local and os.path.exists(os.path.join(REPO, here, m + ".py"))]
        if hits:
            sibling[os.path.basename(p)] = sorted(hits)

    citations = {}
    for m in tracked("*.md"):
        with open(os.path.join(REPO, m), errors="ignore") as fh:
            hits = sorted({os.path.basename(h) for h in CITE.findall(fh.read())})
        if hits:
            citations[os.path.basename(m)] = hits

    json.dump({"n_py": len(pys), "n_md": len(tracked("*.md")),
               "sibling_imports": sibling, "citations": citations},
              sys.stdout, indent=1, sort_keys=True)
    print()


if __name__ == "__main__":
    main()
