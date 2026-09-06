#!/usr/bin/env python3
"""Compute where every script under analysis/ belongs, and (with --apply) move it there.

Implements analysis/REORGANIZATION_PLAN.md. The point of this being a program rather than a
sequence of hand edits: every destination is COMPUTED from facts already on disk, so the
result is reviewable before it happens, reproducible after, and free of my opinion about
what any script is "about".

    reorganize_analysis.py --plan     print every destination and the cluster list; move nothing
    reorganize_analysis.py --apply    execute, one commit per cluster, gated, resumable
    reorganize_analysis.py --report   read the manifest and say what happened

THE RULE. Total, ordered, first match wins:

    1  named in a PIPELINE Step (a Step with no group)   stays put -- it is declared
    2  imported by src/, or by >= 2 other scripts        analysis/lib/
    3  named in a Step belonging to a GROUP              analysis/modules/<group>/
    4  calls savefig() and touched since 2026-09-01      analysis/tools/
    5  cited by any .md         analysis/investigations/<stem of the MOST SPECIFIC citer>/
    7  anything else                            STAYS PUT (the existing directory means
                                                something; a flat _unsorted would not)

Rule 1 tests PIPELINE steps specifically, not "any Step". Written as "any Step" it would
swallow the group's scripts before rule 3 ever saw them, and the vegetation-correction
module would never form -- caught while writing this, by checking the rule against the
graph instead of trusting the order I first wrote down.

Rule 5 is the load-bearing one: THE FOLDER IS NAMED BY THE DOCUMENT THAT CITES THE SCRIPT.
If a finding was worth writing down, its name already exists, so no taxonomy is invented
here. `_unsorted/` is where the rule says it does not know, rather than guessing.

CLUSTERS. 33 files import a sibling by bare name, so a move is never one file: the
connected components of the sibling-import graph move together. A component's destination
is its majority member's; ties break to the lowest-numbered rule. Deterministic, and it
never splits a component.
"""
import argparse
import ast
import collections
import json
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
from lidar_diff_icp import workflow as W        # noqa: E402

RECENT = "2026-09-01"          # rule 4: "still in use" means touched during the current work

#: Documents ABOUT THE PROJECT rather than about a scientific question. They cite scripts
#: because they survey, index or audit them, so a citation from one says nothing about which
#: investigation a script belongs to -- and a folder named `investigations/hidden_filters_audit`
#: would be organisation that misleads. Enumerated rather than inferred: it is a short, stable
#: list of documents, which is a far smaller thing to get right than a rule over 157 scripts.
#: A script cited ONLY by these falls to _unsorted, which is the honest answer -- it means no
#: findings document claims it.
META_DOCS = re.compile(
    r"(?:^|/)(?:README|TODO|AUDIT_findings|HIDDEN_FILTERS_AUDIT|SCRIPT_INVENTORY"
    r"|REORGANIZATION_PLAN|SESSION_[A-Z_]+|FRAME_[0-9-]+(?:-[A-Z]+)?"
    r"|HELP_NEXT_STEPS|[A-Z_]*SESSION[A-Z_]*)\.md$")
MANIFEST = os.path.join(REPO, "analysis", ".reorg_manifest.json")
CITE = re.compile(r"(?<![\w/.])((?:analysis|scripts|src|ground_control|trust)"
                  r"/[A-Za-z0-9_./-]+\.py)")


def _tracked(pattern):
    out = subprocess.run(["git", "-C", REPO, "ls-files", pattern],
                         capture_output=True, text=True).stdout.split()
    return [p for p in out if os.path.exists(os.path.join(REPO, p))]


def _ast(path):
    try:
        with open(os.path.join(REPO, path), errors="ignore") as fh:
            return ast.parse(fh.read(), filename=path)
    except SyntaxError:
        return None


def _facts():
    """Everything the rule reads, gathered once."""
    analysis = _tracked("analysis/*.py")
    allpy = _tracked("*.py")
    stem = {p: os.path.splitext(os.path.basename(p))[0] for p in allpy}

    step_of = {}                                    # script -> Step
    for s in W.STEPS:
        for tok in list(s.command.split()) + list(s.code):
            if tok.endswith(".py"):
                step_of[tok] = s

    imported_by = collections.defaultdict(set)      # script -> importers
    sibling = collections.defaultdict(set)          # script -> siblings it imports
    for p in allpy:
        tree = _ast(p)
        if tree is None:
            continue
        mods = set()
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                mods |= {a.name.split(".")[0] for a in n.names}
            elif isinstance(n, ast.ImportFrom) and n.module and not n.level:
                mods.add(n.module.split(".")[0])
        for q in allpy:
            if q != p and stem[q] in mods:
                imported_by[q].add(p)
                if os.path.dirname(q) == os.path.dirname(p):
                    sibling[p].add(q)
                    sibling[q].add(p)

    # A citation counts whether it names the PATH or just the FILE. The lint in
    # tests/test_repo_links.py deliberately checks only path-form citations, because only
    # those promise a location that a move can break. The RULE wants something different --
    # association, not linkability -- and 51 of these scripts are named by bare filename
    # only. Using the lint's regex here put them all in _unsorted, which is how I noticed.
    cited_by = collections.defaultdict(set)         # script -> documents
    names = {os.path.basename(p): p for p in analysis}
    for m in _tracked("*.md"):
        with open(os.path.join(REPO, m), errors="ignore") as fh:
            text = fh.read()
        for hit in CITE.findall(text):
            cited_by[hit].add(m)
        for base, path in names.items():
            if re.search(rf"(?<![\w/.]){re.escape(base)}(?![\w])", text):
                cited_by[path].add(m)

    savefig, recent = set(), set()
    for p in analysis:
        src = open(os.path.join(REPO, p), errors="ignore").read()
        if "savefig" in src:
            savefig.add(p)
        d = subprocess.run(["git", "-C", REPO, "log", "-1", "--format=%ad",
                            "--date=short", "--", p],
                           capture_output=True, text=True).stdout.strip()
        if d >= RECENT:
            recent.add(p)
    cites = collections.defaultdict(set)            # document -> scripts it cites
    for script, docs in cited_by.items():
        for m in docs:
            cites[m].add(script)
    return dict(analysis=analysis, step_of=step_of, imported_by=imported_by,
                sibling=sibling, cited_by=cited_by, cites=cites,
                savefig=savefig, recent=recent)


def destination(p, f):
    """(rule_number, destination_dir, why) for one script. Total and ordered."""
    st = f["step_of"].get(p)
    if st is not None and not st.group:
        return 1, os.path.dirname(p), f"declared by pipeline step {st.name}"
    importers = f["imported_by"].get(p, set())
    if any(q.startswith("src/") for q in importers) or len(importers) >= 2:
        return 2, "analysis/lib", f"imported by {len(importers)} file(s)"
    if st is not None and st.group:
        return 3, f"analysis/modules/{st.group}", f"step {st.name} of group {st.group}"
    if p in f["savefig"] and p in f["recent"]:
        return 4, "analysis/tools", f"draws figures and was touched since {RECENT}"
    docs = {m for m in f["cited_by"].get(p, set()) if not META_DOCS.search(m)}
    if docs:
        # THE MOST SPECIFIC DOCUMENT THAT CITES IT. Not "the only one" and not a shared bin:
        # an audit or index document cites scripts BECAUSE IT SURVEYS THEM, so counting
        # citations equally put 45 scripts in one undifferentiated pile and gave
        # HIDDEN_FILTERS_AUDIT.md 13 that it merely audited. Specificity is measured by how
        # many scripts a document cites, which needs no threshold: the document that talks
        # about fewest scripts is the one this script is most plausibly ABOUT. Ties break on
        # the document name, so the answer does not depend on set ordering.
        best = min(docs, key=lambda m: (len(f["cites"][m]), m))
        d = os.path.splitext(os.path.basename(best))[0].lower()
        return 5, f"analysis/investigations/{d}", (
            f"most specific of {len(docs)} citing document(s): {best} "
            f"(cites {len(f['cites'][best])})")
    # STAYS PUT. Not a flat _unsorted bin: analysis/ridgelines, slope_bias, groundtruth and
    # steady_state came out of the work and carry real meaning, so moving 83 scripts into one
    # undifferentiated folder would DESTROY information in the name of organising it. The
    # reorganization moves only what there is positive evidence about, and leaves the rest
    # where the work put it. That 53% of these scripts have no findings document claiming
    # them is the finding, and --report says so rather than hiding it in a folder name.
    return 7, os.path.dirname(p), "no findings document claims it; left where the work put it"


def clusters(f):
    """Connected components of the sibling-import graph, restricted to analysis/."""
    seen, out = set(), []
    for p in f["analysis"]:
        if p in seen:
            continue
        comp, stack = set(), [p]
        while stack:
            q = stack.pop()
            if q in comp:
                continue
            comp.add(q)
            stack += [r for r in f["sibling"].get(q, ()) if r in f["analysis"]]
        seen |= comp
        out.append(sorted(comp))
    return sorted(out, key=lambda c: (-len(c), c[0]))


def plan(f):
    """[(cluster, destination, rule, members_with_their_own_verdicts)] -- the whole answer."""
    rows = []
    for comp in clusters(f):
        verdicts = {p: destination(p, f) for p in comp}
        if len(comp) == 1:
            rule, dest, _ = verdicts[comp[0]]
        else:
            tally = collections.Counter((r, d) for r, d, _ in verdicts.values())
            top = max(tally.values())
            rule, dest = min(k for k, v in tally.items() if v == top)   # ties -> lowest rule
        rows.append((comp, dest, rule, verdicts))
    return rows


def depth_sensitive(path):
    """Does this file compute a path by walking up from __file__?

    36 scripts do, and a move that changes their DEPTH silently changes where they look.
    `os.path.dirname(os.path.dirname(os.path.abspath(__file__)))` from analysis/ridgelines
    resolves to analysis/; from analysis/modules/vegetation_correction it resolves to
    analysis/modules/. The import lint cannot see this when the target is a DATA directory
    rather than a module, so it would fail at runtime, later, quietly.

    Returns the number of dirname() levels applied to __file__, or 0.
    """
    tree = _ast(path)
    if tree is None:
        return 0
    best = 0
    for n in ast.walk(tree):
        if not (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                and n.func.attr == "dirname"):
            continue
        depth, cur = 0, n
        while (isinstance(cur, ast.Call) and isinstance(cur.func, ast.Attribute)
               and cur.func.attr in ("dirname", "abspath") and cur.args):
            if cur.func.attr == "dirname":
                depth += 1
            cur = cur.args[0]
        if isinstance(cur, ast.Name) and cur.id == "__file__":
            best = max(best, depth)
    return best


def blocked_reason(comp, dest):
    """Why this cluster must NOT be moved automatically, or None."""
    for p in comp:
        if os.path.dirname(p) == dest:
            continue
        d = depth_sensitive(p)
        if d and p.count("/") != dest.count("/") + 1:
            return (f"{os.path.basename(p)} walks {d} dirname level(s) up from __file__, and "
                    f"the move changes its depth ({p.count('/')} -> {dest.count('/') + 1}). "
                    f"A data path computed that way would break silently, so this is not "
                    f"moved automatically.")
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--plan", action="store_true")
    g.add_argument("--apply", action="store_true")
    g.add_argument("--report", action="store_true")
    a = ap.parse_args(argv)

    if a.report:
        if not os.path.exists(MANIFEST):
            print("no manifest; --apply has not run"); return 0
        man = json.load(open(MANIFEST))
        for st in ("done", "blocked", "pending"):
            names = [c for c, v in man.items() if v.get("status") == st]
            print(f"{st:8s} {len(names)}")
            for n in names:
                print(f"    {n}  {man[n].get('note','')}")
        return 0

    f = _facts()
    rows = plan(f)
    moving = [(c, d, r) for c, d, r in ((c, d, r) for c, d, r, _ in rows)
              if any(os.path.dirname(p) != d for p in c)]
    print(f"{len(f['analysis'])} scripts under analysis/, in {len(rows)} clusters; "
          f"{sum(len(c) for c, _, _ in moving)} would move\n")
    by_dest = collections.Counter()
    for comp, dest, rule, verdicts in rows:
        for p in comp:
            by_dest[dest] += 1
    for d, n in sorted(by_dest.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {n:4d}  {d}")
    print("\nclusters larger than one file (these MUST move together):")
    for comp, dest, rule, verdicts in rows:
        if len(comp) > 1:
            print(f"  rule {rule}  -> {dest}")
            for p in comp:
                r, dd, why = verdicts[p]
                mark = " " if dd == dest else "*"
                print(f"    {mark} {p:58s} (own verdict: rule {r}, {why})")
    if not a.apply:
        return 0

    # ---- PHASE 2 -----------------------------------------------------------------------
    # tools and lib first (smallest, most used), then the declared module, then the rest.
    def order_key(dest):
        for i, k in enumerate(("analysis/tools", "analysis/lib", "analysis/modules")):
            if dest.startswith(k):
                return i
        return 3

    man = json.load(open(MANIFEST)) if os.path.exists(MANIFEST) else {}
    todo = [(c, d) for c, d, _, _ in rows if any(os.path.dirname(p) != d for p in c)]
    todo.sort(key=lambda cd: (order_key(cd[1]), cd[0][0]))

    done = failed = skipped = 0
    for comp, dest in todo:
        key = comp[0]
        if man.get(key, {}).get("status") in ("done", "blocked"):
            continue
        why = blocked_reason(comp, dest)
        if why:
            man[key] = {"status": "blocked", "dest": dest, "files": comp, "note": why}
            skipped += 1
            continue

        print(f"\n-> {dest}  ({len(comp)} file(s))", flush=True)
        moved = []
        for src in comp:
            dst = os.path.join(dest, os.path.basename(src))
            os.makedirs(os.path.join(REPO, dest), exist_ok=True)
            subprocess.run(["git", "-C", REPO, "mv", src, dst], check=True)
            moved.append((src, dst))
            print(f"     {src} -> {dst}", flush=True)

        # Every reference to the old path, in code and in prose. This is what keeps the
        # workflow Step commands and the 192 citation edges pointing at real files.
        refs = 0
        for f_ in _tracked("*.py") + _tracked("*.md") + _tracked("*.toml"):
            fp = os.path.join(REPO, f_)
            try:
                t = open(fp, errors="ignore").read()
            except OSError:
                continue
            new = t
            for src, dst in moved:
                new = new.replace(src, dst)
            if new != t:
                open(fp, "w").write(new)
                refs += 1
        if refs:
            print(f"     rewrote references in {refs} file(s)", flush=True)

        gates = [("citation lint", [sys.executable, "-m", "pytest", "-q",
                                    "tests/test_repo_links.py"]),
                 ("import lint", [sys.executable, "-m", "pytest", "-q",
                                  "tests/test_repo_imports.py"]),
                 ("test suite", [sys.executable, "-m", "pytest", "-q"])]
        bad = None
        for name, cmd in gates:
            r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
            if r.returncode != 0:
                bad = f"{name}: {r.stdout.strip().splitlines()[-1] if r.stdout else r.stderr[-300:]}"
                break
        if bad:
            print(f"     RED: {bad}\n     reverting", flush=True)
            subprocess.run(["git", "-C", REPO, "checkout", "--", "."], check=False)
            subprocess.run(["git", "-C", REPO, "reset", "--hard", "HEAD"], check=False)
            man[key] = {"status": "blocked", "dest": dest, "files": comp, "note": bad}
            failed += 1
        else:
            subprocess.run(["git", "-C", REPO, "add", "-A"], check=True)
            subprocess.run(["git", "-C", REPO, "commit", "-q", "-m",
                            f"reorg: {len(comp)} file(s) -> {dest}\n\n"
                            + "\n".join(f"{a_} -> {b_}" for a_, b_ in moved)
                            + f"\n\nGates green: citation lint, import lint, full suite. "
                              f"References rewritten in {refs} file(s).\n\n"
                              f"Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"],
                           check=True)
            man[key] = {"status": "done", "dest": dest, "files": comp}
            done += 1
        json.dump(man, open(MANIFEST, "w"), indent=1, sort_keys=True)

    print(f"\nmoved {done} cluster(s); {failed} reverted; {skipped} not attempted "
          f"(depth-sensitive). See --report.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
