#!/usr/bin/env python3
"""What already exists for this artifact? Read BEFORE building or changing anything.

Written 2026-09-01 after a session in which I repeatedly acted on a product without
reading the work behind it: built a vegetation index that duplicated q2(cover) and lost to
it, ran a fit at defaults that did not match its own documented recipe, and reinvented a
valley cut that was already established at 230 m. Each cost Andy a correction. The common
cause was not ignorance -- the evidence was on disk, and in two cases had already been
printed to my screen -- it was acting without looking.

So this is one command that answers "what is already known about X":

    ./lidar-icp/bin/python analysis/provenance_of.py floodplain_mask.npy
    ./lidar-icp/bin/python analysis/provenance_of.py q2 dod_cover_q2.npy

For each name it reports, in the order worth reading:
    DOCS      .md files that discuss it -- the recipe, the caveats, the choices made
    PRODUCED  code that WRITES it, with the commit that last touched that line
    CONSUMED  code that READS it
    CONSTANTS module-level constants whose name matches, with their comment

It greps; it does not judge. The point is that reading is one command instead of an
intention.
"""
import argparse, os, re, subprocess, sys

ROOT = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True,
                      text=True).stdout.strip() or "."
SKIP = (".git/", "lidar-icp/", "__pycache__", ".trust/runs/")


def tracked():
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True).stdout
    return [f for f in out.splitlines() if not any(s in f for s in SKIP)]


def last_touched(path):
    r = subprocess.run(["git", "log", "-1", "--format=%h %ad %s", "--date=short", "--", path],
                       cwd=ROOT, capture_output=True, text=True).stdout.strip()
    return r[:96] if r else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("names", nargs="+", help="artifact, product or concept, e.g. floodplain_mask.npy")
    ap.add_argument("--context", type=int, default=0, help="lines of context on doc hits")
    a = ap.parse_args()
    files = tracked()

    for name in a.names:
        stem = re.escape(os.path.basename(name))
        pat = re.compile(stem, re.I)
        docs, writes, reads, consts = [], [], [], []
        for f in files:
            p = os.path.join(ROOT, f)
            try:
                text = open(p, encoding="utf-8", errors="ignore").read()
            except (IsADirectoryError, OSError):
                continue
            if not pat.search(text):
                continue
            for i, ln in enumerate(text.splitlines(), 1):
                if not pat.search(ln):
                    continue
                if f.endswith(".md"):
                    docs.append((f, i, ln.strip()))
                elif f.endswith(".py"):
                    if re.search(r"np\.save|savez|to_parquet|\.to_csv|open\([^)]*['\"]w", ln):
                        writes.append((f, i, ln.strip()))
                    elif re.search(r"np\.load|read_table|read_csv|open\(", ln):
                        reads.append((f, i, ln.strip()))
                    elif re.match(r"[A-Z_]{3,}\s*=", ln.strip()):
                        consts.append((f, i, ln.strip()))

        print("=" * 78)
        print(f"{name}")
        print("=" * 78)
        for label, rows, limit in (("DOCS -- read these first", docs, 14),
                                   ("PRODUCED BY", writes, 8),
                                   ("CONSTANTS", consts, 8),
                                   ("CONSUMED BY", reads, 14)):
            if not rows:
                continue
            print(f"\n  {label}")
            seen = set()
            for f, i, ln in rows[:limit]:
                tag = last_touched(f) if f not in seen else ""
                seen.add(f)
                print(f"    {f}:{i}")
                print(f"        {ln[:110]}")
                if tag:
                    print(f"        last commit: {tag}")
            if len(rows) > limit:
                print(f"    ... and {len(rows) - limit} more")
        if not (docs or writes or reads or consts):
            print("\n  nothing found -- it may not exist yet, which is itself the answer")
        print()


if __name__ == "__main__":
    main()
