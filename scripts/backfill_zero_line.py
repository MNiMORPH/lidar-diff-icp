#!/usr/bin/env python3
"""Record each shipped corrections.json's ZERO LINE.

The ZERO LINE is the flight line defined as zero when a tile's swath network is solved.
It therefore carries a zero shift, and its own vertical error becomes the level the
whole tile inherits. It is arbitrary and per-tile, and an absolute datum cancels it
exactly -- but two products are only relatable if each says which line it was.
`difference_dem` writes `zero_line` now; files written before that do not carry it, so
you cannot tell from a product which line it used without re-deriving it, and
`difference_dem`'s guard (which refuses a datum constant measured against a different
zero line) has nothing to check against.

The zero line is NOT assumed from the rule. It is read from the constants themselves:
it is the swath whose shift is exactly zero on all three axes. That it also equals
`min(swaths)` -- the rule in pipeline.py:714 -- is then ASSERTED, so a file that
disagrees is reported rather than quietly relabelled.

Idempotent: a file that already records a matching value is left alone.

    ./lidar-icp/bin/python scripts/backfill_zero_line.py            # report only
    ./lidar-icp/bin/python scripts/backfill_zero_line.py --write
"""
import argparse, glob, json, os

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--glob", default="data/derived/*/corrections*.json")
ap.add_argument("--write", action="store_true", help="write the key (default: report only)")
A = ap.parse_args()

n_ok = n_written = n_skipped = n_problem = 0
for p in sorted(glob.glob(A.glob)):
    try:
        d = json.load(open(p))
    except Exception as e:
        print(f"  {p}: unreadable ({e})"); n_problem += 1; continue
    sw = d.get("per_swath_internal_alignment_dxdydz_m")
    if not sw:
        n_skipped += 1; continue
    items = {int(k): v for k, v in sw.items()}
    zero = [s for s, v in items.items() if all(abs(c) < 1e-12 for c in v[:3])]
    if len(zero) != 1:
        print(f"  {p}: {len(zero)} all-zero swaths {sorted(zero)} -- cannot identify the "
              f"zero line, NOT written")
        n_problem += 1; continue
    zero_line = zero[0]
    if zero_line != min(items):
        print(f"  {p}: zero line {zero_line} is NOT the lowest swath {min(items)} -- pipeline.py:714 "
              f"uses the lowest, so this file was made some other way. NOT written.")
        n_problem += 1; continue
    have = d.get("zero_line")
    if have is not None:
        if int(have) != zero_line:
            print(f"  {p}: records zero_line={have} but the constants say {zero_line}")
            n_problem += 1
        else:
            n_ok += 1
        continue
    print(f"  {p}: zero line {zero_line}  (swaths {sorted(items)})")
    if A.write:
        d["zero_line"] = zero_line
        with open(p, "w") as fh:
            json.dump(d, fh, indent=2)
            fh.write("\n")
        n_written += 1

print(f"\n  {n_ok} already correct, {n_written} written, {n_skipped} without swath constants, "
      f"{n_problem} needing attention")
if not A.write and n_problem == 0:
    print("  re-run with --write to record them")
