#!/usr/bin/env python3
"""Capture (or check) the DoD each route produces, so a refactor can be proven inert.

The correction chain is load-bearing and mature. Any restructuring of it must change
NOTHING numerically, and "the tests still pass" is a weaker claim than "the arrays are
identical". This writes a baseline before a refactor and compares after it.

    python scripts/route_baseline.py --write     # before
    python scripts/route_baseline.py             # after; exits non-zero on any drift
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tests"))
from test_pipeline import _ground, _bump, _write_laz14, BOUNDS, X0, Y0, W  # noqa: E402
from lidar_diff_icp.pipeline import difference_dem                          # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "tests", "fixtures",
                   "route_baseline.npz")


def build(tmp):
    rng = np.random.default_rng(7)
    n = int(3.0 * 120 * W)
    x1 = rng.uniform(X0, X0 + 120, n); y1 = rng.uniform(Y0, Y0 + W, n)
    x2 = rng.uniform(X0 + 80, X0 + W, n); y2 = rng.uniform(Y0, Y0 + W, n)
    xb = np.concatenate([x1, x2]); yb = np.concatenate([y1, y2])
    ps = np.concatenate([np.ones(n), np.full(n, 2)])
    zb = _ground(xb, yb) + rng.normal(0, 0.02, len(xb))
    _write_laz14(os.path.join(tmp, "before.laz"), xb, yb, zb, ps, yb, np.zeros(len(xb)))
    na = int(4.0 * W * W)
    xa = rng.uniform(X0, X0 + W, na); ya = rng.uniform(Y0, Y0 + W, na)
    za = _ground(xa, ya) + _bump(xa, ya) + rng.normal(0, 0.02, na)
    _write_laz14(os.path.join(tmp, "after.laz"), xa, ya, za, np.ones(na), ya, np.zeros(na))
    return os.path.join(tmp, "before.laz"), os.path.join(tmp, "after.laz")


def run(before, after):
    kw = dict(res=5.0, ground_q=0.10, ground="low_q", ground_source="last_return",
              after_ground="last_return", valley_top_m=-1e9)
    out = {}
    for route, extra in (("independent", dict(geoid_datum=(0.0, 0.0, 0.0))),
                         ("delong", {})):
        r = difference_dem(before, after, BOUNDS, route=route, **kw, **extra)
        out[f"{route}__dod"] = r["dod"]
        out[f"{route}__lod"] = r["lod"]
        out[f"{route}__z_after"] = r["z_after"]
        out[f"{route}__stable"] = r["stable"].astype(float)
    # every combination the fork must keep reachable, not just the two named routes
    for cs in (False, True):
        for dr in (False, True):
            r = difference_dem(before, after, BOUNDS, route="independent",
                               correction_surface=cs, along_track_drift=dr,
                               geoid_datum=(0.0, 0.0, 0.0), **kw)
            out[f"cs{int(cs)}_dr{int(dr)}__dod"] = r["dod"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    import tempfile
    tmp = tempfile.mkdtemp()
    got = run(*build(tmp))
    if a.write:
        np.savez_compressed(OUT, **got)
        print(f"wrote {OUT} with {len(got)} arrays")
        return
    ref = np.load(OUT)
    bad = 0
    for k in sorted(got):
        d = got[k] - ref[k]
        m = np.isfinite(d)
        worst = float(np.nanmax(np.abs(d[m]))) if m.any() else 0.0
        same_nan = np.array_equal(np.isnan(got[k]), np.isnan(ref[k]))
        flag = "" if (worst == 0.0 and same_nan) else "   <-- CHANGED"
        bad += bool(flag)
        print(f"  {k:<28} max|diff| {worst:.3e} m   nan-pattern "
              f"{'same' if same_nan else 'DIFFERENT'}{flag}")
    print(f"\n{'IDENTICAL' if not bad else f'{bad} array(s) CHANGED'}")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
