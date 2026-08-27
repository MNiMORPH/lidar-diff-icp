#!/usr/bin/env python3
"""Apply the NGV vegetation correction to gen2 and re-form the DoD.

gen2 was flown 2021-05-01 at green-up and reads HIGH under vegetation; gen1 (Nov 2008,
leaf-off) is taken as the reference (Andy's call, 2026-08-27). Measured against surveyed
control this session: NVA -2.2 mm, VVA -75.7 mm, i.e. gen2 reads high in vegetation by
73.4 mm.

    offset = surveyed_Z - delivered_Z          +ve = the delivered surface reads LOW
    offset = a + b * NGV                       a = datum term, b*NGV = vegetation term
    z_gen2_corrected = z_gen2 + b * NGV        the VEGETATION term only -- never `a`,
                                               which is gen2's own datum level
    DoD = z_gen2 - z_gen1  (+ve = deposition)
    DoD_corrected = DoD + b * NGV

`a` is deliberately excluded: it is a constant datum offset, not a vegetation response, and
folding it in would shift open ground where no vegetation bias exists.

    ./lidar-icp/bin/python analysis/ngv_correct_dod.py --tile data/derived/elba_fulldensity
"""
import argparse, json, os
import numpy as np

# From analysis/ngv.py --marks, refit on the EXACT index (commit 66b0561):
#   free intercept  a = -8.7 +/- 4.2 mm    b = -325.2 +/- 35.4 mm per unit NGV
#   block bootstrap                        b = -324.2 +/- 44.4
B_MM_PER_NGV = -325.2
B_SE_MM = 44.4          # the block-bootstrap SE, which respects the spatial clustering
CALIB_MAX_NGV = None    # measured below from the control table


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", default="data/derived/elba_fulldensity")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    T = a.tile

    ngv = np.load(os.path.join(T, "ngv.npy"))
    dod = np.load(os.path.join(T, "dod.npy"))
    cfg = json.load(open(os.path.join(T, "corrections.json")))
    assert cfg.get("swath_tie") == "intercept", \
        f"this DoD was not built with the adopted swath tie: swath_tie={cfg.get('swath_tie')!r}"

    import pandas as pd
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import control_lowveg_offset as M
    marks = M.load(0.15, 2.0).merge(
        pd.read_csv("data/derived/control_ngv_exact.csv")[["point_id", "ngv"]],
        on="point_id", how="inner")
    cmax = float(marks.ngv.max())

    fp = os.path.join(T, "floodplain_mask.npy")
    flood = np.load(fp).astype(bool) if os.path.exists(fp) else np.zeros(dod.shape, bool)

    print(f"tile {T}   grid {dod.shape}   swath_tie={cfg['swath_tie']}")
    print(f"correction: DoD_corrected = DoD + ({B_MM_PER_NGV:+.1f} mm) * NGV")
    print(f"  b from {len(marks)} control marks, exact NGV; block-bootstrap SE {B_SE_MM:.1f}")
    print(f"  the datum term a is NOT applied (it is not a vegetation response)\n")

    fin = np.isfinite(dod) & np.isfinite(ngv) & ~flood
    print(f"  cells: {dod.size:,} total, {int(np.isfinite(dod).sum()):,} with a DoD, "
          f"{int(flood.sum()):,} floodplain EXCLUDED (standing rule), {int(fin.sum()):,} used")

    print(f"\n  EXTRAPOLATION -- the control marks calibrate NGV only up to {cmax:.3f}")
    over = fin & (ngv > cmax)
    print(f"    cells above that: {int(over.sum()):,} ({100 * over.sum() / fin.sum():.2f}%)"
          f"   tile NGV max {np.nanmax(ngv[fin]):.3f}")
    for q in (50, 75, 90, 99, 100):
        print(f"    tile NGV p{q:<3d} {np.percentile(ngv[fin], q):.3f}"
              f"    marks NGV p{q:<3d} {np.percentile(marks.ngv, q):.3f}")

    dz = (B_MM_PER_NGV / 1000.0) * ngv
    corr = dod + dz

    print(f"\n  CORRECTION SIZE (mm), over used cells")
    print(f"    median {1000 * np.median(dz[fin]):+7.1f}   p90 {1000 * np.percentile(dz[fin], 10):+7.1f}"
          f"   max {1000 * np.min(dz[fin]):+7.1f}")
    print(f"    DoD median  before {1000 * np.median(dod[fin]):+7.1f}   "
          f"after {1000 * np.median(corr[fin]):+7.1f} mm")

    cov_p = os.path.join(T, "canopy_cover_pfs.npy")
    if os.path.exists(cov_p):
        cov = np.load(cov_p)
        print(f"\n  THE TEST -- DoD against canopy cover, before and after.")
        print(f"  If the cover trend was gen2's leaf-on bias it should FLATTEN; if it is")
        print(f"  differential erosion it should survive. cover from canopy_cover_pfs.npy,")
        print(f"  an INDEPENDENT measure (PyForestScan plant-area density), not NGV.")
        print(f"    {'cover':>12} {'n':>8} {'DoD before':>12} {'DoD after':>12} {'NGV med':>9}")
        edges = [0.0, 0.05, 0.15, 0.30, 0.50, 0.70, 1.01]
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = fin & np.isfinite(cov) & (cov >= lo) & (cov < hi)
            if k.sum() < 50:
                print(f"    {lo:.2f}-{hi:.2f} {int(k.sum()):8d}   -- too few")
                continue
            b_, a_ = 1000 * np.median(dod[k]), 1000 * np.median(corr[k])
            rows.append((0.5 * (lo + hi), b_, a_, int(k.sum())))
            print(f"    {lo:.2f}-{hi:.2f} {int(k.sum()):8d} {b_:+12.1f} {a_:+12.1f} "
                  f"{np.median(ngv[k]):9.3f}")
        if len(rows) > 2:
            x = np.array([r[0] for r in rows])
            sb = np.polyfit(x, [r[1] for r in rows], 1)[0]
            sa = np.polyfit(x, [r[2] for r in rows], 1)[0]
            print(f"    slope of DoD on cover:  before {sb:+8.1f}   after {sa:+8.1f} mm per unit cover")
            print(f"    |slope| reduced by {100 * (1 - abs(sa) / abs(sb)):.0f}%")

    for nm in ("core_forest", "core_open"):
        p = os.path.join(T, nm + ".npy")
        if not os.path.exists(p):
            continue
        k = fin & np.load(p).astype(bool)
        if k.sum():
            print(f"\n  {nm:12s} n={int(k.sum()):7,d}  DoD before {1000*np.median(dod[k]):+7.1f}"
                  f"   after {1000*np.median(corr[k]):+7.1f} mm   NGV med {np.median(ngv[k]):.3f}")

    f_p, o_p = os.path.join(T, "core_forest.npy"), os.path.join(T, "core_open.npy")
    if os.path.exists(f_p) and os.path.exists(o_p):
        kf = fin & np.load(f_p).astype(bool); ko = fin & np.load(o_p).astype(bool)
        if kf.sum() and ko.sum():
            gb = 1000 * (np.median(dod[kf]) - np.median(dod[ko]))
            ga = 1000 * (np.median(corr[kf]) - np.median(corr[ko]))
            print(f"\n  FOREST - OPEN gap: before {gb:+.1f} mm   after {ga:+.1f} mm"
                  f"   ({100*(1-abs(ga)/abs(gb)):.0f}% closed)")

    out = a.out or os.path.join(T, "dod_ngv.npy")
    np.save(out, corr)
    np.save(out.replace(".npy", "_dz.npy"), dz)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
