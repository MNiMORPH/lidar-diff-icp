#!/usr/bin/env python3
"""Andy's stable-ground definition: ZERO-CURVATURE, LOW-SLOPE DRAINAGE DIVIDES, not floodplain.

Physical reasoning, which the packaged definition does not have: on a divide with zero
transverse curvature the diffusive flux divergence is zero, so no elevation change is
expected from hillslope transport. Low slope suppresses advective transport as well.
The packaged mask instead takes flat-or-CONVEX ground -- and convex divides are exactly
where diffusion predicts steady lowering, so it calls "stable" the places most expected
to erode.

    divides   = cells where kappa_L20.npy is finite (the S&S divide network)
    curvature = |kappa| < KTOL          transverse convexity over +/-20 m, m^-1-ish
    slope     = slope.npy < SMAX        degrees
    exclude   = floodplain_mask.npy

No threshold is chosen here. The ladder is printed so the cut is Andy's.

    ./lidar-icp/bin/python analysis/stable_flat_divides.py
"""
import argparse, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from detect_change_ngv import stable_mask, clip_stable, nmad


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", default="data/derived/elba_fulldensity")
    ap.add_argument("--save", default=None, help="write the mask for these cuts, 'KTOL,SMAX'")
    a = ap.parse_args()
    T = a.tile
    cfg = json.load(open(f"{T}/corrections.json")); res = float(cfg["res_m"])
    dod = np.load(f"{T}/dod.npy"); corr = np.load(f"{T}/dod_ngv.npy")
    kap = np.load(f"{T}/kappa_L20.npy"); slope = np.load(f"{T}/slope.npy")
    flood = np.load(f"{T}/floodplain_mask.npy").astype(bool)
    Z21 = np.load(f"{T}/z_after.npy")
    ngv = np.load(f"{T}/ngv.npy")

    divides = np.isfinite(kap)
    print(f"divide network (kappa_L20 finite): {int(divides.sum()):,} cells of {kap.size:,}")
    print(f"  of those, floodplain {int((divides & flood).sum()):,}, "
          f"with a DoD {int((divides & np.isfinite(dod)).sum()):,}")
    print(f"  kappa on divides: " + "  ".join(
        f"p{q}={np.percentile(kap[divides], q):+.4f}" for q in (1, 5, 25, 50, 75, 95, 99)))
    print(f"  slope on divides: " + "  ".join(
        f"p{q}={np.percentile(slope[divides], q):.1f}" for q in (5, 25, 50, 75, 95)))
    print(f"  the packaged crest_mask keeps only kappa > 0.004 (CONVEX): "
          f"{int(np.load(f'{T}/crest_mask.npy').sum()):,} cells\n")

    geom = stable_mask(Z21, res)
    st_pipe, _ = clip_stable(geom, dod)
    print("REFERENCE -- the packaged mask, on the same DoDs")
    print(f"  {'mask':34s} {'n':>7} {'med b':>8} {'NMAD b':>8} {'med a':>8} {'NMAD a':>8} {'NGV':>6}")
    k = st_pipe
    print(f"  {'pipeline (flat OR convex, clipped)':34s} {int(k.sum()):7d} "
          f"{1000*np.median(dod[k]):+8.1f} {1000*nmad(dod[k]):8.1f} "
          f"{1000*np.median(corr[k]):+8.1f} {1000*nmad(corr[k]):8.1f} "
          f"{np.nanmedian(ngv[k]):6.3f}")

    print("\nYOUR DEFINITION -- zero-curvature low-slope divides, floodplain excluded")
    print("  b = DoD before the NGV correction, a = after.  NO sigma-clip applied:")
    print("  the point of a physical criterion is that it does not need one.")
    print(f"  {'|kappa| <':>10} {'slope <':>8} {'n':>7} {'med b':>8} {'NMAD b':>8} "
          f"{'med a':>8} {'NMAD a':>8} {'NGV':>6}")
    for ktol in (0.0005, 0.001, 0.002, 0.004):
        for smax in (3.0, 5.0, 10.0):
            m = (divides & (np.abs(kap) < ktol) & (slope < smax) & ~flood
                 & np.isfinite(dod) & np.isfinite(corr))
            if m.sum() < 30:
                print(f"  {ktol:10.4f} {smax:8.1f} {int(m.sum()):7d}   -- too few to report")
                continue
            print(f"  {ktol:10.4f} {smax:8.1f} {int(m.sum()):7d} "
                  f"{1000*np.median(dod[m]):+8.1f} {1000*nmad(dod[m]):8.1f} "
                  f"{1000*np.median(corr[m]):+8.1f} {1000*nmad(corr[m]):8.1f} "
                  f"{np.nanmedian(ngv[m]):6.3f}")

    if a.save:
        ktol, smax = (float(v) for v in a.save.split(","))
        m = divides & (np.abs(kap) < ktol) & (slope < smax) & ~flood
        out = f"{T}/stable_flat_divides.npy"
        np.save(out, m)
        print(f"\nwrote {out}  ({int(m.sum()):,} cells, |kappa| < {ktol:g}, slope < {smax:g} deg)")


if __name__ == "__main__":
    main()
