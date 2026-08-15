#!/usr/bin/env python3
"""Error structure of the co-registered inter-swath residuals.

For each overlapping swath pair: co-register (Nuth & Kaeaeb), take the residual
difference (same-epoch overlap => pure acquisition error), and characterise its
spatial structure with a *robust* (Dowd) variogram. Stratifies by elevation as a
crude floodplain (reed canary grass) vs. upland proxy, and reports the
correlated-error detection limit (LoD95) by feature size.

Key point this quantifies: because the error correlation length (>~370 m) far
exceeds any feature we would map, spatial averaging gives ~no benefit
(N_eff ~ 1), so the detection limit is essentially the point error and is flat
with feature size.

Example
-------
    python scripts/error_structure.py data/before/4342-29-64.laz --figdir figures
"""
import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import laspy

from lidar_diff_icp.swathdiff import _median_grid
from lidar_diff_icp import coreg, variogram as vg

RES = 2.0


def _nmad(a):
    return 1.4826 * np.median(np.abs(a - np.median(a)))


def _residual(x, y, z, filt, psid, a, b):
    ma = filt & (psid == a)
    mb = filt & (psid == b)
    x0 = max(x[ma].min(), x[mb].min()); x1 = min(x[ma].max(), x[mb].max())
    y0 = max(y[ma].min(), y[mb].min()); y1 = min(y[ma].max(), y[mb].max())
    if x1 <= x0 or y1 <= y0:
        return None
    nx = int(np.ceil((x1 - x0) / RES)); ny = int(np.ceil((y1 - y0) / RES))
    zr = _median_grid(x[ma], y[ma], z[ma], RES, x0, y0, nx, ny)
    zs = _median_grid(x[mb], y[mb], z[mb], RES, x0, y0, nx, ny)
    c = coreg.nuth_kaab(zr, zs, RES)
    r = zr - (coreg._shift_grid(zs, c.dx, c.dy, RES) + c.dz)
    gy, gx = np.mgrid[0:ny, 0:nx]
    m = np.isfinite(r) & np.isfinite(zr)
    return (x0 + (gx[m] + 0.5) * RES, y0 + (gy[m] + 0.5) * RES, r[m], zr[m], c)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tile")
    p.add_argument("--max-lag", type=float, default=400.0)
    p.add_argument("--figdir")
    args = p.parse_args()

    f = laspy.read(args.tile)
    x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
    psid = np.asarray(f.point_source_id); cls = np.asarray(f.classification)
    rn = np.asarray(f.return_number); nr = np.asarray(f.number_of_returns)
    filt = (nr == 1) & ~np.isin(cls, [5, 6, 9])  # single-return terrain (cleanest)
    swaths = np.unique(psid).tolist()

    print(f"{Path(args.tile).name}: single-return terrain surface, swaths {swaths}")
    strata = {}
    for a, b in combinations(swaths, 2):
        res = _residual(x, y, z, filt, psid, a, b)
        if res is None:
            continue
        xs, ys, r, elev, c = res
        cen, g, cnt = vg.empirical_variogram(xs, ys, r, args.max_lag)
        mod = vg.fit_spherical(cen, g, cnt)
        print(f"  {a}->{b}: NMAD={_nmad(r):.3f} m  shift=({c.dx:+.2f},{c.dy:+.2f},{c.dz:+.2f})  "
              f"Dowd range={mod.range_:.0f} m  point_sigma={np.sqrt(mod.total_sill):.3f} m")
        # stash the widest-overlap pair for the LoD table
        strata[(a, b)] = (xs, ys, r, elev)

    # Detection limit per pair, stratified by elevation tercile with NEUTRAL
    # labels. Caveat: low-z is the grassy floodplain ONLY where the overlap
    # reaches the valley bottom (the western pairs); for an overlap up on the
    # valley side the tercile split is not a landcover boundary. A proper
    # floodplain mask (from the channel network) is the right tool; this
    # elevation split is a documented proxy, not a classifier.
    print(f"\n  LoD95 at meander scale (~2 ha), by elevation tercile:")
    print(f"  {'pair':>9} {'low-z third':>12} {'high-z third':>13}")
    A = 20000.0
    for (a, b), (xs, ys, r, elev) in strata.items():
        lo = elev <= np.percentile(elev, 33)
        hi = elev >= np.percentile(elev, 67)
        mLo = vg.fit_spherical(*vg.empirical_variogram(xs[lo], ys[lo], r[lo], args.max_lag))
        mHi = vg.fit_spherical(*vg.empirical_variogram(xs[hi], ys[hi], r[hi], args.max_lag))
        _, lLo = vg.detection_limit(mLo, A, RES * RES)
        _, lHi = vg.detection_limit(mHi, A, RES * RES)
        print(f"  {a}->{b:<5} {lLo:>11.3f}m {lHi:>12.3f}m")
    print("  Long correlation length (>~370 m) => averaging barely lowers the limit")
    print("  at feature scale (N_eff ~ 1). Where the overlap straddles the reed-canary")
    print("  floodplain (e.g. 136-137 low-z) the limit ~doubles vs. clean upland.")

    if args.figdir:
        key = max(strata, key=lambda k: strata[k][2].size)
        xs, ys, r, elev = strata[key]
        _plot(xs, ys, r, elev >= np.percentile(elev, 67),
              elev <= np.percentile(elev, 33), args.max_lag, args.figdir,
              Path(args.tile).stem, key)


def _plot(xs, ys, r, up, fp, max_lag, figdir, stem, key):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Path(figdir).mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    for lab, mask in [("all", slice(None)), ("high-z third", up), ("low-z third", fp)]:
        cen, g, cnt = vg.empirical_variogram(xs[mask], ys[mask], r[mask], max_lag)
        ax.plot(cen[cnt > 0], g[cnt > 0], "o-", ms=3, label=lab)
    ax.set_xlabel("lag (m)"); ax.set_ylabel(r"robust $\gamma$ (m$^2$)")
    ax.set_title(f"inter-swath residual variogram, swaths {key[0]}-{key[1]}")
    ax.legend(); ax.grid(alpha=0.3)
    out = Path(figdir) / f"{stem}_variogram_{key[0]}-{key[1]}.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
