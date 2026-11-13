#!/usr/bin/env python3
"""CONTROL: does gen2 show the same per-swath along-track trend that gen1 does?

gen1's per-swath horizontal offset varies along each pass, and the trend appears to flip
sign with flight direction. But at elbaext flight direction alternates PERFECTLY with
flight order, so heading and elapsed time are the same variable and the sign flip proves
nothing on its own.

gen2 (USGS 3DEP 2021) is the control. It is a modern survey with far better GNSS, flown
over the same ground, and it carries its own point_source_id and gps_time. If the same
direction-flipping along-track trend appears in gen2, the effect is METHOD or TERRAIN. If
it is absent or much smaller, it is gen1's trajectory. This is the same wrong-epoch control
logic that settled the vegetation question, where gen1's return spread predicted the DoD
and gen2's did not.

WHAT IS FITTED. For each gen2 swath, against gen2's own merged reference surface
(z_after), the same three-parameter gradient regression used on gen1:

    d = dx * dz/deast + dy * dz/dnorth + c

BIAS THAT MUST BE STATED: the reference INCLUDES the swath being fitted, so each fit is
pulled toward zero by roughly 1/n_swaths_covering_that_cell. That makes this control
CONSERVATIVE -- it can only understate gen2's trend, so if gen2 shows a trend anyway the
finding is safe, and if it shows none the bias is the first thing to blame. Reported, not
corrected.

Streams the cloud; never loads it whole.

    ./lidar-icp/bin/python analysis/gen2_swath_control.py --tile elbaext_delong \
        --gen2 data/after/elbaext_3dep_fd_class2.laz
"""
import argparse
import json

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", required=True)
ap.add_argument("--gen2", required=True)
ap.add_argument("--chunk", type=int, default=6_000_000)
ap.add_argument("--min-cells", type=int, default=2000)
ap.add_argument("--blocks", type=int, default=4)
A = ap.parse_args()

import laspy

D = f"data/derived/{A.tile}"
c = json.load(open(f"{D}/corrections.json"))
b = c["bounds"]; res = float(c["res_m"]); X0, Y0 = b[0], b[1]
Zref = np.load(f"{D}/z_after.npy"); ny, nx = Zref.shape
stable = np.load(f"{D}/stable.npy").astype(bool)
NC = ny * nx

Zf = Zref.copy(); bad = ~np.isfinite(Zf)
if bad.any():
    Zf = Zf[tuple(distance_transform_edt(bad, return_distances=False, return_indices=True))]
Zs = gaussian_filter(Zf, 1.2)
gnorth, geast = np.gradient(Zs, res)

# Pass 1: accumulate per (swath, cell) elevation sums, and per-swath time/northing sums
acc, cnt, tsum, tcnt = {}, {}, {}, {}
ysum, gsum, gysum, gtt = {}, {}, {}, {}
n_in = 0
with laspy.open(A.gen2) as f:
    total = f.header.point_count
    for pts in f.chunk_iterator(A.chunk):
        cl = np.asarray(pts.classification)
        keep = cl == 2
        if not keep.any():
            continue
        x = np.asarray(pts.x)[keep]; y = np.asarray(pts.y)[keep]
        z = np.asarray(pts.z)[keep]
        ps = np.asarray(pts.point_source_id)[keep]
        gt = np.asarray(pts.gps_time)[keep]
        ix = ((x - X0) / res).astype(np.int64); iy = ((y - Y0) / res).astype(np.int64)
        ing = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        if not ing.any():
            continue
        cid = iy[ing] * nx + ix[ing]
        z = z[ing]; ps = ps[ing]; gt = gt[ing]; y = y[ing]
        n_in += cid.size
        for s in np.unique(ps):
            m = ps == s
            s = int(s)
            if s not in acc:
                acc[s] = np.zeros(NC, np.float64); cnt[s] = np.zeros(NC, np.int32)
                tsum[s] = np.zeros(NC, np.float64)
                ysum[s] = 0.0; gsum[s] = 0.0; gysum[s] = 0.0; gtt[s] = 0
            np.add.at(acc[s], cid[m], z[m])
            np.add.at(cnt[s], cid[m], 1)
            np.add.at(tsum[s], cid[m], gt[m])
            ysum[s] += float(y[m].sum()); gsum[s] += float(gt[m].sum())
            gysum[s] += float((y[m] * gt[m]).sum()); gtt[s] += int(m.sum())
print(f"streamed {total:,} points, {n_in:,} class-2 inside the grid, "
      f"{len(acc)} gen2 swaths", flush=True)

st = stable & np.isfinite(Zref)


def fit(cells_mask, s):
    use = cells_mask & (cnt[s].reshape(ny, nx) >= 3)
    n = int(use.sum())
    if n < A.min_cells:
        return None
    G = (acc[s].reshape(ny, nx)[use] / cnt[s].reshape(ny, nx)[use])
    d = G - Zref[use]
    M = np.column_stack([geast[use], gnorth[use], np.ones(n)])
    cf, *_ = np.linalg.lstsq(M, d, rcond=None)
    r = d - M @ cf
    s2 = float(r @ r) / max(n - 3, 1)
    se = np.sqrt(np.diag(s2 * np.linalg.pinv(M.T @ M)))
    return n, cf, se


print(f"\n{'swath':>7}{'cells':>9}{'dx (mm)':>10}{'SE':>7}{'dy (mm)':>10}{'SE':>7}"
      f"{'dz (mm)':>10}{'dir':>6}{'dur (s)':>9}")
rows = {}
for s in sorted(acc):
    out = fit(st, s)
    if out is None:
        continue
    n, cf, se = out
    tc = cnt[s] > 0
    tmean = tsum[s][tc] / cnt[s][tc]
    dur = float(tmean.max() - tmean.min())
    # direction: sign of cov(gps_time, northing) over the swath's points
    nn = gtt[s]
    cov = gysum[s] / nn - (ysum[s] / nn) * (gsum[s] / nn)
    rows[s] = dict(n=n, dx=cf[0], dy=cf[1], dz=cf[2], se=se, dur=dur,
                   direction=1 if cov > 0 else -1)
    print(f"{s:>7}{n:>9,}{1000*cf[0]:>10.0f}{1000*se[0]:>7.0f}{1000*cf[1]:>10.0f}"
          f"{1000*se[1]:>7.0f}{1000*cf[2]:>10.1f}"
          f"{'N' if cov > 0 else 'S':>6}{dur:>9.0f}")

print(f"\nALONG-TRACK TREND per gen2 swath ({A.blocks} blocks by gps_time within the pass)")
print(f"{'swath':>7}{'dir':>5}" + "".join(f"{f'blk{k+1}':>10}" for k in range(A.blocks))
      + f"{'slope/blk':>11}")
slopes = {}
for s in sorted(rows):
    tc = cnt[s] > 0
    tmean = np.full(NC, np.nan); tmean[tc] = tsum[s][tc] / cnt[s][tc]
    tg = tmean.reshape(ny, nx)
    q = np.nanquantile(tg[st & np.isfinite(tg)], np.linspace(0, 1, A.blocks + 1))
    v = []
    for k in range(A.blocks):
        msk = st & np.isfinite(tg) & (tg >= q[k]) & (tg <= q[k + 1])
        o = fit(msk, s)
        v.append(1000 * o[1][1] if o is not None else np.nan)
    if np.isfinite(v).sum() >= 3:
        ok = np.isfinite(v)
        sl = np.polyfit(np.arange(A.blocks)[ok], np.array(v)[ok], 1)[0]
        slopes[s] = (sl, rows[s]["direction"])
        print(f"{s:>7}{'N' if rows[s]['direction'] > 0 else 'S':>5}"
              + "".join(f"{q_:>10.1f}" if np.isfinite(q_) else f"{'-':>10}" for q_ in v)
              + f"{sl:>11.1f}")
N = [v for v, d_ in slopes.values() if d_ > 0]
S = [v for v, d_ in slopes.values() if d_ < 0]
print(f"\n  mean dy slope  NORTH-flying {np.mean(N) if N else float('nan'):+7.1f} "
      f"(n={len(N)})   SOUTH-flying {np.mean(S) if S else float('nan'):+7.1f} (n={len(S)})")
print("  gen1 for comparison: NORTH -21.1, SOUTH +12.1 mm/block (n=3 each)")
print("\n  If gen2 shows the same flip, the effect is METHOD or TERRAIN, not gen1's GNSS.")
print("  Remember the self-reference bias: this control can only UNDERSTATE gen2's trend.")
