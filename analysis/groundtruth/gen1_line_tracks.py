"""Flight-line ground tracks for the 2008 gen1 acquisition, from the tiles on disk.

Why this exists
---------------
Two different things are needed from the flight lines and they must not be confused:

  * a **targeting** model -- where does line N run, so which tiles hold marks that line
    could have seen?  That is what this file builds: one straight ground track per
    ``point_source_id`` pass, fitted to the nadir track and extrapolated.
  * the **assignment** -- which line actually illuminated a given mark?  That is NOT read
    off this model.  It is read off the ``point_source_id`` of the ground returns at the
    mark itself (``gen1_more_marks_tie.py``).  Distance to a fitted centreline is a proxy
    and mislabels once the search radius approaches the ~1 km line spacing.

Method
------
The nadir track is taken from **near-nadir returns only** (``|scan_angle_rank| <=
NADIR_DEG``).  Using the mean of all returns in a gps_time bin, as a first cut here did,
is biased by up to a swath half-width wherever a tile boundary clips the swath
asymmetrically -- that produced 100-450 m "track residuals" that were tile edges, not
flight.  Near-nadir returns are either in the tile or not; they are never off-centre.

``point_source_id`` is **reused across missions** in this acquisition (psid 151 spans
337 057 s of gps_time, ~3.9 days).  Each psid is therefore split into passes at gps_time
gaps of more than ``GAP_S``, and every pass is fitted and reported separately.  A gap is
not by itself evidence of a second pass -- it is equally a stretch of line under tiles we
do not hold -- so the printed table gives the gap length in seconds and in flight-path
kilometres beside it, and the caller judges.
"""
from __future__ import annotations

import glob
import json
import math
import os

import numpy as np

SCRATCH = os.environ.get("SCRATCH", ".")
STRIDE = 7            # subsample stride; MINE -- only thins the near-nadir strip
NADIR_DEG = 2.0       # |scan_angle_rank| kept for the nadir track; MINE
BIN_S = 1.0           # gps_time bin width, seconds
GAP_S = 120.0         # split a psid into passes at gps_time gaps longer than this; MINE


def main():
    tiles = [t for t in sorted(glob.glob("data/before/*.laz")) if "merged" not in t]
    print(f"tiles on disk: {len(tiles)}")
    print(f"params: STRIDE={STRIDE} (MINE), NADIR_DEG={NADIR_DEG} (MINE), "
          f"BIN_S={BIN_S}, GAP_S={GAP_S} (MINE)")
    import laspy
    per_line: dict[int, list] = {}
    tile_lines: dict[str, dict] = {}
    for i, p in enumerate(tiles):
        f = laspy.read(p)
        sa = np.abs(np.asarray(f.scan_angle_rank).astype(float))
        keep = sa <= NADIR_DEG
        x = np.asarray(f.x)[keep][::STRIDE]
        y = np.asarray(f.y)[keep][::STRIDE]
        t = np.asarray(f.gps_time)[keep][::STRIDE]
        s = np.asarray(f.point_source_id)[keep][::STRIDE].astype(int)
        sall = np.asarray(f.point_source_id).astype(int)
        tall = np.asarray(f.gps_time)
        name = os.path.basename(p)[:-4]
        tl = {}
        for ln in np.unique(sall):
            ma = sall == ln
            tl[int(ln)] = dict(n=int(ma.sum()), t0=float(tall[ma].min()), t1=float(tall[ma].max()))
        for ln in np.unique(s):
            m = s == ln
            tt, xx, yy = t[m], x[m], y[m]
            k = np.round(tt / BIN_S).astype(np.int64)
            order = np.argsort(k)
            k, xx, yy, tt = k[order], xx[order], yy[order], tt[order]
            uk, idx, cnt = np.unique(k, return_index=True, return_counts=True)
            cx = np.add.reduceat(xx, idx) / cnt
            cy = np.add.reduceat(yy, idx) / cnt
            ct = np.add.reduceat(tt, idx) / cnt
            keep2 = cnt >= 20
            if keep2.sum():
                per_line.setdefault(int(ln), []).append(
                    np.stack([ct[keep2], cx[keep2], cy[keep2]], 1))
        tile_lines[name] = tl
        print(f"  [{i+1:2d}/{len(tiles)}] {name}  psids {sorted(tl)}", flush=True)
        del f, x, y, t, s, sall, tall

    print("\n=== fitted ground tracks, one row per pass ===")
    print(f"{'psid':>5} {'pass':>4} {'bins':>5} {'gps_t0':>10} {'gps_t1':>10} {'dur_s':>7} "
          f"{'speed':>6} {'head':>6} {'len_km':>7} {'res_med_m':>9} {'res_p95_m':>9} "
          f"{'gap_before_s':>12} {'gap_km':>7}")
    model = {}
    for ln in sorted(per_line):
        a = np.concatenate(per_line[ln], 0)
        a = a[np.argsort(a[:, 0])]
        uk, idx, cnt = np.unique(np.round(a[:, 0] / BIN_S).astype(np.int64),
                                 return_index=True, return_counts=True)
        ct = np.add.reduceat(a[:, 0], idx) / cnt
        cx = np.add.reduceat(a[:, 1], idx) / cnt
        cy = np.add.reduceat(a[:, 2], idx) / cnt
        brk = np.flatnonzero(np.diff(ct) > GAP_S) + 1
        segs = np.split(np.arange(ct.size), brk)
        prev_end = None
        for si, sidx in enumerate(segs):
            if sidx.size < 4:
                continue
            T, X, Y = ct[sidx], cx[sidx], cy[sidx]
            px = np.polyfit(T, X, 1)
            py = np.polyfit(T, Y, 1)
            r = np.hypot(X - np.polyval(px, T), Y - np.polyval(py, T))
            sp = math.hypot(px[0], py[0])
            head = math.degrees(math.atan2(px[0], py[0])) % 360.0
            gap = (T.min() - prev_end) if prev_end is not None else float("nan")
            prev_end = T.max()
            key = f"{ln}.{si}"
            model[key] = dict(psid=int(ln), seg=si, px=px.tolist(), py=py.tolist(),
                              t0=float(T.min()), t1=float(T.max()), speed=sp, heading=head,
                              resid_med=float(np.median(r)), resid_p95=float(np.percentile(r, 95)),
                              n_bins=int(T.size),
                              p0=[float(np.polyval(px, T.min())), float(np.polyval(py, T.min()))],
                              p1=[float(np.polyval(px, T.max())), float(np.polyval(py, T.max()))])
            d = model[key]
            print(f"{ln:>5} {si:>4} {d['n_bins']:>5} {d['t0']:>10.1f} {d['t1']:>10.1f} "
                  f"{d['t1']-d['t0']:>7.1f} {sp:>6.1f} {head:>6.1f} "
                  f"{(d['t1']-d['t0'])*sp/1000:>7.1f} {d['resid_med']:>9.1f} {d['resid_p95']:>9.1f} "
                  f"{gap:>12.1f} {gap*sp/1000 if gap == gap else float('nan'):>7.1f}")

    # A single straight fit per psid, over every near-nadir bin it has, in every tile.
    # This is the long-baseline targeting track: the per-pass fits above are 3-38 km long
    # and extrapolate to 0.03-1.6 km of miss at the far end of the acquisition, while a
    # fit that spans the whole line does not extrapolate at all where marks actually sit.
    print("\n=== one straight fit per psid, over every near-nadir bin ===")
    print(f"{'psid':>5} {'bins':>5} {'span_km':>8} {'head':>6} {'res_med_m':>9} {'res_p95_m':>9} {'res_max_m':>9}")
    whole = {}
    for ln in sorted(per_line):
        a = np.concatenate(per_line[ln], 0)
        a = a[np.argsort(a[:, 0])]
        uk, idx, cnt = np.unique(np.round(a[:, 0] / BIN_S).astype(np.int64),
                                 return_index=True, return_counts=True)
        ct = np.add.reduceat(a[:, 0], idx) / cnt
        cx = np.add.reduceat(a[:, 1], idx) / cnt
        cy = np.add.reduceat(a[:, 2], idx) / cnt
        if ct.size < 4:
            continue
        # fit x,y to arclength along the dominant direction rather than to time, so a
        # gps_time gap under tiles we do not hold does not lever the fit
        u, v = cx - cx.mean(), cy - cy.mean()
        w = np.polyfit(v, u, 1) if np.ptp(v) >= np.ptp(u) else None
        if w is not None:
            r = np.abs(u - np.polyval(w, v)) / math.hypot(1.0, w[0])
            axis = "N"
        else:
            w = np.polyfit(u, v, 1)
            r = np.abs(v - np.polyval(w, u)) / math.hypot(1.0, w[0])
            axis = "E"
        span = math.hypot(np.ptp(cx), np.ptp(cy)) / 1000.0
        head = math.degrees(math.atan2(cx[-1] - cx[0], cy[-1] - cy[0])) % 360.0
        whole[int(ln)] = dict(axis=axis, w=list(map(float, w)),
                              xm=float(cx.mean()), ym=float(cy.mean()),
                              n_bins=int(ct.size), span_km=span, heading=head,
                              resid_med=float(np.median(r)), resid_p95=float(np.percentile(r, 95)),
                              resid_max=float(r.max()),
                              tmin=float(ct.min()), tmax=float(ct.max()),
                              nmin=float(cy.min()), nmax=float(cy.max()),
                              emin=float(cx.min()), emax=float(cx.max()))
        d = whole[int(ln)]
        print(f"{ln:>5} {d['n_bins']:>5} {span:>8.1f} {head:>6.1f} {d['resid_med']:>9.1f} "
              f"{d['resid_p95']:>9.1f} {d['resid_max']:>9.1f}")
        np.savez(f"{SCRATCH}/linebins_{ln}.npz", t=ct, x=cx, y=cy)

    with open(f"{SCRATCH}/line_tracks.json", "w") as fh:
        json.dump(dict(model=model, whole=whole, tile_lines=tile_lines, stride=STRIDE, bin_s=BIN_S,
                       nadir_deg=NADIR_DEG, gap_s=GAP_S), fh)
    print(f"\nwrote {SCRATCH}/line_tracks.json  ({len(model)} passes)")


if __name__ == "__main__":
    main()
