#!/usr/bin/env python3
"""Do any 2008 gen1 flight lines cross the N-S pattern -- i.e. are there QA tie lines?

The per-line across-track coefficients ``c_s`` measured in
``analysis/SWATH_ACROSS_TRACK_TEST.md`` are recoverable from overlaps only up to an
ALTERNATING vector, because consecutive lines are flown there-and-back so every overlap
constrains a pair SUM ``(c_A+c_B)/2``.  A line flown PERPENDICULAR to the pattern would
cut every N-S line at a different across-track position on each and break that outright.

This run answers, from the tiles already on disk and nothing fetched:

  1. every ``point_source_id`` present, its fitted ground-track HEADING, and whether any
     heading departs from the N-S boustrophedon;
  2. the nadir-track easting and swath half-width of every line, hence the line spacing
     and whether any SECOND-NEIGHBOUR (n to n+2) overlap exists;
  3. the vendor classification histogram, in particular class 9 (water) -- a lake is a
     level surface and would calibrate one line's across-track ramp on its own.

Method for the track fit is the one already used in ``analysis/ELBAEXT2_SCOPE.md`` Sec 2:
returns with |scan_angle_rank| <= 1 sit under the aircraft, and x and y are regressed on
gps_time.  Nothing is fetched; every LAZ is streamed in chunks.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/degeneracy_flightline_inventory.py
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import laspy

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from trust.provenance import Run

REF_N = 4_884_126.0     # northing at which track eastings are reported (ELBAEXT2_SCOPE Sec 2)


def scan_tile(path, nadir_deg, chunk):
    """Streamed pass: per-PSID sufficient statistics, plus the tile's class histogram."""
    acc = {}
    classes = np.zeros(32, np.int64)
    with laspy.open(path) as f:
        hdr = f.header
        for pts in f.chunk_iterator(chunk):
            psid = np.asarray(pts.point_source_id)
            sar = np.asarray(pts.scan_angle_rank).astype(np.float64)
            x = np.asarray(pts.x)
            y = np.asarray(pts.y)
            t = np.asarray(pts.gps_time)
            cl = np.asarray(pts.classification).astype(np.int64)
            classes += np.bincount(np.clip(cl, 0, 31), minlength=32)
            for s in np.unique(psid):
                m = psid == s
                a = acc.setdefault(int(s), dict(n=0, nn=0, tmin=np.inf, tmax=-np.inf,
                                                xmin=np.inf, xmax=-np.inf,
                                                ymin=np.inf, ymax=-np.inf,
                                                amax=0.0, S=np.zeros(7)))
                a["n"] += int(m.sum())
                a["tmin"] = min(a["tmin"], float(t[m].min()))
                a["tmax"] = max(a["tmax"], float(t[m].max()))
                a["xmin"] = min(a["xmin"], float(x[m].min()))
                a["xmax"] = max(a["xmax"], float(x[m].max()))
                a["ymin"] = min(a["ymin"], float(y[m].min()))
                a["ymax"] = max(a["ymax"], float(y[m].max()))
                a["amax"] = max(a["amax"], float(np.abs(sar[m]).max()))
                nd = m & (np.abs(sar) <= nadir_deg)
                if nd.any():
                    tt, xx, yy = t[nd], x[nd], y[nd]
                    a["nn"] += int(nd.sum())
                    a["S"] += np.array([tt.size, tt.sum(), (tt * tt).sum(), xx.sum(),
                                        (tt * xx).sum(), yy.sum(), (tt * yy).sum()])
    return acc, classes, hdr


def fit_track(S):
    """Linear fit of x and y on gps_time from the nadir sufficient statistics."""
    n, st, stt, sx, stx, sy, sty = S
    if n < 50:
        return None
    den = n * stt - st * st
    if den <= 0:
        return None
    bx = (n * stx - st * sx) / den
    by = (n * sty - st * sy) / den
    ax = (sx - bx * st) / n
    ay = (sy - by * st) / n
    speed = float(np.hypot(bx, by))
    heading = float(np.degrees(np.arctan2(bx, by)) % 360.0)   # compass bearing of motion
    if abs(by) < 1e-9:
        x_ref = float("nan")
    else:
        t_ref = (REF_N - ay) / by
        x_ref = ax + bx * t_ref
    return dict(bx=bx, by=by, ax=ax, ay=ay, speed=speed, heading=heading, x_ref=x_ref)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--glob", default="data/before/4*.laz")
    ap.add_argument("--nadir-deg", type=float, default=1.0)
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--out", default="")
    A = ap.parse_args()

    R = Run("do any 2008 gen1 flight lines run perpendicular to the N-S pattern (a QA tie "
            "line), and is there any second-neighbour swath overlap?")
    paths = sorted(glob.glob(A.glob))
    for p in paths:
        R.input(p, role="raw 2008 MN DNR gen1 LAZ tile, vendor classification, "
                        "point_source_id = flight line")
    R.param("nadir_deg", A.nadir_deg, src="repo",
            why="|scan_angle_rank| <= 1 deg selects returns under the aircraft; the track-fit "
                "method of analysis/ELBAEXT2_SCOPE.md Sec 2")
    R.param("ref_northing", REF_N, src="repo",
            why="the northing at which ELBAEXT2_SCOPE.md Sec 2 reports every track easting")
    R.param("chunk", A.chunk, src="MINE",
            why="points per streamed read; a memory ceiling only -- it changes no result, and "
                "the sufficient statistics are exact regardless of how the file is split")

    for c, d in [
        ("tile", "gen1 LAZ tile the line was measured in"),
        ("psid", "point_source_id, i.e. the flight line"),
        ("n", "returns of that line in that tile"),
        ("n_nadir", "returns with |scan_angle_rank| <= 1 deg, used for the track fit"),
        ("heading", "compass bearing of the fitted nadir ground track, degrees"),
        ("speed", "along-track speed of the fitted track, m/s"),
        ("x_ref", "easting where the fitted track crosses N 4 884 126, m (UTM 15N)"),
        ("halfwidth", "99.5th percentile of |cross-track distance| over all that line's "
                      "returns in the tile, m"),
        ("amax", "maximum |scan_angle_rank| of that line in the tile, degrees"),
        ("t0", "minimum gps_time of that line in the tile, GPS seconds of week"),
        ("t1", "maximum gps_time of that line in the tile, GPS seconds of week"),
        ("dx_span", "easting span of that line's returns in the tile, m"),
        ("dy_span", "northing span of that line's returns in the tile, m"),
        ("cls", "vendor LAS classification code"),
        ("count", "returns carrying that classification, summed over all tiles read"),
        ("psid_a", "lower-numbered flight line of the pair"),
        ("psid_b", "higher-numbered flight line of the pair"),
        ("gap", "psid_b - psid_a, i.e. 1 = adjacent, 2 = second neighbour"),
        ("spacing", "|x_ref(b) - x_ref(a)|, the nadir-track separation, m"),
        ("hw_sum", "halfwidth(a) + halfwidth(b), the separation at which the swaths just touch, m"),
        ("sidelap", "hw_sum - spacing: positive = the two swaths overlap by that width, m"),
    ]:
        R.column(c, d)

    R.banner()
    print()

    rows, cls_tot = [], np.zeros(32, np.int64)
    lines = {}
    for p in paths:
        acc, classes, hdr = scan_tile(p, A.nadir_deg, A.chunk)
        cls_tot += classes
        # second pass for half-width: cross-track distance to the fitted track
        fits = {s: fit_track(a["S"]) for s, a in acc.items()}
        hw = {s: [] for s in acc}
        with laspy.open(p) as f:
            for pts in f.chunk_iterator(A.chunk):
                psid = np.asarray(pts.point_source_id)
                x, y = np.asarray(pts.x), np.asarray(pts.y)
                for s, fi in fits.items():
                    if fi is None:
                        continue
                    m = psid == s
                    if not m.any():
                        continue
                    ux, uy = fi["bx"] / fi["speed"], fi["by"] / fi["speed"]
                    # perpendicular offset from the track line through (ax,ay) at t: use the
                    # point-to-line distance in the plane, signed by the left normal (-uy, ux)
                    dx = x[m] - fi["ax"]
                    dy = y[m] - fi["ay"]
                    # remove the along-track component measured from the fit's own origin
                    ct = -uy * dx + ux * dy
                    hw[s].append(np.abs(ct).astype(np.float32))
        base = os.path.basename(p).replace(".laz", "")
        for s in sorted(acc):
            a, fi = acc[s], fits[s]
            h = np.concatenate(hw[s]) if hw[s] else np.array([np.nan], np.float32)
            hwv = float(np.percentile(h, 99.5))
            rows.append([base, s, a["n"], a["nn"],
                         "n/a" if fi is None else f"{fi['heading']:.2f}",
                         "n/a" if fi is None else f"{fi['speed']:.1f}",
                         "n/a" if fi is None else f"{fi['x_ref']:.0f}",
                         f"{hwv:.0f}", f"{a['amax']:.0f}",
                         f"{a['tmin']:.0f}", f"{a['tmax']:.0f}",
                         f"{a['xmax']-a['xmin']:.0f}", f"{a['ymax']-a['ymin']:.0f}"])
            if fi is not None:
                lines.setdefault(s, []).append(dict(tile=base, n=a["n"], hw=hwv, **fi))
        del acc, fits, hw

    print("## 1. every flight line present, per tile")
    R.table(["tile", "psid", "n", "n_nadir", "heading", "speed", "x_ref", "halfwidth",
             "amax", "t0", "t1", "dx_span", "dy_span"], rows)
    print()

    # ---- consolidated per-line geometry, weighted by returns
    print("## 2. consolidated per line (returns-weighted mean of the tiles it appears in)")
    con = []
    for s in sorted(lines):
        e = lines[s]
        w = np.array([d["n"] for d in e], float)
        hd = np.array([d["heading"] for d in e])
        hd = np.where(hd > 180, hd - 360, hd)               # unwrap for averaging
        xr = np.array([d["x_ref"] for d in e])
        hws = np.array([d["hw"] for d in e])
        con.append(dict(psid=s, n=int(w.sum()), heading=float((hd * w).sum() / w.sum()),
                        x_ref=float((xr * w).sum() / w.sum()), hw=float(hws.max()),
                        ntile=len(e), xspread=float(xr.max() - xr.min())))
    R.table(["psid", "n", "heading", "x_ref", "halfwidth"],
            [[c["psid"], c["n"], f"{c['heading'] % 360:.2f}", f"{c['x_ref']:.0f}",
              f"{c['hw']:.0f}"] for c in con])
    print()

    hd_all = np.array([c["heading"] for c in con])
    nsish = np.abs(((hd_all + 90) % 180) - 90) <= 15.0      # within 15 deg of due N or due S
    print(f"lines fitted                          : {len(con)}")
    print(f"headings within 15 deg of the N-S axis: {int(nsish.sum())} of {len(con)}")
    print(f"heading range, unwrapped to [-90,+90] : "
          f"{np.min(((hd_all + 90) % 180) - 90):+.2f} .. {np.max(((hd_all + 90) % 180) - 90):+.2f} deg")
    off = [c for c, k in zip(con, nsish) if not k]
    print(f"CROSS / TIE LINE CANDIDATES (>15 deg off the N-S axis): "
          f"{[c['psid'] for c in off] if off else 'NONE'}")
    print()

    print("## 3. line spacing and second-neighbour overlap")
    ps = [c for c in con if c["n"] >= 100_000]
    ps.sort(key=lambda c: c["x_ref"])
    ov = []
    for i, a in enumerate(ps):
        for j in (i + 1, i + 2):
            if j < len(ps):
                b = ps[j]
                sp = abs(b["x_ref"] - a["x_ref"])
                hs = a["hw"] + b["hw"]
                ov.append([min(a["psid"], b["psid"]), max(a["psid"], b["psid"]),
                           j - i, f"{sp:.0f}", f"{hs:.0f}", f"{hs - sp:+.0f}"])
    R.table(["psid_a", "psid_b", "gap", "spacing", "hw_sum", "sidelap"], ov)
    second = [r for r in ov if r[2] == 2 and float(r[5]) > 0]
    print(f"\nSECOND-NEIGHBOUR OVERLAPS (gap=2 with positive sidelap): "
          f"{[(r[0], r[1], r[5]) for r in second] if second else 'NONE'}")
    n1 = [float(r[3]) for r in ov if r[2] == 1]
    print(f"adjacent nadir-track spacing: min {min(n1):.0f} m, median {np.median(n1):.0f} m, "
          f"max {max(n1):.0f} m over {len(n1)} adjacent pairs")
    print(f"spacing at which a second neighbour would just touch: 2*halfwidth = "
          f"{2*np.median([c['hw'] for c in ps]):.0f} m (median halfwidth "
          f"{np.median([c['hw'] for c in ps]):.0f} m)")
    print()

    print("## 4. vendor classification histogram, all tiles read")
    R.table(["cls", "count"], [[i, int(cls_tot[i])] for i in range(32) if cls_tot[i]])
    print(f"\nclass 9 (water) returns: {int(cls_tot[9]):,} "
          f"({100.0*cls_tot[9]/cls_tot.sum():.4f}% of {int(cls_tot.sum()):,})")

    if A.out:
        json.dump(dict(rows=rows, con=con, ov=ov, cls=cls_tot.tolist()), open(A.out, "w"), indent=1)
        print(f"\nwrote {A.out}")
    R.done(headline="flight-line inventory: headings, spacings, second-neighbour overlap, water")


if __name__ == "__main__":
    main()
