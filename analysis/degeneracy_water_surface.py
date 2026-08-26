#!/usr/bin/env python3
"""Does a level water surface measure a single gen1 flight line's across-track ramp?

The overlap network cannot separate a per-line across-track coefficient ``c_s`` from a
global cross-track tilt of the mosaic (``analysis/degeneracy_identifiability.py`` Sec 1
verifies the two are the same null direction).  Breaking that needs a reference that is
level in the GRAVITY sense, not merely shared between two lines.

``analysis/degeneracy_flightline_inventory.py`` found 562,664 vendor class-9 (water)
returns in the local gen1 tiles, and 458,239 of them are in one tile: ``4358-26-03``, at a
modal elevation near 201 m over N 4,893,254-4,896,747 -- the Mississippi River Pool 5
backwaters (Weaver Bottoms).  A navigation pool above a lock and dam is a level surface.

Three flight lines cross it, each over its FULL scan range.  Within ONE line, regressing
the water-return elevation on ``tan(scan angle)`` measures that line's own ``c_s``, with no
other line, no chain and no gen2 involved.  The three lines fly alternating directions, so
a residual hydraulic slope of the pool itself would show up as a coefficient that ALTERNATES
in sign with the line's body-fixed side, while a per-line roll would not -- printed below as
``c_over_h``, which is the implied cross-track tilt and must be common to all three lines if
the signal is the pool rather than the instrument.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/degeneracy_water_surface.py
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import laspy

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from trust.provenance import Run


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", default="data/before/4358-26-03.laz")
    ap.add_argument("--dz-m", type=float, default=0.30)
    ap.add_argument("--block-m", type=float, default=50.0)
    ap.add_argument("--h-m", type=float, default=2562.0)
    ap.add_argument("--chunk", type=int, default=2_000_000)
    A = ap.parse_args()

    R = Run("does the level Pool 5 water surface measure one gen1 flight line's across-track "
            "coefficient on its own, without any other line?")
    R.input(A.tile, role="raw 2008 gen1 LAZ tile holding the Mississippi Pool 5 backwater, "
                         "vendor classification (9 = water)")
    R.param("dz_m", A.dz_m, src="MINE",
            why="keeps class-9 returns within this distance of the line's modal water "
                "elevation, to exclude a second water body or a mis-classified bank; the "
                "effect of 0.15/0.30/0.60 m is measured and printed below")
    R.param("block_m", A.block_m, src="repo",
            why="the 50 m cluster-robust block size used by swath_across_track_test.py")
    R.param("h_m", A.h_m, src="repo", why="flying height, SWATH_ACROSS_TRACK_TEST.md Sec 0")
    R.param("chunk", A.chunk, src="MINE",
            why="points per streamed read, a memory ceiling only; the reduction is done "
                "after the whole tile is assembled so it changes no result")

    for c, d in [
        ("psid", "point_source_id, i.e. the flight line, dimensionless integer"),
        ("dz_win", "half-window about the modal water elevation, m"),
        ("n", "class-9 water returns of that line inside the window"),
        ("blocks", f"independent {A.block_m:g} m spatial blocks they span"),
        ("x_span", "easting span of those returns, m"),
        ("tan_lo", "minimum tan(scan_angle_rank) over them"),
        ("tan_hi", "maximum of the same"),
        ("z_sd", "standard deviation of the water-return elevation, m"),
        ("c_mm", "OLS slope of elevation on tan(scan angle) within the line, mm per unit tangent "
                 "-- the line's own across-track coefficient if the surface is level"),
        ("c_se", "its cluster-robust standard error, mm per unit tangent"),
        ("h_sign", "sign of the line's body-fixed side, from its heading (N-bound = +1)"),
        ("c_over_h", "c_mm / (h_sign * h): the cross-track TILT the coefficient implies, mm/km. "
                     "Common to all lines = the pool's own slope; differing = per-line roll"),
        ("dzdx", "OLS slope of the water elevation on EASTING, mm/km: the pool surface slope "
                 "as the lidar sees it, pooled over the line"),
        ("dzdy", "the same on NORTHING, mm/km"),
        ("pair", "the two flight lines differenced, dimensionless ids"),
        ("c_A_hat", "the lower-numbered line's raw water-fit coefficient, mm per unit tangent"),
        ("c_B_hat", "the higher-numbered line's raw water-fit coefficient, mm per unit tangent"),
        ("q", "(c_A_hat - c_B_hat)/2: the antisymmetric combination the N-S overlaps cannot "
              "see, mm per unit tangent, with the common water response cancelled"),
        ("q_se", "its standard error from the two independent fits, mm per unit tangent"),
        ("q_per_poolslope", "mm per unit tangent of systematic error in q per 1 mm/km of true "
                            "pool surface slope in easting (= h in km)"),
    ]:
        R.column(c, d)
    R.banner()
    print()

    X, Y, Z, S, P = [], [], [], [], []
    with laspy.open(A.tile) as f:
        for pts in f.chunk_iterator(A.chunk):
            m = np.asarray(pts.classification) == 9
            if m.any():
                X.append(np.asarray(pts.x)[m]); Y.append(np.asarray(pts.y)[m])
                Z.append(np.asarray(pts.z)[m])
                S.append(np.asarray(pts.scan_angle_rank)[m].astype(float))
                P.append(np.asarray(pts.point_source_id)[m])
    x = np.concatenate(X); y = np.concatenate(Y); z = np.concatenate(Z)
    s = np.concatenate(S); p = np.concatenate(P)
    print(f"class-9 returns: {x.size:,}   modal elevation {np.median(z):.2f} m   "
          f"extent E {x.min():.0f}-{x.max():.0f}  N {y.min():.0f}-{y.max():.0f}")
    # heading sign per line, from the nadir returns of the whole tile
    hs = {}
    with laspy.open(A.tile) as f:
        for pts in f.chunk_iterator(A.chunk):
            sr = np.asarray(pts.scan_angle_rank).astype(float)
            nd = np.abs(sr) <= 1
            if not nd.any():
                continue
            pp = np.asarray(pts.point_source_id)[nd]
            yy = np.asarray(pts.y)[nd]; tt = np.asarray(pts.gps_time)[nd]
            for u in np.unique(pp):
                m = pp == u
                a = hs.setdefault(int(u), np.zeros(5))
                a += [m.sum(), tt[m].sum(), (tt[m] ** 2).sum(), yy[m].sum(), (tt[m] * yy[m]).sum()]
    sign = {}
    for u, a in hs.items():
        n, st, stt, sy, sty = a
        sign[u] = 1.0 if (n * sty - st * sy) > 0 else -1.0

    def cluster_se(t, zz, blk):
        Xm = np.c_[np.ones_like(t), t]
        XtXi = np.linalg.pinv(Xm.T @ Xm)
        beta = XtXi @ (Xm.T @ zz)
        e = zz - Xm @ beta
        o = np.argsort(blk, kind="stable")
        g = blk[o]; Xs = Xm[o]; es = e[o]
        st = np.flatnonzero(np.r_[True, g[1:] != g[:-1]]); en = np.r_[st[1:], len(g)]
        meat = np.zeros((2, 2))
        for i, j in zip(st, en):
            sc = Xs[i:j].T @ es[i:j]
            meat += np.outer(sc, sc)
        G = len(st)
        corr = (G / max(G - 1, 1)) * ((len(zz) - 1) / max(len(zz) - 2, 1))
        V = corr * (XtXi @ meat @ XtXi)
        return beta, np.sqrt(np.diag(V)), G

    rows = []
    for dz in (0.15, A.dz_m, 0.60):
        for u in sorted(np.unique(p)):
            m = p == u
            if m.sum() < 5000:
                continue
            zm = np.median(z[m])
            w = m & (np.abs(z - zm) < dz)
            t = np.tan(np.radians(s[w]))
            blk = (np.floor(y[w] / A.block_m).astype(np.int64) << 20) + \
                  np.floor(x[w] / A.block_m).astype(np.int64)
            beta, se, G = cluster_se(t, z[w] * 1000.0, blk)
            sg = sign.get(int(u), np.nan)
            bx = np.polyfit(x[w], z[w] * 1e6, 1)[0]      # mm per km
            by = np.polyfit(y[w], z[w] * 1e6, 1)[0]
            rows.append([int(u), f"{dz:.2f}", int(w.sum()), G, f"{np.ptp(x[w]):.0f}",
                         f"{t.min():+.3f}", f"{t.max():+.3f}", f"{z[w].std():.3f}",
                         f"{beta[1]:+.1f}", f"{se[1]:.1f}", f"{sg:+.0f}",
                         f"{beta[1]/(sg*A.h_m)*1000:+.2f}", f"{bx:+.1f}", f"{by:+.1f}"])
    R.table(["psid", "dz_win", "n", "blocks", "x_span", "tan_lo", "tan_hi", "z_sd",
             "c_mm", "c_se", "h_sign", "c_over_h", "dzdx", "dzdy"], rows)
    print()
    print("## What survives the water's own scan-angle response")
    print("A within-line water fit returns c_s + w + g_pool*h_s, where w is the water surface's")
    print("OWN response to view angle (specular/wave-facet geometry, identical in the body frame")
    print("on every line) and g_pool is the pool's hydraulic slope in easting. Two lines of")
    print("OPPOSITE heading therefore give: c_A - c_B = (c_A_hat - c_B_hat) - 2*g_pool*h.")
    prow = []
    for dz in (0.15, A.dz_m, 0.60):
        r = {int(x[0]): x for x in rows if float(x[1]) == dz}
        if 144 in r and 145 in r:
            ca, sa = float(r[144][8]), float(r[144][9])
            cb, sb = float(r[145][8]), float(r[145][9])
            q = 0.5 * (ca - cb)
            sq = 0.5 * np.hypot(sa, sb)
            prow.append([f"{dz:.2f}", "144-145", f"{ca:+.1f}", f"{cb:+.1f}", f"{q:+.1f}",
                         f"{sq:.1f}", f"{A.h_m/1000.0:.2f}"])
    R.table(["dz_win", "pair", "c_A_hat", "c_B_hat", "q", "q_se", "q_per_poolslope"], prow)
    print("q_per_poolslope: each 1 mm/km of true pool surface slope in EASTING adds this many")
    print("mm per unit tangent to q, so a pool level to 10 mm/km leaves a 26 mm/unit-tangent")
    print("systematic beside the quoted standard error. Compare the N-S overlaps' q_se of")
    print("30-93 mm per unit tangent (degeneracy_identifiability.py Sec 2).")
    R.done(headline="per-line across-track coefficient from a level Pool 5 water surface")


if __name__ == "__main__":
    main()
