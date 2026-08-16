#!/usr/bin/env python3
"""Convex-prior warp correction.

SUPERSEDED / DO NOT USE FOR PRODUCTS. This IDW-on-hillslopes warp was found to be
OVERFIT (it hugs the hillslope control and fits ~100 m terrain-correlated noise;
low out-of-sample skill). The real, physical form of the residual is the smooth
per-swath along-track GNSS drift -- see coreg.fit_along_track_drift and
scripts/along_track_drift.py. Kept only as a record of the diagnostic.

The co-registration leaves an undulating, position-dependent vertical warp of
~0.04 m at ~100 m scale (the correlated part of the residual). It is invisible on
flats but the geomorphic prior makes it legible on convex hillslopes, which must
erode or hold (never gain) yet show coherent depositional patches. DeLong et al.
(2022) hit the same warp and handled it with local ICP; here we use the prior as
a constraint instead.

Assumption: over 13 years convex hillslopes change at sub-LoD magnitude, so their
observed M3C2 value approximates the warp. We fit a smooth warp field by inverse-
distance interpolation over CONTROL = flat stable uplands + convex hillslopes
(excluding concave slopes, the floodplain, and |change| > thresh so real large
features are not fit), then subtract it from the M3C2 change. Cells with no
control within the radius are left uncorrected (warp = 0), preserving real
features (banks, channel) far from control.

Caveat: real convex erosion up to `change_thresh` is absorbed as warp; on
diffusive convex slopes that is well below what 13 years can produce, but it is
the price of the constraint.
"""
import argparse
from pathlib import Path
import numpy as np
import laspy
import pandas as pd
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter, uniform_filter, distance_transform_edt as edt

from lidar_diff_icp import coreg


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("change_laz", help="M3C2 change product (median) with m3c2, lod")
    ap.add_argument("after_last_laz", help="2021 last-return cloud (for terrain)")
    ap.add_argument("--bounds", nargs=4, type=float, required=True)
    ap.add_argument("--radius", type=float, default=150.0)
    ap.add_argument("--curv-scale", type=float, default=50.0, help="curvature length (m)")
    ap.add_argument("--change-thresh", type=float, default=0.15)
    ap.add_argument("--out", default="data/derived/change_core2008_convexcorr.laz")
    a = ap.parse_args()
    X0, Y0, X1, Y1 = a.bounds
    res = 5.0; nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))

    def grid(x, y, v, q=None):
        ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & np.isfinite(v)
        gb = pd.Series(v[ok]).groupby(iy[ok] * nx + ix[ok])
        s = gb.quantile(q) if q is not None else gb.median()
        out = np.full(nx * ny, np.nan); out[s.index.values] = s.values
        return out.reshape(ny, nx)

    # terrain from 2021 ground (low-10%)
    g = laspy.read(a.after_last_laz)
    gx = np.asarray(g.x); gy = np.asarray(g.y); gz = np.asarray(g.z)
    m = (gx >= X0) & (gx < X1) & (gy >= Y0) & (gy < Y1)
    Z = grid(gx[m], gy[m], gz[m], q=0.10)
    Zf = Z.copy(); nanm = np.isnan(Zf)
    if nanm.any():
        Zf = Zf[tuple(edt(nanm, return_distances=False, return_indices=True))]
    tpi = Z - uniform_filter(Zf, size=int(2 * 300 / res), mode="nearest")
    sdeg = np.degrees(coreg.slope_aspect(gaussian_filter(Zf, 2.0), res)[0])
    Zs = gaussian_filter(Zf, sigma=a.curv_scale / res / 2.0)
    lap = (np.gradient(np.gradient(Zs, res, axis=0), res, axis=0)
           + np.gradient(np.gradient(Zs, res, axis=1), res, axis=1))
    convex = (sdeg > 5) & (sdeg < 35) & (tpi > -2) & (lap < 0)
    flat_up = (sdeg < 3) & (tpi > -2)

    # sample terrain classes at each core point
    h = laspy.read(a.change_laz)
    hx = np.asarray(h.x); hy = np.asarray(h.y); d = np.asarray(h.m3c2)
    ci = np.clip(((hx - X0) / res).astype(int), 0, nx - 1)
    ri = np.clip(((hy - Y0) / res).astype(int), 0, ny - 1)
    is_ctrl = ((convex | flat_up)[ri, ci]) & np.isfinite(d) & (np.abs(d) < a.change_thresh)
    is_convex = convex[ri, ci] & np.isfinite(d)
    print(f"core points {len(d):,}; control {is_ctrl.sum():,} "
          f"(convex {(is_convex & is_ctrl).sum():,} + flat {(flat_up[ri, ci] & is_ctrl).sum():,})", flush=True)

    # IDW warp from control, evaluated at all core points
    tree = cKDTree(np.c_[hx[is_ctrl], hy[is_ctrl]])
    dist, idx = tree.query(np.c_[hx, hy], k=32)
    w = 1.0 / np.maximum(dist, 1e-6) ** 2; w[dist > a.radius] = 0.0
    wsum = w.sum(1)
    warp = np.where(wsum > 0, (w * d[is_ctrl][idx]).sum(1) / np.where(wsum == 0, 1.0, wsum), np.nan)
    corr = d - np.where(np.isfinite(warp), warp, 0.0)

    # --- validation ---
    def nmad(v): return 1.4826 * np.median(np.abs(v - np.median(v)))
    # block CV skill on control
    Xs, Ys, Vs = hx[is_ctrl], hy[is_ctrl], d[is_ctrl]
    blk = (np.floor((Xs - X0) / 150) * 1e5 + np.floor((Ys - Y0) / 150)).astype(int)
    oos = np.full(Vs.size, np.nan)
    for b in np.unique(blk):
        te = blk == b; tr = ~te
        if tr.sum() < 32: continue
        t = cKDTree(np.c_[Xs[tr], Ys[tr]]); dd, ii = t.query(np.c_[Xs[te], Ys[te]], k=32)
        ww = 1 / np.maximum(dd, 1e-6) ** 2; ww[dd > a.radius] = 0; az = ww.sum(1) == 0; ww[az, 0] = 1
        oos[te] = (ww * Vs[tr][ii]).sum(1) / ww.sum(1)
    print(f"warp CV skill on control: NMAD {nmad(Vs):.3f} -> {nmad(Vs - oos):.3f}  "
          f"var explained {1 - np.var(Vs - oos) / np.var(Vs):.2f}", flush=True)

    def coh_dep(val):
        bx = (((hx - X0) // 100)).astype(int); by = (((hy - Y0) // 100)).astype(int)
        mk = is_convex & np.isfinite(val)
        df = pd.DataFrame({"b": by[mk] * 1000 + bx[mk], "d": val[mk]})
        bm = df.groupby("b")["d"].median(); bn = df.groupby("b")["d"].size(); big = bn >= 20
        return 100 * np.average(bm[big] > 0.05, weights=bn[big])
    print(f"convex area coherently depositional: {coh_dep(d):.0f}% -> {coh_dep(corr):.0f}%", flush=True)
    big = np.isfinite(d) & (np.abs(d) > 0.3)
    print(f"large real features (|change|>0.3, n={big.sum():,}): median |change| "
          f"{np.median(np.abs(d[big])):.2f} -> {np.median(np.abs(corr[big])):.2f} (preserved if ~unchanged)", flush=True)
    print(f"warp: median {np.nanmedian(warp):+.3f}, NMAD {nmad(warp[np.isfinite(warp)]):.3f}; "
          f"applied to {100 * np.mean(np.isfinite(warp)):.0f}% of points", flush=True)

    # write
    hdr = laspy.LasHeader(point_format=3, version="1.2")
    core = np.c_[hx, hy, np.asarray(h.z)]
    hdr.offsets = core.min(0); hdr.scales = [0.01, 0.01, 0.01]
    for nm in ("m3c2", "m3c2_raw", "warp", "lod"):
        hdr.add_extra_dim(laspy.ExtraBytesParams(name=nm, type=np.float32))
    las = laspy.LasData(hdr)
    las.x, las.y, las.z = hx, hy, np.asarray(h.z)
    las.m3c2 = corr.astype(np.float32); las.m3c2_raw = d.astype(np.float32)
    las.warp = np.where(np.isfinite(warp), warp, 0.0).astype(np.float32)
    las.lod = np.asarray(h.lod).astype(np.float32)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True); las.write(a.out)
    print(f"wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()
