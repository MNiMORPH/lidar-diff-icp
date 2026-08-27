#!/usr/bin/env python3
"""NGV -- the near-ground vegetation index -- defined once, on raw return heights.

    NGV(p) = ( returns with NGV_LO < h <= NGV_HI ) / ( ALL returns within RADIUS of p )

    h = (z - S(x, y)) / sqrt(1 + gx^2 + gy^2)      slope-normal height

with S the order-2 least-squares surface through the CLASS-2 returns within RADIUS of p.
Order 2 removes local slope and curvature, so neither enters the index. Because S is fitted
from the same neighbourhood's own returns, NGV is invariant to any vertical shift of the
cloud: it structurally cannot carry offset information.

WHY THIS MODULE EXISTS. The control-mark table read NGV off STORED 20 mm histograms using a
bin-CENTRE test, so the (0.14, 0.16] bin was counted whole -- that metric is really
"fraction above 0.14 m". Verified 2026-08-27 on mark 1030_2022_MN: identical histograms
(7960 vs 7957 returns, per-bin agreement within LAS quantisation) but numerator 520 exact
against 705 binned. The coefficient fitted against the binned metric therefore cannot be
applied to an exactly-computed one. This module defines the exact form and refits it, so
that one definition runs from the control marks through to the DEM cells.

    ./lidar-icp/bin/python analysis/ngv.py --marks
"""
import argparse, glob, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import control_lowveg_offset as M

NGV_LO, NGV_HI, RADIUS = 0.15, 4.0, 7.5
BOXES = "data/derived/control_boxes"
STRUCT = "data/derived/control_structure"


def design(u, v, order=2):
    """Column order must match the fitted coefficients exactly: [1, u, v, u^2, v^2, uv]."""
    cols = [np.ones_like(u), u, v]
    if order >= 2:
        cols += [u * u, v * v, u * v]
    return np.column_stack(cols)


def ngv(x, y, z, px, py, coef, radius=RADIUS, lo=NGV_LO, hi=NGV_HI):
    """Exact NGV at (px, py) given a pre-fitted order-2 surface `coef`."""
    m = (x - px) ** 2 + (y - py) ** 2 <= radius ** 2
    n = int(m.sum())
    if n == 0:
        return np.nan, 0
    h = (z[m] - design(x[m] - px, y[m] - py) @ coef) / np.sqrt(1 + coef[1] ** 2 + coef[2] ** 2)
    return float(((h > lo) & (h <= hi)).sum()) / n, n


def at_marks(limit=None):
    """Recompute NGV exactly at every control mark, from the kept box clouds."""
    import laspy
    out = []
    boxes = sorted(glob.glob(os.path.join(BOXES, "gen2_2021_control__*.laz")))
    for i, b in enumerate(boxes if limit is None else boxes[:limit]):
        pid = os.path.basename(b).split("__")[1][:-4]
        sf = os.path.join(STRUCT, f"gen2_2021_control__{pid}.npz")
        if not os.path.exists(sf):
            continue
        z_ = np.load(sf)
        if "surface_coef" not in z_.files:
            continue
        f = laspy.read(b)
        X, Y, Z = np.asarray(f.x), np.asarray(f.y), np.asarray(f.z)
        E, N, R = float(z_["easting"]), float(z_["northing"]), float(z_["struct_radius"])
        v, n = ngv(X, Y, Z, E, N, z_["surface_coef"], R)
        # how far outside the stored histograms' support do returns actually reach?
        m = (X - E) ** 2 + (Y - N) ** 2 <= R ** 2
        h = (Z[m] - design(X[m] - E, Y[m] - N) @ z_["surface_coef"]) / \
            np.sqrt(1 + z_["surface_coef"][1] ** 2 + z_["surface_coef"][2] ** 2)
        out.append(dict(point_id=pid, ngv=v, n_disc=n, radius=R,
                        n_below_m1=int((h <= -1.0).sum()), n_above_45=int((h > 45.0).sum())))
    return pd.DataFrame(out)


MIN_CLASS2 = 20      # order-2 needs 6; below this the surface is not meaningfully constrained


def _cell_index(x, y, x0, y0, nx, ny, res):
    """CSR-like index of points into res-sized bins, for O(1) neighbourhood gather."""
    j = np.floor((x - x0) / res).astype(np.int64)
    i = np.floor((y - y0) / res).astype(np.int64)
    ok = (j >= 0) & (j < nx) & (i >= 0) & (i < ny)
    flat = np.where(ok, i * nx + j, -1)
    order = np.argsort(flat, kind="stable")
    flat_s = flat[order]
    start = np.searchsorted(flat_s, np.arange(nx * ny), side="left")
    stop = np.searchsorted(flat_s, np.arange(nx * ny), side="right")
    return order, start, stop


def tile(tile_dir, copc, band_rows=60, out=None, progress=True):
    """Per-cell NGV over a tile's grid, banded so memory stays flat."""
    import json, pdal, time
    cfg = json.load(open(os.path.join(tile_dir, "corrections.json")))
    X0, Y0, X1, Y1 = cfg["bounds"]; res = float(cfg["res_m"])
    nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))
    NG = np.full((ny, nx), np.nan)
    NPTS = np.zeros((ny, nx), np.int32)
    NC2 = np.zeros((ny, nx), np.int32)
    halo = RADIUS + res
    print(f"  grid {ny} x {nx} at {res:g} m, bounds {X0:.1f} {Y0:.1f} {X1:.1f} {Y1:.1f}")
    print(f"  radius {RADIUS} m, window ({NGV_LO}, {NGV_HI}] m, min class-2 {MIN_CLASS2}")
    t0 = time.time()
    for r0 in range(0, ny, band_rows):
        r1 = min(r0 + band_rows, ny)
        yb0, yb1 = Y0 + r0 * res - halo, Y0 + r1 * res + halo
        pl = pdal.Pipeline(json.dumps({"pipeline": [
            {"type": "readers.copc", "filename": copc,
             "bounds": f"([{X0 - halo},{X1 + halo}],[{yb0},{yb1}])"}]}))
        n = pl.execute()
        if n == 0:
            continue
        a = pl.arrays[0]
        px = a["X"].astype(np.float64); py = a["Y"].astype(np.float64)
        pz = a["Z"].astype(np.float64); cl = a["Classification"]
        # local bin grid covering the band plus halo
        bx0 = X0 - halo; by0 = yb0
        bnx = int(np.ceil((X1 + halo - bx0) / res)); bny = int(np.ceil((yb1 - by0) / res))
        order, start, stop = _cell_index(px, py, bx0, by0, bnx, bny, res)
        for i in range(r0, r1):
            cy = Y0 + (i + 0.5) * res
            bi = int((cy - by0) // res)
            for j in range(nx):
                cx = X0 + (j + 0.5) * res
                bj = int((cx - bx0) // res)
                idx = []
                for di in range(bi - 2, bi + 3):
                    if di < 0 or di >= bny:
                        continue
                    base = di * bnx
                    lo_j = max(bj - 2, 0); hi_j = min(bj + 2, bnx - 1)
                    s_ = start[base + lo_j]; e_ = stop[base + hi_j]
                    if e_ > s_:
                        idx.append(order[s_:e_])
                if not idx:
                    continue
                sel = np.concatenate(idx)
                dx = px[sel] - cx; dy = py[sel] - cy
                k = dx * dx + dy * dy <= RADIUS * RADIUS
                if not k.any():
                    continue
                sel = sel[k]; dx = dx[k]; dy = dy[k]
                NPTS[i, j] = sel.size
                g = cl[sel] == 2
                NC2[i, j] = int(g.sum())
                if NC2[i, j] < MIN_CLASS2:
                    continue
                A = design(dx[g], dy[g])
                try:
                    coef = np.linalg.solve(A.T @ A, A.T @ pz[sel][g])
                except np.linalg.LinAlgError:
                    continue
                h = (pz[sel] - design(dx, dy) @ coef) / np.sqrt(1 + coef[1] ** 2 + coef[2] ** 2)
                NG[i, j] = ((h > NGV_LO) & (h <= NGV_HI)).sum() / sel.size
        if progress:
            done = (r1 - 0) / ny
            el = time.time() - t0
            print(f"    rows {r0:4d}-{r1:4d}  {n:9,d} pts read  "
                  f"{el:6.1f} s elapsed, ~{el / done - el:5.1f} s left", flush=True)
    out = out or os.path.join(tile_dir, "ngv.npy")
    np.save(out, NG)
    np.save(out.replace(".npy", "_npts.npy"), NPTS)
    np.save(out.replace(".npy", "_nclass2.npy"), NC2)
    fin = np.isfinite(NG)
    print(f"\n  cells {NG.size:,}   with NGV {int(fin.sum()):,} ({100*fin.mean():.1f}%)")
    print(f"  BLANK cells and why: {int((NPTS == 0).sum()):,} had no returns in the disc; "
          f"{int(((NPTS > 0) & (NC2 < MIN_CLASS2)).sum()):,} had < {MIN_CLASS2} class-2 returns")
    print(f"  NGV percentiles: " + "  ".join(
        f"p{q}={np.percentile(NG[fin], q):.3f}" for q in (0, 10, 25, 50, 75, 90, 99, 100)))
    print(f"  returns per disc: median {np.median(NPTS[fin]):.0f}, "
          f"class-2 median {np.median(NC2[fin]):.0f}")
    print(f"wrote {out}")
    return NG


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--marks", action="store_true",
                    help="recompute NGV exactly at the control marks and refit the offset")
    ap.add_argument("--tile", default=None, help="tile directory, e.g. data/derived/elba_fulldensity")
    ap.add_argument("--copc", default=None, help="gen2 COPC to read")
    ap.add_argument("--band-rows", type=int, default=60)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default="data/derived/control_ngv_exact.csv")
    a = ap.parse_args()
    if not (a.marks or a.tile):
        ap.error("pick a section")
    if a.tile:
        if not a.copc:
            ap.error("--tile needs --copc")
        tile(a.tile, a.copc, a.band_rows)
        return

    print(f"NGV = returns in ({NGV_LO}, {NGV_HI}] m / ALL returns within {RADIUS} m, "
          f"slope-normal above an order-2 class-2 surface")
    t = at_marks(a.limit)
    print(f"  recomputed at {len(t)} marks from the kept boxes")
    print(f"  returns outside the stored histograms' support, summed over marks: "
          f"{int(t.n_below_m1.sum())} below -1 m, {int(t.n_above_45.sum())} above 45 m "
          f"-- so 'all returns in the disc' and 'the histogram support' are the same "
          f"population here, and the exact denominator needs no window")

    m = M.load(0.15, 2.0).merge(t[["point_id", "ngv", "n_disc"]], on="point_id", how="inner")
    print(f"  merged with the offset table: {len(m)} marks\n")

    from scipy import stats as st
    print(f"  NGV exact vs the BINNED metric the old coefficient was fitted against:")
    print(f"    Pearson {st.pearsonr(m.ngv, m.lowveg)[0]:+.4f}   "
          f"exact - binned: median {np.median(m.ngv - m.lowveg):+.4f}, "
          f"max |d| {np.abs(m.ngv - m.lowveg).max():.4f}")

    y = m.resid_mm.to_numpy(float); x = m.ngv.to_numpy(float)
    lr = st.linregress(x, y)
    bo, so = M.fit_origin(x, y)
    print(f"\n  REFIT on the exact index, per-mark:")
    print(f"    free intercept : a {lr.intercept:+7.1f} +/- {lr.intercept_stderr:.1f}   "
          f"b {lr.slope:+8.1f} +/- {lr.stderr:.1f}   (p {lr.pvalue:.1e})")
    print(f"    through origin : b {bo:+8.1f} +/- {so:.1f}")

    blk = M._blocks(m.easting.to_numpy(float), m.northing.to_numpy(float), 10.0)
    ub = np.unique(blk); rng = np.random.default_rng(0); sl = []
    for _ in range(2000):
        idx = np.concatenate([np.where(blk == k)[0] for k in rng.choice(ub, len(ub), True)])
        if len(np.unique(x[idx])) > 2:
            sl.append(np.polyfit(x[idx], y[idx], 1)[0])
    sl = np.array(sl)
    print(f"    block bootstrap, {len(ub)} blocks of 10 km: b = {sl.mean():+.1f} +/- {sl.std(ddof=1):.1f}")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    t.to_csv(a.out, index=False)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
