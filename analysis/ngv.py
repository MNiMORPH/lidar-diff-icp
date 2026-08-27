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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--marks", action="store_true",
                    help="recompute NGV exactly at the control marks and refit the offset")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default="data/derived/control_ngv_exact.csv")
    a = ap.parse_args()
    if not a.marks:
        ap.error("pick a section")

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
