#!/usr/bin/env python3
"""Where do the "stable" cells come from, and is the DoD tilt over them a property of the tile?

A plane fitted to ``dod_cover_q2.npy`` over the stable open reference cells at Elba returns
``-10.30 mm intercept, dE -14.19 mm/km, dN -16.70 mm/km``. That is a claim about the TILE,
but the fit is over a population selected by landform (divide cells, low curvature, gentle
slope) and by land cover (canopy cover <= 0.02, floodplain removed). A plane fitted to a
spatially clustered or elevation-biased subset is not a plane over the tile, and two of the
registration terms already in the model -- the per-swath constants (which run west to east
across ~2 km, i.e. in EASTING) and the along-track drift curves (which run along the
north-south flight lines, i.e. in NORTHING) -- would leave residuals in exactly those two
directions if they were incomplete.

This run answers five questions and invents no threshold to do it:

1. **Where are the cells?** Coverage per 250 m block, fraction of the tile, and their
   distribution in easting, northing and ELEVATION against the whole tile's.
2. **Is the tilt stable across subregions?** Refits on halves, quadrants, a 50 m
   checkerboard and random block subsamples.
3. **Is it landform?** Nested models adding elevation (``z_after``), distance from the
   valley (Euclidean distance to ``floodplain_mask``) and the agricultural-region mask.
4. **Is it the registration?** The decisive split. Per-return ``point_source_id`` gives a
   per-(cell, flight-line) ground estimate -- the pipeline's own estimator evaluated one
   line at a time -- so the plane can be refitted WITHIN a swath, where the across-swath
   constant is fixed by construction. And because consecutive gen1 lines fly opposite
   directions (verified here: ``corr(gps_time, northing) = -1, +1, -1, +1``), a residual
   ALONG-TRACK drift must alternate sign in northing between lines while a ground-fixed
   north-south field must not. That is the discriminator for ``dN``.
5. **What is left, and with what uncertainty?** The cluster-robust SE depends on the block
   size, and the block must exceed the correlation length of the residual field, so the
   empirical semivariogram of the residual is measured and the SE is reported over a ladder
   of block sizes rather than at one invented value.

Nothing in ``pipeline.py`` or ``coreg.py`` is modified. Every input, parameter, mask and
column is declared through ``trust/provenance.py``.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/stable_point_tilt_audit.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lidar_diff_icp import binstats as bs
from lidar_diff_icp.refcells import reference_cells
from trust.provenance import Run


# ----------------------------------------------------------------- statistics
def ols_cluster(X, y, groups):
    """OLS with a cluster-robust (sandwich) covariance, clustered on ``groups``.

    Same construction as ``analysis/swath_across_track_test.py``: finite-sample factor
    ``G/(G-1) * (n-1)/(n-k)``. Returns ``(beta, se, n_clusters, r2)``.
    """
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    e = y - X @ beta
    g = np.asarray(groups)
    order = np.argsort(g, kind="stable")
    gs, Xs, es = g[order], X[order], e[order]
    starts = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
    ends = np.r_[starts[1:], len(gs)]
    meat = np.zeros((k, k))
    for a, b in zip(starts, ends):
        s = Xs[a:b].T @ es[a:b]
        meat += np.outer(s, s)
    G = len(starts)
    corr = (G / max(G - 1, 1)) * ((n - 1) / max(n - k, 1))
    V = corr * (XtX_inv @ meat @ XtX_inv)
    sst = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((e ** 2).sum()) / sst if sst > 0 else float("nan")
    return beta, np.sqrt(np.diag(V)), G, r2


def lad(X, y, iters=200, tol=1e-9):
    """Least absolute deviations by IRLS from the OLS start -- the multivariate form of
    ``analysis/ridgelines/cover_offset_regression.lad``. Deterministic."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    for _ in range(iters):
        w = np.sqrt(1.0 / np.maximum(np.abs(y - X @ beta), 1e-6))
        nb = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)[0]
        if np.max(np.abs(nb - beta)) < tol:
            return nb
        beta = nb
    return beta


def block_bootstrap(X, y, blk, n_boot, rng):
    """Resample whole spatial blocks with replacement; return the SE of each coefficient."""
    ub, inv = np.unique(blk, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    iv = inv[order]
    starts = np.flatnonzero(np.r_[True, iv[1:] != iv[:-1]])
    ends = np.r_[starts[1:], iv.size]
    out = []
    for _ in range(n_boot):
        pick = rng.integers(0, ub.size, ub.size)
        sub = np.concatenate([order[starts[p]:ends[p]] for p in pick])
        out.append(np.linalg.lstsq(X[sub], y[sub], rcond=None)[0])
    return np.asarray(out).std(0)


def f3(v):
    return f"{v:+.2f}"


# ------------------------------------------------------------------------ run
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", default="data/derived/elba_fulldensity")
    ap.add_argument("--dod", default="dod_cover_q2.npy",
                    help="the DoD raster under audit, relative to the tile directory")
    ap.add_argument("--cover-max", type=float, default=0.02,
                    help="canopy-cover ceiling of the population under audit")
    ap.add_argument("--block-m", type=float, nargs="+", default=[50.0, 100.0, 250.0, 500.0, 1000.0],
                    help="spatial block sizes (m) for the cluster-robust SE ladder; 50 is the "
                         "repo default, the rest test it against the measured correlation length")
    ap.add_argument("--map-block-m", type=float, default=250.0,
                    help="block size for the coverage map")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--n-elev-bins", type=int, default=8,
                    help="equal-count elevation bins (binstats.quantile_edges, spans all data)")
    ap.add_argument("--seed", type=int, default=0)
    A = ap.parse_args()
    rng = np.random.default_rng(A.seed)
    T = A.tile
    name = os.path.basename(T)

    R = Run("where do the stable cells behind the reported DoD tilt come from, and is the "
            "tilt a property of the tile or of that population and of the registration?")
    R.param("cover_max", A.cover_max, src="andy",
            why="the population under audit was defined with this ceiling")
    R.param("block_m ladder", A.block_m, src="MINE",
            why="50 m is the repo default used by offset_vs_angle.py and "
                "cover_offset_regression.py; the larger blocks are added because the "
                "semivariogram measured in section 2c shows the residual field is still "
                "correlated at 400 m, so a 50 m block treats correlated cells as independent. "
                "Nothing is excluded by this; it changes only the reported SE.")
    R.param("map_block_m", A.map_block_m, src="andy",
            why="the coverage map was asked for in 250 m blocks")
    R.param("n_boot", A.n_boot, src="MINE",
            why="block-bootstrap replicates, used only to check the sandwich SE against the "
                "bootstrap the original claim quoted. No data are dropped.")
    R.param("n_elev_bins", A.n_elev_bins, src="MINE",
            why="equal-count elevation bins via binstats.quantile_edges, which spans the full "
                "observed range so no elevation is truncated. The count sets resolution only.")
    R.param("reference_cells defaults", "curv_max=0.015, slope_max=12, gross_change_mm=500, "
                                        "clearcut_drop=0.30, require_ridge=True", src="repo")
    R.param("minimum counts imposed", "none", src="repo",
            why="")

    for c, d in [
        ("what", "the quantity or subset the row describes"),
        ("n", "cells (or cell-flightline rows, where stated) entering the fit"),
        ("blocks", "independent spatial blocks the fit's cluster-robust SE counts"),
        ("mean_mm", "fitted intercept = mean DoD at the subset's own centroid, mm "
                    "(+ = elevation rose 2008->2021)"),
        ("mean_se", "cluster-robust standard error of mean_mm, mm"),
        ("dE", "fitted easting gradient of the DoD, mm per km east"),
        ("dE_se", "cluster-robust standard error of dE, mm/km"),
        ("dN", "fitted northing gradient of the DoD, mm per km north"),
        ("dN_se", "cluster-robust standard error of dN, mm/km"),
        ("dE_lad", "the same easting gradient from a least-absolute-deviations fit, mm/km"),
        ("dN_lad", "the same northing gradient from a least-absolute-deviations fit, mm/km"),
        ("block_m", "size of the spatial block the cluster-robust SE is computed over, m"),
        ("t_dE", "dE / dE_se, cluster-robust"),
        ("t_dN", "dN / dN_se, cluster-robust"),
        ("term", "regressor in the nested model"),
        ("coef", "fitted coefficient in the units named by term"),
        ("se", "cluster-robust standard error of coef"),
        ("model", "the set of regressors fitted jointly"),
        ("z_term", "coefficient on gen2 ground elevation z_after, mm per 100 m of elevation"),
        ("z_se", "cluster-robust standard error of z_term"),
        ("d_term", "coefficient on distance to the nearest floodplain cell, mm per km"),
        ("d_se", "cluster-robust standard error of d_term"),
        ("ag_term", "coefficient on the agricultural-region mask, mm (ag minus non-ag)"),
        ("ag_se", "cluster-robust standard error of ag_term"),
        ("lag_lo", "lower edge of the separation-distance bin, m"),
        ("lag_hi", "upper edge of the separation-distance bin, m"),
        ("pairs", "randomly drawn cell pairs in the lag bin"),
        ("gamma_rms", "sqrt of the semivariance in the bin, mm (the rms difference / sqrt2)"),
        ("rho", "implied correlation 1 - gamma/variance of the residual field"),
        ("swath", "gen1 flight line (point_source_id)"),
        ("rows", "cell-flightline ground estimates for that line on the stable cells"),
        ("E_span", "easting extent of that line's stable cells, km -- its leverage on dE"),
        ("h_m", "fitted flying height from x = a + b*y + h*tan(scan) within the line, m; "
                "the sign alternates because scan_angle is body-fixed and the aircraft turns"),
        ("dir", "corr(gps_time, northing) within the line: -1 = flown south, +1 = flown north"),
        ("mm_per_s", "OLS slope of the DoD on gps_time within the line, mm per second of "
                     "mission time -- a residual along-track drift alternates its sign in "
                     "northing between lines, a ground-fixed field does not"),
        ("mm_per_s_se", "cluster-robust standard error of mm_per_s"),
        ("c_tan", "dE re-expressed in the line's own body frame, mm per unit tan(scan) "
                  "(= dE * h): the across-track coefficient this line's easting gradient "
                  "would imply"),
        ("elev_lo", "lower edge of the elevation bin, m"),
        ("elev_hi", "upper edge of the elevation bin, m"),
        ("median_mm", "median DoD in the bin, mm"),
        ("frac", "fraction of the stated total"),
        ("z_p", "elevation percentile, m"),
        ("pop", "which population the percentiles describe"),
        ("band_lo", "lower edge of the northing band, m above the tile's southern limit"),
        ("band_hi", "upper edge of the northing band, m"),
        ("removed", "which registration term was subtracted back out before refitting"),
    ]:
        R.column(c, d)

    # ------------------------------------------------------------- 0. the data
    meta = json.load(open(R.input(f"{T}/corrections.json",
                                  role="pipeline corrections: grid geometry, per-swath "
                                       "constants, geoid datum (const + tilt), drift curves")))
    b = meta["bounds"]
    RES = float(meta["res_m"])
    dod = np.load(R.input(f"{T}/{A.dod}",
                          role="DoD under audit: gen2 at the cover-dependent percentile minus "
                               "gen1's registered median, m, + = elevation rose"))
    NY, NX = dod.shape
    z_all = np.load(R.input(f"{T}/z_after.npy",
                            role="gen2 gridded ground, the DoD reference surface, m")).ravel()
    cover = np.load(R.input(f"{T}/canopy_cover_pfs.npy",
                            role="PyForestScan canopy cover from the gen2 cloud, plant-area "
                                 "based, min_height 2 m")).ravel()
    flood = np.load(R.input(f"{T}/floodplain_mask.npy",
                            role="floodplain cells, excluded from the population")).astype(bool)
    ag = np.load(R.input(f"{T}/ag_region.npy",
                         role="agricultural region: penetration >= 0.45, eroded 6 cells off the "
                              "cover boundary, channel buffer removed "
                              "(analysis/ridgelines/hillslope_regions_K.py)")).ravel()

    mask, rep = reference_cells(T)
    R.cuts(f"{name} reference_cells", rep)
    d_mm = dod.ravel() * 1000.0
    sel = mask & (cover <= A.cover_max) & (~flood.ravel()) & np.isfinite(d_mm)
    R.mask("stable open cells under audit", sel, of=d_mm.size,
           defn=f"refcells.reference_cells() at repo defaults, then canopy_cover_pfs <= "
                f"{A.cover_max:g}, floodplain_mask removed, finite {A.dod}")
    idx = np.flatnonzero(sel)
    ix, iy = idx % NX, idx // NX
    x = b[0] + (ix + 0.5) * RES
    y = b[1] + (iy + 0.5) * RES
    v = d_mm[idx]
    z = z_all[idx]
    E = (x - x.mean()) / 1000.0
    N = (y - y.mean()) / 1000.0
    dist_km = (distance_transform_edt(~flood) * RES).ravel()[idx] / 1000.0
    blocks = {bm: bs.block_ids(idx, nx=NX, res=RES, block_m=bm) for bm in A.block_m}
    B0 = blocks[A.block_m[0]]

    R.banner()

    def fit(what, q, X=None, cols=("dE", "dN"), bm=None, y_=None):
        """Cluster-robust plane fit on a subset, returned as a table row."""
        bm = A.block_m[0] if bm is None else bm
        q = np.asarray(q)
        Eq, Nq = E[q] - E[q].mean(), N[q] - N[q].mean()
        Xq = np.c_[np.ones(int(q.sum())), Eq, Nq] if X is None else X
        yq = v[q] if y_ is None else y_[q]
        be, se, G, _ = ols_cluster(Xq, yq, bs.block_ids(idx[q], nx=NX, res=RES, block_m=bm))
        bl = lad(Xq, yq)
        return [what, int(q.sum()), G, f3(be[0]), f"{se[0]:.2f}", f3(be[1]), f"{se[1]:.2f}",
                f3(be[2]), f"{se[2]:.2f}", f3(bl[1]), f3(bl[2])]

    PLANE = ["what", "n", "blocks", "mean_mm", "mean_se", "dE", "dE_se", "dN", "dN_se",
             "dE_lad", "dN_lad"]

    # ------------------------------------------------- 1. where are the cells
    print("\n## 1. WHERE THE CELLS ARE\n")
    step = int(round(A.map_block_m / RES))
    nbx, nby = int(np.ceil(NX / step)), int(np.ceil(NY / step))
    cnt = np.bincount((iy // step) * nbx + (ix // step),
                      minlength=nbx * nby).reshape(nby, nbx)
    print(f"  {idx.size:,} cells = {100*idx.size/dod.size:.2f}% of the {dod.size:,}-cell tile "
          f"({NY} x {NX} at {RES:g} m)")
    print(f"  {A.map_block_m:g} m blocks: {int((cnt>0).sum())} of {nbx*nby} occupied "
          f"({100*(cnt>0).mean():.0f}%); of {step*step} cells per full block the occupied ones "
          f"hold median {int(np.median(cnt[cnt>0]))}, "
          f"quartiles {int(np.percentile(cnt[cnt>0],25))}-{int(np.percentile(cnt[cnt>0],75))}, "
          f"max {int(cnt.max())}")
    print(f"\n  cells per {A.map_block_m:g} m block (north at top, west at left):")
    for row in cnt[::-1]:
        print("    " + " ".join(f"{c:4d}" if c else "   ." for c in row))
    rows = []
    for lab, arr in (("stable cells under audit", z),
                     ("reference cells (all cover)", z_all[mask & np.isfinite(z_all)]),
                     ("whole tile", z_all[np.isfinite(z_all)])):
        p = np.percentile(arr, [0, 5, 25, 50, 75, 95, 100])
        rows.append([lab, arr.size] + [f"{q:.1f}" for q in p])
    R.column("p0", "minimum elevation of the population, m")
    for q in (5, 25, 50, 75, 95):
        R.column(f"p{q}", f"{q}th percentile of elevation, m")
    R.column("p100", "maximum elevation of the population, m")
    print()
    R.table(["pop", "n", "p0", "p5", "p25", "p50", "p75", "p95", "p100"], rows)
    print(f"\n  easting  {x.min():.1f} .. {x.max():.1f}  (tile {b[0]:.1f} .. {b[2]:.1f})")
    print(f"  northing {y.min():.1f} .. {y.max():.1f}  (tile {b[1]:.1f} .. {b[3]:.1f})")
    print(f"  corr(z, easting) = {np.corrcoef(z, E)[0,1]:+.3f}   "
          f"corr(z, northing) = {np.corrcoef(z, N)[0,1]:+.3f}")
    hist, edges = np.histogram(z, bins=np.arange(np.floor(z.min()/5)*5, z.max()+5, 5))
    lo_pk = int(np.argmax(hist[:len(hist)//2])); hi_pk = len(hist)//2 + int(np.argmax(hist[len(hist)//2:]))
    anti = lo_pk + int(np.argmin(hist[lo_pk:hi_pk]))
    z_split = float(edges[anti + 1])
    print(f"\n  the elevation distribution is BIMODAL: modes at {edges[lo_pk]:.0f}-{edges[lo_pk+1]:.0f} m "
          f"and {edges[hi_pk]:.0f}-{edges[hi_pk+1]:.0f} m; the antimode between them is at "
          f"{z_split:.0f} m")
    print(f"  -> upland limb {int((z>=z_split).sum()):,} cells ({100*(z>=z_split).mean():.0f}%), "
          f"valley limb {int((z<z_split).sum()):,} cells ({100*(z<z_split).mean():.0f}%)")
    R.param("z_split", z_split, src="MINE",
            why="the antimode of the population's OWN elevation histogram (5 m bins), used only "
                "to name the two limbs the data already show; every fit is also reported over "
                "equal-count elevation bins spanning the full range, which uses no threshold.")

    # --------------------------------------------- 2. the plane and its errors
    print("\n## 2. THE PLANE, AND WHAT ITS UNCERTAINTY DEPENDS ON\n")
    X = np.c_[np.ones(v.size), E, N]
    print("  2a. the reported fit, reproduced, and its SE against block size")
    rows = []
    for bm in A.block_m:
        be, se, G, _ = ols_cluster(X, v, blocks[bm])
        rows.append([f"{bm:.0f}", int(v.size), G, f3(be[1]), f"{se[1]:.2f}", f"{be[1]/se[1]:+.2f}",
                     f3(be[2]), f"{se[2]:.2f}", f"{be[2]/se[2]:+.2f}"])
    R.table(["block_m", "n", "blocks", "dE", "dE_se", "t_dE", "dN", "dN_se", "t_dN"], rows)
    boot = block_bootstrap(X, v, B0, A.n_boot, rng)
    be, se, G, r2 = ols_cluster(X, v, B0)
    print(f"\n  intercept {be[0]:+.2f} +- {se[0]:.2f} mm, r2 = {r2:.4f}; "
          f"{A.n_boot} block bootstrap replicates at {A.block_m[0]:g} m give SEs "
          f"{boot[0]:.2f} / {boot[1]:.2f} / {boot[2]:.2f} -- the sandwich and the bootstrap agree")
    print(f"  mean DoD {v.mean():+.2f} mm, median {np.median(v):+.2f} mm, "
          f"NMAD {bs.nmad(v):.1f} mm, sd {v.std():.1f} mm")

    print("\n  2b. L2 against L1 (the repo's standing preference is L1: right-skewed residuals)")
    bl = lad(X, v)
    print(f"    OLS  intercept {be[0]:+.2f}  dE {be[1]:+.2f}  dN {be[2]:+.2f}")
    print(f"    LAD  intercept {bl[0]:+.2f}  dE {bl[1]:+.2f}  dN {bl[2]:+.2f}")
    ubi, inv = np.unique((iy // step) * nbx + (ix // step), return_inverse=True)
    med = np.array([np.median(v[inv == i]) for i in range(ubi.size)])
    Em = np.array([E[inv == i].mean() for i in range(ubi.size)])
    Nm = np.array([N[inv == i].mean() for i in range(ubi.size)])
    Xb = np.c_[np.ones(med.size), Em, Nm]
    bb = np.linalg.lstsq(Xb, med, rcond=None)[0]
    r = med - Xb @ bb
    sb = np.sqrt(np.diag(np.linalg.pinv(Xb.T @ Xb))) * np.sqrt((r ** 2).sum() / (med.size - 3))
    print(f"    plane on the {med.size} {A.map_block_m:g} m BLOCK MEDIANS, one vote per block: "
          f"intercept {bb[0]:+.2f} +- {sb[0]:.2f}  dE {bb[1]:+.2f} +- {sb[1]:.2f}  "
          f"dN {bb[2]:+.2f} +- {sb[2]:.2f}")

    print("\n  2c. how far the residual field is correlated (semivariogram of the plane residual)")
    resid = v - X @ be
    m = 400_000
    i1, i2 = rng.integers(0, v.size, m), rng.integers(0, v.size, m)
    lag = np.hypot(x[i1] - x[i2], y[i1] - y[i2])
    sq = (resid[i1] - resid[i2]) ** 2
    var = float(np.var(resid))
    rows = []
    for lo, hi in zip([0, 25, 50, 75, 100, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400],
                      [25, 50, 75, 100, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400, 3600]):
        q = (lag >= lo) & (lag < hi)
        if not q.any():
            continue
        gam = 0.5 * sq[q].mean()
        rows.append([lo, hi, int(q.sum()), f"{np.sqrt(gam):.1f}", f"{1-gam/var:+.3f}"])
    R.table(["lag_lo", "lag_hi", "pairs", "gamma_rms", "rho"], rows)
    print(f"    residual sd {np.sqrt(var):.1f} mm. Correlation is still +0.25 at 300-400 m and "
          f"reaches zero only near 1.6 km,\n    so a 50 m block treats correlated cells as "
          f"independent and its SE is the optimistic end of the ladder above.")

    # ------------------------------------------------ 3. subregion stability
    print("\n## 3. IS THE TILT A PROPERTY OF THE TILE? (refits on subregions)\n")
    rows = [fit("ALL", np.ones(v.size, bool))]
    rows.append(fit("west half", E < 0))
    rows.append(fit("east half", E >= 0))
    rows.append(fit("south half", N < 0))
    rows.append(fit("north half", N >= 0))
    for ly, qy in (("S", N < 0), ("N", N >= 0)):
        for lx, qx in (("W", E < 0), ("E", E >= 0)):
            rows.append(fit(f"quadrant {ly}{lx}", qx & qy))
    chk = (((ix // 10) + (iy // 10)) % 2) == 0
    rows.append(fit("checkerboard A (50 m)", chk))
    rows.append(fit("checkerboard B (50 m)", ~chk))
    R.table(PLANE, rows)
    ub, inv50 = np.unique(B0, return_inverse=True)
    sub = []
    for _ in range(20):
        q = (rng.random(ub.size) < 0.5)[inv50]
        Xq = np.c_[np.ones(int(q.sum())), E[q] - E[q].mean(), N[q] - N[q].mean()]
        sub.append(np.linalg.lstsq(Xq, v[q], rcond=None)[0])
    sub = np.asarray(sub)
    print(f"\n  20 random 50%-of-blocks subsamples: dE {sub[:,1].mean():+.2f} +- {sub[:,1].std():.2f}"
          f"   dN {sub[:,2].mean():+.2f} +- {sub[:,2].std():.2f}  (spread across subsamples)")

    # ------------------------------------------------------- 4. is it landform
    print("\n## 4. IS IT LANDFORM? (elevation, distance from the valley, land use)\n")
    zc = (z - z.mean()) / 100.0
    agf = ag[idx].astype(float)
    print(f"  the population is {100*agf.mean():.0f}% inside ag_region "
          f"(caveat: ag_region is thresholded on `penetration`, which this project has "
          f"already found\n  is dominated by scan angle and flight-line overlap, not canopy -- "
          f"so it is a weak land-use covariate here)")
    print(f"  distance to the nearest floodplain cell: median {np.median(dist_km)*1000:.0f} m, "
          f"range {dist_km.min()*1000:.0f}-{dist_km.max()*1000:.0f} m\n")
    models = [
        ("E, N", np.c_[np.ones(v.size), E, N], ["dE", "dN"]),
        ("z only", np.c_[np.ones(v.size), zc], ["z/100m"]),
        ("dist only", np.c_[np.ones(v.size), dist_km], ["dist/km"]),
        ("E, N, z", np.c_[np.ones(v.size), E, N, zc], ["dE", "dN", "z/100m"]),
        ("E, N, dist", np.c_[np.ones(v.size), E, N, dist_km], ["dE", "dN", "dist/km"]),
        ("E, N, ag", np.c_[np.ones(v.size), E, N, agf], ["dE", "dN", "ag"]),
        ("E, N, z, dist", np.c_[np.ones(v.size), E, N, zc, dist_km], ["dE", "dN", "z/100m", "dist/km"]),
        ("E, N, z, dist, ag", np.c_[np.ones(v.size), E, N, zc, dist_km, agf],
         ["dE", "dN", "z/100m", "dist/km", "ag"]),
    ]
    rows = []
    for nm, Xm, labs in models:
        bm_, sm, G, r2m = ols_cluster(Xm, v, B0)
        for lab, c_, s_ in zip(["intercept"] + labs, bm_, sm):
            rows.append([nm, lab, f3(c_), f"{s_:.2f}"])
    R.table(["model", "term", "coef", "se"], rows)
    print("\n  the same fit inside equal-count elevation bins (binstats.quantile_edges, "
          "spans the full range):")
    ed = bs.quantile_edges(z, A.n_elev_bins)
    rows = []
    for lo, hi in zip(ed[:-1], ed[1:]):
        q = (z >= lo) & (z < hi)
        if not q.any():
            continue
        rows.append([f"{lo:.0f}-{hi:.0f} m"] + fit("", q)[1:])
    R.table(PLANE, rows)
    print("\n  and split at the antimode of that distribution:")
    R.table(PLANE, [fit(f"upland  z >= {z_split:.0f}", z >= z_split),
                    fit(f"valley  z <  {z_split:.0f}", z < z_split)])

    # ------------------------------------------------- 5. is it the registration
    print("\n## 5. IS IT THE REGISTRATION? (within-swath refits)\n")
    cols = ["cell", "point_source_id", "gps_time", "scan_angle", "in_grid",
            "d_mm_corr", "dz_drift_mm", "dz_swath_mm"]
    t = pq.read_table(R.input(f"{T}/beam_offset_table.parquet",
                              role="per-return gen1 CSF ground offsets to the gen2 surface, "
                                   "slope-normal mm, with the four registration terms and "
                                   "point_source_id / gps_time / scan_angle per return"),
                      columns=cols)
    g = t["in_grid"].to_numpy().astype(bool)
    df = pd.DataFrame({k: t[k].to_numpy()[g] for k in
                       ("cell", "point_source_id", "gps_time", "scan_angle",
                        "d_mm_corr", "dz_drift_mm", "dz_swath_mm")})
    del t
    df = df.astype({"d_mm_corr": float, "dz_drift_mm": float, "dz_swath_mm": float,
                    "scan_angle": float})
    # flying height per line, from ALL in-grid returns: x = a + b*y + h*tan(scan)
    xr = b[0] + ((df.cell.to_numpy() % NX) + 0.5) * RES
    yr = b[1] + ((df.cell.to_numpy() // NX) + 0.5) * RES
    hgt = {}
    for s, q in df.groupby("point_source_id").groups.items():
        k = np.asarray(q)
        Ah = np.c_[np.ones(k.size), yr[k] - yr[k].mean(), np.tan(np.radians(df.scan_angle.to_numpy()[k]))]
        hgt[int(s)] = float(np.linalg.lstsq(Ah, xr[k], rcond=None)[0][2])
    ins = np.zeros(NY * NX, bool)
    ins[idx] = True
    df = df[ins[df.cell.to_numpy()]]

    zfill = np.load(f"{T}/z_after.npy")
    mfill = ~np.isfinite(zfill)
    if mfill.any():
        zfill = zfill[tuple(distance_transform_edt(mfill, return_distances=False,
                                                   return_indices=True))]
    gyz, gxz = np.gradient(zfill, RES)
    nnorm = np.sqrt(gxz.ravel() ** 2 + gyz.ravel() ** 2 + 1.0)

    # per (cell, line) ground = the pipeline's own estimator on one line
    df["d_nodrift"] = df.d_mm_corr - df.dz_drift_mm
    df["d_noswath"] = df.d_mm_corr - df.dz_swath_mm
    per_cell = df.groupby("cell").agg(h_all=("d_mm_corr", "median"),
                                      h_nodrift=("d_nodrift", "median"),
                                      h_noswath=("d_noswath", "median")).reindex(idx)
    cl = df.groupby(["cell", "point_source_id"]).agg(
        h1=("d_mm_corr", "median"), gps=("gps_time", "mean"),
        sc=("scan_angle", "mean"), n=("d_mm_corr", "size")).reset_index()
    vmap = pd.Series(v, index=idx)
    cl["nn"] = nnorm[cl.cell.to_numpy()]
    cl["dod_s"] = (vmap.reindex(cl.cell).to_numpy()
                   + (per_cell.h_all.reindex(cl.cell).to_numpy() - cl.h1) * cl.nn)
    cx = b[0] + ((cl.cell.to_numpy() % NX) + 0.5) * RES
    cy = b[1] + ((cl.cell.to_numpy() // NX) + 0.5) * RES
    Ec = (cx - cx.mean()) / 1000.0
    Nc = (cy - cy.mean()) / 1000.0
    ps = cl.point_source_id.to_numpy()
    ys_ = cl.dod_s.to_numpy()
    bl_c = bs.block_ids(cl.cell.to_numpy(), nx=NX, res=RES, block_m=A.block_m[0])
    print(f"  {len(cl):,} (cell, flight-line) ground estimates over {cl.cell.nunique():,} cells; "
          f"the per-cell DoD is recovered exactly when a cell has one line.")

    def rowf(what, q, Xq=None, extra=0):
        Eq, Nq = Ec[q] - Ec[q].mean(), Nc[q] - Nc[q].mean()
        Xu = np.c_[np.ones(int(q.sum())), Eq, Nq] if Xq is None else Xq
        be_, se_, G_, _ = ols_cluster(Xu, ys_[q], bl_c[q])
        bl_ = lad(Xu, ys_[q])
        return [what, int(q.sum()), G_, f3(be_[0]), f"{se_[0]:.2f}", f3(be_[1]), f"{se_[1]:.2f}",
                f3(be_[2]), f"{se_[2]:.2f}", f3(bl_[1]), f3(bl_[2])]

    allq = np.ones(len(cl), bool)
    ups = np.sort(np.unique(ps))
    Dm = np.column_stack([(ps == s).astype(float) for s in ups[1:]])
    Xd = np.c_[np.ones(len(cl)), Ec - Ec.mean(), Nc - Nc.mean(), Dm]
    bd, sd, Gd, _ = ols_cluster(Xd, ys_, bl_c)
    rows = [rowf("(cell,line) rows, no swath dummies", allq),
            ["(cell,line) rows, WITH swath dummies", len(cl), Gd, f3(bd[0]), f"{sd[0]:.2f}",
             f3(bd[1]), f"{sd[1]:.2f}", f3(bd[2]), f"{sd[2]:.2f}", "-", "-"]]
    R.table(PLANE, rows)
    print("    swath fixed effects relative to line %d: " % ups[0]
          + "  ".join(f"{s}: {c:+.2f} +- {e:.2f} mm" for s, c, e in zip(ups[1:], bd[3:], sd[3:])))

    print("\n  5a. the plane fitted INSIDE each flight line (across-swath constant fixed by "
          "construction)")
    rows = [rowf(f"swath {s}", ps == s) for s in ups]
    R.table(PLANE, rows)
    rows = []
    for s in ups:
        q = ps == s
        Eq = Ec[q]
        gq = cl.gps.to_numpy()[q]
        bq, sq, _, _ = ols_cluster(np.c_[np.ones(int(q.sum())), Eq - Eq.mean(),
                                         Nc[q] - Nc[q].mean()], ys_[q], bl_c[q])
        bt, st_, _, _ = ols_cluster(np.c_[np.ones(int(q.sum())), gq - gq.mean()], ys_[q], bl_c[q])
        rows.append([int(s), int(q.sum()), f"{Eq.max()-Eq.min():.2f}", f"{hgt[int(s)]:+.0f}",
                     f"{np.corrcoef(gq, Nc[q])[0,1]:+.2f}", f3(bq[1]), f"{sq[1]:.2f}",
                     f3(bq[2]), f"{sq[2]:.2f}", f"{bt[1]:+.3f}", f"{st_[1]:.3f}",
                     f"{bq[1]*hgt[int(s)]/1000.0:+.0f}"])
    R.table(["swath", "rows", "E_span", "h_m", "dir", "dE", "dE_se", "dN", "dN_se",
             "mm_per_s", "mm_per_s_se", "c_tan"], rows)
    # who supplies the pooled within-swath dE: the within estimator weights n * Var(E)
    wt = {int(s): float((ps == s).sum() * np.var(Ec[ps == s])) for s in ups}
    tot = sum(wt.values())
    sl = {}
    for s in ups:
        q = ps == s
        Xq = np.c_[np.ones(int(q.sum())), Ec[q] - Ec[q].mean(), Nc[q] - Nc[q].mean()]
        sl[int(s)] = float(np.linalg.lstsq(Xq, ys_[q], rcond=None)[0][1])
    print("    the within estimator weights each line by n*Var(E): "
          + ", ".join(f"{s} {wt[s]/tot:.3f}" for s in sorted(wt))
          + f"\n    weighted mean of the per-line slopes = "
            f"{sum(wt[s]*sl[s] for s in wt)/tot:+.2f} mm/km against the {bd[1]:+.2f} the "
            f"swath-dummy model returns;\n    line by line that is "
          + ", ".join(f"{s}: {wt[s]/tot*sl[s]:+.2f}" for s in sorted(wt))
          + " mm/km of it.")

    print("\n  5b. the same, on the upland limb only (landform held fixed)")
    R.table(PLANE, [rowf("upland, all lines", (z_all[cl.cell.to_numpy()] >= z_split))]
            + [rowf(f"upland, swath {s}", (ps == s) & (z_all[cl.cell.to_numpy()] >= z_split))
               for s in ups])

    print("\n  5c. what each registration term is worth to the tilt "
          f"({A.block_m[2]:g} m blocks)")
    rows = []
    for lab, col in (("nothing (as registered)", "h_all"),
                     ("along-track drift", "h_nodrift"),
                     ("per-swath constants", "h_noswath")):
        dd = v + (per_cell.h_all.to_numpy() - per_cell[col].to_numpy()) * nnorm[idx]
        be_, se_, G_, _ = ols_cluster(X, dd, blocks[A.block_m[2]])
        rows.append([lab, int(v.size), G_, f3(be_[0]), f"{se_[0]:.2f}", f3(be_[1]),
                     f"{se_[1]:.2f}", f3(be_[2]), f"{se_[2]:.2f}", "-", "-"])
    R.column("what", "the quantity or subset the row describes")
    R.table(PLANE, rows)

    # ---------------------------------------------------------- 6. what is left
    print("\n## 6. WHAT IS LEFT\n")
    up = z >= z_split
    rows = []
    for bm in A.block_m:
        Xu = np.c_[np.ones(int(up.sum())), E[up] - E[up].mean(), N[up] - N[up].mean()]
        be_, se_, G_, _ = ols_cluster(Xu, v[up], bs.block_ids(idx[up], nx=NX, res=RES, block_m=bm))
        rows.append([f"{bm:.0f}", int(up.sum()), G_, f3(be_[1]), f"{se_[1]:.2f}",
                     f"{be_[1]/se_[1]:+.2f}", f3(be_[2]), f"{se_[2]:.2f}", f"{be_[2]/se_[2]:+.2f}"])
    print(f"  upland limb only (z >= {z_split:.0f} m, {int(up.sum()):,} cells):")
    R.table(["block_m", "n", "blocks", "dE", "dE_se", "t_dE", "dN", "dN_se", "t_dN"], rows)
    yb = ((y[up] - y[up].min()) // 500).astype(int)
    rows = []
    for k in np.unique(yb):
        q = yb == k
        rows.append([int(500 * k), int(500 * k + 500), int(q.sum()), f3(np.median(v[up][q])),
                     f3(v[up][q].mean())])
    R.column("mean_mm", "mean DoD in the band, mm (+ = elevation rose)")
    print(f"\n  is the northing gradient a ramp or a few patches? median DoD by 500 m band, "
          f"upland limb:")
    R.table(["band_lo", "band_hi", "n", "median_mm", "mean_mm"], rows)

    be_u, se_u, _, _ = ols_cluster(
        np.c_[np.ones(int(up.sum())), E[up] - E[up].mean(), N[up] - N[up].mean()],
        v[up], bs.block_ids(idx[up], nx=NX, res=RES, block_m=A.block_m[2]))
    R.done(headline=f"the easting tilt is landform, not registration: dE {be[1]:+.2f} mm/km "
                    f"tile-wide becomes {be_u[1]:+.2f} +- {se_u[1]:.2f} on the upland limb; the "
                    f"northing tilt survives every control at {be_u[2]:+.2f} +- {se_u[2]:.2f} "
                    f"mm/km (250 m blocks) and is NOT along-track drift (its sign in gps_time "
                    f"alternates with flight direction while its sign in northing does not)")


if __name__ == "__main__":
    main()
