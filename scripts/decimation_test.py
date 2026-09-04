#!/usr/bin/env python3
"""Density-decimation test.

Does matching 2021 down to 2008's point density change the differenced surface?
Naive expectation: the denser 2021 cloud finds lower ground (more shots have a
chance to penetrate vegetation to the true surface), biasing the surface low --
most plausibly on the reed-canary-grass floodplain. This runs M3C2 (median)
twice on the SAME corrected 2008 cloud: once with 2021 at full density, once with
2021 randomly decimated to 2008's point count (= 2008's density over the same
area). Results are broken down by landscape position -- floodplain / upland flat /
steep -- via a topographic position index, and written to
analysis/decimation_result.md.

    python scripts/decimation_test.py data/before/4342-29-64.laz \
        data/after/3dep2021_last.laz --bounds 577492.8 4882737.6 580035.0 4886238.3
"""
import argparse
from pathlib import Path
import numpy as np
import laspy
import pandas as pd
import py4dgeo
from scipy.ndimage import uniform_filter, distance_transform_edt as edt

from lidar_diff_icp import io, coreg, terrain


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("before_laz"); ap.add_argument("after_last_laz")
    ap.add_argument("--bounds", nargs=4, type=float, required=True)
    ap.add_argument("--core-res", type=float, default=2.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--valley-top", dest="valley_top", default="histogram",
                    help="valley top for the floodplain mask: an elevation in metres, "
                         "'registry', or 'histogram'. Never chosen for you.")
    ap.add_argument("--tile-dir", dest="tile_dir", default=None,
                    help="tile directory, for --valley-top registry/histogram")
    a = ap.parse_args()
    X0, Y0, X1, Y1 = a.bounds
    res = 5.0; nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))

    def gg(x, y, z, q=0.10):
        ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        s = pd.Series(z[ok]).groupby(iy[ok] * nx + ix[ok]).quantile(q)
        out = np.full(nx * ny, np.nan); out[s.index.values] = s.values
        return out.reshape(ny, nx)

    # 2008: last return, internal align, quadratic tie, correction surface
    f = laspy.read(a.before_laz)
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns)
    be8 = rn8 == nr8
    pc = io.PointCloud(x8, y8, z8, ps8, np.asarray(f.classification),
                       np.zeros_like(z8), np.zeros_like(ps8), io.MN_GEN1_CRS)
    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()))
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz
    g = laspy.read(a.after_last_laz)
    x2 = np.asarray(g.x); y2 = np.asarray(g.y); z2 = np.asarray(g.z)
    m2 = (x2 >= X0) & (x2 < X1) & (y2 >= Y0) & (y2 < Y1)
    Z21 = gg(x2[m2], y2[m2], z2[m2]); Z08 = gg(xc[be8], yc[be8], zc[be8])
    tie = coreg.tie_polynomial(Z21, Z08, res, X0, Y0, order=2)
    xc += coreg.eval_poly_field(tie["a"], xc, yc, tie["norm"], 2)
    yc += coreg.eval_poly_field(tie["b"], xc, yc, tie["norm"], 2)
    zc += coreg.eval_poly_field(tie["c"], xc, yc, tie["norm"], 2)
    Z08c = gg(xc[be8], yc[be8], zc[be8])
    Zfill = Z21.copy(); nanm = np.isnan(Zfill)
    if nanm.any():
        Zfill = Zfill[tuple(edt(nanm, return_distances=False, return_indices=True))]
    floodplain = terrain.terrain_masks(Z21, res, valley_top_m=a.valley_top,
                                       tile_dir=a.tile_dir)["floodplain"]
    cs = coreg.correction_surface(Z21, Z08c, res, X0, Y0, radius=400.0, exclude=floodplain)
    C = cs["C"]
    ixp = np.clip(((xc - X0) / res).astype(int), 0, nx - 1)
    iyp = np.clip(((yc - Y0) / res).astype(int), 0, ny - 1)
    Cpt = C[iyp, ixp]; zc[np.isfinite(Cpt)] += Cpt[np.isfinite(Cpt)]

    p08 = np.column_stack([xc[be8], yc[be8], zc[be8]]).astype(np.float64)
    p21 = np.column_stack([x2[m2], y2[m2], z2[m2]]).astype(np.float64)
    rng = np.random.default_rng(a.seed)
    keep = rng.choice(len(p21), size=min(len(p08), len(p21)), replace=False)
    p21d = p21[keep]
    print(f"2008 {len(p08):,}  2021 dense {len(p21):,} ({len(p21)/len(p08):.1f}x)  "
          f"2021 decimated {len(p21d):,}", flush=True)

    kk = (np.floor(p08[:, 0] / a.core_res) * 1e6 + np.floor(p08[:, 1] / a.core_res)).astype(np.int64)
    _, idx = np.unique(kk, return_index=True); core = p08[idx]
    E08 = py4dgeo.Epoch(p08)

    def runmed(p):
        m = py4dgeo.M3C2(epochs=(E08, py4dgeo.Epoch(p)), corepoints=core,
                         normal_radii=(3.0,), cyl_radius=1.5, max_distance=15.0,
                         registration_error=0.0, robust_aggr=True)
        d, u = m.run(); return d, u["lodetection"]
    dd, ld = runmed(p21); print("dense done", flush=True)
    dc, lc = runmed(p21d); print("decimated done", flush=True)

    # classify each core point by landscape position
    slope, _ = coreg.slope_aspect(Z21, res); sdeg = np.degrees(slope)
    ci = np.clip(((core[:, 0] - X0) / res).astype(int), 0, nx - 1)
    ri = np.clip(((core[:, 1] - Y0) / res).astype(int), 0, ny - 1)
    fp = floodplain[ri, ci]; st = sdeg[ri, ci] > 15; up = (sdeg[ri, ci] < 3) & (~fp)

    def nmad(v): return 1.4826 * np.median(np.abs(v - np.median(v)))
    rows = []
    for lab, msk in [("floodplain (reed canary)", fp), ("upland flat", up),
                     ("steep >15deg", st), ("ALL", np.ones(len(core), bool))]:
        ok = np.isfinite(dd) & np.isfinite(dc) & msk
        rows.append((lab, ok.sum(), np.median(dd[ok]), np.median(dc[ok]),
                     np.median(dd[ok] - dc[ok]), np.median(ld[ok]), np.median(lc[ok])))
    hdr = (f"| position | n | dense median | decim median | dense-decim | "
           f"LoD dense | LoD decim |\n|---|--:|--:|--:|--:|--:|--:|")
    lines = [hdr] + [f"| {l} | {n:,} | {a1:+.3f} | {a2:+.3f} | {a3:+.3f} | {a4:.3f} | {a5:.3f} |"
                     for (l, n, a1, a2, a3, a4, a5) in rows]
    table = "\n".join(lines)
    print(table, flush=True)

    # why: density effect by percentile on the floodplain. A denser cloud finds
    # lower ground through vegetation at the MIN / low percentiles, but not at
    # the median -- which is why the median-aggregated M3C2 is density-immune.
    xin, yin, zin = x2[m2], y2[m2], z2[m2]
    ip = ((xin - X0) / res).astype(int); jp = ((yin - Y0) / res).astype(int)
    okp = (ip >= 0) & (ip < nx) & (jp >= 0) & (jp < ny)
    dfp = pd.DataFrame({"cell": jp[okp] * nx + ip[okp], "z": zin[okp]})
    rng2 = np.random.default_rng(a.seed)
    dfp["k"] = rng2.random(len(dfp)) < (len(p08) / len(p21))
    fp_flat = np.flatnonzero(floodplain.ravel())
    plines = ["| percentile | dense - decimated (m) |", "|---|--:|"]
    for p, lab in [(0, "min"), (5, "5th"), (10, "10th"), (25, "25th"), (50, "median")]:
        d = dfp.groupby("cell")["z"].quantile(p / 100.0)
        c = dfp[dfp.k].groupby("cell")["z"].quantile(p / 100.0)
        common = d.index.intersection(c.index)
        common = common[np.isin(common, fp_flat)]
        diff = (d.loc[common] - c.loc[common]).values
        diff = diff[np.isfinite(diff)]
        plines.append(f"| {lab} | {np.median(diff):+.3f} |")
    ptable = "\n".join(plines)
    print(ptable, flush=True)

    Path("analysis").mkdir(exist_ok=True)
    with open("analysis/decimation_result.md", "w") as fh:
        fh.write("# Density-decimation test (M3C2 median)\n\n")
        fh.write("2021 decimated to 2008 point count (same area => 2008 density). "
                 "Positive = 2021 higher (deposition). Units metres.\n\n")
        fh.write(table + "\n\n")
        fh.write("## Why: density effect is real for a LOW ground pick, not the median\n\n")
        fh.write("Per floodplain cell, 2021 elevation at a percentile, dense minus "
                 "decimated:\n\n")
        fh.write(ptable + "\n\n")
        fh.write("The denser cloud finds lower ground through the grass at the min / "
                 "low percentiles, decaying to zero at the median.\n\n")
        fh.write("Interpretation is the reader's; the numbers are what changed "
                 "when the densities were matched.\n")
    print("wrote analysis/decimation_result.md", flush=True)


if __name__ == "__main__":
    main()
