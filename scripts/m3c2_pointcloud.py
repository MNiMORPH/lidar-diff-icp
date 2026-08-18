#!/usr/bin/env python3
"""Point-cloud change: M3C2 between the two epochs' LAST-RETURN (bare-earth)
clouds. Output stays a point cloud (core points with change + LoD), no gridding.

2008 is internally aligned (per-flight-line) and tied to 2021 with the spatially
varying (quadratic) fit; the corrections are applied to the points, then M3C2
runs directly on the last-return clouds.
"""
import argparse
from pathlib import Path
import numpy as np
import laspy
import py4dgeo

from lidar_diff_icp import io, coreg
from lidar_diff_icp.swathdiff import _median_grid


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("before_laz"); ap.add_argument("after_last_laz")
    ap.add_argument("--bounds", nargs=4, type=float, required=True)
    ap.add_argument("--core-res", type=float, default=2.0, help="core-point spacing (m)")
    ap.add_argument("--core-from", choices=("2021", "2008"), default="2021",
                    help="which cloud to hang the core points (and thus the "
                         "cylinders) on; the OTHER cloud is searched inside them")
    ap.add_argument("--correction-surface", action="store_true",
                    help="after the quadratic tie, add a DeLong 400 m vertical "
                         "correction surface fit on upland stable ground (TPI "
                         "floodplain buffer), applied to the 2008 points")
    ap.add_argument("--robust", action="store_true",
                    help="M3C2 median (robust_aggr) instead of mean -- robust to "
                         "above-ground last returns on steep/rough slopes")
    ap.add_argument("--normal-radius", type=float, default=3.0)
    ap.add_argument("--cyl-radius", type=float, default=1.5)
    ap.add_argument("--along-track-drift", action="store_true",
                    help="correct each 2008 point for per-swath along-track GNSS "
                         "drift (coreg.fit_along_track_drift) BEFORE the difference")
    ap.add_argument("--corrections-json", help="write the fitted corrections here")
    ap.add_argument("--out", default="data/derived/change_pointcloud.laz")
    a = ap.parse_args()
    X0, Y0, X1, Y1 = a.bounds
    res = 3.0; nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))

    # 2008: last return, internal align, then quadratic tie to 2021 ground
    f = laspy.read(a.before_laz)
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); cl8 = np.asarray(f.classification)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns)
    gt8 = np.asarray(f.gps_time)
    be8 = rn8 == nr8
    pc = io.PointCloud(x8, y8, z8, ps8, cl8, np.zeros_like(z8), np.zeros_like(ps8), io.MN_GEN1_CRS)
    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()))
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz
    print(f"2008 last-return pts: {be8.sum():,}", flush=True)

    g = laspy.read(a.after_last_laz)
    x2 = np.asarray(g.x); y2 = np.asarray(g.y); z2 = np.asarray(g.z)
    m2 = (x2 >= X0) & (x2 < X1) & (y2 >= Y0) & (y2 < Y1)
    print(f"2021 last-return pts (in bounds): {m2.sum():,}", flush=True)

    # tie on low-percentile (ground) surfaces, apply to 2008 points
    def pg(x, y, z, q):
        import pandas as pd
        ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny); flat = iy[ok] * nx + ix[ok]
        s = pd.Series(z[ok]).groupby(flat).quantile(q)
        out = np.full(nx * ny, np.nan); out[s.index.values] = s.values; return out.reshape(ny, nx)
    mi8 = be8 & (xc >= X0) & (xc < X1) & (yc >= Y0) & (yc < Y1)
    tie = coreg.tie_polynomial(pg(x2[m2], y2[m2], z2[m2], 0.10),
                               pg(xc[mi8], yc[mi8], zc[mi8], 0.10), res, X0, Y0, order=2)
    xc += coreg.eval_poly_field(tie["a"], xc, yc, tie["norm"], 2)
    yc += coreg.eval_poly_field(tie["b"], xc, yc, tie["norm"], 2)
    zc += coreg.eval_poly_field(tie["c"], xc, yc, tie["norm"], 2)
    print(f"tie dx range {np.ptp(tie['dx_field']):.2f} m, residual NMAD {tie['nmad_after']:.3f} m", flush=True)

    # optional DeLong 400 m vertical correction surface (the tie's residual at
    # finer scale). Fit on UPLAND stable ground only -- a topographic position
    # index buffers out the valley floor, since flow accumulation cannot
    # reliably place the channel in a flat floodplain. Applied to 2008 points.
    if a.correction_surface:
        from scipy.ndimage import uniform_filter, distance_transform_edt as edt
        mi8b = be8 & (xc >= X0) & (xc < X1) & (yc >= Y0) & (yc < Y1)
        Z21g = pg(x2[m2], y2[m2], z2[m2], 0.10)
        Z08g = pg(xc[mi8b], yc[mi8b], zc[mi8b], 0.10)
        Zfill = Z21g.copy(); nanm = np.isnan(Zfill)
        if nanm.any():
            Zfill = Zfill[tuple(edt(nanm, return_distances=False, return_indices=True))]
        tpi = Z21g - uniform_filter(Zfill, size=int(2 * 300 / res), mode="nearest")
        floodplain = np.isfinite(Z21g) & (tpi < -2.0)
        cs = coreg.correction_surface(Z21g, Z08g, res, X0, Y0, radius=400.0,
                                      exclude=floodplain)
        C = cs["C"]
        ixp = np.clip(((xc - X0) / res).astype(int), 0, nx - 1)
        iyp = np.clip(((yc - Y0) / res).astype(int), 0, ny - 1)
        Cpt = C[iyp, ixp]; ap_ = np.isfinite(Cpt)
        zc[ap_] += Cpt[ap_]
        print(f"correction surface: floodplain (TPI) {100*floodplain.mean():.0f}% "
              f"excluded from fit; |C| median {np.nanmedian(np.abs(C)):.3f} m; "
              f"applied to {100*ap_.mean():.0f}% of 2008 pts", flush=True)

    # per-swath along-track GNSS-drift correction, applied to EACH 2008 point
    # (via its own gps_time) BEFORE it enters the difference -- the physical,
    # deterministic form of the 2008 registration error, not a post-hoc subtract.
    drift_curves = {}
    if a.along_track_drift:
        from scipy.ndimage import (gaussian_filter, uniform_filter,
                                   distance_transform_edt as edt)
        Z21g = pg(x2[m2], y2[m2], z2[m2], 0.10)
        mi8c = be8 & (xc >= X0) & (xc < X1) & (yc >= Y0) & (yc < Y1)
        resid = Z21g - pg(xc[mi8c], yc[mi8c], zc[mi8c], 0.10)   # ground change (~0 stable)
        Zfill = Z21g.copy(); nanm = np.isnan(Zfill)
        if nanm.any():
            Zfill = Zfill[tuple(edt(nanm, return_distances=False, return_indices=True))]
        tpi = Z21g - uniform_filter(Zfill, size=int(2 * 300 / res), mode="nearest")
        sdeg = np.degrees(coreg.slope_aspect(gaussian_filter(Zfill, 2.0), res)[0])
        Zsm = gaussian_filter(Zfill, 10.0)
        lap = (np.gradient(np.gradient(Zsm, res, axis=0), res, axis=0)
               + np.gradient(np.gradient(Zsm, res, axis=1), res, axis=1))
        stable = ((sdeg < 3) & (tpi > -2)) | ((sdeg > 5) & (sdeg < 35) & (tpi > -2) & (lap < 0))
        ixp = np.clip(((xc - X0) / res).astype(int), 0, nx - 1)
        iyp = np.clip(((yc - Y0) / res).astype(int), 0, ny - 1)
        chg_pt = resid[iyp, ixp]
        stab_pt = be8 & stable[iyp, ixp] & np.isfinite(chg_pt) & (np.abs(chg_pt) < 0.15)
        drift, curves = coreg.fit_along_track_drift(gt8, chg_pt, stab_pt, ps8)
        zc += drift                                            # per-point, pre-difference
        drift_curves = {p: (c[0].tolist(), c[1].tolist()) for p, c in curves.items()}
        print(f"along-track drift: {len(curves)} swaths, drift NMAD "
              f"{1.4826*np.median(np.abs(drift-np.median(drift))):.3f} m; per-swath "
              f"amplitude " + ", ".join(f"{p}:{(c[1].max()-c[1].min()):.3f}" for p, c in curves.items()),
              flush=True)

    p08 = np.column_stack([xc[be8], yc[be8], zc[be8]]).astype(np.float64)
    p21 = np.column_stack([x2[m2], y2[m2], z2[m2]]).astype(np.float64)

    # core points: subsample the chosen cloud onto a grid. The core cloud is
    # where the cylinders (and normals) are hung; the OTHER cloud is searched
    # inside them. Sign is unchanged -- epochs stay (2008, 2021), so positive =
    # 2021 higher regardless of which cloud carries the core points.
    src = p21 if a.core_from == "2021" else p08
    k = (np.floor(src[:, 0] / a.core_res) * 1e6 + np.floor(src[:, 1] / a.core_res)).astype(np.int64)
    _, idx = np.unique(k, return_index=True); core = src[idx]
    print(f"core points ({a.core_res} m) from {a.core_from}: {len(core):,}", flush=True)

    # robust_aggr=True -> median (not mean) of the cylinder's points along the
    # normal. On steep/rough slopes the mean is pulled up by above-ground last
    # returns (scatter in the upper tail); the median ignores that tail up to
    # ~50% contamination, so it tracks the ground surface better.
    m3 = py4dgeo.M3C2(epochs=(py4dgeo.Epoch(p08), py4dgeo.Epoch(p21)), corepoints=core,
                      normal_radii=(a.normal_radius,), cyl_radius=a.cyl_radius,
                      max_distance=15.0, registration_error=0.0, robust_aggr=a.robust)
    dist, unc = m3.run(); lod = unc["lodetection"]
    ok = np.isfinite(dist)
    print(f"M3C2 done: {ok.sum():,} core points with a value; median LoD {np.nanmedian(lod):.3f} m", flush=True)

    # write core points + change + LoD as a point cloud (extra dims)
    h = laspy.LasHeader(point_format=3, version="1.2")
    h.offsets = core.min(0); h.scales = [0.01, 0.01, 0.01]
    h.add_extra_dim(laspy.ExtraBytesParams(name="m3c2", type=np.float32))
    h.add_extra_dim(laspy.ExtraBytesParams(name="lod", type=np.float32))
    las = laspy.LasData(h)
    las.x, las.y, las.z = core[:, 0], core[:, 1], core[:, 2]
    las.m3c2 = dist.astype(np.float32); las.lod = lod.astype(np.float32)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    las.write(a.out)
    print(f"wrote {a.out} ({len(core):,} points, dims: x y z m3c2 lod)", flush=True)

    if a.corrections_json:
        import json
        out = {
            "epochs": "2021 - 2008 (positive = 2021 higher)",
            "crs": "EPSG:26915", "bounds": list(a.bounds),
            "per_swath_internal_alignment_dxdydz_m":
                {str(k): [round(float(v), 4) for v in val] for k, val in corr.items()},
            "cross_epoch_tie": {"order": 2,
                "dx_coef": [round(float(v), 6) for v in tie["a"]],
                "dy_coef": [round(float(v), 6) for v in tie["b"]],
                "dz_coef": [round(float(v), 6) for v in tie["c"]],
                "normalization_xm_xhr_ym_yhr": [round(float(v), 3) for v in tie["norm"]]},
            "along_track_drift_curves_gpsTime_to_metres":
                {str(p): {"gps_time": [round(t, 3) for t in c[0]],
                          "drift_m": [round(d, 4) for d in c[1]]}
                 for p, c in drift_curves.items()},
        }
        with open(a.corrections_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"wrote corrections to {a.corrections_json}", flush=True)


if __name__ == "__main__":
    main()
