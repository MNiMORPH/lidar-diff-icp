#!/usr/bin/env python3
"""Final differenced product: gridded low-10% bare-earth DEM of Difference.

Chosen over point-based M3C2-median because a controlled comparison (same grid,
mask, and corrections) showed M3C2-median carries a coherent slope bias (~16% of
convex hillslope falsely depositional) while the gridded low-10% ground surface
does not (~4%). For coherent change detection the coherent bias is the dangerous
error -- it masquerades as change and does not average out -- whereas the ground
DoD's higher point noise is incoherent and is carried by the LoD.

Ground = 10th-percentile last-return elevation per cell (consistent between
epochs, so intra-cell relief cancels in the difference). 2008 is corrected in the
acquisition-honest order: per-swath internal alignment -> spatially varying
quadratic tie -> DeLong 400 m correction surface (flats, TPI floodplain buffer)
-> per-swath along-track GNSS-drift spline f(gps_time), applied to each point
before gridding. DoD = 2021 - 2008 (positive = deposition). Cell size is set by
the sparse 2008 density (~0.8 pts/m^2 -> ~20 pts per 5 m cell for a stable
percentile). Outputs GeoTIFFs (dod, lod), a corrections JSON, and a figure.

    env PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal \
      python scripts/gridded_ground_dod.py data/before/4342-29-64.laz \
      data/after/3dep2021_last.laz --bounds 577492.8 4882737.6 580035.0 4886238.3
"""
import argparse, json
from pathlib import Path
import numpy as np
import laspy
import pandas as pd
from scipy.ndimage import gaussian_filter, uniform_filter, distance_transform_edt as edt

from lidar_diff_icp import io, coreg


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("before_laz"); ap.add_argument("after_last_laz")
    ap.add_argument("--bounds", nargs=4, type=float, required=True)
    ap.add_argument("--res", type=float, default=5.0)
    ap.add_argument("--ground-q", type=float, default=0.10)
    ap.add_argument("--outdir", default="data/derived/final")
    ap.add_argument("--figdir", default="figures")
    a = ap.parse_args()
    X0, Y0, X1, Y1 = a.bounds
    res = a.res; nx = int(round((X1 - X0) / res)); ny = int(round((Y1 - Y0) / res))

    def cellstat(x, y, z, how):
        ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        gb = pd.Series(z[ok]).groupby(iy[ok] * nx + ix[ok])
        if how == "ground": s = gb.quantile(a.ground_q)
        elif how == "spread": s = 1.4826 * (gb.quantile(0.75) - gb.quantile(0.25)) / 1.349
        elif how == "count": s = gb.size()
        out = np.full(nx * ny, np.nan); out[s.index.values] = s.values
        return out.reshape(ny, nx)

    # --- 2021 reference ground + its within-cell spread/count ---
    g = laspy.read(a.after_last_laz)
    x2 = np.asarray(g.x); y2 = np.asarray(g.y); z2 = np.asarray(g.z)
    m2 = (x2 >= X0) & (x2 < X1) & (y2 >= Y0) & (y2 < Y1)
    Z21 = cellstat(x2[m2], y2[m2], z2[m2], "ground")
    s21 = cellstat(x2[m2], y2[m2], z2[m2], "spread"); n21 = cellstat(x2[m2], y2[m2], z2[m2], "count")

    # terrain masks (from 2021 ground)
    Zf = Z21.copy(); nanm = np.isnan(Zf)
    if nanm.any():
        Zf = Zf[tuple(edt(nanm, return_distances=False, return_indices=True))]
    tpi = Z21 - uniform_filter(Zf, size=int(2 * 300 / res), mode="nearest")
    sdeg = np.degrees(coreg.slope_aspect(gaussian_filter(Zf, 2.0), res)[0])
    Zsm = gaussian_filter(Zf, 50 / res / 2)
    lap = (np.gradient(np.gradient(Zsm, res, axis=0), res, axis=0)
           + np.gradient(np.gradient(Zsm, res, axis=1), res, axis=1))
    convex = (sdeg > 5) & (sdeg < 35) & (tpi > -2) & (lap < 0)
    stable = ((sdeg < 3) & (tpi > -2)) | convex
    floodplain = np.isfinite(Z21) & (tpi < -2)

    # --- 2008: internal align -> tie -> correction surface -> along-track drift ---
    f = laspy.read(a.before_laz)
    x8 = np.asarray(f.x); y8 = np.asarray(f.y); z8 = np.asarray(f.z)
    ps8 = np.asarray(f.point_source_id); gt8 = np.asarray(f.gps_time)
    rn8 = np.asarray(f.return_number); nr8 = np.asarray(f.number_of_returns); be = rn8 == nr8
    pc = io.PointCloud(x8, y8, z8, ps8, np.asarray(f.classification),
                       np.zeros_like(z8), np.zeros_like(ps8), io.MN_2008_CRS)
    corr, _, _ = coreg.align_swaths(pc, ref=int(ps8.min()))
    xc, yc, zc = x8.copy(), y8.copy(), z8.copy()
    for s, (dx, dy, dz) in corr.items():
        m = ps8 == s; xc[m] += dx; yc[m] += dy; zc[m] += dz
    Z08 = cellstat(xc[be], yc[be], zc[be], "ground")
    tie = coreg.tie_polynomial(Z21, Z08, res, X0, Y0, order=2)
    xc += coreg.eval_poly_field(tie["a"], xc, yc, tie["norm"], 2)
    yc += coreg.eval_poly_field(tie["b"], xc, yc, tie["norm"], 2)
    zc += coreg.eval_poly_field(tie["c"], xc, yc, tie["norm"], 2)
    cs = coreg.correction_surface(Z21, cellstat(xc[be], yc[be], zc[be], "ground"),
                                  res, X0, Y0, radius=400.0, exclude=floodplain)
    ixp = np.clip(((xc - X0) / res).astype(int), 0, nx - 1)
    iyp = np.clip(((yc - Y0) / res).astype(int), 0, ny - 1)
    Cpt = cs["C"][iyp, ixp]; zc[np.isfinite(Cpt)] += Cpt[np.isfinite(Cpt)]
    resid = Z21 - cellstat(xc[be], yc[be], zc[be], "ground")
    chg_pt = resid[iyp, ixp]
    stab_pt = be & stable[iyp, ixp] & np.isfinite(chg_pt) & (np.abs(chg_pt) < 0.15)
    drift, curves = coreg.fit_along_track_drift(gt8, chg_pt, stab_pt, ps8)
    zc += drift

    # --- final gridded ground DoD + per-cell LoD ---
    Z08c = cellstat(xc[be], yc[be], zc[be], "ground")
    s08 = cellstat(xc[be], yc[be], zc[be], "spread"); n08 = cellstat(xc[be], yc[be], zc[be], "count")
    dod = Z21 - Z08c
    lod = 1.96 * np.sqrt(np.nan_to_num(s08**2 / np.maximum(n08, 1))
                         + np.nan_to_num(s21**2 / np.maximum(n21, 1)))
    r = dod[stable & np.isfinite(dod)]
    emp = 1.4826 * np.median(np.abs(r - np.median(r)))
    print(f"gridded-ground DoD @ {res:.0f} m: stable NMAD {emp:.3f} m (empirical 1-sigma); "
          f"median LoD {np.nanmedian(lod):.3f} m; {100*np.mean(np.abs(dod[np.isfinite(dod)])>lod[np.isfinite(dod)]):.0f}% of cells exceed LoD", flush=True)

    # --- write GeoTIFFs + corrections + figure ---
    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    _tif(dod, res, X0, Y0, ny, f"{a.outdir}/dod.tif")
    _tif(lod, res, X0, Y0, ny, f"{a.outdir}/lod.tif")
    with open(f"{a.outdir}/corrections.json", "w") as fh:
        json.dump({"epochs": "2021 - 2008 (positive = deposition)", "crs": "EPSG:26915",
                   "res_m": res, "bounds": list(a.bounds), "stable_1sigma_m": round(float(emp), 4),
                   "per_swath_internal_alignment_dxdydz_m":
                       {str(k): [round(float(v), 4) for v in val] for k, val in corr.items()},
                   "cross_epoch_tie_order2_coef": {"dx": [round(float(v), 6) for v in tie["a"]],
                       "dy": [round(float(v), 6) for v in tie["b"]], "dz": [round(float(v), 6) for v in tie["c"]],
                       "norm_xm_xhr_ym_yhr": [round(float(v), 3) for v in tie["norm"]]},
                   "along_track_drift_gpsTime_to_m":
                       {str(p): {"gps_time": [round(t, 3) for t in c[0]], "drift_m": [round(d, 4) for d in c[1]]}
                        for p, c in curves.items()}}, fh, indent=2)
    print(f"wrote {a.outdir}/dod.tif, lod.tif, corrections.json", flush=True)
    _fig(Z21, dod, lod, res, X0, Y0, nx, ny, a.figdir)


def _tif(arr, res, x0, y0, ny, out):
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(out, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype="float32", crs="EPSG:26915", nodata=np.nan,
                       transform=from_origin(x0, y0 + ny * res, res, res)) as d:
        d.write(np.flipud(arr).astype("float32"), 1)


def _fig(Z21, dod, lod, res, X0, Y0, nx, ny, figdir):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    Path(figdir).mkdir(exist_ok=True)
    ls = LightSource(azdeg=315, altdeg=45)
    q = np.nan_to_num(Z21, nan=np.nanmin(Z21))
    hs = np.flipud(ls.hillshade(np.flipud(q), vert_exag=2, dx=res, dy=res))  # standard NW
    ext = (X0, X0 + nx * res, Y0, Y0 + ny * res); v = 0.3
    fig, ax = plt.subplots(1, 2, figsize=(15, 9))
    ax[0].imshow(hs, extent=ext, origin="lower", cmap="gray", alpha=0.6)
    im0 = ax[0].imshow(dod, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v)
    ax[0].set_title("DEM of Difference (gridded ground), 2021 - 2008 (m)\nred = erosion, blue = deposition")
    fig.colorbar(im0, ax=ax[0], shrink=0.6, extend="both")
    im1 = ax[1].imshow(lod, extent=ext, origin="lower", cmap="viridis", vmin=0, vmax=0.2)
    ax[1].set_title("level of detection (m)")
    fig.colorbar(im1, ax=ax[1], shrink=0.6, extend="max")
    for a in ax: a.set_xlabel("Easting (m)"); a.set_ylabel("Northing (m)")
    fig.savefig(f"{figdir}/final_dod.png", dpi=130, bbox_inches="tight")
    print(f"wrote {figdir}/final_dod.png", flush=True)


if __name__ == "__main__":
    main()
