#!/usr/bin/env python3
"""Re-run every validated site with gen2 bare earth from 3DEP's OWN ASPRS ground
classification (``after_ground="class2"``) instead of the last-return heuristic,
and emit the two standard per-site figures:

  A. DEM of Difference (gen2 - gen1) + the level-of-detection raster
  B. the DEM hillshade with the robustly-detected DoD cells drawn at 70% opacity

Every product is saved (dod/lod/change GeoTIFF, .npy arrays, regions + corrections
JSON) so a run is never throwaway. gen1 always uses CSF; gen2 uses class 2 with a
region-level CSF fallback (handled in the pipeline).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python scripts/run_all_sites.py
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python scripts/run_all_sites.py --only battlecreek
"""
import argparse, json, time
from pathlib import Path

import numpy as np
import laspy
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lidar_diff_icp.pipeline import difference_dem
from lidar_diff_icp.detect import detect_change_standard
from lidar_diff_icp.viz import hillshade


def header_bounds(before, res):
    """Frame from the gen1 tile header, ceil/floor to the grid (as the MNRV /
    Whitewater drivers do), for sites without a saved GeoTIFF frame."""
    h = laspy.open(before).header
    X0 = np.ceil(h.mins[0] / res) * res; Y0 = np.ceil(h.mins[1] / res) * res
    X1 = np.floor(h.maxs[0] / res) * res; Y1 = np.floor(h.maxs[1] / res) * res
    return (X0, Y0, X1, Y1)


# name -> (before, after (FULL classified 3DEP cloud), bounds | None=from header, stream)
SITES = {
    "elba": ("data/before/4342-29-64.laz", "data/after/3dep2021_fulltile.laz",
             (577492.8, 4882737.6, 580032.8, 4886237.6), False),
    "whitewater": ("data/before/4358-26-03.laz", "data/after/3dep_4358_fulltile.laz",
                   None, True),
    "mnrv": ("data/before_mnrv/4342-23-01.laz", "data/after_mnrv/mnrv_3dep2021.laz",
             None, True),
    "cook": ("data/before_ne/1158-31-59.laz", "data/after_ne/ne_3dep.laz",
             (709531.0, 5323589.0, 711986.0, 5327144.0), True),
    "carlton": ("data/before_carlton/2742-12-53.laz", "data/after_carlton/carlton_3dep.laz",
                (547805.0, 5163676.0, 550225.0, 5167166.0), True),
    "battlecreek": ("data/before_battlecreek/4342-03-32_b_a.laz",
                    "data/after_battlecreek/battlecreek_3dep.laz",
                    (498750.0, 4975136.0, 499365.0, 4976006.0), False),
}


def _tif(arr, res, x0, y0, ny, out):
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(out, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype="float32", crs="EPSG:26915", nodata=np.nan,
                       transform=from_origin(x0, y0 + ny * res, res, res)) as d:
        d.write(np.flipud(arr).astype("float32"), 1)


def fig_dod_lod(name, Z21, dod, lod, res, X0, Y0, nx, ny, figdir):
    hs = hillshade(Z21, res, X0, Y0, fill_gaps=False)  # nodata -> white, consistent with the LoD panel
    ext = (X0, X0 + nx * res, Y0, Y0 + ny * res); v = 0.3
    fig, ax = plt.subplots(1, 2, figsize=(15, 9))
    ax[0].imshow(hs, extent=ext, origin="lower", cmap="gray", alpha=0.6)
    im0 = ax[0].imshow(dod, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v)
    ax[0].set_title(f"{name}: DEM of Difference (gridded ground): gen2 - gen1 (m)\n"
                    "red = erosion, blue = deposition")
    fig.colorbar(im0, ax=ax[0], shrink=0.6, extend="both")
    im1 = ax[1].imshow(lod, extent=ext, origin="lower", cmap="viridis", vmin=0, vmax=0.2)
    ax[1].set_title("level of detection (m)")
    fig.colorbar(im1, ax=ax[1], shrink=0.6, extend="max")
    for a in ax: a.set_xlabel("Easting (m)"); a.set_ylabel("Northing (m)")
    out = f"{figdir}/{name}_dod_lod.png"
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    return out


def fig_dem_change(name, Z21, dod, change, regions, res, X0, Y0, nx, ny, figdir):
    """DEM hillshade with the robustly-detected DoD cells at 70% opacity."""
    hs = hillshade(Z21, res, X0, Y0, fill_gaps=True)  # gap-filled backdrop, no white holes
    ext = (X0, X0 + nx * res, Y0, Y0 + ny * res); v = 0.3
    over = np.where(change, dod, np.nan)
    net = sum(r["volume_m3"] for r in regions)
    fig, ax = plt.subplots(figsize=(10, 11))
    ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
    im = ax.imshow(over, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v, alpha=0.7)
    ax.set_title(f"{name}: robustly-detected change over the DEM\n"
                 f"{len(regions)} regions, net {net:+,.0f} m3 "
                 "(red = erosion, blue = deposition)")
    ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
    fig.colorbar(im, ax=ax, shrink=0.6, extend="both", label="detected DoD (m)")
    out = f"{figdir}/{name}_change.png"
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    return out


def run_site(name, figdir="figures/rerun_class2"):
    before, after, bounds, stream = SITES[name]
    res = 5.0
    if bounds is None:
        bounds = header_bounds(before, res)
    outdir = Path(f"data/derived/{name}")
    outdir.mkdir(parents=True, exist_ok=True)
    Path(figdir).mkdir(parents=True, exist_ok=True)

    t = time.time()
    print(f"[{name}] difference_dem  bounds={tuple(round(b,1) for b in bounds)} stream={stream}", flush=True)
    r = difference_dem(before, after, bounds, res=res, ground="slope_normal",
                       ground_source="csf", after_ground="class2", stream=stream,
                       robust_stable=True, csf_cache=f"data/csf_cache/{name}.las")
    dod, lod, Z21, stable = r["dod"], r["lod"], r["z_after"], r["stable"]
    nx, ny = r["nx"], r["ny"]; X0, Y0 = r["bounds"][0], r["bounds"][1]
    ex = np.isfinite(dod)
    print(f"[{name}] done in {time.time()-t:.0f}s  stable_sigma={r['stable_sigma']:.3f} m  "
          f"median LoD={np.nanmedian(lod):.3f} m  "
          f"{100*np.mean(np.abs(dod[ex])>lod[ex]):.0f}% cells exceed LoD", flush=True)
    n_null = int((~np.isfinite(Z21)).sum())
    print(f"[{name}] gen2 null cells: {n_null} ({100*n_null/Z21.size:.2f}%) -- water/"
          "dropouts, interpolated in the shaded-relief backdrop ONLY, never on the map "
          "or in the DoD (recorded in regions.json)", flush=True)

    det = detect_change_standard(dod, lod, stable, res)
    change = det["change"]; regions = det["regions"]
    print(f"[{name}] detected {len(regions)} robust regions  "
          f"tau_sys={det['tau_sys_m']:.3f} m  L={det['corr_length_m']:.0f} m", flush=True)

    # persist EVERY output
    _tif(dod, res, X0, Y0, ny, f"{outdir}/dod.tif")
    _tif(lod, res, X0, Y0, ny, f"{outdir}/lod.tif")
    _tif(change.astype("float32"), res, X0, Y0, ny, f"{outdir}/change.tif")
    np.save(f"{outdir}/dod.npy", dod); np.save(f"{outdir}/lod.npy", lod)
    np.save(f"{outdir}/z_after.npy", Z21); np.save(f"{outdir}/change.npy", change)
    with open(f"{outdir}/corrections.json", "w") as fh:
        json.dump(r["corrections"], fh, indent=2)
    with open(f"{outdir}/regions.json", "w") as fh:
        json.dump({**{k: det[k] for k in ("regions", "sigma", "corr_length_m",
                                          "tau_sys_m", "method")},
                   "gen2_null_cells": int((~np.isfinite(Z21)).sum())}, fh, indent=2)

    fa = fig_dod_lod(name, Z21, dod, lod, res, X0, Y0, nx, ny, figdir)
    fb = fig_dem_change(name, Z21, dod, change, regions, res, X0, Y0, nx, ny, figdir)
    print(f"[{name}] wrote {fa}  and  {fb}", flush=True)
    return dict(name=name, sigma=r["stable_sigma"], med_lod=float(np.nanmedian(lod)),
                n_regions=len(regions))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", help="run only these site names")
    ap.add_argument("--figdir", default="figures/rerun_class2")
    a = ap.parse_args()
    names = a.only if a.only else list(SITES)
    summary = []
    for nm in names:
        try:
            summary.append(run_site(nm, figdir=a.figdir))
        except Exception as exc:
            import traceback; traceback.print_exc()
            print(f"[{nm}] FAILED: {exc}", flush=True)
    print("\n=== SUMMARY ===")
    for s in summary:
        print(f"  {s['name']:12s} sigma={s['sigma']:.3f} m  medLoD={s['med_lod']:.3f} m  "
              f"regions={s['n_regions']}")


if __name__ == "__main__":
    main()
