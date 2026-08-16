#!/usr/bin/env python3
"""CLI for the final gridded-ground DEM of Difference (2008 MN lidar vs 3DEP).

Thin wrapper over ``lidar_diff_icp.pipeline.difference_dem`` (which holds the
validated workflow and its lessons). Writes GeoTIFFs (dod, lod), a corrections
JSON, and a standard NW-hillshade figure. Run with PROJ_DATA UNSET so the pip
rasterio uses its bundled PROJ:

    env -u PROJ_DATA -u GDAL_DATA python scripts/gridded_ground_dod.py \
      data/before/4342-29-64.laz data/after/3dep2021_last.laz \
      --bounds 577492.8 4882737.6 580035.0 4886238.3
"""
import argparse, json
from pathlib import Path
import numpy as np

from lidar_diff_icp.pipeline import difference_dem


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("before_laz"); ap.add_argument("after_last_laz")
    ap.add_argument("--bounds", nargs=4, type=float, required=True)
    ap.add_argument("--res", type=float, default=5.0)
    ap.add_argument("--ground-q", type=float, default=0.10)
    ap.add_argument("--no-correction-surface", action="store_true")
    ap.add_argument("--no-drift", action="store_true")
    ap.add_argument("--outdir", default="data/derived/final")
    ap.add_argument("--figdir", default="figures")
    a = ap.parse_args()

    r = difference_dem(a.before_laz, a.after_last_laz, a.bounds, res=a.res,
                       ground_q=a.ground_q,
                       correction_surface=not a.no_correction_surface,
                       along_track_drift=not a.no_drift)
    dod, lod = r["dod"], r["lod"]
    ex = np.isfinite(dod)
    print(f"gridded-ground DoD @ {a.res:.0f} m: stable 1-sigma {r['stable_sigma']:.3f} m; "
          f"median LoD {np.nanmedian(lod):.3f} m; "
          f"{100*np.mean(np.abs(dod[ex]) > lod[ex]):.0f}% of cells exceed LoD", flush=True)

    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    _tif(dod, r["res"], *r["bounds"][:2], r["ny"], f"{a.outdir}/dod.tif")
    _tif(lod, r["res"], *r["bounds"][:2], r["ny"], f"{a.outdir}/lod.tif")
    with open(f"{a.outdir}/corrections.json", "w") as fh:
        json.dump(r["corrections"], fh, indent=2)
    print(f"wrote {a.outdir}/dod.tif, lod.tif, corrections.json", flush=True)
    _fig(r["z_after"], dod, lod, r["res"], r["bounds"][0], r["bounds"][1], r["nx"], r["ny"], a.figdir)


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
    ax[0].set_title("DEM of Difference (gridded ground), after - before (m)\nred = erosion, blue = deposition")
    fig.colorbar(im0, ax=ax[0], shrink=0.6, extend="both")
    im1 = ax[1].imshow(lod, extent=ext, origin="lower", cmap="viridis", vmin=0, vmax=0.2)
    ax[1].set_title("level of detection (m)")
    fig.colorbar(im1, ax=ax[1], shrink=0.6, extend="max")
    for a in ax: a.set_xlabel("Easting (m)"); a.set_ylabel("Northing (m)")
    fig.savefig(f"{figdir}/final_dod.png", dpi=130, bbox_inches="tight")
    print(f"wrote {figdir}/final_dod.png", flush=True)


if __name__ == "__main__":
    main()
