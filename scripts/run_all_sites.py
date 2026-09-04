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
from lidar_diff_icp import figures, refcells
from lidar_diff_icp.viz import hillshade


def header_bounds(before, res):
    """Frame from the gen1 tile header, ceil/floor to the grid (as the MNRV /
    Whitewater drivers do), for sites without a saved GeoTIFF frame."""
    h = laspy.open(before).header
    X0 = np.ceil(h.mins[0] / res) * res; Y0 = np.ceil(h.mins[1] / res) * res
    X1 = np.floor(h.maxs[0] / res) * res; Y1 = np.floor(h.maxs[1] / res) * res
    return (X0, Y0, X1, Y1)


# name -> (before, after (FULL classified 3DEP cloud), bounds | None=from header, stream)
#: Where the two standard per-site figures go. It was "figures/rerun_class2" -- the name of
#: the 2026-08-19 experiment that first ran every site with ASPRS class-2 gen2 ground
#: (cd14c35). class2 has been the default ever since, so the name recorded a question
#: settled two weeks ago while the directory held the CURRENT products for all six sites.
#: Named for its content now.
FIGDIR = figures.DEFAULT_FIGDIR

SITES = {
    # gen2 was 3dep2021_fulltile.laz until 2026-09-02: 16,785,971 points, 1.89 pts/m2,
    # against 182,923,322 (20.55 pts/m2) sitting beside it in the same extent -- a 10.90x
    # thinning, and 19.27x on the ground class (0.30 vs 5.78 pts/m2). Every other site
    # already read its full cloud, so Elba alone was differenced from something 6.03-27.53x
    # thinner than its peers: a cross-site method artifact before it is an uncertainty
    # problem. It also put 109,894 of 339,950 DoD cells (32.33%) below cell_plane_roughness
    # min_n=6, so they lost their LoD entirely. stream=True because 182.9M points will not
    # fit in memory whole; it is the same path the other four large sites use.
    "elba": ("data/before/4342-29-64.laz", "data/after/3dep2021_fulldensity.laz",
             (577492.8, 4882737.6, 580032.8, 4886237.6), True),
    # 3dep_4358_fulltile.laz was a TRUNCATED fetch: 5.52 returns/m2 east of easting 586362
    # against 15.45 west, same five flight lines both sides. Replaced 2026-09-04 by an
    # uncapped re-fetch, 148,050,625 points, west/east density ratio 1.06.
    "whitewater": ("data/before/4358-26-03.laz", "data/after/3dep_4358_fulldensity.laz",
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


#: Where each site's valley top comes from. STATED per site, never defaulted: "registry"
#: uses refcells.VALLEY_TOP_M (an established, cited elevation); "histogram" computes it
#: from the landscape's pooled elevations. Elba's landscape has a cited 230.0 m, so it uses
#: it; the others have none, and say so by asking for the computation.
VALLEY_TOP_SOURCE = {
    "elba": "registry",
    "whitewater": "histogram",
    "mnrv": "histogram",
    "cook": "histogram",
    "carlton": "histogram",
    "battlecreek": "histogram",
}


def run_site(name, figdir=FIGDIR, *, skip_penetration=False):
    before, after, bounds, stream = SITES[name]
    res = 5.0
    if bounds is None:
        bounds = header_bounds(before, res)
    outdir = Path(f"data/derived/{name}")
    outdir.mkdir(parents=True, exist_ok=True)
    Path(figdir).mkdir(parents=True, exist_ok=True)

    t = time.time()
    print(f"[{name}] difference_dem  bounds={tuple(round(b,1) for b in bounds)} stream={stream}", flush=True)
    # The valley cut is STATED here, per site, and never chosen by the pipeline:
    # "registry" for a landscape with an established, cited elevation; "histogram" to
    # compute it from the landscape's pooled elevations. Andy, 2026-09-04: we must always
    # tell it which.
    vt = VALLEY_TOP_SOURCE[name]
    r = difference_dem(before, after, bounds, res=res, ground="slope_normal",
                       ground_source="csf", after_ground="class2", stream=stream,
                       robust_stable=True, csf_cache=f"data/csf_cache/{name}.las",
                       valley_top_m=(refcells.VALLEY_TOP_M[name] if vt == "registry" else vt),
                       tile_dir_for_landscape=str(outdir))
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

    # gen2 ground penetration. OPTIONAL, and it drives nothing -- it is saved for
    # inspection only.
    #
    # The leaf-on/forest-slope flag that used to be derived here was RETIRED on 2026-09-02
    # (Andy). It combined two undefended thresholds, penetration < 0.25 AND slope > 12 deg,
    # and once its LoD widening was switched off it changed no output: nothing in the
    # library, the scripts or analysis/ read leafon_flag.npy, and the detector takes no
    # vegetation term. It was 98,143 cells (28%) at Elba that looked like a filter and was
    # not one. The mechanism it rested on (gen2 ground STARVED under leaf-on canopy) is
    # also contradicted by later work here: gen2 carries ~12x MORE ground than gen1, and
    # the offset tracks the 2008 canopy rather than 2021 cover. See canopy.py's status
    # section. Vegetation belongs in the cover-matched q2 chain, not in a flag on the LoD.
    if skip_penetration:
        print(f"[{name}] penetration: SKIPPED as stated (--no-penetration)", flush=True)
    else:
        from lidar_diff_icp.canopy import ground_penetration
        pen = ground_penetration(after, r["bounds"], res, nx, ny)
        print(f"[{name}] gen2 ground penetration: "
              f"{int((~np.isfinite(pen)).sum())} cells have no gen2 returns and are NaN, "
              f"not zero", flush=True)

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

    # Written HERE, after the grid products, not where they are computed. penetration is
    # built on this run's grid (bounds/res -> nx, ny) and the workflow graph says so, but
    # saving it before corrections.json/z_after.npy made it one second OLDER than its own
    # declared inputs -- so `lidar-diff-workflow` called it STALE the moment it was made,
    # permanently, at every site. The dependency is real; only the write order was wrong.
    # (Its VALUES do not come from z_after: make_penetration.py loads that array solely to
    # assert the shape, and ground_penetration reads the gen2 cloud.) The computation stays
    # above so that `pen` is computed while the gen2 cloud context is still open.
    if not skip_penetration:
        np.save(f"{outdir}/penetration.npy", pen)

    # The figures read the products just saved above, via the library, so either can
    # be rebuilt later without re-running any point-cloud work.
    fa = figures.dod_lod_figure(str(outdir), figdir, name)
    fb = figures.change_figure(str(outdir), figdir, name)
    print(f"[{name}] wrote {fa}  and  {fb}", flush=True)
    return dict(name=name, sigma=r["stable_sigma"], med_lod=float(np.nanmedian(lod)),
                n_regions=len(regions))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", help="run only these site names")
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--no-penetration", action="store_true",
                    help="skip the gen2 ground-penetration layer entirely. It is flagged by "
                         "analysis/ridgelines/AUDIT_findings.md as a gen2-derived variable "
                         "that should not drive gen1-internal conclusions, and "
                         "canopy_cover_pfs is the cover measure.")
    a = ap.parse_args()
    names = a.only if a.only else list(SITES)
    summary = []
    for nm in names:
        try:
            summary.append(run_site(nm, figdir=a.figdir,
                                    skip_penetration=a.no_penetration))
        except Exception as exc:
            import traceback; traceback.print_exc()
            print(f"[{nm}] FAILED: {exc}", flush=True)
    print("\n=== SUMMARY ===")
    for s in summary:
        print(f"  {s['name']:12s} sigma={s['sigma']:.3f} m  medLoD={s['med_lod']:.3f} m  "
              f"regions={s['n_regions']}")


if __name__ == "__main__":
    main()
