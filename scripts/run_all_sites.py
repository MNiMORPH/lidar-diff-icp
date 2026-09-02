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
from lidar_diff_icp import figures
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


def run_site(name, figdir=FIGDIR, *, skip_penetration=False,
             leafon_lod_factor=None):
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

    # gen2 ground penetration, and the leaf-on/forest-slope flag derived from it.
    #
    # BOTH ARE OPTIONAL, and the LoD widening they used to drive is now OFF unless asked
    # for. It was applied at every site with factor 2.0, a number with no derivation
    # anywhere in this repo, on 5.8-40% of each tile -- so it silently set what counts as
    # detected change over a large fraction of every product. The mechanism it rested on
    # (gen2 ground STARVED under leaf-on canopy) is also contradicted by later work in this
    # repo: gen2 carries ~12x MORE ground than gen1, and the offset tracks the 2008 canopy
    # rather than 2021 cover. See canopy.py's status section.
    #
    # The flag is still computed and SAVED when penetration is available, because it is a
    # reasonable thing to look at. It just no longer changes the LoD behind the reader's
    # back.
    if skip_penetration:
        print(f"[{name}] penetration: SKIPPED as stated (--no-penetration); no leaf-on "
              f"flag, and the LoD is the heteroscedastic one alone", flush=True)
    else:
        from lidar_diff_icp.canopy import (ground_penetration, leafon_slope_flag,
                                           inflate_lod)
        _Zf = Z21.copy(); _nm = ~np.isfinite(_Zf)
        if _nm.any():
            from scipy.ndimage import distance_transform_edt as _edt
            _Zf = _Zf[tuple(_edt(_nm, return_distances=False, return_indices=True))]
        _sl = np.degrees(np.arctan(np.hypot(*np.gradient(_Zf, res)[::-1])))
        # Neither threshold is derived: 12 deg is the repo-wide gentle-ground cut and 0.25
        # mirrors the forest strata cut. Stated in the run's own output, not left implicit.
        PEN_MAX, SLOPE_MIN = 0.25, 12.0
        pen = ground_penetration(after, r["bounds"], res, nx, ny)
        leafon = leafon_slope_flag(pen, _sl, min_penetration=PEN_MAX, min_slope=SLOPE_MIN)
        msg = (f"[{name}] leaf-on/forest-slope flag: {int(leafon.sum())} cells "
               f"({100*leafon.mean():.0f}%), flagged where penetration < {PEN_MAX:g} AND "
               f"slope > {SLOPE_MIN:g} deg (conventions, not derivations); "
               f"{int((~np.isfinite(pen)).sum())} cells have no gen2 returns and are NaN, "
               f"not flagged")
        if leafon_lod_factor is None:
            msg += " -- LoD NOT widened (pass --leafon-lod-factor to widen it)"
        else:
            lod = inflate_lod(lod, leafon, factor=leafon_lod_factor)
            msg += f" -- LoD x{leafon_lod_factor:g} there, AS REQUESTED on the command line"
        print(msg, flush=True)

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
    # above because --leafon-lod-factor must inflate `lod` before `lod` is persisted.
    if not skip_penetration:
        np.save(f"{outdir}/leafon_flag.npy", leafon)
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
                         "canopy_cover_pfs is the cover measure. Without it there is no "
                         "leaf-on flag.")
    ap.add_argument("--leafon-lod-factor", type=float, default=None,
                    help="multiply the LoD by this on leaf-on/forest-slope cells. NOT "
                         "applied unless given. It used to default to 2.0, which has no "
                         "derivation in this repo and silently set the detection bar on "
                         "5.8-40%% of every tile.")
    a = ap.parse_args()
    names = a.only if a.only else list(SITES)
    summary = []
    for nm in names:
        try:
            summary.append(run_site(nm, figdir=a.figdir,
                                    skip_penetration=a.no_penetration,
                                    leafon_lod_factor=a.leafon_lod_factor))
        except Exception as exc:
            import traceback; traceback.print_exc()
            print(f"[{nm}] FAILED: {exc}", flush=True)
    print("\n=== SUMMARY ===")
    for s in summary:
        print(f"  {s['name']:12s} sigma={s['sigma']:.3f} m  medLoD={s['med_lod']:.3f} m  "
              f"regions={s['n_regions']}")


if __name__ == "__main__":
    main()
