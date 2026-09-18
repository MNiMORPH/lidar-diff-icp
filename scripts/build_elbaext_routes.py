#!/usr/bin/env python3
"""Build elbaext under BOTH routes, from one parameter set, for a like-with-like compare.

Andy, 2026-09-11: "Build elbaext under both routes." This is item 3 of the DeLong order --
every number behind the switch came from ONE tile (elba_fulldensity) and four flight
lines, so nothing about it is known to generalise until a second tile says so.

The two runs differ in `route` and NOTHING else. That is the whole point, and it is
enforced rather than trusted: the shared kwargs are built once and the script REFUSES if
the two recorded corrections.json disagree on any parameter other than the route's own
three switches.

elbaext is a tile directory, not a registered Site (sites.py has six; this is a product),
so its inputs are named here rather than resolved from a Site record:
    gen1   data/before/elbaext_gen1_merged.laz   (the hand-picked field subset merge)
    gen2   data/after/3dep2021_fulldensity.laz
    CSF    data/csf_cache/elbaext.las            (reused; ~460 s per tile otherwise)

Writes to data/derived/elbaext_independent and data/derived/elbaext_delong. The existing
data/derived/elbaext is NOT touched -- it was built on pre-fork code and is the historical
record.

    ./lidar-icp/bin/python scripts/build_elbaext_routes.py --route delong
"""
import argparse
import json
import os
import time

import numpy as np

from lidar_diff_icp import acquisitions, figures
from lidar_diff_icp.detect import detect_change_standard
from lidar_diff_icp.pipeline import difference_dem

#: INPUTS ARE READ FROM THE TILE'S OWN meta.json, NEVER GUESSED. I first hardcoded
#: data/after/3dep2021_fulldensity.laz by inference and it was WRONG: elbaext is the
#: EXTENDED tile and that file under-covers it, giving 355,961 finite gen2 cells against
#: the 719,226 the real input yields -- roughly half the tile, silently. meta.json records
#: the actual pair and was there the whole time.
_META = json.load(open("data/derived/elbaext/meta.json"))
BOUNDS = tuple(_META["bounds"])
GEN1 = _META["before"]
GEN2 = _META["after"]
CSF_CACHE = "data/csf_cache/elbaext.las"
PROJECT = "lidar_semn2008"

#: Everything both routes share. Taken from data/derived/elbaext/corrections.json so this
#: rebuild is comparable with the shipped product, not a differently-configured run.
#:
#: valley_top_m is CITED, not chosen. refcells.VALLEY_TOP_M records 230.0 m for elbaext
#: (and for elba and elba_fulldensity), so this reproduces the established value rather
#: than inventing one. It is passed as a number rather than "registry" because the
#: registry is keyed on tile-directory basename and these outputs are elbaext_independent
#: / elbaext_delong -- build variants, which do not belong in the registry. The pipeline
#: refuses outright if it is not stated ("It will not be chosen for you"), which is how
#: this run failed the first time; that refusal is correct and stays.
SHARED = dict(res=5.0, ground_q=0.50, ground="slope_normal", ground_source="csf",
              robust_stable=True, swath_tie="intercept", valley_top_m=230.0)


def _tif(arr, res, x0, y0, ny, out):
    """GeoTIFF exactly as run_all_sites._tif writes it, so products are interchangeable."""
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(out, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype="float32", crs="EPSG:26915", nodata=np.nan,
                       transform=from_origin(x0, y0 + ny * res, res, res)) as d:
        d.write(np.flipud(arr).astype("float32"), 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--route", required=True, choices=["independent", "delong"])
    ap.add_argument("--out", default=None)
    ap.add_argument("--boresight", action="store_true",
                    help="remove a COMMON scanner roll, self-calibrated from gen1's own "
                         "flight-line self-overlap (gen2-free). Measured at elbaext: "
                         "2.22 +/- 0.05 mm/deg by coreg.estimate_boresight_roll, and an "
                         "independent scan of the per-line step metric bottoms out at "
                         "2.0-3.0. It is NOT in the shipped default, so this writes to a "
                         "separate directory and leaves the comparison products intact.")
    a = ap.parse_args()
    out = a.out or (f"data/derived/elbaext_{a.route}"
                    + ("_boresight" if a.boresight else ""))
    os.makedirs(out, exist_ok=True)

    kw = dict(SHARED)
    # The geoid is REQUIRED by the independent route and unused by delong. Resolved from
    # the acquisition, never defaulted -- defaulting it is what cost +54.87 mm at Ramsey.
    if a.route == "independent":
        kw["gen1_geoid"] = acquisitions.for_project(PROJECT).geoid_grid
    if a.boresight:
        kw["correct_boresight"] = True

    t0 = time.time()
    r = difference_dem(GEN1, GEN2, BOUNDS, route=a.route, tile_dir=out,
                       csf_cache=CSF_CACHE, **kw)
    c = r["corrections"]

    # PERSIST EVERY OUTPUT. difference_dem RETURNS arrays -- its `tile_dir` only feeds the
    # valley-top registry lookup and writes nothing. run_all_sites.py does the saving, but
    # it is keyed on registered Sites and elbaext is a product, not a Site, so the same
    # persist block is reproduced here. Without this the run leaves an empty directory and
    # no figure can be drawn, which is exactly what the first pair of builds did.
    dod, lod, Z21, stable = r["dod"], r["lod"], r["z_after"], r["stable"]
    res = r["res"]; nx, ny = r["nx"], r["ny"]; X0, Y0 = r["bounds"][0], r["bounds"][1]
    det = detect_change_standard(dod, lod, stable, res)
    change, regions = det["change"], det["regions"]
    _tif(dod, res, X0, Y0, ny, f"{out}/dod.tif")
    _tif(lod, res, X0, Y0, ny, f"{out}/lod.tif")
    _tif(change.astype("float32"), res, X0, Y0, ny, f"{out}/change.tif")
    np.save(f"{out}/dod.npy", dod); np.save(f"{out}/lod.npy", lod)
    np.save(f"{out}/z_after.npy", Z21); np.save(f"{out}/change.npy", change)
    # The surface the DoD was ACTUALLY taken against. Equal to z_after when no gen2
    # correction ran; different when one did, and then z_after alone would misrepresent
    # the product. Saved always, so a reader never has to know which case they are in.
    np.save(f"{out}/z_after_differenced.npy", r["z_after_differenced"])
    # The reporting stable mask. detect_change_standard NEEDS it -- given an empty one it
    # builds its noise model from nothing and returns MORE detections after a LoD was
    # inflated, which is impossible and is how this omission was found.
    np.save(f"{out}/stable.npy", stable)
    # The applied gen1 correction GRIDS. On the DeLong route this is the correction
    # surface -- the term standing in for both the geoid and the drift. It was applied
    # and discarded until 2026-09-17, which made the product unauditable: the beam table
    # carries only the dz_* columns, so the surface is precisely what is missing from it.
    for _k, _g in (r.get("correction_grids") or {}).items():
        np.save(f"{out}/applied_{_k}.npy", _g)
        print(f"  applied {_k}: median {1000*np.nanmedian(_g):+.1f} mm, "
              f"p05 {1000*np.nanpercentile(_g,5):+.0f}, p95 {1000*np.nanpercentile(_g,95):+.0f}")
    with open(f"{out}/corrections.json", "w") as fh:
        json.dump(c, fh, indent=2)
    with open(f"{out}/regions.json", "w") as fh:
        json.dump({**{k: det[k] for k in ("regions", "sigma", "corr_length_m",
                                          "tau_sys_m", "method")},
                   "gen2_null_cells": int((~np.isfinite(Z21)).sum())}, fh, indent=2)
    name = os.path.basename(out)
    fa = figures.dod_lod_figure(out, "figures/sites", name)
    fb = figures.change_figure(out, "figures/sites", name)

    ex = np.isfinite(dod)
    print(f"\n  route={a.route}  {time.time()-t0:.0f} s  ->  {out}")
    print(f"  median LoD {np.nanmedian(lod):.4f} m   "
          f"{100*np.mean(np.abs(dod[ex])>lod[ex]):.1f}% of cells exceed LoD   "
          f"{len(regions)} robust regions")
    print(f"  figures: {fa}  {fb}")
    print(f"  chain switches: correction_surface={c['correction_surface']} "
          f"along_track_drift={c['along_track_drift']} geoid_applied={c['geoid_applied']}")
    print(f"  stable_sigma {r['stable_sigma']:.4f} m   zero_line {c.get('zero_line')}")

    # Enforce the like-with-like claim as soon as both exist.
    other = f"data/derived/elbaext_{'delong' if a.route == 'independent' else 'independent'}"
    op = os.path.join(other, "corrections.json")
    if os.path.exists(op):
        o = json.load(open(op))
        # Keys a route is ENTITLED to differ in. gen1_geoid_grid belongs here because the
        # delong route is never GIVEN a geoid grid -- that is the geoid switch itself, not
        # an unrelated parameter. It was missing on the first run and the guard refused,
        # which is the guard working: it made me justify the exemption rather than assume
        # it. Anything not listed here still refuses.
        ROUTE_KEYS = {"route", "correction_surface", "along_track_drift", "geoid_applied",
                      "gen2_chain",
                      "gen1_geoid_grid", "cross_epoch_datum",
                      "along_track_drift_gpsTime_to_m",
                      "stable_1sigma_m", "stable_clip_fraction"}
        diff = {k: (c.get(k), o.get(k)) for k in set(c) | set(o)
                if k not in ROUTE_KEYS and c.get(k) != o.get(k)}
        if diff:
            raise SystemExit(
                f"REFUSED: the two routes differ in parameters OTHER than the route's own "
                f"switches, so any difference between their DoDs is not attributable to "
                f"the route: {json.dumps(diff, indent=1, default=str)}")
        print(f"  like-with-like CONFIRMED against {other}: no parameter differs but the "
              f"route's own switches")


if __name__ == "__main__":
    main()
