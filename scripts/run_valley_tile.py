#!/usr/bin/env python3
"""Phase B: one Whitewater-valley tile, end to end, the way the Elba product was built.

Andy, 2026-10-20: "Run all of your analyses through the entire valley. This might require
also downloading lidar tiles. Do it." and "I want you to be as autonomous as possible."
Decisions in force (analysis/PLAN_VALLEY_FIELD_SITES.md): D1 the valley floor from the
official 3DEP DEM by Clubb et al. 2017 (the 1500 m-corridor run at 10 m); D2 class-2 gen2
for the DoD; D6 public-land tiles only (29); D7 time the first tile before any loop.

PER TILE, resumable (every stage is skipped when its product exists):
  1. gen1   data/before/<tile>.laz from MnGeo (scripts/fetch_tile.py; county resolved from
            the three county listings, data/derived/valley/public_tiles.json).
  2. bounds from the gen1 header, ceil/floor to the 5 m grid, exactly as run_all_sites does.
  3. gen2   the FULL-density 3DEP cloud for the bounds (scripts/fetch_3dep_curl.py --auto,
            every node at every depth), written to data/after/valley/<tile>_3dep_full.laz;
            then a class-2-only copy <tile>_3dep_class2.laz (laspy, chunked) for the DoD.
            Full density is kept because the lowveg grid needs ALL near-ground returns
            (scripts/make_lowveg.py); that is my reason for keeping more than D2 asked.
            The per-tile EPT node cache is deleted after the LAS is written.
  4. completeness record (area-weighted EPT point count inside the bounds vs the file),
            the same arithmetic as analysis/ept_coverage_check.py, recorded never judged.
  5. base   pipeline.difference_dem with the Elba SHARED kwargs, route delong, swath_frame
            star, correction_surface False (= data/derived/elbaext_delong_star_nosurf),
            CSF cached at data/csf_cache/<tile>.las; detect_change_standard; every
            product saved as run_all_sites saves it. VALLEY TOP per tile: the 95th
            percentile of the official DEM over the Clubb corridor floodplain cells inside
            the tile (MINE, D1 detail; Elba's registry value is 230.0 m, and the rule is
            printed beside it for the check). Falls back to the tile's DEM 20th percentile
            if the tile holds no Clubb floodplain (printed when it happens).
  6. steps  slope, ridge_mask, convexity (with that valley top), curvature, gen1_angles,
            beam_table, then gen1 span + penetration failure + the vegetation LoD
            (scripts/figures_with_vegetation_lod.py, Elba's 450 mm floor and 2 deg),
            and the gen2 lowveg grid (scripts/make_lowveg.py --gen2 <full cloud>).
            pfs_cover is skipped (plan section 3: needs an untwine COPC and ~14 GB scratch).
  7. clubb  floodplain_mask_clubb.npy: the D1 floodplain on the tile grid (nearest), beside
            the pipeline's own convexity mask, never over it.
Stage times go to data/derived/<tile>/valley_run.json.

    ./lidar-icp/bin/python scripts/run_valley_tile.py --tile 4342-28-64
    ./lidar-icp/bin/python scripts/run_valley_tile.py --tile 4342-28-64 --stop-after base
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request

import laspy
import numpy as np
import rasterio
from pyproj import Transformer

sys.path.insert(0, "scripts")
from run_all_sites import header_bounds, _tif  # noqa: E402
from lidar_diff_icp import threedep, completeness, figures  # noqa: E402
from lidar_diff_icp.detect import detect_change_standard  # noqa: E402
from lidar_diff_icp.pipeline import difference_dem  # noqa: E402

PY = "env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python"
RES = 5.0
SHARED = dict(res=RES, ground_q=0.50, ground="slope_normal", ground_source="csf",
              robust_stable=True, swath_tie="intercept")      # build_elbaext_routes.SHARED minus valley_top_m
ROUTE = dict(route="delong", swath_frame="star", correction_surface=False)   # = elbaext_delong_star_nosurf
CLUBB_FP = "data/derived/valley/clubb/corr1500_10m/corr1500_10m_valley.bil"
CLUBB_DEM = "data/derived/valley/clubb/valley10m/valley10m.bil"
PUBLIC = json.load(open("data/derived/valley/public_tiles.json"))

ap = argparse.ArgumentParser()
ap.add_argument("--tile", required=True, help="AAAA-BB-CC")
ap.add_argument("--stop-after", default=None, choices=["gen1", "gen2", "base", "steps"])
ap.add_argument("--workers", type=int, default=8, help="3DEP download workers (16 is the script default; 8 is gentler on the shared laptop)")
A = ap.parse_args()
T = A.tile; D = f"data/derived/{T}"; os.makedirs(D, exist_ok=True); os.makedirs("data/after/valley", exist_ok=True)
GEN1 = f"data/before/{T}.laz"; FULL = f"data/after/valley/{T}_3dep_full.laz"; CLS2 = f"data/after/valley/{T}_3dep_class2.laz"; CSF = f"data/csf_cache/{T}.las"
LOG = f"{D}/valley_run.json"; log = json.load(open(LOG)) if os.path.exists(LOG) else {"tile": T, "stages": {}}
def stamp(stage, t0, **kw):
    log["stages"][stage] = dict(seconds=round(time.time() - t0, 1), when=time.strftime("%Y-%m-%dT%H:%M:%S"), **kw); json.dump(log, open(LOG, "w"), indent=1)
    print(f"[{T}] {stage}: {time.time() - t0:.0f} s {kw if kw else ''}", flush=True)
def sh(cmd):
    print(f"[{T}] $ {cmd}", flush=True); r = subprocess.run(cmd, shell=True)
    if r.returncode != 0: sys.exit(f"[{T}] FAILED ({r.returncode}): {cmd}")

# ---- 1. gen1
if not os.path.exists(GEN1):
    t0 = time.time(); county = PUBLIC[T]["counties"][0]
    sh(f"{PY} scripts/fetch_tile.py --tile {T} --county {county} --out data/before"); stamp("gen1_fetch", t0, county=county, bytes=os.path.getsize(GEN1))
else: print(f"[{T}] gen1 present: {GEN1}")
B = header_bounds(GEN1, RES); print(f"[{T}] bounds from the gen1 header: {tuple(round(b, 1) for b in B)}")
if A.stop_after == "gen1": sys.exit(0)

# ---- 3. gen2 full, then class 2
if not os.path.exists(FULL):
    t0 = time.time(); cache = f"data/after/valley/_ept_{T}"
    sh(f"PROJ_DATA=/usr/share/proj GDAL_DATA=/usr/share/gdal ./lidar-icp/bin/python scripts/fetch_3dep_curl.py --auto --bounds {B[0]} {B[1]} {B[2]} {B[3]} --workers {A.workers} --out {FULL} --ept-cache data/ept_index_cache.json")
    stamp("gen2_fetch", t0, bytes=os.path.getsize(FULL))
    for c in glob.glob("data/after/valley/_ept_*") + glob.glob("data/after/_ept_tiles_*"): shutil.rmtree(c, ignore_errors=True)
else: print(f"[{T}] gen2 full present: {FULL}")
if not os.path.exists(CLS2):
    t0 = time.time(); n_in = n_out = 0
    with laspy.open(FULL) as f:
        w = laspy.open(CLS2, mode="w", header=f.header)
        for pts in f.chunk_iterator(5_000_000):
            keep = np.asarray(pts.classification) == 2; n_in += len(pts); n_out += int(keep.sum())
            if keep.any(): w.write_points(pts[keep])
        w.close()
    stamp("gen2_class2", t0, points_full=n_in, points_class2=n_out)
else: print(f"[{T}] gen2 class-2 present: {CLS2}")

# ---- 4. completeness (gen2), the ept_coverage_check arithmetic
if "completeness" not in log["stages"]:
    t0 = time.time()
    def available(bounds):
        to_ll = Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True); get = lambda u: json.load(urllib.request.urlopen(u, timeout=60))
        pts = [to_ll.transform(cx, cy) for cx in bounds[0::2] for cy in bounds[1::2]]; lons = [p[0] for p in pts]; lats = [p[1] for p in pts]
        ref = threedep.resolve_reference((min(lons) + max(lons)) / 2, (min(lats) + max(lats)) / 2, (min(lons), min(lats), max(lons), max(lats)), cache="data/ept_index_cache.json")
        base = ref["url"].rsplit("/ept.json", 1)[0]; ept = get(f"{base}/ept.json"); b = ept["bounds"]; side = b[3] - b[0]
        t = Transformer.from_crs("EPSG:26915", "EPSG:3857", always_xy=True); xs, ys = zip(*[t.transform(cx, cy) for cx in bounds[0::2] for cy in bounds[1::2]]); bx0, bx1, by0, by1 = min(xs), max(xs), min(ys), max(ys)
        def node(d, x, y): s = side / (2 ** d); return (b[0] + x * s, b[1] + y * s, b[0] + (x + 1) * s, b[1] + (y + 1) * s, s)
        def ov(d, x, y): n = node(d, x, y); return not (n[2] < bx0 or n[0] > bx1 or n[3] < by0 or n[1] > by1)
        def frac_in(d, x, y): n = node(d, x, y); return max(0.0, min(n[2], bx1) - max(n[0], bx0)) * max(0.0, min(n[3], by1) - max(n[1], by0)) / (n[4] * n[4])
        frontier = [get(f"{base}/ept-hierarchy/0-0-0-0.json")]; seen = set(); inb = 0.0; avail = 0
        while frontier:
            H = frontier.pop()
            for key, cnt in H.items():
                d, x, y, z = map(int, key.split("-"))
                if not ov(d, x, y): continue
                if cnt == -1:
                    if key not in seen: seen.add(key); frontier.append(get(f"{base}/ept-hierarchy/{key}.json"))
                else: avail += cnt; inb += cnt * frac_in(d, x, y)
        return ref["name"], avail, inb
    proj, avail, inb = available(B)
    with laspy.open(FULL) as f: have = f.header.point_count
    completeness.write(D, epoch="gen2", cloud=FULL, points_in_file=have, points_available=inb, measured_by=f"scripts/run_valley_tile.py against {proj} (ept_coverage_check arithmetic)", note="area-weighted share of the EPT node set inside the bbox; nodes clip the box, so a raw node-count sum overstates a complete fetch")
    stamp("completeness", t0, project=proj, points_in_file=int(have), points_available=round(inb), ratio=round(have / max(inb, 1), 4))
if A.stop_after == "gen2": sys.exit(0)

# ---- 5. base: valley top from the Clubb corridor floodplain, then difference_dem
if not os.path.exists(f"{D}/dod.npy"):
    t0 = time.time()
    with rasterio.open(CLUBB_FP) as r: fp = r.read(1); tr = r.transform
    with rasterio.open(CLUBB_DEM) as r: dem = r.read(1); dem = np.where(dem == r.nodata, np.nan, dem)
    ny_, nx_ = fp.shape; jy, jx = np.mgrid[0:ny_, 0:nx_]; X, Y = tr * (jx + 0.5, jy + 0.5)
    inside = (X >= B[0]) & (X < B[2]) & (Y >= B[1]) & (Y < B[3])
    # only the Clubb floodplain WITHIN 500 m OF THE WHITEWATER CENTRELINES counts: the
    # corridor run still holds flat pockets on the plateau inside a bluff tile, and p95 over
    # them put 4342-29-63's "valley top" at 351 m (checked 2026-10-20). 500 m is mine.
    import shapely
    from shapely.geometry import shape
    riv = shapely.union_all([shape(f["geometry"]) for f in json.load(open("data/derived/valley/public_land/rivers_whitewater.geojson"))["features"]]).buffer(500.0)
    near = inside.copy(); near[inside] = shapely.contains_xy(riv, X[inside], Y[inside])
    fpc = near & (fp == 1) & np.isfinite(dem)
    if fpc.sum() >= 50:
        vt = float(np.nanpercentile(dem[fpc], 95)); vsrc = f"p95 of the official DEM over {int(fpc.sum())} Clubb floodplain cells within 500 m of the Whitewater centrelines in the tile (MINE)"
    else:
        vt = float(np.nanpercentile(dem[inside & np.isfinite(dem)], 5)) + 15.0; vsrc = f"only {int(fpc.sum())} Clubb floodplain cells near the river in the tile: p05 of the tile's DEM + 15 m (fallback, MINE)"
    print(f"[{T}] valley top {vt:.1f} m: {vsrc}; Elba's registry value is 230.0 m", flush=True)
    r = difference_dem(GEN1, CLS2, B, tile_dir=D, csf_cache=CSF, valley_top_m=vt, **SHARED, **ROUTE)
    c = r["corrections"]; dod, lod, Z21, stable = r["dod"], r["lod"], r["z_after"], r["stable"]; nx, ny = r["nx"], r["ny"]; X0, Y0 = r["bounds"][0], r["bounds"][1]
    det = detect_change_standard(dod, lod, stable, RES); change, regions = det["change"], det["regions"]
    _tif(dod, RES, X0, Y0, ny, f"{D}/dod.tif"); _tif(lod, RES, X0, Y0, ny, f"{D}/lod.tif"); _tif(change.astype("float32"), RES, X0, Y0, ny, f"{D}/change.tif")
    np.save(f"{D}/dod.npy", dod); np.save(f"{D}/lod.npy", lod); np.save(f"{D}/z_after.npy", Z21); np.save(f"{D}/change.npy", change); np.save(f"{D}/z_after_differenced.npy", r["z_after_differenced"]); np.save(f"{D}/stable.npy", stable)
    for _k, _g in (r.get("correction_grids") or {}).items(): np.save(f"{D}/applied_{_k}.npy", _g)
    c["valley_top_rule"] = vsrc; json.dump(c, open(f"{D}/corrections.json", "w"), indent=2)
    json.dump({**{k: det[k] for k in ("regions", "sigma", "corr_length_m", "tau_sys_m", "method")}, "gen2_null_cells": int((~np.isfinite(Z21)).sum())}, open(f"{D}/regions.json", "w"), indent=2)
    json.dump({"bounds": list(B), "before": GEN1, "after": CLS2, "after_full": FULL, "csf_cache": CSF, "valley_top_m": vt, "valley_top_rule": vsrc, "built_by": "scripts/run_valley_tile.py", "like": "data/derived/elbaext_delong_star_nosurf"}, open(f"{D}/meta.json", "w"), indent=2)
    ex = np.isfinite(dod)
    stamp("base", t0, valley_top_m=round(vt, 2), stable_sigma_m=round(float(r["stable_sigma"]), 4), median_lod_mm=round(float(np.nanmedian(lod)) * 1000, 1), pct_exceed_lod=round(100 * float(np.mean(np.abs(dod[ex]) > lod[ex])), 1), regions=len(regions), zero_line=c.get("zero_line"))
    try: figures.dod_lod_figure(D, "figures/sites", T); figures.change_figure(D, "figures/sites", T)
    except Exception as e: print(f"[{T}] figures skipped: {e}")
else: print(f"[{T}] base present")
if A.stop_after == "base": sys.exit(0)

# ---- 6. the derived steps, each skipped when its product exists
vt = json.load(open(f"{D}/meta.json"))["valley_top_m"]
steps = [("slope", "slope.npy", f"{PY} scripts/make_slope.py --tile {T}"),
         ("ridge_mask", "ridge_mask.npy", f"{PY} -m lidar_diff_icp.steps.trace_ridgelines {D} --out ridge_mask.npy --no-fig"),
         ("convexity", "floodplain_mask.npy", f"{PY} -m lidar_diff_icp.steps.convexity_dod_landcover --tile {T} --dod {D}/dod.npy --without cover --valley-top {vt}"),
         ("curvature", "curv_laplacian.npy", f"{PY} -m lidar_diff_icp.steps.curvature_diffusion --tile {T}"),
         ("gen1_angles", "gen1_csf_angles.npz", f"{PY} analysis/modules/vegetation_correction/gen1_save_angles_slope.py {T} {CSF}"),
         ("beam_table", "beam_offset_table.parquet", f"{PY} analysis/modules/vegetation_correction/beam_offset_table.py {D} {CSF}"),
         ("vegetation_lod", "penetration_failure.npy", f"{PY} scripts/figures_with_vegetation_lod.py --tile {D} --failure-floor-mm 450 --slope-max-deg 2"),
         ("lowveg", "lowveg_grid.npy", f"{PY} scripts/make_lowveg.py --tile {T} --gen2 {FULL}")]
for name, product, cmd in steps:
    if os.path.exists(f"{D}/{product}"): print(f"[{T}] {name} present"); continue
    t0 = time.time(); sh(cmd); stamp(name, t0)

# ---- 7. the D1 floodplain on the tile grid
if not os.path.exists(f"{D}/floodplain_mask_clubb.npy"):
    with rasterio.open(CLUBB_FP) as r: fp = r.read(1); tr = r.transform
    z = np.load(f"{D}/z_after.npy"); ny, nx = z.shape; jy, jx = np.mgrid[0:ny, 0:nx]; X = B[0] + (jx + 0.5) * RES; Y = B[1] + (jy + 0.5) * RES
    c, rr = ~tr * (X.ravel(), Y.ravel()); rr = np.clip(rr.astype(int), 0, fp.shape[0] - 1); c = np.clip(c.astype(int), 0, fp.shape[1] - 1)
    m = (fp[rr, c] == 1).reshape(ny, nx); np.save(f"{D}/floodplain_mask_clubb.npy", m)
    own = np.load(f"{D}/floodplain_mask.npy").astype(bool)
    print(f"[{T}] floodplain_mask_clubb: {m.sum() * 25 / 1e4:.0f} ha; pipeline convexity mask {own.sum() * 25 / 1e4:.0f} ha; both {np.sum(m & own) * 25 / 1e4:.0f} ha", flush=True)
    log["clubb_vs_pipeline_ha"] = [round(m.sum() * 25 / 1e4, 1), round(own.sum() * 25 / 1e4, 1), round(np.sum(m & own) * 25 / 1e4, 1)]; json.dump(log, open(LOG, "w"), indent=1)
print(f"[{T}] DONE; stages: " + ", ".join(f"{k} {v['seconds']:.0f}s" for k, v in log["stages"].items()), flush=True)
