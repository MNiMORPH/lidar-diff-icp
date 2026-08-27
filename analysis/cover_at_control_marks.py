#!/usr/bin/env python3
"""Canopy cover at EVERY control mark, both epochs, from the gen2 3DEP EPT store.

WHY GEN2 FOR BOTH EPOCHS (Andy's call, 2026-08-27): the forest mostly did not change in the
13 years between acquisitions, so ONE source, ONE estimator and ONE scale beats epoch-matched
processing. It also matches `canopy_cover_pfs.npy`, which is already a gen2 quantity, and --
via the EPT -- it reaches every mark in the state rather than the 88/6 that happen to sit
under clouds on local disk.

WHAT THE COVARIATE IS: 2021 LEAF-ON structure. That is the right thing for ranking stands by
density. It is NOT an absolute optical-penetration quantity for gen1's 2008 leaf-off beams, so
a coefficient fitted against gen1 residuals absorbs the leaf-state conversion. Stands harvested
between the epochs are a real exception and are NOT screened here -- flag them downstream
(`refcells.py` carries a one-sided clear-cut test for exactly this).

COMPARABILITY, MEASURED NOT ASSUMED. Settings are Elba's: res=5.0, voxel_z=1.0, min_height=2.0,
DTM-HAG (not Delaunay). The one substitution is the DTM SOURCE -- Elba used the pipeline's
`z_after`; a scattered mark has no such grid, so the DTM is built from the box's own class-2
returns. Checked against `data/derived/elba_fulldensity/canopy_cover_pfs.npy` on 40 cells
stratified across the cover range:

    mine r=2.5 (single cell) vs canopy_cover_pfs: bias +0.020  sd 0.037  r=+0.989
    mine r=7.5               vs canopy_cover_pfs: bias -0.002  sd 0.090  r=+0.928

RADIUS LADDER, NOT A CHOSEN BOX. Read cost is latency-bound and FLAT in box size (0.9-1.4 s
from 5 m to 100 m), so one generous read serves every radius. Cover is reported at 2.5, 5,
7.5, 10, 15 and 25 m. **r=7.5 m is primary**: it matches the tie's own report radius
(1.5 x res), so cover describes the ground the residual actually integrated.
`cover_ladder_spread` is a per-mark QUALITY column, not a filter -- stable across radii means a
homogeneous stand and a meaningful value; swinging means a canopy edge. (This is the opposite
of the radius-spread SCREEN on ties, which was measured useless and must not be revived.)

PARAMETERS I CHOSE, flagged: `k=10` (a 105 m box, 21x21 cells of 5 m) and the MARK-CENTRED
grid. Mark-centred is deliberate: snapping to an arbitrary global grid would make the smallest
radius depend on where a mark happens to fall inside its cell. `k` only needs to exceed the
largest radius plus the wrapper's 25 m halo.

Resumable: one row appended per mark, marks already in the output are skipped.

    PROJ_DATA=$CENV/share/proj GDAL_DATA=$CENV/share/gdal \
    $CENV/bin/python analysis/cover_at_control_marks.py --out data/derived/control_cover.csv
"""
import argparse, csv, json, os, re, sys, tempfile, time
import numpy as np, pdal, pandas as pd
from pyproj import Transformer
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from forest_metrics_pfs import canopy_cover_raster, _write_dtm

# The 2021 acquisition is served as FIVE EPT blocks, not one. A mark outside all five is
# reported as uncovered -- never silently served from a DIFFERENT acquisition (IA/WI projects
# and MN_WasecaCo_2010 / MN_BlueEarth_2011 all overlap this footprint), which would destroy
# the one-source premise this whole approach rests on.
EPT_FAMILY = re.compile(r"^MN_SEDriftless_\d_2021$")
EPT_URL = ("https://s3-us-west-2.amazonaws.com/usgs-lidar-public/{name}/ept.json")
BOUNDARIES = "data/3dep_boundaries.geojson"
RADII = (2.5, 5.0, 7.5, 10.0, 15.0, 25.0)


def _ept_index():
    """(STRtree, polygons, names) over the 2021 acquisition's five EPT blocks."""
    from shapely.geometry import shape
    from shapely.strtree import STRtree
    gj = json.load(open(BOUNDARIES))
    polys, names = [], []
    for f in gj["features"]:
        nm = f["properties"].get("name", "")
        if EPT_FAMILY.match(nm):
            polys.append(shape(f["geometry"])); names.append(nm)
    if not polys:
        raise SystemExit(f"no MN_SEDriftless_*_2021 blocks in {BOUNDARIES}")
    return STRtree(polys), polys, names


def _block_for(tree, polys, names, lon, lat):
    from shapely.geometry import Point
    pt = Point(lon, lat)
    for i in tree.query(pt):
        if polys[i].contains(pt):
            return names[i]
    return None
G1 = "src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv"
G2 = "src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv"
_TR = Transformer.from_crs(26915, 3857, always_xy=True)
_TR4326 = Transformer.from_crs(26915, 4326, always_xy=True)


def cover_at(E, N, *, res, k, voxel_z, min_height, td, ept_url):
    half = (k + 0.5) * res
    b = (E - half, N - half, E + half, N + half)
    x0, y0 = _TR.transform(b[0], b[1]); x1, y1 = _TR.transform(b[2], b[3])
    laz = os.path.join(td, "box.laz")
    p = pdal.Pipeline(json.dumps({"pipeline": [
        {"type": "readers.ept", "filename": ept_url, "bounds": f"([{x0},{x1}],[{y0},{y1}])"},
        {"type": "filters.reprojection", "out_srs": "EPSG:26915"},
        {"type": "writers.las", "filename": laz, "compression": "laszip", "a_srs": "EPSG:26915"}]}))
    n = p.execute()
    if n == 0:
        return {"status": "no points in the EPT box"}
    a = p.arrays[0]; g = a[a["Classification"] == 2]
    nx = ny = 2 * k + 1
    dtm = np.full((ny, nx), np.nan)
    if len(g):
        ix = ((g["X"] - b[0]) // res).astype(int); iy = ((g["Y"] - b[1]) // res).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        for i, j, z in zip(ix[ok], iy[ok], g["Z"][ok]):
            dtm[j, i] = z if not np.isfinite(dtm[j, i]) else min(dtm[j, i], z)
    m = ~np.isfinite(dtm)
    if m.all():
        os.unlink(laz); return {"status": "no class-2 ground in the box", "n_points": int(n)}
    if m.any():
        dtm = dtm[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
    tif = os.path.join(td, "dtm.tif")
    _write_dtm(dtm, b[0], b[1], res, "EPSG:26915", tif)
    cc, pai = canopy_cover_raster(laz, b, res, dtm_path=tif, tile_m=400.0, halo_m=25.0,
                                  voxel_z=voxel_z, min_height=min_height)
    for f in (laz, tif):
        if os.path.exists(f):
            os.unlink(f)
    yy, xx = np.mgrid[0:ny, 0:nx]
    d = np.hypot(b[0] + (xx + 0.5) * res - E, b[1] + (yy + 0.5) * res - N)
    out = {"status": "ok", "n_points": int(n), "n_ground": int(len(g))}
    for r in RADII:
        s = (d <= r) & np.isfinite(cc)
        out[f"cover_r{r:g}"] = round(float(np.nanmean(cc[s])), 6) if s.any() else ""
        out[f"pai_r{r:g}"] = round(float(np.nanmean(pai[s])), 6) if s.any() else ""
        out[f"ncell_r{r:g}"] = int(s.sum())
    v = [out[f"cover_r{r:g}"] for r in RADII if out[f"cover_r{r:g}"] != ""]
    out["cover_ladder_spread"] = round(max(v) - min(v), 6) if len(v) > 1 else ""
    return out


def marks():
    g1 = pd.read_csv(G1).drop_duplicates(subset=["easting", "northing", "elevation"])
    for _, r in g1.iterrows():
        yield ("gen1_2008_control", r.point_id, float(r.easting), float(r.northing),
               str(r.point_type))
    g2 = pd.read_csv(G2)
    for _, r in g2.iterrows():
        yield ("gen2_2021_control", r.point_id, float(r.easting), float(r.northing),
               str(r.point_type))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--res", type=float, default=5.0, help="Elba's grid resolution")
    ap.add_argument("--k", type=int, default=10, help="half-width in cells; box = (2k+1)*res")
    ap.add_argument("--voxel-z", type=float, default=1.0)
    ap.add_argument("--min-height", type=float, default=2.0)
    ap.add_argument("--sleep", type=float, default=0.0, help="pause between EPT reads, s")
    a = ap.parse_args()

    cols = (["set", "point_id", "easting", "northing", "point_type", "ept_block", "status",
             "n_points", "n_ground"]
            + [f"{p}_r{r:g}" for r in RADII for p in ("cover", "pai", "ncell")]
            + ["cover_ladder_spread"])
    done = set()
    if os.path.exists(a.out):
        d = pd.read_csv(a.out)
        done = set(zip(d["set"], d.point_id))
        print(f"resuming: {len(done)} marks already done", flush=True)
    new = not os.path.exists(a.out)
    todo = [m for m in marks() if (m[0], m[1]) not in done]
    print(f"{len(todo)} marks to do   res={a.res} voxel_z={a.voxel_z} min_height={a.min_height} "
          f"k={a.k} (box {(2*a.k+1)*a.res:.0f} m)", flush=True)
    tree, polys, names = _ept_index()
    print(f"EPT blocks: {', '.join(sorted(names))}", flush=True)
    td = tempfile.mkdtemp(prefix="cover_marks_")
    t0 = time.time()
    with open(a.out, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        if new:
            w.writeheader()
        for i, (setname, pid, E, N, ptype) in enumerate(todo, 1):
            lon, lat = _TR4326.transform(E, N)
            blk = _block_for(tree, polys, names, lon, lat)
            if blk is None:
                o = {"status": "outside the MN_SEDriftless_*_2021 acquisition"}
            else:
                try:
                    o = cover_at(E, N, res=a.res, k=a.k, voxel_z=a.voxel_z,
                                 min_height=a.min_height, td=td,
                                 ept_url=EPT_URL.format(name=blk))
                except Exception as ex:
                    o = {"status": f"FAIL {type(ex).__name__}: {ex}"[:120]}
            w.writerow({"set": setname, "point_id": pid, "easting": E, "northing": N,
                        "point_type": ptype, "ept_block": blk or "", **o})
            fh.flush()
            if i % 25 == 0 or i == len(todo):
                el = time.time() - t0
                print(f"  {i}/{len(todo)}  {el/i:.2f} s/mark  eta {(len(todo)-i)*el/i/60:.0f} min",
                      flush=True)
            if a.sleep:
                time.sleep(a.sleep)
    print(f"wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()
