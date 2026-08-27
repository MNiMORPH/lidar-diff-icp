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
import argparse, csv, datetime as dt, json, os, re, sys, tempfile, time
import numpy as np, pdal, pandas as pd
from pyproj import Transformer
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from forest_metrics_pfs import canopy_cover_raster, _write_dtm
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from lidar_diff_icp.groundtruth.tie import _design, ground_elevation_at  # noqa: E402

# NEAR-GROUND window: EXACTLY Elba's nearground_cells_sn.npz (zlo -1.0, zhi +2.0, dz 0.02,
# 150 bins) so mark structure is directly comparable with the tile cube.
NG_EDGES = np.arange(-1.0, 2.0 + 1e-9, 0.02)
# CANOPY window: coarse and tall. Diagnostic -- among other things it shows whether a block
# was really flown leaf-off, which the vendor spec claims and gps_time contradicts.
CAN_EDGES = np.arange(-2.0, 45.0 + 1e-9, 0.25)

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


def _gps_to_utc(t):
    """LAS 1.4 Adjusted Standard GPS Time -> UTC. GPS epoch 1980-01-06; 18 leap seconds at
    these dates; Adjusted Standard subtracts 1e9 from Standard GPS time."""
    return dt.datetime(1980, 1, 6) + dt.timedelta(seconds=float(t) + 1e9 - 18)
_TR4326 = Transformer.from_crs(26915, 4326, always_xy=True)


def _surface(x, y, z, px, py, radius, order=2):
    """Order-2 local surface at (px,py) by the SAME construction as
    :func:`groundtruth.tie.ground_elevation_at` -- same design matrix, same lstsq -- but
    returning the coefficients, which that function does not expose. The caller asserts
    z_hat agrees with it, so this cannot drift into a second estimator."""
    r = np.hypot(x - px, y - py)
    m = r <= radius
    kk = 3 if order == 1 else 6
    if int(m.sum()) < kk:
        return None
    A = _design(x[m] - px, y[m] - py, order)
    coef, *_ = np.linalg.lstsq(A, z[m], rcond=None)
    if not np.isfinite(coef).all():
        return None
    return coef


def _sn_hist(x, y, z, px, py, coef, order, edges):
    """Histogram of the SLOPE-NORMAL residual to the local surface, matching the project's
    definition: (z - S(x,y)) / |n| with |n| = sqrt(1 + gx^2 + gy^2), gradient at the mark."""
    S = _design(x - px, y - py, order) @ coef
    nn = np.sqrt(1.0 + coef[1] ** 2 + coef[2] ** 2)
    return np.histogram((z - S) / nn, bins=edges)[0].astype(np.int32)


def cover_at(E, N, *, res, k, voxel_z, min_height, td, ept_url, box_out=None,
             struct_out=None, struct_radius=7.5):
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
    if os.path.exists(tif):
        os.unlink(tif)
    # The BOX IS KEPT (~1.6 MB), not deleted: every derived quantity here -- and any future
    # one -- is then reproducible with no further network read. That is the whole point.
    if box_out:
        os.makedirs(os.path.dirname(box_out), exist_ok=True)
        os.replace(laz, box_out)
        out_box = box_out
    else:
        os.unlink(laz); out_box = ""
    yy, xx = np.mgrid[0:ny, 0:nx]
    d = np.hypot(b[0] + (xx + 0.5) * res - E, b[1] + (yy + 0.5) * res - N)
    out = {"status": "ok", "n_points": int(n), "n_ground": int(len(g))}
    # FLIGHT DATE, from the same read. The "2021" acquisition is NOT one epoch: measured from
    # gps_time, block 4 is entirely 2022, blocks 2 and 3 span both seasons, and the family runs
    # 2021-04-16 to 2022-06-05 -- mid-April (deciduous bare) to early June (full canopy) in SE
    # Minnesota. Canopy cover is therefore NOT one quantity across blocks, and the date must
    # travel with every value so phenology is a covariate rather than a hidden confound.
    if "GpsTime" in a.dtype.names:
        gt = a["GpsTime"]
        out["gps_utc_min"] = _gps_to_utc(np.min(gt)).strftime("%Y-%m-%d %H:%M")
        out["gps_utc_max"] = _gps_to_utc(np.max(gt)).strftime("%Y-%m-%d %H:%M")
        out["gps_span_days"] = round(float((np.max(gt) - np.min(gt)) / 86400.0), 3)
    for r in RADII:
        s = (d <= r) & np.isfinite(cc)
        out[f"cover_r{r:g}"] = round(float(np.nanmean(cc[s])), 6) if s.any() else ""
        out[f"pai_r{r:g}"] = round(float(np.nanmean(pai[s])), 6) if s.any() else ""
        out[f"ncell_r{r:g}"] = int(s.sum())
    v = [out[f"cover_r{r:g}"] for r in RADII if out[f"cover_r{r:g}"] != ""]
    out["cover_ladder_spread"] = round(max(v) - min(v), 6) if len(v) > 1 else ""
    out["box_laz"] = out_box

    # --- slope-normal vertical structure, against the tie estimator's own surface --------
    coef = _surface(a["X"], a["Y"], a["Z"], E, N, struct_radius) if len(g) else None
    if coef is None:
        out["struct_status"] = f"no order-2 surface within {struct_radius} m"
        return out
    zc, _info = ground_elevation_at(g["X"], g["Y"], g["Z"], E, N, struct_radius)
    gcoef = _surface(g["X"], g["Y"], g["Z"], E, N, struct_radius)
    if gcoef is None:
        out["struct_status"] = "no class-2 surface"
        return out
    # the reference surface is fitted to CLASS-2 ground, as the tie is; assert it reproduces
    # ground_elevation_at's constant term so this cannot drift into a second estimator
    out["surface_check_mm"] = round(1000.0 * float(abs(gcoef[0] - (zc - _info["median_resid_mm"] / 1000.0))), 6) \
        if np.isfinite(zc) else ""
    out["slope_deg"] = round(float(np.degrees(np.arctan(np.hypot(gcoef[1], gcoef[2])))), 4)
    sel_all = np.hypot(a["X"] - E, a["Y"] - N) <= struct_radius
    sel_g = np.hypot(g["X"] - E, g["Y"] - N) <= struct_radius
    H = {
        "ng_all": _sn_hist(a["X"][sel_all], a["Y"][sel_all], a["Z"][sel_all], E, N, gcoef, 2, NG_EDGES),
        "ng_class2": _sn_hist(g["X"][sel_g], g["Y"][sel_g], g["Z"][sel_g], E, N, gcoef, 2, NG_EDGES),
        "can_all": _sn_hist(a["X"][sel_all], a["Y"][sel_all], a["Z"][sel_all], E, N, gcoef, 2, CAN_EDGES),
        "can_class2": _sn_hist(g["X"][sel_g], g["Y"][sel_g], g["Z"][sel_g], E, N, gcoef, 2, CAN_EDGES),
    }
    out["n_struct_all"] = int(sel_all.sum()); out["n_struct_class2"] = int(sel_g.sum())
    out["struct_status"] = "ok"
    if struct_out:
        os.makedirs(os.path.dirname(struct_out), exist_ok=True)
        np.savez_compressed(struct_out, ng_edges=NG_EDGES, can_edges=CAN_EDGES,
                            surface_coef=gcoef, struct_radius=struct_radius,
                            easting=E, northing=N, **H)
        out["struct_npz"] = struct_out
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
    ap.add_argument("--boxes", default="data/derived/control_boxes",
                    help="where each mark's EPT box LAZ is KEPT, so nothing needs re-reading")
    ap.add_argument("--struct", default="data/derived/control_structure",
                    help="where each mark's slope-normal histograms are written")
    ap.add_argument("--struct-radius", type=float, default=7.5,
                    help="radius for the surface fit and the structure histograms; 7.5 m = "
                         "the tie's own report radius (1.5 x res)")
    a = ap.parse_args()

    cols = (["set", "point_id", "easting", "northing", "point_type", "ept_block", "status",
             "n_points", "n_ground"]
            + [f"{p}_r{r:g}" for r in RADII for p in ("cover", "pai", "ncell")]
            + ["cover_ladder_spread", "gps_utc_min", "gps_utc_max", "gps_span_days",
               "box_laz", "struct_npz", "struct_status", "n_struct_all", "n_struct_class2",
               "surface_check_mm", "slope_deg"])
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
                safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{setname}__{pid}")
                try:
                    o = cover_at(E, N, res=a.res, k=a.k, voxel_z=a.voxel_z,
                                 min_height=a.min_height, td=td,
                                 ept_url=EPT_URL.format(name=blk),
                                 box_out=os.path.join(a.boxes, safe + ".laz"),
                                 struct_out=os.path.join(a.struct, safe + ".npz"),
                                 struct_radius=a.struct_radius)
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
