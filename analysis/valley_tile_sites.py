#!/usr/bin/env python3
"""Candidate RTK-GNSS survey sites on public land for ANY valley tile: the Elba screener
(analysis/field_sites_rtk.py) re-run on a data/derived/<tile>/ directory written by
scripts/run_valley_tile.py.

    ./lidar-icp/bin/python analysis/valley_tile_sites.py --tile data/derived/4342-28-64
    ./lidar-icp/bin/python analysis/valley_tile_sites.py --tile data/derived/elbaext_delong_star_nosurf \\
        --lowveg data/derived/elbaext_independent/lowveg_grid.npy --compare data/derived/field_sites_rtk.csv
    ./lidar-icp/bin/python analysis/valley_tile_sites.py --all

CLASSES. F1, F2, F3 and H2 are the Elba criteria VERBATIM (field_sites_rtk.py lines 76-84):

    flat = fp & (sl <= 2.0) & np.isfinite(dod)
    F1 sedge-correction test            flat & public & (span_pool >= 400) & (dod < -lod) & (D > 0) & ~pf
    F2 closed canopy (zero penetration) flat & public & pf
    F3 bare floodplain control          flat & public & (sp < 140) & (np.abs(dod) < 30)
    H2 steep slope, gen1-low            fin & public & (sl > 25) & (dod < -lod)

where dod, lod, sp (gen1_spread_mm) are in mm, sl in degrees, pf is penetration_failure,
span_pool is the pooled 3x3 gen1 span rebuilt from beam_offset_table.parquet exactly as
field_sites_rtk.py builds it (lines 67-74: 20 mm bins over +/-5000 mm, 3x3 uniform_filter
x9, p90-p10 of the pooled histogram), and D is the sedge-correction delta
floodplain_vegcorr_delta.npy. Valley tiles carry NO sedge correction, so the (D > 0) term
is applied only when that layer exists and its absence is printed. If the beam table is
absent, span_pool falls back to the per-cell gen1_spread_mm (printed; NOT the Elba
criterion). The cluster rule is Elba's: 8-connected components, "cluster minimum 40 cells
(0.1 ha) so a site holds a few dozen RTK points; top 4 per class by size".

Elba's H1 (forest, PyForestScan forest_pfs) and H3 (open, open_pfs) cannot be built here:
run_valley_tile.py does not produce the PyForestScan layers. In their place, the gen2
near-ground composition lowveg_grid.npy (scripts/make_lowveg.py) stands in, following the
inline analysis in analysis/FIELD_SITES_RTK.md ("gen2 lowveg on the Elba hillslopes"):

    H3 low-lowveg hillslope control  fin & public & (sl >= 10) & (sl <= 35) & (lowveg <= Q1) & (|dod| < 0.5 lod)
    H4 gen2-high in low vegetation   fin & public & (sl >= 10) & (sl <= 35) & (lowveg >= Q3) & (dod > lod)

Q1 and Q3 are the 25th and 75th percentiles of lowveg over the tile's OWN public hillslope
cells (10-35 deg, outside the floodplain domain, finite DoD and lowveg), the population
the inline analysis used (Elba: 169,697 cells, quartiles 0.068 / 0.161 / 0.339). The
quartile edges and that population are MINE (a PROPOSAL), and are printed as such.
H3 AND H4 ARE DEFINED DIFFERENTLY FROM ELBA'S H3 AND ARE NOT COMPARABLE TO IT: Elba's H3 is
"open cover" by the PyForestScan voxel metric; this H3 is "lowest-quartile lowveg", a
different quantity (canopy cover vs lowveg r 0.035 on the Elba hillslopes).

FLOODPLAIN DOMAIN: floodplain_mask_clubb.npy when present (decision D1, the Clubb et al.
2017 valley floor), else the pipeline's convexity floodplain_mask.npy; the one used and the
hectares of each are printed.

PUBLIC LAND: WMA units and state parks (data/derived/valley/public_land/{wma_units,
state_parks}.geojson, shapely contains_xy on cell centres). Distances to WMA parking lots
(facility_subgroup == "Parking Lot" in wma_parking.geojson) and WMA roads (wma_roads.geojson)
are EXACT shapely distances from the cluster's cell centres (nearest cell), not a raster
transform. Elba's CSV measured them from the cluster centroid; --compare reports both.

OUTPUT: data/derived/<tile>/sites_<tile>.csv with field_sites_rtk.csv's columns plus
tile, lowveg (cluster median), fingerprint_score (cluster median; NaN when
fingerprint_score.npy is absent) and domain (which floodplain mask). Every input file,
criterion and cell count is printed above the table. --all loops over every
data/derived/AAAA-BB-CC directory holding dod.npy, skips (and lists) the ones with missing
layers, and concatenates into data/derived/valley/sites_valley_lidar.csv.
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd
import pyproj
import shapely
from scipy.ndimage import binary_dilation, label, uniform_filter
from shapely.geometry import shape

MIN_CELLS = 40; TOP = 4            # Elba's, verbatim (field_sites_rtk.py line 49)
PUBLIC_DEFAULT = "data/derived/valley/public_land"
REQUIRED = ["dod.npy", "lod.npy", "slope.npy", "gen1_spread_mm.npy", "penetration_failure.npy", "corrections.json"]
OPTIONAL = ["floodplain_mask_clubb.npy", "floodplain_mask.npy", "beam_offset_table.parquet", "floodplain_vegcorr_delta.npy",
            "dod_floodplain_vegcorr.npy", "canopy_cover_pfs.npy", "fingerprint_score.npy", "forest_pfs.npy", "open_pfs.npy"]


def pooled_span(parquet, ny, nx, rows_per_chunk=100):
    """field_sites_rtk.py lines 67-74, computed in row chunks (1-row halo) so the
    ny*nx*500-bin histogram never sits in memory at once; identical cell values."""
    import pyarrow.parquet as pq
    LO, HI, DZ = -5000.0, 5000.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
    t = pq.read_table(parquet, columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool)
    cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
    cell = cell[ok]; b = np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1); del t, d, g, ok
    row = cell // nx; order = np.argsort(row, kind="stable"); cell = cell[order]; b = b[order]; row = row[order]
    starts = np.searchsorted(row, np.arange(ny + 1)); out = np.full((ny, nx), np.nan)
    for r0 in range(0, ny, rows_per_chunk):
        r1 = min(ny, r0 + rows_per_chunk); h0 = max(0, r0 - 1); h1 = min(ny, r1 + 1); s = slice(starts[h0], starts[h1])
        H = np.bincount((cell[s] - h0 * nx) * NB + b[s], minlength=(h1 - h0) * nx * NB).astype(np.int32).reshape(h1 - h0, nx, NB)
        Hp = np.rint(uniform_filter(H.astype(np.float32), size=(3, 3, 1), mode="nearest") * 9); C = np.cumsum(Hp, 2); tot = C[..., -1]
        q = lambda p: ctr[np.argmax(C >= p * tot[..., None], axis=2)]
        sp = np.where(tot > 0, q(0.90) - q(0.10), np.nan); out[r0:r1] = sp[r0 - h0:r0 - h0 + (r1 - r0)]
    return out, cell.size


def load_public(pub):
    gj = lambda n: json.load(open(f"{pub}/{n}.geojson"))["features"]
    wma_f = gj("wma_units"); park_f = gj("state_parks"); fac = gj("wma_parking"); rd = gj("wma_roads")
    polys = [shape(f["geometry"]) for f in wma_f] + [shape(f["geometry"]) for f in park_f]
    names = [f["properties"].get("name") for f in wma_f] + [f"state park: {f['properties'].get('unit_name')}" for f in park_f]
    parking = [shape(f["geometry"]) for f in fac if f["properties"].get("facility_subgroup") == "Parking Lot"]
    roads = [shape(f["geometry"]) for f in rd]
    return shapely.union_all(polys), names, parking, roads, len(fac)


def screen(tile_dir, pub=PUBLIC_DEFAULT, lowveg_path=None, compare=None, out=None):
    tile = os.path.basename(os.path.normpath(tile_dir)); P = lambda n: f"{tile_dir}/{n}"
    lowveg_path = lowveg_path or P("lowveg_grid.npy")
    missing = [n for n in REQUIRED if not os.path.exists(P(n))] + ([] if os.path.exists(lowveg_path) else ["lowveg_grid.npy"])
    if not (os.path.exists(P("floodplain_mask_clubb.npy")) or os.path.exists(P("floodplain_mask.npy"))): missing.append("floodplain_mask(_clubb).npy")
    if missing: return None, missing
    print(f"\n{'=' * 100}\nTILE {tile}  ({tile_dir})\n{'=' * 100}")
    print("INPUT FILES AND ROLES")
    c = json.load(open(P("corrections.json"))); X0, Y0, X1, Y1 = c["bounds"]; res = c["res_m"]
    dod = np.load(P("dod.npy")) * 1e3; lod = np.load(P("lod.npy")) * 1e3; sl = np.load(P("slope.npy"))
    sp = np.load(P("gen1_spread_mm.npy")); pf = np.load(P("penetration_failure.npy")).astype(bool); lv = np.load(lowveg_path)
    ny, nx = dod.shape; yy, xx = np.mgrid[0:ny, 0:nx]; E = X0 + (xx + 0.5) * res; N = Y0 + (yy + 0.5) * res
    for n, a, role in (("dod.npy", dod, "DoD = gen2 - gen1, m -> mm"), ("lod.npy", lod, "level of detection, m -> mm"), ("slope.npy", sl, "slope, deg"),
                       ("gen1_spread_mm.npy", sp, "per-cell gen1 p90-p10 span, mm (F3's sp)"), ("penetration_failure.npy", pf, "no gen1 return at the gen2 surface (pf)")):
        print(f"  {P(n):70s} {role}; shape {a.shape}; finite {np.isfinite(a).sum() if a.dtype.kind == 'f' else a.sum():,}")
    print(f"  {lowveg_path:70s} gen2 lowveg (near-ground composition); shape {lv.shape}; finite {np.isfinite(lv).sum():,}")
    assert lv.shape == dod.shape, f"lowveg grid {lv.shape} is not on the tile grid {dod.shape}"
    print(f"  {P('corrections.json'):70s} grid bounds {c['bounds']}, res {res} m -> {ny} x {nx} cells, {ny * nx * res * res / 1e4:.0f} ha")
    # floodplain domain (D1)
    own = np.load(P("floodplain_mask.npy")).astype(bool) if os.path.exists(P("floodplain_mask.npy")) else None
    clubb = np.load(P("floodplain_mask_clubb.npy")).astype(bool) if os.path.exists(P("floodplain_mask_clubb.npy")) else None
    ha = lambda m: m.sum() * res * res / 1e4
    print("FLOODPLAIN DOMAIN")
    if own is not None: print(f"  {P('floodplain_mask.npy'):70s} pipeline convexity mask: {ha(own):.1f} ha")
    if clubb is not None: print(f"  {P('floodplain_mask_clubb.npy'):70s} Clubb et al. 2017 valley floor (D1): {ha(clubb):.1f} ha")
    fp, domain = (clubb, "clubb") if clubb is not None else (own, "pipeline")
    print(f"  USED: {domain} ({'floodplain_mask_clubb.npy, decision D1' if domain == 'clubb' else 'floodplain_mask.npy; no Clubb layer in this directory'}), {ha(fp):.1f} ha")
    # optional layers
    have = {n: os.path.exists(P(n)) for n in OPTIONAL}
    D = np.load(P("floodplain_vegcorr_delta.npy")) * 1e3 if have["floodplain_vegcorr_delta.npy"] else None
    new = np.load(P("dod_floodplain_vegcorr.npy")) * 1e3 if have["dod_floodplain_vegcorr.npy"] else None
    cc = np.load(P("canopy_cover_pfs.npy")) if have["canopy_cover_pfs.npy"] else None
    fs = np.load(P("fingerprint_score.npy")) if have["fingerprint_score.npy"] else None
    print("OPTIONAL LAYERS")
    if have["beam_offset_table.parquet"]:
        span_pool, nret = pooled_span(P("beam_offset_table.parquet"), ny, nx)
        print(f"  {P('beam_offset_table.parquet'):70s} pooled 3x3 gen1 span (F1's span_pool), Elba recipe, {nret:,} in-grid returns; finite cells {np.isfinite(span_pool).sum():,}")
    else:
        span_pool = sp; print(f"  beam_offset_table.parquet ABSENT: F1's span_pool = per-cell gen1_spread_mm (NOT the Elba criterion)")
    print(f"  floodplain_vegcorr_delta.npy {'present: F1 keeps the (D > 0) term' if D is not None else 'ABSENT: F1 has no (D > 0) term (valley tiles carry no sedge correction)'}")
    print(f"  dod_floodplain_vegcorr.npy   {'present: corrected_dod_mm filled' if new is not None else 'ABSENT: corrected_dod_mm is NaN'}")
    print(f"  canopy_cover_pfs.npy         {'present: canopy_cover filled' if cc is not None else 'ABSENT (no PyForestScan on valley tiles): canopy_cover is NaN'}")
    print(f"  fingerprint_score.npy        {'present' if fs is not None else 'ABSENT: fingerprint_score is NaN'}")
    # public land
    wma_u, names, parking, roads, nfac = load_public(pub)
    public = shapely.contains_xy(wma_u, E.ravel(), N.ravel()).reshape(ny, nx)
    bx = shapely.box(X0, Y0, X1, Y1); park_in = [p for p in parking if bx.contains(p)]; roads_in = [r for r in roads if bx.intersects(r)]
    print("PUBLIC LAND")
    print(f"  {pub}/wma_units.geojson + state_parks.geojson: {len(names)} polygons ({names}); public cells on the grid {public.sum():,} of {ny * nx:,} ({100 * public.mean():.1f}%)")
    print(f"  {pub}/wma_parking.geojson: {len(parking)} of {nfac} facilities are facility_subgroup == 'Parking Lot'; {len(park_in)} inside the tile bounds (distances use ALL {len(parking)})")
    print(f"  {pub}/wma_roads.geojson: {len(roads)} roads, {len(roads_in)} touching the tile bounds (distances use all)")
    # hillslope lowveg quartiles (PROPOSAL)
    fin = np.isfinite(dod); hs = fin & public & (sl >= 10) & (sl <= 35) & ~fp & np.isfinite(lv)
    q1, q2, q3 = np.percentile(lv[hs], [25, 50, 75]) if hs.sum() else (np.nan,) * 3
    print("LOWVEG QUARTILES (PROPOSAL, mine)")
    print(f"  population: public hillslope cells, 10 <= slope <= 35 deg, outside the {domain} floodplain, finite DoD and lowveg: {hs.sum():,} cells")
    print(f"  edges Q1 / median / Q3 = {q1:.3f} / {q2:.3f} / {q3:.3f}  (Elba's inline analysis: 169,697 cells, 0.068 / 0.161 / 0.339)")
    hs_all = fin & (sl >= 10) & (sl <= 35) & ~fp & np.isfinite(lv)
    if hs_all.sum(): print(f"  for reference, the same quartiles over ALL hillslope cells in the tile ({hs_all.sum():,}): {' / '.join(f'{v:.3f}' for v in np.percentile(lv[hs_all], [25, 50, 75]))}")
    # classes
    flat = fp & (sl <= 2.0) & fin
    f1 = flat & public & (span_pool >= 400) & (dod < -lod) & ~pf
    if D is not None: nD = int((f1 & ~(D > 0)).sum()); f1 = f1 & (D > 0); print(f"  F1's (D > 0) term removes {nD:,} of {f1.sum() + nD:,} cells on this tile")
    classes = {
        "F1 sedge-correction test": (f1, "flat & public & (span_pool >= 400) & (dod < -lod)" + (" & (D > 0)" if D is not None else "") + " & ~pf"),
        "F2 closed canopy (zero penetration)": (flat & public & pf, "flat & public & pf"),
        "F3 bare floodplain control": (flat & public & (sp < 140) & (np.abs(dod) < 30), "flat & public & (sp < 140) & (|dod| < 30)"),
        "H2 steep slope, gen1-low": (fin & public & (sl > 25) & (dod < -lod), "fin & public & (sl > 25) & (dod < -lod)"),
        "H3 low-lowveg hillslope control": (fin & public & (sl >= 10) & (sl <= 35) & (lv <= q1) & (np.abs(dod) < 0.5 * lod), f"fin & public & (10 <= sl <= 35) & (lowveg <= Q1 {q1:.3f}) & (|dod| < 0.5 lod)"),
        "H4 gen2-high in low vegetation": (fin & public & (sl >= 10) & (sl <= 35) & (lv >= q3) & (dod > lod), f"fin & public & (10 <= sl <= 35) & (lowveg >= Q3 {q3:.3f}) & (dod > lod)"),
    }
    print("CRITERIA  (flat = floodplain domain & slope <= 2 deg & finite DoD: %s cells; fin = finite DoD: %s cells)" % (f"{flat.sum():,}", f"{fin.sum():,}"))
    print("  H3 AND H4 USE gen2 lowveg, NOT PyForestScan cover: they are DEFINED DIFFERENTLY from Elba's H1/H3 and are NOT COMPARABLE to them.")
    to_ll = pyproj.Transformer.from_crs(c.get("crs", "EPSG:26915"), "EPSG:4326", always_xy=True)
    park_u = shapely.union_all(parking) if parking else None; road_u = shapely.union_all(roads) if roads else None
    rows = []
    for cname, (m, crit) in classes.items():
        lab, n = label(m, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); order = [k for k in np.argsort(sizes[1:])[::-1] + 1 if sizes[k] >= MIN_CELLS][:TOP]
        print(f"\n{cname}: {crit}\n  {m.sum():,} cells; clusters >= {MIN_CELLS} cells: {int((sizes[1:] >= MIN_CELLS).sum())}; top {TOP} listed"
              + (f"; largest cluster {sizes[1:].max() if n else 0} cells" if len(order) < TOP else ""))
        for r, k in enumerate(order):
            mm = lab == k; e0, n0 = E[mm].mean(), N[mm].mean(); lon, lat = to_ll.transform(e0, n0)
            pts = shapely.MultiPoint(np.column_stack([E[mm], N[mm]]))
            dpark = shapely.distance(pts, park_u) if park_u is not None else np.nan; droad = shapely.distance(pts, road_u) if road_u is not None else np.nan
            edge = (binary_dilation(mm, structure=np.ones((3, 3))) & ~mm & flat & ~pf).sum() if cname.startswith("F2") else np.nan
            rows.append(dict(site=f"{cname.split()[0]}-{r + 1}", cls=cname, cells=int(mm.sum()), ha=round(mm.sum() * res * res / 1e4, 2),
                             easting=round(e0, 1), northing=round(n0, 1), lat=round(lat, 6), lon=round(lon, 6),
                             raw_dod_mm=round(float(np.median(dod[mm])), 0), corrected_dod_mm=round(float(np.median(new[mm])), 0) if (new is not None and cname.startswith("F")) else np.nan,
                             gen1_span_mm=round(float(np.nanmedian(sp[mm])), 0), slope_deg=round(float(np.median(sl[mm])), 1),
                             canopy_cover=round(float(np.nanmedian(cc[mm])), 2) if cc is not None else np.nan, penfail_pct=round(100 * float(pf[mm].mean()), 0),
                             edge_cells=edge, dist_parking_m=round(dpark, 0), dist_wma_road_m=round(droad, 0),
                             tile=tile, lowveg=round(float(np.nanmedian(lv[mm])), 3), fingerprint_score=round(float(np.nanmedian(fs[mm])), 3) if fs is not None else np.nan, domain=domain))
    cols = ["site", "cls", "cells", "ha", "easting", "northing", "lat", "lon", "raw_dod_mm", "corrected_dod_mm", "gen1_span_mm", "slope_deg", "canopy_cover",
            "penfail_pct", "edge_cells", "dist_parking_m", "dist_wma_road_m", "tile", "lowveg", "fingerprint_score", "domain"]
    F = pd.DataFrame(rows, columns=cols); out = out or P(f"sites_{tile}.csv"); F.to_csv(out, index=False)
    pd.set_option("display.width", 300); print(f"\nCOLUMNS: dist_* = exact shapely distance from the cluster's nearest CELL CENTRE; lowveg, fingerprint_score = cluster medians")
    print(F.drop(columns=["cls", "tile"]).to_string(index=False) if len(F) else "(no clusters)"); print(f"table: {out}")
    if compare: reproduce(compare, classes, E, N, have, P, fin, public, sl, dod, lod, park_u, road_u)
    return F, []


def reproduce(ref_csv, classes, E, N, have, P, fin, public, sl, dod, lod, park_u, road_u):
    """Elba check: same counts and centroids as field_sites_rtk.csv for F1/F2/F3/H2 (and
    Elba's PFS-defined H1/H3 when forest_pfs/open_pfs are present, computed here for the
    check only)."""
    R = pd.read_csv(ref_csv); print(f"\nREPRODUCTION CHECK against {ref_csv}")
    extra = {}
    if have["forest_pfs.npy"] and have["open_pfs.npy"]:
        forest = np.load(P("forest_pfs.npy")).astype(bool); openc = np.load(P("open_pfs.npy")).astype(bool)
        extra["H1"] = fin & public & (sl >= 10) & (sl <= 35) & forest & (dod > lod)
        extra["H3"] = fin & public & (sl >= 10) & (sl <= 35) & openc & (np.abs(dod) < 0.5 * lod)
        print("  forest_pfs/open_pfs present: Elba's H1/H3 (PyForestScan) computed for the check only, not written")
    ours = {}
    for cname, (m, _) in classes.items():
        k = cname.split()[0]
        if k in ("F1", "F2", "F3", "H2"): ours[k] = m
    ours.update(extra)
    worst = 0.0; nbad = 0
    for k, m in ours.items():
        lab, n = label(m, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); order = [j for j in np.argsort(sizes[1:])[::-1] + 1 if sizes[j] >= MIN_CELLS][:TOP]
        ref = R[R.site.str.startswith(k + "-")].reset_index(drop=True)
        print(f"  {k}: ours {m.sum():,} cells, {len(order)} clusters listed; Elba CSV lists {len(ref)}")
        for r, j in enumerate(order):
            mm = lab == j; e0, n0 = E[mm].mean(), N[mm].mean()
            if r >= len(ref): print(f"    {k}-{r + 1}: ours {mm.sum()} cells E {e0:.1f} N {n0:.1f}; NOT IN Elba CSV"); nbad += 1; continue
            de, dn = e0 - ref.easting[r], n0 - ref.northing[r]; dc = int(mm.sum()) - int(ref.cells[r]); worst = max(worst, abs(de), abs(dn)); nbad += (dc != 0) or (abs(de) > 0.05) or (abs(dn) > 0.05)
            pt = shapely.Point(e0, n0); dp = (pt.distance(park_u) if park_u is not None else np.nan, pt.distance(road_u) if road_u is not None else np.nan)
            print(f"    {k}-{r + 1}: cells {mm.sum()} vs {ref.cells[r]} (d {dc:+d}); centroid dE {de:+.2f} m dN {dn:+.2f} m; centroid dist to lot {dp[0]:.0f} vs {ref.dist_parking_m[r]:.0f}, road {dp[1]:.0f} vs {ref.dist_wma_road_m[r]:.0f}")
        if len(order) < len(ref): print(f"    Elba CSV has {len(ref) - len(order)} more {k} cluster(s) than we list"); nbad += len(ref) - len(order)
    print(f"  RESULT: {'REPRODUCED' if nbad == 0 else f'{nbad} MISMATCH(ES)'}; largest centroid difference {worst:.2f} m")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", help="tile directory data/derived/<tile>")
    ap.add_argument("--public", default=PUBLIC_DEFAULT); ap.add_argument("--lowveg", help="lowveg_grid.npy override (default <tile>/lowveg_grid.npy)")
    ap.add_argument("--compare", help="reference CSV (data/derived/field_sites_rtk.csv) for the Elba reproduction check")
    ap.add_argument("--out", help="output CSV (default <tile>/sites_<tile>.csv)")
    ap.add_argument("--all", action="store_true", help="every data/derived/AAAA-BB-CC with dod.npy -> data/derived/valley/sites_valley_lidar.csv")
    a = ap.parse_args()
    if a.all:
        dirs = sorted(d for d in glob.glob("data/derived/*") if re.fullmatch(r"\d{4}-\d{2}-\d{2}", os.path.basename(d)) and os.path.exists(f"{d}/dod.npy"))
        print(f"--all: {len(dirs)} tile directories with dod.npy: {[os.path.basename(d) for d in dirs]}")
        frames = []; skipped = {}
        for d in dirs:
            F, missing = screen(d, a.public)
            if missing: skipped[os.path.basename(d)] = missing
            elif len(F): frames.append(F)
        out = "data/derived/valley/sites_valley_lidar.csv"; os.makedirs(os.path.dirname(out), exist_ok=True)
        A = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(); A.to_csv(out, index=False)
        print(f"\n--all: {len(frames)} tiles screened, {len(A)} sites -> {out}")
        if skipped: print("--all: SKIPPED (incomplete layers):"); [print(f"  {t}: missing {m}") for t, m in skipped.items()]
    else:
        if not a.tile: ap.error("--tile or --all")
        F, missing = screen(a.tile, a.public, a.lowveg, a.compare, a.out)
        if missing: print(f"{a.tile}: INCOMPLETE, missing {missing}"); sys.exit(1)


if __name__ == "__main__":
    main()
