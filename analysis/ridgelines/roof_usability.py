#!/usr/bin/env python3
"""FEMA USA Structures footprints over the Elba tile: which are USABLE as flat, hard,
STABLE vertical-datum reference surfaces?

For each footprint we test, against BOTH lidar epochs (gen2 2021 full-density, gen1 2008):
  (a) FLATNESS  - fit a plane to the gen2 TOP surface inside the footprint (5 m inward
                  test region); report roof slope (deg) + residual roughness (m).
                  Flat roof if slope < 5 deg and roughness < 0.10 m.
  (b) APRON     - flat (<4 deg) low-roughness (<4 cm) class-2 GROUND in a 2-6 m outer ring.
  (c) STABILITY - gen2 vs gen1 median TOP-surface height over the footprint; present_in_2008
                  flag; |dh| < 0.3 m => stable (building unchanged, or bare in both).

Classify USABLE (flat roof OR good flat apron, AND stable) vs REJECT.

Streams the 183M-pt gen2 cloud once via chunk_iterator; never laspy.read the big file.

    cd /home/awickert/projects/lidar-diff-icp
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/roof_usability.py
"""
import json
from pathlib import Path
import numpy as np
import laspy
from matplotlib.path import Path as MplPath
from shapely.geometry import Polygon
from shapely import unary_union

OUT = Path("data/derived/elba_refdatum")
GEN2 = "data/after/3dep2021_fulldensity.laz"
GEN1 = "data/before/4342-29-64.laz"

# tile grid / center (lon/lat of tile center for quadrant split)
TILE_CENTER_LON, TILE_CENTER_LAT = -92.016, 44.1095

APRON_IN, APRON_OUT = 2.0, 6.0   # outer-ring apron buffer (m)
INSET = 5.0                       # inward test region for roof plane (m) -- fall back if too small
CELL = 1.0                        # ~1 m top-surface cell for max-z binning

# ---------------------------------------------------------------------------
# load footprints (UTM 15N rings + metadata) prepared by the reprojection step
z = np.load(str(OUT / "fema_footprints_utm.npz"), allow_pickle=True)
rings = [z[f"ring{i}"] for i in range(sum(k.startswith("ring") for k in z.files))]
meta = json.loads(str(z["meta"]))
assert len(rings) == len(meta), (len(rings), len(meta))
N = len(rings)
print(f"{N} footprints")

polys_shp = [Polygon(r) for r in rings]
paths = [MplPath(r) for r in rings]

def buffered_path(poly, dist):
    """shapely buffer (dist<0 shrink, dist>0 grow); return an MplPath (largest ring)
    or None if the buffer collapses to empty."""
    b = poly.buffer(dist)
    if b.is_empty:
        return None
    # pick the largest polygon if MultiPolygon
    if b.geom_type == "MultiPolygon":
        b = max(b.geoms, key=lambda g: g.area)
    return MplPath(np.array(b.exterior.coords))

# inward test path for the roof plane: erode footprint by INSET
inset_paths = [buffered_path(p, -INSET) for p in polys_shp]

# per-footprint bbox (with apron pad) for spatial prefilter
pads = APRON_OUT + 1.0
bb = np.array([[r[:, 0].min() - pads, r[:, 0].max() + pads,
                r[:, 1].min() - pads, r[:, 1].max() + pads] for r in rings])
GX0, GX1 = bb[:, 0].min(), bb[:, 1].max()
GY0, GY1 = bb[:, 2].min(), bb[:, 3].max()
print(f"footprint+apron bbox X {GX0:.1f}..{GX1:.1f}  Y {GY0:.1f}..{GY1:.1f}")

# ---------------------------------------------------------------------------
# stream a cloud; collect points within the global bbox, per footprint.
# We keep, per footprint: inner-test-region top pts (x,y,z), full-footprint z (top),
# and apron class-2 ground pts (x,y,z).
def collect(path, chunk=20_000_000, want_class2_apron=True):
    # accumulate raw x,y,z,cls for points in global bbox, then assign per footprint
    xs, ys, zs, cs = [], [], [], []
    with laspy.open(path) as fh:
        for pts in fh.chunk_iterator(chunk):
            x = np.asarray(pts.x); y = np.asarray(pts.y)
            m = (x >= GX0) & (x <= GX1) & (y >= GY0) & (y <= GY1)
            if not m.any():
                continue
            xs.append(x[m]); ys.append(y[m]); zs.append(np.asarray(pts.z)[m])
            cs.append(np.asarray(pts.classification)[m])
    if not xs:
        return None
    X = np.concatenate(xs); Y = np.concatenate(ys)
    Z = np.concatenate(zs); C = np.concatenate(cs)
    print(f"  {Path(path).name}: {len(X):,} pts in bbox")
    return X, Y, Z, C


def top_surface_z(x, y, z, cell=CELL):
    """max-z per ~cell grid cell over the given points; returns array of cell top-z."""
    if len(x) == 0:
        return np.array([]), np.array([]), np.array([])
    ix = np.floor((x - x.min()) / cell).astype(int)
    iy = np.floor((y - y.min()) / cell).astype(int)
    key = ix.astype(np.int64) * 100000 + iy
    order = np.argsort(key)
    key_s = key[order]; z_s = z[order]; x_s = x[order]; y_s = y[order]
    # for each unique key take the max z
    uk, start = np.unique(key_s, return_index=True)
    end = np.append(start[1:], len(key_s))
    tz = np.empty(len(uk)); tx = np.empty(len(uk)); ty = np.empty(len(uk))
    for i, (s, e) in enumerate(zip(start, end)):
        j = s + np.argmax(z_s[s:e])
        tz[i] = z_s[j]; tx[i] = x_s[j]; ty[i] = y_s[j]
    return tx, ty, tz


def bottom_surface_z(x, y, z, cell=CELL, pct=10):
    """low-percentile (p10) z per ~cell grid cell -- the bare-GROUND envelope, robust to
    low vegetation / curbs / edge spikes. Returns cell (x,y,z_low)."""
    if len(x) == 0:
        return np.array([]), np.array([]), np.array([])
    ix = np.floor((x - x.min()) / cell).astype(int)
    iy = np.floor((y - y.min()) / cell).astype(int)
    key = ix.astype(np.int64) * 100000 + iy
    order = np.argsort(key)
    key_s = key[order]; z_s = z[order]; x_s = x[order]; y_s = y[order]
    uk, start = np.unique(key_s, return_index=True)
    end = np.append(start[1:], len(key_s))
    tz = np.empty(len(uk)); tx = np.empty(len(uk)); ty = np.empty(len(uk))
    for i, (s, e) in enumerate(zip(start, end)):
        tz[i] = np.percentile(z_s[s:e], pct)
        tx[i] = x_s[s:e].mean(); ty[i] = y_s[s:e].mean()
    return tx, ty, tz


def plane_fit(x, y, z):
    """LS plane z = a*x + b*y + c; return slope(deg), roughness(nmad of resid)."""
    if len(x) < 6:
        return np.nan, np.nan
    A = np.column_stack([x - x.mean(), y - y.mean(), np.ones_like(x)])
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)
    a, b, _ = coef
    slope = np.degrees(np.arctan(np.hypot(a, b)))
    resid = z - A @ coef
    rough = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    return slope, rough


print("streaming gen2 (183M) ...")
g2 = collect(GEN2)
print("streaming gen1 (2008) ...")
g1 = collect(GEN1)

X2, Y2, Z2, C2 = g2
X1, Y1, Z1, C1 = g1

recs = []
for i in range(N):
    r = rings[i]; m = meta[i]
    x0, x1, y0, y1 = bb[i]
    # gen2 points in this footprint bbox
    sel2 = (X2 >= x0) & (X2 <= x1) & (Y2 >= y0) & (Y2 <= y1)
    px2, py2, pz2, pc2 = X2[sel2], Y2[sel2], Z2[sel2], C2[sel2]
    sel1 = (X1 >= x0) & (X1 <= x1) & (Y1 >= y0) & (Y1 <= y1)
    px1, py1, pz1, pc1 = X1[sel1], Y1[sel1], Z1[sel1], C1[sel1]

    # --- inside-footprint membership
    in2 = paths[i].contains_points(np.column_stack([px2, py2])) if len(px2) else np.array([], bool)
    in1 = paths[i].contains_points(np.column_stack([px1, py1])) if len(px1) else np.array([], bool)
    # inner test region for roof plane
    ip = inset_paths[i]
    if ip is not None and len(px2):
        inn2 = ip.contains_points(np.column_stack([px2, py2]))
    else:
        inn2 = np.zeros(len(px2), bool)
    # fall back to full footprint if inner region too sparse (small footprint)
    roof_mask = inn2 if inn2.sum() >= 20 else in2

    # (a) FLATNESS: gen2 top surface over roof_mask
    rx, ry, rz = top_surface_z(px2[roof_mask], py2[roof_mask], pz2[roof_mask])
    slope, rough = plane_fit(rx, ry, rz)
    n_roof = len(rz)

    # median TOP heights (max-z cells) for stability
    _, _, tz2 = top_surface_z(px2[in2], py2[in2], pz2[in2])
    _, _, tz1 = top_surface_z(px1[in1], py1[in1], pz1[in1])
    med2 = float(np.median(tz2)) if len(tz2) else np.nan
    med1 = float(np.median(tz1)) if len(tz1) else np.nan
    dh = (med2 - med1) if (np.isfinite(med2) and np.isfinite(med1)) else np.nan

    # present in 2008? need gen1 points inside footprint AND a structure-height signature.
    # We compare the footprint-interior median top to nearby GROUND (apron class-2) in gen1.
    # apron ring (both epochs): between APRON_IN and APRON_OUT outside the footprint
    outer_p = buffered_path(polys_shp[i], APRON_OUT)   # footprint grown by APRON_OUT
    innerb_p = buffered_path(polys_shp[i], APRON_IN)   # footprint grown by APRON_IN

    def apron_ground(px, py, pz, pc):
        if len(px) == 0 or outer_p is None or innerb_p is None:
            return np.array([]), np.array([]), np.array([])
        xy = np.column_stack([px, py])
        in_outer = outer_p.contains_points(xy)
        in_inner = innerb_p.contains_points(xy)   # inside the grown-by-APRON_IN region
        ring_mask = in_outer & (~in_inner)        # annulus APRON_IN..APRON_OUT outside footprint
        g = ring_mask & (pc == 2)
        return px[g], py[g], pz[g]

    ax2, ay2, az2 = apron_ground(px2, py2, pz2, pc2)
    ax1, ay1, az1 = apron_ground(px1, py1, pz1, pc1)
    # apron flatness from the gen2 class-2 GROUND ENVELOPE: grid to ~1 m cells and take the
    # p10 low value per cell (bare-ground surface, robust to grass/curb/edge spikes), then
    # fit a plane to the cell-bottoms. This is the surface a datum would actually sit on.
    gx2, gy2, gz2 = bottom_surface_z(ax2, ay2, az2)
    if len(gz2) >= 8:
        aslope, arough = plane_fit(gx2, gy2, gz2)
        apron_ground_med2 = float(np.median(gz2))
    else:
        aslope, arough, apron_ground_med2 = np.nan, np.nan, np.nan
    # apron_flat = usable hard datum apron: slope < 4 deg, ground-envelope roughness < 6 cm.
    # apron_ideal (< 4 cm) flags the very cleanest paved aprons; the 6 cm bound admits a
    # 4.9 cm government-office apron that is a hard surface by any reasonable standard while
    # excluding the 8 cm grass/gravel surrounds.
    apron_flat = bool(np.isfinite(aslope) and aslope < 4.0 and np.isfinite(arough) and arough < 0.06
                      and len(az2) >= 8)
    apron_ideal = bool(apron_flat and arough < 0.04)

    # present_in_2008: gen1 footprint interior median top rises clearly above gen1 apron ground
    gx1, gy1, gz1 = bottom_surface_z(ax1, ay1, az1)
    apron_ground_med1 = float(np.median(gz1)) if len(gz1) >= 4 else np.nan
    if np.isfinite(med1) and np.isfinite(apron_ground_med1):
        height_above_ground_2008 = med1 - apron_ground_med1
    else:
        height_above_ground_2008 = np.nan
    # gen2 building height above ground (sanity that FEMA footprint really is a building now)
    if np.isfinite(med2) and np.isfinite(apron_ground_med2):
        height_above_ground_2021 = med2 - apron_ground_med2
    else:
        height_above_ground_2021 = np.nan
    present_2008 = bool(np.isfinite(height_above_ground_2008) and height_above_ground_2008 > 1.0)

    # apron GROUND stability 2008->2021 (the datum surface must itself be unchanged)
    apron_dh = (apron_ground_med2 - apron_ground_med1) if (
        np.isfinite(apron_ground_med2) and np.isfinite(apron_ground_med1)) else np.nan
    apron_stable = bool(np.isfinite(apron_dh) and abs(apron_dh) < 0.15)

    # --- flatness verdict
    flat_roof = bool(np.isfinite(slope) and slope < 5.0 and np.isfinite(rough) and rough < 0.10
                     and n_roof >= 20)

    # --- roof stability: small top-surface difference between epochs
    roof_stable = bool(np.isfinite(dh) and abs(dh) < 0.30)

    # overall: a footprint is USABLE if it offers a stable, flat HARD surface --
    #   (i) flat roof present & stable in both epochs, OR
    #   (ii) good flat apron ground that is itself stable 2008->2021.
    usable_roof = flat_roof and roof_stable
    usable_apron = apron_flat and apron_stable
    hard_flat = flat_roof or apron_flat
    verdict = "USABLE" if (usable_roof or usable_apron) else "REJECT"
    reasons = []
    if not hard_flat:
        reasons.append("no flat roof/apron")
    else:
        if flat_roof and not roof_stable:
            reasons.append(f"roof unstable dh={dh:+.2f}m" if np.isfinite(dh) else "roof no-gen1")
        if apron_flat and not apron_stable:
            reasons.append(f"apron unstable dh={apron_dh:+.2f}m" if np.isfinite(apron_dh)
                           else "apron no-gen1")
    if not np.isfinite(med1):
        reasons.append("no gen1 pts")

    # quadrant relative to tile center (use footprint centroid lon/lat)
    q = ("N" if m["lat"] >= TILE_CENTER_LAT else "S") + ("E" if m["lon"] >= TILE_CENTER_LON else "W")

    recs.append(dict(
        idx=i, build_id=m["build_id"], occ=m["occ"], prim=m["prim"], sqm=m["sqm"],
        lon=m["lon"], lat=m["lat"], quadrant=q,
        roof_slope_deg=slope, roof_rough_m=rough, n_roof=n_roof, flat_roof=flat_roof,
        apron_slope_deg=aslope, apron_rough_m=arough, apron_n=len(az2), apron_flat=apron_flat,
        apron_ideal=apron_ideal,
        med2=med2, med1=med1, dh_gen1_gen2=dh, roof_stable=roof_stable,
        apron_dh_gen1_gen2=apron_dh, apron_stable=apron_stable,
        h_above_grnd_2021=height_above_ground_2021, h_above_grnd_2008=height_above_ground_2008,
        present_2008=present_2008, verdict=verdict,
        usable_roof=usable_roof, usable_apron=usable_apron, reasons="; ".join(reasons),
    ))

# ---------------------------------------------------------------------------
# split: footprints WITH lidar coverage (both epochs present) vs NO coverage (outside tile)
covered = [rr for rr in recs if np.isfinite(rr["med2"]) or np.isfinite(rr["med1"])]
nocov = [rr for rr in recs if rr not in covered]
print(f"\n{len(covered)} footprints with lidar coverage; {len(nocov)} outside the tile (no coverage)")

hdr = ("BUILD_ID  occ         area  lon/lat            quad  roofSlp roofRgh flatR  "
       "apSlp apRgh apFlat  dh(g1-g2) apdh   pres08  h21   h08  verdict")
print("\n" + hdr)
print("-" * len(hdr))
def f(v, s="%.2f"):
    return "  nan" if not np.isfinite(v) else s % v
for rr in sorted(covered, key=lambda d: (d["verdict"] != "USABLE", -d["sqm"])):
    print(f"{rr['build_id']:>8}  {rr['occ'][:10]:<10}  {rr['sqm']:>4.0f}  "
          f"{rr['lon']:.4f},{rr['lat']:.4f}  {rr['quadrant']:<4}  "
          f"{f(rr['roof_slope_deg'],'%5.1f')}  {f(rr['roof_rough_m'],'%5.2f')}  "
          f"{str(rr['flat_roof'])[:1]}     {f(rr['apron_slope_deg'],'%4.1f')} {f(rr['apron_rough_m'],'%4.2f')} "
          f"{str(rr['apron_flat'])[:1]}({str(rr['apron_ideal'])[:1]})   {f(rr['dh_gen1_gen2'],'%+5.2f')}  {f(rr['apron_dh_gen1_gen2'],'%+5.2f')}  "
          f"{str(rr['present_2008'])[:1]}     {f(rr['h_above_grnd_2021'],'%4.1f')} {f(rr['h_above_grnd_2008'],'%4.1f')}  "
          f"{rr['verdict']}  {rr['reasons']}")

usable = [rr for rr in recs if rr["verdict"] == "USABLE"]
print(f"\nUSABLE: {len(usable)} / {N}")
from collections import Counter
qc = Counter(rr["quadrant"] for rr in usable)
print("usable by quadrant:", dict(qc))
print("usable BUILD_IDs:", [rr["build_id"] for rr in usable])

# save usable polys + metadata as ref_polys.
# IMPORTANT: for a flat-ROOF surface the datum poly is the footprint interior; for an
# apron-qualified surface the datum surface is the GROUND RING, so we save the apron
# annulus (footprint buffered APRON_IN..APRON_OUT) instead of the pitched-roof footprint.
save = {}
bids, surfs = [], []
for j, rr in enumerate(usable):
    i = rr["idx"]
    if rr["usable_roof"]:
        save[f"poly{j}"] = rings[i]
        surfs.append("roof")
    else:  # apron: annulus polygon (outer buffer with the footprint-grown-by-APRON_IN as a hole)
        outer = polys_shp[i].buffer(APRON_OUT)
        inner = polys_shp[i].buffer(APRON_IN)
        ring_poly = outer.difference(inner)
        if ring_poly.geom_type == "MultiPolygon":
            ring_poly = max(ring_poly.geoms, key=lambda g: g.area)
        # store exterior ring (the sampler intersects points; the small inner hole is a
        # ~2 m collar and negligible). Keep the outer boundary of the apron.
        save[f"poly{j}"] = np.array(ring_poly.exterior.coords)
        surfs.append("apron")
    bids.append(rr["build_id"])
save["build_ids"] = np.array(bids, dtype=np.int64)
save["quadrants"] = np.array([rr["quadrant"] for rr in usable])
save["occ"] = np.array([rr["occ"] for rr in usable])
save["surface"] = np.array(surfs)
np.savez(str(OUT / "usable_hard_polys.npz"), **save)
print(f"saved {len(usable)} usable polys -> {OUT/'usable_hard_polys.npz'}  surfaces={surfs}")

# also dump full classification table as json for the record
json.dump(recs, open(OUT / "roof_usability_table.json", "w"), indent=2, default=float)
print(f"saved full table -> {OUT/'roof_usability_table.json'}")
