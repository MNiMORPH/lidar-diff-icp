#!/usr/bin/env python3
"""Test Andy's decomposition: is flat-forest ~0 actually FIELD-EDGE EROSION (negative)
averaging against a real FOREST-INTERIOR / CONVEX-HILLTOP INCREASE (positive)?

Two independent claims:
  (1) field-margin erosion: forest DoD is negative right at the field edge, distinct from
      interior -> the flat-forest ~0 is a mix, not a null.
  (2) convex-hilltop increase: strongly CONVEX, LOW-SLOPE, forested cells FAR from any
      field edge show POSITIVE DoD -> a convex top cannot deposit, so + = gen2-high artifact
      independent of slope -> revives a cover/canopy offset.

If (2) holds, my "no canopy offset at flat" is WRONG. Decisive.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/edge_vs_convex.py
"""
import numpy as np, laspy, os
from scipy.ndimage import uniform_filter, distance_transform_edt

RES = 5.0; X0, Y0 = 577492.8, 4882737.6
z = np.load("data/derived/elba_fulldensity/z_after.npy")
slope = np.load("data/derived/elba_fulldensity/slope.npy")
dod = np.load("data/derived/elba_refdatum/dod_geoid.npy")
ny, nx = z.shape
zf = z.copy(); nm = ~np.isfinite(zf)
if nm.any():
    zf = zf[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]

pf = "data/derived/elba_fulldensity/penetration.npy"
if os.path.exists(pf):
    pen = np.load(pf)
else:
    g = np.zeros(ny*nx, np.int64); t = np.zeros(ny*nx, np.int64)
    with laspy.open("data/after/3dep2021_fulldensity.laz") as fh:
        for pts in fh.chunk_iterator(20_000_000):
            cl = np.asarray(pts.classification); keep = cl != 7
            ix = ((np.asarray(pts.x)[keep]-X0)/RES).astype(np.int64)
            iy = ((np.asarray(pts.y)[keep]-Y0)/RES).astype(np.int64)
            ok=(ix>=0)&(ix<nx)&(iy>=0)&(iy<ny); c=iy[ok]*nx+ix[ok]
            t+=np.bincount(c,minlength=ny*nx); g+=np.bincount(c[(cl[keep]==2)[ok]],minlength=ny*nx)
    pen = np.divide(g, np.maximum(t,1)).reshape(ny, nx); np.save(pf, pen)

# convexity: z minus local mean (+ = convex/local high); plan-ish via mid window
tpi_mid   = zf - uniform_filter(zf, size=9)      # ~45 m: + convex/nose, - concave/hollow
tpi_large = zf - uniform_filter(zf, size=61)
upland = tpi_large >= 0.0
forest = pen < 0.25
field  = (pen >= 0.45) & (slope <= 5) & upland
dist_field = distance_transform_edt(~field) * RES
fin = np.isfinite(dod)

def stat(m):
    d = dod[m]*1000
    return (m.sum(), np.median(d), 1.4826*np.median(np.abs(d-np.median(d))),
            1.2533*1.4826*np.median(np.abs(d-np.median(d)))/np.sqrt(max(m.sum(),1)))

print("=== CLAIM 1: field-margin. Forest DoD vs distance from field edge, slope<6 deg ===")
print(f"{'dist(m)':>10} {'n':>6} {'medDoD':>8} {'NMAD':>6} {'SE':>5}")
for lo,hi in [(0,10),(10,20),(20,40),(40,80),(80,1e9)]:
    m = forest & fin & (slope<6) & (dist_field>lo) & (dist_field<=hi)
    n,md,nd,se = stat(m); lab=f"{lo:.0f}-{hi:.0f}" if hi<1e9 else f">{lo:.0f}"
    print(f"{lab:>10} {n:>6} {md:>+8.1f} {nd:>6.0f} {se:>5.1f}")

print("\n=== CLAIM 2: convex forested hilltops, LOW slope, FAR from field (>40 m) ===")
print(f"{'convexity(TPI_mid)':>18} {'n':>6} {'medDoD':>8} {'NMAD':>6} {'SE':>5}  (slope<6, forest, dist>40)")
base = forest & fin & (slope<6) & (dist_field>40)
for lo,hi,lab in [(-9,-0.3,"concave <-0.3"),(-0.3,0.3,"flat -0.3..0.3"),
                  (0.3,0.8,"convex 0.3..0.8"),(0.8,99,"convex >0.8")]:
    m = base & (tpi_mid>=lo) & (tpi_mid<hi)
    n,md,nd,se = stat(m); print(f"{lab:>18} {n:>6} {md:>+8.1f} {nd:>6.0f} {se:>5.1f}")

print("\n=== cross: strongly CONVEX (TPI_mid>0.3) forest, by slope, far from field ===")
cx = forest & fin & (dist_field>40) & (tpi_mid>0.3)
for lo,hi in [(0,3),(3,6),(6,10),(10,20)]:
    m = cx & (slope>=lo)&(slope<hi); n,md,nd,se=stat(m)
    print(f"  slope {lo:>2}-{hi:<2}: n={n:>5}  medDoD={md:>+7.1f} mm  SE={se:.1f}")

print("\nVERDICT: if convex low-slope forest FAR from field is clearly POSITIVE (>~ +5 mm, "
      ">2*SE), a convex top can't deposit -> gen2-high artifact independent of slope -> "
      "canopy offset REVIVED, and the flat-forest ~0 was edge-erosion masking it.")
