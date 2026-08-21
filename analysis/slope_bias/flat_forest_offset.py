#!/usr/bin/env python3
"""FLAT-FOREST OFFSET (Task 4) — is there a DoD offset on NEAR-ZERO-SLOPE forested cells?

If the +slope bias is really a leaf-on canopy lift (gen2 ground reads high under May
green-up), it should appear as a CONSTANT positive offset on FLAT forest, where every
slope mechanism (registration x slope, within-cell relief, cloth-on-slope) is ~0 and the
flat-hard datum pins OPEN flat ground to zero. Flat forest is NOT used for the datum, so
it is a free readout of the pure canopy effect.

Uses the full-density DoD (gen1-CSF vs gen2-class2), slope, and gen2 penetration on the
AOI grid. Splits by canopy closure (penetration) for a dose-response.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/flat_forest_offset.py
"""
import numpy as np, laspy

RES = 5.0; X0, Y0 = 577492.8, 4882737.6
dod = np.load("data/derived/elba_fulldensity/dod.npy")
slope = np.load("data/derived/elba_fulldensity/slope.npy")
ny, nx = dod.shape

# gen2 penetration (ground/total) on the AOI grid
g = np.zeros(ny*nx, np.int64); t = np.zeros(ny*nx, np.int64)
with laspy.open("data/after/3dep2021_fulldensity.laz") as fh:
    for pts in fh.chunk_iterator(20_000_000):
        cl = np.asarray(pts.classification); keep = cl != 7
        ix = ((np.asarray(pts.x)[keep]-X0)/RES).astype(np.int64)
        iy = ((np.asarray(pts.y)[keep]-Y0)/RES).astype(np.int64)
        ok = (ix>=0)&(ix<nx)&(iy>=0)&(iy<ny)
        c = iy[ok]*nx+ix[ok]
        t += np.bincount(c, minlength=ny*nx)
        g += np.bincount(c[(cl[keep]==2)[ok]], minlength=ny*nx)
pen = np.divide(g, np.maximum(t,1)).reshape(ny, nx)

fin = np.isfinite(dod)
def med(m):
    return (m.sum(), np.median(dod[m])*1000 if m.any() else np.nan,
            1.4826*np.median(np.abs(dod[m]-np.median(dod[m])))*1000 if m.any() else np.nan)

print("=== median DoD (mm) by slope, FOREST (pen<0.25) vs OPEN (pen>=0.40) ===")
print("  flat-open should be ~0 (datum); a nonzero flat-FOREST value = the canopy offset")
print(f"{'slope':>8} | {'forest n':>8} {'forest mm':>9} {'NMAD':>6} | {'open n':>7} {'open mm':>8}")
for lo, hi in [(0,1),(1,2),(2,3),(3,5),(5,10),(10,15)]:
    s = fin & (slope>=lo) & (slope<hi)
    fn, fm, fnm = med(s & (pen<0.25)); on, om, _ = med(s & (pen>=0.40))
    print(f"{lo:>3}-{hi:<3} | {fn:>8} {fm:>+9.1f} {fnm:>6.0f} | {on:>7} {om:>+8.1f}")

print("\n=== flat forest (slope<3) by canopy closure (penetration) — dose-response ===")
flat = fin & (slope<3)
print(f"{'penetration':>14} {'n':>7} {'medDoD mm':>10}")
for lo, hi in [(0.0,0.10),(0.10,0.20),(0.20,0.30),(0.30,0.45),(0.45,1.01)]:
    m = flat & (pen>=lo) & (pen<hi)
    if m.any():
        print(f"  {lo:.2f}-{hi:.2f}   {m.sum():>7} {np.median(dod[m])*1000:>+10.1f}")
print("\nRead: a positive flat-forest offset that GROWS as penetration drops (denser canopy) "
      "= a real leaf-on canopy lift, present at ZERO slope -> the bias is not purely slope.")
