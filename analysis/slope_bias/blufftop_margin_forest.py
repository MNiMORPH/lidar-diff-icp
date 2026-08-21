#!/usr/bin/env python3
"""BLUFFTOP-MARGIN FOREST selection (Task 4, clean canopy-offset readout).

Pick the first 30 m of forest beyond the field/forest edge on the blufftop RIM, excluding
drainages -- i.e. upland forest fringe, NOT floodplain and NOT re-entrant hollows. This is
the clean flat/low-slope forest the floodplain contaminated in flat_forest_offset.py.

Definitions on the 5 m AOI grid (gen2 DEM = z_after):
  field   = OPEN (penetration>=0.45) AND gently sloping (<=5 deg) AND UPLAND
            (large-scale TPI >= 0, i.e. plateau not valley bottom)
  forest  = canopy (penetration < 0.25)
  margin  = forest cell within 30 m (<=6 cells) of a field cell  -> the first 30 m in
  drainage= concave hollow (mid-scale TPI < -0.75 m)  -> EXCLUDED
  SELECT  = forest AND within 30 m of field AND not drainage AND upland

Emits counts, a verification figure, and the median DoD on the selection by distance band.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/blufftop_margin_forest.py
"""
import numpy as np, laspy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter, distance_transform_edt
from lidar_diff_icp.viz import hillshade

RES = 5.0; X0, Y0 = 577492.8, 4882737.6
z = np.load("data/derived/elba_fulldensity/z_after.npy")
slope = np.load("data/derived/elba_fulldensity/slope.npy")
dod = np.load("data/derived/elba_fulldensity/dod.npy")
ny, nx = z.shape

# fill DEM holes (nearest) for TPI/hillshade
zf = z.copy(); nm = ~np.isfinite(zf)
if nm.any():
    zf = zf[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]

# penetration (canopy) on the AOI grid
g = np.zeros(ny*nx, np.int64); t = np.zeros(ny*nx, np.int64)
with laspy.open("data/after/3dep2021_fulldensity.laz") as fh:
    for pts in fh.chunk_iterator(20_000_000):
        cl = np.asarray(pts.classification); keep = cl != 7
        ix = ((np.asarray(pts.x)[keep]-X0)/RES).astype(np.int64)
        iy = ((np.asarray(pts.y)[keep]-Y0)/RES).astype(np.int64)
        ok = (ix>=0)&(ix<nx)&(iy>=0)&(iy<ny); c = iy[ok]*nx+ix[ok]
        t += np.bincount(c, minlength=ny*nx); g += np.bincount(c[(cl[keep]==2)[ok]], minlength=ny*nx)
pen = np.divide(g, np.maximum(t,1)).reshape(ny, nx)

# topographic position: large scale (plateau vs valley) and mid scale (nose vs hollow)
tpi_large = zf - uniform_filter(zf, size=61)     # ~305 m: + = upland, - = valley
tpi_mid   = zf - uniform_filter(zf, size=15)     # ~75 m:  + = nose/rim (SHEDS), - = hollow
upland    = tpi_large >= 0.0
shedding  = tpi_mid >= 0.0                        # divergent local high = least material creep-in
convergent = tpi_mid < -0.75                      # hollows/drainages (for the figure)

field  = (pen >= 0.45) & (slope <= 5) & upland
forest = pen < 0.25
dist_field = distance_transform_edt(~field) * RES            # m to nearest field cell
margin = forest & (dist_field <= 30) & upland & shedding

print(f"field cells:      {field.sum():>7}")
print(f"forest cells:     {forest.sum():>7}")
print(f"convergent cells: {convergent.sum():>7}")
print(f"SELECTED blufftop-margin forest (<=30 m, SHEDDING, upland): {margin.sum()}")

fin = np.isfinite(dod)
print("\n=== median DoD (mm), by distance from field edge (slope RISES inward) ===")
for lo, hi in [(0,10),(10,20),(20,30)]:
    m = margin & fin & (dist_field > lo) & (dist_field <= hi)
    if m.any():
        d = dod[m]*1000
        print(f"  {lo:>2}-{hi:<2} m: n={m.sum():>5}  medDoD={np.median(d):+6.1f} mm  "
              f"NMAD={1.4826*np.median(np.abs(d-np.median(d))):5.0f}  medslope={np.median(slope[m]):.1f} deg")

print("\n=== median DoD (mm) by SLOPE on the selection (lowest-slope = cleanest) ===")
for lo, hi in [(0,3),(3,6),(6,10),(10,15),(15,90)]:
    m = margin & fin & (slope>=lo) & (slope<hi)
    if m.any():
        d = dod[m]*1000
        print(f"  {lo:>2}-{hi:<2} deg: n={m.sum():>5}  medDoD={np.median(d):+6.1f} mm  "
              f"NMAD={1.4826*np.median(np.abs(d-np.median(d))):5.0f}")
m = margin & fin
print(f"\n  ALL: n={m.sum()}  medDoD={np.median(dod[m]*1000):+.1f} mm  "
      f"medslope={np.median(slope[m]):.1f} deg  medpen={np.median(pen[m]):.2f}")

# verification figure
hs = hillshade(zf, RES, X0, Y0, fill_gaps=True)
ext = (X0, X0+nx*RES, Y0, Y0+ny*RES)
fig, ax = plt.subplots(figsize=(11, 13))
ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
ov = np.full((ny, nx, 4), 0.0)
ov[field]   = (0.9, 0.8, 0.1, 0.35)     # fields  yellow
ov[forest & ~margin] = (0.1, 0.5, 0.1, 0.18)  # forest  faint green
ov[convergent] = (0.2, 0.4, 0.9, 0.25)  # hollows/drainages blue
ov[margin]  = (0.9, 0.1, 0.1, 0.9)      # SELECTION red
ax.imshow(ov, extent=ext, origin="lower")
ax.set_title("blufftop-margin forest (red) = first 30 m of forest past the field rim,\n"
             "shedding/divergent only.  fields=yellow, other forest=green, hollows=blue")
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
out = "figures/blufftop_margin_forest.png"
import os; os.makedirs("figures", exist_ok=True)
fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")
np.save("data/derived/elba_fulldensity/blufftop_margin_mask.npy", margin)
