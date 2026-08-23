#!/usr/bin/env python3
"""Forest-cover density map from the 2021 cloud + test: is DENSITY (not slope) the driver
of the gen2-high canopy artifact on forested ridge crests? (Andy's prediction.)

Forest-cover density = FRACTIONAL CANOPY COVER = share of FIRST returns whose height above
the gen2 ground surface exceeds 2 m (standard lidar fractional cover). Also a low-veg
fraction (0.5-2 m) for context. Then, on forested crest pixels:
  (2) DoD vs slope in FINE bins -> locate the rollover where steeper slope yields LESS
      apparent increase;
  (3) DoD vs canopy cover; cover-vs-slope covariation; and the 2-D control DoD(cover x
      slope) + partial correlations to see which variable actually drives DoD.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/forest_density_driver.py
"""
import numpy as np, laspy, os
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RES = 5.0; X0, Y0 = 577492.8, 4882737.6
z = np.load("data/derived/elba_fulldensity/z_after.npy")
slope = np.load("data/derived/elba_fulldensity/slope.npy")
dod = np.load("data/derived/elba_refdatum/dod_geoid.npy")
pen = np.load("data/derived/elba_fulldensity/penetration.npy")
ny, nx = z.shape
zf = z.copy(); nm = ~np.isfinite(zf)
if nm.any():
    zf = zf[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]
zg = zf.ravel()

# ---- fractional canopy cover from the 2021 cloud (height above the gen2 ground) ---------
cf = "data/derived/elba_fulldensity/canopy_cover.npy"
if os.path.exists(cf):
    cover = np.load(cf)
else:
    nf = np.zeros(ny*nx); nfc = np.zeros(ny*nx); nlow = np.zeros(ny*nx); nall = np.zeros(ny*nx)
    with laspy.open("data/after/3dep2021_fulldensity.laz") as fh:
        for pts in fh.chunk_iterator(20_000_000):
            cl = np.asarray(pts.classification); keep = cl != 7
            x = np.asarray(pts.x)[keep]; y = np.asarray(pts.y)[keep]
            zz = np.asarray(pts.z)[keep]; rn = np.asarray(pts.return_number)[keep]
            ix = ((x-X0)/RES).astype(np.int64); iy = ((y-Y0)/RES).astype(np.int64)
            ok = (ix>=0)&(ix<nx)&(iy>=0)&(iy<ny)
            c = iy[ok]*nx+ix[ok]; h = zz[ok]-zg[c]; first = rn[ok]==1
            nall += np.bincount(c, minlength=ny*nx)
            nf  += np.bincount(c[first], minlength=ny*nx)
            nfc += np.bincount(c[first & (h>2.0)], minlength=ny*nx)
            nlow+= np.bincount(c[(h>0.5)&(h<=2.0)], minlength=ny*nx)
    cover = np.divide(nfc, np.maximum(nf,1)).reshape(ny, nx)     # fractional canopy cover 0-1
    np.save(cf, cover)
    np.save("data/derived/elba_fulldensity/lowveg_frac.npy",
            np.divide(nlow, np.maximum(nall,1)).reshape(ny, nx))

# map
import matplotlib
ext = (X0, X0+nx*RES, Y0, Y0+ny*RES)
fig, ax = plt.subplots(figsize=(10, 12))
im = ax.imshow(cover, extent=ext, origin="lower", cmap="YlGn", vmin=0, vmax=1)
fig.colorbar(im, ax=ax, shrink=0.6, label="fractional canopy cover (first returns >2 m)")
ax.set_title("2021 forest-cover density (fractional canopy cover)")
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
fig.savefig("figures/forest_cover_density.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote figures/forest_cover_density.png")

# ---- crest pixels + covariation tests --------------------------------------------------
crest = np.load("data/derived/elba_fulldensity/crest_mask.npy")
fin = np.isfinite(dod)
cm = crest & fin
forest = pen < 0.25
S = slope; C = cover; D = dod*1000

print("\n=== (2) FORESTED-crest DoD vs slope (fine bins) — locate the rollover ===")
print(f"{'slope':>9} {'n':>5} {'medDoD':>8} {'med cover':>9}")
for lo,hi in [(0,2),(2,4),(4,6),(6,8),(8,10),(10,13),(13,17),(17,22),(22,30),(30,90)]:
    m = cm & forest & (S>=lo)&(S<hi)
    if m.sum()>=20:
        print(f"{lo:>4}-{hi:<4} {m.sum():>5} {np.median(D[m]):>+8.1f} {np.median(C[m]):>9.2f}")

print("\n=== cover-vs-slope covariation on forested crests ===")
mf = cm & forest
print(f"  corr(slope, cover) = {np.corrcoef(S[mf], C[mf])[0,1]:+.2f}   "
      f"(n={mf.sum()})   median cover={np.median(C[mf]):.2f}")

print("\n=== (3) DoD by CANOPY COVER quintiles (forested crests), and vs slope ===")
qs = np.quantile(C[mf], [0,.2,.4,.6,.8,1.0])
print(f"{'cover bin':>14} {'n':>5} {'medDoD':>8} {'med slope':>9}")
for i in range(5):
    m = mf & (C>=qs[i]) & (C<=qs[i+1] if i==4 else C<qs[i+1])
    if m.sum()>=20:
        print(f"  {qs[i]:.2f}-{qs[i+1]:.2f}   {m.sum():>5} {np.median(D[m]):>+8.1f} {np.median(S[m]):>9.1f}")

print("\n=== 2-D control: median DoD (mm) by SLOPE (rows) x COVER (cols) ===")
sbin = [0,5,10,20,90]; cbin = [0,.5,.75,.9,1.01]
hdr = "slope\\cover | " + " | ".join(f"{cbin[j]:.2f}-{cbin[j+1]:.2f}" for j in range(4))
print(hdr)
for i in range(4):
    row = f"  {sbin[i]:>2}-{sbin[i+1]:<3}   | "
    for j in range(4):
        m = mf & (S>=sbin[i])&(S<sbin[i+1]) & (C>=cbin[j])&(C<cbin[j+1])
        row += f"{np.median(D[m]):>+6.1f}({m.sum():>4})" if m.sum()>=15 else f"{'  --':>11}"
        row += " | "
    print(row)

# partial correlations (does DoD track cover controlling slope, and vice versa?)
def partial(a, b, ctrl):
    ra = a - np.polyval(np.polyfit(ctrl, a, 1), ctrl)
    rb = b - np.polyval(np.polyfit(ctrl, b, 1), ctrl)
    return np.corrcoef(ra, rb)[0,1]
print(f"\n  corr(DoD, cover)={np.corrcoef(D[mf],C[mf])[0,1]:+.2f}  "
      f"corr(DoD, slope)={np.corrcoef(D[mf],S[mf])[0,1]:+.2f}")
print(f"  partial corr(DoD, cover | slope)={partial(D[mf],C[mf],S[mf]):+.2f}  "
      f"partial corr(DoD, slope | cover)={partial(D[mf],S[mf],C[mf]):+.2f}")
print("\nRead: if DoD rises with COVER at fixed slope (rows vary across cols) but is flat "
      "with slope at fixed cover (cols constant down rows), and partial(DoD,cover|slope) >> "
      "partial(DoD,slope|cover), then forest DENSITY is the driver, as predicted.")
