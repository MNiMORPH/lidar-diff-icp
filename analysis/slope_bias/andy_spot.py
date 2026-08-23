#!/usr/bin/env python3
"""Zoom on Andy's ground-known flat blufftop forest: 44.1144712, -92.0109470
(UTM15N E=579143.9 N=4885062.5; AOI grid row 464 col 330). Report DoD/slope/canopy in a
window and whether it matches the clean-forest finding (no canopy offset at low slope).
"""
import numpy as np, laspy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import distance_transform_edt
from lidar_diff_icp.viz import hillshade

RES = 5.0; X0, Y0 = 577492.8, 4882737.6
E, N = 579143.9, 4885062.5
z = np.load("data/derived/elba_fulldensity/z_after.npy")
slope = np.load("data/derived/elba_fulldensity/slope.npy")
dod = np.load("data/derived/elba_refdatum/dod_geoid.npy")
margin = np.load("data/derived/elba_fulldensity/blufftop_margin_mask.npy")
ny, nx = z.shape
cc, cr = int((E-X0)/RES), int((N-Y0)/RES)
HW = 40                                                 # +/- 200 m window
r0, r1, c0, c1 = cr-HW, cr+HW, cc-HW, cc+HW
Eb0, Eb1 = X0+c0*RES, X0+c1*RES; Nb0, Nb1 = Y0+r0*RES, Y0+r1*RES

# penetration in the window (crop cloud)
g = np.zeros(ny*nx, np.int64); t = np.zeros(ny*nx, np.int64)
with laspy.open("data/after/3dep2021_fulldensity.laz") as fh:
    for pts in fh.chunk_iterator(20_000_000):
        x = np.asarray(pts.x); y = np.asarray(pts.y); cl = np.asarray(pts.classification)
        keep = (cl!=7)&(x>=Eb0)&(x<Eb1)&(y>=Nb0)&(y<Nb1)
        if not keep.any(): continue
        ix=((x[keep]-X0)/RES).astype(np.int64); iy=((y[keep]-Y0)/RES).astype(np.int64)
        c=iy*nx+ix; t+=np.bincount(c,minlength=ny*nx); g+=np.bincount(c[cl[keep]==2],minlength=ny*nx)
pen = np.divide(g, np.maximum(t,1)).reshape(ny, nx)

W = (slice(r0,r1), slice(c0,c1))
dw, sw, pw, mw = dod[W], slope[W], pen[W], margin[W]
fin = np.isfinite(dw); forest = pw < 0.25
print(f"window {2*HW*RES:.0f} m around Andy's spot (row {cr} col {cc})")
print(f"  cell AT the point: DoD={dod[cr,cc]*1000:+.1f} mm  slope={slope[cr,cc]:.1f} deg  pen={pen[cr,cc]:.2f}")
print(f"  window forest cells: {int((forest&fin).sum())}  median slope={np.nanmedian(sw[forest&fin]):.1f} deg  "
      f"median pen={np.nanmedian(pw[forest&fin]):.2f}")
print(f"  in my blufftop-margin selection: {int((mw).sum())} cells")
print("\n  median DoD (mm) by slope, forest cells in this window:")
for lo,hi in [(0,3),(3,6),(6,10),(10,15),(15,90)]:
    m = forest & fin & (sw>=lo) & (sw<hi)
    if m.any(): print(f"    {lo:>2}-{hi:<2} deg: n={m.sum():>4}  medDoD={np.median(dw[m]*1000):+6.1f} mm")
mm = forest & fin & (sw<6)
if mm.any(): print(f"  FLAT forest (<6 deg) here: n={mm.sum()}  medDoD={np.median(dw[mm]*1000):+.1f} mm")

# zoom figure: hillshade + DoD, with the point and selection
zf = z.copy(); nm=~np.isfinite(zf)
if nm.any(): zf = zf[tuple(distance_transform_edt(nm,return_distances=False,return_indices=True))]
hs = hillshade(zf, RES, X0, Y0, fill_gaps=True)[W]
ext = (Eb0, Eb1, Nb0, Nb1)
fig, ax = plt.subplots(1,2, figsize=(16,8))
for a in ax: a.imshow(hs, extent=ext, origin="lower", cmap="gray")
im = ax[0].imshow(np.where(fin,dw,np.nan), extent=ext, origin="lower", cmap="RdBu", vmin=-0.1, vmax=0.1, alpha=0.7)
ax[0].set_title("DoD gen2-gen1 (m)"); fig.colorbar(im, ax=ax[0], shrink=0.6)
ov=np.zeros((r1-r0,c1-c0,4)); ov[mw]=(0.9,0.1,0.1,0.8); ov[forest&~mw]=(0.1,0.5,0.1,0.2)
ax[1].imshow(ov, extent=ext, origin="lower"); ax[1].set_title("forest (green), my selection (red)")
for a in ax:
    a.plot(E, N, "c*", ms=20, mec="k"); a.set_xlabel("Easting (m)"); a.set_ylabel("Northing (m)")
fig.suptitle("Andy's flat blufftop forest spot (cyan star): 44.11447, -92.01095")
import os; os.makedirs("figures", exist_ok=True)
fig.savefig("figures/andy_spot.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/andy_spot.png")
