#!/usr/bin/env python3
"""Hillslope diffusion over ALL agricultural / ALL forested cells (not just crests), with
buffers around channels (fluvial transport, not diffusion) and around the forest/ag edge
(no cross-contamination). Full curvature range (concave hollows -> convex noses) pins K.

  dz/dt = K * Laplacian(z) + offset      offset: c_ag (const) / beta*cover (forest)

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/hillslope_regions_K.py
"""
import sys, numpy as np
sys.path.insert(0, "/home/awickert/dataanalysis/r.fluvial")
from rivernetworkx import dreich as Dr
Dr._RichDEMVersion = lambda: "local"
from scipy.ndimage import binary_erosion, binary_dilation, distance_transform_edt, gaussian_filter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lidar_diff_icp.viz import hillshade

RES = 5.0; X0, Y0 = 577492.8, 4882737.6; DT = 12.44
D = "data/derived/elba_fulldensity/"
z = np.load(D + "z_after.npy"); ny, nx = z.shape
zf = z.copy(); m = ~np.isfinite(zf)
if m.any(): zf = zf[tuple(distance_transform_edt(m, return_distances=False, return_indices=True))]
lap = np.load(D + "curv_laplacian.npy")               # d2z/dx2 + d2z/dy2 (1/m)
pen = np.load(D + "penetration.npy"); cover = np.load(D + "canopy_cover.npy")
dod = np.load("data/derived/elba_refdatum/dod_osm.npy"); rate = dod / DT

# --- channels: flow accumulation on the (smoothed) DEM, then buffer --------------------
filled = Dr.fill(gaussian_filter(zf, 1.5).astype(np.float64), nodata=-9999.0, cellsize=RES)
fi = Dr.build_flowinfo(filled, nodata=-9999.0, cellsize=RES); Dr.contributing_area(fi)
nc = np.zeros((ny, nx)); nc[fi["row_of"], fi["col_of"]] = fi["ncontrib"]
channel = nc > 1500                                   # ~3.75 ha contributing -> valley network
CHAN_BUF, EDGE_BUF = 10, 6                             # 50 m around channels, 30 m off cover edge
chan_zone = binary_dilation(channel, iterations=CHAN_BUF)

# --- land-cover cores (eroded off their mutual boundary; channel zone removed) ---------
ag = pen >= 0.45; forest = pen < 0.25
ag_core = binary_erosion(ag, iterations=EDGE_BUF) & ~chan_zone
for_core = binary_erosion(forest, iterations=EDGE_BUF) & ~chan_zone
fin = np.isfinite(dod) & np.isfinite(lap)
agc = ag_core & fin; forc = for_core & fin
print(f"ag cells={agc.sum()}  forest cells={forc.sum()}  (channel-zone {chan_zone.mean()*100:.0f}% of tile)")

def fit(y, X):
    c, *_ = np.linalg.lstsq(X, y, rcond=None); p = X @ c
    return c, 1 - np.sum((y-p)**2)/np.sum((y-np.mean(y))**2)

(K_ag, c_ag), r2a = fit(rate[agc], np.c_[lap[agc], np.ones(agc.sum())])
(K_f, beta), r2f = fit(rate[forc], np.c_[lap[forc], cover[forc]])
(K_fj, a_fj, b_fj), r2fj = fit(rate[forc], np.c_[lap[forc], np.ones(forc.sum()), cover[forc]])
print(f"\nFARMLAND: dz/dt = K*Lap + c   K_ag={K_ag:.4f} m2/yr  c_ag={c_ag*1000:+.2f} mm/yr  R2={r2a:.3f}")
print(f"FOREST  : dz/dt = K*Lap + b*cover   K_for={K_f:.4f} m2/yr  beta={beta*1000:+.2f} mm/yr/cover  R2={r2f:.3f}")
print(f"FOREST  : dz/dt = K*Lap + a + b*cover   K_for={K_fj:.4f}  a={a_fj*1000:+.2f}  b={b_fj*1000:+.2f}  R2={r2fj:.3f}")

# --- plots: region map + binned dz/dt vs curvature ------------------------------------
ext = (X0, X0+nx*RES, Y0, Y0+ny*RES); hs = hillshade(zf, RES, X0, Y0, fill_gaps=True)
fig, ax = plt.subplots(figsize=(10, 11)); ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
ov = np.zeros((ny, nx, 4)); ov[ag_core] = (0.9, 0.8, 0.1, 0.5); ov[for_core] = (0.1, 0.6, 0.1, 0.5)
ov[chan_zone] = (0.2, 0.4, 0.9, 0.25)
ax.imshow(ov, extent=ext, origin="lower")
ax.set_title("hillslope-diffusion regions: farmland (yellow), forest (green),\n"
             "channel buffer (blue, excluded); forest/ag edges buffered off")
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
fig.savefig("figures/refdatum/hillslope_regions.png", dpi=130, bbox_inches="tight"); plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 6))
for tag, mask, K, c, col in [("farmland", agc, K_ag, c_ag, "goldenrod"),
                             ("forest", forc, K_fj, a_fj+b_fj*np.median(cover[forc]), "green")]:
    L = lap[mask]; y = rate[mask]*1000
    edges = np.quantile(L, np.linspace(0, 1, 16)); mid = 0.5*(edges[:-1]+edges[1:])
    med = [np.median(y[(L>=edges[i])&(L<edges[i+1])]) for i in range(len(edges)-1)]
    ax.plot(mid, med, "o-", color=col, label=f"{tag}: K={K:.3f} m²/yr")
    ax.plot(mid, (K*mid + c)*1000, "--", color=col, alpha=0.5)
ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
ax.set_xlabel("Laplacian curvature ∇²z (1/m)  [convex<0  concave>0]")
ax.set_ylabel("dz/dt (mm/yr)  [erosion<0  deposition>0]")
ax.set_title("dz/dt vs curvature — hillslope diffusion (binned medians)"); ax.legend(); ax.grid(alpha=0.3)
fig.savefig("figures/refdatum/dzdt_vs_curvature.png", dpi=130, bbox_inches="tight"); plt.close(fig)
np.save(D + "ag_region.npy", ag_core); np.save(D + "forest_region.npy", for_core)
print("wrote figures/refdatum/hillslope_regions.png, dzdt_vs_curvature.png")
