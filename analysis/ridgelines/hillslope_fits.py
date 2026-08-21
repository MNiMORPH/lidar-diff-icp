#!/usr/bin/env python3
"""Hillslope-diffusion + lidar-offset fits on the CLEAN OSM-datum DoD, on ridge crests
(no incoming material). Curvature = d2z/dx2 (per Andy's equation); dt from flight dates.

  FARMLAND (open crests):  dz/dt = K_ag   * d2z/dx2 + c_ag            (c_ag = error term)
  FOREST   (forest crests): dz/dt = K_for * d2z/dx2 + f(cover)        (f = lidar offset)

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/hillslope_fits.py
"""
import numpy as np
DT = 12.44                                              # yr, 2008-11-21 -> 2021-05-01
D = "data/derived/elba_fulldensity/"
R = np.load(D + "ridgecrest_pixels.npz", allow_pickle=True)
row, col = R["row"].astype(int), R["col"].astype(int)
kap = R["curv_xx"].astype(float)                        # d2z/dx2 (1/m; convex-up negative)
lap = R["curv_laplacian"].astype(float)
slope = R["slope_deg"].astype(float); pen = R["penetration"].astype(float)
dod = np.load("data/derived/elba_refdatum/dod_osm.npy")[row, col]
cover = np.load(D + "canopy_cover.npy")[row, col]
elev = np.load(D + "z_after.npy")[row, col]
rate = dod / DT                                         # m/yr
fin = np.isfinite(dod) & np.isfinite(kap)

def fit(y, X):                                          # least squares, return coefs + R2
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ c; ss = 1 - np.sum((y-pred)**2)/np.sum((y-np.mean(y))**2)
    return c, ss

print(f"dt = {DT} yr;  curvature = d2z/dx2 (curv_xx)\n")
# ---- FARMLAND: dz/dt = K*kap + c -------------------------------------------------------
ag = fin & (pen >= 0.45)
(K_ag, c_ag), r2 = fit(rate[ag], np.c_[kap[ag], np.ones(ag.sum())])
print(f"FARMLAND (open crests, n={ag.sum()}):")
print(f"  dz/dt = K*d2z/dx2 + c   ->  K_ag = {K_ag:.4f} m^2/yr   c_ag(error) = {c_ag*1000:+.2f} mm/yr   R2={r2:.3f}")
# also with Laplacian
(K_agL, c_agL), _ = fit(rate[ag], np.c_[lap[ag], np.ones(ag.sum())])
print(f"  [Laplacian form]        K_ag = {K_agL:.4f} m^2/yr   c_ag = {c_agL*1000:+.2f} mm/yr")

# ---- FOREST: dz/dt = K*kap + f(cover) --------------------------------------------------
fo = fin & (pen < 0.25)
# (a) linear offset f(cover) = a + b*cover, shared K
(K_f, a_f, b_f), r2f = fit(rate[fo], np.c_[kap[fo], np.ones(fo.sum()), cover[fo]])
print(f"\nFOREST (forest crests, n={fo.sum()}):")
print(f"  dz/dt = K*d2z/dx2 + (a + b*cover)")
print(f"  K_for = {K_f:.4f} m^2/yr   a = {a_f*1000:+.2f} mm/yr   b = {b_f*1000:+.2f} mm/yr per unit cover   R2={r2f:.3f}")
print(f"  -> f(cover) offset at cover=0.5: {(a_f+b_f*0.5)*1000:+.1f}  at 0.9: {(a_f+b_f*0.9)*1000:+.1f} mm/yr")
# (b) shared K + per-cover-bin offset (nonparametric f)
edges = np.quantile(cover[fo], [0,.2,.4,.6,.8,1.0]); nb = 5
binid = np.clip(np.digitize(cover[fo], edges[1:-1]), 0, nb-1)
X = np.zeros((fo.sum(), 1+nb)); X[:,0] = kap[fo]
for b in range(nb): X[binid==b, 1+b] = 1.0
c, r2b = fit(rate[fo], X)
print(f"  [shared K + per-cover-bin f]  K_for = {c[0]:.4f} m^2/yr   R2={r2b:.3f}")
for b in range(nb):
    cc = cover[fo][binid==b]
    print(f"    cover {edges[b]:.2f}-{edges[b+1]:.2f} (med {np.median(cc):.2f}): f = {c[1+b]*1000:+.1f} mm/yr  (n={ (binid==b).sum() })")

# ---- K_forest vs elevation (dolostone caprock check; single K reported above) ----------
print("\nK_for by elevation band (dolostone note; curvature-only slope of dz/dt vs kap):")
for lo,hi in [(210,260),(260,300),(300,360)]:
    m = fo & (elev>=lo) & (elev<hi)
    if m.sum()>50:
        (k,_a,_b),_ = fit(rate[m], np.c_[kap[m], np.ones(m.sum()), cover[m]])
        print(f"  elev {lo}-{hi} m: K_for = {k:.4f} m^2/yr  (n={m.sum()})")

print(f"\nCompare: K_ag = {K_ag:.4f}  vs  K_for = {K_f:.4f} m^2/yr")
