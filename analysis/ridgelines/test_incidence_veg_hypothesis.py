#!/usr/bin/env python3
"""Test Andy's hypothesis (per-return, gen1 CSF ground, from gen1_csf_angles.npz):
 P-B: at FIXED slope, floor elevation d INCREASES with incidence angle (slope-perpendicular
      beams read deeper; slope-oblique read higher because longer path through near-ground veg).
 P-A: at fixed slope, more-vegetated cells retain lower-incidence (more perpendicular) ground
      returns (oblique beams blocked by canopy) -> median ground-return incidence drops with veg.
Veg proxy: canopy_height_p95 (gen2 canopy_struct) -- a height, less scan-geometry-confounded
than a return fraction.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/test_incidence_veg_hypothesis.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D="data/derived/elba_fulldensity/"
Z=np.load(D+"gen1_csf_angles.npz")
inc=Z["incidence"]; slp=Z["slope"]; d=Z["d_mm"]; cell=Z["cell"]; strat=Z["stratum"]
p95=np.load(D+"canopy_struct.npz")["canopy_height_p95"].ravel()
veg=np.load(D+"canopy_struct.npz")["understory_frac"].ravel()[cell]  # near-ground veg (0.5-2m)
F=(strat==1)&np.isfinite(inc)&np.isfinite(d)   # forest

# ---- P-B: floor d vs incidence, in slope bands ----
print("P-B: floor d (mm) vs incidence, by slope band (forest):")
for slo,shi in [(8,16),(16,24),(24,32)]:
    m=F&(slp>=slo)&(slp<shi)
    print(f"  slope {slo}-{shi}deg (n={m.sum():,}):")
    e=np.array([0,5,10,15,20,25,30,40,60])
    xs=[];ys=[]
    for i in range(len(e)-1):
        b=m&(inc>=e[i])&(inc<e[i+1])
        if b.sum()<500: continue
        xs.append((e[i]+e[i+1])/2); ys.append(np.median(d[b]))
        print(f"     inc {e[i]:2d}-{e[i+1]:2d}: d {np.median(d[b]):+7.1f}  n={b.sum():,}")
    if len(xs)>2: print(f"     -> corr(d,inc | this slope band) = {np.corrcoef(d[m],inc[m])[0,1]:+.3f}")

# ---- P-A: median ground-return incidence vs veg, fixed slope band 16-24 ----
print("\nP-A: median ground-return incidence vs vegetation (slope 16-24deg forest):")
m=F&(slp>=16)&(slp<24)&np.isfinite(veg)
vv=veg[m]; ii=inc[m]; dd=d[m]
e=np.quantile(vv,np.linspace(0,1,7))
for i in range(len(e)-1):
    b=(vv>=e[i])&(vv<e[i+1] if i<len(e)-2 else vv<=e[i+1])
    if b.sum()<500: continue
    print(f"  understory {e[i]:4.2f}-{e[i+1]:4.2f}: median incidence {np.median(ii[b]):4.1f}deg  floor d {np.median(dd[b]):+6.1f} mm  n={b.sum():,}")
print(f"  corr(incidence, veg | slope 16-24) = {np.corrcoef(ii,vv)[0,1]:+.3f}   (P-A predicts NEGATIVE)")
print(f"  corr(floor d, veg | slope 16-24)   = {np.corrcoef(dd,vv)[0,1]:+.3f}")

# ---- figure ----
fig,ax=plt.subplots(1,2,figsize=(14,6))
for slo,shi,c in [(8,16,"C0"),(16,24,"C1"),(24,32,"C3")]:
    m=F&(slp>=slo)&(slp<shi); e=np.array([0,5,10,15,20,25,30,40,60]); xs=[];ys=[]
    for i in range(len(e)-1):
        b=m&(inc>=e[i])&(inc<e[i+1])
        if b.sum()<500: continue
        xs.append((e[i]+e[i+1])/2); ys.append(np.median(d[b]))
    ax[0].plot(xs,ys,"o-",color=c,label=f"slope {slo}-{shi}°")
ax[0].set_xlabel("incidence angle to surface (deg)"); ax[0].set_ylabel("floor elevation d (mm)")
ax[0].set_title("P-B: floor vs incidence, by slope band"); ax[0].legend(); ax[0].grid(alpha=.3)
m=F&(slp>=16)&(slp<24)&np.isfinite(veg); vv=veg[m]; ii=inc[m]
e=np.quantile(vv,np.linspace(0,1,8)); xs=[];ys=[]
for i in range(len(e)-1):
    b=(vv>=e[i])&(vv<e[i+1] if i<len(e)-2 else vv<=e[i+1])
    if b.sum()<500: continue
    xs.append(np.median(vv[b])); ys.append(np.median(ii[b]))
ax[1].plot(xs,ys,"C2o-"); ax[1].set_xlabel("understory fraction (0.5-2m)"); ax[1].set_ylabel("median ground-return incidence (deg)")
ax[1].set_title("P-A: ground-return incidence vs veg (slope 16-24°)"); ax[1].grid(alpha=.3)
fig.suptitle("Testing incidence/veg hypothesis for gen1 forest floor",y=1.0)
fig.savefig("figures/refdatum/test_incidence_veg.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/test_incidence_veg.png")
