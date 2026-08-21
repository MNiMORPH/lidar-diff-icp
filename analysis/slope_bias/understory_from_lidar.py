"""Lidar-derived understory/canopy structure from the gen2 (2021) cloud, and test
whether the tan(S) DoD bias tracks understory. Height above ground = z_return −
z_after(cell). Bands: ground |h|<0.15; LOW understory 0.15–1.0 m (the brush/forb layer
gen2's ground can sit on); brush 1.0–2.5 m; canopy >2.5 m. Classifies WHERE understory
is present and tests DoD vs understory at fixed slope (the mechanism check)."""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt as edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
X0,Y0,X1,Y1=577492.8,4882737.6,580032.8,4886237.6; res=5.0; ext=(X0,X1,Y0,Y1)
nx=int(round((X1-X0)/res)); ny=int(round((Y1-Y0)/res))
Z=np.load("data/derived/elba/z_after.npy"); dod=np.load("data/derived/elba/dod_refdatum.npy")
Zf=Z.copy(); nmz=~np.isfinite(Zf)
if nmz.any(): Zf=Zf[tuple(edt(nmz,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,res); slope=np.degrees(np.arctan(np.hypot(gx,gy)))
with laspy.open("data/after/3dep2021_fulltile.laz") as fh: p=fh.read()
x=np.array(p.x); y=np.array(p.y); z=np.array(p.z); cl=np.array(p.classification)
m=(x>=X0)&(x<X1)&(y>=Y0)&(y<Y1)&(cl!=7)
ix=((x[m]-X0)/res).astype(int); iy=((y[m]-Y0)/res).astype(int); cid=iy*nx+ix
h=z[m]-Zf[iy,ix]                                   # height above ground DEM
N=nx*ny
tot=np.bincount(cid,minlength=N).astype(float)
low =np.bincount(cid[(h>0.15)&(h<=1.0)],minlength=N).astype(float)   # forb/brush base
brush=np.bincount(cid[(h>1.0)&(h<=2.5)],minlength=N).astype(float)
can =np.bincount(cid[h>2.5],minlength=N).astype(float)
understory_frac=np.where(tot>0,(low+brush)/np.maximum(tot,1),np.nan).reshape(ny,nx)
lowveg_frac=np.where(tot>0,low/np.maximum(tot,1),np.nan).reshape(ny,nx)
np.save("data/derived/elba/understory_frac.npy",understory_frac)
fin=np.isfinite(dod)&np.isfinite(Z)&np.isfinite(understory_frac)
print(f"understory (0.15-2.5m) fraction: tile median {np.nanmedian(understory_frac[fin]):.2f}, "
      f"cells >20% understory: {100*np.mean(understory_frac[fin]>0.2):.0f}%")
print("\nDoD (mm) by slope, split by UNDERSTORY fraction (does the bias track understory?):")
print(f"  {'slope':8s} {'LOW und (<10%)':>16s} {'HIGH und (>30%)':>16s}")
for lo,hi in [(2,10),(10,20),(20,30),(30,90)]:
    band=fin&(slope>=lo)&(slope<hi)
    lu=band&(understory_frac<0.10); hu=band&(understory_frac>0.30)
    ls=f"{1000*np.median(dod[lu]):+.1f} (n={lu.sum()})" if lu.sum()>30 else f"-- ({lu.sum()})"
    hs=f"{1000*np.median(dod[hu]):+.1f} (n={hu.sum()})" if hu.sum()>30 else f"-- ({hu.sum()})"
    print(f"  {lo:2d}-{hi:2d} deg  {ls:>16s} {hs:>16s}")
# DoD vs understory at fixed slope (20-30) -- the cleanest mechanism check
band=fin&(slope>=20)&(slope<30)
print("\nDoD (mm) vs understory fraction, slope 20-30 deg (mechanism check):")
for a,b in [(0,.1),(.1,.2),(.2,.35),(.35,.6),(.6,1.01)]:
    mm=band&(understory_frac>=a)&(understory_frac<b)
    if mm.sum()>30: print(f"  understory {a:.2f}-{b:.2f}:  DoD {1000*np.median(dod[mm]):+6.1f} mm  (n={mm.sum()})")
# map
fig,ax=plt.subplots(1,2,figsize=(15,7))
im0=ax[0].imshow(understory_frac,extent=ext,origin="lower",cmap="YlGn",vmin=0,vmax=0.6)
ax[0].set_title("gen2 understory fraction (0.15-2.5 m returns / total)"); fig.colorbar(im0,ax=ax[0],shrink=.6)
im1=ax[1].imshow(np.where(fin,dod,np.nan),extent=ext,origin="lower",cmap="RdBu",vmin=-.3,vmax=.3)
ax[1].set_title("DoD (gen2-gen1)"); fig.colorbar(im1,ax=ax[1],shrink=.6)
fig.savefig("figures/understory_map.png",dpi=120,bbox_inches="tight"); print("\nwrote figures/understory_map.png")
