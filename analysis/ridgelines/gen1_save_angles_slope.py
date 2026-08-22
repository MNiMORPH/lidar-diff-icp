#!/usr/bin/env python3
"""Save per-return ANGLE information for gen1 CSF cloth ground returns (incidence to local
surface normal, signed scan angle, local slope, forest-floor elevation d, cell, flight line,
stratum), and plot the SLOPE DEPENDENCY of the forest-floor elevation.

Incidence reconstruction validated in incidence_angle.py (flat ground -> |scan angle|).
Saves data/derived/elba_fulldensity/gen1_csf_angles.npz.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_save_angles_slope.py
"""
import numpy as np, laspy, math
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
CSF="data/csf_cache/elba.las"
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES); slope_deg=np.degrees(np.arctan(np.hypot(gx,gy)))
gxf=gx.ravel(); gyf=gy.ravel(); nnorm=np.sqrt(gxf**2+gyf**2+1.0); Zflat=Zf.ravel()
g2pen=np.load("data/derived/elba_fulldensity/penetration.npy").ravel()
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool).ravel()
core=np.load("data/derived/elba_fulldensity/core_forest.npy").ravel()
copen=np.load("data/derived/elba_fulldensity/core_open.npy").ravel()

las=laspy.read(CSF)
x=np.asarray(las.x,np.float64); y=np.asarray(las.y,np.float64); z=np.asarray(las.z,np.float64)
sa=np.asarray(las.scan_angle).astype(float)*0.006
psid=np.asarray(las.point_source_id); gt=np.asarray(las.gps_time)
ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
ing=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY); cell=np.where(ing,iy*NX+ix,0)

bhx=np.zeros(len(x)); bhy=np.zeros(len(x))
for pl in np.unique(psid):
    m=psid==pl
    vx=np.polyfit(gt[m],x[m],1)[0]; vy=np.polyfit(gt[m],y[m],1)[0]; H=math.atan2(vy,vx)
    cxu=np.array([-math.sin(H),math.cos(H)])
    cross=(x[m]-x[m].mean())*cxu[0]+(y[m]-y[m].mean())*cxu[1]
    sgn=np.sign(np.corrcoef(cross,sa[m])[0,1])
    bhx[m]=-np.sign(sa[m])*sgn*cxu[0]; bhy[m]=-np.sign(sa[m])*sgn*cxu[1]
th=np.radians(np.abs(sa)); bx=np.sin(th)*bhx; by=np.sin(th)*bhy; bz=np.cos(th)
inc=np.degrees(np.arccos(np.clip((bx*(-gxf[cell])+by*(-gyf[cell])+bz)/nnorm[cell],-1,1)))
slp=slope_deg.ravel()[cell]
xc=X0+((cell%NX)+0.5)*RES; yc=Y0+((cell//NX)+0.5)*RES
d=(z-(Zflat[cell]+gxf[cell]*(x-xc)+gyf[cell]*(y-yc)))*(1.0/nnorm[cell])*1000  # mm

# stratum code: 1 forest, 2 farmland(open), 0 other ; plus core flags
strat=np.zeros(len(x),np.int8)
strat[ing&((g2pen[cell]<0.25)&~fld[cell])]=1
strat[ing&((g2pen[cell]>=0.45)&~fld[cell])]=2
np.savez_compressed("data/derived/elba_fulldensity/gen1_csf_angles.npz",
    incidence=inc.astype(np.float32), scan_angle=sa.astype(np.float32), slope=slp.astype(np.float32),
    d_mm=d.astype(np.float32), cell=cell.astype(np.int32), point_source_id=psid.astype(np.int32),
    stratum=strat, core_forest=core[cell]&ing, core_open=copen[cell]&ing, in_grid=ing)
print("saved data/derived/elba_fulldensity/gen1_csf_angles.npz  (n=%d returns)"%len(x))

# --- SLOPE DEPENDENCY of forest-floor elevation (forest, CSF ground) ---
F=ing&(strat==1); O=ing&(strat==2)
def bin_med(mask,xv,nb=10,lo=0,hi=40):
    e=np.linspace(lo,hi,nb+1); mx=[];my=[];mc=[]
    for i in range(nb):
        b=mask&(xv>=e[i])&(xv<e[i+1])
        if b.sum()<300: continue
        mx.append((e[i]+e[i+1])/2); my.append(np.median(d[b])); mc.append(b.sum())
    return np.array(mx),np.array(my),mc
print("\nforest-floor d vs SLOPE (forest CSF ground):")
sx,sy,sc=bin_med(F,slp)
for a,b,c in zip(sx,sy,sc): print(f"  slope {a:4.0f}deg: d {b:+7.1f} mm  n={c:,}")
print(f"corr(d, slope) forest = {np.corrcoef(d[F],slp[F])[0,1]:+.3f}")

fig,ax=plt.subplots(1,2,figsize=(14,6))
sxo,syo,_=bin_med(O,slp,nb=8,lo=0,hi=20)
ax[0].plot(sx,sy,"C0o-",label="forest"); ax[0].plot(sxo,syo,"C2s-",label="farmland")
ax[0].set_xlabel("slope (deg)"); ax[0].set_ylabel("forest-floor elevation d (mm)")
ax[0].set_title("floor elevation vs SLOPE"); ax[0].legend(); ax[0].grid(alpha=.3)
ix2,iy2,_=bin_med(F,inc,nb=12,lo=0,hi=45)
ax[1].plot(ix2,iy2,"C3o-",label="forest")
ax[1].set_xlabel("incidence angle to surface (deg)"); ax[1].set_ylabel("forest-floor elevation d (mm)")
ax[1].set_title("floor elevation vs INCIDENCE (beam vs surface)"); ax[1].legend(); ax[1].grid(alpha=.3)
fig.suptitle("gen1 CSF forest-floor elevation: slope and incidence dependency",y=1.0)
fig.savefig("figures/refdatum/gen1_floor_vs_slope_incidence.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_floor_vs_slope_incidence.png")
