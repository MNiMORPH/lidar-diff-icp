#!/usr/bin/env python3
"""Compute per-return INCIDENCE ANGLE to the local surface normal (beam vs ground), so we can
select "near-surface-perpendicular" pulses instead of "near-nadir". On flat ground incidence ==
|scan angle|; on a slope it depends on scan angle, slope, and the beam azimuth relative to aspect.

Beam azimuth is reconstructed per flight line (point_source_id): heading from (x,y) vs gps_time;
the cross-track side of +scan_angle from corr(cross-track pos, scan_angle) (|corr|~0.99 here).
Surface normal from the gen2 bare-earth (z_after) gradient.

Validation (printed): on flat farmland incidence must reduce to |scan angle|.
Then: redo the all-forest vs core-forest floor histogram using NEAR-PERPENDICULAR (incidence<=2°)
CSF cloth ground, to compare against the near-nadir version.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/incidence_angle.py
"""
import numpy as np, laspy, math
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
CSF="data/csf_cache/elba.las"
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES)                     # gy=d/dNorth, gx=d/dEast (m/m)
slope=np.degrees(np.arctan(np.hypot(gx,gy)))
gxf=gx.ravel(); gyf=gy.ravel(); nnorm=np.sqrt(gxf**2+gyf**2+1.0)
g2pen=np.load("data/derived/elba_fulldensity/penetration.npy").ravel()
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool).ravel()
core=np.load("data/derived/elba_fulldensity/core_forest.npy").ravel()

las=laspy.read(CSF)
x=np.asarray(las.x,np.float64); y=np.asarray(las.y,np.float64); z=np.asarray(las.z,np.float64)
sa=np.asarray(las.scan_angle).astype(float)*0.006          # signed degrees
psid=np.asarray(las.point_source_id); gt=np.asarray(las.gps_time)
ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
ing=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)
cell=np.where(ing,iy*NX+ix,0)

# --- per flight line: heading + cross-track side sign ---
bhx=np.zeros(len(x)); bhy=np.zeros(len(x))   # beam horizontal unit (ground->sensor)
for pl in np.unique(psid):
    m=psid==pl
    vx=np.polyfit(gt[m],x[m],1)[0]; vy=np.polyfit(gt[m],y[m],1)[0]
    H=math.atan2(vy,vx); cxu=np.array([-math.sin(H),math.cos(H)])   # cross-track unit (left of heading)
    cross=(x[m]-x[m].mean())*cxu[0]+(y[m]-y[m].mean())*cxu[1]
    sgn=np.sign(np.corrcoef(cross,sa[m])[0,1])                      # +scan_angle on sgn*cxu side
    # ground point offset dir = sign(sa)*sgn*cxu ; beam points OPPOSITE (toward sensor)
    hx=-np.sign(sa[m])*sgn*cxu[0]; hy=-np.sign(sa[m])*sgn*cxu[1]
    bhx[m]=hx; bhy[m]=hy
    print(f"  line {pl}: heading {math.degrees(H):+.0f}deg  corr(cross,scan)={np.corrcoef(cross,sa[m])[0,1]:+.2f}  n={m.sum():,}")

th=np.radians(np.abs(sa))
bx=np.sin(th)*bhx; by=np.sin(th)*bhy; bz=np.cos(th)          # beam unit (ground->sensor)
# surface normal (upward) = (-gx,-gy,1)/norm ; incidence cos = b . n
nn=nnorm[cell]
cosi=(bx*(-gxf[cell])+by*(-gyf[cell])+bz*1.0)/nn
inc=np.degrees(np.arccos(np.clip(cosi,-1,1)))

# --- VALIDATION: flat farmland incidence should equal |scan angle| ---
flat=ing&(g2pen[cell]>=0.45)&(~fld[cell])&(slope.ravel()[cell]<2.0)
print(f"\nVALIDATION on flat farmland (slope<2deg), incidence vs |scan angle|:")
for lo,hi in [(0,2),(4,6),(8,10),(12,16)]:
    b=flat&(np.abs(sa)>=lo)&(np.abs(sa)<hi)
    if b.sum()<500: continue
    print(f"  |scan| {lo}-{hi}deg: median incidence {np.median(inc[b]):.1f}deg  (should ~ {(lo+hi)/2})  n={b.sum():,}")
steep=ing&(g2pen[cell]<0.25)&(~fld[cell])&(slope.ravel()[cell]>20)
print(f"\nOn STEEP forest (slope>20deg): median |scan angle| {np.median(np.abs(sa)[steep]):.1f}deg  "
      f"vs median incidence {np.median(inc[steep]):.1f}deg  (incidence != scan angle on slopes)  n={steep.sum():,}")

# --- forest floor histogram: NEAR-PERPENDICULAR (incidence<=2) all vs core forest ---
def dof(sel):
    c=cell[sel]; xc=X0+((c%NX)+0.5)*RES; yc=Y0+((c//NX)+0.5)*RES
    return (z[sel]-(Zf.ravel()[c]+gxf[c]*(x[sel]-xc)+gyf[c]*(y[sel]-yc)))*(1.0/nnorm[c])*1000
allf=(g2pen<0.25)&(~fld)&np.isfinite(g2pen)
perp=ing&(inc<=2.0)
da=dof(perp&allf[cell]); dc=dof(perp&core[cell])
print(f"\nNEAR-PERPENDICULAR (incidence<=2deg) CSF cloth forest floor:")
print(f"  all forest median {np.median(da):+.1f} (n={da.size:,})  core {np.median(dc):+.1f} (n={dc.size:,})  shift {np.median(dc)-np.median(da):+.1f} mm")

fig,ax=plt.subplots(figsize=(10,6)); b=np.arange(-350,251,10)
for d,col,l in [(da,"C2","ALL forest"),(dc,"C1","CORE forest")]:
    ax.hist(d,bins=b,density=True,histtype="step",lw=2,color=col,label=f"{l} (n={d.size:,}, med {np.median(d):+.0f} mm)")
    ax.axvline(np.median(d),color=col,ls=":",lw=1.5)
ax.set_xlabel("forest-floor elevation d (mm)"); ax.set_ylabel("density")
ax.set_title("NEAR-SURFACE-PERPENDICULAR (incidence<=2°) CSF cloth ground: all vs core forest")
ax.legend(); ax.grid(alpha=.3)
fig.savefig("figures/refdatum/incidence_perp_all_vs_core.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/incidence_perp_all_vs_core.png")
