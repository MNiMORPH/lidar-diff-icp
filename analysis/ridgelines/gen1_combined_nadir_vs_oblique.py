#!/usr/bin/env python3
"""FOREST + FARMLAND combined, gen1 ground returns: relative elevation vs intensity, split by
scan geometry — NEAR-NADIR (|angle|<=2 deg) vs OBLIQUE (|angle|>=8 deg).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_combined_nadir_vs_oblique.py
"""
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
Zg=np.load("data/derived/elba_fulldensity/z_after.npy"); Zf=Zg.copy(); m=~np.isfinite(Zf)
if m.any(): Zf=Zf[tuple(distance_transform_edt(m,return_distances=False,return_indices=True))]
gy,gx=np.gradient(Zf,RES); cosd=1.0/np.sqrt(1.0+gx*gx+gy*gy)
Zff=Zf.ravel(); gxf=gx.ravel(); gyf=gy.ravel(); cosf=cosd.ravel()
g2pen=np.load("data/derived/elba_fulldensity/penetration.npy").ravel()
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool).ravel()
mask=((g2pen<0.25)|(g2pen>=0.45))&(~fld)&np.isfinite(g2pen)   # forest OR farmland (open)

Dn=[];In=[];Do=[];Io=[]
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(CHUNK):
        cl=np.asarray(pts.classification); ang=np.abs(np.asarray(pts.scan_angle_rank)).astype(float)
        base=cl==2
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64); z=np.asarray(pts.z,np.float64)
        inten=np.asarray(pts.intensity).astype(float)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        ok=base&(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)
        cell=np.where(ok,iy*NX+ix,0)
        cm=ok&mask[cell]
        for sel,DD,II in [(cm&(ang<=2.0),Dn,In),(cm&(ang>=8.0),Do,Io)]:
            if not sel.any(): continue
            c=cell[sel]; xc=X0+((c%NX)+0.5)*RES; yc=Y0+((c//NX)+0.5)*RES
            d=(z[sel]-(Zff[c]+gxf[c]*(x[sel]-xc)+gyf[c]*(y[sel]-yc)))*cosf[c]
            DD.append(d*1000); II.append(inten[sel])
dn=np.concatenate(Dn); iN=np.concatenate(In); do=np.concatenate(Do); iO=np.concatenate(Io)

def summarize(d,I,lbl):
    print(f"\n{lbl}: n={d.size:,}  median d={np.median(d):+.1f} mm  corr(d,I)={np.corrcoef(d,I)[0,1]:+.3f}")
    e=np.quantile(I,np.linspace(0,1,11)); mx=[];my=[]
    for i in range(len(e)-1):
        b=(I>=e[i])&(I<e[i+1] if i<len(e)-2 else I<=e[i+1])
        if b.sum()<200: continue
        mx.append(np.median(I[b])); my.append(np.median(d[b]))
        print(f"  intensity med {np.median(I[b]):4.0f}: d {np.median(d[b]):+7.1f} mm  n={b.sum()}")
    return mx,my
mxn,myn=summarize(dn,iN,"NEAR-NADIR (|angle|<=2)")
mxo,myo=summarize(do,iO,"OBLIQUE (|angle|>=8)")

fig,ax=plt.subplots(1,2,figsize=(15,6),sharey=True)
for a,(d,I,mx,my,lbl) in zip(ax,[(dn,iN,mxn,myn,"NEAR-NADIR (|angle|<=2deg)"),(do,iO,mxo,myo,"OBLIQUE (|angle|>=8deg)")]):
    hb=a.hexbin(I,d,gridsize=45,bins="log",cmap="viridis",mincnt=1,extent=(0,80,-400,300))
    a.plot(mx,my,"r.-",lw=2.5,ms=10,label="median d")
    a.set_xlim(0,80); a.set_ylim(-350,250); a.set_xlabel("ground-return intensity")
    a.set_title(f"{lbl}  (n={d.size:,}, med d={np.median(d):+.0f} mm, r={np.corrcoef(d,I)[0,1]:+.2f})"); a.legend(fontsize=8)
    fig.colorbar(hb,ax=a,label="log10 count")
ax[0].set_ylabel("relative ground elevation d (mm)")
fig.suptitle("gen1 FOREST+FARMLAND ground: elevation vs intensity — near-nadir vs oblique pulses",y=1.0)
fig.savefig("figures/refdatum/gen1_combined_nadir_vs_oblique.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_combined_nadir_vs_oblique.png")
