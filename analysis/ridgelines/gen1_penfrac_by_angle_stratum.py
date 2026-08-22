#!/usr/bin/env python3
"""Ground-return fraction (penetration = class-2 / non-noise) split by scan geometry
(NEAR-NADIR |angle|<=2 vs OBLIQUE |angle|>=8), for three strata: CORE FARMLAND, FOREST (all),
CORE FOREST. Shows whether oblique pulses under-penetrate dense canopy (the mechanism behind
the all->core forest-floor sink). Full gen1 cloud (has all returns).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/gen1_penfrac_by_angle_stratum.py
"""
import numpy as np, laspy
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NY,NX=700,508; X0,Y0=577492.8,4882737.6; RES=5.0
GEN1="data/before/4342-29-64.laz"; CHUNK=5_000_000
g2pen=np.load("data/derived/elba_fulldensity/penetration.npy").ravel()
fld=np.load("data/derived/elba_fulldensity/floodplain_mask.npy").astype(bool).ravel()
strata={"core farmland":np.load("data/derived/elba_fulldensity/core_open.npy").ravel(),
        "forest (all)":(g2pen<0.25)&(~fld)&np.isfinite(g2pen),
        "core forest":np.load("data/derived/elba_fulldensity/core_forest.npy").ravel()}

# counts[stratum][angleclass] = [total, ground]
C={s:{"nadir":np.zeros(2),"oblique":np.zeros(2)} for s in strata}
with laspy.open(GEN1) as f:
    for pts in f.chunk_iterator(CHUNK):
        cl=np.asarray(pts.classification); aa=np.abs(np.asarray(pts.scan_angle_rank)).astype(float)
        x=np.asarray(pts.x,np.float64); y=np.asarray(pts.y,np.float64)
        ix=((x-X0)/RES).astype(np.int64); iy=((y-Y0)/RES).astype(np.int64)
        ok=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(cl!=7); cell=np.where(ok,iy*NX+ix,0)
        g=(cl==2)
        for s,mask in strata.items():
            ins=ok&mask[cell]
            for acl,asel in [("nadir",aa<=2.0),("oblique",aa>=8.0)]:
                sel=ins&asel
                C[s][acl][0]+=sel.sum(); C[s][acl][1]+=(sel&g).sum()

print(f"{'stratum':16s} {'nadir pen':>12s} {'oblique pen':>12s}   (ground/total; low=under-penetration)")
res={}
for s in strata:
    pn=C[s]["nadir"][1]/max(C[s]["nadir"][0],1); po=C[s]["oblique"][1]/max(C[s]["oblique"][0],1)
    res[s]=(pn,po)
    print(f"{s:16s} {pn:>10.3f}   {po:>10.3f}   nadir/oblique ratio {pn/max(po,1e-6):.2f}   "
          f"(n nadir {int(C[s]['nadir'][0]):,}, oblique {int(C[s]['oblique'][0]):,})")

fig,ax=plt.subplots(figsize=(9,6))
xs=np.arange(len(strata)); w=0.35
ax.bar(xs-w/2,[res[s][0] for s in strata],w,color="C0",label="near-nadir (|angle|<=2°)")
ax.bar(xs+w/2,[res[s][1] for s in strata],w,color="C3",label="oblique (|angle|>=8°)")
for i,s in enumerate(strata):
    ax.text(i-w/2,res[s][0]+.01,f"{res[s][0]:.2f}",ha="center",fontsize=9)
    ax.text(i+w/2,res[s][1]+.01,f"{res[s][1]:.2f}",ha="center",fontsize=9)
ax.set_xticks(xs); ax.set_xticklabels(list(strata)); ax.set_ylabel("ground-return fraction (penetration)")
ax.set_title("gen1 ground-return fraction: near-nadir vs oblique, by stratum")
ax.legend(); ax.grid(alpha=.3,axis="y")
fig.savefig("figures/refdatum/gen1_penfrac_by_angle_stratum.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/gen1_penfrac_by_angle_stratum.png")
