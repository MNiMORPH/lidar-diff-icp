#!/usr/bin/env python3
"""gen1 vs ITSELF: gen1's ground-return median sits lower in the core (densest) forest than
in the full forest.  What drives gen1's own downward shift — canopy density or slope?
Bin gen1 forest cells (per-cell median, from slope_normal_returns.npz) by penetration
(=ground-return fraction; lower=denser canopy) and by slope, incl. slope-controlled.

Frame caveat: d is measured against gen2 bare earth, so this is gen1 relative to gen2's
surface; a purely gen1 vs gen2-independent datum is not available in forest.  Slope is the
one geometry covariate independent of both classifiers.

    ./lidar-icp/bin/python analysis/ridgelines/gen1_sink_vs_density.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D="data/derived/elba_fulldensity/"
Z=np.load(D+"slope_normal_returns.npz")
g1=Z["gen1_ground_median_d"]*1000; n1=Z["gen1_n_ground"]; n2=Z["gen2_n_ground"]
pen=np.load(D+"penetration.npy"); slp=np.load(D+"slope.npy"); fld=np.load(D+"floodplain_mask.npy").astype(bool)
cf=np.load(D+"core_forest.npy")
forest=(pen<0.25)&~fld&np.isfinite(pen)&(n1>=8)&np.isfinite(g1)
p=pen[forest]; a=g1[forest]; s=slp[forest]

def binned(x,y,xe):
    mid=[];md=[]
    for i in range(len(xe)-1):
        m=(x>=xe[i])&(x<xe[i+1] if i<len(xe)-2 else x<=xe[i+1])
        if m.sum()<50: continue
        mid.append(np.median(x[m])); md.append(np.median(y[m]))
    return np.array(mid),np.array(md)

e=np.quantile(p,np.linspace(0,1,9)); mp,ma=binned(p,a,e)
sc=(s>=15)&(s<25); e2=np.quantile(p[sc],np.linspace(0,1,7)); mp2,ma2=binned(p[sc],a[sc],e2)
core_pen=np.median(pen[cf&forest]); core_med=np.median(g1[cf&forest])
full_med=np.median(a)

fig,ax=plt.subplots(figsize=(9,6))
ax.plot(mp,ma,"C0o-",lw=1.8,label="gen1 median vs canopy density (all forest)")
ax.plot(mp2,ma2,"C2s--",lw=1.6,label="gen1 median, slope-controlled (15-25°)")
ax.axvline(core_pen,color="0.4",ls=":",lw=1.2)
ax.plot([core_pen],[core_med],"r*",ms=16,label=f"core forest (pen {core_pen:.2f}, {core_med:+.0f} mm)")
ax.annotate("← denser canopy",(0.02,0.05),xycoords="axes fraction",fontsize=9)
ax.annotate("sparser canopy →",(0.70,0.05),xycoords="axes fraction",fontsize=9)
ax.set_xlabel("penetration = ground-return fraction  (LEFT = denser canopy)")
ax.set_ylabel("gen1 ground median, slope-normal d (mm)")
ax.set_title("gen1 vs itself: its own ground median sinks as canopy thickens (driver = density, not slope)")
ax.legend(fontsize=9); ax.grid(alpha=.3)
fig.savefig("figures/refdatum/gen1_sink_vs_density.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print(f"full forest gen1 med {full_med:+.0f} mm ; core forest {core_med:+.0f} mm")
print("wrote figures/refdatum/gen1_sink_vs_density.png")
