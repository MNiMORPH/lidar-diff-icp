#!/usr/bin/env python3
"""Core forest: PDF of near-nadir ground-return elevations vs PDF of ALL ground-return
elevations, as a function of slope (one panel per slope band).
near-nadir = |scan angle| <= 2 deg. Elevation = d_mm (slope-normal to gen2 bare earth; contains
the ~67 mm constant datum offset). From gen1_csf_angles.npz (gen1 CSF cloth ground).

    ./lidar-icp/bin/python analysis/ridgelines/nadir_vs_all_pdf_by_slope.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Z=np.load("data/derived/elba_fulldensity/gen1_csf_angles.npz")
d=Z["d_mm"]; sa=Z["scan_angle"]; slp=Z["slope"]; cf=Z["core_forest"]; ing=Z["in_grid"]
base=cf&ing&np.isfinite(d)&np.isfinite(sa)&np.isfinite(slp)
nadir=base&(np.abs(sa)<=2.0)

bands=[(0,8),(8,16),(16,24),(24,32),(32,45)]
bins=np.arange(-350,251,10)
fig,axes=plt.subplots(1,len(bands),figsize=(2.9*len(bands),4.6),sharey=True)
for ax,(lo,hi) in zip(axes,bands):
    ma=base&(slp>=lo)&(slp<hi); mn=nadir&(slp>=lo)&(slp<hi)
    if ma.sum()>50:
        ax.hist(d[ma],bins=bins,density=True,histtype="step",lw=2,color="0.4",
                label=f"all (n={ma.sum():,}, med {np.median(d[ma]):+.0f})")
        ax.axvline(np.median(d[ma]),color="0.4",ls=":",lw=1)
    if mn.sum()>50:
        ax.hist(d[mn],bins=bins,density=True,histtype="step",lw=2,color="C0",
                label=f"near-nadir (n={mn.sum():,}, med {np.median(d[mn]):+.0f})")
        ax.axvline(np.median(d[mn]),color="C0",ls=":",lw=1)
    ax.set_title(f"slope {lo}-{hi}°"); ax.set_xlabel("ground-return elevation d (mm)")
    ax.legend(fontsize=7,loc="upper left"); ax.grid(alpha=.3)
axes[0].set_ylabel("PDF (density)")
fig.suptitle("Core forest: near-nadir vs all ground-return elevation PDFs, by slope",y=1.02)
fig.savefig("figures/refdatum/nadir_vs_all_pdf_by_slope.png",dpi=100,bbox_inches="tight"); plt.close(fig)

print("core forest, ground-return elevation medians (mm):")
print(f"{'slope':>10} {'all n':>10} {'all med':>9} {'nadir n':>9} {'nadir med':>10}")
for lo,hi in bands:
    ma=base&(slp>=lo)&(slp<hi); mn=nadir&(slp>=lo)&(slp<hi)
    am=np.median(d[ma]) if ma.sum()>50 else float("nan")
    nm=np.median(d[mn]) if mn.sum()>50 else float("nan")
    print(f"{lo:>4}-{hi:<4} {ma.sum():>10,} {am:>+9.1f} {mn.sum():>9,} {nm:>+10.1f}")
print("\nwrote figures/refdatum/nadir_vs_all_pdf_by_slope.png")
