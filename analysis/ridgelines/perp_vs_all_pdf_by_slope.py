#!/usr/bin/env python3
"""Forest (PyForestScan cover>=0.5): PDF of near-slope-PERPENDICULAR ground-return elevations vs PDF of ALL
ground-return elevations, as a function of slope. near-perpendicular = incidence to surface
<= 2 deg. Elevation = d_mm (contains ~67 mm constant datum). gen1 CSF cloth ground.
(Expect the perpendicular sample to dwindle as slope steepens: scanner +-17deg can't reach
perpendicular on steep slopes.)

    ./lidar-icp/bin/python analysis/ridgelines/perp_vs_all_pdf_by_slope.py
"""
import argparse, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
A_ = ap.parse_args()
TILE = os.path.basename(A_.tile.rstrip("/"))
TAG = "" if TILE == "elba_fulldensity" else f"_{TILE}"
TITLE_TAG = "" if not TAG else f"  ({TILE})"
Z=np.load(f"{A_.tile}/gen1_csf_angles.npz")
d=Z["d_mm"]; inc=Z["incidence"]; slp=Z["slope"]; cf=Z["pfs_forest"]; ing=Z["in_grid"]
base=cf&ing&np.isfinite(d)&np.isfinite(inc)&np.isfinite(slp)
perp=base&(inc<=2.0)

bands=[(0,8),(8,16),(16,24),(24,32),(32,45)]
bins=np.arange(-350,251,10)
fig,axes=plt.subplots(1,len(bands),figsize=(2.9*len(bands),4.6),sharey=True)
for ax,(lo,hi) in zip(axes,bands):
    ma=base&(slp>=lo)&(slp<hi); mp=perp&(slp>=lo)&(slp<hi)
    if ma.sum()>50:
        ax.hist(d[ma],bins=bins,density=True,histtype="step",lw=2,color="0.4",
                label=f"all (n={ma.sum():,}, med {np.median(d[ma]):+.0f})")
        ax.axvline(np.median(d[ma]),color="0.4",ls=":",lw=1)
    if mp.sum()>50:
        ax.hist(d[mp],bins=bins,density=True,histtype="step",lw=2,color="C3",
                label=f"perp (n={mp.sum():,}, med {np.median(d[mp]):+.0f})")
        ax.axvline(np.median(d[mp]),color="C3",ls=":",lw=1)
    else:
        ax.text(0.5,0.5,f"perp n={mp.sum()}\n(too few)",transform=ax.transAxes,ha="center",color="C3")
    ax.set_title(f"slope {lo}-{hi}°"); ax.set_xlabel("ground-return elevation d (mm)")
    ax.legend(fontsize=7,loc="upper left"); ax.grid(alpha=.3)
axes[0].set_ylabel("PDF (density)")
fig.suptitle(f"Forest (PyForestScan cover>=0.5): near-slope-perpendicular vs all ground-return elevation PDFs, by slope{TITLE_TAG}",y=1.02)
fig.savefig(f"figures/refdatum/perp_vs_all_pdf_by_slope{TAG}.png",dpi=100,bbox_inches="tight"); plt.close(fig)

print("forest (PFS cover>=0.5), ground-return elevation medians (mm); perp = incidence<=2deg:")
print(f"{'slope':>10} {'all n':>10} {'all med':>9} {'perp n':>9} {'perp med':>10}")
for lo,hi in bands:
    ma=base&(slp>=lo)&(slp<hi); mp=perp&(slp>=lo)&(slp<hi)
    am=np.median(d[ma]) if ma.sum()>50 else float("nan")
    pm=np.median(d[mp]) if mp.sum()>50 else float("nan")
    print(f"{lo:>4}-{hi:<4} {ma.sum():>10,} {am:>+9.1f} {mp.sum():>9,} {pm:>+10.1f}")
print(f"\nwrote figures/refdatum/perp_vs_all_pdf_by_slope{TAG}.png")
