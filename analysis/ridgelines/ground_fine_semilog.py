#!/usr/bin/env python3
"""Replot the 1 cm ground-class distributions with a LOG density axis (semi-log), to
expose the tails: gen1's deeper low tail (leaf-off penetration) vs the near-coincident
upper tails. Reuses ground_fine_pooled.npz (no re-streaming).

    ./lidar-icp/bin/python analysis/ridgelines/ground_fine_semilog.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Z = np.load("data/derived/elba_fulldensity/ground_fine_pooled.npz")
fc = Z["fc"]; FW = fc[1]-fc[0]
def dens(gen,s):
    c = Z[f"{gen}_{s}"].astype(float); return c/c.sum()/FW
def med(gen,s):
    c = Z[f"{gen}_{s}"].astype(float); cdf=np.cumsum(c)/c.sum(); return np.interp(.5,cdf,fc)

fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
for ax, s, lbl in [(axes[0], 1, "FOREST"), (axes[1], 2, "OPEN")]:
    for gen, col in [("gen1", "C0"), ("gen2", "C3")]:
        d = dens(gen, s)
        ax.semilogx(np.where(d > 0, d, np.nan), fc, col, lw=1.7, label=f"{gen} ground-class")
        ax.axhline(med(gen, s), color=col, ls=":", lw=1)
    ax.axhline(0, color="k", lw=.5)
    ax.set_ylim(-0.7, 0.9); ax.set_xlim(1e-2, 3e1)
    ax.set_xlabel("density (1/m, log)"); ax.set_title(f"{lbl}: fine 1 cm ground-class (semi-log)")
    ax.grid(alpha=.3, which="both"); ax.legend(fontsize=9)
axes[0].set_ylabel("slope-normal d (m)  [plane = gen2 bare earth]")
fig.suptitle("Ground-class return density, log scale — gen1 deeper low tail; upper tails ~coincide")
fig.savefig("figures/refdatum/ground_fine_semilog.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote figures/refdatum/ground_fine_semilog.png")
