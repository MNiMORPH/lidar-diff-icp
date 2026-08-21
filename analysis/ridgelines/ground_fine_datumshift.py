#!/usr/bin/env python3
"""Replot the 1 cm ground-class distributions (linear + log density) WITH the vertical
datum applied: the pipeline adds const_m = +0.067 m to gen1 (the GEOID03->GEOID18 geoid
difference), so in the d-frame gen1 shifts up by +67 mm.  Plotting gen1 at fc+DATUM is
exact (no re-binning).  After the shift the bulk 'gen2 everywhere higher' offset is gone
and only the leaf-state residual (forest > open) remains.

Note: this pooled-by-return view lands open near 52-67 = -15 mm; the finished per-cell DoD
grid gives open -6.4 mm.  The ~10 mm gap is method (pooled vs per-cell + lateral/tilt), not
a second datum term.

    ./lidar-icp/bin/python analysis/ridgelines/ground_fine_datumshift.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

DATUM = 0.067   # m, added to gen1 (geoid GEOID03 - GEOID18)
Z = np.load("data/derived/elba_fulldensity/ground_fine_pooled.npz")
fc = Z["fc"]; FW = fc[1]-fc[0]
def dens(gen, s):
    c = Z[f"{gen}_{s}"].astype(float); return c/c.sum()/FW
def med(gen, s, shift=0.0):
    c = Z[f"{gen}_{s}"].astype(float); cdf = np.cumsum(c)/c.sum()
    return np.interp(.5, cdf, fc) + shift

# gen1 y-coordinate shifted up by the datum; gen2 unchanged
yc = {"gen1": fc + DATUM, "gen2": fc}
sh = {"gen1": DATUM, "gen2": 0.0}

fig, axes = plt.subplots(2, 2, figsize=(13, 11), sharey=True)
for j, (s, strat) in enumerate([(1, "FOREST"), (2, "OPEN")]):
    for i, (scale, tag) in enumerate([("linear", "linear density"), ("log", "log density")]):
        ax = axes[i, j]
        for gen, col in [("gen1", "C0"), ("gen2", "C3")]:
            d = dens(gen, s)
            if scale == "log":
                ax.semilogx(np.where(d > 0, d, np.nan), yc[gen], col, lw=1.7,
                            label=f"{gen} (datum-shifted)" if gen == "gen1" else gen)
            else:
                ax.plot(d, yc[gen], col, lw=1.7,
                        label=f"{gen} (+{DATUM*1000:.0f} mm datum)" if gen == "gen1" else gen)
            ax.axhline(med(gen, s, sh[gen]), color=col, ls=":", lw=1)
        rf = (med("gen2", s) - med("gen1", s, DATUM)) * 1000
        ax.axhline(0, color="k", lw=.5); ax.set_ylim(-0.5, 0.7)
        if scale == "log": ax.set_xlim(1e-2, 3e1); ax.grid(alpha=.3, which="both")
        else: ax.grid(alpha=.3)
        ax.set_xlabel(f"density (1/m{', log' if scale=='log' else ''})")
        ax.set_title(f"{strat}: {tag}  —  residual gen2-gen1 = {rf:+.0f} mm")
        ax.legend(fontsize=8)
axes[0, 0].set_ylabel("slope-normal d (m), datum-corrected")
axes[1, 0].set_ylabel("slope-normal d (m), datum-corrected")
fig.suptitle(f"1 cm ground-class with vertical datum applied (gen1 +{DATUM*1000:.0f} mm): "
             "bulk offset removed, leaf-state residual remains", y=0.995)
fig.savefig("figures/refdatum/ground_fine_datumshift.png", dpi=130, bbox_inches="tight"); plt.close(fig)

print(f"datum applied to gen1: +{DATUM*1000:.0f} mm")
for s, strat in [(1, "forest"), (2, "open")]:
    print(f"  {strat}: gen1 median {med('gen1',s):+.3f} -> {med('gen1',s,DATUM):+.3f} m ;"
          f"  gen2 {med('gen2',s):+.3f} m ;  residual gen2-gen1 = "
          f"{(med('gen2',s)-med('gen1',s,DATUM))*1000:+.0f} mm")
print("wrote figures/refdatum/ground_fine_datumshift.png")
