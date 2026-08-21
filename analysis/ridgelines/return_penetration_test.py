#!/usr/bin/env python3
"""Classifier-INDEPENDENT test of the forest ground rise.

The ground-class comparison uses each vendor's own ground classifier (gen1 2008 vendor
vs gen2 2021 3DEP), which differ.  To separate a real physical obscuration from a
classifier/point-selection artifact, compare the DEEPEST echoes (low percentile of ALL
returns) — the physical floor the pulses reached, independent of classification.

Logic:
  * If gen2's deepest returns reach AS LOW as gen1's ground in forest, but gen2's ground
    CLASS sits ~30 mm higher, the rise is a classifier/selection difference (gen2 ground
    surface chosen higher than the reachable floor).
  * If gen2's deepest returns CANNOT reach as low as gen1's (a floor gen2 never sees),
    the ground is physically obscured in 2021 (leaf-on groundcover/litter) — a real
    surface difference, not a classifier choice.

All heights are slope-normal d (m); d=0 is the tilted reference plane.  gen2 has ~24x the
return count, so we compare per-epoch PERCENTILES (shape), never raw counts.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/return_penetration_test.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D = "data/derived/elba_fulldensity/"
Z = np.load(D + "slope_normal_returns.npz")
edges = Z["edges"]; ctr = 0.5*(edges[:-1]+edges[1:])
pen = np.load(D + "penetration.npy")
fld = np.load(D + "floodplain_mask.npy").astype(bool)
g1ng = Z["gen1_n_ground"]; g2ng = Z["gen2_n_ground"]
g1med = Z["gen1_ground_median_d"]; g2med = Z["gen2_ground_median_d"]

have = (g1ng >= 4) & (g2ng >= 4) & ~fld & np.isfinite(pen)
forest = have & (pen < 0.25); open_ = have & (pen >= 0.45)

# materialize once — npz indexing decompresses the whole (700,508,164) array each access
G1ALL = Z["gen1_all"]; G2ALL = Z["gen2_all"]
G1GND = Z["gen1_ground"]; G2GND = Z["gen2_ground"]
HALL = {"gen1": G1ALL, "gen2": G2ALL}; HGND = {"gen1": G1GND, "gen2": G2GND}

def pooled_q(hist3d, mask, qs):
    h = hist3d[mask].sum(0).astype(float); cdf = np.cumsum(h)/h.sum()
    return {q: float(np.interp(q, cdf, ctr)) for q in qs}

qs = [0.005, 0.01, 0.02, 0.05, 0.10, 0.50]
print("POOLED low percentiles of ALL returns (deepest echoes), slope-normal d (m):")
print("  (p0.5 = the physical floor the pulses reached)")
for mask, lbl in [(forest,"forest"),(open_,"open")]:
    q1 = pooled_q(G1ALL, mask, qs)
    q2 = pooled_q(G2ALL, mask, qs)
    print(f"  {lbl}:")
    print("    gen1  " + "  ".join(f"p{q*100:g}={q1[q]:+.3f}" for q in qs))
    print("    gen2  " + "  ".join(f"p{q*100:g}={q2[q]:+.3f}" for q in qs))
    print("    d(gen2-gen1) " + "  ".join(f"p{q*100:g}={(q2[q]-q1[q])*1000:+.0f}mm" for q in qs))

# Ground-class median for reference (what the DoD uses)
print("\nGround-CLASS pooled median d (what the DoD differences):")
for mask, lbl in [(forest,"forest"),(open_,"open")]:
    m1 = pooled_q(G1GND, mask, [0.5])[0.5]
    m2 = pooled_q(G2GND, mask, [0.5])[0.5]
    print(f"  {lbl}: gen1={m1*1000:+.0f}mm  gen2={m2*1000:+.0f}mm  rise={ (m2-m1)*1000:+.0f}mm")

# --- decisive per-cell test: can gen2 reach gen1's ground floor? ----------------------
# For each forest cell: does gen2_all have appreciable mass BELOW gen1's ground median?
# If yes, gen2 physically reaches that low but its ground CLASS sits higher -> selection.
print("\nPer-cell: fraction of gen2 ALL-return mass at or below gen1 ground median")
edL = edges[:-1]; edR = edges[1:]
def mass_below(hist3d_row, thr):
    # fraction of a cell's histogram mass below height thr (linear within bin)
    h = hist3d_row.astype(float); tot = h.sum()
    if tot == 0: return np.nan
    frac_bin = np.clip((thr - edL)/(edR - edL), 0, 1)
    return (h*frac_bin).sum()/tot

G2ALL = Z["gen2_all"]   # materialize once (npz access decompresses the whole array)
for mask, lbl in [(forest,"forest"),(open_,"open")]:
    idx = np.argwhere(mask)
    # sample up to 8000 cells for speed
    if len(idx) > 8000:
        sel = idx[np.linspace(0, len(idx)-1, 8000).astype(int)]
    else:
        sel = idx
    fr = []
    for iy, ix in sel:
        thr = g1med[iy, ix]
        if not np.isfinite(thr): continue
        fr.append(mass_below(G2ALL[iy, ix], thr))
    fr = np.array(fr); fr = fr[np.isfinite(fr)]
    print(f"  {lbl}: median {np.median(fr)*100:.1f}%  mean {np.mean(fr)*100:.1f}%  "
          f"(frac of gen2 returns that DO reach below gen1 ground median)  n={fr.size}")

# --- figure: deep-tail CDF of ALL returns, forest, zoomed near ground -----------------
fig, axes = plt.subplots(1, 2, figsize=(13,6), sharey=True)
for ax, mask, lbl in [(axes[0],forest,"FOREST"),(axes[1],open_,"OPEN")]:
    for gen, col in [("gen1","C0"),("gen2","C3")]:
        h = HALL[gen][mask].sum(0).astype(float); cdf = np.cumsum(h)/h.sum()
        ax.plot(cdf*100, ctr, col, lw=1.8, label=f"{gen} all")
    for gen, col, ls in [("gen1","C0",":"),("gen2","C3",":")]:
        m = pooled_q(HGND[gen], mask, [0.5])[0.5]
        ax.axhline(m, color=col, ls=ls, lw=1.2, alpha=.8)
    ax.set_ylim(-0.8, 1.5); ax.set_xlim(0, 20)
    ax.set_xlabel("cumulative % of all returns (from bottom)")
    ax.set_title(f"{lbl}: deep-echo CDF (dotted = ground-class median)")
    ax.grid(alpha=.3); ax.legend()
axes[0].set_ylabel("slope-normal height d (m)")
fig.suptitle("Deepest echoes vs ground-class median — does gen2 reach gen1's floor?")
fig.savefig("figures/refdatum/return_penetration.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("\nwrote figures/refdatum/return_penetration.png")
