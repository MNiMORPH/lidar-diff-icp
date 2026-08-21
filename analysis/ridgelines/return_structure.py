#!/usr/bin/env python3
"""Vertical structure of lidar returns (slope-normal height d), gen1 vs gen2, to find
what drives the FALSE forest elevation rise (gen2 ground surface sits above gen1's).

Data: data/derived/elba_fulldensity/slope_normal_returns.npz
  gen{1,2}_all[iy,ix,:]    per-cell histogram of ALL returns over d-bins (edges, 0.25 m)
  gen{1,2}_ground[iy,ix,:] per-cell histogram of GROUND-classified returns over d-bins
  d = slope-normal distance above the tilted ground plane; d=0 is the ground plane.

The DoD grids on the per-cell MEDIAN of ground-class d (ground_q=0.50), so the forest
rise is  median(gen2_ground_d) - median(gen1_ground_d) > 0.  This script asks whether
that is (a) real ground motion, (b) gen2 understory leaking into the "ground" class and
pulling its median up, or (c) gen1 ground returns behaving differently (sparser/lower).

Layers we try to separate, by slope-normal height d:
  ground      d in [-0.25, +0.25]   (the return that defines the surface)
  understory  d in [+0.5, +3.0]     (shrubs, forbs, tree boles low)
  canopy      d > +3.0

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/return_structure.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D = "data/derived/elba_fulldensity/"
Z = np.load(D + "slope_normal_returns.npz")
edges = Z["edges"]; ctr = 0.5*(edges[:-1]+edges[1:]); bw = np.diff(edges)
ny, nx, nb = Z["gen1_all"].shape

pen   = np.load(D + "penetration.npy")          # ground-return fraction proxy; low=forest
veg   = np.load(D + "canopy_struct.npz")["veg_frac"]
under = np.load(D + "canopy_struct.npz")["understory_frac"]
fld   = np.load(D + "floodplain_mask.npy").astype(bool)
g1med = Z["gen1_ground_median_d"]; g2med = Z["gen2_ground_median_d"]
g1ng  = Z["gen1_n_ground"]; g2ng = Z["gen2_n_ground"]

# --- strata: forest / open, hillslope only (drop floodplain) -------------------------
have = (g1ng >= 4) & (g2ng >= 4) & ~fld & np.isfinite(pen)
forest = have & (pen < 0.25)
open_  = have & (pen >= 0.45)
print(f"cells: forest={forest.sum()}  open={open_.sum()}  (hillslope, both epochs >=4 ground)")

def stack_profile(hist3d, mask):
    """Sum per-cell histograms over a mask, return density (per m, area-normalized)."""
    h = hist3d[mask].sum(0).astype(float)
    dens = h / h.sum() / bw
    return dens, h

# ============ FIG 1: full vertical structure, all returns, forest vs open ============
fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
for ax, mask, lbl in [(axes[0], forest, "FOREST (pen<0.25)"), (axes[1], open_, "OPEN (pen>=0.45)")]:
    d1, _ = stack_profile(Z["gen1_all"], mask)
    d2, _ = stack_profile(Z["gen2_all"], mask)
    ax.plot(d1, ctr, color="C0", lw=1.6, label="gen1 (2008) all")
    ax.plot(d2, ctr, color="C3", lw=1.6, label="gen2 (2021) all")
    ax.axhline(0, color="k", lw=.5); ax.set_ylim(-1, 20)
    ax.set_xscale("log"); ax.set_xlabel("return density (1/m, area-norm)")
    ax.set_title(lbl); ax.grid(alpha=.3); ax.legend()
axes[0].set_ylabel("slope-normal height d (m)")
fig.suptitle("Vertical structure of ALL returns — layered ground/understory/canopy")
fig.savefig("figures/refdatum/return_structure_all.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ============ FIG 2: near-ground GROUND-CLASS structure (the median driver) ==========
fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
for ax, mask, lbl in [(axes[0], forest, "FOREST"), (axes[1], open_, "OPEN")]:
    d1, _ = stack_profile(Z["gen1_ground"], mask)
    d2, _ = stack_profile(Z["gen2_ground"], mask)
    ax.plot(d1, ctr, color="C0", lw=1.8, label="gen1 ground")
    ax.plot(d2, ctr, color="C3", lw=1.8, label="gen2 ground")
    ax.axhline(0, color="k", lw=.5); ax.set_ylim(-0.8, 3.0)
    ax.set_xlabel("ground-return density (1/m, area-norm)")
    ax.set_title(lbl+": ground-class returns"); ax.grid(alpha=.3); ax.legend()
axes[0].set_ylabel("slope-normal height d (m)")
fig.suptitle("Ground-classified returns near the surface — does gen2 have an upper tail?")
fig.savefig("figures/refdatum/return_structure_ground.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# --- quantitative: pooled ground-class quantiles + upper-tail fraction ---------------
def pooled_quantiles(hist3d, mask, qs=(0.10,0.25,0.50,0.75,0.90)):
    h = hist3d[mask].sum(0).astype(float); cdf = np.cumsum(h)/h.sum()
    return {q: np.interp(q, cdf, ctr) for q in qs}, h
def frac_above(h, thr):
    return h[ctr > thr].sum()/h.sum()

print("\nPOOLED ground-class slope-normal quantiles (m):")
for mask, lbl in [(forest,"forest"),(open_,"open")]:
    q1,h1 = pooled_quantiles(Z["gen1_ground"], mask)
    q2,h2 = pooled_quantiles(Z["gen2_ground"], mask)
    print(f"  {lbl}:")
    print(f"    gen1  " + "  ".join(f"p{int(q*100):02d}={v:+.3f}" for q,v in q1.items())
          + f"   frac(d>0.5m)={frac_above(h1,0.5)*100:.1f}%  frac(d>1m)={frac_above(h1,1.0)*100:.2f}%")
    print(f"    gen2  " + "  ".join(f"p{int(q*100):02d}={v:+.3f}" for q,v in q2.items())
          + f"   frac(d>0.5m)={frac_above(h2,0.5)*100:.1f}%  frac(d>1m)={frac_above(h2,1.0)*100:.2f}%")

# --- per-cell median difference (what the DoD actually differences) ------------------
dmed = (g2med - g1med)*1000   # mm; slope-normal, before datum. +ve = gen2 ground higher
print("\nPer-cell ground-median difference gen2-gen1 (mm, slope-normal, pre-datum):")
for mask, lbl in [(forest,"forest"),(open_,"open")]:
    v = dmed[mask]; v = v[np.isfinite(v)]
    print(f"  {lbl}: median={np.median(v):+.1f}  mean={np.mean(v):+.1f}  "
          f"p25={np.percentile(v,25):+.1f}  p75={np.percentile(v,75):+.1f}  n={v.size}")

# --- correlate the rise with understory presence ------------------------------------
print("\nForest ground-median rise vs understory fraction (binned):")
m = forest & np.isfinite(under) & np.isfinite(dmed)
u = under[m]; dv = dmed[m]
e = np.quantile(u, np.linspace(0,1,7))
for i in range(len(e)-1):
    b = (u>=e[i])&(u<e[i+1] if i<len(e)-2 else u<=e[i+1])
    if b.sum()<50: continue
    print(f"  understory {e[i]:.2f}-{e[i+1]:.2f}: median rise {np.median(dv[b]):+.1f} mm  n={b.sum()}")
r = np.corrcoef(u, dv)[0,1]
print(f"  corr(understory_frac, rise) = {r:+.3f}")

print("\nwrote figures/refdatum/return_structure_all.png, return_structure_ground.png")
