#!/usr/bin/env python3
"""Test Andy's model of the forest ground rise:
    gen2 ground-class column = gen1's true-ground return (still detected) + an added
    UPWARD shoulder from leaf-on herbaceous GROUNDCOVER (~0.1-0.6 m above ground).

Distinct from a bulk shift. Predictions, if the model holds:
  (P1) gen2 retains mass at/below gen1's ground mode (true ground still seen), and gains
       mass mainly in the groundcover band above it -> asymmetric, one-sided upper excess.
  (P2) the per-cell median RISE grows with a GROUNDCOVER proxy taken from the RIGHT band
       (d in 0.1-0.6 m) -- unlike understory_frac (0.5-2 m shrubs), which missed it.
  (P3) THE PAYOFF: estimating gen2 ground by a LOWER percentile (recovering the true ground
       under the groundcover) shrinks the forest-specific rise toward the open value.

Everything slope-normal d (m), common gen2-bare-earth frame. gen2 has ~24x the returns;
compare per-epoch shape/percentiles, never raw counts.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/groundcover_decomp.py
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D = "data/derived/elba_fulldensity/"
Z = np.load(D + "slope_normal_returns.npz")
edges = Z["edges"]; ctr = 0.5*(edges[:-1]+edges[1:]); edL = edges[:-1]; edR = edges[1:]
pen = np.load(D + "penetration.npy")
fld = np.load(D + "floodplain_mask.npy").astype(bool)
g1ng = Z["gen1_n_ground"]; g2ng = Z["gen2_n_ground"]
g1med = Z["gen1_ground_median_d"]; g2med = Z["gen2_ground_median_d"]
g1p10 = Z["gen1_ground_p10_d"];    g2p10 = Z["gen2_ground_p10_d"]
G1G = Z["gen1_ground"]; G2G = Z["gen2_ground"]     # materialize once

have = (g1ng >= 8) & (g2ng >= 8) & ~fld & np.isfinite(pen)   # >=8 ground for stable percentiles
forest = have & (pen < 0.25); open_ = have & (pen >= 0.45)
print(f"cells: forest={forest.sum()} open={open_.sum()} (>=8 ground each, hillslope)")

# ---------- P1: pooled shape, gen1 vs gen2, matched on the DEEP flank ----------------
# Scale so the two agree on the low flank (d in [-0.6,-0.3], only-true-ground region);
# then any gen2 EXCESS above ground is the added component, not an artifact of area-norm.
def pooled(hist3d, mask):
    return hist3d[mask].sum(0).astype(float)
h1 = pooled(G1G, forest); h2 = pooled(G2G, forest)
d1 = h1/h1.sum(); d2 = h2/h2.sum()                 # unit-area densities (per bin), NO scaling
mode1 = ctr[np.argmax(d1)]; mode2 = ctr[np.argmax(d2)]
mode = 0.5*(mode1+mode2)                            # common ground mode
# Decompose each distribution's mass about the common mode: shoulder (above) vs tail (below)
def split(d):
    return d[ctr>mode].sum(), d[ctr<mode].sum()     # (upper shoulder frac, lower tail frac)
up1,lo1 = split(d1); up2,lo2 = split(d2)
print(f"\nP1 pooled forest ground-class (unit area, no scaling):")
print(f"  ground MODE: gen1 {mode1:+.3f}  gen2 {mode2:+.3f} m  (shift {(mode2-mode1)*1000:+.0f} mm) -> common ground return preserved")
print(f"  mass ABOVE mode (upper shoulder): gen1 {up1*100:.1f}%  gen2 {up2*100:.1f}%  "
      f"-> gen2 adds {(up2-up1)*100:+.1f} pts (groundcover shoulder)")
print(f"  mass BELOW mode (deeper tail):    gen1 {lo1*100:.1f}%  gen2 {lo2*100:.1f}%  "
      f"-> gen1 has {(lo1-lo2)*100:+.1f} pts more (leaf-off penetration)")
diff = d2 - d1
fig, ax = plt.subplots(1, 2, figsize=(13,6))
ax[0].plot(d1, ctr, "C0", lw=1.8, label="gen1 ground (2008 leaf-off)")
ax[0].plot(d2, ctr, "C3", lw=1.8, label="gen2 ground (2021 leaf-on)")
ax[0].axhline(mode,color="0.4",ls="--",lw=1,label=f"common mode {mode:+.2f} m")
ax[0].axhspan(0.0,0.6, color="#d8f0d8", alpha=.5, label="groundcover band (0-0.6 m)")
ax[0].axhline(0,color="k",lw=.5); ax[0].set_ylim(-0.8,1.2)
ax[0].set_xlabel("unit-area density"); ax[0].set_ylabel("slope-normal d (m)")
ax[0].set_title("Forest ground-class: SAME mode, two-sided difference"); ax[0].legend(fontsize=8)
ax[1].fill_betweenx(ctr, 0, diff, where=diff>0, color="C2", alpha=.6, label="gen2 excess (groundcover up)")
ax[1].fill_betweenx(ctr, 0, diff, where=diff<0, color="C1", alpha=.6, label="gen1 excess (deeper penetration)")
ax[1].axhline(mode,color="0.4",ls="--",lw=1); ax[1].axhline(0,color="k",lw=.5); ax[1].axvline(0,color="k",lw=.5)
ax[1].set_ylim(-0.8,1.2); ax[1].set_xlabel("gen2 - gen1 unit-area density")
ax[1].set_title("Two mechanisms, opposite signs of DoD"); ax[1].legend(fontsize=8)
fig.savefig("figures/refdatum/groundcover_decomp.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ---------- P2: per-cell median rise vs a RIGHT-BAND groundcover proxy ---------------
# groundcover proxy = fraction of gen2 ground-class returns in d in (0.1, 0.6] (herb layer)
def band_frac(hist3d, mask, lo, hi):
    h = hist3d.astype(float)
    w = np.clip((np.minimum(edR,hi)-np.maximum(edL,lo))/(edR-edL),0,1)   # partial-bin overlap
    num = (h*w).sum(-1); den = h.sum(-1)
    out = np.full(den.shape, np.nan); m = den>0; out[m]=num[m]/den[m]; return out
gc2 = band_frac(G2G, None, 0.1, 0.6)     # gen2 groundcover fraction, per cell
rise = (g2med - g1med)*1000
mf = forest & np.isfinite(gc2) & np.isfinite(rise)
g = gc2[mf]; rv = rise[mf]
print("\nP2 forest median-rise vs gen2 groundcover fraction d(0.1-0.6 m):")
e = np.quantile(g, np.linspace(0,1,7))
for i in range(len(e)-1):
    b = (g>=e[i]) & (g<=e[i+1] if i==len(e)-2 else g<e[i+1])
    if b.sum()<50: continue
    print(f"  gc {e[i]:.3f}-{e[i+1]:.3f}: median rise {np.median(rv[b]):+.1f} mm  n={b.sum()}")
print(f"  corr(gc2, rise) = {np.corrcoef(g,rv)[0,1]:+.3f}   "
      f"(cf. understory_frac r=+0.05 = wrong band)")

# ---------- P3: does a lower percentile cancel the forest-specific rise? -------------
# per-cell p25 from the ground-class histograms (cumsum interp)
def per_cell_q(hist3d, q):
    h = hist3d.astype(float); c = np.cumsum(h,axis=-1); tot = c[...,-1:]
    out = np.full(h.shape[:-1], np.nan); m = tot[...,0]>0
    cn = np.where(tot>0, c/np.where(tot>0,tot,1), np.nan)
    # vectorized interp per cell is awkward; loop only over masked cells
    idx = np.argwhere(m)
    for iy,ix in idx:
        out[iy,ix] = np.interp(q, cn[iy,ix], ctr)
    return out
print("\nP3 forest-specific rise by ground estimator (per-cell, median across cells):")
print(f"  {'estimator':10s} {'forest rise':>12s} {'open rise':>11s} {'forest-open':>12s}")
for name, e1, e2 in [("p50 (now)", g1med, g2med), ("p10 (saved)", g1p10, g2p10)]:
    r = (e2-e1)*1000
    rf = np.median(r[forest & np.isfinite(r)]); ro = np.median(r[open_ & np.isfinite(r)])
    print(f"  {name:10s} {rf:>10.1f}mm {ro:>9.1f}mm {rf-ro:>10.1f}mm")
# p25 computed on the fly
g1p25 = per_cell_q(G1G, 0.25); g2p25 = per_cell_q(G2G, 0.25)
r = (g2p25-g1p25)*1000
rf = np.median(r[forest & np.isfinite(r)]); ro = np.median(r[open_ & np.isfinite(r)])
print(f"  {'p25 (calc)':10s} {rf:>10.1f}mm {ro:>9.1f}mm {rf-ro:>10.1f}mm")
print("\n  (forest-open = the anomaly that survives the open-tied datum; smaller = estimator "
      "recovers true ground under groundcover)")
print("\nwrote figures/refdatum/groundcover_decomp.png")
