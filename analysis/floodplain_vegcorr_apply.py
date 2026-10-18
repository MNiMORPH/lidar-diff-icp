#!/usr/bin/env python3
"""Apply the 3x3 logistic transfer pct0(span) as a vegetation correction to gen1's ground on
the flat floodplain, erosion-side only, tapered -- and measure whether the vegetation
"erosion" patches fade into the DoD's background.

Andy, 2026-10-18: accept the mat hypothesis as likely correct; it gives the logistic a
reason. "Could we use the 3x3 logistic fit to correct the observed floodplain elevations?
Deposition is possible, so correct only those that would show erosion; taper in a
reasonable way instead of making a hard cutoff. The low R^2 does not mean failure: mostly
we are seeing real topographic + vegetation noise, but there is a coherent global signal
behind it. The test of success is that the patches of 'erosion' due to vegetation are
removed and fade into the background ups and downs of the lidar."

WHAT THIS IS, stated once. A cross-epoch calibration: gen1's ground is moved to the
percentile of its own pooled cloud at which gen2's surface sits on average at that
vegetation span. The mean floodplain DoD at a given span therefore goes to ~0 by
construction, and real erosion in vegetated cells would be absorbed. Andy reaffirmed with
that stated. The fitted parameters are site-specific; the METHOD travels, the numbers do
not.

THE CORRECTION, per flat-floodplain cell (floodplain mask AND terrain slope <= 2 deg):
  1. pool the cell's 3x3 neighbourhood of gen1 slope-normal residuals d (mm vs gen2);
  2. span = p90 - p10 of the pooled cloud; pct = logistic(span), REFITTED HERE on the same
     histograms it is applied to (same bins, same window, same population);
  3. gen1's ground estimate becomes the pct-th percentile of the pooled cloud; the DoD
     change is  delta = (median of the cell's own returns) - (that percentile),
     the same bookkeeping as gen1_coarsened_ground.py;
  4. EROSION-SIDE ONLY, TAPERED (both MINE, as proposals): only delta > 0 (raising the DoD)
     is applied, with weight w = clip(-DoD_raw / T, 0, 1) where T is the NMAD of the DoD on
     stable cells -- a measured noise scale, so a cell already reading deposition gets
     nothing, a cell at -T or below gets the full correction, and between them it tapers.

SUCCESS TESTS, before vs after, all on the flat floodplain:
  * DoD-vs-span coefficient and the quartile medians (the artifact's signature);
  * NMAD of the DoD on the floodplain against NMAD on stable ground (the background);
  * share of floodplain cells below -LoD, and the size of their connected patches;
  * maps: raw DoD, corrected DoD, the applied correction.

Outputs (next to the product, clearly named): dod_floodplain_vegcorr.npy,
floodplain_vegcorr_delta.npy, floodplain_vegcorr.json, figures/floodplain_vegcorr_*.png.

    ./lidar-icp/bin/python analysis/floodplain_vegcorr_apply.py
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq
from scipy.ndimage import label, uniform_filter
from scipy.optimize import least_squares

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext_delong_star_nosurf")
ap.add_argument("--pool", type=int, default=3)
ap.add_argument("--slope-max", type=float, default=2.0, help="deg; the flat-floodplain cut used throughout")
ap.add_argument("--figdir", default="figures")
A = ap.parse_args()
S = A.tile; P = A.pool
LO, HI, DZ = -5000.0, 5000.0, 20.0
edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size

dod = np.load(f"{S}/dod.npy") * 1000.0
sl = np.load(f"{S}/slope.npy"); fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
st = np.load(f"{S}/stable.npy").astype(bool); lod = np.load(f"{S}/lod.npy") * 1000.0
ny, nx = dod.shape
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float)
ok = np.isfinite(d) & (d > LO) & (d < HI); cell, d = cell[ok], d[ok]
H = np.zeros((ny * nx, NB), np.int32)
np.add.at(H, (cell, np.clip(((d - LO) / DZ).astype(int), 0, NB - 1)), 1)
H = H.reshape(ny, nx, NB); n_own = H.sum(2)
flat = fp & (sl <= A.slope_max) & np.isfinite(dod) & (n_own >= 3)
Hp = np.rint(uniform_filter(H.astype(np.float32), size=(P, P, 1), mode="nearest") * P * P)


def qtl(Hx, q):
    """Per-cell quantile q (scalar or per-cell array) read off histogram rows."""
    c = np.cumsum(Hx, -1); tot = c[..., -1]
    qq = np.broadcast_to(np.asarray(q, float), tot.shape)
    out = np.full(tot.shape, np.nan); m = tot > 0
    out[m] = ctr[(c[m] >= qq[m][:, None] * tot[m][:, None]).argmax(-1)]
    return out


Hf = Hp.reshape(-1, NB)[flat.ravel()]
span = qtl(Hf, 0.90) - qtl(Hf, 0.10)
pct0 = 100.0 * Hf[:, ctr <= 0].sum(1) / Hf.sum(1)
logistic = lambda p, x: p[1] + (p[0] - p[1]) / (1 + np.exp(np.clip((x - p[2]) / p[3], -50, 50)))
par = least_squares(lambda p: logistic(p, span) - pct0, x0=[65, 20, 250, 60], loss="soft_l1", f_scale=10,
                    bounds=([0, 0, 0, 1], [100, 100, 5000, 5000])).x
T = 1.4826 * np.median(np.abs(dod[st] - np.median(dod[st])))
print("INPUTS: beam_offset_table d_mm_corr (gen1 slope-normal residual vs gen2, mm); dod, slope, floodplain_mask, stable, lod")
print(f"POPULATION corrected: flat floodplain (mask AND slope <= {A.slope_max:g} deg): {flat.sum():,} cells; pooling {P}x{P}")
print(f"LOGISTIC refitted here: p_hi {par[0]:.1f}  p_lo {par[1]:.1f}  s0 {par[2]:.0f}  w {par[3]:.0f}")
print(f"TAPER scale T = NMAD of DoD on stable cells = {T:.1f} mm (measured); weight = clip(-DoD_raw/T, 0, 1); erosion side only")

est_own = qtl(H.reshape(-1, NB)[flat.ravel()], 0.50)                 # the product's gen1 ground, per cell
pct = logistic(par, span) / 100.0
est_new = qtl(Hf, pct)
delta_full = est_own - est_new                                        # mm added to the DoD if applied in full
raw = dod[flat]
w = np.clip(-raw / T, 0, 1)
delta = np.where(delta_full > 0, w * delta_full, 0.0)
dod_new = dod.copy(); dod_new[flat] = raw + delta
D = np.zeros_like(dod); D[flat] = delta
np.save(f"{S}/dod_floodplain_vegcorr.npy", dod_new / 1000.0)
np.save(f"{S}/floodplain_vegcorr_delta.npy", D / 1000.0)

# ----------------------------------------------------------------- success tests
sp_own = np.load(f"{S}/gen1_spread_mm.npy")[flat]
qs = np.percentile(span, [25, 50, 75])
qm = [(span >= lo) & (span < hi) for lo, hi in zip([0] + list(qs), list(qs) + [1e9])]
nmad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))
print(f"\nSPAN QUARTILES (pooled {P}x{P} span, mm): edges {qs.round(0)}")
print(f"  {'DoD on flat floodplain (mm)':<34}{'Q1':>8}{'Q2':>8}{'Q3':>8}{'Q4':>8}{'coeff':>8}{'NMAD':>7}{'median':>8}")
for lab, v in (("raw", raw), ("corrected", raw + delta), ("full correction, no taper", raw + np.maximum(delta_full, 0)),
               ("full correction, both signs", raw + delta_full)):
    row = [np.median(v[m]) for m in qm]; c = np.polyfit(span, v, 1)[0]
    print(f"  {lab:<34}" + "".join(f"{x:>8.1f}" for x in row) + f"{c:>8.3f}{nmad(v):>7.1f}{np.median(v):>8.1f}")
print(f"  {'stable ground (background)':<34}{'':>32}{'':>8}{nmad(dod[st]):>7.1f}{np.median(dod[st]):>8.1f}")
print(f"\nCORRECTION applied: cells with delta > 0: {np.mean(delta > 0)*100:.1f}%; median delta where applied {np.median(delta[delta > 0]):.1f} mm; "
      f"p90 {np.percentile(delta[delta > 0], 90):.1f}; max {delta.max():.1f}")

print("\nEROSION PATCHES below -LoD on the flat floodplain (8-connected):")
lodf = lod
for lab, field in (("raw", dod), ("corrected", dod_new)):
    below = flat & np.isfinite(field) & (field < -lodf)
    lab_, n = label(below, structure=np.ones((3, 3)))
    sizes = np.bincount(lab_.ravel())[1:] if n else np.array([])
    print(f"  {lab:<11} cells below -LoD {below.sum():>6,} ({100*below.sum()/flat.sum():.1f}% of flat floodplain); "
          f"patches {n:>5,}; cells in patches >= 20 cells (0.05 ha): {int(sizes[sizes >= 20].sum()) if n else 0:>6,}; largest {int(sizes.max()) if n else 0:,}")
above = lambda f: (flat & np.isfinite(f) & (f > lodf)).sum()
print(f"  deposition cells above +LoD: raw {above(dod):,}  corrected {above(dod_new):,}  (must be unchanged: erosion side only)")

json.dump(dict(tile=S, pool=P, slope_max_deg=A.slope_max, logistic=dict(p_hi=float(par[0]), p_lo=float(par[1]), s0=float(par[2]), w=float(par[3])),
               taper_T_mm=float(T), erosion_side_only=True, n_cells=int(flat.sum()),
               note="cross-epoch calibration; floodplain DoD at a given span -> ~0 by construction; site-fitted numbers"),
          open(f"{S}/floodplain_vegcorr.json", "w"), indent=1)

# ----------------------------------------------------------------- maps
iy, ix = np.nonzero(fp); r0, r1, c0, c1 = iy.min(), iy.max() + 1, ix.min(), ix.max() + 1
fig, axes = plt.subplots(1, 3, figsize=(18, 6.5))
v = 150.0
for ax, (lab, field, cm, lim) in zip(axes, (("raw DoD", dod, "RdBu", v), ("corrected DoD", dod_new, "RdBu", v), ("correction added (mm)", D, "Blues", 150))):
    im = ax.imshow(np.where(fp, field, np.nan)[r0:r1, c0:c1], cmap=cm, vmin=(-lim if cm == "RdBu" else 0), vmax=lim, origin="lower", interpolation="nearest")
    ax.set_title(f"{lab}, floodplain mask (mm)", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
fig.suptitle("Flat-floodplain vegetation correction from the 3x3 logistic, erosion side only, tapered over the stable NMAD", fontsize=11)
fig.tight_layout(); fig.savefig(f"{A.figdir}/floodplain_vegcorr_maps.png", dpi=110); plt.close(fig)
fig, ax = plt.subplots(figsize=(8, 5))
bi = np.minimum((np.argsort(np.argsort(span)) * 20) // span.size, 19)
for lab, v_, col in (("raw", raw, "#2a78d6"), ("corrected", raw + delta, "#eb6834")):
    xm = [np.median(span[bi == b]) for b in range(20)]
    ax.plot(xm, [np.median(v_[bi == b]) for b in range(20)], color=col, lw=2, marker="o", ms=4, label=lab)
    ax.fill_between(xm, [np.percentile(v_[bi == b], 25) for b in range(20)], [np.percentile(v_[bi == b], 75) for b in range(20)], color=col, alpha=0.15, lw=0)
ax.axhline(0, color="#52514e", lw=0.8, ls=":"); ax.set_xlabel("gen1 pooled return span (mm)"); ax.set_ylabel("DoD on flat floodplain (mm), median and 25-75% band")
ax.legend(frameon=False); [ax.spines[s_].set_visible(False) for s_ in ("top", "right")]
fig.tight_layout(); fig.savefig(f"{A.figdir}/floodplain_vegcorr_vs_span.png", dpi=140)
print(f"\nwrote {S}/dod_floodplain_vegcorr.npy, floodplain_vegcorr_delta.npy, floodplain_vegcorr.json; figures/floodplain_vegcorr_maps.png, _vs_span.png")
