#!/usr/bin/env python3
"""Fit a rigid translation (dE, dN, dZ) of the UNTIED hillslope RTK points (B2 base, 25 Oct 2026)
that minimizes their misfit to the 2021 lidar ground (tile 4342-29-64 z_after.npy, 5 m, bilinear).
The B2 tie is left out on purpose: the hillslope then estimates the whole translation on its own,
and the pole-pair tie (dE +0.802, dN -1.123, dZ -0.473 m) is the independent check.
Residual r_i = (navd88_untied_i + dZ) - Z2021(E_i + dE, N_i + dN). For each (dE, dN) on a grid the
best dZ is closed-form: median(r) for the L1 misfit, mean(r) after a 3-NMAD clip for the RMS one.
Grid +/-2.5 m in 0.05 m steps (mine). Uncertainty by bootstrap over points (500 draws, mine).
Writes data/derived/4342-29-64/bluff_check/translation_fit.csv (the grid) and
figures/rtk_2026-10-25_hillslope_translation.png.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates

T = "4342-29-64"; b = json.load(open(f"data/derived/{T}/corrections.json"))["bounds"]; res = 5.0
Z = np.load(f"data/derived/{T}/z_after.npy")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
H = P[P.status.str.startswith("FIX") & P.group.str.startswith("hillslope")].copy()
E, N, h = H.easting.to_numpy(), H.northing.to_numpy(), H.navd88_m.to_numpy()   # UNTIED
print(f"points: {len(H)} hillslope FIX, untied (B2 as exported)")
PAIR = (0.802, -1.123, -0.473)

def z2021(x, y):
    return map_coordinates(Z, [(y - b[1]) / res - 0.5, (x - b[0]) / res - 0.5], order=1, mode="nearest")

def misfits(dE, dN, idx=slice(None)):
    r = h[idx] - z2021(E[idx] + dE, N[idx] + dN)
    dz1 = -np.median(r); l1 = np.mean(np.abs(r + dz1))
    rr = r.copy()
    for _ in range(5):
        m = np.mean(rr); nm = 1.4826 * np.median(np.abs(rr - np.median(rr))); keep = np.abs(rr - m) < 3 * max(nm, 1e-3)
        if keep.all(): break
        rr = rr[keep]
    dz2 = -np.mean(rr); rms = np.sqrt(np.mean((rr + dz2) ** 2)); nmad = 1.4826 * np.median(np.abs(r + dz1))
    return dz1, l1, dz2, rms, nmad, len(rr)

g = np.arange(-2.5, 2.5001, 0.05); rows = []
for dE in g:
    for dN in g:
        dz1, l1, dz2, rms, nmad, n = misfits(dE, dN); rows.append((dE, dN, dz1, l1, dz2, rms, nmad, n))
G = pd.DataFrame(rows, columns=["dE", "dN", "dZ_L1", "L1_m", "dZ_rms", "rms_m", "nmad_m", "n_kept"])
os.makedirs(f"data/derived/{T}/bluff_check", exist_ok=True); G.to_csv(f"data/derived/{T}/bluff_check/translation_fit.csv", index=False)
best1 = G.loc[G.L1_m.idxmin()]; best2 = G.loc[G.rms_m.idxmin()]
print(f"L1 minimum : dE {best1.dE:+.2f} dN {best1.dN:+.2f} dZ {best1.dZ_L1:+.3f} m, mean|r| {best1.L1_m*1000:.0f} mm, NMAD {best1.nmad_m*1000:.0f} mm")
print(f"RMS minimum: dE {best2.dE:+.2f} dN {best2.dN:+.2f} dZ {best2.dZ_rms:+.3f} m, RMS {best2.rms_m*1000:.0f} mm after 3-NMAD clip (n kept {int(best2.n_kept)})")
p = misfits(*PAIR[:2]); print(f"at the pair tie (dE {PAIR[0]:+.3f}, dN {PAIR[1]:+.3f}): best dZ_L1 {p[0]:+.3f} (pair says {PAIR[2]:+.3f}), mean|r| {p[1]*1000:.0f} mm, NMAD {p[4]*1000:.0f} mm")
z0 = misfits(0, 0); print(f"no shift: dZ_L1 {z0[0]:+.3f}, mean|r| {z0[1]*1000:.0f} mm, NMAD {z0[4]*1000:.0f} mm")

# bootstrap over points: refit on a coarser local grid around the L1 minimum
rng = np.random.default_rng(0); gl = np.arange(-1.0, 1.0001, 0.1); boots = []
for k in range(500):
    idx = rng.integers(0, len(H), len(H)); bestv = (np.inf, None)
    for dE in best1.dE + gl:
        for dN in best1.dN + gl:
            dz1, l1, *_ = misfits(dE, dN, idx)
            if l1 < bestv[0]: bestv = (l1, (dE, dN, dz1))
    boots.append(bestv[1])
B = np.array(boots); q = lambda a: np.percentile(a, [16, 50, 84])
print("bootstrap (500 draws, L1, 0.1 m local grid) p16/p50/p84: dE", q(B[:, 0]).round(2), " dN", q(B[:, 1]).round(2), " dZ", q(B[:, 2]).round(3))
print("bootstrap corr(dE, dZ) =", np.corrcoef(B[:, 0], B[:, 2])[0, 1].round(2), " corr(dN, dZ) =", np.corrcoef(B[:, 1], B[:, 2])[0, 1].round(2))
# the trade-off slope: dZ per m of dE along the valley floor
sel = G[np.isclose(G.dN, best1.dN)]; c = np.polyfit(sel.dE, sel.dZ_L1, 1); print(f"along dN={best1.dN:+.2f}: dZ_L1 changes by {c[0]*1000:.0f} mm per m of dE (linear fit over the +/-2.5 m row)")

fig, axs = plt.subplots(1, 3, figsize=(16, 5))
M = G.pivot(index="dN", columns="dE", values="L1_m") * 1000
im = axs[0].pcolormesh(M.columns, M.index, M.values, cmap="viridis", vmin=M.values.min(), vmax=min(M.values.max(), M.values.min() * 4), shading="nearest")
axs[0].plot(best1.dE, best1.dN, "r*", ms=14, label=f"L1 min ({best1.dE:+.2f}, {best1.dN:+.2f})"); axs[0].plot(PAIR[0], PAIR[1], "wo", mec="k", ms=9, label="pole-pair tie (+0.80, -1.12)"); axs[0].plot(0, 0, "k+", ms=12, mew=2, label="no shift")
axs[0].set_xlabel("dE (m)"); axs[0].set_ylabel("dN (m)"); axs[0].set_title("mean |residual| (mm) vs horizontal shift, dZ optimal at each", fontsize=9); axs[0].legend(fontsize=8, loc="lower left"); fig.colorbar(im, ax=axs[0], shrink=0.8)
axs[1].plot(sel.dE, sel.L1_m * 1000, "k-", label="mean |r|"); axs[1].plot(sel.dE, sel.nmad_m * 1000, "b--", label="NMAD"); axs[1].axvline(PAIR[0], color="gray", ls=":", label="pair dE"); axs[1].set_xlabel("dE (m)"); axs[1].set_ylabel("mm"); axs[1].set_title(f"misfit along dN = {best1.dN:+.2f}", fontsize=9); axs[1].legend(fontsize=8); axs[1].grid(alpha=0.3)
axs[2].plot(sel.dE, sel.dZ_L1, "k-"); axs[2].axvline(PAIR[0], color="gray", ls=":"); axs[2].axhline(PAIR[2], color="gray", ls=":"); axs[2].plot(PAIR[0], PAIR[2], "wo", mec="k", ms=9, label="pair tie"); axs[2].plot(best1.dE, best1.dZ_L1, "r*", ms=14, label="L1 min")
axs[2].scatter(B[:, 0], B[:, 2], s=6, c="orange", alpha=0.5, label="bootstrap"); axs[2].set_xlabel("dE (m)"); axs[2].set_ylabel("best dZ (m)"); axs[2].set_title("the dE-dZ trade-off on an east-facing bluff", fontsize=9); axs[2].legend(fontsize=8); axs[2].grid(alpha=0.3)
fig.suptitle("Hillslope RTK (39 FIX, B2 untied) vs 2021 lidar ground: rigid translation fit", fontsize=10); fig.tight_layout()
fig.savefig("figures/rtk_2026-10-25_hillslope_translation.png", dpi=120); print("figures/rtk_2026-10-25_hillslope_translation.png")
