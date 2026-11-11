#!/usr/bin/env python3
"""One or more translations (Andy, 2026-11-11): for each hillslope group, for both pooled, and for every FIX point, the
lateral shift (dE, dN) of the TIED RTK points that minimizes the spread of lidar minus RTK against
the 2021 lidar ground (tile 4342-29-64 z_after.npy, 5 m, bilinear, the same surface the table's
gen2_minus_rtk_m reads), and the vertical offset of the lidar at that shift.
  r_i(dE, dN) = Z2021(E_i + dE, N_i + dN) - navd88_tied_i        (lidar minus RTK, m)
  spread      = NMAD(r)  (the criterion that picks dE, dN; the mean-square after a 3-NMAD clip is
                printed beside it)
  offset      = median(r) at the chosen (dE, dN)                  (what dZ would remove)
Shifts are relative to the tied positions, i.e. residual to the DR9925 tie. Grid +/-1.5 m in 0.05 m
steps (mine). Uncertainty: bootstrap over points, 500 draws (mine). Nothing is cut: every FIX
hillslope point enters, the -1.6 m point included, which is why the NMAD and not the RMS chooses.
Writes data/derived/4342-29-64/bluff_check/offsets_by_group.csv (the three grids) and
figures/rtk_2026-10-25_hillslope_offsets_by_group.png (NMAD over dE, dN; one panel per fit).
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
H = P[P.status.str.startswith("FIX") & P.gen2_minus_rtk_m.notna()].copy()   # every FIX point; the fits select
FITS = {"RS3 (T3a)": H.group == "hillslope_RS3", "T4a": H.group == "hillslope_T4a", "both": H.group.str.startswith("hillslope"), "all FIX points": H.group.notna()}

def z2021(x, y):
    return map_coordinates(Z, [(y - b[1]) / res - 0.5, (x - b[0]) / res - 0.5], order=1, mode="nearest")

def stats(E, N, h, dE, dN):
    r = z2021(E + dE, N + dN) - h; med = np.median(r); nmad = 1.4826 * np.median(np.abs(r - med)); l1 = np.mean(np.abs(r - med))
    keep = np.abs(r - med) < 3 * max(nmad, 1e-3); rms = np.sqrt(np.mean((r[keep] - r[keep].mean()) ** 2))
    return med, nmad, rms, int(keep.sum()), l1

g = np.round(np.arange(-1.5, 1.5001, 0.05), 2); rows = []; out = []
for name, sel in FITS.items():
    s = H[sel]; E, N, h = s.easting_tied.to_numpy(), s.northing_tied.to_numpy(), s.navd88_tied_m.to_numpy()
    for dE in g:
        for dN in g:
            med, nmad, rms, n, l1 = stats(E, N, h, dE, dN); rows.append((name, dE, dN, med, nmad, rms, n, l1))
    G = pd.DataFrame(rows, columns=["fit", "dE", "dN", "offset_m", "nmad_m", "rms_clip_m", "n_kept", "l1_m"]); G = G[G.fit == name]
    k = G.nmad_m.idxmin(); best = G.loc[k]; k2 = G.rms_clip_m.idxmin(); best2 = G.loc[k2]; zero = G[(G.dE == 0) & (G.dN == 0)].iloc[0]
    rng = np.random.default_rng(0); B = []
    for _ in range(500):
        i = rng.integers(0, len(s), len(s)); loc = []
        for dE in np.arange(best.dE - 0.3, best.dE + 0.3001, 0.1):
            for dN in np.arange(best.dN - 0.3, best.dN + 0.3001, 0.1):
                med, nmad, rms, n, l1 = stats(E[i], N[i], h[i], dE, dN); loc.append((dE, dN, med, nmad))
        loc = np.array(loc); B.append(loc[loc[:, 3].argmin()])
    B = np.array(B); p = lambda j: np.percentile(B[:, j], [16, 50, 84])
    print(f"\n{name}: n {len(s)}")
    print(f"  no shift (tied as is):      offset {1000 * zero.offset_m:+.0f} mm, NMAD {1000 * zero.nmad_m:.0f}, clipped rms {1000 * zero.rms_clip_m:.0f} mm (n kept {zero.n_kept})")
    print(f"  min-NMAD shift:  dE {best.dE:+.2f} dN {best.dN:+.2f} m -> offset {1000 * best.offset_m:+.0f} mm, NMAD {1000 * best.nmad_m:.0f}, clipped rms {1000 * best.rms_clip_m:.0f} mm (n kept {best.n_kept})")
    print(f"  min-clipped-rms shift: dE {best2.dE:+.2f} dN {best2.dN:+.2f} m -> offset {1000 * best2.offset_m:+.0f} mm, NMAD {1000 * best2.nmad_m:.0f}, clipped rms {1000 * best2.rms_clip_m:.0f} mm")
    best3 = G.loc[G.l1_m.idxmin()]
    print(f"  BEST FIT, min mean|r| with the vertical closed-form (the combined shift): dE {best3.dE:+.2f} dN {best3.dN:+.2f} m, lidar minus RTK {1000 * best3.offset_m:+.0f} mm (i.e. RTK + {1000 * best3.offset_m:+.0f} mm matches the lidar there); mean|r| {1000 * zero.l1_m:.0f} -> {1000 * best3.l1_m:.0f} mm, NMAD {1000 * best3.nmad_m:.0f}")
    print(f"  bootstrap p16/p50/p84 (min-NMAD, 0.1 m local grid): dE {p(0).round(2)}  dN {p(1).round(2)}  offset mm {(1000 * p(2)).round(0)}")
    out.append(dict(fit=name, n=len(s), dE=best.dE, dN=best.dN, offset_m=best.offset_m, nmad_m=best.nmad_m, rms_clip_m=best.rms_clip_m, nmad0_m=zero.nmad_m, offset0_m=zero.offset_m))
G = pd.DataFrame(rows, columns=["fit", "dE", "dN", "offset_m", "nmad_m", "rms_clip_m", "n_kept", "l1_m"])
os.makedirs(f"data/derived/{T}/bluff_check", exist_ok=True); G.to_csv(f"data/derived/{T}/bluff_check/offsets_by_group.csv", index=False)
pd.DataFrame(out).to_csv(f"data/derived/{T}/bluff_check/offsets_by_group_best.csv", index=False)

fig, axs = plt.subplots(1, 4, figsize=(17, 4.2), sharey=True)
for ax, (name, sel) in zip(axs, FITS.items()):
    Gi = G[G.fit == name].pivot(index="dN", columns="dE", values="nmad_m") * 1000
    im = ax.imshow(Gi.values, origin="lower", extent=[g[0], g[-1], g[0], g[-1]], cmap="viridis", vmin=0, vmax=300)
    bi = [o for o in out if o["fit"] == name][0]; ax.plot(bi["dE"], bi["dN"], "w*", ms=12, mec="k"); ax.plot(0, 0, "wo", ms=6, mec="k")
    ax.set_xlabel("dE (m)"); ax.text(0.03, 0.95, name, transform=ax.transAxes, color="white", va="top", fontsize=9)
axs[0].set_ylabel("dN (m)"); cb = fig.colorbar(im, ax=axs, shrink=0.85, pad=0.02); cb.set_label("NMAD of lidar minus RTK (mm)")
fig.savefig("figures/rtk_2026-10-25_hillslope_offsets_by_group.png", dpi=120, bbox_inches="tight"); print("\nfigures/rtk_2026-10-25_hillslope_offsets_by_group.png")
