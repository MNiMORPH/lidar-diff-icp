#!/usr/bin/env python3
"""A rigid shift in all three axes (Andy, 2026-11-11: "I am asking E, N, and z"): the (dE, dN, dZ)
added to the TIED RTK points that best fits the 2021 lidar ground (tile 4342-29-64 z_after.npy,
5 m, bilinear, the surface the table's gen2_minus_rtk_m reads).
  r_i(dE, dN, dZ) = Z2021(E_i + dE, N_i + dN) - (navd88_tied_i + dZ)       (lidar minus shifted RTK, m)
  criteria: RMS(r) (the fit), with mean|r| printed beside it
Full grid search, nothing profiled out: dE, dN in +/-1.5 m by 0.05, dZ in -0.2 .. +0.6 m by 0.01
(grid is mine). Nothing is cut. Fits: hillslope both groups, every FIX point, RS3 alone, T4a alone.
Writes data/derived/4342-29-64/bluff_check/shift_3d_best.csv and, for the first two fits,
figures/rtk_2026-10-25_shift_3d.png: RMS over (dE, dZ) at the best dN, (dN, dZ) at the best dE,
(dE, dN) at the best dZ; star = the 3-D minimum, circle = no shift.
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
H = P[P.status.str.startswith("FIX") & P.gen2_minus_rtk_m.notna()].copy()
FITS = {"hillslope, both groups": H.group.str.startswith("hillslope"), "all FIX points": H.group.notna(), "RS3 (T3a)": H.group == "hillslope_RS3", "T4a": H.group == "hillslope_T4a"}
gE = np.round(np.arange(-1.5, 1.5001, 0.05), 2); gZ = np.round(np.arange(-0.2, 0.6001, 0.01), 2)

def z2021(x, y):
    return map_coordinates(Z, [(y - b[1]) / res - 0.5, (x - b[0]) / res - 0.5], order=1, mode="nearest")

out = []; cubes = {}
for name, sel in FITS.items():
    s = H[sel]; E, N, h = s.easting_tied.to_numpy(), s.northing_tied.to_numpy(), s.navd88_tied_m.to_numpy()
    RMS = np.empty((len(gE), len(gE), len(gZ))); L1 = np.empty_like(RMS)
    for i, dE in enumerate(gE):
        for j, dN in enumerate(gE):
            r0 = z2021(E + dE, N + dN) - h                      # lidar minus RTK before the vertical shift
            r = r0[:, None] - gZ[None, :]                       # minus dZ
            RMS[i, j] = np.sqrt(np.mean(r ** 2, axis=0)); L1[i, j] = np.mean(np.abs(r), axis=0)
    cubes[name] = RMS
    k = np.unravel_index(RMS.argmin(), RMS.shape); k1 = np.unravel_index(L1.argmin(), L1.shape); z0 = (np.argmin(np.abs(gE)), np.argmin(np.abs(gE)), np.argmin(np.abs(gZ)))
    print(f"\n{name}: n {len(s)}")
    print(f"  no shift:                 RMS {1000 * RMS[z0]:.0f} mm, mean|r| {1000 * L1[z0]:.0f} mm")
    print(f"  dZ only (dE = dN = 0):    dZ {gZ[RMS[z0[0], z0[1]].argmin()]:+.2f} m -> RMS {1000 * RMS[z0[0], z0[1]].min():.0f} mm")
    print(f"  3-D RMS minimum:          dE {gE[k[0]]:+.2f} dN {gE[k[1]]:+.2f} dZ {gZ[k[2]]:+.2f} m -> RMS {1000 * RMS[k]:.0f} mm, mean|r| {1000 * L1[k]:.0f} mm")
    print(f"  3-D mean|r| minimum:      dE {gE[k1[0]]:+.2f} dN {gE[k1[1]]:+.2f} dZ {gZ[k1[2]]:+.2f} m -> mean|r| {1000 * L1[k1]:.0f} mm, RMS {1000 * RMS[k1]:.0f} mm")
    # how sharp: the dZ range over which RMS stays within 5 % of its minimum along the dE-dZ ridge (the trade-off), at the best dN
    sl = RMS[:, k[1], :]; near = np.argwhere(sl <= 1.05 * RMS[k])
    print(f"  within 5 % of the minimum at the best dN: dE {gE[near[:, 0].min()]:+.2f} .. {gE[near[:, 0].max()]:+.2f}, dZ {gZ[near[:, 1].min()]:+.2f} .. {gZ[near[:, 1].max()]:+.2f} m")
    out.append(dict(fit=name, n=len(s), dE=gE[k[0]], dN=gE[k[1]], dZ=gZ[k[2]], rms_m=RMS[k], l1_m=L1[k], rms0_m=RMS[z0], l1_0_m=L1[z0], dZ_only=gZ[RMS[z0[0], z0[1]].argmin()], rms_dZ_only=RMS[z0[0], z0[1]].min()))
os.makedirs(f"data/derived/{T}/bluff_check", exist_ok=True); pd.DataFrame(out).to_csv(f"data/derived/{T}/bluff_check/shift_3d_best.csv", index=False)

fig, axs = plt.subplots(2, 3, figsize=(14, 8))
for row, name in enumerate(list(FITS)[:2]):
    RMS = cubes[name] * 1000; k = np.unravel_index(RMS.argmin(), RMS.shape); vmax = min(RMS[k] * 3, 600)
    panels = ((RMS[:, k[1], :].T, gE, gZ, "dE (m)", "dZ (m)", (gE[k[0]], gZ[k[2]])), (RMS[k[0], :, :].T, gE, gZ, "dN (m)", "dZ (m)", (gE[k[1]], gZ[k[2]])), (RMS[:, :, k[2]].T, gE, gE, "dE (m)", "dN (m)", (gE[k[0]], gE[k[1]])))
    for ax, (A, gx, gy, xl, yl, best) in zip(axs[row], panels):
        im = ax.imshow(A, origin="lower", extent=[gx[0], gx[-1], gy[0], gy[-1]], aspect="auto", cmap="viridis", vmin=RMS[k], vmax=vmax)
        ax.plot(*best, "w*", ms=12, mec="k"); ax.plot(0, 0, "wo", ms=6, mec="k"); ax.set_xlabel(xl); ax.set_ylabel(yl)
    axs[row, 0].text(0.03, 0.95, name, transform=axs[row, 0].transAxes, color="white", va="top", fontsize=9)
    fig.colorbar(im, ax=axs[row], shrink=0.9, pad=0.02).set_label("RMS of lidar minus shifted RTK (mm)")
fig.savefig("figures/rtk_2026-10-25_shift_3d.png", dpi=120, bbox_inches="tight"); print("\nfigures/rtk_2026-10-25_shift_3d.png")
