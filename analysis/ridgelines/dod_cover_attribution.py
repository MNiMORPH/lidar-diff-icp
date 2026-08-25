#!/usr/bin/env python3
"""How much of the DoD's apparent hillslope AGGRADATION is the canopy measurement effect?

Earlier DoDs showed net elevation increase across the steeper hillslopes, including their
bottoms. Two candidate explanations were never separated: real deposition, or the fact that
the leaf-on 2021 epoch reads the forest floor high and the forest lives on the slopes.

cover_offset_reference.py calibrates the second one on ground where the first cannot
operate (low-gradient upland, TPI>0): offset = a + k*cover, with k ~ -49 mm per unit cover
on both tiles. In DoD sign (gen2 - gen1, + = gain) that is an apparent GAIN of |k|*cover.
This subtracts that predicted gain and asks whether the slope and hillslope-position
signals survive. Whatever remains is the part canopy cannot explain.

The calibration is deliberately fitted OFF the slopes it is used to correct, so this is a
prediction being tested, not a fit being re-reported.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/dod_cover_attribution.py --tile data/derived/elba_fulldensity
"""
import argparse, json, os
import numpy as np
from scipy.ndimage import uniform_filter, distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--dod", default="dod.npy")
ap.add_argument("--k-cover", type=float, default=None,
                help="mm of apparent DoD GAIN per unit cover; default = the tile's own "
                     "low-gradient-upland calibration (49.6 elba, 48.4 elbaext)")
ap.add_argument("--tpi-window", type=float, default=500.0)
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
K = A.k_cover if A.k_cover is not None else (48.4 if TILE == "elbaext" else 49.6)


def grid_meta(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            return b[0], b[1], r
    raise SystemExit(f"no grid meta in {tile}")


X0, Y0, RES = grid_meta(A.tile)
dod = np.load(f"{A.tile}/{A.dod}") * 1000.0                  # m -> mm, + = gain
slope = np.load(f"{A.tile}/slope.npy")
cover = np.load(f"{A.tile}/canopy_cover_pfs.npy")
z = np.load(f"{A.tile}/z_after.npy")
zf = z.copy(); miss = ~np.isfinite(zf)
if miss.any():
    zf = zf[tuple(distance_transform_edt(miss, return_distances=False, return_indices=True))]
k = max(3, int(round(A.tpi_window / RES)) | 1)
tpi = zf - uniform_filter(zf, size=k, mode="nearest")

fin = np.isfinite(dod) & np.isfinite(slope) & np.isfinite(cover)
pred = K * np.where(np.isfinite(cover), cover, 0.0)          # predicted canopy-driven gain
resid = dod - pred
print("=" * 92)
print(f"DoD ATTRIBUTION  [{TILE}, {A.dod}]   canopy model: apparent gain = {K:.1f} mm x cover")
print(f"{fin.sum():,} cells;  median DoD {np.median(dod[fin]):+.1f} mm  ->  "
      f"after removing the canopy term {np.median(resid[fin]):+.1f} mm")
print("=" * 92)

SB = [(0, 3), (3, 6), (6, 10), (10, 15), (15, 20), (20, 25), (25, 30), (30, 60)]
print(f"\n{'slope band':>12s} {'n':>9s} {'med cover':>10s} {'DoD (mm)':>10s} "
      f"{'canopy pred':>12s} {'RESIDUAL':>10s}")
for lo, hi in SB:
    m = fin & (slope >= lo) & (slope < hi)
    if m.sum() < 200: print(f"{lo:5.0f}-{hi:<5.0f} {m.sum():>9,}  (sparse)"); continue
    print(f"{lo:5.0f}-{hi:<5.0f} {m.sum():>9,} {np.median(cover[m]):>10.2f} "
          f"{np.median(dod[m]):>+10.1f} {np.median(pred[m]):>+12.1f} {np.median(resid[m]):>+10.1f}")

print(f"\nby HILLSLOPE POSITION (TPI {A.tpi_window:g} m quintile, low = valley bottom):")
qs = np.percentile(tpi[fin], [0, 20, 40, 60, 80, 100])
print(f"{'TPI band (m)':>18s} {'n':>9s} {'med slope':>10s} {'med cover':>10s} "
      f"{'DoD (mm)':>10s} {'RESIDUAL':>10s}")
for lo, hi in zip(qs[:-1], qs[1:]):
    m = fin & (tpi >= lo) & (tpi < hi)
    if m.sum() < 200: continue
    print(f"{lo:>+8.1f}..{hi:>+7.1f} {m.sum():>9,} {np.median(slope[m]):>10.1f} "
          f"{np.median(cover[m]):>10.2f} {np.median(dod[m]):>+10.1f} {np.median(resid[m]):>+10.1f}")

fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2), dpi=130)
cen = [0.5*(a+b) for a, b in SB]
for arr, nm, col in ((dod, "DoD as measured", "C3"), (resid, "residual after canopy term", "C0")):
    med = [np.median(arr[fin & (slope >= a) & (slope < b)])
           if (fin & (slope >= a) & (slope < b)).sum() >= 200 else np.nan for a, b in SB]
    ax[0].plot(cen, med, "o-", color=col, lw=1.6, ms=5, label=nm)
ax[0].axhline(0, color="k", lw=.7); ax[0].set_xlabel("surface slope (deg)")
ax[0].set_ylabel("median DoD (mm; + = gain)")
ax[0].set_title("apparent hillslope gain, before and after\nremoving the canopy term", fontsize=10)
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
qc = [0.5*(a+b) for a, b in zip(qs[:-1], qs[1:])]
for arr, nm, col in ((dod, "DoD as measured", "C3"), (resid, "residual after canopy term", "C0")):
    med = [np.median(arr[fin & (tpi >= a) & (tpi < b)]) for a, b in zip(qs[:-1], qs[1:])]
    ax[1].plot(qc, med, "o-", color=col, lw=1.6, ms=5, label=nm)
ax[1].axhline(0, color="k", lw=.7)
ax[1].set_xlabel(f"hillslope position: TPI over {A.tpi_window:g} m (m; - = valley bottom)")
ax[1].set_ylabel("median DoD (mm)")
ax[1].set_title("gain by hillslope position", fontsize=10); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
fig.suptitle(f"how much of the DoD is canopy? — {TILE} (canopy term = {K:.1f} mm x cover)", y=1.02)
out = f"figures/refdatum/dod_cover_attribution_{TILE}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")
