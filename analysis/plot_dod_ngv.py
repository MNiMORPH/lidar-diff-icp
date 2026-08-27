#!/usr/bin/env python3
"""Plot the DoD before and after the NGV vegetation correction, and the correction itself."""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--out", default="figures/dod_ngv_corrected.png")
ap.add_argument("--clip", type=float, default=200.0, help="colour limit, mm")
a = ap.parse_args()

T = a.tile
dod = np.load(os.path.join(T, "dod.npy")) * 1000.0
corr = np.load(os.path.join(T, "dod_ngv.npy")) * 1000.0
dz = np.load(os.path.join(T, "dod_ngv_dz.npy")) * 1000.0
cfg = json.load(open(os.path.join(T, "corrections.json")))
X0, Y0, X1, Y1 = cfg["bounds"]
ext = (X0, X1, Y0, Y1)

fig, ax = plt.subplots(1, 3, figsize=(16.5, 7.2), sharex=True, sharey=True)
for A, D, ttl, cmap, lim in (
        (ax[0], dod, "DoD before", "RdBu_r", a.clip),
        (ax[1], corr, "DoD after NGV correction", "RdBu_r", a.clip),
        (ax[2], dz, "the correction applied  ($-325.2\\,$mm $\\times$ NGV)", "viridis", None)):
    if lim is not None:
        im = A.imshow(D, origin="lower", extent=ext, cmap=cmap, vmin=-lim, vmax=lim,
                      interpolation="nearest")
    else:
        im = A.imshow(D, origin="lower", extent=ext, cmap=cmap,
                      vmin=np.nanpercentile(dz, 0.5), vmax=0, interpolation="nearest")
    A.set_title(ttl, fontsize=11)
    A.set_xlabel("easting (m)")
    fig.colorbar(im, ax=A, fraction=0.046, pad=0.03, label="mm")
ax[0].set_ylabel("northing (m)")
fin = np.isfinite(dod) & np.isfinite(corr)
fig.suptitle(
    f"Elba DoD (gen2 $-$ gen1, $+$ve = deposition), before and after correcting gen2 for "
    f"leaf-on vegetation\nmedian {np.median(dod[fin]):+.1f} $\\to$ {np.median(corr[fin]):+.1f} mm   "
    f"|   correction median {np.median(dz[fin]):+.1f} mm, most negative {np.nanmin(dz):+.1f} mm",
    fontsize=11.5)
fig.tight_layout()
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
fig.savefig(a.out, dpi=150)
print(f"  DoD before : median {np.median(dod[fin]):+7.1f}  NMAD {1.4826*np.median(np.abs(dod[fin]-np.median(dod[fin]))):7.1f} mm")
print(f"  DoD after  : median {np.median(corr[fin]):+7.1f}  NMAD {1.4826*np.median(np.abs(corr[fin]-np.median(corr[fin]))):7.1f} mm")
print(f"  correction : median {np.median(dz[fin]):+7.1f}  min {np.nanmin(dz):+7.1f}  max {np.nanmax(dz):+7.1f} mm")
print(f"wrote {a.out}")
