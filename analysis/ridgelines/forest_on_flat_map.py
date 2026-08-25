#!/usr/bin/env python3
"""WHERE is the forest that sits on gentle ground?

The slope x cover analyses are confounded because steep ground is preferentially forested.
The converse matters just as much: the forest that does occupy FLAT ground is what carries
the "cover effect at low slope", so if it sits systematically at hilltops (or in valley
bottoms) then that term is entangled with hillslope position, not with cover alone.

This maps the forest-on-gentle-ground cells and scores their hillslope position with the
topographic position index (TPI = elevation minus the mean over a disc), at two scales:
positive = local high (crest, plateau), negative = local low (valley floor, hollow). The
test is not "where are they" by eye but whether their TPI distribution differs from that of
ALL gentle cells -- the honest null being that flat forest is spread like flat ground is.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/forest_on_flat_map.py --tile data/derived/elbaext
"""
import argparse, json, os
import numpy as np
from scipy.ndimage import uniform_filter, distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from lidar_diff_icp.viz import hillshade

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext")
ap.add_argument("--slope-max", type=float, default=12.0, help="gentle-ground cut (deg)")
ap.add_argument("--cover-min", type=float, default=0.50, help="forest cut (PFS cover fraction)")
ap.add_argument("--tpi-small", type=float, default=150.0, help="small TPI window (m)")
ap.add_argument("--tpi-large", type=float, default=500.0, help="large TPI window (m)")
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))


def grid_meta(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            return b[0], b[1], r
    raise SystemExit(f"no grid meta in {tile}")


X0, Y0, RES = grid_meta(A.tile)
z = np.load(f"{A.tile}/z_after.npy")
slope = np.load(f"{A.tile}/slope.npy")
cover = np.load(f"{A.tile}/canopy_cover_pfs.npy")

zf = z.copy(); miss = ~np.isfinite(zf)
if miss.any():
    zf = zf[tuple(distance_transform_edt(miss, return_distances=False, return_indices=True))]


def tpi(win_m):
    """Elevation minus the mean over a square window of `win_m` (NaN-safe via the fill)."""
    k = max(3, int(round(win_m / RES)) | 1)
    return zf - uniform_filter(zf, size=k, mode="nearest")


tpi_s, tpi_l = tpi(A.tpi_small), tpi(A.tpi_large)
fin = np.isfinite(z) & np.isfinite(slope) & np.isfinite(cover)
gentle = fin & (slope < A.slope_max)
forest = fin & (cover > A.cover_min)
ff = gentle & forest                                   # forest on gentle ground

print("=" * 84)
print(f"FOREST ON GENTLE GROUND  [{TILE}: slope < {A.slope_max:g} deg, cover > {A.cover_min:g}]")
print("=" * 84)
print(f"  finite cells      {fin.sum():>9,}")
print(f"  gentle (<{A.slope_max:g} deg)  {gentle.sum():>9,}  ({100*gentle.sum()/fin.sum():4.1f}% of tile)")
print(f"  forest            {forest.sum():>9,}  ({100*forest.sum()/fin.sum():4.1f}%)")
print(f"  forest AND gentle {ff.sum():>9,}  ({100*ff.sum()/gentle.sum():4.1f}% of gentle ground,"
      f" {100*ff.sum()/max(forest.sum(),1):4.1f}% of all forest)")

for name, T in (("TPI %g m" % A.tpi_small, tpi_s), ("TPI %g m" % A.tpi_large, tpi_l)):
    a, b = T[ff], T[gentle]
    print(f"\n  {name}:  forest-on-gentle   median {np.median(a):+6.2f} m   "
          f"IQR {np.percentile(a,25):+.2f}..{np.percentile(a,75):+.2f}")
    print(f"  {'':>{len(name)}}   ALL gentle ground  median {np.median(b):+6.2f} m   "
          f"IQR {np.percentile(b,25):+.2f}..{np.percentile(b,75):+.2f}")
    print(f"  {'':>{len(name)}}   shift {np.median(a)-np.median(b):+.2f} m "
          f"({'HIGH: crests/plateaux' if np.median(a) > np.median(b) else 'LOW: valley bottoms'})")
    # the cleanest statement: P(forest | gentle) as a function of position
    qs = np.percentile(b, [0, 20, 40, 60, 80, 100])
    print(f"  {'':>{len(name)}}   P(forest | gentle) by {name} quintile (low -> high position):")
    for lo, hi in zip(qs[:-1], qs[1:]):
        m = gentle & (T >= lo) & (T < hi)
        print(f"  {'':>{len(name)}}      {lo:+6.2f}..{hi:+6.2f} m : "
              f"{100*(m & ff).sum()/max(m.sum(),1):5.2f}%   (n={m.sum():,})")

hs = hillshade(zf, RES, X0, Y0)
fig, ax = plt.subplots(1, 3, figsize=(19, 6.2), dpi=130)
ax[0].imshow(hs, origin="lower", cmap="gray", vmin=0, vmax=1)
ov = np.full(z.shape + (4,), 0.0)
ov[gentle & ~forest] = (0.30, 0.55, 0.95, 0.30)        # gentle, not forest
ov[forest & ~gentle] = (0.95, 0.60, 0.15, 0.35)        # forest, steep
ov[ff] = (0.05, 0.65, 0.10, 0.95)                      # the subject
ax[0].imshow(ov, origin="lower")
ax[0].set_title(f"green = forest on slope < {A.slope_max:g}°\n"
                f"blue = gentle & open, orange = forest & steep", fontsize=10)
ax[0].set_xticks([]); ax[0].set_yticks([])

for axi, (T, nm) in zip(ax[1:], ((tpi_s, f"{A.tpi_small:g} m"), (tpi_l, f"{A.tpi_large:g} m"))):
    lo, hi = np.percentile(T[gentle], [0.5, 99.5])
    bins = np.linspace(lo, hi, 60)
    axi.hist(T[gentle], bins=bins, density=True, histtype="step", lw=2, color="0.4",
             label=f"all gentle ground (n={gentle.sum():,})")
    axi.hist(T[ff], bins=bins, density=True, histtype="step", lw=2, color="C2",
             label=f"forest on gentle ground (n={ff.sum():,})")
    axi.axvline(0, color="k", lw=.7)
    axi.axvline(np.median(T[gentle]), color="0.4", ls=":", lw=1.2)
    axi.axvline(np.median(T[ff]), color="C2", ls=":", lw=1.2)
    axi.set_xlabel(f"TPI over {nm}  (m; + = local high, - = local low)")
    axi.set_ylabel("density"); axi.legend(fontsize=8); axi.grid(alpha=.3)
    axi.set_title(f"hillslope position of flat forest ({nm})", fontsize=10)
fig.suptitle(f"Where is the forest on gentle ground? — {TILE} "
             f"(slope < {A.slope_max:g}°, PFS cover > {A.cover_min:g})", y=1.02)
out = f"figures/refdatum/forest_on_flat_{TILE}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")
