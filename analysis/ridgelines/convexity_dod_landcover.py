#!/usr/bin/env python3
"""STEPS 2-4 — transverse convexity on the ridgeline network, DoD along the crests, and
the DoD-vs-land-cover comparison.

Step 2 (convexity): at each divide (candidate ridge) cell, orient PERPENDICULAR to the
local ridge (cross-ridge = eigenvector of the smoothed-DEM Hessian with the most-negative
eigenvalue). Sample the ORIGINAL DEM along +/-L, fit z = a + b*s + c*s^2, convexity
kappa = -2c (>0 convex-up), cross-slope b (=0 at a true crest). Sweep L in {10,20,30} m.
A CREST cell = divide & convex (kappa>kmin) & near-crest (|b|<bmax) & on a topographic high.
Furrow-immune: tillage furrows are flat at +/-20 m -> kappa~0 -> rejected.

Step 3 (DoD): a convex crest sheds -> real change <= 0, so any POSITIVE DoD is an artifact.
Step 4 (land cover): forest (pen<0.25) vs open (pen>=0.45) crests, matched on slope.

Grid geometry comes from the tile's own corrections.json, so nothing here is tied to Elba.
`penetration.npy` is REQUIRED for step 4 unless you state `--without penetration`, in which
case steps 2-3 and the masks are still written and step 4 is skipped and said to be skipped.
Steps 2-3 need no cover layer at all, so `floodplain_mask.npy` and `crest_mask.npy` are
produced for any tile with z_after + slope + ridge_mask.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/convexity_dod_landcover.py \
        --tile elba_fulldensity --dod data/derived/elba_refdatum/dod_geoid.npy
"""
import argparse, json, os
import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter, distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lidar_diff_icp.viz import hillshade

_ap = argparse.ArgumentParser()
_ap.add_argument("--tile", default="elba_fulldensity",
                 help="tile name under data/derived/; grid geometry is read from its "
                      "corrections.json, so no origin is hardcoded")
_ap.add_argument("--dod", required=True,
                 help="path to the DoD grid to read (m, + = elevation rose). REQUIRED and "
                      "not defaulted: Elba's shipped run used data/derived/elba_refdatum/"
                      "dod_geoid.npy, which is NOT inside the tile directory, so no rule "
                      "of the form <tile>/dod_*.npy reproduces it. Name the grid you mean.")
_ap.add_argument("--without", default="",
                 help="comma-separated optional layers to run without, stated explicitly; "
                      "only 'penetration' is optional here (it drives step 4 alone)")
ARGS = _ap.parse_args()

TILE = ARGS.tile
D = f"data/derived/{TILE}"
OPTIONAL_STRATA = ("penetration",)
WITHOUT = {t.strip() for t in ARGS.without.split(",") if t.strip()}
_bad = WITHOUT - set(OPTIONAL_STRATA)
if _bad:
    raise SystemExit(f"--without names layers that are not optional here: {sorted(_bad)}; "
                     f"choose from {list(OPTIONAL_STRATA)}")
if WITHOUT:
    print(f"  running WITHOUT, as stated: {sorted(WITHOUT)}", flush=True)


def _grid(tile):                                    # (X0,Y0,RES) from the tile's own meta
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"data/derived/{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]
            return b[0], b[1], float(j.get("res") or j.get("res_m"))
    raise SystemExit(f"no grid meta for {tile}: none of meta.json, corrections_geoid.json, "
                     f"corrections.json exists under data/derived/{tile}")


def _opt(name):
    """REQUIRED unless named in --without. Missing and unstated -> refuse, because a layer
    that is absent must not read as one that was measured and came back empty."""
    if name in WITHOUT:
        return None
    p = f"{D}/{name}.npy"
    if not os.path.exists(p):
        raise SystemExit(f"{p} is missing. It is required for step 4 (the forest/open crest "
                         f"split). Produce it, or state that you are running without it: "
                         f"--without {name}")
    return np.load(p)


X0, Y0, RES = _grid(TILE)
z = np.load(f"{D}/z_after.npy")
dod = np.load(ARGS.dod)
slope = np.load(f"{D}/slope.npy")
pen = _opt("penetration")
ridge = np.load(f"{D}/ridge_mask.npy")
if dod.shape != z.shape:
    raise SystemExit(f"--dod {ARGS.dod} is {dod.shape} but {TILE}'s grid is {z.shape}; "
                     f"these are different grids and must not be combined")
ny, nx = z.shape
zf = z.copy(); nm = ~np.isfinite(zf)
if nm.any():
    zf = zf[tuple(distance_transform_edt(nm, return_distances=False, return_indices=True))]
# TPI is kept ONLY as a recorded per-cell covariate below (rec["tpi"]), never as a cut.
# It defines no population here any more: the floodplain is by elevation and the crest gate
# is convexity + a valley top. Andy, 2026-09-04: TPI is not used for floodplain anywhere.
tpi_large = zf - uniform_filter(zf, size=61)

# --- cross-ridge direction from the Hessian of a lightly smoothed DEM --------------------
zs = gaussian_filter(zf, sigma=3.0)                 # ~15 m, stable orientation
zx = np.gradient(zs, RES, axis=1); zy = np.gradient(zs, RES, axis=0)
zxx = np.gradient(zx, RES, axis=1); zyy = np.gradient(zy, RES, axis=0)
zxy = np.gradient(zx, RES, axis=0)
# principal-axis angle of the larger eigenvalue; cross-ridge (most negative) is +90 deg
phi = 0.5 * np.arctan2(2 * zxy, (zxx - zyy))
cross = phi + np.pi / 2
dx = np.cos(cross); dy = np.sin(cross)              # unit cross-ridge dir (east, north)

rr, cc = np.where(ridge)
Ec = X0 + (cc + 0.5) * RES; Nc = Y0 + (rr + 0.5) * RES
ddx = dx[rr, cc]; ddy = dy[rr, cc]

def bilinear(E, N):
    col = (E - X0) / RES - 0.5; row = (N - Y0) / RES - 0.5
    c0 = np.clip(np.floor(col).astype(int), 0, nx - 2); r0 = np.clip(np.floor(row).astype(int), 0, ny - 2)
    fc = np.clip(col - c0, 0, 1); fr = np.clip(row - r0, 0, 1)
    return ((1-fr)*(1-fc)*zf[r0, c0] + (1-fr)*fc*zf[r0, c0+1]
            + fr*(1-fc)*zf[r0+1, c0] + fr*fc*zf[r0+1, c0+1])

def convexity(L):
    s = np.arange(-L, L + 1e-6, RES)                # profile stations
    G = np.vstack([np.ones_like(s), s, s**2]).T
    pinv = np.linalg.pinv(G)                          # (3, len(s))
    Z = np.stack([bilinear(Ec + si*ddx, Nc + si*ddy) for si in s], axis=1)   # (ncells, ns)
    coef = (pinv @ Z.T)                               # (3, ncells): a, b, c
    return -2 * coef[2], coef[1]                      # kappa, cross-slope b

print("kappa (convexity, 1/m) distribution on divide cells, by profile half-length L:")
kap = {}
for L in [10, 20, 30]:
    k, b = convexity(L); kap[L] = (k, b)
    print(f"  L=+/-{L:>2} m: median kappa={np.median(k):+.4f}  "
          f"frac convex(kappa>0)={np.mean(k>0):.0%}  p75={np.percentile(k,75):+.4f}")

# map convexity kappa (each L) and cross-slope b back to the grid (for per-pixel records)
for L in (10, 20, 30):
    G = np.full((ny, nx), np.nan); G[rr, cc] = kap[L][0]
    np.save(f"{D}/kappa_L{L}.npy", G)
L0 = 20; kappa, bslope = kap[L0]
kappa_g = np.full((ny, nx), np.nan); b_g = np.full((ny, nx), np.nan)
kappa_g[rr, cc] = kappa; b_g[rr, cc] = bslope

# THE FLOODPLAIN MASK, BY ELEVATION (Andy, 2026-09-04). It was `tpi_large < -2.0`, and
# this file is the PRODUCER of floodplain_mask.npy, so that one line set the cut for every
# consumer in the project. Measured at whitewater, the TPI version missed 150,873 cells of
# valley floor and removed 16,218 upland hollows at a median 305.5 m -- wrong in both
# directions. The cut is now the first local minimum above the tile's dominant elevation
# mode; see refcells.valley_top_from_histogram for what that assumes and where it fails.
from lidar_diff_icp.refcells import VALLEY_TOP_M, valley_top_from_histogram
_vt = VALLEY_TOP_M.get(TILE) or valley_top_from_histogram(zf)
if _vt is None:
    floodplain = np.zeros_like(tpi_large, bool)
    print(f"NO valley cut for {TILE}: its elevation histogram has no minimum above the "
          f"dominant mode, so floodplain_mask.npy is EMPTY and every consumer of it "
          f"applies no floodplain cut on this tile.")
else:
    floodplain = np.isfinite(zf) & (zf < float(_vt))
    print(f"floodplain by ELEVATION < {float(_vt):.1f} m "
          f"({'registry' if VALLEY_TOP_M.get(TILE) else 'histogram'}): "
          f"{int(floodplain.sum()):,} cells ({100*floodplain.mean():.1f}% of the grid)")
np.save(f"{D}/floodplain_mask.npy", floodplain)

# crest = convex, near-crest, on a topographic high, NOT floodplain
KMIN = 0.004                                          # ~ >=0.8 m of convex relief over +/-20 m
crest_sel = (kappa > KMIN) & (np.abs(bslope) < 0.15) & (~floodplain[rr, cc])
crest = np.zeros((ny, nx), bool); crest[rr[crest_sel], cc[crest_sel]] = True
np.save(f"{D}/crest_mask.npy", crest)
print(f"\ncrest cells (convex, near-crest, above the valley top): {int(crest.sum())} "
      f"of {len(rr)} divide cells ({int(floodplain[rr, cc].sum())} divide cells were below it)")

# per-pixel record at every ridgecrest cell: slope AND curvature (+ DoD, land cover, coords)
cr, cco = np.where(crest)
rec = dict(row=cr, col=cco, E=X0+(cco+0.5)*RES, N=Y0+(cr+0.5)*RES,
           slope_deg=slope[cr, cco], curvature_kappa=kappa_g[cr, cco], cross_slope=b_g[cr, cco],
           dod_m=dod[cr, cco])
# column ORDER is preserved from the single-tile version; the two cover columns are simply
# absent when penetration is, rather than present and filled with a stand-in value
if pen is not None:
    rec["penetration"] = pen[cr, cco]
rec["tpi"] = tpi_large[cr, cco]
if pen is not None:
    rec["landcover"] = np.where(pen[cr, cco] < 0.25, "forest",
                                np.where(pen[cr, cco] >= 0.45, "open", "mixed"))
np.savez(f"{D}/ridgecrest_pixels.npz", **rec)
import csv
with open(f"{D}/ridgecrest_pixels.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(list(rec))
    for i in range(len(cr)):
        w.writerow([rec[k][i] for k in rec])
print(f"saved per-pixel slope+curvature table: ridgecrest_pixels.npz/.csv ({len(cr)} pixels, "
      f"cols: {', '.join(rec)})")

# ---- STEP 3: DoD along the crests ------------------------------------------------------
fin = np.isfinite(dod)
cm = crest & fin
print(f"\n=== STEP 3: DoD (gen2-gen1, mm) on convex crests (should be <=0; + = artifact) ===")
d = dod[cm]*1000
print(f"  all crests: n={cm.sum()}  median={np.median(d):+.1f}  "
      f"frac>0={np.mean(d>0):.0%}  NMAD={1.4826*np.median(np.abs(d-np.median(d))):.0f}")
print("  by slope:")
for lo,hi in [(0,3),(3,6),(6,10),(10,15),(15,90)]:
    m = cm & (slope>=lo)&(slope<hi)
    if m.any(): print(f"    {lo:>2}-{hi:<2} deg: n={m.sum():>4}  medDoD={np.median(dod[m]*1000):+6.1f} mm")

# ---- STEP 4: DoD on crests vs land cover ------------------------------------------------
if pen is None:
    print(f"\n=== STEP 4 SKIPPED: no penetration layer, stated via --without penetration ===")
    print("  the forest/open crest split is NOT reported for this tile; it is absent, not null")
else:
    print(f"\n=== STEP 4: convex-crest DoD by land cover (forest vs open), matched on slope ===")
    print(f"{'slope':>8} | {'forest n':>8} {'forest mm':>9} | {'open n':>7} {'open mm':>8}")
    forest = pen < 0.25; openc = pen >= 0.45
    for lo,hi in [(0,3),(3,6),(6,10),(10,15)]:
        s = cm & (slope>=lo)&(slope<hi)
        fn = s & forest; on = s & openc
        fm = np.median(dod[fn]*1000) if fn.any() else np.nan
        om = np.median(dod[on]*1000) if on.any() else np.nan
        print(f"{lo:>3}-{hi:<3} | {int(fn.sum()):>8} {fm:>+9.1f} | {int(on.sum()):>7} {om:>+8.1f}")
    mf = cm & forest; mo = cm & openc
    print(f"  ALL crests: forest n={int(mf.sum())} med={np.median(dod[mf]*1000):+.1f} mm  |  "
          f"open n={int(mo.sum())} med={np.median(dod[mo]*1000):+.1f} mm")

# ---- figure: crests colored by DoD -----------------------------------------------------
hs = hillshade(zf, RES, X0, Y0, fill_gaps=True)
ext = (X0, X0+nx*RES, Y0, Y0+ny*RES)
fig, ax = plt.subplots(figsize=(11, 13))
ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
cd = np.where(crest & fin, dod, np.nan)
im = ax.imshow(cd, extent=ext, origin="lower", cmap="RdBu", vmin=-0.08, vmax=0.08)
fig.colorbar(im, ax=ax, shrink=0.6, label="DoD gen2-gen1 (m) on convex crests")
ax.set_title(f"Steps 2-4: convex ridge crests (n={int(crest.sum())}), DoD.\n"
             "convex crest sheds -> blue (increase) = artifact")
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
_fig = ("figures/ridgeline_crests_dod.png" if TILE == "elba_fulldensity"
        else f"figures/ridgeline_crests_dod_{TILE}.png")
fig.savefig(_fig, dpi=130, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {_fig}")
