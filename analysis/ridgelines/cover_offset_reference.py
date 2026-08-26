#!/usr/bin/env python3
"""Calibrate VERTICAL OFFSET against FOREST DENSITY on ground that should not be eroding.

The problem with fitting an offset-vs-cover relationship anywhere else is that cover and
erosion are both organised by topography, so a cover term fitted across a whole tile
absorbs real geomorphic change. This picks a reference population where the change term is
as close to zero as the landscape allows and only the measurement effect should remain:

    LOW-GRADIENT UPLAND      slope < --slope-max AND large-scale TPI > 0
                             gentle enough that overland-flow erosion is negligible, and
                             ABOVE the local mean so it is not receiving deposition either.
    (optionally) DIVIDES     --ridge adds the Scherler & Schwanghart divide network, which
                             has zero contributing area by construction.

On that population the offset is regressed on PyForestScan canopy cover. Reported as binned
medians with a robust standard error (1.2533*NMAD/sqrt(n)) so the reader can see which bins
carry the fit, and candidate FORMS ARE COMPARED, not assumed: linear, quadratic, and a Beer-Lambert optical
depth -ln(1-cover), which is the physically motivated one if cover behaves like 1 minus a gap
fraction. Selection is by AIC on the weighted binned medians. The binned medians themselves
are the model-free answer; the fitted forms are descriptions of them.

Uses the REGISTRATION-CORRECTED offset by default: on raw d_mm the per-swath misalignment
would be read as a cover effect wherever cover and flight-line geometry covary.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/cover_offset_reference.py --tile data/derived/elbaext
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy.ndimage import uniform_filter, distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext")
ap.add_argument("--offset", default="corr", choices=("raw", "corr"))
ap.add_argument("--slope-max", type=float, default=12.0)
ap.add_argument("--tpi-window", type=float, default=500.0, help="upland test window (m)")
ap.add_argument("--ridge", action="store_true", help="ALSO include divide cells of any slope")
ap.add_argument("--ref", default="upland", choices=("upland", "divides"),
                help="upland = low-gradient, TPI>0 (default); divides = the S&S divide "
                     "network ONLY, which has zero contributing area by construction")
ap.add_argument("--inc-max", type=float, default=None,
                help="keep only returns below this beam incidence (deg). Fixing incidence "
                     "removes the beam-geometry term, leaving cover as the only variable -- "
                     "but note incidence ~ slope for near-nadir beams, so this also selects "
                     "near-flat ground and cannot speak to the cover effect on steep slopes")
ap.add_argument("--curv-max", type=float, default=None)
ap.add_argument("--min-n", type=int, default=1,
                help="minimum returns per cover bin. 1 is the definitional floor (a median "
                     "needs a return). It removes nothing at 200 on these tiles either, but "
                     "the top cover bin is only n=273 (elba) / 400 (elbaext) and carries the "
                     "-140 / -130 mm point that is the whole finding, so any threshold here "
                     "is one tile away from deleting the result")
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
TAG = ("" if TILE == "elba_fulldensity" else f"_{TILE}") + ("" if A.offset == "corr" else "_raw")
TAG += "_ridge" if A.ridge else ""
TAG += "" if A.ref == "upland" else "_divides"
TAG += "" if A.inc_max is None else f"_inc{A.inc_max:g}"
DCOL = "d_mm_corr" if A.offset == "corr" else "d_mm"


def grid_meta(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            return b[0], b[1], r
    raise SystemExit(f"no grid meta in {tile}")


X0, Y0, RES = grid_meta(A.tile)
z = np.load(f"{A.tile}/z_after.npy")
zf = z.copy(); miss = ~np.isfinite(zf)
if miss.any():
    zf = zf[tuple(distance_transform_edt(miss, return_distances=False, return_indices=True))]
k = max(3, int(round(A.tpi_window / RES)) | 1)
tpi = (zf - uniform_filter(zf, size=k, mode="nearest")).ravel()

df = pd.read_parquet(f"{A.tile}/beam_offset_table.parquet",
                     columns=["cell", DCOL, "slope", "canopy_cover", "curv_laplacian",
                              "incidence", "in_grid"])
df = df[df.in_grid.values].copy()
if A.curv_max is not None:
    df = df[(df.curv_laplacian.abs() <= A.curv_max).to_numpy()].copy()
cell = df.cell.to_numpy()
if A.ref == "divides":
    ref = np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()[cell]
    label = "divides (zero contributing area)"
else:
    ref = (df.slope.to_numpy() < A.slope_max) & (tpi[cell] > 0)
    label = f"low-gradient upland (slope<{A.slope_max:g} deg, TPI{A.tpi_window:g}>0)"
if A.inc_max is not None:
    ref &= df.incidence.to_numpy() < A.inc_max
    label += f", incidence<{A.inc_max:g} deg"
if A.ridge:
    rm = np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()
    ref = ref | rm[cell]
    label += " + divides"
ref &= np.isfinite(df[DCOL].to_numpy()) & np.isfinite(df.canopy_cover.to_numpy())
d = df[DCOL].to_numpy(float)[ref]; cov = df.canopy_cover.to_numpy(float)[ref]
sl = df.slope.to_numpy(float)[ref]

print("=" * 86)
print(f"OFFSET vs FOREST DENSITY on non-eroding reference ground  [{TILE}]")
print(f"reference = {label};  offset = {DCOL}")
print(f"{ref.sum():,} returns of {len(df):,} in-grid  ({100*ref.sum()/len(df):.1f}%)")
print(f"median slope of the reference population: {np.median(df.slope.to_numpy()[ref]):.1f} deg;"
      f"  median incidence {np.nanmedian(df.incidence.to_numpy()[ref]):.1f} deg")
print("=" * 86)

EDGES = np.array([0, .02, .05, .10, .15, .20, .30, .40, .50, .65, 1.01])
nmad = lambda a: 1.4826 * np.median(np.abs(a - np.median(a)))
cc, mm, se, nn = [], [], [], []
print(f"{'cover bin':>14s} {'median d (mm)':>14s} {'robust SE':>10s} {'NMAD':>8s} {'n':>10s} {'med slope':>10s}")
for lo, hi in zip(EDGES[:-1], EDGES[1:]):
    m = (cov >= lo) & (cov < hi)
    if m.sum() < max(1, A.min_n):
        print(f"{lo:6.2f}-{hi:<6.2f} {'--':>14s} {'':>10s} {'':>8s} {m.sum():>10,}"); continue
    v = d[m]; s = nmad(v); e = 1.2533 * s / np.sqrt(m.sum())
    cc.append(cov[m].mean()); mm.append(np.median(v)); se.append(e); nn.append(int(m.sum()))
    print(f"{lo:6.2f}-{hi:<6.2f} {np.median(v):>14.1f} {e:>10.2f} {s:>8.1f} {m.sum():>10,} {np.median(sl[m]):>10.1f}")
cc, mm, se, nn = map(np.array, (cc, mm, se, nn))

if cc.size >= 3:
    w = 1.0 / se**2
    # --- compare candidate forms rather than assume linearity ---
    # optical depth: if canopy cover behaves like 1 - gap fraction, Beer-Lambert makes the
    # path length through canopy proportional to -ln(1 - cover), which ACCELERATES as cover
    # approaches 1 -- the opposite of saturating, and what the high-cover bins actually do.
    FORMS = {
        "linear         d = a + b*cover":        lambda c: np.c_[np.ones_like(c), c],
        "quadratic      d = a + b*c + e*c^2":    lambda c: np.c_[np.ones_like(c), c, c**2],
        "optical depth  d = a + b*(-ln(1-c))":   lambda c: np.c_[np.ones_like(c),
                                                                -np.log(1 - np.clip(c, 0, 0.98))],
    }
    print(f"\n{'form':>36s} {'chi2_red':>9s} {'AIC':>9s}   coefficients")
    bestf = None
    for nm, dsg in FORMS.items():
        X = dsg(cc); rw = np.sqrt(w)[:, None]
        bb, *_ = np.linalg.lstsq(X*rw, mm*np.sqrt(w), rcond=None)
        chi2 = float(np.sum(w*(mm - X@bb)**2)); k = X.shape[1]
        print(f"{nm:>36s} {chi2/max(cc.size-k,1):>9.1f} {chi2+2*k:>9.1f}   "
              + ", ".join(f"{v:+.4g}" for v in bb))
        if bestf is None or chi2+2*k < bestf[1]: bestf = (nm, chi2+2*k, bb, dsg)
    print(f"   SELECTED: {bestf[0]}   coefficients "
          + ", ".join(f"{v:+.5g}" for v in bestf[2]))
    print(f"   predicted offset at cover 0.1/0.3/0.5/0.7: "
          + ", ".join(f"{float((bestf[3](np.array([c])) @ bestf[2])[0]):+.0f}" for c in (.1,.3,.5,.7)) + " mm")
    b, pop, f = bestf[2], None, bestf[3]
else:
    b = None; pop = None; f = None; bestf = None
    print("\ntoo few populated cover bins to fit")

fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4), dpi=130)
ax[0].errorbar(cc, mm, yerr=se, fmt="o-", ms=5, lw=1.4, capsize=3, color="C2",
               label="binned median ± robust SE")
if bestf is not None:
    xs = np.linspace(0, min(max(cc.max(), 0.7), 0.95), 200)
    for nm, dsg in FORMS.items():
        X = dsg(cc); rw = np.sqrt(w)[:, None]
        bb, *_ = np.linalg.lstsq(X*rw, mm*np.sqrt(w), rcond=None)
        sel = nm == bestf[0]
        ax[0].plot(xs, dsg(xs) @ bb, "-" if sel else "--", lw=1.8 if sel else 0.9,
                   alpha=1.0 if sel else .55,
                   label=nm.split()[0] + (" (selected)" if sel else ""))
ax[0].axhline(0, color="k", lw=.7)
ax[0].set_xlabel("PyForestScan canopy cover (fraction)")
ax[0].set_ylabel(f"median offset {DCOL} (mm)   [gen1 \u2212 gen2; + = lower in 2021]")
ax[0].set_title(f"offset vs forest density on non-eroding ground\n{label}", fontsize=10)
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].bar(cc, nn, width=0.03, color="0.6")
ax[1].set_yscale("log"); ax[1].set_xlabel("canopy cover (fraction)")
ax[1].set_ylabel("returns per bin"); ax[1].set_title("sample size behind each point", fontsize=10)
ax[1].grid(alpha=.3)
fig.suptitle(f"forest density vs vertical offset — {TILE}, {DCOL}", y=1.02)
out = f"figures/refdatum/cover_offset_reference{TAG}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")
