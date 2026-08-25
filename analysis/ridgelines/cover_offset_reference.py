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
carry the fit, plus a straight line and a saturating form -- canopy cover is bounded and an
effect that grows without limit in cover is not physical.

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
ap.add_argument("--curv-max", type=float, default=None)
ap.add_argument("--min-n", type=int, default=200, help="minimum returns per cover bin")
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
TAG = ("" if TILE == "elba_fulldensity" else f"_{TILE}") + ("" if A.offset == "corr" else "_raw")
TAG += "_ridge" if A.ridge else ""
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
                     columns=["cell", DCOL, "slope", "canopy_cover", "curv_laplacian", "in_grid"])
df = df[df.in_grid.values].copy()
if A.curv_max is not None:
    df = df[(df.curv_laplacian.abs() <= A.curv_max).to_numpy()].copy()
cell = df.cell.to_numpy()
ref = (df.slope.to_numpy() < A.slope_max) & (tpi[cell] > 0)
label = f"low-gradient upland (slope<{A.slope_max:g} deg, TPI{A.tpi_window:g}>0)"
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
print("=" * 86)

EDGES = np.array([0, .02, .05, .10, .15, .20, .30, .40, .50, .65, 1.01])
nmad = lambda a: 1.4826 * np.median(np.abs(a - np.median(a)))
cc, mm, se, nn = [], [], [], []
print(f"{'cover bin':>14s} {'median d (mm)':>14s} {'robust SE':>10s} {'NMAD':>8s} {'n':>10s} {'med slope':>10s}")
for lo, hi in zip(EDGES[:-1], EDGES[1:]):
    m = (cov >= lo) & (cov < hi)
    if m.sum() < A.min_n:
        print(f"{lo:6.2f}-{hi:<6.2f} {'--':>14s} {'':>10s} {'':>8s} {m.sum():>10,}"); continue
    v = d[m]; s = nmad(v); e = 1.2533 * s / np.sqrt(m.sum())
    cc.append(cov[m].mean()); mm.append(np.median(v)); se.append(e); nn.append(int(m.sum()))
    print(f"{lo:6.2f}-{hi:<6.2f} {np.median(v):>14.1f} {e:>10.2f} {s:>8.1f} {m.sum():>10,} {np.median(sl[m]):>10.1f}")
cc, mm, se, nn = map(np.array, (cc, mm, se, nn))

if cc.size >= 3:
    w = 1.0 / se**2
    b = np.polyfit(cc, mm, 1, w=np.sqrt(w))
    pred = np.polyval(b, cc)
    chi2 = np.sum(w * (mm - pred)**2) / max(cc.size - 2, 1)
    print(f"\nLINEAR   d(mm) = {b[1]:+.1f} {b[0]:+.1f} * cover        "
          f"(SE-weighted; reduced chi2 {chi2:.1f})")
    print(f"   open (cover 0) -> {b[1]:+.1f} mm ;  full canopy (cover 1) -> {b[1]+b[0]:+.1f} mm")
    # saturating: cover cannot grow an effect without bound
    from scipy.optimize import curve_fit
    f = lambda x, a, c, k: a + c * (1 - np.exp(-x / k))
    try:
        p0 = [mm[0], mm[-1] - mm[0], 0.2]
        pop, _ = curve_fit(f, cc, mm, p0=p0, sigma=se, absolute_sigma=True, maxfev=20000)
        chi2s = np.sum(w * (mm - f(cc, *pop))**2) / max(cc.size - 3, 1)
        print(f"SATURATING d(mm) = {pop[0]:+.1f} {pop[1]:+.1f} * (1 - exp(-cover/{pop[2]:.3f}))  "
              f"(reduced chi2 {chi2s:.1f})")
    except Exception as e:                                    # noqa: BLE001
        pop = None; print(f"saturating fit did not converge: {e}")
else:
    b = None; pop = None; print("\ntoo few populated cover bins to fit")

fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4), dpi=130)
ax[0].errorbar(cc, mm, yerr=se, fmt="o-", ms=5, lw=1.4, capsize=3, color="C2",
               label="binned median ± robust SE")
if b is not None:
    xs = np.linspace(0, max(cc.max(), 0.7), 100)
    ax[0].plot(xs, np.polyval(b, xs), "--", color="0.35", lw=1.2,
               label=f"linear {b[1]:+.0f} {b[0]:+.0f}·cover")
    if pop is not None:
        ax[0].plot(xs, f(xs, *pop), "-", color="C1", lw=1.2, alpha=.8, label="saturating")
ax[0].axhline(0, color="k", lw=.7)
ax[0].set_xlabel("PyForestScan canopy cover (fraction)")
ax[0].set_ylabel(f"median offset {DCOL} (mm)")
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
