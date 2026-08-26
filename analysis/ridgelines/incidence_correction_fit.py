#!/usr/bin/env python3
"""Elevation change vs BEAM INCIDENCE, anchored on near-nadir beams declared correct.

The declaration: on low-curvature ridgetops, returns at incidence below --ref-max (default
0-5 deg) are taken to have measured the ground truly. Divides have zero contributing area,
so there is no overland-flow erosion and nothing upslope to deposit, and near-planar cells
give hillslope diffusion no purchase -- so on this ground the epochs SHOULD agree, and any
systematic departure at larger incidence is measurement geometry rather than change.

All canopy covers are included: the anchor is a beam-geometry band, not a land-cover class.
Because of that, the median canopy cover is printed per incidence bin -- if cover drifts
with incidence, part of the fitted curve is the canopy effect riding along, and the
per-cover-band fits at the end show whether the SHAPE is cover-independent.

Candidate forms are compared rather than assumed, since they encode different physics:
linear in angle, tan (a lateral/footprint displacement projected onto a slope), and
sec-1 = 1/cos(theta)-1 (a path-length or penetration term). Selected by AIC on n-weighted
binned medians with robust standard errors.

No incidence truncation. Every populated bin out to the largest observed incidence is
reported and fitted, with its own robust error bar; the sparse tail beyond ~35 deg is the
rare regime and carries the largest deltas, so deleting it changes the fitted shape (it
manufactures a turnover near 28 deg). The tail is also the highest-cover part of the range
-- corr(incidence, cover) is about +0.9 -- which is a reason to plot the cover alongside,
as the figure does, not a reason to remove the bins.

CAVEAT carried from the figures: the axis is incidence + delta with delta an unknown nadir
offset, so the angle origin is not guaranteed to be true nadir.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/incidence_correction_fit.py --tile data/derived/elbaext
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext")
ap.add_argument("--offset", default="corr", choices=("raw", "corr"))
ap.add_argument("--ref-max", type=float, default=5.0, help="incidence band declared correct (deg)")
ap.add_argument("--curv-max", type=float, default=0.015)
ap.add_argument("--no-ridge", action="store_true", help="drop the divide restriction")
ap.add_argument("--bin", type=float, default=2.0)
ap.add_argument("--min-n", type=int, default=1,
                help="returns needed to report a bin. 1 is the definitional floor -- a "
                     "median needs one point. Raising it DELETES the sparse high-incidence "
                     "bins, which carry the largest deltas; do not.")
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
TAG = ("" if TILE == "elba_fulldensity" else f"_{TILE}") + ("_noridge" if A.no_ridge else "")
DCOL = "d_mm_corr" if A.offset == "corr" else "d_mm"

df = pd.read_parquet(f"{A.tile}/beam_offset_table.parquet",
                     columns=["cell", DCOL, "incidence", "canopy_cover", "curv_laplacian", "in_grid"])
df = df[df.in_grid.values]
sel = ((df.curv_laplacian.abs().to_numpy() <= A.curv_max)
       & np.isfinite(df[DCOL].to_numpy()) & np.isfinite(df.incidence.to_numpy()))
if not A.no_ridge:
    sel &= np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()[df.cell.to_numpy()]
d = df[DCOL].to_numpy(float)[sel]; th = df.incidence.to_numpy(float)[sel]
cov = df.canopy_cover.to_numpy(float)[sel]
print("=" * 92)
print(f"ELEVATION CHANGE vs BEAM INCIDENCE  [{TILE}, {DCOL}]")
print(f"anchor: incidence < {A.ref_max:g} deg declared correct;  ground: "
      f"{'divides, ' if not A.no_ridge else ''}|Laplacian|<={A.curv_max:g};  ALL canopy covers")
print(f"{sel.sum():,} returns of {len(df):,} in-grid ({100*sel.sum()/len(df):.2f}%)")
print("=" * 92)

nmad = lambda a: 1.4826 * np.median(np.abs(a - np.median(a)))
ref_m = th < A.ref_max
REF = float(np.median(d[ref_m]))
REF_SE = 1.2533 * nmad(d[ref_m]) / np.sqrt(ref_m.sum())
print(f"\nANCHOR  incidence < {A.ref_max:g} deg: median d = {REF:+.2f} ± {REF_SE:.2f} mm "
      f"(n={ref_m.sum():,}, median cover {np.nanmedian(cov[ref_m]):.2f})")
print(f"        -> declared to be zero change; everything below is measured RELATIVE to it")

# Bin to the largest incidence actually observed. There is no upper truncation: the
# high-incidence tail is sparse but it holds the largest deltas, so it is shown with its
# (large) error bars rather than deleted.
edges = np.arange(0, np.ceil(th.max() / A.bin) * A.bin + A.bin, A.bin)
x, y, se, nn, cb = [], [], [], [], []
print(f"\n{'incidence':>12s} {'delta d (mm)':>13s} {'robust SE':>10s} {'med cover':>10s} {'n':>10s}")
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (th >= lo) & (th < hi)
    if m.sum() < A.min_n:
        print(f"{lo:5.0f}-{hi:<6.0f} {'--':>13s} {'':>10s} {'':>10s} {m.sum():>10,}"); continue
    v = d[m]
    x.append(th[m].mean()); y.append(np.median(v) - REF)
    # Robust SE of the bin median. A bin of one point (or of identical points) has a
    # degenerate NMAD of 0, which would give it infinite weight; fall back to the pooled
    # spread in that case so a sparse bin is uncertain rather than authoritative.
    s_bin = 1.2533 * nmad(v) / np.sqrt(m.sum())
    if not np.isfinite(s_bin) or s_bin <= 0:
        s_bin = 1.2533 * nmad(d) / np.sqrt(m.sum())
    se.append(np.hypot(s_bin, REF_SE))
    nn.append(int(m.sum())); cb.append(float(np.nanmedian(cov[m])))
    print(f"{lo:5.0f}-{hi:<6.0f} {y[-1]:>13.1f} {se[-1]:>10.2f} {cb[-1]:>10.2f} {m.sum():>10,}")
x, y, se, nn, cb = map(np.array, (x, y, se, nn, cb))
w = 1.0 / se**2
print(f"\n  CONFOUND CHECK  median cover across these bins: {cb.min():.2f} to {cb.max():.2f} "
      f"(drift {cb.max()-cb.min():+.2f}); corr(incidence, cover) = "
      f"{np.corrcoef(x, cb)[0,1]:+.2f}")

MODELS = {
    "linear      dz = b*theta":            lambda t: np.c_[t],
    "tan         dz = b*tan(theta)":       lambda t: np.c_[np.tan(np.radians(t))],
    "sec-1       dz = b*(1/cos(theta)-1)": lambda t: np.c_[1/np.cos(np.radians(t)) - 1],
    "linear+int  dz = a + b*theta":        lambda t: np.c_[np.ones_like(t), t],
    "quadratic   dz = b*th + c*th^2":      lambda t: np.c_[t, t**2],
}
print(f"\n{'model':>38s} {'chi2_red':>9s} {'AIC':>9s}   coefficients")
best = None
for name, design in MODELS.items():
    X = design(x); rw = np.sqrt(w)[:, None]
    b, *_ = np.linalg.lstsq(X*rw, y*np.sqrt(w), rcond=None)
    chi2 = float(np.sum(w*(y - X@b)**2)); k = X.shape[1]
    print(f"{name:>38s} {chi2/max(len(x)-k,1):>9.2f} {chi2+2*k:>9.1f}   "
          + ", ".join(f"{c:+.5g}" for c in b))
    if best is None or chi2+2*k < best[1]: best = (name, chi2+2*k, b, design)
name, aic, b, design = best
print(f"\nSELECTED (lowest AIC): {name}")
print(f"   dz(theta) = " + " + ".join(f"{c:+.5g}*term{i}" for i, c in enumerate(b)))
print(f"\n   apparent elevation change vs near-nadir, and the correction to apply:")
print(f"   {'theta':>7s} {'dz (mm)':>10s} {'correction (mm)':>17s}")
for t in (5, 10, 15, 20, 25, 30, 35):
    v = float((design(np.array([float(t)])) @ b)[0])
    print(f"   {t:>5d}   {v:>10.1f} {-v:>17.1f}")

# --- is the SHAPE cover-independent, or is the curve partly the canopy effect? ---
print(f"\n  SHAPE BY COVER BAND (each anchored on its OWN incidence<{A.ref_max:g} median):")
print(f"   {'cover band':>14s} {'n':>10s}" + "".join(f"{t:>8d}°" for t in (10, 15, 20, 25, 30)))
print(f"   {'':>14s} {'':>10s}   (n of each anchor and cell printed beneath the deltas)")
for lo, hi in ((0, .05), (.05, .20), (.20, .35), (.35, 1.01)):
    bm = np.isfinite(cov) & (cov >= lo) & (cov < hi)
    r = bm & (th < A.ref_max)
    if r.sum() == 0:          # definitional: the band needs an anchor to be measured from
        print(f"   {100*lo:5.0f}-{100*hi:<5.0f}% {bm.sum():>10,}   "
              f"(no return below {A.ref_max:g} deg to anchor on)"); continue
    r0 = float(np.median(d[r]))
    row = f"   {100*lo:5.0f}-{100*min(hi,1.0):<5.0f}% {bm.sum():>10,}"
    cnt = f"   {'anchor n=' + format(int(r.sum()), ',') :>14s} {'':>10s}"
    for t in (10, 15, 20, 25, 30):
        mm = bm & (th >= t-2) & (th < t+2)
        # every cell with at least one return is reported, with its n underneath;
        # a sparse cell is read against its count, not deleted for being sparse
        row += f"{np.median(d[mm])-r0:>+9.1f}" if mm.sum() else f"{'--':>9s}"
        cnt += f"{mm.sum():>9,}"
    print(row); print(cnt)

# Two panels on the SAME points: the left spans every populated bin so the sparse
# high-incidence tail is visible with its (large) error bars; the right zooms on the
# densely-sampled region so the model comparison stays readable. Nothing is cut from
# either -- the zoom is a viewport, not a filter, and its x-limit is stated on the axis.
ZOOM_MAX = 36.0
fig, (axf, ax) = plt.subplots(1, 2, figsize=(13.2, 5.6), dpi=130)
tt_full = np.linspace(0, x.max(), 400)
tt = np.linspace(0, ZOOM_MAX, 200)
fits = {}
for nm, dsg in MODELS.items():
    X = dsg(x); rw = np.sqrt(w)[:, None]
    fits[nm], *_ = np.linalg.lstsq(X*rw, y*np.sqrt(w), rcond=None)

for a, tgrid, ttl in ((axf, tt_full, f"ALL {len(x)} populated bins, to {x.max():.0f}°"),
                      (ax, tt, f"detail: 0–{ZOOM_MAX:g}° (same points, no bins removed)")):
    a.errorbar(x, y, yerr=se, fmt="o", ms=5, capsize=3, color="C0",
               label=f"binned median − anchor (n={sel.sum():,})", zorder=3)
    for nm, dsg in MODELS.items():
        a.plot(tgrid, dsg(tgrid) @ fits[nm], "-" if nm == name else "--",
               lw=2.0 if nm == name else 0.9, alpha=1.0 if nm == name else .5,
               label=nm.split()[0] + (" (selected)" if nm == name else ""))
    a.axhline(0, color="k", lw=.7)
    a.axvspan(0, A.ref_max, color="0.85", zorder=0)
    a.set_xlabel("incidence + δ (deg; δ = unknown nadir offset)")
    a.set_title(ttl, fontsize=9)
    a.grid(alpha=.3)
axf.set_ylabel("apparent elevation change vs near-nadir (mm)")
axf.legend(fontsize=8)
ax.set_xlim(-1, ZOOM_MAX)
_in = x <= ZOOM_MAX
ax.set_ylim(min(0, (y[_in]-se[_in]).min()) - 5, (y[_in]+se[_in]).max() + 5)
# the confound the tail rides on: median canopy cover per bin, same x
axc = ax.twinx()
axc.plot(x, cb, ":", color="C3", lw=1.4, label="median canopy cover")
axc.set_ylabel("median canopy cover in bin", color="C3", fontsize=9)
axc.tick_params(axis="y", labelcolor="C3", labelsize=8)
axc.legend(fontsize=8, loc="lower right")
fig.suptitle(f"elevation change vs beam incidence — {TILE}   "
             f"[anchor incidence < {A.ref_max:g}°, "
             f"{'divides, ' if not A.no_ridge else ''}|Laplacian|≤{A.curv_max:g}, all covers; "
             f"corr(incidence, cover) = {np.corrcoef(x, cb)[0,1]:+.2f}]", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = f"figures/refdatum/incidence_correction{TAG}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")
