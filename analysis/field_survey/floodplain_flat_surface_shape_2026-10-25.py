#!/usr/bin/env python
"""What shape does this sensor give on a hard flat surface? (Andy, 2026-11-21: research the shape returns should take.)
The road through the field site, in the same 2021 flight: road cells = 1 m cells with mean intensity below 2200, not
tree, more than 10 m from any tree cell, and lying on the road corridor (a hand-drawn polyline buffer of 4 m; cuts
MINE). Returns in those cells are detrended by a least-squares plane per 5 x 5 m window (windows with >= 40 returns and
plane RMS below 0.05 m kept; MINE) and the residuals pooled. Reported: histogram shape (log slopes), sigma, NMAD,
kurtosis, Laplace vs Gaussian likelihood, and the residual against intensity (range walk would show as a trend).
The sward's own core is treated the same way for comparison: region cells, same plane detrending.
Outputs  data/field/2026-10-25/intensity/flat_surface_shape_log.txt; figures/rtk_2026-10-25_floodplain_flat_surface_shape.png
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from scipy import stats
from scipy.ndimage import distance_transform_edt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/flat_surface_shape_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579780.0, 580160.0, 4885280.0, 4885760.0; RES = 1.0; nx, ny = 380, 480
with rasterio.open(f"{OUT}/intensity_2021_1m_mean.tif") as r: MI = r.read(1).astype(float); MI[MI == r.nodata] = np.nan
with rasterio.open(f"{OUT}/intensity_2021_1m_tree.tif") as r: TREE = r.read(1) == 1
with rasterio.open(f"{OUT}/region_sward_mask.tif") as r: REG = r.read(1) == 1
dist = distance_transform_edt(~TREE) * RES
# road corridor: polyline read off the intensity map (E, N), buffered 4 m
road = np.array([(579893, 4885760), (579900, 4885700), (579912, 4885640), (579925, 4885580), (579937, 4885520), (579947, 4885460), (579955, 4885400), (579962, 4885340), (579968, 4885280)], float)
cy = N1 - (np.arange(ny) + 0.5) * RES; cx = E0 + (np.arange(nx) + 0.5) * RES; CX, CY = np.meshgrid(cx, cy)
dmin = np.full(CX.shape, np.inf)
for (x0, y0), (x1, y1) in zip(road[:-1], road[1:]):
    t = np.clip(((CX - x0) * (x1 - x0) + (CY - y0) * (y1 - y0)) / ((x1 - x0) ** 2 + (y1 - y0) ** 2), 0, 1); dmin = np.minimum(dmin, np.hypot(CX - (x0 + t * (x1 - x0)), CY - (y0 + t * (y1 - y0))))
ROAD = (dmin <= 4) & (MI < 2200) & ~TREE & (dist > 10)
say(f"road cells: {int(ROAD.sum())}; intensity median {np.nanmedian(MI[ROAD]):.0f}")
parts = []
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0) & (x < E1) & (y >= N0) & (y < N1)
            if t == "4358-29-01": m &= x > 580030.0
            if not m.any(): continue
            i = np.clip(((N1 - y[m]) / RES).astype(int), 0, ny - 1); j = np.clip(((x[m] - E0) / RES).astype(int), 0, nx - 1); k = ROAD[i, j] | REG[i, j]
            parts.append(pd.DataFrame(dict(x=x[m][k], y=y[m][k], z=np.asarray(ch.z)[m][k], inten=np.asarray(ch.intensity)[m][k].astype(float), nr=np.asarray(ch.number_of_returns)[m][k], road=ROAD[i, j][k])))
R = pd.concat(parts, ignore_index=True)
def detrend(df, win=5.0):
    df = df.copy(); df["wx"] = np.floor((df.x - E0) / win).astype(int); df["wy"] = np.floor((N1 - df.y) / win).astype(int); out = []
    for _, g in df.groupby(["wx", "wy"]):
        if len(g) < 40: continue
        A = np.c_[g.x - g.x.mean(), g.y - g.y.mean(), np.ones(len(g))]; c, *_ = np.linalg.lstsq(A, g.z, rcond=None); res = g.z - A @ c
        if np.sqrt(np.mean(res ** 2)) < 0.05: out.append(pd.DataFrame(dict(res=res, inten=g.inten, nr=g.nr)))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["res", "inten", "nr"])
sets = {"road (hard, flat)": detrend(R[R.road]), "sward region": detrend(R[~R.road])}
W = 0.005; edges = np.arange(-0.15, 0.15 + W / 2, W); mids = 0.5 * (edges[1:] + edges[:-1])
fig, axs = plt.subplots(1, 3, figsize=(18, 5))
for c, (name, d) in zip(("k", "#0072b2"), sets.items()):
    v = d.res.to_numpy(); sd = v.std(); nmad = 1.4826 * np.median(np.abs(v - np.median(v))); kurt = stats.kurtosis(v)
    llg = stats.norm.logpdf(v, v.mean(), sd).sum(); b = np.mean(np.abs(v - np.median(v))); lll = stats.laplace.logpdf(v, np.median(v), b).sum()
    say(f"\n{name}: {len(v):,} detrended returns; sd {sd:.4f} m, NMAD {nmad:.4f}, excess kurtosis {kurt:.2f}; single-return share {100 * (d.nr == 1).mean():.1f} %; log-likelihood Laplace minus Gaussian per return {(lll - llg) / len(v):+.3f} (positive favors exponential flanks)")
    dens = np.histogram(v, bins=edges)[0] / (len(v) * W); ld = np.log10(np.where(dens > 0, dens, np.nan)); sl = np.gradient(ld, W)
    say("   log10-density slope (decades/m) at |res| = 1/2/3/4/6/8 cm, lower side: " + " ".join(f"{np.interp(-h, mids, sl):6.1f}" for h in (0.01, 0.02, 0.03, 0.04, 0.06, 0.08)) + " | upper side: " + " ".join(f"{np.interp(h, mids, sl):6.1f}" for h in (0.01, 0.02, 0.03, 0.04, 0.06, 0.08)))
    say("   (a Gaussian's slope grows linearly with |res|: at sigma, 2 sigma, 3 sigma it is 0.43/sigma, 0.87/sigma, 1.30/sigma; an exponential's is constant)")
    # residual vs intensity (range walk)
    q = d[(d.nr == 1)]; bins = np.percentile(q.inten, np.linspace(0, 100, 11)); lab = np.digitize(q.inten, bins[1:-1]); med = q.groupby(lab).agg(inten=("inten", "median"), res=("res", "median"), n=("res", "size"))
    say("   single returns, residual median by intensity decile: " + "; ".join(f"{int(r.inten)}: {1000 * r.res:+.0f} mm" for _, r in med.iterrows()) + f" | Spearman(res, intensity) {q[['res', 'inten']].corr(method='spearman').iloc[0, 1]:+.3f}")
    axs[0].plot(mids, dens, color=c, lw=1.5, label=f"{name} (n {len(v):,}, sd {100 * sd:.1f} cm)"); xx = np.linspace(-0.15, 0.15, 400); axs[0].plot(xx, stats.norm.pdf(xx, 0, sd), color=c, lw=0.8, ls="--"); axs[0].plot(xx, stats.laplace.pdf(xx, 0, b), color=c, lw=0.8, ls=":")
    axs[1].plot(mids, sl, color=c, lw=1.3, label=name); axs[2].plot(med.inten, 1000 * med.res, "o-", color=c, label=name)
axs[0].set_yscale("log"); axs[0].set_ylim(1e-2, 60); axs[0].set_xlabel("return minus local 5 m plane (m)"); axs[0].set_ylabel("density (1/m), log"); axs[0].legend(frameon=False, fontsize=8); axs[0].set_title("dashed: Gaussian of the same sd; dotted: Laplace of the same mean |res|", fontsize=9)
axs[1].set_xlim(-0.12, 0.12); axs[1].set_ylim(-60, 60); axs[1].set_xlabel("return minus local plane (m)"); axs[1].set_ylabel("slope of log10 density (decades/m)"); axs[1].legend(frameon=False, fontsize=8)
axs[2].set_xlabel("intensity (decile medians, single returns)"); axs[2].set_ylabel("median residual (mm)"); axs[2].legend(frameon=False, fontsize=8)
for a in axs: a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_flat_surface_shape.png", dpi=int(os.environ.get("FIGDPI", 150))); say("\nwrote figures/rtk_2026-10-25_floodplain_flat_surface_shape.png")

# scale-mixture test (theory note, 2026-11-21): is each intensity class Gaussian with its own width while the pool has exponential flanks?
d = sets["sward region"]; q = d[d.nr == 1]; bins = np.percentile(q.inten, np.linspace(0, 100, 11)); lab = np.digitize(q.inten, bins[1:-1])
say("\nsward region, single returns, by intensity decile: n, sd, NMAD, excess kurtosis, Laplace minus Gaussian log-likelihood per return (negative = Gaussian fits better):")
for k in range(10):
    v = q.res[lab == k].to_numpy(); v = v - np.median(v); sd = v.std(); b = np.mean(np.abs(v))
    say(f"   decile {k + 1:2d} intensity {bins[k]:5.0f}-{bins[k + 1]:5.0f}: n {len(v):6d}  sd {100 * sd:5.2f} cm  NMAD {100 * 1.4826 * np.median(np.abs(v)):5.2f}  kurtosis {stats.kurtosis(v):+5.2f}  Laplace-Gaussian {(stats.laplace.logpdf(v, 0, b).sum() - stats.norm.logpdf(v, 0, sd).sum()) / len(v):+.3f}")
v = q.res.to_numpy() - np.median(q.res); say(f"   all single returns pooled: n {len(v)} sd {100 * v.std():.2f} cm kurtosis {stats.kurtosis(v):+.2f} Laplace-Gaussian {(stats.laplace.logpdf(v, 0, np.mean(np.abs(v))).sum() - stats.norm.logpdf(v, 0, v.std()).sum()) / len(v):+.3f}")
LOG.close()
