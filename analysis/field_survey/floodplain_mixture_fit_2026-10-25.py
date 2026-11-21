#!/usr/bin/env python
"""Fit the stacked return-height distribution (returns minus the official 1 m DEM) as a mixture and ask where the
soil component sits, without the RTK (Andy, 2026-11-21: "do we model it as a sum of two Gaussians?").
Models, fitted by Poisson maximum likelihood on 1 cm bins over -0.45..+0.12 m (the upper range stops before the tall
herb and canopy; MINE):
  M1  one Gaussian (sward top)
  M2  two Gaussians (sward top + soil)
  M3  two Gaussians + a uniform bridge between the soil and sward means (merged echoes; Straatsma's layer)
  M4  M3 + an exponential tail below the soil mean (hollows and sub-cell relief)
Each on the sward region and the three distance-from-crown bands. Reported: parameters, log-likelihood, BIC, and the
soil mean against the RTK value (DEM minus soil: region 0.130, crown edge 0.120, near 0.145, open 0.121 m).
Outputs  data/field/2026-10-25/intensity/mixture_fit_log.txt; figures/rtk_2026-10-25_floodplain_mixture_fit.png
"""
import os
import numpy as np, pandas as pd, laspy, rasterio
from scipy.ndimage import distance_transform_edt
from scipy.optimize import minimize
from scipy.stats import norm
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/mixture_fit_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579780.0, 580160.0, 4885280.0, 4885760.0; RES = 1.0; nx, ny = 380, 480
with rasterio.open(f"{OUT}/intensity_2021_1m_tree.tif") as r: TREE = r.read(1) == 1
with rasterio.open(f"{OUT}/region_sward_mask.tif") as r: REG = r.read(1) == 1
with rasterio.open("data/field/2026-10-25/dem3dep/usgs_1m_dem_field.tif") as r: DEM = r.read(1).astype(float); DEM[DEM == r.nodata] = np.nan
z10 = np.load("data/field/2026-10-25/imagery/floodplain_model_features.npz", allow_pickle=True); xs, ys, fpm = z10["xs"], z10["ys"], z10["fpm"].astype(bool)
cy = N1 - (np.arange(ny) + 0.5) * RES; cx = E0 + (np.arange(nx) + 0.5) * RES
FP = fpm[np.clip(np.round((ys[0] - cy) / 10).astype(int), 0, len(ys) - 1)[:, None], np.clip(np.round((cx - xs[0]) / 10).astype(int), 0, len(xs) - 1)[None, :]]
dist = distance_transform_edt(~TREE) * RES; BAND = np.full((ny, nx), -1); BAND[(dist > 0) & (dist <= 3)] = 0; BAND[(dist > 3) & (dist <= 10)] = 1; BAND[dist > 10] = 2; BAND[~FP | TREE] = -1
parts = []
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0) & (x < E1) & (y >= N0) & (y < N1)
            if t == "4358-29-01": m &= x > 580030.0
            if not m.any(): continue
            i = np.clip(((N1 - y[m]) / RES).astype(int), 0, ny - 1); j = np.clip(((x[m] - E0) / RES).astype(int), 0, nx - 1); b = BAND[i, j]; k = (b >= 0) | REG[i, j]
            parts.append(pd.DataFrame(dict(dz=np.asarray(ch.z)[m][k] - DEM[i[k], j[k]], band=b[k], reg=REG[i, j][k])))
R = pd.concat(parts, ignore_index=True); R = R[np.isfinite(R.dz)]
sets = {"sward region": (R[R.reg].dz, 0.130), "crown edge 0-3 m": (R[R.band == 0].dz, 0.120), "near trees 3-10 m": (R[R.band == 1].dz, 0.145), "open > 10 m": (R[R.band == 2].dz, 0.121)}
LO, HI, W = -0.45, 0.12, 0.01; edges = np.arange(LO, HI + W / 2, W); mids = 0.5 * (edges[1:] + edges[:-1])
def comp_pdf(theta, model, x):
    # theta: mu1, s1, [w2, mu2, s2, [w3, [w4, lam]]]; weights are fractions of the fitted range's total; returns density on x
    mu1, s1 = theta[0], abs(theta[1]); pdf = norm.pdf(x, mu1, s1); parts = [pdf]
    if model >= 2:
        w2, mu2, s2 = theta[2], theta[3], abs(theta[4]); g2 = norm.pdf(x, mu2, s2); parts = [(1 - w2) * pdf, w2 * g2]
    if model >= 3:
        w3 = theta[5]; lo_, hi_ = min(mu1, mu2), max(mu1, mu2); u = np.where((x >= lo_) & (x <= hi_), 1.0 / max(hi_ - lo_, 1e-3), 0.0); parts = [p * (1 - w3) for p in parts] + [w3 * u]
    if model >= 4:
        w4, lam = theta[6], abs(theta[7]); e = np.where(x <= mu2, np.exp((x - mu2) / max(lam, 1e-3)) / max(lam, 1e-3), 0.0); parts = [p * (1 - w4) for p in parts] + [w4 * e]
    return parts
def nll(theta, model, cnt, n):
    parts = comp_pdf(theta, model, mids); dens = np.sum(parts, axis=0); mu = n * W * dens + 1e-12
    pen = 0.0
    if model >= 2 and not (0 < theta[2] < 0.5): pen += 1e6
    if model >= 3 and not (0 <= theta[5] < 0.5): pen += 1e6
    if model >= 4 and not (0 <= theta[6] < 0.5): pen += 1e6
    return float(np.sum(mu - cnt * np.log(mu))) + pen
results = {}
for name, (v, rtk) in sets.items():
    v = v[(v >= LO) & (v < HI)].to_numpy(); cnt = np.histogram(v, bins=edges)[0]; n = len(v); results[name] = (cnt, n, {})
    say(f"\n{name}: {n:,} returns in {LO}..{HI} m")
    inits = {1: [0.0, 0.02], 2: [0.0, 0.02, 0.05, -0.13, 0.04], 3: [0.0, 0.02, 0.05, -0.13, 0.04, 0.05], 4: [0.0, 0.02, 0.04, -0.13, 0.03, 0.04, 0.02, 0.06]}
    for model in (1, 2, 3, 4):
        best = None
        for mu2_0 in (-0.08, -0.13, -0.18):
            th0 = list(inits[model]); 
            if model >= 2: th0[3] = mu2_0
            r_ = minimize(nll, th0, args=(model, cnt, n), method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-6, fatol=1e-6))
            if best is None or r_.fun < best.fun: best = r_
        th = best.x; k = len(th); ll = -best.fun; bic = k * np.log(len(cnt)) + 2 * best.fun; results[name][2][model] = (th, ll, bic)
        desc = f"M{model}: sward mu {th[0]:+.3f} sigma {abs(th[1]):.3f}"
        if model >= 2: desc += f" | soil w {th[2]:.3f} mu {th[3]:+.3f} sigma {abs(th[4]):.3f}  -> DEM-soil {-th[3]:.3f} vs RTK {rtk:.3f} (diff {-th[3] - rtk:+.3f})"
        if model >= 3: desc += f" | bridge w {th[5]:.3f}"
        if model >= 4: desc += f" | tail w {th[6]:.3f} e-fold {abs(th[7]):.3f}"
        say(f"   {desc} | BIC {bic:.0f}")
fig, axs = plt.subplots(2, 2, figsize=(14, 9))
for ax, (name, (cnt, n, fits)) in zip(axs.ravel(), results.items()):
    ax.bar(mids, cnt / (n * W), width=W, color="0.8", label="data"); xx = np.linspace(LO, HI, 600)
    for model, c, ls in ((2, "#0072b2", "--"), (4, "#d55e00", "-")):
        th = fits[model][0]; parts = comp_pdf(th, model, xx); ax.plot(xx, np.sum(parts, axis=0), color=c, ls=ls, lw=1.5, label=f"M{model}, soil at {th[3]:+.3f} m")
        if model == 4:
            for p_, lab in zip(parts, ("sward top", "soil", "bridge", "tail")): ax.plot(xx, p_, color="#d55e00", lw=0.7, alpha=0.6)
    ax.axvline(-sets[name][1], color="k", lw=0.8, ls=":", label=f"RTK soil {-sets[name][1]:+.3f}"); ax.set_yscale("log"); ax.set_ylim(5e-3, 40); ax.set_xlim(LO, HI); ax.set_title(name, fontsize=10); ax.set_xlabel("return minus official DEM (m)"); ax.set_ylabel("density (1/m), log"); ax.legend(frameon=False, fontsize=8); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_mixture_fit.png", dpi=int(os.environ.get("FIGDPI", 150))); say("\nwrote figures/rtk_2026-10-25_floodplain_mixture_fit.png (thin lines: the M4 components; dotted: RTK soil, for the score only)")
LOG.close()
