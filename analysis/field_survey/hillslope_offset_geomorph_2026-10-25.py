#!/usr/bin/env python
"""Do the bluff offsets (2021 lidar floor above the 2026 surveyed ground) follow geomorphic drivers? (Andy, 2026-11-21:
soil creep via slope and curvature, overland-flow erosion via topographic convergence; the alternative is a duff layer.)
Points: the 39 hillslope FIX points minus RS3 10, RS3 24 and T4a 23 (Andy's call: one displaced 1.4 m, one degraded fix,
one footing), 36 points. Offset = official 1 m DEM minus RTK, slope-normal (hillslope_points_summary.csv).
Covariates from a wider official 1 m DEM window (E 579500-580050, N 4885200-4885900, fetched here into
dem3dep/usgs_1m_dem_bluff_geomorph.tif so that flow routing is not cut at the bluff top; the window is MINE):
  slope (deg) and aspect; Laplacian curvature (positive = concave, hollow) of the DEM smoothed with Gaussians of 2.5 and
  5 m sigma (scales MINE), plan and profile curvature at the 5 m scale; D8 upslope contributing area (cells drained in
  descending elevation order to the steepest lower neighbor, pits left as sinks) and the topographic wetness index
  ln(area / tan slope); and, from the point stacks, canopy share and ground share.
Spearman rank correlations of the offset with each covariate, with a permutation p-value; an OLS on standardized
slope, Laplacian (5 m) and log area with bootstrap intervals; the same for the first traverse alone. n = 36, so the
test can see only strong structure.
Outputs  data/field/2026-10-25/intensity/hillslope_offset_geomorph.csv, _log.txt; figures/rtk_2026-10-25_hillslope_offset_geomorph.png
"""
import os
import numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, TwoSlopeNorm
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; LOG = open(f"{OUT}/hillslope_offset_geomorph_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = 579500.0, 580050.0, 4885200.0, 4885900.0; DEMF = "data/field/2026-10-25/dem3dep/usgs_1m_dem_bluff_geomorph.tif"
if not os.path.exists(DEMF):
    arrs = []
    for u in ("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1m/Projects/MN_SE_Driftless_2021_B21/TIFF/USGS_1M_15_x57y489_MN_SE_Driftless_2021_B21.tif", "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1m/Projects/MN_SE_Driftless_2021_B21/TIFF/USGS_1M_15_x58y489_MN_SE_Driftless_2021_B21.tif"):
        r = rasterio.open("/vsicurl/" + u); w = from_bounds(max(E0, r.bounds.left), max(N0, r.bounds.bottom), min(E1, r.bounds.right), min(N1, r.bounds.top), transform=r.transform); a = r.read(1, window=w); t = r.window_transform(w); arrs.append((a, t, r.nodata, r.crs))
    full = np.full((int(N1 - N0), int(E1 - E0)), np.nan, np.float32)
    for a, t, nod, crs in arrs:
        r0 = int(round(N1 - t.f)); c0 = int(round(t.c - E0)); a = a.astype(np.float32); a[a == nod] = np.nan; full[r0:r0 + a.shape[0], c0:c0 + a.shape[1]] = np.where(np.isnan(full[r0:r0 + a.shape[0], c0:c0 + a.shape[1]]), a, full[r0:r0 + a.shape[0], c0:c0 + a.shape[1]])
    with rasterio.open(DEMF, "w", driver="GTiff", height=full.shape[0], width=full.shape[1], count=1, dtype="float32", crs=arrs[0][3], transform=rasterio.transform.from_origin(E0, N1, 1.0, 1.0), nodata=-9999.0, compress="deflate") as dst: dst.write(np.where(np.isfinite(full), full, -9999.0), 1)
with rasterio.open(DEMF) as r: Z = r.read(1).astype(float); Z[Z == r.nodata] = np.nan
Zf = np.where(np.isfinite(Z), Z, np.nanmedian(Z))
def grads(a): gy, gx = np.gradient(a, 1.0); return gx, -gy
gx, gy = grads(gaussian_filter(Zf, 1.0)); SLOPE = np.degrees(np.arctan(np.hypot(gx, gy))); ASPECT = (np.degrees(np.arctan2(-gx, -gy))) % 360
def curvs(sig):
    a = gaussian_filter(Zf, sig); zx, zy = grads(a); zxx = np.gradient(zx, 1.0, axis=1); zyy = -np.gradient(zy, 1.0, axis=0); zxy = -np.gradient(zx, 1.0, axis=0)
    p2 = zx ** 2 + zy ** 2 + 1e-9; lap = zxx + zyy
    prof = -(zxx * zx ** 2 + 2 * zxy * zx * zy + zyy * zy ** 2) / (p2 * (1 + p2) ** 1.5); plan = -(zxx * zy ** 2 - 2 * zxy * zx * zy + zyy * zx ** 2) / (p2 ** 1.5)
    return lap, prof, plan
LAP25, _, _ = curvs(2.5); LAP5, PROF5, PLAN5 = curvs(5.0)
# D8 accumulation
ny, nx = Zf.shape; order = np.argsort(-Zf.ravel()); acc = np.ones(ny * nx); nbr = [(-1, -1, np.sqrt(2)), (-1, 0, 1), (-1, 1, np.sqrt(2)), (0, -1, 1), (0, 1, 1), (1, -1, np.sqrt(2)), (1, 0, 1), (1, 1, np.sqrt(2))]
zr = Zf.ravel()
for idx in order:
    i, j = divmod(idx, nx); best = None; bs = 0
    for di, dj, dd in nbr:
        ii, jj = i + di, j + dj
        if 0 <= ii < ny and 0 <= jj < nx:
            s_ = (zr[idx] - zr[ii * nx + jj]) / dd
            if s_ > bs: bs = s_; best = ii * nx + jj
    if best is not None: acc[best] += acc[idx]
ACC = acc.reshape(ny, nx); TWI = np.log(ACC / np.maximum(np.tan(np.radians(SLOPE)), 1e-3))
def at(a, x, y): return map_coordinates(a, [((N1 - y) - 0.5), ((x - E0) - 0.5)], order=1, mode="nearest")
S = pd.read_csv(f"{OUT}/hillslope_points_summary.csv"); excl = {("hillslope_RS3", 10), ("hillslope_RS3", 24), ("hillslope_T4a", 23)}
S = S[~S.apply(lambda r: (r.group, int(r.point)) in excl, axis=1)].copy(); S["off"] = S.dem_minus_rtk_n
for nm, a in (("slope_dem", SLOPE), ("aspect", ASPECT), ("lap_2.5m", LAP25), ("lap_5m", LAP5), ("prof_5m", PROF5), ("plan_5m", PLAN5), ("acc_m2", ACC), ("twi", TWI)): S[nm] = at(a, S.E.values, S.N.values)
S["log_acc"] = np.log10(S.acc_m2); S.to_csv(f"{OUT}/hillslope_offset_geomorph.csv", index=False, float_format="%.4f")
say(f"{len(S)} points after excluding RS3 10, RS3 24, T4a 23; offset (DEM minus RTK, slope-normal) median {S.off.median():+.3f}, IQR {S.off.quantile(.25):+.3f}..{S.off.quantile(.75):+.3f}")
say("covariates at the points, median (p10-p90): " + "; ".join(f"{c} {S[c].median():.3g} ({S[c].quantile(.1):.3g}-{S[c].quantile(.9):.3g})" for c in ("slope_dem", "lap_2.5m", "lap_5m", "prof_5m", "plan_5m", "acc_m2", "twi")))
def perm_p(x, y, n=20000, seed=0):
    r = np.random.default_rng(seed); obs = stats.spearmanr(x, y).statistic; return obs, float(np.mean([abs(stats.spearmanr(x, r.permutation(y)).statistic) >= abs(obs) for _ in range(n)]))
say("\nSpearman rank correlation of the offset with each covariate (permutation p, 20,000 draws), all 36 | first traverse (16) | second (20):")
for c in ("slope_dem", "lap_2.5m", "lap_5m", "prof_5m", "plan_5m", "log_acc", "twi", "canopy_share", "ground_share", "inten_near"):
    line = f"   {c:13s}"
    for lab, q in (("all", S), ("first", S[S.group == "hillslope_RS3"]), ("second", S[S.group == "hillslope_T4a"])):
        rho, p = perm_p(q[c].values, q.off.values); line += f"  {lab} {rho:+.2f} (p {p:.2f})"
    say(line)
# OLS on standardized covariates
X = S[["slope_dem", "lap_5m", "log_acc"]]; Xs = (X - X.mean()) / X.std(); A = np.c_[np.ones(len(S)), Xs.values]; y = S.off.values
b = np.linalg.lstsq(A, y, rcond=None)[0]; res = y - A @ b; r2 = 1 - res.var() / y.var()
rng = np.random.default_rng(0); bb = np.array([np.linalg.lstsq(A[k], y[k], rcond=None)[0] for k in (rng.integers(0, len(y), len(y)) for _ in range(5000))])
say(f"\nOLS offset ~ slope + Laplacian(5 m) + log area (standardized), n {len(S)}: R^2 {r2:.2f}; coefficients (m per sd) with bootstrap 95 %: " + ", ".join(f"{nm} {b[k + 1]:+.3f} ({np.percentile(bb[:, k + 1], 2.5):+.3f}..{np.percentile(bb[:, k + 1], 97.5):+.3f})" for k, nm in enumerate(["slope", "Laplacian", "log area"])))
say("sign conventions: Laplacian > 0 concave (hollow, creep deposition expected), < 0 convex (nose, creep erosion); larger log area = more convergent, overland-flow erosion expected to lower the ground, i.e. a POSITIVE offset if erosion since 2021")
say(f"expected magnitudes for reference (hypotheses, not measured here): soil creep surface velocities of mm/yr give vertical change over 5.5 years far below 1 cm except at sharp curvature; overland-flow rilling can remove cm in a few events in convergent hollows")
fig, axs = plt.subplots(1, 4, figsize=(21, 5), gridspec_kw=dict(width_ratios=[1, 1, 1, 1.5]))
for ax, c, lab in zip(axs[:3], ("slope_dem", "lap_5m", "log_acc"), ("slope (deg, 1 m DEM)", "Laplacian curvature, 5 m scale (1/m; > 0 concave)", "log10 upslope area (m^2)")):
    for g, col in (("hillslope_RS3", "#0072b2"), ("hillslope_T4a", "#d55e00")):
        q = S[S.group == g]; ax.plot(q[c], q.off, "o", color=col, ms=6, label=g.replace("hillslope_", "RS3 "))
    ax.axhline(0, color="0.4", lw=0.8, ls=":"); ax.set_xlabel(lab); ax.set_ylabel("DEM minus RTK, slope-normal (m)"); ax.spines[["top", "right"]].set_visible(False)
axs[0].legend(frameon=False, fontsize=8)
ax = axs[3]; ls = LightSource(315, 45); pad = 25; e0, e1, n0, n1 = S.E.min() - pad, S.E.max() + pad, S.N.min() - pad, S.N.max() + pad
ax.imshow(ls.hillshade(Zf, vert_exag=1, dx=1, dy=1), cmap="gray", extent=(E0, E1, N0, N1)); ax.imshow(np.clip(LAP5, -0.05, 0.05), cmap="RdBu", vmin=-0.05, vmax=0.05, extent=(E0, E1, N0, N1), alpha=0.45)
sc = ax.scatter(S.E, S.N, c=S.off, cmap="PuOr_r", norm=TwoSlopeNorm(0, -0.15, 0.25), s=45, edgecolors="k", linewidths=0.5); plt.colorbar(sc, ax=ax, shrink=0.7, label="DEM minus RTK (m)")
ax.set_xlim(e0, e1); ax.set_ylim(n0, n1); ax.set_aspect("equal"); ax.set_title("offsets on hillshade; blue = concave, red = convex (Laplacian, 5 m)", fontsize=9); ax.tick_params(labelsize=7)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_hillslope_offset_geomorph.png", dpi=int(os.environ.get("FIGDPI", 150))); say("wrote figures/rtk_2026-10-25_hillslope_offset_geomorph.png")
LOG.close()
