#!/usr/bin/env python3
"""Weighted inversion of the rigid shift (dE, dN, dZ) between the 25 Oct 2026 RTK survey and the 2021
lidar ground, every logged point included (FIX, FLOAT, SINGLE; Andy, 2026-11-11: "Instead of
classifying FIX points, perform the inversion taking all uncertainties into account").
  r_i = Z2021(E_i + dE, N_i + dN) - (navd88_tied_i + dZ)                 (lidar minus shifted RTK)
  sigma_i^2 = rms_v,i^2  +  (tan(slope_i) * rms_h,i)^2  +  s_lidar,i^2
     rms_v,i    the rover's reported vertical RMS for that point (FLOAT and SINGLE carry their own, large, values)
     rms_h,i    the rover's reported horizontal RMS (hypot of E and N; T4a's table reports the lateral RMS directly),
                turned into height through the local lidar slope, since a horizontal error on a slope is a height error
     s_lidar,i  the rms of the 2021 class-2 returns within 3 m of the point about their best-fit plane, from the
                tile's class-2 cloud: the lidar surface's own roughness at the scale the 5 m grid averages
  chi2(dE, dN, dZ) = sum_i (r_i / sigma_i)^2, full grid search; 1-sigma = where the profile chi2 rises by 1.
Choices that are mine: the 3 m disk, the grid (dE, dN +/-1.5 m by 0.05; dZ -0.2..+0.6 m by 0.01), and the
error model above. Nothing is cut; a point's weight is its 1/sigma^2. A point whose 3 m disk holds no
class-2 return has no measured lidar surface to compare against (the grid there is filled from
neighbours) and is listed and left out: that is a missing datum, not a cut.
TWO SOLUTIONS ARE PRINTED. (1) The sigmas as reported. On this survey T1a's FLOAT points sit 4-24 m
off the lidar with reported vertical RMS of 0.05-0.24 m, so the reported FLOAT RMS is not an uncertainty
and that solution is theirs. (2) Variance components by status: each status's sigmas are scaled by the
robust spread of its residuals over sigma (1.4826 * median|r/sigma|) from solution 1, so the data say
what each solution type's errors really are; the factors are printed. The 1-sigma ranges are where the
profile chi2 rises by 1 from the minimum.
Writes data/derived/4342-29-64/bluff_check/gen2_local_rms_all_points.csv (cached),
       data/derived/4342-29-64/bluff_check/weighted_shift.json,
       figures/rtk_2026-10-25_shift_weighted.png (chi2 slices through the minimum; star = minimum, circle = no shift).
"""
import json, os, sys
import laspy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

T = "4342-29-64"; D = f"data/derived/{T}/bluff_check"; os.makedirs(D, exist_ok=True)
b = json.load(open(f"data/derived/{T}/corrections.json"))["bounds"]; res = 5.0
Z = np.load(f"data/derived/{T}/z_after.npy")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")
RADIUS = 3.0
SUBSET = sys.argv[1] if len(sys.argv) > 1 else "all"   # "all" or "hillslope" (Andy, 2026-11-11: "the hillslope lidar alone in the inversion")
SUF = "" if SUBSET == "all" else f"_{SUBSET}"

# --- lidar roughness at every point, from the class-2 cloud (cached)
cache = f"{D}/gen2_local_rms_all_points.csv"
if os.path.exists(cache):
    R = pd.read_csv(cache)
else:
    x0, x1 = P.easting_tied.min() - 10, P.easting_tied.max() + 10; y0, y1 = P.northing_tied.min() - 10, P.northing_tied.max() + 10
    xs, ys, zs = [], [], []
    with laspy.open(f"data/after/valley/{T}_3dep_class2.laz") as F:
        for ch in F.chunk_iterator(2_000_000):
            x, y = ch.x, ch.y; m = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
            if m.any(): xs.append(x[m]); ys.append(y[m]); zs.append(ch.z[m])
    x, y, z = np.concatenate(xs), np.concatenate(ys), np.concatenate(zs); tree = cKDTree(np.c_[x, y]); rows = []
    for _, p in P.iterrows():
        idx = tree.query_ball_point([p.easting_tied, p.northing_tied], RADIUS); n = len(idx)
        if n >= 4:
            A = np.c_[x[idx] - p.easting_tied, y[idx] - p.northing_tied, np.ones(n)]; c, *_ = np.linalg.lstsq(A, z[idx], rcond=None)
            rms = float(np.sqrt(np.mean((z[idx] - A @ c) ** 2)))
        else: rms = np.nan
        rows.append(dict(group=p.group, point=p.point, gen2_n_3m=n, gen2_plane_rms_m=rms))
    R = pd.DataFrame(rows); R.to_csv(cache, index=False); print(f"class-2 returns kept in the survey box: {len(x)}; cached {cache}")
P = P.merge(R, on=["group", "point"], how="left")
if SUBSET == "hillslope": P = P[P.group.str.startswith("hillslope")].reset_index(drop=True)
print(f"SUBSET: {SUBSET}, {len(P)} logged points")

# --- the error model
rms_h = np.where(P.rms_n.isna(), P.rms_e, np.hypot(P.rms_e, P.rms_n))
sig_rtk_h = np.tan(np.radians(P.slope_deg.to_numpy())) * rms_h
sig = np.sqrt(P.rms_v.to_numpy() ** 2 + sig_rtk_h ** 2 + P.gen2_plane_rms_m.to_numpy() ** 2)
ok = np.isfinite(sig) & P.navd88_tied_m.notna().to_numpy()
P["sigma_m"] = sig; P["sigma_rtk_v_m"] = P.rms_v; P["sigma_rtk_h_as_height_m"] = sig_rtk_h
left = P[~ok]; print(f"points: {int(ok.sum())} of {len(P)} enter. Left out, no class-2 return within {RADIUS} m (the grid there is filled, not measured): " + "; ".join(f"{g} {len(d)} pts ({', '.join(str(int(q)) for q in d.point)})" for g, d in left.groupby("group")))
print("by status (sigmas as reported):")
print(P[ok].groupby("status").agg(n=("point", "size"), rms_v_med=("rms_v", "median"), h_as_height_med=("sigma_rtk_h_as_height_m", "median"), lidar_rms_med=("gen2_plane_rms_m", "median"), sigma_med=("sigma_m", "median"), sigma_max=("sigma_m", "max")).round(3).to_string())
E, N, h, w = P.easting_tied.to_numpy()[ok], P.northing_tied.to_numpy()[ok], P.navd88_tied_m.to_numpy()[ok], 1.0 / sig[ok] ** 2
print(f"share of the total weight: floodplain {w[P.group[ok].str.startswith('floodplain').to_numpy()].sum() / w.sum():.2f}, hillslope {w[P.group[ok].str.startswith('hillslope').to_numpy()].sum() / w.sum():.2f}, benchmark {w[(P.group[ok] == 'benchmark_DR9925').to_numpy()].sum() / w.sum():.2f}; FIX {w[P.status[ok].str.startswith('FIX').to_numpy()].sum() / w.sum():.2f}")

def z2021(x, y):
    return map_coordinates(Z, [(y - b[1]) / res - 0.5, (x - b[0]) / res - 0.5], order=1, mode="nearest")

gE = np.round(np.arange(-1.5, 1.5001, 0.05), 2); gZ = np.round(np.arange(-0.2, 0.6001, 0.01), 2)
status = P.status[ok].str.split(" ").str[0].to_numpy()   # FIX / FLOAT / SINGLE

def solve(w, label):
    C = np.empty((len(gE), len(gE), len(gZ)))
    for i, dE in enumerate(gE):
        for j, dN in enumerate(gE):
            r = (z2021(E + dE, N + dN) - h)[:, None] - gZ[None, :]; C[i, j] = np.sum(w[:, None] * r ** 2, axis=0)
    k = np.unravel_index(C.argmin(), C.shape); nu = len(w) - 3; chi2nu = C[k] / nu; z0 = (np.argmin(np.abs(gE)), np.argmin(np.abs(gE)), np.argmin(np.abs(gZ)))
    def one_sigma(axis, grid):
        others = tuple(a for a in range(3) if a != axis); p = C.min(axis=others); inside = grid[p <= C[k] + 1.0]; return float(inside.min()), float(inside.max())
    r_best = z2021(E + gE[k[0]], N + gE[k[1]]) - h - gZ[k[2]]
    print(f"\n{label}: dE {gE[k[0]]:+.2f} dN {gE[k[1]]:+.2f} dZ {gZ[k[2]]:+.2f} m   chi2 {C[k]:.1f} on nu = {nu} -> reduced chi2 {chi2nu:.2f}  (no shift: chi2 {C[z0]:.1f})")
    res = dict(label=label, dE=float(gE[k[0]]), dN=float(gE[k[1]]), dZ=float(gZ[k[2]]), chi2=float(C[k]), nu=nu, chi2_nu=float(chi2nu), n=len(w), one_sigma={})
    for a, g_, nm in ((0, gE, "dE"), (1, gE, "dN"), (2, gZ, "dZ")):
        lo, hi = one_sigma(a, g_); res["one_sigma"][nm] = (lo, hi); print(f"  {nm}: 1-sigma (delta chi2 = 1) {lo:+.2f} .. {hi:+.2f} m")
    sg = 1 / np.sqrt(w); print(f"  residuals at the minimum: unweighted NMAD {1.4826 * np.median(np.abs(r_best - np.median(r_best))) * 1000:.0f} mm; |r/sigma| > 3 at {int((np.abs(r_best) / sg > 3).sum())} of {len(w)} points; weight share floodplain {w[grp_f].sum() / w.sum():.2f} hillslope {w[grp_h].sum() / w.sum():.2f}")
    return res, C, k, r_best

grp_f = P.group[ok].str.startswith("floodplain").to_numpy(); grp_h = P.group[ok].str.startswith("hillslope").to_numpy()
w0 = 1.0 / sig[ok] ** 2
res1, C1, k1, r1 = solve(w0, "SOLUTION 1, sigmas as reported")
# variance components by status, from the residuals of solution 1 (one pass: reported sigma scaled by the robust spread of r/sigma within each status)
sg = sig[ok].copy(); factors = {}
for st in np.unique(status):
    m = status == st; f = 1.4826 * np.median(np.abs(r1[m] / sg[m])); factors[st] = f; sg[m] *= f
print("\nvariance-component factors (reported sigma x factor = sigma the residuals show): " + ", ".join(f"{st} x{f:.1f} (n {int((status == st).sum())})" for st, f in factors.items()))
res2, C2, k2, r2 = solve(1.0 / sg ** 2, "SOLUTION 2, sigmas calibrated by status")
P.loc[ok, "sigma_calibrated_m"] = sg; P.loc[ok, "r_at_solution2_m"] = r2
json.dump(dict(radius_m=RADIUS, error_model="sigma^2 = rms_v^2 + (tan(slope)*rms_h)^2 + gen2_plane_rms_3m^2; chi2 = sum (r/sigma)^2; r = lidar - (RTK + dZ) at the shifted position", left_out=left[["group", "point", "status"]].to_dict("records"),
               solution_1_reported_sigmas=res1, variance_component_factors=factors, solution_2_calibrated_sigmas=res2), open(f"{D}/weighted_shift{SUF}.json", "w"), indent=1)
P.to_csv(f"{D}/weighted_shift_points{SUF}.csv", index=False)
C, k = C2, k2

fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
D2 = C - C[k]; vmax = 50 if SUBSET == "all" else 20
panels = ((D2[:, k[1], :].T, gE, gZ, "dE (m)", "dZ (m)", (gE[k[0]], gZ[k[2]])), (D2[k[0], :, :].T, gE, gZ, "dN (m)", "dZ (m)", (gE[k[1]], gZ[k[2]])), (D2[:, :, k[2]].T, gE, gE, "dE (m)", "dN (m)", (gE[k[0]], gE[k[1]])))
for ax, (A, gx, gy, xl, yl, best) in zip(axs, panels):
    im = ax.imshow(A, origin="lower", extent=[gx[0], gx[-1], gy[0], gy[-1]], aspect="auto", cmap="viridis", vmin=0, vmax=vmax)
    ax.contour(gx, gy, A, levels=[1.0, 4.0], colors=["white", "0.7"], linewidths=0.8); ax.plot(*best, "w*", ms=12, mec="k"); ax.plot(0, 0, "wo", ms=6, mec="k"); ax.set_xlabel(xl); ax.set_ylabel(yl)
fig.colorbar(im, ax=axs, shrink=0.9, pad=0.02).set_label("$\\chi^2 - \\chi^2_{min}$ (calibrated sigmas)")
fig.savefig(f"figures/rtk_2026-10-25_shift_weighted{SUF}.png", dpi=120, bbox_inches="tight"); print(f"\nfigures/rtk_2026-10-25_shift_weighted{SUF}.png")
