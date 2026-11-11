#!/usr/bin/env python3
"""Nuth & Kaab (2011) test of HORIZONTAL misregistration of the 2008 (gen1) cloud against
the 2021 (gen2) ground grid on tile 4342-29-64, run on the RAW gen1 last returns.

Hypothesis under test: the ~0.4 m gen1-below-gen2 offset on the 30 deg east-facing bluff
(RTK traverse window rows 483:550, cols 440:469) is a horizontal shift of gen1 relative
to gen2 (aspect-dependent, dh/tan(slope) = a cos(aspect - b) + c), not a slope-dependent
incidence error (dh = f(slope), aspect-independent).

Conventions (also in NUTH_KAAB_RAW.md):
  dh      = z_gen1 - z_gen2(bilinear at the return's x, y)      [m]
  slope, aspect from np.gradient of the reference grid; aspect = DOWNSLOPE direction,
  degrees clockwise from north.
  A gen1 return truly at p but recorded at p + s carries z(p - s), so
      dh = -s . grad z = tan(slope) * (s_x sin(aspect) + s_y cos(aspect))
  => s_x = a sin(b), s_y = a cos(b) is the displacement of gen1 RELATIVE TO gen2, and
  the CORRECTION to add to gen1 coordinates is (-s_x, -s_y). Checked numerically below.

Writes every printed block to the report (append) and tables to bluff_check/.
Shared laptop: the cloud is loaded once; no multiprocessing.
"""
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import laspy  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402

ROOT = "/home/awickert/projects/lidar-diff-icp"
TILE = "4342-29-64"
DER = f"{ROOT}/data/derived/{TILE}"
OUT = f"{DER}/bluff_check"
REPORT = f"{ROOT}/analysis/field_survey/bluff/NUTH_KAAB_RAW.md"
RAW = f"{ROOT}/data/before/{TILE}.laz"
CSF = f"{ROOT}/data/csf_cache/{TILE}.las"

SLOPE_MIN_DEG = 15.0          # Andy's band edge (brief); his to change
NBINS = 36                    # brief
SLOPE_BANDS = [0, 5, 10, 15, 20, 25, 30, 35, 45]   # brief
WINDOW = (slice(483, 550), slice(440, 469))         # brief: RTK traverse
NBOOT = 1000
RNG = np.random.default_rng(20261111)
QUADS = {"E 45-135": (45, 135), "S 135-225": (135, 225),
         "W 225-315": (225, 315), "N 315-45": (315, 45)}

os.makedirs(OUT, exist_ok=True)
_rep = open(REPORT, "a")
_rep.write("\n## Script output (verbatim, appended as produced)\n\n```\n")
_rep.flush()


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    _rep.write(s + "\n")
    _rep.flush()


say(f"# nuth_kaab_raw_{TILE}.py  run {time.strftime('%Y-%m-%d %H:%M:%S')}")

# ----------------------------------------------------------------------------- inputs
corr = json.load(open(f"{DER}/corrections.json"))
X0, Y0, X1, Y1 = corr["bounds"]
RES = float(corr["res_m"])
STAR = {int(k): np.array(v, float) for k, v in corr["per_swath_star_frame_dxdydz_m"].items()}
ALIGN = {int(k): np.array(v, float)
         for k, v in corr["per_swath_internal_alignment_dxdydz_m"].items()}
LAT = np.array(corr["cross_epoch_datum"]["horizontal_shift_m"], float)
say(f"corrections.json: swath_frame={corr['swath_frame']} zero_line={corr['zero_line']}")
say(f"  per_swath_star_frame_dxdydz_m (APPLIED, added to x,y,z): "
    f"{ {k: v.tolist() for k, v in STAR.items()} }")
say(f"  per_swath_internal_alignment_dxdydz_m (SOLVED, NOT applied under star): "
    f"{ {k: v.tolist() for k, v in ALIGN.items()} }")
say(f"  cross_epoch_datum.horizontal_shift_m (APPLIED after star, added to x,y): {LAT.tolist()}")

FP = np.load(f"{DER}/floodplain_mask_clubb.npy").astype(bool)
GRIDS = {
    "Zref(z_after_differenced, slope-normal median)": np.load(f"{DER}/z_after_differenced.npy"),
    "Z21(z_after, horizontal median)": np.load(f"{DER}/z_after.npy"),
}
NY, NX = FP.shape
say(f"grid {NY}x{NX} res {RES} bounds {[X0, Y0, X1, Y1]}; floodplain cells {int(FP.sum())}")


class Ref:
    """Reference grid + gradients, with bilinear sampling at points."""

    def __init__(self, name, Z):
        self.name = name
        self.Z = Z
        gy, gx = np.gradient(Z, RES)            # axis 0 = north (rows), axis 1 = east
        self.gx, self.gy = gx, gy
        self.slope = np.degrees(np.arctan(np.hypot(gx, gy)))
        self.aspect = np.degrees(np.arctan2(-gx, -gy)) % 360.0

    def sample(self, x, y):
        r = (y - Y0) / RES - 0.5
        c = (x - X0) / RES - 0.5
        coords = np.vstack([r, c])
        z = map_coordinates(self.Z, coords, order=1, mode="constant", cval=np.nan)
        gx = map_coordinates(self.gx, coords, order=1, mode="constant", cval=np.nan)
        gy = map_coordinates(self.gy, coords, order=1, mode="constant", cval=np.nan)
        slope = np.degrees(np.arctan(np.hypot(gx, gy)))
        aspect = np.degrees(np.arctan2(-gx, -gy)) % 360.0
        return z, slope, aspect


def cell_index(x, y):
    ix = np.floor((x - X0) / RES).astype(int)
    iy = np.floor((y - Y0) / RES).astype(int)
    ok = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
    return iy, ix, ok


def nmad(v):
    v = v[np.isfinite(v)]
    return 1.4826 * np.median(np.abs(v - np.median(v)))


# ----------------------------------------------------------------------------- the clouds
t0 = time.time()
with laspy.open(RAW) as f:
    P = f.read()
x_raw = np.asarray(P.x, float); y_raw = np.asarray(P.y, float); z_raw = np.asarray(P.z, float)
psid_raw = np.asarray(P.point_source_id).astype(int)
last = np.asarray(P.return_number) == np.asarray(P.number_of_returns)
say(f"raw cloud {RAW}: {len(x_raw)} returns, {int(last.sum())} last returns, "
    f"read in {time.time()-t0:.1f} s")
x_raw, y_raw, z_raw, psid_raw = x_raw[last], y_raw[last], z_raw[last], psid_raw[last]
del P, last
say(f"  last returns per swath: { {int(k): int(v) for k, v in zip(*np.unique(psid_raw, return_counts=True))} }")

t0 = time.time()
with laspy.open(CSF) as f:
    P = f.read()
x_csf = np.asarray(P.x, float); y_csf = np.asarray(P.y, float); z_csf = np.asarray(P.z, float)
psid_csf = np.asarray(P.point_source_id).astype(int)
say(f"csf cache {CSF}: {len(x_csf)} points (all class 2), read in {time.time()-t0:.1f} s")
del P


def shifted(x, y, z, psid, star=False, align=False, lateral=False):
    """Apply the pipeline's shifts, ADDING the recorded values (chain.StarFrame,
    chain.LateralShift conventions)."""
    x = x.copy(); y = y.copy(); z = z.copy()
    for s in np.unique(psid):
        m = psid == s
        d = np.zeros(3)
        if align:
            d += ALIGN[int(s)]
        if star:
            d += STAR[int(s)]
        x[m] += d[0]; y[m] += d[1]; z[m] += d[2]
    if lateral:
        x += LAT[0]; y += LAT[1]
    return x, y, z


# ----------------------------------------------------------------------------- the fit
def nk_fit(dh, slope, aspect, label, save=None):
    """Nuth & Kaab on 36 aspect-bin medians of dh/tan(slope) for slope > SLOPE_MIN_DEG.
    Returns dict with a, b, c, s_x, s_y, their sigmas, NMAD before/after."""
    ok = np.isfinite(dh) & np.isfinite(slope) & (slope > SLOPE_MIN_DEG)
    y = dh[ok] / np.tan(np.radians(slope[ok]))
    asp = aspect[ok]
    edges = np.linspace(0, 360, NBINS + 1)
    ib = np.clip(np.digitize(asp, edges) - 1, 0, NBINS - 1)
    cnt = np.bincount(ib, minlength=NBINS)
    med = np.array([np.median(y[ib == k]) if cnt[k] else np.nan for k in range(NBINS)])
    ctr = np.radians(0.5 * (edges[:-1] + edges[1:]))
    use = np.isfinite(med)
    G = np.column_stack([np.cos(ctr[use]), np.sin(ctr[use]), np.ones(use.sum())])
    coef, *_ = np.linalg.lstsq(G, med[use], rcond=None)
    A, B, c = coef                       # A = a cos b = s_y ; B = a sin b = s_x
    resid_bins = med[use] - G @ coef
    dof = max(use.sum() - 3, 1)
    s2 = (resid_bins ** 2).sum() / dof
    cov = s2 * np.linalg.inv(G.T @ G)
    sA, sB, sc = np.sqrt(np.diag(cov))
    a = np.hypot(A, B); b = np.degrees(np.arctan2(B, A)) % 360.0
    # bootstrap over bins
    boots = []
    idx_all = np.nonzero(use)[0]
    for _ in range(NBOOT):
        pick = RNG.choice(idx_all, size=len(idx_all), replace=True)
        Gb = np.column_stack([np.cos(ctr[pick]), np.sin(ctr[pick]), np.ones(len(pick))])
        try:
            cb, *_ = np.linalg.lstsq(Gb, med[pick], rcond=None)
        except np.linalg.LinAlgError:
            continue
        boots.append(cb)
    boots = np.array(boots)
    bsx, bsy = boots[:, 1], boots[:, 0]
    ba = np.hypot(boots[:, 0], boots[:, 1])
    # NMAD of dh/tan(slope) before and after removing the fit (evaluated per point)
    fit_pt = A * np.cos(np.radians(asp)) + B * np.sin(np.radians(asp)) + c
    nm_before = nmad(y); nm_after = nmad(y - fit_pt)
    # NMAD of dh itself, before and after removing only the aspect term
    dh_ok = dh[ok]
    nm_dh_before = nmad(dh_ok)
    nm_dh_after = nmad(dh_ok - np.tan(np.radians(slope[ok])) * (A * np.cos(np.radians(asp)) + B * np.sin(np.radians(asp))))
    r = dict(case=label, n=int(ok.sum()), n_bins=int(use.sum()), min_bin_n=int(cnt[use].min()),
             a_m=a, b_deg=b, c_m=c, s_x=B, s_y=A,
             sig_sx_lsq=sB, sig_sy_lsq=sA, sig_c_lsq=sc,
             sig_sx_boot=bsx.std(ddof=1), sig_sy_boot=bsy.std(ddof=1), sig_a_boot=ba.std(ddof=1),
             a_boot_p2_5=np.percentile(ba, 2.5), a_boot_p97_5=np.percentile(ba, 97.5),
             corr_dx_to_gen1=-B, corr_dy_to_gen1=-A,
             nmad_y_before=nm_before, nmad_y_after=nm_after,
             nmad_dh_before=nm_dh_before, nmad_dh_after=nm_dh_after,
             median_dh=float(np.median(dh_ok)), median_slope=float(np.median(slope[ok])),
             bin_rms_resid=float(np.sqrt(s2)))
    say(f"[{label}] n={r['n']} bins={r['n_bins']} (min bin n={r['min_bin_n']}) "
        f"a={a:.4f} m b={b:.1f} deg c={c:.4f} m | s_x={B:+.4f}±{sB:.4f}(lsq)/{bsx.std(ddof=1):.4f}(boot) "
        f"s_y={A:+.4f}±{sA:.4f}/{bsy.std(ddof=1):.4f} | a boot 95% [{r['a_boot_p2_5']:.4f},{r['a_boot_p97_5']:.4f}] "
        f"| correction to add to gen1 (dx,dy)=({-B:+.4f},{-A:+.4f}) "
        f"| NMAD dh/tan: {nm_before:.4f} -> {nm_after:.4f} m | NMAD dh: {nm_dh_before:.4f} -> {nm_dh_after:.4f} m "
        f"| median dh {r['median_dh']:+.4f} m, median slope {r['median_slope']:.1f}, bin rms resid {r['bin_rms_resid']:.4f}")
    if save:
        pd.DataFrame(dict(aspect_bin_center_deg=np.degrees(ctr), n=cnt, median_dh_over_tan=med,
                          fit=A * np.cos(ctr) + B * np.sin(ctr) + c)).to_csv(save, index=False)
    return r


def band_table(dh, slope, aspect, label, fp_ok):
    """Median dh (m) in slope bands, pooled and per aspect quadrant. fp_ok = non-floodplain."""
    rows = []
    for lo, hi in zip(SLOPE_BANDS[:-1], SLOPE_BANDS[1:]):
        m = fp_ok & np.isfinite(dh) & (slope >= lo) & (slope < hi)
        row = dict(case=label, band=f"{lo}-{hi}", n=int(m.sum()),
                   pooled=float(np.median(dh[m])) if m.any() else np.nan)
        for q, (qa, qb) in QUADS.items():
            if qa < qb:
                mq = m & (aspect >= qa) & (aspect < qb)
            else:
                mq = m & ((aspect >= qa) | (aspect < qb))
            row[q] = float(np.median(dh[mq])) if mq.any() else np.nan
            row["n " + q] = int(mq.sum())
        rows.append(row)
    df = pd.DataFrame(rows)
    say(f"-- median dh (m) by gen2 slope band, {label}")
    say("   band      n_pooled   pooled   E(45-135)   S(135-225)   W(225-315)   N(315-45)     nE      nS      nW      nN")
    for _, r in df.iterrows():
        say(f"   {r['band']:<7} {r['n']:>9d}  {r['pooled']:+.4f}    {r['E 45-135']:+.4f}      {r['S 135-225']:+.4f}       "
            f"{r['W 225-315']:+.4f}      {r['N 315-45']:+.4f}  {r['n E 45-135']:>6d}  {r['n S 135-225']:>6d}  {r['n W 225-315']:>6d}  {r['n N 315-45']:>6d}")
    return df


def evaluate(ref, x, y, z):
    z2, sl, asp = ref.sample(x, y)
    iy, ix, ok = cell_index(x, y)
    fp = np.ones(len(x), bool)
    fp[ok] = FP[iy[ok], ix[ok]]
    good = ok & ~fp & np.isfinite(z2)
    return z - z2, sl, asp, good


# ============================================================================= step 2
say("\n=== STEP 2: Nuth & Kaab fits (slope > %.0f deg, non-floodplain) ===" % SLOPE_MIN_DEG)
REF = Ref(*list(GRIDS.items())[0])
REF21 = Ref(*list(GRIDS.items())[1])
say(f"primary reference: {REF.name}; sensitivity reference: {REF21.name}")
fits = []
band_dfs = []

# (i) raw pooled
dh, sl, asp, good = evaluate(REF, x_raw, y_raw, z_raw)
fits.append(nk_fit(dh[good], sl[good], asp[good], "i raw pooled", save=f"{OUT}/bins_i_raw_pooled.csv"))
# sign check: apply the correction and refit
sx, sy = fits[-1]["s_x"], fits[-1]["s_y"]
dh_c, sl_c, asp_c, good_c = evaluate(REF, x_raw - sx, y_raw - sy, z_raw)
chk = nk_fit(dh_c[good_c], sl_c[good_c], asp_c[good_c], "i raw pooled, after adding (-s_x,-s_y) to gen1 x,y [SIGN CHECK]")
say(f"  sign check: a went {fits[-1]['a_m']:.4f} -> {chk['a_m']:.4f} m; a residual a near zero confirms the sign convention")
# (ii) raw per swath
for s in sorted(STAR):
    m = good & (psid_raw == s)
    fits.append(nk_fit(dh[m], sl[m], asp[m], f"ii raw swath {s}", save=f"{OUT}/bins_ii_raw_swath{s}.csv"))

# (iii-a) registered AS THE PIPELINE APPLIED IT: star + lateral
xr, yr, zr = shifted(x_raw, y_raw, z_raw, psid_raw, star=True, lateral=True)
dhr, slr, aspr, goodr = evaluate(REF, xr, yr, zr)
fits.append(nk_fit(dhr[goodr], slr[goodr], aspr[goodr], "iii-a registered(star+lateral, AS APPLIED) pooled",
                   save=f"{OUT}/bins_iiia_registered_pooled.csv"))
fit_reg_pooled = fits[-1]
fit_reg_swath = {}
for s in sorted(STAR):
    m = goodr & (psid_raw == s)
    fit_reg_swath[s] = nk_fit(dhr[m], slr[m], aspr[m], f"iii-a registered(star+lateral) swath {s}",
                              save=f"{OUT}/bins_iiia_registered_swath{s}.csv")
    fits.append(fit_reg_swath[s])

# (iii-b) the brief's description: internal alignment + star (NOT what the pipeline applied)
xb, yb, zb = shifted(x_raw, y_raw, z_raw, psid_raw, align=True, star=True, lateral=False)
dhb, slb, aspb, goodb = evaluate(REF, xb, yb, zb)
fits.append(nk_fit(dhb[goodb], slb[goodb], aspb[goodb], "iii-b internal+star (HYPOTHETICAL, not applied) pooled"))
for s in sorted(STAR):
    m = goodb & (psid_raw == s)
    fits.append(nk_fit(dhb[m], slb[m], aspb[m], f"iii-b internal+star (hypothetical) swath {s}"))
del xb, yb, zb, dhb, slb, aspb, goodb

# (iv) CSF cache, registered (star + lateral), pooled
xc, yc, zc = shifted(x_csf, y_csf, z_csf, psid_csf, star=True, lateral=True)
dhc, slc, aspc, goodc = evaluate(REF, xc, yc, zc)
fits.append(nk_fit(dhc[goodc], slc[goodc], aspc[goodc], "iv CSF cache registered(star+lateral) pooled",
                   save=f"{OUT}/bins_iv_csf_registered_pooled.csv"))
for s in sorted(STAR):
    m = goodc & (psid_csf == s)
    fits.append(nk_fit(dhc[m], slc[m], aspc[m], f"iv CSF cache registered swath {s}"))

# sensitivity: z_after.npy (horizontal median) as reference
dh21, sl21, asp21, good21 = evaluate(REF21, x_raw, y_raw, z_raw)
fits.append(nk_fit(dh21[good21], sl21[good21], asp21[good21], "sens: i raw pooled vs z_after.npy (Z21)"))
dh21r, sl21r, asp21r, good21r = evaluate(REF21, xr, yr, zr)
fits.append(nk_fit(dh21r[good21r], sl21r[good21r], asp21r[good21r], "sens: iii-a registered pooled vs z_after.npy (Z21)"))
del dh21, sl21, asp21, good21, dh21r, sl21r, asp21r, good21r

FT = pd.DataFrame(fits)
FT.to_csv(f"{OUT}/nuth_kaab_fits.csv", index=False)
say("\n-- FIT TABLE (m, deg). s_x,s_y = displacement of gen1 relative to gen2; correction = -s")
say("   case | n | a | b | c | s_x(dx) | s_y(dy) | sig_sx(lsq/boot) | sig_sy(lsq/boot) | NMAD dh/tan before->after | NMAD dh before->after")
for _, r in FT.iterrows():
    say(f"   {r['case']} | {r['n']} | {r['a_m']:.4f} | {r['b_deg']:.1f} | {r['c_m']:+.4f} | {r['s_x']:+.4f} | {r['s_y']:+.4f} | "
        f"{r['sig_sx_lsq']:.4f}/{r['sig_sx_boot']:.4f} | {r['sig_sy_lsq']:.4f}/{r['sig_sy_boot']:.4f} | "
        f"{r['nmad_y_before']:.4f}->{r['nmad_y_after']:.4f} | {r['nmad_dh_before']:.4f}->{r['nmad_dh_after']:.4f}")

# ============================================================================= step 3
say("\n=== STEP 3: median dh by slope band, pooled and by aspect quadrant (non-floodplain) ===")
say("'after fit' = gen1 x,y moved by the fitted correction (-s_x,-s_y) and gen2 RESAMPLED (not a linear removal)")
band_dfs.append(band_table(dh, sl, asp, "raw, no shift", good))
band_dfs.append(band_table(dhr, slr, aspr, "registered (star+lateral), before fit", goodr))
# after pooled fit, resampled
sx, sy = fit_reg_pooled["s_x"], fit_reg_pooled["s_y"]
dha, sla, aspa, gooda = evaluate(REF, xr - sx, yr - sy, zr)
band_dfs.append(band_table(dha, sla, aspa, f"registered + pooled fitted shift ({-sx:+.3f},{-sy:+.3f}) resampled", gooda))
# linear removal, for comparison with the resampled version
dh_lin = dhr - np.tan(np.radians(slr)) * (sy * np.cos(np.radians(aspr)) + sx * np.sin(np.radians(aspr)))
band_dfs.append(band_table(dh_lin, slr, aspr, "registered + pooled fit removed LINEARLY (check)", goodr))
# after per-swath fits, resampled
xs, ys = xr.copy(), yr.copy()
for s, f in fit_reg_swath.items():
    m = psid_raw == s
    xs[m] -= f["s_x"]; ys[m] -= f["s_y"]
dhs, sls, asps, goods = evaluate(REF, xs, ys, zr)
band_dfs.append(band_table(dhs, sls, asps, "registered + PER-SWATH fitted shifts resampled", goods))
# also the pooled residual after removing c*tan(slope) is not reported: c absorbs the vertical offset
pd.concat(band_dfs).to_csv(f"{OUT}/slope_band_medians.csv", index=False)

# the residual-vs-slope test, stated as numbers: pooled median dh in 0-5 vs 30-35 vs 35-45, after fit
for df in band_dfs:
    lab = df["case"].iloc[0]
    flat = df.loc[df.band == "0-5", "pooled"].iloc[0]
    say(f"   {lab}: pooled median dh 0-5 deg {flat:+.4f}; minus that: "
        + ", ".join(f"{b} {v - flat:+.4f}" for b, v in zip(df.band, df.pooled)))
    for q in QUADS:
        flatq = df.loc[df.band == "0-5", q].iloc[0]
        say(f"      {q}: " + ", ".join(f"{b} {v - flatq:+.4f}" for b, v in zip(df.band, df[q])))

# ============================================================================= step 4
say("\n=== STEP 4: the bluff window rows %s, cols %s ===" % (f"{WINDOW[0].start}:{WINDOW[0].stop}", f"{WINDOW[1].start}:{WINDOW[1].stop}"))
wx0 = X0 + WINDOW[1].start * RES; wx1 = X0 + WINDOW[1].stop * RES
wy0 = Y0 + WINDOW[0].start * RES; wy1 = Y0 + WINDOW[0].stop * RES
say(f"window x [{wx0}, {wx1}) y [{wy0}, {wy1}) = {WINDOW[0].stop-WINDOW[0].start}x{WINDOW[1].stop-WINDOW[1].start} cells")
Wsl = REF.slope[WINDOW]; Wasp = REF.aspect[WINDOW]; Wfp = FP[WINDOW]
say(f"  gen2 slope in window: median {np.nanmedian(Wsl):.1f} deg, p10/p90 {np.nanpercentile(Wsl,10):.1f}/{np.nanpercentile(Wsl,90):.1f}; "
    f"aspect median {np.nanmedian(Wasp):.1f} deg; floodplain cells {int(Wfp.sum())} of {Wfp.size}")
steep = Wsl > 25
say(f"  cells with slope > 25 deg: {int(np.nansum(steep))}, their median aspect {np.nanmedian(Wasp[steep]):.1f} deg")


def in_window(x, y):
    iy, ix, ok = cell_index(x, y)
    return ok & (iy >= WINDOW[0].start) & (iy < WINDOW[0].stop) & (ix >= WINDOW[1].start) & (ix < WINDOW[1].stop)


wm = in_window(x_raw, y_raw)
say(f"  gen1 raw last returns in window: {int(wm.sum())}; per swath "
    f"{ {int(k): int(v) for k, v in zip(*np.unique(psid_raw[wm], return_counts=True))} }")
wmc = in_window(x_csf, y_csf)
say(f"  CSF-cache points in window: {int(wmc.sum())}; per swath "
    f"{ {int(k): int(v) for k, v in zip(*np.unique(psid_csf[wmc], return_counts=True))} }")
# swaths per cell in the window (raw last returns)
iy, ix, ok = cell_index(x_raw, y_raw)
cells_by_swath = {}
for s in sorted(STAR):
    m = wm & (psid_raw == s)
    cells_by_swath[s] = int(len(np.unique(iy[m] * NX + ix[m])))
say(f"  window cells (of {Wfp.size}) holding at least one raw last return, per swath: {cells_by_swath}")
# steep window cells by swath
steep_cells = set((np.nonzero(steep)[0] + WINDOW[0].start) * NX + np.nonzero(steep)[1] + WINDOW[1].start)
for s in sorted(STAR):
    m = wm & (psid_raw == s)
    say(f"    swath {s}: raw last returns on window cells steeper than 25 deg: "
        f"{int(np.isin(iy[m] * NX + ix[m], list(steep_cells)).sum())}")


def window_bands(dh_, sl_, label, mask):
    rows = []
    out = []
    for lo, hi in zip(SLOPE_BANDS[:-1], SLOPE_BANDS[1:]):
        m = mask & np.isfinite(dh_) & (sl_ >= lo) & (sl_ < hi)
        v = float(np.median(dh_[m])) if m.any() else np.nan
        rows.append(dict(case=label, band=f"{lo}-{hi}", n=int(m.sum()), median_dh=v))
        out.append(f"{lo}-{hi}: {v:+.4f} (n={int(m.sum())})")
    say(f"-- window median dh (m) by gen2 slope band, {label}:  " + " | ".join(out))
    return pd.DataFrame(rows)


wrows = []
wrows.append(window_bands(dh, sl, "raw, no shift", wm & good))
wrows.append(window_bands(dhr, slr, "registered (star+lateral)", wm & goodr))
wrows.append(window_bands(dha, sla, f"registered + pooled fitted shift ({-sx:+.3f},{-sy:+.3f})", wm & gooda))
wrows.append(window_bands(dhs, sls, "registered + per-swath fitted shifts", wm & goods))
# per-swath inside the window, registered
for s in sorted(STAR):
    m = wm & goodr & (psid_raw == s)
    if m.any():
        wrows.append(window_bands(dhr, slr, f"registered, swath {s} only", m))
# raw + raw-pooled fit
sx_i, sy_i = fits[0]["s_x"], fits[0]["s_y"]
wrows.append(window_bands(dh_c, sl_c, f"raw + raw-pooled fitted shift ({-sx_i:+.3f},{-sy_i:+.3f})", wm & good_c))
pd.concat(wrows).to_csv(f"{OUT}/window_slope_band_medians.csv", index=False)

# the window-implied east shift if the whole steep-cell offset were horizontal
m = wm & goodr & (slr > 25)
yv = dhr[m] / np.tan(np.radians(slr[m]))
say(f"  window, registered, slope > 25 deg: n={int(m.sum())}, median dh {np.median(dhr[m]):+.4f} m, "
    f"median aspect {np.median(aspr[m]):.1f} deg, median dh/tan(slope) {np.median(yv):+.4f} m "
    f"(= implied s_x if purely horizontal and east-facing; correction would be {-np.median(yv):+.4f} m east)")
m0 = wm & goodr & (slr < 10)
say(f"  window, registered, slope < 10 deg: n={int(m0.sum())}, median dh {np.median(dhr[m0]):+.4f} m (the local vertical offset)")

# DoD proxy in the window: per-cell median over gen1 last returns of (gen2_bilinear - gen1) = -dh
def dod_proxy(dh_, mask):
    f = iy[mask] * NX + ix[mask]
    s = pd.Series(-dh_[mask]).groupby(f).median()
    g = np.full(NY * NX, np.nan); g[s.index.values] = s.values
    return g.reshape(NY, NX)[WINDOW]


DOD = np.load(f"{DER}/dod.npy")[WINDOW]
say("-- window DoD (gen2 - gen1, m): pipeline dod.npy vs the per-cell proxy (median over raw last returns of gen2_bilinear - gen1)")
say(f"   pipeline dod.npy: n finite {int(np.isfinite(DOD).sum())}, median {np.nanmedian(DOD):+.4f}, "
    f"steep(>25) cells median {np.nanmedian(DOD[steep]):+.4f} (n {int(np.isfinite(DOD[steep]).sum())}), "
    f"gentle(<10) median {np.nanmedian(DOD[Wsl < 10]):+.4f}")
proxies = {
    "raw, no shift": dod_proxy(dh, wm & good),
    "registered (star+lateral)": dod_proxy(dhr, wm & goodr),
    f"registered + pooled fitted shift": dod_proxy(dha, wm & gooda),
    "registered + per-swath fitted shifts": dod_proxy(dhs, wm & goods),
}
for k, g in proxies.items():
    say(f"   proxy {k}: n finite {int(np.isfinite(g).sum())}, median {np.nanmedian(g):+.4f}, "
        f"steep(>25) median {np.nanmedian(g[steep]):+.4f} (n {int(np.isfinite(g[steep]).sum())}), "
        f"gentle(<10) median {np.nanmedian(g[Wsl < 10]):+.4f}, "
        f"steep minus gentle {np.nanmedian(g[steep]) - np.nanmedian(g[Wsl < 10]):+.4f}")
    np.save(f"{OUT}/window_dod_proxy_{k.split(',')[0].split(' ')[0]}_{'fit' if 'fitted' in k else 'nofit'}{'_swath' if 'PER' in k.upper() and 'fitted' in k else ''}.npy", g)
np.save(f"{OUT}/window_dod_pipeline.npy", DOD)
np.save(f"{OUT}/window_slope_deg.npy", Wsl)
# both proxies and pipeline on the SAME cells (where both finite), steep
both = np.isfinite(DOD) & np.isfinite(proxies["registered (star+lateral)"]) & steep
say(f"   on the {int(both.sum())} steep cells where BOTH pipeline dod and registered proxy are finite: "
    f"pipeline median {np.nanmedian(DOD[both]):+.4f}, proxy median {np.nanmedian(proxies['registered (star+lateral)'][both]):+.4f}, "
    f"proxy after pooled fit {np.nanmedian(proxies['registered + pooled fitted shift'][both]):+.4f}, "
    f"after per-swath fit {np.nanmedian(proxies['registered + per-swath fitted shifts'][both]):+.4f}")

# ============================================================================= step 5 numbers
say("\n=== STEP 5: the deciding numbers ===")
r = fit_reg_pooled
say(f"registered pooled fit amplitude a = {r['a_m']:.4f} m (boot 95% [{r['a_boot_p2_5']:.4f}, {r['a_boot_p97_5']:.4f}]), "
    f"direction b = {r['b_deg']:.1f} deg, i.e. gen1 displaced ({r['s_x']:+.4f}, {r['s_y']:+.4f}) m relative to gen2")
say(f"a horizontal shift explaining the window's steep offset alone would need |s_x| ~ {abs(np.median(yv)):.3f} m")
say(f"window steep median dh, registered: {np.median(dhr[m]):+.4f}; after pooled fitted shift: "
    f"{np.median(dha[wm & gooda & (sla > 25)]):+.4f}; after per-swath shifts: {np.median(dhs[wm & goods & (sls > 25)]):+.4f}")
say(f"done {time.strftime('%Y-%m-%d %H:%M:%S')}")
_rep.write("```\n")
_rep.close()
