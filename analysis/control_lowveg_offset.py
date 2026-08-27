#!/usr/bin/env python3
"""gen2's ground-surface offset against surveyed control, as a function of near-ground vegetation.

Reproduces every number in `analysis/CONTROL_LOWVEG_OFFSET.md` from LOCAL data only -- no
network. Inputs are the products of `analysis/cover_at_control_marks.py`
(`control_cover.csv`, `control_structure/*.npz`) plus the bundled control CSV.

THE TWO QUANTITIES ARE INDEPENDENT BY CONSTRUCTION, which is the point of the design:

  density  from the CLOUD'S SHAPE -- the fraction of returns in a near-ground band above the
           order-2 least-squares surface fitted to class-2 returns within the mark's radius.
           Because that surface is fitted from the box's own returns, the histogram is
           INVARIANT to any vertical shift of the cloud, so this metric structurally cannot
           carry offset information.
  offset   from OUTSIDE the cloud -- USGS's published surveyed_Z minus delivered_LAZ_Z at the
           mark. Computed by an estimator we had no hand in. `+ve = the surface reads LOW`.

USE THE PUBLISHED RESIDUAL, NOT OUR OWN. Our least-squares surface diverges from the vendor's
by ~62 mm in vegetation and ~2 mm in the open -- i.e. AS A FUNCTION OF VEGETATION -- so using
our own offset would put a vegetation-dependent term on both axes.

THE METRIC IS ORDINAL, NOT ABSOLUTE, and its definition must travel with any coefficient
fitted from it. Moving the band's lower edge over +/-0.10 m changes the metric's value by a
factor of ~50 while leaving the rank correlation between -0.30 and -0.37: the ordering of
marks is stable, the scale is not. `--band-lo` is swept by `--sweep` for exactly this reason.

LCPs ARE EXCLUDED. The 143 LiDAR Control Points calibrated the acquisition; the NVA/VVA
checkpoints were held out. Checking gen2 against its own calibration points would be circular.
They carry `role=calibration` and no residual, so they drop out naturally -- asserted here.

    ./lidar-icp/bin/python analysis/control_lowveg_offset.py --out data/derived/control_lowveg_offset.csv
"""
import argparse, os, sys
import numpy as np, pandas as pd
from scipy import stats

COVER = "data/derived/control_cover.csv"
STRUCT = "data/derived/control_structure"
G2 = "src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv"


def lowveg(point_id, lo, hi, setname="gen2_2021_control"):
    """Fraction of ALL returns with slope-normal height in (lo, hi] above the local surface."""
    f = os.path.join(STRUCT, f"{setname}__{point_id}.npz")
    if not os.path.exists(f):
        return np.nan
    z = np.load(f)
    e = z["ng_edges"]; mid = 0.5 * (e[:-1] + e[1:])
    h = z["ng_all"].astype(float)
    t = h.sum()
    return h[(mid > lo) & (mid <= hi)].sum() / t if t else np.nan


def load(band_lo, band_hi):
    missing = [p for p in (COVER, STRUCT) if not os.path.exists(p)]
    if missing:
        raise SystemExit(
            f"missing {', '.join(missing)}.\n"
            "These are the products of the acquisition step, which reads the gen2 3DEP EPT\n"
            "over the network (~55 min for 1,497 marks, ~1.6 GB of kept boxes):\n"
            "    PROJ_DATA=$CENV/share/proj GDAL_DATA=$CENV/share/gdal $CENV/bin/python \\\n"
            "        analysis/cover_at_control_marks.py --out data/derived/control_cover.csv\n"
            "It is resumable: rerun it and it skips marks already present.")
    cov = pd.read_csv(COVER)
    cov = cov[(cov.set == "gen2_2021_control") & (cov.status == "ok")]
    g2 = pd.read_csv(G2)
    g2["resid_mm"] = g2[["usgs_ql1_laz_error_m", "usgs_ql0_laz_error_m"]].mean(axis=1, skipna=True) * 1000
    assert g2[g2.role == "calibration"].resid_mm.notna().sum() == 0, \
        "an LCP carries a residual: the held-out/calibration split is broken"
    m = cov.merge(g2[["point_id", "point_type", "resid_mm", "role", "elevation"]],
                  on="point_id", how="inner", suffixes=("", "_g2"))
    m["lowveg"] = [lowveg(p, band_lo, band_hi) for p in m.point_id]
    return m.dropna(subset=["lowveg", "resid_mm"]).copy()


def near_ground_stats(point_id, setname="gen2_2021_control"):
    """Spread of the near-ground returns at a mark: IQR, p90-p10 and NMAD, mm, for class-2
    and for all returns. Read off the stored histogram, so these are QUANTISED at the 20 mm
    bin width -- treat the ladder as ordinal. The boxes resolve it finer if ever needed."""
    f = os.path.join(STRUCT, f"{setname}__{point_id}.npz")
    if not os.path.exists(f):
        return None
    z = np.load(f); e = z["ng_edges"]; mid = 0.5 * (e[:-1] + e[1:])
    out = {}
    for key, tag in (("ng_class2", "c2"), ("ng_all", "all")):
        h = z[key].astype(float); t = h.sum()
        if t < 20:
            return None
        c = np.cumsum(h) / t
        q = lambda pp: mid[np.searchsorted(c, pp)]
        med = q(0.5)
        order = np.argsort(np.abs(mid - med)); cw = np.cumsum(h[order]) / t
        out[f"iqr_{tag}"] = (q(0.75) - q(0.25)) * 1000
        out[f"p9010_{tag}"] = (q(0.90) - q(0.10)) * 1000
        out[f"nmad_{tag}"] = 1.4826 * np.abs(mid - med)[order][np.searchsorted(cw, 0.5)] * 1000
    return out


def canopy_frac(point_id, lo, hi, setname="gen2_2021_control"):
    """Fraction of ALL returns in (lo, hi] m above the local surface, from the TALL window
    (-2..+45 m, 0.25 m bins) -- so the upper edge can be swept past the near-ground cube."""
    f = os.path.join(STRUCT, f"{setname}__{point_id}.npz")
    if not os.path.exists(f):
        return np.nan
    z = np.load(f); e = z["can_edges"]; mid = 0.5 * (e[:-1] + e[1:])
    h = z["can_all"].astype(float); t = h.sum()
    return h[(mid > lo) & (mid <= hi)].sum() / t if t else np.nan


def _partial(x, y, z_):
    """Pearson correlation of x and y after linearly removing z_ from both."""
    rx = x - np.poly1d(np.polyfit(z_, x, 1))(z_)
    ry = y - np.poly1d(np.polyfit(z_, y, 1))(z_)
    return stats.pearsonr(rx, ry)


LOWVEG_DEFINITION = """\
lowveg, EXACT DEFINITION -- it must travel with any coefficient fitted from it
--------------------------------------------------------------------------------
1. At the mark, read gen2's OWN acquisition from the covering MN_SEDriftless_*_2021
   EPT block: a 105 m box, every return, every class.
2. Fit an ORDER-2 least-squares surface S to the CLASS-2 returns within 7.5 m of the
   mark (7.5 m = 1.5 x the pipeline's 5 m grid, the tie estimator's own report radius).
   Order 2 removes slope AND curvature as a trend.
3. For every return within 7.5 m, take the SLOPE-NORMAL height above that surface:
       h = (z - S(x,y)) / sqrt(1 + gx^2 + gy^2),  gradient at the mark.
4. lowveg = the fraction of those returns with 0.15 m < h <= 2.00 m.
   Lower edge 0.15 m: ~2.5x the bare-ground class-2 NMAD (59.3 mm), so above the
   surface's own noise. Upper edge 2.00 m: below tree crowns. NEITHER IS PHYSICAL --
   see --sweep and --strata; the metric is ORDINAL, its scale moves ~50x with the
   lower edge while the rank correlation moves 0.07.
Because S comes from the box's own returns, lowveg is INVARIANT to any vertical shift
of the cloud: it structurally cannot carry the offset it is used to predict."""


def fit_origin(x, y, w=None):
    """Weighted least squares THROUGH THE ORIGIN: y = b*x. Returns (b, se)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    w = np.ones_like(x) if w is None else np.asarray(w, float)
    b = (w * x * y).sum() / (w * x * x).sum()
    r = y - b * x
    dof = max(len(x) - 1, 1)
    se = np.sqrt((w * r ** 2).sum() / dof / (w * x * x).sum())
    return b, se


def wls(x, y, w):
    X = np.c_[np.ones(len(x)), x]; W = np.diag(w)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
    r = y - X @ beta
    s2 = (w * r ** 2).sum() / max(len(x) - 2, 1)
    return beta, np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ W @ X)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band-lo", type=float, default=0.15,
                    help="lower edge of the near-ground band, m. NOT a physical constant: see "
                         "the module docstring and --sweep")
    ap.add_argument("--band-hi", type=float, default=2.0)
    ap.add_argument("--bin-width", type=float, default=0.06, help="UNIFORM bins in lowveg")
    ap.add_argument("--block-km", type=float, default=10.0, help="spatial block for the bootstrap")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--sweep", action="store_true", help="report the band-edge sensitivity")
    ap.add_argument("--scatter", action="store_true",
                    help="near-ground spread vs vegetation vs offset, with partial correlations")
    ap.add_argument("--slope-check", action="store_true",
                    help="size the vertical-vs-slope-normal mismatch and test slope as a confound")
    ap.add_argument("--strata", action="store_true",
                    help="sweep the band's UPPER edge and isolate each height stratum")
    ap.add_argument("--all", action="store_true", help="every section")
    ap.add_argument("--plot", default=None, help="write the regression figure here")
    ap.add_argument("--out", default=None, help="write the per-mark table")
    a = ap.parse_args()

    if a.all:
        a.sweep = a.scatter = a.slope_check = a.strata = True
    m = load(a.band_lo, a.band_hi)
    print(f"n = {len(m)} gen2 checkpoints (held-out NVA/VVA; LCPs excluded and asserted)")
    print(f"lowveg = fraction of returns in ({a.band_lo}, {a.band_hi}] m above the local surface")
    print(f"offset = USGS surveyed_Z - delivered_LAZ_Z, +ve = surface reads LOW\n")
    print(LOWVEG_DEFINITION + "\n")

    E = np.arange(0, m.lowveg.max() + a.bin_width, a.bin_width)
    rows = []
    for lo, hi in zip(E[:-1], E[1:]):
        s = m[(m.lowveg >= lo) & (m.lowveg < hi)]
        se = s.resid_mm.std(ddof=1) / np.sqrt(len(s)) if len(s) > 1 else np.nan
        rows.append((lo, hi, len(s), s.resid_mm.median() if len(s) else np.nan,
                     s.resid_mm.mean() if len(s) else np.nan, se))
    b = pd.DataFrame(rows, columns=["lo", "hi", "n", "median", "mean", "se"])
    print(f"{'bin':>13} {'n':>4} {'median':>8} {'mean':>8} {'SE':>7}")
    for _, r in b.iterrows():
        print(f"  {r.lo:.2f}-{r.hi:.2f} {int(r.n):4d} " +
              ("      -        -       -" if r.n == 0 else
               f"{r['median']:8.1f} {r['mean']:8.1f} {r.se:7.1f}"))

    g = b[(b.n > 1) & np.isfinite(b.se)].copy(); g["x"] = 0.5 * (g.lo + g.hi)
    # THROUGH-ORIGIN fit. Justified, not assumed: with no vegetation there is no
    # vegetation-induced bias, and both the free intercept here (-5.7 +/- 4.3 mm) and gen2's
    # own open-ground level against held-out control (NVA -2.22 +/- 2.35 mm) are consistent
    # with zero. Forcing the origin would ABSORB a real datum offset into the slope, so the
    # free-intercept fit is always printed beside it as the check.
    bo_bin, so_bin = fit_origin(g.x.values, g["mean"].values, 1 / g.se.values ** 2)
    bo_mk, so_mk = fit_origin(m.lowveg.values, m.resid_mm.values)
    print(f"\nTHROUGH-ORIGIN (offset = b * lowveg), the form to carry back to the DEM:")
    print(f"  binned, 1/SE^2 weighted : b = {bo_bin:+8.1f} +/- {so_bin:.1f} mm per unit lowveg")
    print(f"  per-mark, unweighted    : b = {bo_mk:+8.1f} +/- {so_mk:.1f} mm per unit lowveg")
    bd, sd_ = wls(g.x.values, g["mean"].values, 1 / g.se.values ** 2)
    ba, sa = wls(g.x.values, g["mean"].values, g.n.values.astype(float))
    print(f"\nfits on the binned means ({len(g)} bins with n>1):")
    print(f"  DESIGN-weighted  1/SE^2 : intercept {bd[0]:+7.1f} +/- {sd_[0]:.1f}   "
          f"slope {bd[1]:+8.1f} +/- {sd_[1]:.1f} mm per unit")
    print(f"  ABUNDANCE-weighted by n : intercept {ba[0]:+7.1f} +/- {sa[0]:.1f}   "
          f"slope {ba[1]:+8.1f} +/- {sa[1]:.1f} mm per unit")

    B = a.block_km * 1000.0
    m["blk"] = (m.easting // B).astype(int).astype(str) + "_" + (m.northing // B).astype(int).astype(str)
    ub = m.blk.unique(); rng = np.random.default_rng(0); sl = []
    for _ in range(a.n_boot):
        s = pd.concat([m[m.blk == k] for k in rng.choice(ub, size=len(ub), replace=True)])
        if s.lowveg.nunique() > 2:
            sl.append(np.polyfit(s.lowveg, s.resid_mm, 1)[0])
    sl = np.array(sl); nv = stats.linregress(m.lowveg, m.resid_mm)
    print(f"\nper-mark slope: naive {nv.slope:+.1f} +/- {nv.stderr:.1f} mm per unit (p {nv.pvalue:.1e})")
    print(f"  block bootstrap on {len(ub)} blocks of {a.block_km:.0f} km: "
          f"{np.mean(sl):+.1f} +/- {np.std(sl, ddof=1):.1f}   "
          f"-> SE inflated {np.std(sl, ddof=1)/nv.stderr:.2f}x")

    print(f"\nconfound check -- slope within each EPT block (they differ in FLIGHT DATE):")
    for k, s in m.groupby("ept_block"):
        if len(s) < 25 or s.lowveg.nunique() < 5:
            print(f"  {k:24s} n={len(s):4d}  too few"); continue
        r = stats.linregress(s.lowveg, s.resid_mm)
        print(f"  {k:24s} n={len(s):4d}  slope {r.slope:+8.1f} +/- {r.stderr:6.1f}  p {r.pvalue:.1e}")

    if a.sweep:
        print(f"\nband-edge sensitivity (the metric is ORDINAL; its scale is not meaningful):")
        print(f"{'lower edge':>11} {'median lowveg':>14} {'rho vs offset':>15} {'p':>10}")
        for d in (-0.10, -0.06, -0.02, 0.0, 0.02, 0.06, 0.10):
            v = np.array([lowveg(p, a.band_lo + d, a.band_hi) for p in m.point_id])
            ok = np.isfinite(v)
            r, p = stats.spearmanr(v[ok], m.resid_mm[ok])
            print(f"  {a.band_lo+d:9.2f} m {np.nanmedian(v):14.4f} {r:15.3f} {p:10.2e}")

    if a.scatter:
        rec = []
        for _, r in m.iterrows():
            st = near_ground_stats(r.point_id)
            if st:
                rec.append({**r.to_dict(), **st})
        t = pd.DataFrame(rec)
        print(f"\nNEAR-GROUND SCATTER (mm; quantised at the 20 mm bin width), n={len(t)}:")
        print(f"{'lowveg bin':>13} {'n':>4} {'iqr_c2':>8} {'nmad_c2':>8} {'iqr_all':>8} {'offset':>9}")
        E2 = np.arange(0, t.lowveg.max() + 0.09, 0.09)
        for lo, hi in zip(E2[:-1], E2[1:]):
            s_ = t[(t.lowveg >= lo) & (t.lowveg < hi)]
            if not len(s_):
                continue
            print(f"  {lo:.2f}-{hi:.2f} {len(s_):4d} {s_.iqr_c2.median():8.1f} "
                  f"{s_.nmad_c2.median():8.1f} {s_.iqr_all.median():8.1f} {s_.resid_mm.median():9.1f}")
        print(f"\n{'':22s} {'vs lowveg':>12} {'vs offset':>12}")
        for c in ("iqr_c2", "nmad_c2", "p9010_c2", "iqr_all", "nmad_all"):
            r1, _ = stats.spearmanr(t[c], t.lowveg); r2, p2 = stats.spearmanr(t[c], t.resid_mm)
            print(f"  {c:20s} {r1:+12.3f} {r2:+12.3f}   (p {p2:.1e})")
        r1, _ = stats.spearmanr(t.lowveg, t.resid_mm)
        print(f"  {'lowveg':20s} {'':12s} {r1:+12.3f}")
        pa, ppa = _partial(t.nmad_c2.values, t.resid_mm.values, t.lowveg.values)
        pb, ppb = _partial(t.lowveg.values, t.resid_mm.values, t.nmad_c2.values)
        print(f"\n  class-2 NMAD | controlling for lowveg : r {pa:+.3f}  p {ppa:.2e}")
        print(f"  lowveg       | controlling for NMAD   : r {pb:+.3f}  p {ppb:.2e}")
        print("  -> scatter is a SYMPTOM of vegetation, not an independent driver of the offset")

    if a.slope_check:
        sl = m.slope_deg.dropna()
        f = 1 / np.cos(np.radians(sl))
        print(f"\nSLOPE. The scatter is slope-normal by construction; the OFFSET is a VERTICAL")
        print(f"difference (USGS surveyed_Z - delivered_Z) and is not converted. Size of that:")
        print(f"  slope: median {sl.median():.2f} deg  p90 {sl.quantile(.9):.2f}  max {sl.max():.2f}")
        print(f"  1/cos(slope): median {f.median():.4f}  max {f.max():.4f}  "
              f"-> {100*(f.median()-1):.2f}% typical, {100*(f.max()-1):.1f}% worst")
        rv, _ = stats.spearmanr(m.lowveg, m.resid_mm)
        rn, _ = stats.spearmanr(m.lowveg, m.resid_mm / np.cos(np.radians(m.slope_deg)))
        print(f"  lowveg vs offset:  vertical rho {rv:+.3f}   slope-normal rho {rn:+.3f}")
        for c, lab in (("lowveg", "lowveg"), ("resid_mm", "offset")):
            r, p = stats.spearmanr(m.slope_deg, m[c])
            print(f"  slope vs {lab:8s} rho {r:+.3f}  p {p:.1e}")
        r, p = _partial(m.lowveg.values, m.resid_mm.values, m.slope_deg.values)
        print(f"  lowveg vs offset, controlling for slope: r {r:+.3f}  p {p:.2e}")

    if a.strata:
        print(f"\nUPPER-EDGE sweep (lower edge 0.25 m; tall window, 0.25 m bins):")
        print(f"{'band':>14} {'median frac':>12} {'rho vs offset':>15} {'p':>11}")
        for hi in (0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 45.0):
            v = np.array([canopy_frac(p_, 0.25, hi) for p_ in m.point_id]); ok = np.isfinite(v)
            r, p = stats.spearmanr(v[ok], m.resid_mm[ok])
            print(f"  0.25-{hi:5.1f} m {np.nanmedian(v):12.4f} {r:15.3f} {p:11.2e}")
        print(f"\nSTRATA in isolation -- which layer carries the signal:")
        print(f"{'band':>14} {'median frac':>12} {'rho vs offset':>15} {'p':>11}")
        for lo, hi in ((0.25, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 45.0)):
            v = np.array([canopy_frac(p_, lo, hi) for p_ in m.point_id]); ok = np.isfinite(v)
            r, p = stats.spearmanr(v[ok], m.resid_mm[ok])
            print(f"  {lo:5.2f}-{hi:5.1f} m {np.nanmedian(v):12.4f} {r:15.3f} {p:11.2e}")
        print("  NOTE: the tall strata have a median fraction of 0.0000 -- surveyors do not put")
        print("  marks under closed canopy -- so this cannot test tall canopy, only show that the")
        print("  near-ground layer alone reproduces the full signal.")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.6, 5.4))
        ax.axhline(0, color="0.8", lw=0.8, zorder=0)
        ax.scatter(m.lowveg, m.resid_mm, s=12, c="0.65", edgecolor="none", zorder=1,
                   label=f"{len(m)} checkpoints")
        ax.errorbar(g.x, g["mean"], yerr=g.se, fmt="o", ms=7, color="C0", capsize=3, zorder=3,
                    label="uniform bins, mean $\\pm$ SE")
        for _, r in g.iterrows():
            ax.annotate(f"{int(r.n)}", (r.x, r["mean"]), textcoords="offset points",
                        xytext=(7, 5), fontsize=7, color="C0")
        xs = np.linspace(0, m.lowveg.max() * 1.02, 50)
        ax.plot(xs, bo_bin * xs, "C3-", lw=2, zorder=4,
                label=f"through origin: ${bo_bin:.0f} \\pm {so_bin:.0f}$ mm per unit")
        ax.plot(xs, bd[0] + bd[1] * xs, "C3--", lw=1.2, zorder=4,
                label=f"free intercept: ${bd[0]:+.1f} \\pm {sd_[0]:.1f}$ mm")
        ax.set_xlabel(f"lowveg  =  fraction of returns {a.band_lo:.2f}–{a.band_hi:.2f} m "
                      "above the local ground surface")
        ax.set_ylabel("surveyed $-$ delivered lidar elevation  (mm)")
        ax.set_title("gen2 ground surface floats above true ground in proportion to low vegetation\n"
                     f"{len(m)} held-out NVA/VVA checkpoints, single epoch", fontsize=10)
        ax.text(0.985, 0.04, "negative = lidar reads HIGH", transform=ax.transAxes,
                ha="right", fontsize=8, color="0.35")
        ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
        fig.tight_layout(); os.makedirs("figures", exist_ok=True)
        fig.savefig(a.plot, dpi=160)
        print(f"\nwrote {a.plot}")

    if a.out:
        cols = ["point_id", "point_type_g2", "easting", "northing", "ept_block",
                "gps_utc_min", "lowveg", "resid_mm", "cover_r7.5", "n_ground", "slope_deg"]
        m[[c for c in cols if c in m.columns]].to_csv(a.out, index=False)
        print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
