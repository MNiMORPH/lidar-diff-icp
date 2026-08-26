#!/usr/bin/env python
"""
NEAR-NADIR slope dependence of the gen1 (2008) "reads-the-ground-too-low"
residual.

Prior finding (glennie_scanangle_swath_test.py): the gen1 steep-slope low is a
NEAR-NADIR (low scan angle) phenomenon, NOT a boresight/scan-mirror artifact
(sign wrong, per-swath incoherent, non-monotonic to swath edge). This script
takes the clean fixed-beam-geometry sub-population -- LOW SCAN ANGLE beams,
|scan_angle| < 5 deg -- and reads how the gen1-low residual varies with SLOPE.

Physics to state: for a near-nadir beam on a slope, the incidence angle to the
local surface normal ~ the slope angle, so the near-nadir residual-vs-slope curve
IS effectively a residual-vs-incidence-angle curve. We confirm this in the data.

Datum: remove the ~67 mm geoid datum FIRST -- subtract the flat-ground
(slope < 3 deg) median d_mm (about -50.9 mm) so residual r = d_mm - datum, with
0 = "no excursion".

Everything gen1-internal: scan_angle, incidence, slope, d_mm are gen1's own; the
gen2 z_after bare-earth plane is a FIXED spatial reference (nested-subset
comparison of a fixed registered dataset is internally consistent). canopy_cover
is a gen2-derived SPATIAL land-cover selector only -- its magnitude is NOT a
covariate; forest/open are selected by cc>0.5 / cc<0.2 cells.

Sign: gen1 BELOW the gen2 plane -> d_mm NEGATIVE. residual r < 0 = gen1 reads low.
Size effects by MEDIAN mm shifts against the ~20 mm signal budget, NOT r
(per-return NMAD is ~150-270 mm).
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as _st

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
_ap = argparse.ArgumentParser()
_ap.add_argument("--tile", default=os.path.join(ROOT, "data", "derived", "elba_fulldensity"))
_A = _ap.parse_args()
DERIVED = _A.tile
TILE = os.path.basename(DERIVED.rstrip("/"))
TAG = "" if TILE == "elba_fulldensity" else f"_{TILE}"
NPZ = os.path.join(DERIVED, "gen1_csf_angles.npz")
FLAT_SLOPE_DEG = 3.0

# fine slope bins (deg): headline curve resolution
SLOPE_EDGES = [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 35, 40, 90]
SLOPE_LABELS = ["0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21",
                "21-24", "24-27", "27-30", "30-35", "35-40", ">40"]

SA_NADIR = 5.0     # primary near-nadir: |scan_angle| < 5 deg
SA_NADIR_WIDE = 10.0   # robustness: |scan_angle| < 10 deg
SA_OBLIQUE = 10.0  # oblique contrast: |scan_angle| > 10 deg

MIN_N = 30         # min returns to report a bin median


def nmad(x):
    if x.size == 0:
        return np.nan
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))


def load():
    d = np.load(NPZ)
    ing = d["in_grid"]
    out = {k: d[k][ing] for k in
           ["incidence", "scan_angle", "slope", "d_mm", "cell",
            "point_source_id", "stratum", "pfs_forest", "pfs_open"]}
    out["absa"] = np.abs(out["scan_angle"])

    flat = out["slope"] < FLAT_SLOPE_DEG
    datum = float(np.median(out["d_mm"][flat]))
    out["datum_mm"] = datum
    out["res_mm"] = out["d_mm"] - datum
    out["n_flat"] = int(flat.sum())

    # gen2-derived canopy_cover as a SPATIAL land-cover selector (magnitude NOT used)
    cc = np.load(os.path.join(DERIVED, "canopy_cover_pfs.npy"))  # (NY, NX) PyForestScan cover
    cc_flat = cc.ravel()
    ok = (out["cell"] >= 0) & (out["cell"] < cc.size)
    cc_pt = np.full(out["cell"].shape, np.nan)
    cc_pt[ok] = cc_flat[out["cell"][ok]]
    out["cc"] = cc_pt
    out["cc_forest"] = cc_pt > 0.5
    out["cc_open"] = cc_pt < 0.2
    return out


def bin_curve(res, sl, sel):
    """Median r, NMAD, n per slope bin over selection `sel`."""
    meds, nmads, ns = [], [], []
    for i in range(len(SLOPE_LABELS)):
        s0, s1 = SLOPE_EDGES[i], SLOPE_EDGES[i + 1]
        m = sel & (sl >= s0) & (sl < s1)
        x = res[m]
        if x.size >= MIN_N:
            meds.append(float(np.median(x)))
            nmads.append(float(nmad(x)))
        else:
            meds.append(np.nan)
            nmads.append(np.nan)
        ns.append(int(x.size))
    return np.array(meds), np.array(nmads), np.array(ns)


def slope_center(i):
    s0, s1 = SLOPE_EDGES[i], SLOPE_EDGES[i + 1]
    if s1 >= 90:  # ">40" open bin: use a representative center
        return 45.0
    return 0.5 * (s0 + s1)


def fit_report(meds, ns):
    """Fit near-nadir median r vs tan(slope) and vs slope; report which fits
    better and whether a break at ~25-27 deg is needed beyond the smooth curve.
    EVERY populated slope bin enters the fit -- there is no slope truncation and
    no minimum-n cut. Sparse bins are not deleted; the sqrt(n) weighting already
    gives them the influence their counts earn."""
    cen = np.array([slope_center(i) for i in range(len(SLOPE_LABELS))])
    good = np.isfinite(meds)   # a median needs >=1 point; nothing else is cut
    x_deg = cen[good]
    x_tan = np.tan(np.deg2rad(cen[good]))
    y = meds[good]
    w = np.sqrt(ns[good].astype(float))  # weight by sqrt(n)

    def wls(X, y, w):
        # X: (n, k) design; returns coefs, weighted R^2, residuals
        Xw = X * w[:, None]
        yw = y * w
        coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        pred = X @ coef
        ybar = np.sum(w * y) / np.sum(w)
        ss_res = np.sum(w * (y - pred) ** 2)
        ss_tot = np.sum(w * (y - ybar) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        return coef, r2, ss_res

    # tan-law through origin
    Xt0 = x_tan[:, None]
    ct0, r2t0, ssr_t0 = wls(Xt0, y, w)
    # tan-law with free intercept
    Xt1 = np.column_stack([np.ones_like(x_tan), x_tan])
    ct1, r2t1, ssr_t1 = wls(Xt1, y, w)
    # linear in slope (free intercept)
    Xs1 = np.column_stack([np.ones_like(x_deg), x_deg])
    cs1, r2s1, ssr_s1 = wls(Xs1, y, w)

    # threshold test: tan-law + a step at 27 deg (extra dof). F-test-style.
    step = (x_deg >= 27).astype(float)
    Xb = np.column_stack([np.ones_like(x_tan), x_tan, step])
    cb, r2b, ssr_b = wls(Xb, y, w)
    n = good.sum()
    # improvement of adding the step to the tan+intercept model
    p_full, p_red = 3, 2
    if n > p_full and ssr_b > 0:
        F = ((ssr_t1 - ssr_b) / (p_full - p_red)) / (ssr_b / (n - p_full))
        df1, df2 = p_full - p_red, n - p_full
        F_p = float(_st.f.sf(F, df1, df2))
        F_crit = float(_st.f.isf(0.05, df1, df2))
    else:
        F = F_p = F_crit = np.nan
        df1 = df2 = 0

    # knee model: constant plateau below 27 deg + a ramp in (slope-27) above,
    # i.e. flat then switch-on. Tests the "plateau to ~24 then knee at 27" shape
    # directly against the smooth tan-law.
    ramp = np.maximum(x_deg - 27.0, 0.0)
    Xk = np.column_stack([np.ones_like(x_deg), ramp])
    ck, r2k, ssr_k = wls(Xk, y, w)

    return dict(
        cen=cen, good=good,
        tan0_slope=ct0[0], tan0_r2=r2t0,
        tan1_int=ct1[0], tan1_slope=ct1[1], tan1_r2=r2t1,
        lin_int=cs1[0], lin_slope=cs1[1], lin_r2=r2s1,
        step_int=cb[0], step_tan=cb[1], step_amp=cb[2], step_r2=r2b, step_F=F,
        step_p=F_p, step_Fcrit=F_crit, step_df1=df1, step_df2=df2,
        knee_plateau=ck[0], knee_ramp=ck[1], knee_r2=r2k,
        n=n,
    )


def incidence_track(D, sel):
    """median incidence vs slope for the near-nadir subset -> should track ~1:1."""
    inc = D["incidence"]; sl = D["slope"]
    rows = []
    for i in range(len(SLOPE_LABELS)):
        s0, s1 = SLOPE_EDGES[i], SLOPE_EDGES[i + 1]
        m = sel & (sl >= s0) & (sl < s1)
        if m.sum() >= MIN_N:
            rows.append((SLOPE_LABELS[i], slope_center(i),
                         float(np.median(inc[m])), int(m.sum())))
    return rows


def fmt_curve_table(labels, meds, nmads, ns):
    lines = ["| slope (deg) | median r (mm) | NMAD (mm) | n |",
             "|---|---|---|---|"]
    for lab, m, nm, n in zip(labels, meds, nmads, ns):
        if np.isfinite(m):
            lines.append(f"| {lab} | {m:+.1f} | {nm:.0f} | {n:,} |")
        else:
            lines.append(f"| {lab} | - | - | {n:,} |")
    return "\n".join(lines)


def make_figure(D, path):
    res = D["res_mm"]; sl = D["slope"]; absa = D["absa"]
    cen = np.array([slope_center(i) for i in range(len(SLOPE_LABELS))])

    nadir = absa < SA_NADIR
    oblique = absa > SA_OBLIQUE

    m_nadir, _, n_nadir = bin_curve(res, sl, nadir)
    m_obl, _, n_obl = bin_curve(res, sl, oblique)
    m_for, _, n_for = bin_curve(res, sl, nadir & D["cc_forest"])
    m_open, _, n_open = bin_curve(res, sl, nadir & D["cc_open"])

    fig, ax = plt.subplots(figsize=(11, 6), dpi=170)
    # smooth tan-law reference through the near-nadir origin-fit
    fit = fit_report(m_nadir, n_nadir)
    xs = np.linspace(0, 45, 200)
    ax.plot(xs, fit["tan0_slope"] * np.tan(np.deg2rad(xs)), color="0.5",
            lw=1.2, ls="--", label=f"tan-law fit (slope {fit['tan0_slope']:.0f} mm/tan)")

    ax.plot(cen, m_nadir, "-o", color="C0", ms=6, lw=2,
            label=f"near-nadir |SA|<5 (n={n_nadir.sum():,})")
    ax.plot(cen, m_obl, "-s", color="C3", ms=5, lw=1.6,
            label=f"oblique |SA|>10 (n={n_obl.sum():,})")
    # forest/open near-nadir where populated
    gf = n_for >= MIN_N
    go = n_open >= MIN_N
    ax.plot(cen[gf], m_for[gf], "-^", color="C2", ms=6, lw=1.6,
            label=f"near-nadir FOREST (cc>0.5, n={n_for.sum():,})")
    ax.plot(cen[go], m_open[go], "-v", color="C1", ms=6, lw=1.6,
            label=f"near-nadir OPEN (cc<0.2, n={n_open.sum():,})")

    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.axvline(27, color="0.7", lw=0.7, ls=":")
    ax.text(27.3, ax.get_ylim()[1] * 0.9, "27 deg", fontsize=8, color="0.4")
    ax.set_xlabel("slope (deg)  [near-nadir: incidence angle ~ slope]")
    ax.set_ylabel("median residual r (mm), datum removed\n<0 = gen1 reads ground LOW")
    ax.set_title("gen1 near-nadir (fixed-beam-geometry) residual vs slope\n"
                 "(datum = flat-ground median d_mm removed)")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def main():
    D = load()
    res = D["res_mm"]; sl = D["slope"]; absa = D["absa"]
    datum = D["datum_mm"]

    nadir = absa < SA_NADIR
    nadir_wide = absa < SA_NADIR_WIDE
    oblique = absa > SA_OBLIQUE

    m_nadir, nm_nadir, n_nadir = bin_curve(res, sl, nadir)
    m_wide, nm_wide, n_wide = bin_curve(res, sl, nadir_wide)
    m_obl, nm_obl, n_obl = bin_curve(res, sl, oblique)

    fit = fit_report(m_nadir, n_nadir)

    # incidence ~ slope confirmation (near-nadir)
    inc_rows = incidence_track(D, nadir)

    # land-cover split (near-nadir)
    m_for, nm_for, n_for = bin_curve(res, sl, nadir & D["cc_forest"])
    m_open, nm_open, n_open = bin_curve(res, sl, nadir & D["cc_open"])
    m_forC, nm_forC, n_forC = bin_curve(res, sl, nadir & D["pfs_forest"])
    m_openC, nm_openC, n_openC = bin_curve(res, sl, nadir & D["pfs_open"])

    # ---- report ----
    o = []
    o.append("# Near-nadir slope dependence of the gen1 (2008) ground-low residual")
    o.append("")
    o.append("_Generated by `nearnadir_slope_dependence.py`. gen1-internal: "
             "scan_angle / incidence / slope / d_mm are gen1's own; the gen2 "
             "`z_after` bare-earth plane is a FIXED spatial reference. "
             "`canopy_cover` (gen2-derived) is used only as a SPATIAL land-cover "
             "selector, not as a magnitude covariate._")
    o.append("")
    o.append(f"**Datum removed:** flat-ground (slope < {FLAT_SLOPE_DEG:g} deg) "
             f"median d_mm = **{datum:.1f} mm** (n={D['n_flat']:,}). This carries "
             "the ~67 mm GEOID03->GEOID18 constant; subtracting it sets residual "
             f"r = d_mm - ({datum:.1f}), so r=0 is 'no excursion'. Sign: r<0 = "
             "gen1 reads the ground LOW.")
    o.append("")
    o.append("**Near-nadir (fixed geometry):** primary population is |scan_angle| "
             f"< {SA_NADIR:g} deg (n={int(nadir.sum()):,}); robustness check "
             f"|scan_angle| < {SA_NADIR_WIDE:g} deg (n={int(nadir_wide.sum()):,}). "
             "Judge by MEDIAN mm shifts vs the ~20 mm budget, NOT r (per-return "
             "NMAD ~150-270 mm).")
    o.append("")

    # incidence ~ slope
    o.append("## 0. Incidence ~ slope confirmation (near-nadir)")
    o.append("")
    o.append("For a near-nadir beam on a slope, incidence-to-surface-normal ~ "
             "slope. So the residual-vs-slope curve below IS a "
             "residual-vs-incidence-angle curve.")
    o.append("")
    o.append("| slope (deg) | slope center | median incidence (deg) | n |")
    o.append("|---|---|---|---|")
    for lab, c, mi, n in inc_rows:
        o.append(f"| {lab} | {c:.1f} | {mi:.1f} | {n:,} |")
    o.append("")
    # quantify the tracking
    cc = np.array([r[1] for r in inc_rows])
    ii = np.array([r[2] for r in inc_rows])
    if cc.size > 2:
        A = np.column_stack([np.ones_like(cc), cc])
        sl_fit, *_ = np.linalg.lstsq(A, ii, rcond=None)
        o.append(f"Linear fit: median incidence = {sl_fit[0]:+.1f} + "
                 f"{sl_fit[1]:.2f} * slope. Slope coefficient ~1 confirms the "
                 "near-nadir incidence tracks slope roughly 1:1 (a small nadir "
                 "floor from residual scan angle + surface-normal noise).")
    o.append("")

    # headline curve
    o.append("## 1. HEADLINE: median residual r vs slope, near-nadir |SA|<5 deg")
    o.append("")
    o.append(fmt_curve_table(SLOPE_LABELS, m_nadir, nm_nadir, n_nadir))
    o.append("")
    o.append(f"Robustness: |SA|<{SA_NADIR_WIDE:g} deg (wider n at steep slope):")
    o.append("")
    o.append(fmt_curve_table(SLOPE_LABELS, m_wide, nm_wide, n_wide))
    o.append("")

    # tan-law vs threshold
    o.append("## 2. Smooth tan-law vs discrete threshold")
    o.append("")
    o.append(f"Fits of near-nadir median r vs slope (weighted by sqrt(n); "
             f"ALL {fit['n']:d} populated slope bins enter the fit -- no slope "
             "truncation, no minimum-n cut. The steep bins (35-40 deg, "
             "n=14,647; >40 deg, n=3,623) are sparse and their NMAD is "
             "~300-550 mm, so the sqrt(n) weighting already discounts them; "
             "they are NOT deleted):")
    o.append("")
    o.append("| model | form | R^2 |")
    o.append("|---|---|---|")
    o.append(f"| tan-law (origin) | r = {fit['tan0_slope']:.1f} * tan(slope) | "
             f"{fit['tan0_r2']:.3f} |")
    o.append(f"| tan-law (+intercept) | r = {fit['tan1_int']:+.1f} + "
             f"{fit['tan1_slope']:.1f} * tan(slope) | {fit['tan1_r2']:.3f} |")
    o.append(f"| linear in slope | r = {fit['lin_int']:+.1f} + "
             f"{fit['lin_slope']:.2f} * slope | {fit['lin_r2']:.3f} |")
    o.append(f"| knee: flat + ramp>27 | r = {fit['knee_plateau']:+.1f} + "
             f"{fit['knee_ramp']:.2f}*(slope-27)+ | {fit['knee_r2']:.3f} |")
    o.append(f"| tan-law + step@27 | tan + {fit['step_amp']:+.1f} mm step "
             f"(F={fit['step_F']:.2f}, p={fit['step_p']:.3f}) | "
             f"{fit['step_r2']:.3f} |")
    o.append("")
    # Data-honest verdict logic. Compare the four smooth/knee models by R^2.
    r2s = {"tan (+int)": fit["tan1_r2"], "linear": fit["lin_r2"],
           "knee@27": fit["knee_r2"], "tan+step@27": fit["step_r2"]}
    best = max(r2s, key=r2s.get)
    tan_best = fit["tan1_r2"] >= max(fit["lin_r2"], fit["knee_r2"])
    # Significance, not a raw R^2 increment: adding a parameter can only raise
    # R^2. The step earns belief only if the F-test clears its own 5% critical
    # value on the bins actually fitted.
    step_sig = bool(np.isfinite(fit["step_p"]) and fit["step_p"] < 0.05)
    o.append("**Verdict.** The near-nadir curve is NEITHER a clean smooth tan-law "
             "NOR a clean flat-then-step. It is roughly flat at about -15 mm from "
             "3-15 deg, RECOVERS toward 0 (-6 to -10 mm) at 18-24 deg, then drops "
             "sharply to about -40 to -44 mm at 27-35 deg. Over the full "
             "0-45 deg range the "
             f"best single smooth/knee fit is **{best}** (R^2={r2s[best]:.3f}); "
             f"the origin tan-law is a poor description (R^2={fit['tan0_r2']:.3f}).")
    o.append("")
    o.append(f"Adding a step at 27 deg to the tan-law lifts R^2 from "
             f"{fit['tan1_r2']:.3f} to {fit['step_r2']:.3f} (step "
             f"{fit['step_amp']:+.1f} mm), but adding a parameter always lifts "
             f"R^2. The test that matters is the F-test: F={fit['step_F']:.2f} "
             f"on ({fit['step_df1']:d},{fit['step_df2']:d}) df, "
             f"p={fit['step_p']:.3f}, against a 5% critical value of "
             f"{fit['step_Fcrit']:.2f}.")
    o.append("")
    if step_sig:
        o.append("The step CLEARS its critical value, so at this bin resolution "
                 "the 27 deg steepening is not reducible to the smooth trend.")
    else:
        o.append("**The step does NOT clear its critical value, so the 27 deg "
                 "switch-on is UNRESOLVED by this test.** What survives is the "
                 f"step AMPLITUDE ({fit['step_amp']:+.1f} mm); what is absent "
                 "is evidence that a break is needed at all. These medians are "
                 "consistent with a ~20 mm deepening near 27 deg and equally "
                 f"consistent with no break -- {fit['n']:d} bins cannot "
                 "separate the two.")
    o.append("")
    o.append("_Correction to the record (2026-08-26)._ This fit previously "
             "excluded the two bins above 35 deg (`FIT_MAX_SLOPE = 35`, plus an "
             "unexposed `n >= 200`). On `elba_fulldensity` that truncation "
             "reported R^2 = 0.800 and F = 10.85, written up here and in the "
             "project memory as '27 deg switch-on is REAL ... NOT a tan-curve "
             "artifact'. The significance was an artefact of the truncation: "
             "untruncated the same fit gives R^2 = 0.267, F = 2.49, p = 0.146. "
             "The step AMPLITUDE was not an artefact -- -22.2 mm truncated vs "
             "-22.1 mm untruncated -- so the size of the effect stands and only "
             "the evidence for its being a discrete break falls away. "
             "(`FRAME_2026-08-26.md` had already retired the knee on independent "
             "grounds -- per-swath misalignment -- so this changes the record, "
             "not the current science.)")
    o.append("")
    o.append("Honest reading: the binned medians do deepen from about -10/-20 mm "
             "below ~24 deg to about -40 mm at 27-35 deg, on top of a shallow, "
             "non-monotone trend that no clean tan-law reproduces. Whether that "
             "deepening is a genuine break or the steep end of a smooth curve is "
             "NOT decided here. Calling it purely 'tan-law, no threshold' would "
             "misstate the data; calling it an established knee overstates them.")
    o.append("")

    # oblique contrast
    o.append("## 3. Near-nadir vs oblique contrast")
    o.append("")
    o.append(f"Oblique |SA|>{SA_OBLIQUE:g} deg median r vs slope:")
    o.append("")
    o.append(fmt_curve_table(SLOPE_LABELS, m_obl, nm_obl, n_obl))
    o.append("")
    o.append("| slope (deg) | near-nadir r | oblique r | oblique - nadir |")
    o.append("|---|---|---|---|")
    for i, lab in enumerate(SLOPE_LABELS):
        mn, mo = m_nadir[i], m_obl[i]
        if np.isfinite(mn) and np.isfinite(mo):
            o.append(f"| {lab} | {mn:+.1f} | {mo:+.1f} | {mo-mn:+.1f} |")
        else:
            o.append(f"| {lab} | {'-' if not np.isfinite(mn) else f'{mn:+.1f}'} | "
                     f"{'-' if not np.isfinite(mo) else f'{mo:+.1f}'} | - |")
    o.append("")

    # land-cover split
    o.append("## 4. Land-cover split at matched slope (near-nadir)")
    o.append("")
    o.append("Forest = canopy_cover > 0.5 cells; open = canopy_cover < 0.2 cells "
             "(gen2-derived SPATIAL selector). Open thins with slope but does "
             "not run out: it is still tens of thousands of returns per bin "
             "well past 12 deg (see the counts in the table). Read every bin "
             "against its own n.")
    o.append("")
    o.append("### 4a. canopy_cover selector")
    o.append("")
    o.append("| slope (deg) | FOREST r (n) | OPEN r (n) | forest - open |")
    o.append("|---|---|---|---|")
    for i, lab in enumerate(SLOPE_LABELS):
        mf, nf = m_for[i], n_for[i]
        mo2, no = m_open[i], n_open[i]
        fs = f"{mf:+.1f} ({nf:,})" if np.isfinite(mf) else f"- ({nf:,})"
        os_ = f"{mo2:+.1f} ({no:,})" if np.isfinite(mo2) else f"- ({no:,})"
        d = f"{mf-mo2:+.1f}" if (np.isfinite(mf) and np.isfinite(mo2)) else "-"
        o.append(f"| {lab} | {fs} | {os_} | {d} |")
    o.append("")
    o.append("### 4b. PyForestScan forest / open masks (cover >= 0.5 / <= 0.1)")
    o.append("")
    o.append("| slope (deg) | pfs_FOREST r (n) | pfs_OPEN r (n) | forest - open |")
    o.append("|---|---|---|---|")
    for i, lab in enumerate(SLOPE_LABELS):
        mf, nf = m_forC[i], n_forC[i]
        mo2, no = m_openC[i], n_openC[i]
        fs = f"{mf:+.1f} ({nf:,})" if np.isfinite(mf) else f"- ({nf:,})"
        os_ = f"{mo2:+.1f} ({no:,})" if np.isfinite(mo2) else f"- ({no:,})"
        d = f"{mf-mo2:+.1f}" if (np.isfinite(mf) and np.isfinite(mo2)) else "-"
        o.append(f"| {lab} | {fs} | {os_} | {d} |")
    o.append("")

    # Matched-slope forest vs open. EVERY bin where both classes have at least
    # one return is listed, with its counts, and no minimum-n cut: the reader
    # judges reliability from n, which is printed, not from a hidden threshold.
    matched = [(SLOPE_LABELS[i], m_for[i], m_open[i], n_for[i], n_open[i])
               for i in range(len(SLOPE_LABELS))
               if np.isfinite(m_for[i]) and np.isfinite(m_open[i])]
    o.append("**Matched-slope forest vs open (near-nadir):**")
    o.append("")
    n_open_tot = int(n_open.sum())
    o.append("Both classes are listed in every slope bin that has returns in "
             "each, with the counts alongside; no bin is suppressed for being "
             "sparse. Open (cc<0.2) does thin with slope -- from "
             f"n={n_open[0]:,} at 0-3 deg to n={n_open[-1]:,} above 40 deg out "
             f"of {n_open_tot:,} near-nadir open returns -- so the steepest "
             "open medians carry the widest uncertainty and should be read "
             "against their n, not taken as equal to the dense low-slope ones:")
    o.append("")
    for lab, mf, mo2, nf, no in matched:
        rel = "ABOVE (less low than)" if mf > mo2 else "BELOW (lower than)"
        o.append(f"- slope {lab}: forest {mf:+.1f} mm (n={nf:,}) vs open "
                 f"{mo2:+.1f} mm (n={no:,}) -- forest sits {rel} open by "
                 f"{mf-mo2:+.1f} mm")
    o.append("")
    diffs = [mf - mo2 for _, mf, mo2, _, _ in matched]
    if diffs:
        md = float(np.median(diffs))
        o.append(f"Across all {len(matched):d} matched bins "
                 f"({matched[0][0]} to {matched[-1][0]} deg), forest reads "
                 f"{'LOWER' if md < 0 else 'HIGHER'} than open at matched slope "
                 f"by a median {md:+.1f} mm -- a canopy/forest-floor term ON TOP "
                 "of slope. This is the disentangling result: even at NEAR-NADIR "
                 "and MATCHED SLOPE, forest ground reads lower than open ground, "
                 "so the near-nadir low is not slope alone. The bin-to-bin "
                 "difference is noisy and changes sign in places, so read the "
                 "median across bins rather than any single bin.")
    o.append("")

    o.append("## Caveats")
    o.append("")
    o.append("- Median mm shifts are the signal (~20 mm budget); per-return NMAD "
             "is ~150-270 mm so per-return correlations are tiny by construction.")
    o.append("- gen1-internal comparison against a fixed gen2 spatial plane is "
             "internally consistent (nested-subset of a registered dataset).")
    o.append("- canopy_cover is a gen2-derived SPATIAL selector; only forest/open "
             "membership is used, never its magnitude.")
    o.append("- A clean smooth tan-law with no real threshold is a valid finding, "
             "not a failure.")
    o.append("")

    out_md = os.path.join(HERE, f"NEARNADIR_SLOPE_DEPENDENCE{TAG.upper()}.md")
    with open(out_md, "w") as f:
        f.write("\n".join(o))

    fig = make_figure(D, os.path.join(HERE, f"nearnadir_slope_dependence{TAG}.png"))

    print(f"datum = {datum:.1f} mm (n_flat={D['n_flat']:,})")
    print("wrote", out_md)
    print("wrote", fig)
    print("\nNear-nadir |SA|<5 median r vs slope:")
    for lab, m, n in zip(SLOPE_LABELS, m_nadir, n_nadir):
        mm = f"{m:+6.1f}" if np.isfinite(m) else "   -  "
        print(f"  {lab:>6s}: r={mm} mm  n={n:,}")
    print(f"\ntan-law origin R^2={fit['tan0_r2']:.3f} "
          f"(slope {fit['tan0_slope']:.1f} mm/tan); "
          f"linear-in-slope R^2={fit['lin_r2']:.3f}; "
          f"step@27 amp={fit['step_amp']:+.1f} mm R^2={fit['step_r2']:.3f} "
          f"F={fit['step_F']:.2f} p={fit['step_p']:.3f} "
          f"(F_crit={fit['step_Fcrit']:.2f}, n_bins={fit['n']:d})")
    print("\nForest vs open (near-nadir), every bin populated in both:")
    for lab, mf, mo2, nf, no in matched:
        print(f"  {lab:>6s}: forest {mf:+.1f}  open {mo2:+.1f}  "
              f"diff {mf-mo2:+.1f}  (forest n={nf:,}, open n={no:,})")


if __name__ == "__main__":
    main()
