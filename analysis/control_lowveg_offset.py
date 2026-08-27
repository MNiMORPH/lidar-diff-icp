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
    ap.add_argument("--out", default=None, help="write the per-mark table")
    a = ap.parse_args()

    m = load(a.band_lo, a.band_hi)
    print(f"n = {len(m)} gen2 checkpoints (held-out NVA/VVA; LCPs excluded and asserted)")
    print(f"lowveg = fraction of returns in ({a.band_lo}, {a.band_hi}] m above the local surface")
    print(f"offset = USGS surveyed_Z - delivered_LAZ_Z, +ve = surface reads LOW\n")

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

    if a.out:
        cols = ["point_id", "point_type_g2", "easting", "northing", "ept_block",
                "gps_utc_min", "lowveg", "resid_mm", "cover_r7.5", "n_ground", "slope_deg"]
        m[[c for c in cols if c in m.columns]].to_csv(a.out, index=False)
        print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
