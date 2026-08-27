#!/usr/bin/env python3
"""What height window, and what weighting inside it, best turns near-ground returns into a
predictor of gen2's ground-surface offset at surveyed control?

`analysis/control_lowveg_offset.py` fixes the predictor at a BOXCAR over (0.15, 2.0] m and
sweeps only its edges. This script sweeps the whole weight FUNCTION w(h), and repeats the
sweep within strata, because the optimum need not be the same in the open and under
vegetation.

    metric(mark) = sum_h w(h) * counts(h) / sum_h counts(h)

Height axis is SPLICED: the near-ground cube's 20 mm bins from -1 m to +2 m, then the tall
window's 250 mm bins above +2 m. Both are raw counts from the same box, so the ratio is a
genuine fraction of returns; the bin width changes at 2 m but w(h) is evaluated at bin
centres, so a wide window is not biased by the change.

TWO SCORINGS ARE REPORTED AND THEY ARE NOT COMPARABLE TO EACH OTHER:
  alone   the metric is the only predictor. This is the number quoted in
          `analysis/CONTROL_LOWVEG_OFFSET.md` (boxcar 0.15-2.0 scores +0.170).
  +block  the metric PLUS one intercept per EPT block. Each block is a separate acquisition
          with its own unpublished vendor bias, so this asks what the metric adds once that
          is absorbed. It scores ~0.02 higher throughout. Compare within a column only.

Scored on HELD-OUT 10 km spatial blocks, 5 folds, 20 fold seeds; mean and sd over seeds.

    ./lidar-icp/bin/python analysis/control_band_window_sweep.py --all
"""
import argparse, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import control_lowveg_offset as M

# Values recovered from the inline sweeps of 2026-08-27T18:15-18:19Z, kept so that any change
# to the loaders or the splice shows up as a mismatch instead of passing silently.
REFERENCE = {
    ("shape", "boxcar 0.15-2"): 0.1844, ("shape", "boxcar 0.15-3"): 0.1942,
    ("shape", "boxcar 0.15-5"): 0.1989, ("shape", "sqrt-h weight (0.15-3)"): 0.2018,
    ("shape", "exp decay lam=1.2"): 0.1751, ("shape", "gauss h0=0.25 sd=0.15"): 0.1167,
    ("upper", 2): 0.1844, ("upper", 4): 0.1993, ("upper", 8): 0.2007, ("upper", 20): 0.1910,
    ("lower", 0.05): 0.1249, ("lower", 0.15): 0.1989, ("lower", 0.30): 0.1456,
}


def build(min_returns=200):
    """Spliced raw-count profiles, the offset, coordinates, EPT block and stratifiers."""
    m = M.load(0.15, 2.0)
    NG, CAN, keep = [], [], []
    for i, r in m.iterrows():
        f = os.path.join(M.STRUCT, f"gen2_2021_control__{r.point_id}.npz")
        if not os.path.exists(f):
            continue
        z = np.load(f)
        ng = z["ng_all"].astype(float); ca = z["can_all"].astype(float)
        if ng.sum() < min_returns or ca.sum() < min_returns:
            continue
        NG.append(ng); CAN.append(ca); keep.append(i)
        ngm = 0.5 * (z["ng_edges"][:-1] + z["ng_edges"][1:])
        cam = 0.5 * (z["can_edges"][:-1] + z["can_edges"][1:])
    m = m.loc[keep].reset_index(drop=True)
    NG = np.array(NG); CAN = np.array(CAN)
    hs = np.concatenate([ngm[ngm > -1.0], cam[cam > 2.0]])
    CNT = np.concatenate([NG[:, ngm > -1.0], CAN[:, cam > 2.0]], axis=1)
    return m, hs, CNT, np.maximum(CNT.sum(1), 1.0)


def weights(hs):
    """The candidate weight functions, in the order they were first tried."""
    W = {}
    for hi in (2.0, 3.0, 5.0):
        W[f"boxcar 0.15-{hi:g}"] = ((hs > 0.15) & (hs <= hi)).astype(float)
    for lam in (0.15, 0.3, 0.6, 1.2):
        W[f"exp decay lam={lam}"] = np.where(hs > 0.15, np.exp(-(hs - 0.15) / lam), 0.0)
    for h0, sg in ((0.25, 0.15), (0.35, 0.25), (0.5, 0.4)):
        W[f"gauss h0={h0} sd={sg}"] = np.where(hs > 0.05, np.exp(-0.5 * ((hs - h0) / sg) ** 2), 0.0)
    W["soft sigmoid edge, k=40"] = 1 / (1 + np.exp(-40 * (hs - 0.15))) * np.where(hs <= 5, 1, 0)
    W["ramp up then decay"] = np.where(
        hs > 0.05, np.clip((hs - 0.05) / 0.2, 0, 1) * np.exp(-np.maximum(hs - 0.25, 0) / 0.8), 0.0)
    W["sqrt-h weight (0.15-3)"] = np.where((hs > 0.15) & (hs <= 3.0), np.sqrt(np.maximum(hs, 0)), 0.0)
    return W


class Scorer:
    def __init__(self, m, seeds=20):
        self.y = m.resid_mm.to_numpy(float)
        self.blk = M._blocks(m.easting.to_numpy(float), m.northing.to_numpy(float), 10.0)
        ub = sorted(set(m.ept_block))
        self.D = np.column_stack([(m.ept_block == b).to_numpy(float) for b in ub[1:]]) \
            if len(ub) > 1 else np.zeros((len(m), 0))
        self.seeds = seeds
        self.n_blocks = len(np.unique(self.blk))

    def cv(self, x, with_block):
        X = np.c_[x, self.D] if with_block and self.D.shape[1] else np.asarray(x).reshape(-1, 1)
        out = []
        for s in range(self.seeds):
            ub = np.unique(self.blk); r = np.random.default_rng(s); r.shuffle(ub)
            out.append(M._cv_r2(X, self.y, self.blk, np.array_split(ub, 5)))
        return float(np.mean(out)), float(np.std(out, ddof=1))


def metric(CNT, TOT, w):
    return (CNT * w).sum(1) / TOT


def _ref(row, obs):
    r = REFERENCE.get(row)
    return "" if r is None else f"   [2026-08-27 run {r:+.4f}, d {obs - r:+.4f}]"


# Candidate rules for --adaptive. These grids are MY proposal, not measured: a
# coarse ladder chosen to bracket the optima the stratified sweep found (4 m pooled,
# 8-12 m vegetated) without being fine enough to fit noise.
EDGES_CONST = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0)
FIRSTPASS = (0.15, 2.0)          # the cheap statistic a first scan of the cell yields
THRESHOLDS = (0.03, 0.06, 0.09, 0.15)
EDGES_LO = (1.0, 2.0, 3.0, 4.0)
EDGES_HI = (4.0, 6.0, 8.0, 12.0)


def _edge_metrics(hs, CNT, TOT, edges):
    """metric for every candidate upper edge, precomputed once."""
    return {e: metric(CNT, TOT, ((hs > 0.15) & (hs <= e)).astype(float)) for e in edges}


def _fit_pred(x_tr, y_tr, x_te):
    D = np.c_[np.ones(len(x_tr)), x_tr]
    b, *_ = np.linalg.lstsq(D, y_tr, rcond=None)
    return np.c_[np.ones(len(x_te)), x_te] @ b


def _inner_select(idx, y, blk, cand, rng):
    """Pick the rule with the best INNER-CV score, using training marks only."""
    ub = np.unique(blk[idx]); rng.shuffle(ub)
    folds = np.array_split(ub, 4)
    best, bestscore = None, -np.inf
    for nm, x in cand.items():
        pred = np.full(len(idx), np.nan)
        for fo in folds:
            te = np.isin(blk[idx], fo); tr = ~te
            if te.sum() < 2 or tr.sum() < 15:
                continue
            pred[te] = _fit_pred(x[idx][tr], y[idx][tr], x[idx][te])
        ok = np.isfinite(pred)
        if ok.sum() < 10:
            continue
        sc = 1 - np.sum((y[idx][ok] - pred[ok]) ** 2) / np.sum((y[idx][ok] - y[idx][ok].mean()) ** 2)
        if sc > bestscore:
            best, bestscore = nm, sc
    return best


def adaptive(m, hs, CNT, TOT, seeds):
    """Nested CV: does a window chosen from a first scan of each cell beat a constant?"""
    y = m.resid_mm.to_numpy(float)
    blk = M._blocks(m.easting.to_numpy(float), m.northing.to_numpy(float), 10.0)
    EM = _edge_metrics(hs, CNT, TOT, sorted(set(EDGES_CONST) | set(EDGES_LO) | set(EDGES_HI)))
    v = metric(CNT, TOT, ((hs > FIRSTPASS[0]) & (hs <= FIRSTPASS[1])).astype(float))

    const = {f"const {e:g} m": EM[e] for e in EDGES_CONST}
    adapt = {}
    for t in THRESHOLDS:
        for lo in EDGES_LO:
            for hi in EDGES_HI:
                if hi <= lo:
                    continue
                adapt[f"v<={t:g} -> {lo:g} m, else {hi:g} m"] = np.where(v <= t, EM[lo], EM[hi])

    print("=" * 84)
    print("IS THE WINDOW ALLOWED TO DEPEND ON A FIRST SCAN OF THE CELL?")
    print("=" * 84)
    print(f"  first-pass statistic v = fraction of returns in ({FIRSTPASS[0]:g}, {FIRSTPASS[1]:g}] m")
    print(f"  -- computed from the cell's OWN returns, so the rule needs no external label.")
    print(f"  NESTED CV: the rule is chosen by an inner 4-fold on the training blocks only,")
    print(f"  then scored on the held-out block. Both families select the same way, so the")
    print(f"  adaptive family is not rewarded for having more knobs.")
    print(f"  {len(const)} constant rules, {len(adapt)} adaptive rules.\n")

    rows = []
    for nm, fam in (("constant only", const), ("adaptive only", adapt),
                    ("both families", {**const, **adapt})):
        scores, picks = [], {}
        for s in range(seeds):
            ub = np.unique(blk); r = np.random.default_rng(s); r.shuffle(ub)
            pred = np.full(len(y), np.nan)
            for fo in np.array_split(ub, 5):
                te = np.isin(blk, fo); tr = np.where(~te)[0]
                if te.sum() < 2 or len(tr) < 15:
                    continue
                pick = _inner_select(tr, y, blk, fam, np.random.default_rng(1000 + s))
                if pick is None:
                    continue
                picks[pick] = picks.get(pick, 0) + 1
                x = fam[pick]
                pred[te] = _fit_pred(x[tr], y[tr], x[te])
            ok = np.isfinite(pred)
            scores.append(1 - np.sum((y[ok] - pred[ok]) ** 2) / np.sum((y[ok] - y[ok].mean()) ** 2))
        rows.append((nm, np.mean(scores), np.std(scores, ddof=1), np.array(scores), picks))

    print(f"  {'family':16s} {'nested CV R2':>13} {'sd':>7}   most-chosen rule (of "
          f"{seeds * 5} fold fits)")
    for nm, mu, sd, _, picks in rows:
        top = sorted(picks.items(), key=lambda kv: -kv[1])[:2]
        t = "; ".join(f"{k} x{n}" for k, n in top)
        print(f"  {nm:16s} {mu:+13.4f} {sd:7.4f}   {t}")

    a = dict((r[0], r[3]) for r in rows)
    d = a["adaptive only"] - a["constant only"]
    print(f"\n  paired over {seeds} seeds:  adaptive - constant = {d.mean():+.4f} "
          f"+- {d.std(ddof=1) / np.sqrt(seeds):.4f} (SE), wins {int((d > 0).sum())}/{seeds}")
    print(f"\n  For reference, the same data with the edge FIXED, no selection at all:")
    for e in EDGES_CONST:
        sc = Scorer(m, seeds)
        print(f"    fixed {e:4g} m   CV R2 {sc.cv(EM[e], False)[0]:+.4f}")


NGV_LO, NGV_HI = 0.15, 4.0

NGV_DEFINITION = """\
NGV -- near-ground vegetation fraction.  Name is a proposal; the DEFINITION is what
must travel with any coefficient fitted from it.

  1. Fit an order-2 least-squares surface to the CLASS-2 returns within 7.5 m of the mark.
     Order 2 removes local slope AND curvature, so neither enters the index.
  2. For every return in that radius -- every class, not just ground -- take its
     SLOPE-NORMAL height above that surface:  h = (z - S(x,y)) / sqrt(1 + gx^2 + gy^2).
  3. NGV = ( number of returns with 0.15 < h <= 4.00 m ) / ( ALL returns in the radius ).

  Range 0 to 1.  Denominator is EVERY return, so NGV falls if canopy above 4 m grows
  while understory stays fixed -- it is a fraction, not a density.

  The window is (0.15, 4.0] m and both edges are chosen, not arbitrary:
    lower 0.15 m -- the ground peak's own tail leaks in below this; 0.05 m costs 0.074 CV R2.
    upper 4.0  m -- the top of a plateau flat from 3 to 8 m (within 0.005), so the exact
                    value is not delicate. See --edges and --adaptive.

  Because the surface is fitted from the box's OWN returns, NGV is invariant to any
  vertical shift of the cloud. It structurally cannot carry offset information, which is
  what makes the regression below a real test rather than a tautology."""


def index_regression(m, hs, CNT, TOT, block_km, n_boot, bin_width=0.06):
    from scipy import stats as st
    y = m.resid_mm.to_numpy(float)
    ngv = metric(CNT, TOT, ((hs > NGV_LO) & (hs <= NGV_HI)).astype(float))
    old = metric(CNT, TOT, ((hs > 0.15) & (hs <= 2.0)).astype(float))

    print("=" * 84)
    print("THE VEGETATION INDEX")
    print("=" * 84)
    print(NGV_DEFINITION)
    print(f"\n  n = {len(y)} marks.  NGV percentiles:")
    print("   ", "  ".join(f"p{q}={np.percentile(ngv, q):.3f}" for q in (0, 10, 25, 50, 75, 90, 99, 100)))
    print(f"    NGV vs the (0.15, 2.0] incumbent: Pearson {st.pearsonr(ngv, old)[0]:+.4f}, "
          f"NGV/old ratio median {np.median(ngv / np.maximum(old, 1e-9)):.3f}")
    print(f"  offset = USGS surveyed_Z - delivered_LAZ_Z, mm, +ve = the surface reads LOW")
    print(f"    sd {y.std(ddof=1):.1f} mm, median {np.median(y):+.1f}\n")

    E = np.arange(0, ngv.max() + bin_width, bin_width)
    rows = []
    for lo, hi in zip(E[:-1], E[1:]):
        k = (ngv >= lo) & (ngv < hi)
        se = y[k].std(ddof=1) / np.sqrt(k.sum()) if k.sum() > 1 else np.nan
        rows.append((lo, hi, int(k.sum()),
                     np.median(y[k]) if k.sum() else np.nan,
                     y[k].mean() if k.sum() else np.nan, se))
    b = pd.DataFrame(rows, columns=["lo", "hi", "n", "median", "mean", "se"])
    print(f"  EVERY bin at its true span, counts visible -- no bin dropped, none merged.")
    print(f"  {'bin':>13} {'n':>4} {'median':>8} {'mean':>8} {'SE':>7}")
    for _, r in b.iterrows():
        print(f"  {r.lo:.2f}-{r.hi:.2f} {int(r.n):4d} " +
              ("      -        -       -" if r.n == 0 else
               f"{r['median']:8.1f} {r['mean']:8.1f} {r.se:7.1f}"))

    g = b[(b.n > 1) & np.isfinite(b.se)].copy(); g["x"] = 0.5 * (g.lo + g.hi)
    w = 1 / g.se.values ** 2
    print(f"\n  REGRESSION  offset = a + b * NGV      ({len(g)} bins with n>1 of {len(b)})")
    bd, sd_ = M.wls(g.x.values, g["mean"].values, w)
    ba, sa = M.wls(g.x.values, g["mean"].values, g.n.values.astype(float))
    print(f"    binned, DESIGN-weighted 1/SE^2 : a {bd[0]:+7.1f} +/- {sd_[0]:.1f}   "
          f"b {bd[1]:+8.1f} +/- {sd_[1]:.1f} mm per unit NGV")
    print(f"    binned, ABUNDANCE-weighted by n: a {ba[0]:+7.1f} +/- {sa[0]:.1f}   "
          f"b {ba[1]:+8.1f} +/- {sa[1]:.1f} mm per unit NGV")
    lr = st.linregress(ngv, y)
    print(f"    per-mark, unweighted           : a {lr.intercept:+7.1f} +/- {lr.intercept_stderr:.1f}   "
          f"b {lr.slope:+8.1f} +/- {lr.stderr:.1f}   (p {lr.pvalue:.1e})")

    bo_b, so_b = M.fit_origin(g.x.values, g["mean"].values, w)
    bo_m, so_m = M.fit_origin(ngv, y)
    print(f"\n  THROUGH THE ORIGIN  offset = b * NGV   -- the form to carry back to the DEM")
    print(f"    binned, 1/SE^2 weighted        : b = {bo_b:+8.1f} +/- {so_b:.1f} mm per unit NGV")
    print(f"    per-mark, unweighted           : b = {bo_m:+8.1f} +/- {so_m:.1f} mm per unit NGV")
    print(f"    the origin is a CHECK, not an assumption: the free intercept above is "
          f"{bd[0]:+.1f} +/- {sd_[0]:.1f} mm.")

    B = block_km * 1000.0
    blk = M._blocks(m.easting.to_numpy(float), m.northing.to_numpy(float), block_km)
    ub = np.unique(blk); rng = np.random.default_rng(0); sl = []
    for _ in range(n_boot):
        pick = rng.choice(ub, size=len(ub), replace=True)
        idx = np.concatenate([np.where(blk == k)[0] for k in pick])
        if len(np.unique(ngv[idx])) > 2:
            sl.append(np.polyfit(ngv[idx], y[idx], 1)[0])
    sl = np.array(sl)
    print(f"\n  UNCERTAINTY that respects the spatial clustering")
    print(f"    naive per-mark SE              : {lr.stderr:.1f}")
    print(f"    block bootstrap, {len(ub)} blocks of {block_km:.0f} km, {len(sl)} draws: "
          f"b = {sl.mean():+.1f} +/- {sl.std(ddof=1):.1f}  "
          f"-> SE inflated {sl.std(ddof=1)/lr.stderr:.2f}x")

    sc = Scorer(m, 20)
    a_ngv = sc.cv(ngv, False); a_old = sc.cv(old, False)
    print(f"\n  HELD-OUT (10 km blocks, 5 folds, 20 seeds)")
    print(f"    NGV  (0.15, 4.0]               : CV R2 {a_ngv[0]:+.4f} +/- {a_ngv[1]:.4f}   "
          f"RMSE {y.std(ddof=1)*np.sqrt(1-a_ngv[0]):.1f} mm")
    print(f"    incumbent (0.15, 2.0]          : CV R2 {a_old[0]:+.4f} +/- {a_old[1]:.4f}   "
          f"RMSE {y.std(ddof=1)*np.sqrt(1-a_old[0]):.1f} mm")
    print(f"    no model                       : {'':22s} RMSE {y.std(ddof=1):.1f} mm")
    return ngv


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shape", action="store_true", help="bake-off of weight functions w(h)")
    ap.add_argument("--edges", action="store_true", help="sweep the upper and lower edges")
    ap.add_argument("--strata", action="store_true",
                    help="repeat the edge sweep within NVA/VVA and lowveg terciles")
    ap.add_argument("--index", action="store_true",
                    help="define the vegetation index on the chosen window and regress it")
    ap.add_argument("--block-km", type=float, default=10.0)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--adaptive", action="store_true",
                    help="does a window chosen per-cell beat a constant one?")
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.all:
        a.shape = a.edges = a.strata = a.adaptive = a.index = True
    if not (a.shape or a.edges or a.strata or a.adaptive or a.index):
        ap.error("pick a section, or --all")

    m, hs, CNT, TOT = build()
    sc = Scorer(m, a.seeds)
    print(f"n = {len(m)} marks   {sc.n_blocks} spatial blocks of 10 km   "
          f"{m.ept_block.nunique()} EPT blocks   offset sd {sc.y.std(ddof=1):.1f} mm")
    print(f"height axis: {hs.min():+.3f} to {hs.max():+.2f} m, {len(hs)} bins "
          f"({int((hs <= 2).sum())} at 20 mm, {int((hs > 2).sum())} at 250 mm)")
    print(f"CV: held-out 10 km blocks, 5 folds, mean +/- sd over {a.seeds} fold seeds\n")

    if a.shape:
        print("=" * 84)
        print("WEIGHT FUNCTION w(h), ranked by the +block score")
        print("=" * 84)
        print(f"  {'w(h)':30s} {'alone':>9} {'sd':>7} {'+block':>9} {'sd':>7}")
        rows = []
        for nm, w in weights(hs).items():
            x = metric(CNT, TOT, w)
            rows.append((sc.cv(x, True), sc.cv(x, False), nm))
        for (b, sb), (al, sa), nm in sorted(rows, reverse=True):
            print(f"  {nm:30s} {al:+9.4f} {sa:7.4f} {b:+9.4f} {sb:7.4f}"
                  + _ref(("shape", nm), b))
        print()

    if a.edges:
        print("=" * 84)
        print("UPPER EDGE, lower fixed at 0.15 m")
        print("=" * 84)
        print(f"  {'upper (m)':>9} {'boxcar alone':>14} {'boxcar +block':>15} {'sqrt-h +block':>15}")
        for hi in (2, 3, 4, 5, 6, 8, 10, 14, 20):
            w = ((hs > 0.15) & (hs <= hi)).astype(float)
            w2 = np.where((hs > 0.15) & (hs <= hi), np.sqrt(np.maximum(hs, 0)), 0.0)
            al, _ = sc.cv(metric(CNT, TOT, w), False)
            b, _ = sc.cv(metric(CNT, TOT, w), True)
            b2, _ = sc.cv(metric(CNT, TOT, w2), True)
            print(f"  {hi:9d} {al:+14.4f} {b:+15.4f} {b2:+15.4f}" + _ref(("upper", hi), b))
        print(f"\n  LOWER EDGE, upper fixed at 5 m (boxcar)")
        for lo in (0.05, 0.10, 0.13, 0.15, 0.18, 0.22, 0.30):
            w = ((hs > lo) & (hs <= 5.0)).astype(float)
            al, _ = sc.cv(metric(CNT, TOT, w), False)
            b, _ = sc.cv(metric(CNT, TOT, w), True)
            print(f"  lower {lo:.2f} m {al:+12.4f} {b:+15.4f}" + _ref(("lower", lo), b))
        print()

    if a.strata:
        print("=" * 84)
        print("WITHIN STRATA -- does the best window differ under vegetation?")
        print("=" * 84)
        print("  NVA/VVA is the SURVEY's own non-vegetated / vegetated designation, not a cut")
        print("  of mine. The lowveg terciles ARE my cut, chosen only to give three equal groups.")
        print("  Scored ALONE (a stratum can hold too few EPT blocks for block intercepts).\n")
        q = m.lowveg.quantile([1 / 3, 2 / 3]).to_numpy()
        strata = [("all marks", np.ones(len(m), bool)),
                  ("NVA (open)", (m.point_type_g2 == "NVA").to_numpy()),
                  ("VVA (vegetated)", (m.point_type_g2 == "VVA").to_numpy()),
                  (f"lowveg <= {q[0]:.3f}", (m.lowveg <= q[0]).to_numpy()),
                  (f"lowveg {q[0]:.3f}-{q[1]:.3f}", ((m.lowveg > q[0]) & (m.lowveg <= q[1])).to_numpy()),
                  (f"lowveg > {q[1]:.3f}", (m.lowveg > q[1]).to_numpy())]
        uppers = (1, 2, 4, 8, 12, 16, 24, 40)
        print(f"  {'stratum':22s} {'n':>4} {'blk':>4} {'sd_mm':>7} " +
              " ".join(f"{u:>7d}" for u in uppers) + "   best")
        for nm, k in strata:
            if k.sum() < 30:
                print(f"  {nm:22s} {k.sum():4d}   -- too few marks to cross-validate")
                continue
            s = Scorer(m[k].reset_index(drop=True), a.seeds)
            vals = []
            for hi in uppers:
                w = ((hs > 0.15) & (hs <= hi)).astype(float)
                vals.append(s.cv(metric(CNT[k], TOT[k], w), False)[0])
            best = uppers[int(np.argmax(vals))]
            print(f"  {nm:22s} {k.sum():4d} {s.n_blocks:4d} {s.y.std(ddof=1):7.1f} " +
                  " ".join(f"{v:+7.3f}" for v in vals) + f"   {best} m")
        print("\n  LIMIT OF THIS CONTROL SET: it does not sample full forest. Canopy cover")
        print("  (PyForestScan, r=10 m) is 0.000 at the median and reaches 0.30 at only "
              f"{int((m['cover_r10'] >= 0.30).sum())} marks;")
        print(f"  lowveg reaches 0.30 at {int((m.lowveg >= 0.30).sum())} marks, 0.40 at "
              f"{int((m.lowveg >= 0.40).sum())}. Marks are sited for sky view.")
        print("  Any window chosen here is calibrated on understory, NOT on closed canopy.")
        print()

    if a.index:
        index_regression(m, hs, CNT, TOT, a.block_km, a.n_boot)
        print()

    if a.adaptive:
        adaptive(m, hs, CNT, TOT, a.seeds)


if __name__ == "__main__":
    main()
