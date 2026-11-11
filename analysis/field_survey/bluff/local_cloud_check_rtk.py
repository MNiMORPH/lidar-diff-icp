"""Bluff RTK points against the gen2 and gen1 POINT CLOUDS directly (no 5 m grid, no
slope-normal median, no bilinear sampling).

For each FIX hillslope RTK point and radius r (3 m proposed, 5 m also run): fit a plane
z = a + b*(x-x0) + c*(y-y0) by least squares to the cloud points within r of the tied
position, evaluate it at the RTK x,y, and report residual = navd88_tied - plane.

Populations:
  gen2       2021 class-2 ground cloud, as delivered.
  gen1_raw   2008 last returns (return_number == number_of_returns), RAW positions.
  gen1_reg   same returns at REGISTERED positions: raw + per-swath star shift + lateral shift
             (what the pipeline applied under swath_frame="star"; see LOCAL_CLOUD_CHECK.md).
  gen1_is    VARIANT, NOT APPLIED BY THE PRODUCT: raw + internal alignment + star + lateral,
             the composition named in the task brief.
  csf_reg    2008 CSF-selected ground points (what the pipeline gridded), registered.

Every number in LOCAL_CLOUD_CHECK.md is pasted from results.md written by this script.
"""
import json
import os
import sys
import numpy as np
import pandas as pd
import laspy
from scipy.spatial import cKDTree

ROOT = "/home/awickert/projects/lidar-diff-icp"
RTK = f"{ROOT}/data/field/2026-10-25/rtk_points_all.csv"
GEN2 = f"{ROOT}/data/after/valley/4342-29-64_3dep_class2.laz"
GEN1 = f"{ROOT}/data/before/4342-29-64.laz"
CSF = f"{ROOT}/data/csf_cache/4342-29-64.las"
CORR = f"{ROOT}/data/derived/4342-29-64/corrections.json"
OUT = f"{ROOT}/data/derived/4342-29-64/bluff_check/local"
os.makedirs(OUT, exist_ok=True)

RADII = (3.0, 5.0)          # 3 m = Andy's proposal; 5 m also requested
MIN_PTS = 3                 # given in the brief for CSF; used for every plane fit (3 points
                            # define a plane; fewer cannot) -- PROPOSAL where not CSF
SHIFTS = (-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5)
PAD = 14.0                  # bbox pad around the RTK points for the cloud reads (m)

md = []
def P(s=""):
    print(s, flush=True); md.append(s)

def nmad(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return 1.4826 * np.median(np.abs(v - np.median(v))) if v.size else np.nan

def plane_fit(x, y, z, x0, y0):
    """LS plane through points; returns (z at x0,y0, slope_deg, downslope_aspect_deg, rms_m, n)."""
    n = len(x)
    if n < MIN_PTS:
        return np.nan, np.nan, np.nan, np.nan, n
    A = np.column_stack([np.ones(n), x - x0, y - y0])
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)
    a, b, c = coef
    res = z - A @ coef
    slope = np.degrees(np.arctan(np.hypot(b, c)))
    aspect = np.degrees(np.arctan2(-b, -c)) % 360.0   # compass azimuth of the DOWNSLOPE direction
    return a, slope, aspect, float(np.sqrt(np.mean(res**2))), n

# ---------------------------------------------------------------- RTK
d = pd.read_csv(RTK)
h = d[d.status.astype(str).str.startswith("FIX") & d.group.astype(str).str.startswith("hillslope")].copy()
h = h.reset_index(drop=True)
P(f"RTK rows: {len(h)} (status FIX*, group hillslope*)")
E = h.easting_tied.to_numpy(float); N = h.northing_tied.to_numpy(float); Z = h.navd88_tied_m.to_numpy(float)
bbox = (E.min() - PAD, E.max() + PAD, N.min() - PAD, N.max() + PAD)
P(f"RTK bbox E {E.min():.1f}-{E.max():.1f}, N {N.min():.1f}-{N.max():.1f}; read pad {PAD} m")

corr = json.load(open(CORR))
star = {int(k): np.array(v, float) for k, v in corr["per_swath_star_frame_dxdydz_m"].items()}
internal = {int(k): np.array(v, float) for k, v in corr["per_swath_internal_alignment_dxdydz_m"].items()}
lat = np.array(corr["cross_epoch_datum"]["horizontal_shift_m"], float)
P(f"swath_frame={corr['swath_frame']}; star shifts (ADDED to x,y,z) {corr['per_swath_star_frame_dxdydz_m']}")
P(f"internal alignment (NOT applied under star) {corr['per_swath_internal_alignment_dxdydz_m']}")
P(f"lateral shift (ADDED to x,y) {lat.tolist()}")

def in_box(x, y):
    return (x >= bbox[0]) & (x <= bbox[1]) & (y >= bbox[2]) & (y <= bbox[3])

# ---------------------------------------------------------------- gen2 (chunked, bbox only)
g2 = []
with laspy.open(GEN2) as r:
    for ch in r.chunk_iterator(4_000_000):
        x = np.asarray(ch.x, float); y = np.asarray(ch.y, float)
        m = in_box(x, y)
        if m.any():
            g2.append(np.column_stack([x[m], y[m], np.asarray(ch.z, float)[m],
                                       np.asarray(ch.classification)[m]]))
g2 = np.vstack(g2)
P(f"gen2 points in bbox: {len(g2):,}; classification values {np.unique(g2[:,3]).astype(int).tolist()}")
g2 = g2[g2[:, 3] == 2][:, :3]
P(f"gen2 class-2 points in bbox: {len(g2):,}")

# ---------------------------------------------------------------- gen1 raw (whole tile, then bbox)
with laspy.open(GEN1) as r:
    L = r.read()
x = np.asarray(L.x, float); y = np.asarray(L.y, float); z = np.asarray(L.z, float)
rn = np.asarray(L.return_number); nr = np.asarray(L.number_of_returns)
ps = np.asarray(L.point_source_id); cls = np.asarray(L.classification)
m = in_box(x, y)
P(f"gen1 raw returns in bbox: {m.sum():,} of {len(x):,}; point_source_ids {np.unique(ps[m]).tolist()}; "
  f"classification values {np.unique(cls[m]).tolist()}")
x, y, z, rn, nr, ps = x[m], y[m], z[m], rn[m], nr[m], ps[m]
XI = np.asarray(L.X)[m]; YI = np.asarray(L.Y)[m]; ZI = np.asarray(L.Z)[m]   # integer coords for CSF matching
del L
last = rn == nr
P(f"gen1 last returns in bbox: {last.sum():,} ({100*last.mean():.1f}%); single returns "
  f"(nr==1): {100*(nr==1).mean():.1f}% of all returns")
g1 = np.column_stack([x, y, z])[last]; ps1 = ps[last]
XI1, YI1, ZI1 = XI[last], YI[last], ZI[last]
def register(xyz, psid, with_internal=False):
    out = xyz.copy()
    for s in np.unique(psid):
        mm = psid == s
        sh = star[int(s)] + (internal[int(s)] if with_internal else 0.0)
        out[mm] += sh
    out[:, 0] += lat[0]; out[:, 1] += lat[1]
    return out
g1r = register(g1, ps1)
g1i = register(g1, ps1, with_internal=True)

# ---------------------------------------------------------------- CSF (chunked)
cs = []
with laspy.open(CSF) as r:
    for ch in r.chunk_iterator(4_000_000):
        cx_ = np.asarray(ch.x, float); cy_ = np.asarray(ch.y, float)
        mm = in_box(cx_, cy_)
        if mm.any():
            cs.append(np.column_stack([cx_[mm], cy_[mm], np.asarray(ch.z, float)[mm],
                                       np.asarray(ch.point_source_id)[mm],
                                       np.asarray(ch.X)[mm], np.asarray(ch.Y)[mm], np.asarray(ch.Z)[mm]]))
cs = np.vstack(cs)
P(f"gen1 CSF ground points in bbox: {len(cs):,}; point_source_ids {np.unique(cs[:,3]).astype(int).tolist()}")
csf = cs[:, :3]; pscs = cs[:, 3].astype(int)
csf_r = register(csf, pscs)
# which raw last returns are in the CSF selection (match on integer coordinates)
csf_keys = set(map(tuple, cs[:, 4:7].astype(np.int64)))
raw_keys = list(map(tuple, np.column_stack([XI1, YI1, ZI1]).astype(np.int64)))
in_csf = np.array([k in csf_keys for k in raw_keys])
P(f"CSF membership of gen1 last returns in bbox: {in_csf.sum():,} of {len(in_csf):,} "
  f"({100*in_csf.mean():.1f}%) matched by integer X,Y,Z")

clouds = {"gen2": g2, "gen1_raw": g1, "gen1_reg": g1r, "gen1_is": g1i, "csf_reg": csf_r}
trees = {k: cKDTree(v[:, :2]) for k, v in clouds.items()}

# ---------------------------------------------------------------- step 2: per-point planes
rows = []
planes_r3 = {}
for r in RADII:
    for i in range(len(h)):
        rec = dict(point=h.point[i], group=h.group[i], radius_m=r, E=E[i], N=N[i], navd88_tied_m=Z[i],
                   slope_grid_deg=h.slope_deg[i], tied_minus_gen2_grid_m=h.tied_minus_gen2_m[i],
                   tied_minus_gen1_grid_m=h.tied_minus_gen1_m[i])
        for k, c in clouds.items():
            idx = trees[k].query_ball_point([E[i], N[i]], r)
            sub = c[idx]
            a, sl, asp, rms, n = plane_fit(sub[:, 0], sub[:, 1], sub[:, 2], E[i], N[i])
            rec[f"{k}_n"] = n; rec[f"{k}_plane_z_m"] = a; rec[f"{k}_resid_m"] = Z[i] - a
            rec[f"{k}_slope_deg"] = sl; rec[f"{k}_aspect_deg"] = asp; rec[f"{k}_rms_m"] = rms
            if k == "gen2":
                g2coef = None
                if n >= MIN_PTS:
                    A_ = np.column_stack([np.ones(n), sub[:, 0] - E[i], sub[:, 1] - N[i]])
                    g2coef, *_ = np.linalg.lstsq(A_, sub[:, 2], rcond=None)
                if r == RADII[0]:
                    planes_r3[i] = (a, sl, asp, n, sub)
            if k in ("gen1_reg", "csf_reg") and g2coef is not None and n > 0:
                # PROPOSAL (not in the brief): height of each gen1 return above the gen2 plane
                # at the return's own x,y; the median is a one-parameter offset that borrows the
                # 2021 slope and cannot extrapolate the way a 3-point plane does.
                hh = sub[:, 2] - (g2coef[0] + g2coef[1] * (sub[:, 0] - E[i]) + g2coef[2] * (sub[:, 1] - N[i]))
                rec[f"{k}_h_above_gen2_median_m"] = float(np.median(hh))
                rec[f"{k}_h_above_gen2_p25_m"] = float(np.percentile(hh, 25))
                rec[f"{k}_h_above_gen2_min_m"] = float(hh.min())
                rec[f"{k}_resid_gen2slope_m"] = Z[i] - (g2coef[0] + float(np.median(hh)))
            elif k in ("gen1_reg", "csf_reg"):
                for c_ in ("h_above_gen2_median_m", "h_above_gen2_p25_m", "h_above_gen2_min_m", "resid_gen2slope_m"):
                    rec[f"{k}_{c_}"] = np.nan
        rows.append(rec)
tab = pd.DataFrame(rows)
tab.to_csv(f"{OUT}/per_point_planes.csv", index=False)
P(f"\nWrote {OUT}/per_point_planes.csv ({len(tab)} rows = 39 points x {len(RADII)} radii)")

def stat_line(name, v, n=None):
    v = np.asarray(v, float); ok = np.isfinite(v)
    if ok.sum() == 0:
        return f"| {name} | 0 | - | - | - | - |"
    return (f"| {name} | {ok.sum()} | {1000*np.median(v[ok]):.1f} | {1000*nmad(v[ok]):.1f} | "
            f"{1000*np.mean(v[ok]):.1f} | {1000*np.std(v[ok], ddof=1) if ok.sum()>1 else np.nan:.1f} |")

P("\n## Step 2. Residual = navd88_tied - local LS plane, mm (median, NMAD, mean, SD; n = points with a fit)")
for r in RADII:
    t = tab[tab.radius_m == r]
    P(f"\n### r = {r:.0f} m")
    P("| population / group | n | median | NMAD | mean | SD |")
    P("|---|---|---|---|---|---|")
    for k in clouds:
        P(stat_line(f"{k}, ALL", t[f"{k}_resid_m"]))
        for g in ("hillslope_RS3", "hillslope_T4a"):
            P(stat_line(f"{k}, {g}", t[t.group == g][f"{k}_resid_m"]))
    P(f"\nPoint counts within r={r:.0f} m (median [min-max] over the 39 points): " +
      "; ".join(f"{k} {int(np.median(t[f'{k}_n']))} [{int(t[f'{k}_n'].min())}-{int(t[f'{k}_n'].max())}]"
                for k in clouds))
    P(f"Plane-fit RMS about the plane, m (median over points): " +
      "; ".join(f"{k} {np.nanmedian(t[f'{k}_rms_m']):.3f}" for k in clouds))
    P(f"gen2 plane slope deg (median [min-max]): {np.nanmedian(t.gen2_slope_deg):.1f} "
      f"[{np.nanmin(t.gen2_slope_deg):.1f}-{np.nanmax(t.gen2_slope_deg):.1f}]; "
      f"downslope aspect deg (median): {np.nanmedian(t.gen2_aspect_deg):.0f}; "
      f"grid slope_deg column (median): {np.nanmedian(t.slope_grid_deg):.1f}")
    P(f"gen1_reg plane slope deg (median [min-max]): {np.nanmedian(t.gen1_reg_slope_deg):.1f} "
      f"[{np.nanmin(t.gen1_reg_slope_deg):.1f}-{np.nanmax(t.gen1_reg_slope_deg):.1f}]; "
      f"downslope aspect deg (median): {np.nanmedian(t.gen1_reg_aspect_deg):.0f}")
    dg = t.gen2_resid_m - t.tied_minus_gen2_grid_m
    P(f"gen2 cloud residual MINUS grid-based tied_minus_gen2_m: median {1000*np.nanmedian(dg):.1f} mm, "
      f"NMAD {1000*nmad(dg):.1f} mm, n {np.isfinite(dg).sum()}; grid column itself: median "
      f"{1000*np.nanmedian(t.tied_minus_gen2_grid_m):.1f} mm, NMAD {1000*nmad(t.tied_minus_gen2_grid_m):.1f} mm")
    ok = np.isfinite(t.tied_minus_gen1_grid_m)
    if ok.any():
        P(f"Points where the gen1 GRID has a value (n={ok.sum()}): grid tied_minus_gen1 median "
          f"{1000*np.nanmedian(t.tied_minus_gen1_grid_m[ok]):.1f} mm; cloud csf_reg residual at those points "
          f"median {1000*np.nanmedian(t.csf_reg_resid_m[ok]):.1f} mm (n fit {np.isfinite(t.csf_reg_resid_m[ok]).sum()}); "
          f"gen1_reg residual median {1000*np.nanmedian(t.gen1_reg_resid_m[ok]):.1f} mm")
    # gen1 relative to gen2, same points: the cloud-level DoD at the bluff
    for k in ("gen1_raw", "gen1_reg", "gen1_is", "csf_reg"):
        dd = t[f"{k}_plane_z_m"] - t.gen2_plane_z_m
        P(f"{k} plane minus gen2 plane at the RTK x,y: median {1000*np.nanmedian(dd):.1f} mm, "
          f"NMAD {1000*nmad(dd):.1f} mm, n {np.isfinite(dd).sum()}")
    P("PROPOSAL (gen2-slope offset): per-point MEDIAN height of gen1 returns above the gen2 plane, then "
      "residual = navd88_tied - (gen2 plane + that median)")
    for k in ("gen1_reg", "csf_reg"):
        hm = t[f"{k}_h_above_gen2_median_m"]; rr = t[f"{k}_resid_gen2slope_m"]
        P(f"  {k}: median-height-above-gen2 over points: median {1000*np.nanmedian(hm):.1f} mm, NMAD {1000*nmad(hm):.1f} mm, "
          f"n {np.isfinite(hm).sum()}; p25-height median over points {1000*np.nanmedian(t[f'{k}_h_above_gen2_p25_m']):.1f} mm; "
          f"residual RTK minus (gen2 plane + median height): median {1000*np.nanmedian(rr):.1f} mm, NMAD {1000*nmad(rr):.1f} mm")

# per-point table for the report, r=3
t3 = tab[tab.radius_m == RADII[0]]
P(f"\n### Per-point table, r = {RADII[0]:.0f} m (mm; n = points in the fit)")
P("| group | pt | gen2 n | gen2 resid | grid resid | gen2 slope | gen1_reg n | gen1_reg plane resid | gen1_reg-gen2 plane | gen1_reg h>gen2 median | gen1_reg h>gen2 p25 | csf n | csf plane resid | csf h>gen2 median |")
P("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for _, q in t3.iterrows():
    f = lambda v: "-" if not np.isfinite(v) else f"{1000*v:.0f}"
    P(f"| {q.group.replace('hillslope_','')} | {q.point} | {q.gen2_n} | {f(q.gen2_resid_m)} | {f(q.tied_minus_gen2_grid_m)} | "
      f"{q.gen2_slope_deg:.1f} | {q.gen1_reg_n} | {f(q.gen1_reg_resid_m)} | {f(q.gen1_reg_plane_z_m - q.gen2_plane_z_m)} | "
      f"{f(q.gen1_reg_h_above_gen2_median_m)} | {f(q.gen1_reg_h_above_gen2_p25_m)} | {q.csf_reg_n} | {f(q.csf_reg_resid_m)} | "
      f"{f(q.csf_reg_h_above_gen2_median_m)} |")

# ---------------------------------------------------------------- step 3: 2008 return structure
P("\n## Step 3. 2008 last returns within 5 m of any RTK point: height above the LOCAL gen2 plane (r=3 m fit at the nearest RTK point)")
rtk_tree = cKDTree(np.column_stack([E, N]))
def height_above_gen2(xyz):
    dist, j = rtk_tree.query(xyz[:, :2])
    keep = dist <= 5.0
    out = np.full(len(xyz), np.nan)
    for i in np.unique(j[keep]):
        a, sl, asp, n, sub = planes_r3[i]
        if not np.isfinite(a):
            continue
        # re-fit coefficients to evaluate the plane away from the RTK x,y
        A = np.column_stack([np.ones(len(sub)), sub[:, 0] - E[i], sub[:, 1] - N[i]])
        coef, *_ = np.linalg.lstsq(A, sub[:, 2], rcond=None)
        mm = keep & (j == i)
        out[mm] = xyz[mm, 2] - (coef[0] + coef[1] * (xyz[mm, 0] - E[i]) + coef[2] * (xyz[mm, 1] - N[i]))
    return out
def pct_line(name, v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if v.size == 0:
        return f"| {name} | 0 | - | - | - | - | - |"
    p = np.percentile(v, [5, 25, 50, 75, 95])
    return f"| {name} | {v.size} | " + " | ".join(f"{1000*q:.0f}" for q in p) + " |"
P("| population | n | p5 | p25 | p50 | p75 | p95 |  (mm above gen2 plane)")
P("|---|---|---|---|---|---|---|")
hr = height_above_gen2(g1r); hraw = height_above_gen2(g1); hi = height_above_gen2(g1i)
hc = height_above_gen2(csf_r)
P(pct_line("gen1 last returns, REGISTERED positions", hr))
P(pct_line("  ...of which CSF-selected (matched)", hr[in_csf]))
P(pct_line("  ...NOT CSF-selected", hr[~in_csf]))
P(pct_line("gen1 last returns, RAW positions", hraw))
P(pct_line("gen1 last returns, internal+star VARIANT (not applied)", hi))
P(pct_line("CSF ground cloud, registered (direct)", hc))
# gen2 own scatter about its planes, for scale
hg2 = height_above_gen2(g2)
P(pct_line("gen2 class-2 points about the same planes (scale)", hg2))
okr = np.isfinite(hr)
P(f"Fractions of registered gen1 last returns (n={okr.sum()}): below gen2 plane by >0.1 m "
  f"{100*np.mean(hr[okr] < -0.1):.1f}%, within +/-0.1 m {100*np.mean(np.abs(hr[okr]) <= 0.1):.1f}%, "
  f"above by >0.1 m {100*np.mean(hr[okr] > 0.1):.1f}%, above by >1 m {100*np.mean(hr[okr] > 1.0):.1f}%")
P(f"CSF-selected share of those returns: {100*in_csf[okr].mean():.1f}%; "
  f"CSF-selected median {1000*np.nanmedian(hr[okr & in_csf]):.0f} mm vs not-selected median "
  f"{1000*np.nanmedian(hr[okr & ~in_csf]):.0f} mm")
# return-structure of those 2008 returns
near = rtk_tree.query(np.column_stack([x, y]))[0] <= 5.0   # all returns, raw positions
P(f"All 2008 returns within 5 m (raw positions): n {near.sum():,}; number_of_returns histogram "
  f"{dict(zip(*[a.tolist() for a in np.unique(nr[near], return_counts=True)]))}; "
  f"last-return share {100*last[near].mean():.1f}%; "
  f"last returns per RTK 5 m disk: median {np.median([len(trees['gen1_reg'].query_ball_point([E[i],N[i]],5.0)) for i in range(len(h))]):.0f}")
pd.DataFrame(dict(x_reg=g1r[:, 0], y_reg=g1r[:, 1], z=g1r[:, 2], psid=ps1, in_csf=in_csf,
                  h_above_gen2_reg_m=hr, h_above_gen2_raw_m=hraw)).dropna(subset=["h_above_gen2_reg_m"]) \
    .to_csv(f"{OUT}/gen1_last_returns_near_rtk.csv", index=False)
P(f"Wrote {OUT}/gen1_last_returns_near_rtk.csv")

# ---------------------------------------------------------------- step 4: degeneracy
P("\n## Step 4. gen2 cloud residual (r=3 m) under horizontal shifts of the RTK positions, mm")
P("| shift | dE: median | dE: NMAD | dN: median | dN: NMAD |")
P("|---|---|---|---|---|")
shift_rows = []
for s in SHIFTS:
    out = {}
    for ax in ("dE", "dN"):
        res = []
        for i in range(len(h)):
            e = E[i] + (s if ax == "dE" else 0.0); nn = N[i] + (s if ax == "dN" else 0.0)
            idx = trees["gen2"].query_ball_point([e, nn], RADII[0]); sub = g2[idx]
            a, *_ = plane_fit(sub[:, 0], sub[:, 1], sub[:, 2], e, nn)
            res.append(Z[i] - a)
        res = np.array(res); out[ax] = (np.nanmedian(res), nmad(res))
    shift_rows.append(dict(shift_m=s, dE_median_m=out["dE"][0], dE_nmad_m=out["dE"][1],
                           dN_median_m=out["dN"][0], dN_nmad_m=out["dN"][1]))
    P(f"| {s:+.1f} | {1000*out['dE'][0]:.1f} | {1000*out['dE'][1]:.1f} | {1000*out['dN'][0]:.1f} | {1000*out['dN'][1]:.1f} |")
sr = pd.DataFrame(shift_rows); sr.to_csv(f"{OUT}/shift_sensitivity.csv", index=False)
for ax in ("dE", "dN"):
    k = np.polyfit(sr.shift_m, sr[f"{ax}_median_m"], 1)[0]
    P(f"Sensitivity of the median residual to {ax}: {1000*k:.0f} mm per m (LS slope over the 9 shifts)")
sl = np.nanmedian(t3.gen2_slope_deg); asp = np.nanmedian(t3.gen2_aspect_deg)
P(f"Expected from the median gen2 plane (slope {sl:.1f} deg, downslope aspect {asp:.0f} deg): "
  f"d(resid)/dE = +tan(slope)*sin(aspect) = {1000*np.tan(np.radians(sl))*np.sin(np.radians(asp)):.0f} mm/m, "
  f"d(resid)/dN = +tan(slope)*cos(aspect) = {1000*np.tan(np.radians(sl))*np.cos(np.radians(asp)):.0f} mm/m "
  f"(moving the RTK point downslope lowers the plane under it, so the residual rises)")

open(f"{OUT}/results.md", "w").write("\n".join(md) + "\n")
print(f"\nWrote {OUT}/results.md")
