#!/usr/bin/env python3
"""Rebuild the DEM of difference with the canopy-cover correction applied to gen2.

gen2's ground is taken at the cover-dependent percentile of its own near-ground return
column instead of at the median:

    q2(c) = a + b * c               a and b are READ, never assumed (see --relation)

gen1's ground is its median, with the four registration terms applied. Both are per-cell
quantiles of the slope-normal residual to the same plane, so

    DoD = (h2(q2) - h1(0.50)) * |n|          positive = elevation ROSE 2008 -> 2021

Requires one streaming pass over the gen2 cloud to build the class-2 near-ground histogram
for EVERY grid cell (the existing class-split cube covers only the divide reference cells).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \\
        analysis/ridgelines/dod_cover_corrected.py
"""
import argparse, json, os
import numpy as np, laspy, pyarrow.parquet as pq
from lidar_diff_icp import registration as reg
from scipy.ndimage import distance_transform_edt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--gen2", default="data/after/3dep2021_fulldensity.laz")
ap.add_argument("--slope", type=float, default=None,
                help="q2 = intercept + slope * cover. Default: READ from the tile's own "
                     "q2_cover_fit.json (analysis/ridgelines/q2_cover_fit.py), never typed. "
                     "It used to default to -0.1922, which was Elba's number AND stale -- "
                     "the shipped dod_cover_q2.json records -0.1835 -- so on any other "
                     "region it silently applied Elba's correction. Refuses if absent.")
ap.add_argument("--intercept", type=float, default=None,
                help="q2 at zero cover. Required with --slope, because a typed slope has no "
                     "fit to read it from. There is no default: 0.5 was one for years, and "
                     "it silently encoded the since-abandoned assumption that both epochs "
                     "see identical ground on bare earth.")
ap.add_argument("--relation", default=None,
                help="a relation JSON carrying intercept, slope AND the population they were "
                     "fitted on -- e.g. data/derived/control_q_ctrl_fit.json from "
                     "analysis/control_percentile_fit.py. When given, the return column and "
                     "the covariate are BUILT TO THAT SPEC rather than to this script's "
                     "defaults, because a rank measured in one population does not index "
                     "another. Without it the tile's own q2_cover_fit.json is used, on the "
                     "class-2 near-ground column, as before.")
ap.add_argument("--chunk", type=int, default=3_000_000)
ap.add_argument("--zlo", type=float, default=-1.0)
ap.add_argument("--zhi", type=float, default=2.0)
ap.add_argument("--dz", type=float, default=0.02)
A = ap.parse_args()

D = A.tile


def _q2_slope(tile_dir):
    """The q2 slope from THIS tile's own fit. The relation is per-site -- it depends on each
    pair's phenology and undergrowth -- so there is no site-invariant value to default to."""
    p = os.path.join(tile_dir, "q2_cover_fit.json")
    if not os.path.exists(p):
        raise SystemExit(
            f"no {p}. The q2 slope is NOT defaulted: it is per-site, and a value carried "
            f"from another tile would be applied here without saying so. Produce it:\n"
            f"    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python "
            f"analysis/ridgelines/q2_cover_fit.py --tile {tile_dir}\n"
            f"or state your own with --slope.")
    j = json.load(open(p))
    st, pop = j["settings"], j["population"]
    print(f"q2 slope {j['linear_slope']:+.4f}, read from {p}")
    print(f"  fitted on {pop['cells_fitted']:,} cells in {pop['bins_used_in_fit']} of "
          f"{pop['bins']} cover bins; weight={st['weight']}, include_valley="
          f"{st['include_valley']}, valley_top_m={st['valley_top_m']}")
    mt = j.get("inputs_mtime", {})
    if "corrections.json" in mt and "beam_offset_table.parquet" in mt:
        if mt["beam_offset_table.parquet"] < mt["corrections.json"]:
            print(f"  WARNING: that fit read a beam_offset_table ({mt['beam_offset_table.parquet']}) "
                  f"OLDER than corrections.json ({mt['corrections.json']}), so its gen1 "
                  f"offsets carry registration terms that are no longer in force.")
    return float(j["linear_slope"])


def _q2_intercept(tile_dir):
    """The tile fit's intercept. NOT defaulted. It was hardcoded to 0.5 here, which was right
    only while q2_cover_fit PINNED it; since 2026-09-02 it is fitted, and assuming 0.5 would
    apply half of a relation and silently drop the rest. A fit that predates the change does
    not carry the field, and the honest response to that is to refuse -- the old fit's SLOPE
    was estimated under the pin, so pairing it with any other intercept is not that fit."""
    p = os.path.join(tile_dir, "q2_cover_fit.json")
    j = json.load(open(p))
    if "linear_intercept" not in j:
        raise SystemExit(
            f"{p} carries no linear_intercept: it predates 2026-09-02, when the intercept "
            f"stopped being pinned to 0.5. Its slope was fitted UNDER that pin, so it cannot "
            f"be combined with a fitted intercept. Re-run "
            f"analysis/ridgelines/q2_cover_fit.py --tile {tile_dir} (add --pin-intercept to "
            f"reproduce the old relation deliberately).")
    return float(j["linear_intercept"])


# ---- the relation, and the population it obliges us to build -------------------------
REL = json.load(open(A.relation)) if A.relation else None
if REL is not None:
    if A.slope is not None:
        raise SystemExit("--slope and --relation are mutually exclusive: one types a "
                         "coefficient, the other reads it with its population attached.")
    INTERCEPT = float(REL["intercept"]); SLOPE = float(REL["slope"])
    PP = REL["percentile_population"]; CV = REL["covariate"]
    if "ALL" not in PP["classes"].upper():
        raise SystemExit(f"this producer builds an all-return or a class-2 column; the "
                         f"relation asks for {PP['classes']!r}.")
    COL_LO, COL_HI, COL_DZ = PP["window_lo_m"], PP["window_hi_m"], PP["bin_m"]
    COL_CLASSES = "all"
    COVER_NAME = CV["name"]
    print(f"relation: {REL['relation']}  intercept {INTERCEPT:+.4f}  slope {SLOPE:+.4f}")
    print(f"  read from {A.relation}, fitted on {REL['fitted_on']['marks']} "
          f"{REL['fitted_on']['set']} marks, {REL['weighting']}")
    print(f"  percentile population REQUIRED: {PP['classes']}, {COL_LO:+.2f}..{COL_HI:+.2f} m, "
          f"{COL_DZ:g} m bins -> building exactly that")
    print(f"  covariate: {COVER_NAME} = ({CV['band_lo_m']:g}, {CV['band_hi_m']:g}] m over "
          f"({CV['denominator_lo_m']:g}, {CV['denominator_hi_m']:g}] m, {CV['bin_m']:g} m bins, "
          f"{CV['test']}")
    print(f"  KNOWN SCALE DIFFERENCE: {REL['known_scale_difference']}")
    if REL["fitted_on"].get("lowveg_max_observed") is not None:
        print(f"  the relation was MEASURED over {COVER_NAME} 0.."
              f"{REL['fitted_on']['lowveg_max_observed']:.3f}; cells beyond that are "
              f"EXTRAPOLATION and are counted below")
else:
    if A.slope is not None and A.intercept is None:
        raise SystemExit("--slope needs --intercept: the pair defines the relation, and "
                         "there is no defensible default for the intercept alone.")
    SLOPE = A.slope if A.slope is not None else _q2_slope(D)
    INTERCEPT = A.intercept if A.intercept is not None else _q2_intercept(D)
    COL_LO, COL_HI, COL_DZ = A.zlo, A.zhi, A.dz
    COL_CLASSES = "class2"
    COVER_NAME = "canopy_cover_pfs"
    print(f"relation: q2 = {INTERCEPT:.4f} + {SLOPE:+.4f} * cover   "
          f"(tile fit; class-2 near-ground column)")
j = reg.read_corrections(D)   # geoid sidecar wins where a tile carries both;
                              # reading "corrections.json" by name picks elbaext's
                              # obsolete reference_plane product
b = j["bounds"]; RES = float(j["res_m"]); X0, Y0 = b[0], b[1]
zf = np.load(f"{D}/z_after.npy"); NY, NX = zf.shape; NC = zf.size
_zf = zf.copy(); _m = ~np.isfinite(_zf)
if _m.any():
    _zf = _zf[tuple(distance_transform_edt(_m, return_distances=False, return_indices=True))]
gy, gx = np.gradient(_zf, RES)
gxf = gx.ravel(); gyf = gy.ravel(); zflat = _zf.ravel()
nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)
NZ = int(round((COL_HI - COL_LO) / COL_DZ))

# ---- gen2, ONE streaming pass -------------------------------------------------------
# Two accumulators, because the relation's covariate and its rank live in different
# windows: the percentile column at the relation's own resolution, and (when the relation
# needs it) a 0.02 m near-ground column from which lowveg is computed by exactly the
# bin-centre rule scripts/make_lowveg.py uses. Building both here is what lets a lowveg
# relation be applied over the WHOLE grid -- lowveg.npy itself exists only on the divide
# cells of the near-ground cube, which is enough to fit a relation and not to apply one.
H = np.zeros((NC, NZ), np.int32)
LV_EDGES = None
if REL is not None:
    LV_EDGES = np.arange(CV["denominator_lo_m"], CV["denominator_hi_m"] + 0.5 * CV["bin_m"],
                         CV["bin_m"])
    LVH = np.zeros((NC, LV_EDGES.size - 1), np.int32)
n_in = 0
n_lv = 0
with laspy.open(A.gen2) as f:
    for pts in f.chunk_iterator(A.chunk):
        cl = np.asarray(pts.classification)
        keep = np.ones(cl.shape, bool) if COL_CLASSES == "all" else (cl == 2)
        x = np.asarray(pts.x)[keep]; y = np.asarray(pts.y)[keep]; z = np.asarray(pts.z)[keep]
        ix = ((x - X0) / RES).astype(np.int64); iy = ((y - Y0) / RES).astype(np.int64)
        ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
        cc = iy[ing] * NX + ix[ing]
        xc = X0 + ((cc % NX) + 0.5) * RES; yc = Y0 + ((cc // NX) + 0.5) * RES
        h = (z[ing] - (zflat[cc] + gxf[cc] * (x[ing] - xc) + gyf[cc] * (y[ing] - yc))) / nnorm[cc]
        zi = np.floor((h - COL_LO) / COL_DZ).astype(np.int64)
        m = (zi >= 0) & (zi < NZ)
        np.add.at(H, (cc[m], zi[m]), 1)
        n_in += int(m.sum())
        if LV_EDGES is not None:
            li = np.floor((h - LV_EDGES[0]) / CV["bin_m"]).astype(np.int64)
            ml = (li >= 0) & (li < LV_EDGES.size - 1)
            np.add.at(LVH, (cc[ml], li[ml]), 1)
            n_lv += int(ml.sum())
print(f"gen2: {n_in:,} {COL_CLASSES} returns in {COL_LO:+.2f}..{COL_HI:+.2f} m over "
      f"{int((H.sum(1) > 0).sum()):,} cells", flush=True)
if LV_EDGES is not None:
    print(f"      {n_lv:,} returns in the covariate window over "
          f"{int((LVH.sum(1) > 0).sum()):,} cells", flush=True)

# ---- gen1: per-cell median of the registered per-return offsets ----------------------
t = pq.read_table(f"{D}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"])
g = t["in_grid"].to_numpy().astype(bool)
ce = t["cell"].to_numpy()[g]; dc = t["d_mm_corr"].to_numpy()[g].astype(float)
o = np.argsort(ce, kind="stable"); cs = ce[o]; ds = dc[o]
u, st = np.unique(cs, return_index=True); sp = np.r_[st[1:], cs.size]
h1 = np.full(NC, np.nan)
for uu, a, bb in zip(u, st, sp):
    h1[uu] = np.median(ds[a:bb])
print(f"gen1: ground on {int(np.isfinite(h1).sum()):,} cells", flush=True)

# ---- gen2 at the cover-dependent percentile -----------------------------------------
if REL is not None:
    mid = 0.5 * (LV_EDGES[:-1] + LV_EDGES[1:])
    band = (mid > CV["band_lo_m"]) & (mid <= CV["band_hi_m"])   # bin CENTRE, as at the marks
    tot = LVH.sum(1).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        cover = np.where(tot > 0, LVH[:, band].sum(1) / tot, np.nan)
    np.save(f"{D}/{COVER_NAME}_grid.npy", cover.reshape(NY, NX))
    _lvmax = REL["fitted_on"].get("lowveg_max_observed")
    if _lvmax:
        n_ex = int(np.sum(np.isfinite(cover) & (cover > _lvmax)))
        print(f"cover: {int(np.isfinite(cover).sum()):,} cells; {n_ex:,} "
              f"({100*n_ex/max(int(np.isfinite(cover).sum()),1):.1f}%) exceed the "
              f"{_lvmax:.3f} the relation was measured to -- EXTRAPOLATED, not dropped")
else:
    cover = np.load(f"{D}/canopy_cover_pfs.npy").ravel()
q2 = INTERCEPT + SLOPE * np.where(np.isfinite(cover), cover, 0.0)
_oob = int(np.sum((q2 < 0) | (q2 > 1)))
if _oob:
    print(f"q2 outside [0,1] on {_oob:,} cells; clipped to the column's ends, which is a "
          f"FLOOR/CEILING not a fit -- those cells take the extreme return, not a percentile")
q2 = np.clip(q2, 0.0, 1.0)
C = np.cumsum(H, 1).astype(float); ntot = C[:, -1]
have = ntot > 0
idx = np.arange(NC)
r = q2 * ntot
k = (C >= r[:, None]).argmax(1)
below = np.where(k > 0, C[idx, np.maximum(k - 1, 0)], 0.0)
inbin = C[idx, k] - below
frac = np.where(inbin > 0, (r - below) / np.maximum(inbin, 1e-9), 0.0)
h2 = np.where(have, (COL_LO + (k + np.clip(frac, 0, 1)) * COL_DZ) * 1000.0, np.nan)
h2_med = np.where(have, (COL_LO + ((C >= 0.5 * ntot[:, None]).argmax(1) + 0.5) * COL_DZ) * 1000.0,
                  np.nan)

dod_corr = (h2 - h1) / 1000.0 * nnorm            # m, positive = elevation rose
dod_med = (h2_med - h1) / 1000.0 * nnorm         # same but gen2 at its plain median
np.save(f"{D}/dod_cover_q2.npy", dod_corr.reshape(NY, NX))
np.save(f"{D}/dod_gen2_median.npy", dod_med.reshape(NY, NX))
np.save(f"{D}/gen2_q2_used.npy", np.where(have, q2, np.nan).reshape(NY, NX))
json.dump({"relation": f"q2 = {INTERCEPT:.6f} + slope * {COVER_NAME}",
           "intercept": INTERCEPT, "slope": SLOPE,
           "relation_source": A.relation or f"{D}/q2_cover_fit.json",
           "percentile_population": {"classes": COL_CLASSES, "window_lo_m": COL_LO,
                                     "window_hi_m": COL_HI, "bin_m": COL_DZ},
           "source": "analysis/ridgelines/Q2_COVER_RELATION.md",
           "gen1": "beam_offset_table.parquet median of d_mm_corr (4 registration terms)",
           "gen2": f"{A.gen2}, {COL_CLASSES}, column {COL_LO}..{COL_HI} m, {COL_DZ} m bins",
           "sign": "gen2 - gen1, positive = elevation rose",
           "cells_with_both": int((np.isfinite(dod_corr)).sum())},
          open(f"{D}/dod_cover_q2.json", "w"), indent=2)
ok = np.isfinite(dod_corr) & np.isfinite(dod_med)
print(f"\nwrote {D}/dod_cover_q2.npy ({int(ok.sum()):,} cells)")
print(f"      {D}/dod_gen2_median.npy   {D}/gen2_q2_used.npy   {D}/dod_cover_q2.json")
print(f"\nmedian DoD (mm): gen2 median {np.median(dod_med[ok])*1000:+.1f}  ->  "
      f"cover-corrected {np.median(dod_corr[ok])*1000:+.1f}")
old = np.load(f"{D}/dod.npy").ravel()
oo = ok & np.isfinite(old)
print(f"existing dod.npy on the same cells: {np.median(old[oo])*1000:+.1f} mm  "
      f"(n={int(oo.sum()):,})")
