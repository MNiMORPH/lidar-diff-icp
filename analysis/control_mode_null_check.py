#!/usr/bin/env python3
"""NULL CHECK: on open, well-observed ground, does the half-sample mode agree with the median?

Andy, 2026-10-16: the control marks are not on the floodplain, so they cannot validate the
mode for VEGETATION. What they can do is falsify the estimator outright: on open ground the
mode and the median target the same population, so if they disagree there the estimator is
broken everywhere -- including in the coarsened-gen1 floodplain result that already leans
on it. This is a null check, not a validation.

TWO EPOCHS, TWO DATA PATHS, SAME ESTIMATOR.

  gen2  every control mark with a stored near-ground histogram (``control_structure``),
        both control sets. The boxes at the 2008 marks are gen2 clouds too (verified: LAS
        point format 7, GPS time 3.03e8 = 2021), so at those marks the surveyed level is on
        GEOID03 against a GEOID18 cloud and only hsm - median is meaningful. At gen2's own
        marks (GEOID18 both) hsm - surveyed is meaningful too.
  gen1  every 2008 control mark whose gen1 tile is on disk under data/before (29 L1O open
        marks plus 58 of other covers). Returns read from the tile, order-2 local surface
        through vendor class-2 within R, slope-normal residuals of ALL returns, histogram
        at the same 20 mm bins the floodplain work used. gen1 and its 2008 control share
        NAVD88(GEOID03), so hsm - surveyed is meaningful; the known delivered tie (~ +63
        mm reads-low at Elba's open marks) is a LEVEL and is gauge here.

The histogram half-sample mode is the SAME function the coarsened-gen1 experiment used
(``gen1_coarsened_ground.hsm_from_hist``), so a disagreement here is about the estimator,
not about an implementation difference. For gen1 the raw-sample half-sample mode is
reported alongside to show what the 20 mm binning costs.

TWO RADII FOR gen1. R = 7.5 m matches the stored gen2 boxes (``struct_radius``). R = 2.5 m
is one 5 m cell, where gen1 holds ~16 returns -- the count at which the resampling
experiment measured the mode biased +19.7 mm IN HEAVY VEGETATION. On open ground that
bias should not appear; if it does, it is a small-n property of the estimator, not a
vegetation effect.

    ./lidar-icp/bin/python analysis/control_mode_null_check.py
"""
import glob
import os
import sys

import laspy
import numpy as np
import pandas as pd

_REPO = os.path.dirname(os.path.abspath(__file__))
while _REPO != "/" and not os.path.exists(os.path.join(_REPO, "pyproject.toml")):
    _REPO = os.path.dirname(_REPO)
os.chdir(_REPO)
sys.path.insert(0, os.path.join(_REPO, "src"))
from lidar_diff_icp.groundtruth.tie import _design          # noqa: E402
from lidar_diff_icp.tiles import find_tile                  # noqa: E402
sys.path.insert(0, os.path.join(_REPO, "analysis"))

STRUCT = "data/derived/control_structure"
CONTROL = {"gen1_2008_control": "src/lidar_diff_icp/groundtruth/data/mn_dnr_2008_control_semn.csv",
           "gen2_2021_control": "src/lidar_diff_icp/groundtruth/data/mn_se_driftless_2021_control.csv"}
TILES = "data/before"
LO, HI, DZ = -1.0, 2.0, 0.02        # m: the near-ground cube window and its bin (stored gen2 histograms)
RADII = (7.5, 2.5)                  # m: gen2 box radius (struct_radius) and one 5 m cell
OUT = "data/derived/control_mode_null_check.csv"

print("INPUTS AND ROLES")
print(f"  {STRUCT}/*.npz            gen2 near-ground histograms at every mark, both sets; bins ng_edges")
print(f"  {CONTROL['gen1_2008_control']}   2008 marks: point_type, surveyed z (NAVD88 GEOID03)")
print(f"  {CONTROL['gen2_2021_control']}   2021 marks: point_type, surveyed z (NAVD88 GEOID18)")
print(f"  {TILES}/<tile>.laz                 gen1 2008 returns, read whole per tile, all classes but 7")
print("PARAMETERS")
print(f"  window {LO}..{HI} m, bin {DZ*1000:g} mm  -- the stored histograms' own edges; same bin as the floodplain hsm")
print(f"  gen1 radii {RADII} m  -- box radius, and one 5 m cell")
print("COLUMNS (mm unless stated): n = returns in window; med = median; hsm = histogram half-sample mode;")
print("  hsm_raw = half-sample mode on the unbinned returns (gen1 only); *_surv = estimator minus surveyed ground")


from hsm import half_sample_mode as hsm_raw, hsm_from_hist   # noqa: E402  tie-symmetric


def med_from_hist(counts, ctr):
    c = np.cumsum(counts); tot = c[-1]
    return float(ctr[np.searchsorted(c, 0.5 * tot)]) if tot > 0 else np.nan


rows = []
# ------------------------------------------------------------------ gen2: stored histograms
for set_, path in CONTROL.items():
    M = pd.read_csv(path).drop_duplicates(subset=["easting", "northing", "elevation"])
    types = dict(zip(M.point_id, M.point_type)); surv = dict(zip(M.point_id, M.elevation))
    n_ok = 0
    for f in sorted(glob.glob(f"{STRUCT}/{set_}__*.npz")):
        pid = os.path.basename(f)[len(set_) + 2:-4]
        if pid not in surv:
            continue
        z = np.load(f)
        e = z["ng_edges"]; ctr = 0.5 * (e[:-1] + e[1:])
        coef = z["surface_coef"]; nn = np.sqrt(1.0 + coef[1] ** 2 + coef[2] ** 2)
        mu_true = (surv[pid] - coef[0]) / nn
        H = z["ng_all"].astype(int)
        if H.sum() < 20:
            continue
        rows.append(dict(epoch="gen2", set=set_, point_id=pid, type=types[pid], radius=float(z["struct_radius"]),
                         n=int(H.sum()), med=med_from_hist(H, ctr) * 1e3, hsm=hsm_from_hist(H, ctr) * 1e3,
                         hsm_raw=np.nan, surv=mu_true * 1e3))
        n_ok += 1
    print(f"gen2 histograms used: {n_ok} marks of {set_}")

# ------------------------------------------------------------------ gen1: tiles on disk
M1 = pd.read_csv(CONTROL["gen1_2008_control"]).drop_duplicates(subset=["easting", "northing", "elevation"])
by_tile = {}
for t in M1.itertuples():
    tl = find_tile(t.easting, t.northing)
    p = f"{TILES}/{tl}.laz"
    if os.path.exists(p):
        by_tile.setdefault(p, []).append(t)
print(f"gen1 tiles on disk holding 2008 marks: {len(by_tile)}; marks: {sum(len(v) for v in by_tile.values())}")
edges = np.arange(LO, HI + DZ / 2, DZ); ctr1 = 0.5 * (edges[:-1] + edges[1:])
for p, mk in sorted(by_tile.items()):
    f = laspy.read(p)
    x, y, zz = np.asarray(f.x), np.asarray(f.y), np.asarray(f.z)
    cl = np.asarray(f.classification)
    keep = cl != 7
    for t in mk:
        for R in RADII:
            r = np.hypot(x - t.easting, y - t.northing)
            m = (r <= R) & keep
            g = m & (cl == 2)
            if g.sum() < 6:
                continue
            A = _design(x[g] - t.easting, y[g] - t.northing, 2)
            coef, *_ = np.linalg.lstsq(A, zz[g], rcond=None)
            nn = np.sqrt(1.0 + coef[1] ** 2 + coef[2] ** 2)
            h = (zz[m] - _design(x[m] - t.easting, y[m] - t.northing, 2) @ coef) / nn
            h = h[(h > LO) & (h <= HI)]
            if h.size < 5:
                continue
            H = np.histogram(h, bins=edges)[0]
            rows.append(dict(epoch="gen1", set="gen1_2008_control", point_id=t.point_id, type=t.point_type,
                             radius=R, n=int(h.size), med=float(np.median(h)) * 1e3,
                             hsm=hsm_from_hist(H, ctr1) * 1e3, hsm_raw=hsm_raw(h) * 1e3,
                             surv=(t.elevation - coef[0]) / nn * 1e3))
    print(f"  {os.path.basename(p)}: {len(mk)} marks", flush=True)

F = pd.DataFrame(rows)
F["hsm_med"] = F.hsm - F.med
F["med_surv"] = F.med - F.surv
F["hsm_surv"] = F.hsm - F.surv
F["hsmraw_med"] = F.hsm_raw - F.med
F.to_csv(OUT, index=False)


def nmad(v):
    v = np.asarray(v); v = v[np.isfinite(v)]
    return 1.4826 * np.median(np.abs(v - np.median(v))) if v.size else np.nan


def table(G, label, cols):
    print(f"\n{label}")
    print(f"  {'type':<6}{'n marks':>8}{'n ret':>7}" + "".join(f"{c:>14}" for c in cols))
    print(f"  {'':<6}{'':>8}{'median':>7}" + "".join(f"{'med  (NMAD)':>14}" for c in cols))
    order = ["L1O", "NVA", "L2T", "L3B", "L4F", "L5U", "VVA", "LCP", "other"]
    for ty in [o for o in order if o in set(G.type)]:
        g = G[G.type == ty]
        print(f"  {ty:<6}{len(g):>8}{np.median(g.n):>7.0f}"
              + "".join(f"{np.nanmedian(g[c]):>7.1f} ({nmad(g[c]):>4.1f})" for c in cols))


print("\n" + "=" * 90)
print("THE NULL: on OPEN ground (L1O, NVA) hsm - med must be ~0 and hsm - surveyed must equal med - surveyed.")
print("Vegetated types are the contrast: there hsm - med is EXPECTED to be negative (the mode sheds the tail).")
table(F[(F.epoch == "gen2") & (F.set == "gen2_2021_control")],
      "gen2 cloud at gen2's OWN marks (GEOID18 both): hsm - med, and each against surveyed",
      ["hsm_med", "med_surv", "hsm_surv"])
table(F[(F.epoch == "gen2") & (F.set == "gen1_2008_control")],
      "gen2 cloud at the 2008 marks (surveyed on GEOID03 -- its LEVEL is gauge; read hsm_med and the NMADs)",
      ["hsm_med", "med_surv", "hsm_surv"])
for R in RADII:
    table(F[(F.epoch == "gen1") & (F.radius == R)],
          f"gen1 returns at the 2008 marks, R = {R} m (GEOID03 both)",
          ["hsm_med", "hsmraw_med", "med_surv", "hsm_surv"])
print(f"\nwrote {OUT}")
