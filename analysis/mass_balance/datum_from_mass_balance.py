#!/usr/bin/env python3
"""How large a UNIFORM vertical RAISE of gen1 does routed sediment continuity demand?

The DoD is ``gen2 - gen1`` (positive = elevation rose). Raising gen1 by ``delta`` lowers
the DoD by ``delta`` everywhere, so ``V_acc`` -- the volume budget accumulated down the
flow network -- moves by ``-delta * area * N_up(c)`` at every cell, where ``N_up`` is the
routed count of KNOWN upstream cells (``mass_balance``'s ``N_up``; that is exactly its
documented sensitivity to a uniform vertical shift). The error envelope ``sigma_Vacc``
does NOT move: it is built from the per-cell error, the routed weights and the DoD's MAD,
none of which a uniform translation changes. That is asserted in code below, not assumed.

Consequences used here:

* the whole delta sweep is ONE routed budget plus arithmetic, so the curve is exact
  rather than sampled;
* the continuity surplus per unit upstream area,

      S(c) = ( V_acc(c) - V_in_acc(c) - z*sigma_Vacc(c) ) / ( area * N_up(c) )   [m]

  is the shift each cell demands on its own. Cell ``c`` stops being flagged as soon as
  ``delta >= S(c)``. Therefore
      delta needed to leave a fraction ``f`` of cells still flagged  =  quantile(S, 1-f),
  and ``catchment_dod_balance.balancing_offset(keep=q)`` returns ``-max(quantile(S,q),0)``
  -- the same curve, parameterised by ``keep = 1 - f``. Both are printed, and the
  agreement between them is checked, so the reuse is verified rather than claimed.

There is NO single "the" delta: the answer is a curve, and where you read it off is the
caller's choice, not this script's.

Everything about the error envelope, the off-map boundary terms and the variogram follows
``catchment-dod-balance/scripts/validate_site.py`` so the two are the same method.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/mass_balance/datum_from_mass_balance.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.abspath("src"))
sys.path.insert(0, os.path.expanduser("~/projects/catchment-dod-balance"))

from trust.provenance import Run                                       # noqa: E402
from lidar_diff_icp.variogram import empirical_variogram, fit_spherical  # noqa: E402
from catchment_dod_balance import (dinf_proportions, mass_balance,      # noqa: E402
                                   balancing_offset)
from catchment_dod_balance.offmap import offmap_terms, contributing_area  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--dod", default="dod_cover_q2.npy")
ap.add_argument("--lod", default="lod_cover_q2.npy")
ap.add_argument("--dem", default="z_after.npy")
ap.add_argument("--floodplain", default="floodplain_mask.npy")
ap.add_argument("--floodplain-mode", default="routing",
                choices=["routing", "constraint", "none"],
                help="routing: floodplain cells leave the flow graph (hillslope budget); "
                     "constraint: they route but cannot set the datum; none: kept")
ap.add_argument("--envelope", default="correlated", choices=["correlated", "independent"])
ap.add_argument("--channel-area", type=float, default=1.0e4)
ap.add_argument("--max-lag", type=float, default=600.0)
ap.add_argument("--z", type=float, default=1.96)
ap.add_argument("--n-boot", type=int, default=400)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--save", default="")
A = ap.parse_args()

D = A.tile
R = Run("what uniform vertical raise of gen1 does routed sediment continuity demand, "
        "on the hillslopes only?")
dem_p = R.input(f"{D}/{A.dem}", role="DEM routed (gen2 ground surface, m)")
dod_p = R.input(f"{D}/{A.dod}", role="DEM of Difference, gen2 - gen1, m (+ = rose)")
lod_p = R.input(f"{D}/{A.lod}", role="level of detection, 1.96 sigma, m")
fld_p = (R.input(f"{D}/{A.floodplain}", role="floodplain cells excluded from the hillslope "
                 "budget (STANDING PROJECT RULE)")
         if A.floodplain_mode != "none" else None)
meta_p = R.input(f"{D}/corrections.json", role="pipeline metadata: res_m, stable_1sigma_m")

meta = json.load(open(meta_p))
RES = float(meta["res_m"]); AREA = RES * RES
R.param("res_m", RES, src="repo")
R.param("z", A.z, src="repo", why="the shipped surplus test z in mass_balance()")
R.param("channel_area_m2", A.channel_area, src="repo",
        why="catchment-dod-balance offmap default, as used in its own validation run")
R.param("variogram_max_lag_m", A.max_lag, src="repo", why="validate_site.py default")
R.param("floodplain_mode", A.floodplain_mode, src="repo",
        why="standing project rule: floodplain_mask cells out of hillslope mass balance")
R.param("envelope", A.envelope, src="repo",
        why="catchment-dod-balance README: independent is a documented LOWER BOUND")
R.param("n_boot", A.n_boot, src="MINE", why="block-bootstrap replicates for the SE; "
        "raise it if the SE is not stable")

z_dem = np.load(dem_p); dod = np.load(dod_p); lod = np.load(lod_p)
flood = np.load(fld_p).astype(bool) if fld_p else np.zeros(dod.shape, bool)

props, valid_all = dinf_proportions(z_dem, breach=True)
carea_all = contributing_area(props, valid_all, RES)

# ---- error structure: fitted, not assumed (validate_site.py, verbatim method) --------
ny, nx = dod.shape
yy, xx = np.mgrid[0:ny, 0:nx]
stable = np.isfinite(dod) & np.isfinite(lod) & (np.abs(dod) <= lod)
centers, gamma, counts = empirical_variogram(
    (xx[stable] * RES).astype(float), (yy[stable] * RES).astype(float), dod[stable],
    max_lag=A.max_lag, n_lags=25, n_pairs=600_000, estimator="dowd", seed=0)
vgm = fit_spherical(centers, gamma, counts)
perror_total = lod / A.z
perror_nugget = perror_total * np.sqrt(vgm.nugget / vgm.total_sill)

# ---- off-map boundary terms (bounded hillslope area vs unbounded channel) ------------
terms = offmap_terms(z_dem, props, valid_all, RES, channel_area=A.channel_area)
cls = terms["classification"]

# ---- the hillslope routing graph ----------------------------------------------------
valid = valid_all & ~flood if A.floodplain_mode == "routing" else valid_all

R.mask("floodplain", flood, of=flood.size,
       defn=f"{A.floodplain} as shipped (TPI over an 800 m window < -2 m; "
            "analysis/ridgelines/convexity_dod_landcover.py)")
R.mask("routed", valid, of=valid.size,
       defn=("D-infinity graph after breaching, floodplain cells REMOVED"
             if A.floodplain_mode == "routing" else "D-infinity graph after breaching"))
R.mask("variogram_cells", stable, of=stable.size,
       defn="|dod| <= lod (insignificant change), the validate_site.py selection")

kw = dict(unaccounted_area=terms["unaccounted_area"], contaminate=terms["contaminate"])
if A.envelope == "correlated":
    out = mass_balance(dod, perror_nugget, props, valid, RES, z=A.z,
                       corr_sill=float(vgm.sill), corr_range=float(vgm.range_), **kw)
else:
    out = mass_balance(dod, perror_total, props, valid, RES, z=A.z, **kw)

V = out["V_acc"]; SIG = out["sigma_Vacc"]; NUP = out["N_up"]; VIN = out["V_in_acc"]
known = out["known"]; contam = out["contaminated"]

ev = valid & np.isfinite(V) & np.isfinite(SIG) & ~contam & (NUP > 0)
R.mask("evaluable", ev, of=int(valid.sum()),
       defn="routed & finite V_acc & not off-map-contaminated & N_up > 0")

# S(c): the shift, in metres, that cell c on its own demands. delta >= S(c) clears it.
S = np.full(dod.shape, np.nan)
S[ev] = (V[ev] - VIN[ev] - A.z * SIG[ev]) / (AREA * NUP[ev])
s = S[ev]

# ---- assert the two things the whole sweep rests on ---------------------------------
_probe = 0.050
_o2 = (mass_balance(dod - _probe, perror_nugget, props, valid, RES, z=A.z,
                    corr_sill=float(vgm.sill), corr_range=float(vgm.range_), **kw)
       if A.envelope == "correlated" else
       mass_balance(dod - _probe, perror_total, props, valid, RES, z=A.z, **kw))
d_sig = float(np.nanmax(np.abs(_o2["sigma_Vacc"][ev] - SIG[ev])))
d_V = float(np.nanmax(np.abs(_o2["V_acc"][ev] - (V[ev] - _probe * AREA * NUP[ev]))))
_bo = balancing_offset(dod, perror_nugget if A.envelope == "correlated" else perror_total,
                       props, valid, RES, z=A.z,
                       corr_sill=float(vgm.sill) if A.envelope == "correlated" else 0.0,
                       corr_range=float(vgm.range_) if A.envelope == "correlated" else None,
                       exclude=flood if A.floodplain_mode == "constraint" else None,
                       keep=0.99, **kw)

R.column("f_flagged", "fraction of evaluable cells still flagged as unphysical deposition "
                      "(V_acc - V_in > z*sigma) after the raise, dimensionless")
R.column("keep", "1 - f_flagged; the balancing_offset 'keep' quantile that returns this delta")
R.column("delta_mm", "uniform raise of gen1 (mm) = quantile(S, keep); "
                     "= -1000 * balancing_offset(keep)['offset']")
R.column("se_mm", "SE of delta_mm: SD over block-bootstrap replicates of the same "
                  "quantile, blocks of side L (mm)")
R.column("net_vol_m3", "sum over KNOWN evaluable cells of (dod - delta)*area, m^3; "
                       "negative = the hillslope exports sediment")
R.column("L_m", "side of the square spatial block resampled by the bootstrap, m")
R.column("n_blocks", "number of blocks holding at least one evaluable cell")
R.column("band", "drainage-area band of the evaluable cells, ha")
R.column("n", "count of evaluable cells in the band")
R.column("median_S_mm", "median over the band of S = the raise each cell demands, mm")
R.column("f0", "fraction of the band flagged at delta = 0, dimensionless")
R.column("quantity", "name of the reported scalar")
R.column("value", "its value, units given in the name")

R.banner()

print(f"\ngrid {dod.shape}  res {RES:g} m   routed {valid_all.mean():.2%} of cells "
      f"(largest drainage {carea_all.max()/1e6:.2f} km^2)")
print(f"floodplain mask: {flood.mean():.2%} of the grid, mode={A.floodplain_mode}")
print(f"variogram on {stable.mean():.0%} insignificant-change cells: nugget sd "
      f"{np.sqrt(vgm.nugget):.4f} m, sill sd {np.sqrt(vgm.sill):.4f} m, range "
      f"{vgm.range_:.0f} m (total {np.sqrt(vgm.total_sill):.4f} m; pipeline's "
      f"stable_1sigma_m {meta.get('stable_1sigma_m')})")
print(f"boundary: {int(cls['entry'].sum())} crossings, {int(cls['channel_entry'].sum())} "
      f"off-map trunks; hillslope length {cls['hillslope_length']:.0f} m")
print(f"evaluable {int(ev.sum()):,} of {int(valid.sum()):,} routed ({ev.sum()/valid.sum():.2%}); "
      f"known DoD on {int((ev & known).sum()):,} of them")
print(f"\nCHECK  sigma_Vacc invariant under a {1000*_probe:.0f} mm uniform shift: "
      f"max |change| {d_sig:.3e} m^3")
print(f"CHECK  V_acc shift equals -delta*area*N_up: max |residual| {d_V:.3e} m^3")
print(f"CHECK  balancing_offset(keep=0.99) = {-1000*_bo['offset']:+.3f} mm vs "
      f"quantile(S,0.99) = {1000*max(np.quantile(s,0.99),0.0):+.3f} mm; "
      f"coverage {_bo['evaluable_fraction']:.1%}, binding cell {_bo['binding_cell']}")

# ---- the curve ----------------------------------------------------------------------
FS = [0.0, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20,
      0.25, 0.30, 0.40, 0.50]
rng = np.random.default_rng(A.seed)
BLOCKS = [25.0, 50.0, 100.0, 200.0, 400.0, 800.0]
rows_i, cols_i = np.where(ev)


def boot_se(q, L):
    """SD over bootstrap replicates of quantile(S, q), resampling whole LxL blocks."""
    bi = (rows_i * RES / L).astype(np.int64); bj = (cols_i * RES / L).astype(np.int64)
    bid = bi * (int(nx * RES / L) + 2) + bj
    uniq, inv = np.unique(bid, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    starts = np.searchsorted(inv[order], np.arange(uniq.size))
    ends = np.searchsorted(inv[order], np.arange(uniq.size), side="right")
    vals = np.empty(A.n_boot)
    for b in range(A.n_boot):
        pick = rng.integers(0, uniq.size, uniq.size)
        idx = np.concatenate([order[starts[p]:ends[p]] for p in pick])
        vals[b] = np.quantile(s[idx], q)
    return float(np.std(vals, ddof=1)), int(uniq.size)


print("\nCURVE  uniform raise of gen1 vs the fraction of evaluable hillslope cells left "
      "flagged\n       (delta = quantile(S, keep); the SE column uses L = "
      f"{BLOCKS[-2]:.0f} m blocks)")
L_ref = BLOCKS[-2]
rows = []
for f in FS:
    q = 1.0 - f
    dmm = 1000.0 * float(np.quantile(s, q))
    se, nb = boot_se(q, L_ref)
    net = float(np.sum((dod[ev & known] - dmm / 1000.0) * AREA))
    rows.append([f"{f:.4f}", f"{q:.4f}", f"{dmm:+.1f}", f"{1000*se:.1f}", f"{net:+,.0f}"])
R.table(["f_flagged", "keep", "delta_mm", "se_mm", "net_vol_m3"], rows)

print("\nBLOCK SIZE  SE of delta_mm vs the side of the resampled block "
      "(spatial autocorrelation is what this is chasing)")
rows = []
for L in BLOCKS:
    r = [f"{L:.0f}"]
    nb = None
    for f in (0.025, 0.05, 0.10):
        se, nb = boot_se(1.0 - f, L)
        r.append(f"{1000*se:.1f}")
    r.append(str(nb))
    rows.append(r)
R.column("se_f0025_mm", "SE of delta_mm at f_flagged = 0.025, mm")
R.column("se_f005_mm", "SE of delta_mm at f_flagged = 0.05, mm")
R.column("se_f010_mm", "SE of delta_mm at f_flagged = 0.10, mm")
R.table(["L_m", "se_f0025_mm", "se_f005_mm", "se_f010_mm", "n_blocks"], rows)

print("\nBY DRAINAGE AREA  is the demand uniform down the network, or structured?")
edges = np.nanquantile(carea_all[ev], [0, 0.5, 0.8, 0.95, 0.99, 1.0])
rows = []
for i in range(len(edges) - 1):
    lo, hi = edges[i], edges[i + 1]
    last = i == len(edges) - 2
    sel = ev & (carea_all >= lo) & ((carea_all <= hi) if last else (carea_all < hi))
    if sel.sum() < 50:
        continue
    rows.append([f"{lo/1e4:.2f}-{hi/1e4:.2f}", str(int(sel.sum())),
                 f"{1000*np.nanmedian(S[sel]):+.1f}", f"{np.mean(S[sel] > 0):.4f}"])
R.table(["band", "n", "median_S_mm", "f0"], rows)

# ---- the aggregate (net-export) criterion -------------------------------------------
kn = ev & known
mean_dod = float(np.mean(dod[kn]))
print("\nAGGREGATE CRITERION  the whole evaluated hillslope must not gain sediment: "
      "sum (dod - delta) * area <= 0")
rows = [["mean DoD over evaluable+known hillslope cells (mm)", f"{1000*mean_dod:+.2f}"],
        ["delta that zeroes the hillslope net volume (mm)", f"{1000*mean_dod:+.2f}"],
        ["hillslope net volume at delta = 0 (m^3)", f"{np.sum(dod[kn])*AREA:+,.0f}"],
        ["cells in that sum", f"{int(kn.sum()):,}"]]
R.table(["quantity", "value"], rows)

# block SE of the mean DoD (the aggregate criterion's own uncertainty)
print("\n            SE of that mean-DoD delta, by block side (mm)")
rk, ck = np.where(kn)
dk = dod[kn]
rows = []
for L in BLOCKS:
    bid = (rk * RES / L).astype(np.int64) * (int(nx * RES / L) + 2) + (ck * RES / L).astype(np.int64)
    uniq, inv = np.unique(bid, return_inverse=True)
    bm = np.bincount(inv, weights=dk) / np.bincount(inv)
    rows.append([f"{L:.0f}", f"{1000*np.mean(bm):+.2f}",
                 f"{1000*np.std(bm, ddof=1)/np.sqrt(uniq.size):.2f}", str(uniq.size)])
R.column("mean_block_dod_mm", "mean over blocks of the block-mean DoD, mm")
R.column("se_mean_mm", "SE of that mean of block means = SD/sqrt(n_blocks), mm")
R.table(["L_m", "mean_block_dod_mm", "se_mean_mm", "n_blocks"], rows)

if A.save:
    np.savez_compressed(A.save, S=S, evaluable=ev, known=known, N_up=NUP,
                        V_acc=V, sigma_Vacc=SIG, contaminated=contam,
                        carea=carea_all, dod=dod,
                        vgm=np.array([vgm.nugget, vgm.sill, vgm.range_]))
    print(f"\nwrote {A.save}")

R.done(headline=f"delta(f=0.05) = {1000*np.quantile(s,0.95):+.1f} mm on {A.dod}, "
                f"{A.envelope} envelope, floodplain {A.floodplain_mode}")
