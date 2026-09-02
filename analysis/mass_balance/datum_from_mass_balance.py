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
ap.add_argument("--z", type=float, default=1.0,
                help="how many sigma the budget is allowed to sit above zero; "
                     "1.0 = the edge of 1 sigma (Andy 2026-08-26), 1.96 = the "
                     "shipped mass_balance default")
ap.add_argument("--n-boot", type=int, default=400)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--null-real", type=int, default=0,
                help="realisations of the FITTED error model, routed, to measure the flag "
                     "rate the envelope itself produces when there is no real surplus")
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
cover_p = os.path.join(D, "canopy_cover_pfs.npy")
cover = (np.load(R.input(cover_p, role="gen2 canopy cover fraction (PyForestScan), "
                         "for the residual-structure table"))
         if os.path.exists(cover_p) else None)

meta = json.load(open(meta_p))
RES = float(meta["res_m"]); AREA = RES * RES
R.param("res_m", RES, src="repo")
R.param("z", A.z, src="andy", why="the edge of 1 sigma, asked for 2026-08-26; "
        "mass_balance ships 1.96")
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
# An all-False mask is "measured, no floodplain anywhere", which is a claim. The standing
# project rule is that floodplain cells stay OUT of hillslope mass balance, so a run without
# the mask is a DIFFERENT population and must not be produced by accident.
# refcells.reference_cells refuses on exactly this; so does this.
if not fld_p:
    raise SystemExit(
        "no floodplain mask given. An empty one is not the same thing: it asserts that no "
        "cell is floodplain, and the standing rule keeps floodplain cells out of hillslope "
        "mass balance, so the two populations are not comparable. Pass the mask, or state "
        "that you mean to run over ALL cells.")
flood = np.load(fld_p).astype(bool)

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
same_masks = bool((_o2["contaminated"] == contam).all() and (_o2["known"] == known).all()
                  and np.array_equal(_o2["N_up"], NUP))
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
print(f"CHECK  contaminated / known / N_up unchanged by the shift: {same_masks}")
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

print("\nWORST-CELL READING  delta at which the keep-quantile cell reaches +%.2f sigma\n"
      "       (= -balancing_offset(keep)['offset']; binding cell and coverage as the API "
      "reports them)" % A.z)
R.column("delta_api_mm", "the same delta from catchment_dod_balance.balancing_offset, mm")
R.column("binding_cell", "(row, col) of the cell that sets it")
R.column("cover", "evaluable_fraction: share of the routed tile the datum rests on")
rows = []
for q in (1.0, 0.999, 0.99, 0.975, 0.95, 0.90):
    bo = balancing_offset(dod, perror_nugget if A.envelope == "correlated" else perror_total,
                          props, valid, RES, z=A.z, keep=q,
                          corr_sill=float(vgm.sill) if A.envelope == "correlated" else 0.0,
                          corr_range=float(vgm.range_) if A.envelope == "correlated" else None,
                          exclude=flood if A.floodplain_mode == "constraint" else None, **kw)
    se, _ = boot_se(q, L_ref)
    rows.append([f"{q:.4f}", f"{1000*np.quantile(s, q):+.1f}", f"{-1000*bo['offset']:+.1f}",
                 f"{1000*se:.1f}", str(bo["binding_cell"]),
                 f"{bo['evaluable_fraction']:.3f}"])
R.table(["keep", "delta_mm", "delta_api_mm", "se_mm", "binding_cell", "cover"], rows)

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

# ---- what flag rate does the ERROR MODEL alone produce? ------------------------------
# The criterion "no cell is flagged" is blunder-bound, and any other cut is a number
# somebody has to justify. So measure it: draw the DoD error the fitted model says exists
# -- white part at each cell's own perror, correlated part as a spherical Gaussian random
# field at the fitted sill and range, plus one coherent draw on the unaccounted-area term
# exactly as the envelope adds it in quadrature -- route it through the SAME graph, and
# see what fraction of evaluable cells it flags with no real deposition present at all.
# That fraction is the flag rate the answer cannot go below, and it is measured, not set.
if A.null_real > 0:
    src_l, dst_l, frc_l = [], [], []
    Pm = np.asarray(props, float).reshape(ny * nx, 9)
    vf = valid.ravel()
    Ii = np.arange(ny * nx) // nx; Jj = np.arange(ny * nx) % nx
    from catchment_dod_balance.massbalance import _DX as _dx, _DY as _dy
    for k in range(1, 9):
        fk = Pm[:, k]; m = vf & (fk > 0)
        ni = Ii + _dy[k]; nj = Jj + _dx[k]
        ok = m & (ni >= 0) & (ni < ny) & (nj >= 0) & (nj < nx)
        dd = np.where(ok, ni * nx + nj, 0)
        ok = ok & vf[dd]
        src_l.append(np.where(ok)[0]); dst_l.append(dd[ok]); frc_l.append(fk[ok])
    src = np.concatenate(src_l); dst = np.concatenate(dst_l); frc = np.concatenate(frc_l)
    o = np.argsort(src, kind="stable"); src_s = src[o]; dst_s = dst[o]; frc_s = frc[o]
    starts = np.searchsorted(src_s, np.arange(ny * nx))
    ends = np.searchsorted(src_s, np.arange(ny * nx), side="right")
    indeg = np.bincount(dst, minlength=ny * nx)
    from collections import deque
    order = np.empty(int(vf.sum()), np.int64); kk = 0
    ind = indeg.copy()
    dq = deque(np.where(vf & (ind == 0))[0].tolist())
    while dq:
        c = dq.popleft(); order[kk] = c; kk += 1
        for e in range(starts[c], ends[c]):
            ind[dst_s[e]] -= 1
            if ind[dst_s[e]] == 0:
                dq.append(dst_s[e])
    assert kk == int(vf.sum()), "flow graph is not a DAG over the routed cells"

    def _acc(w):
        a = np.asarray(w, float).ravel().copy(); a[~vf] = 0.0
        for c in order:
            av = a[c]
            if av == 0.0:
                continue
            for e in range(starts[c], ends[c]):
                a[dst_s[e]] += frc_s[e] * av
        return a.reshape(ny, nx)

    # verify the cached graph against the shipped accumulator before trusting it
    from catchment_dod_balance.massbalance import weighted_accumulation as _wa
    _ref, _ = _wa(np.where(valid, 1.0, 0.0), props, valid)
    _mine = _acc(np.where(valid, 1.0, 0.0))
    print(f"\nCHECK  cached flow graph vs weighted_accumulation: max |diff| "
          f"{np.nanmax(np.abs(_mine - _ref)):.3e} cells")

    if A.envelope == "correlated":
        sill_n, range_n, pe = float(vgm.sill), float(vgm.range_), perror_nugget
    else:
        sill_n, range_n, pe = 0.0, None, perror_total
    M1 = int(2 ** np.ceil(np.log2(2 * ny))); M2 = int(2 ** np.ceil(np.log2(2 * nx)))
    if sill_n > 0:
        di = np.minimum(np.arange(M1), M1 - np.arange(M1)) * RES
        dj = np.minimum(np.arange(M2), M2 - np.arange(M2)) * RES
        hh = np.hypot(di[:, None], dj[None, :])
        Cov = np.where(hh < range_n,
                       sill_n * (1 - 1.5 * hh / range_n + 0.5 * (hh / range_n) ** 3), 0.0)
        lam = np.real(np.fft.fft2(Cov))
        neg = lam < 0
        clip_frac = float(-lam[neg].sum() / lam[~neg].sum()) if neg.any() else 0.0
        lam = np.maximum(lam, 0.0)
        print(f"CHECK  circulant embedding {M1}x{M2}, negative-eigenvalue mass clipped: "
              f"{clip_frac:.3e}")
    U_acc = out["unaccounted_area_acc"]; s_hole = out["sigma_hole"]
    grng = np.random.default_rng(A.seed + 991)
    fr = np.empty(A.null_real); Tn = np.empty(A.null_real)
    for b in range(A.null_real):
        e = np.where(known, grng.standard_normal(dod.shape) * pe, 0.0)
        if sill_n > 0:
            aa = grng.standard_normal((M1, M2)); bb = grng.standard_normal((M1, M2))
            g = np.real(np.fft.fft2((aa + 1j * bb) * np.sqrt(lam / (M1 * M2))))[:ny, :nx]
            e = e + np.where(known, g, 0.0)
        Vn = _acc(e * AREA) + grng.standard_normal() * s_hole * U_acc
        fr[b] = float(np.mean(Vn[ev] > A.z * SIG[ev]))
        Tn[b] = float(np.sum(np.maximum(Vn[ev] - A.z * SIG[ev], 0.0)))
    f_null = float(np.mean(fr)); f_sd = float(np.std(fr, ddof=1))
    d_null = 1000.0 * float(np.quantile(s, 1.0 - f_null))
    d_lo = 1000.0 * float(np.quantile(s, 1.0 - min(f_null + f_sd, 0.999)))
    d_hi = 1000.0 * float(np.quantile(s, 1.0 - max(f_null - f_sd, 1e-6)))
    se_q, _ = boot_se(1.0 - f_null, L_ref)
    print(f"\nNULL CALIBRATION  {A.null_real} routed realisations of the fitted error model, "
          "no real deposition")
    se_null = f_sd / np.sqrt(A.null_real)
    d_se_lo = 1000.0 * float(np.quantile(s, 1.0 - min(f_null + se_null, 0.999)))
    d_se_hi = 1000.0 * float(np.quantile(s, 1.0 - max(f_null - se_null, 1e-6)))
    R.column("stat", "name of the null-calibration quantity")
    R.column("val", "its value, units in the name")
    R.table(["stat", "val"], [
        ["flag rate produced by error alone, mean over realisations",
         f"{f_null:.4f}"],
        ["   its SD over realisations", f"{f_sd:.4f}"],
        ["   nominal one-sided rate implied by z=%.2f" % A.z,
         f"{0.5 * (1 - __import__("math").erf(A.z / np.sqrt(2))):.4f}"],
        ["delta at which the observed flag rate falls to it (mm)", f"{d_null:+.1f}"],
        ["   +/- from the SD of the null flag rate (mm)",
         f"{d_lo:+.1f} .. {d_hi:+.1f}"],
        ["   block-bootstrap SE of that quantile, L=%.0f m (mm)" % L_ref,
         f"{1000 * se_q:.1f}"],
        ["   +/- from the SE of the mean null rate (mm)",
         f"{d_se_lo:+.1f} .. {d_se_hi:+.1f}"],
        ["observed flag rate at delta = 0", f"{np.mean(s > 0):.4f}"],
    ])

    # Scene-integrated version of the same idea: the TOTAL excess deposition volume,
    # T(delta) = sum over evaluable cells of max(V_acc(delta) - z*sigma, 0), in m^3.
    # Monotone decreasing in delta, so one bisection finds where it meets the mean of
    # its own null distribution. Far less blunder-sensitive than the worst cell, because
    # a single spike contributes its own volume and nothing more.
    Vev = V[ev] - VIN[ev]; Sev = A.z * SIG[ev]; Nev = AREA * NUP[ev]

    def T_of(d):
        return float(np.sum(np.maximum(Vev - d * Nev - Sev, 0.0)))

    T_obs = T_of(0.0); T_null = float(np.mean(Tn)); T_sd = float(np.std(Tn, ddof=1))

    def solve_T(target):
        lo, hi = -5.0, 5.0
        if T_of(lo) < target:
            return np.nan
        if T_of(hi) > target:
            return np.nan
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if T_of(mid) > target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    d_T = solve_T(T_null)
    print("\nSCENE-INTEGRATED, NULL-CALIBRATED  total excess deposition volume "
          "T = sum max(V_acc - z*sigma, 0)")
    R.table(["stat", "val"], [
        ["T observed at delta = 0 (m^3)", f"{T_obs:+,.0f}"],
        ["T from the error model alone, mean over realisations (m^3)", f"{T_null:+,.0f}"],
        ["   its SD over realisations (m^3)", f"{T_sd:+,.0f}"],
        ["delta at which T falls to that mean (mm)", f"{1000*d_T:+.1f}"],
        ["   at mean - 1 SD of the null T (mm)",
         f"{1000*solve_T(max(T_null - T_sd, 0.0)):+.1f}"],
        ["   at mean + 1 SD of the null T (mm)", f"{1000*solve_T(T_null + T_sd):+.1f}"],
    ])

    # Does the residual still have structure once that raise is applied? A uniform raise
    # can only remove a uniform level; anything organised by drainage area or by canopy
    # cover that survives it is not a datum error.
    print(f"\nRESIDUAL STRUCTURE at delta = {d_null:+.1f} mm  (null rate {f_null:.4f})")
    R.column("stratum", "the cells the row is measured over")
    R.column("f_res", "fraction of them still flagged after the raise, dimensionless")
    rr = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        last = i == len(edges) - 2
        sel = ev & (carea_all >= lo) & ((carea_all <= hi) if last else (carea_all < hi))
        if sel.sum() < 50:
            continue
        rr.append([f"drainage {lo/1e4:.2f}-{hi/1e4:.2f} ha", str(int(sel.sum())),
                   f"{np.mean(S[sel] > d_null / 1000.0):.4f}"])
    if cover is not None:
        cov = cover
        cedges = [0.0, 0.05, 0.15, 0.35, 0.60, 1.01]
        for i in range(len(cedges) - 1):
            sel = ev & np.isfinite(cov) & (cov >= cedges[i]) & (cov < cedges[i + 1])
            if sel.sum() < 50:
                continue
            rr.append([f"cover {cedges[i]:.2f}-{cedges[i+1]:.2f}", str(int(sel.sum())),
                       f"{np.mean(S[sel] > d_null / 1000.0):.4f}"])
    R.table(["stratum", "n", "f_res"], rr)

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

# The aggregate budget against ITS OWN correlated 1-sigma envelope, by the same
# Rolstad/Hugonnet convention mass_balance uses internally: nugget over every cell,
# sill over N_eff = n*area/(pi*range^2).
from catchment_dod_balance.massbalance import correlated_variance  # noqa: E402
n_kn = int(kn.sum())
_pe = (perror_nugget if A.envelope == "correlated" else perror_total)[kn]
_nug = (AREA ** 2) * float(np.sum(_pe ** 2))
_var = correlated_variance(np.array([float(n_kn)]), np.array([_nug]), AREA,
                           float(vgm.sill) if A.envelope == "correlated" else 0.0,
                           float(vgm.range_) if A.envelope == "correlated" else None)
sig_tot = float(np.sqrt(_var[0]))
V_tot = float(np.sum(dod[kn]) * AREA)
d_scene = (V_tot - A.z * sig_tot) / (AREA * n_kn)
print(f"\nSCENE-INTEGRATED, ALGEBRAIC  the whole hillslope budget against its own "
      f"{A.z:g}-sigma envelope")
R.table(["quantity", "value"], [
    ["scene net volume, sum over evaluated hillslope cells of DoD*area, m^3",
     f"{V_tot:+,.0f}"],
    [f"its {A.envelope} {A.z:g}-sigma envelope, m^3", f"{A.z*sig_tot:,.0f}"],
    ["N_eff used (cells)",
     f"{max(n_kn*AREA/(np.pi*float(vgm.range_)**2), 1.0):.1f}" if A.envelope == "correlated"
     else f"{n_kn}"],
    [f"delta putting the scene net volume at +{A.z:g} sigma (mm)", f"{1000*d_scene:+.1f}"],
    ["   the envelope's own share of that delta (mm)",
     f"{1000*A.z*sig_tot/(AREA*n_kn):+.1f}"],
    ["   the mean-DoD share of that delta (mm)", f"{1000*V_tot/(AREA*n_kn):+.1f}"],
])

# How much of that rests on the variogram fit? Refit it under different lag windows and
# pair samples and carry each fit all the way through to the delta.
print("\n   variogram sensitivity: refit, then carried through to the same delta")
R.column("fit", "how the variogram was refitted")
R.column("nugget_sd_m", "sqrt(nugget) of that fit, m")
R.column("sill_sd_m", "sqrt(partial sill) of that fit, m")
R.column("range_m", "correlation length of that fit, m")
R.column("d_scene_mm", "delta putting the scene net volume at +z sigma under that fit, mm")
rows = []
for lab, mlag, sd, npair in [("as used", A.max_lag, 0, 600_000),
                             ("max_lag 400 m", 400.0, 0, 600_000),
                             ("max_lag 800 m", 800.0, 0, 600_000),
                             ("seed 1", A.max_lag, 1, 600_000),
                             ("2x pairs", A.max_lag, 0, 1_200_000)]:
    c2, g2, n2 = empirical_variogram(
        (xx[stable] * RES).astype(float), (yy[stable] * RES).astype(float), dod[stable],
        max_lag=mlag, n_lags=25, n_pairs=npair, estimator="dowd", seed=sd)
    v2 = fit_spherical(c2, g2, n2)
    pe2 = (lod / A.z) * np.sqrt(v2.nugget / v2.total_sill)
    var2 = correlated_variance(np.array([float(n_kn)]),
                               np.array([(AREA ** 2) * float(np.sum(pe2[kn] ** 2))]),
                               AREA, float(v2.sill), float(v2.range_))
    d2 = (V_tot - A.z * float(np.sqrt(var2[0]))) / (AREA * n_kn)
    rows.append([lab, f"{np.sqrt(v2.nugget):.4f}", f"{np.sqrt(v2.sill):.4f}",
                 f"{v2.range_:.0f}", f"{1000*d2:+.1f}"])
R.table(["fit", "nugget_sd_m", "sill_sd_m", "range_m", "d_scene_mm"], rows)

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
