#!/usr/bin/env python3
"""Does the gen1-minus-gen2 offset on non-eroding ground come from PHENOLOGY or from DUFF?

On divide cells where erosion cannot operate, the offset d = (2008 ground) - (2021 ground)
grows increasingly NEGATIVE with canopy cover: the 2021 surface sits HIGHER under forest.
Two explanations survive:

    (A) PHENOLOGY -- 2008 was flown leaf-off (Nov) and reached the true forest floor, while
        2021 was flown leaf-on at green-up (May) and stopped short of it.
    (B) DUFF -- 13 years of litter and organic matter really did build the floor up.

Deciduous stands have a huge leaf-off/leaf-on contrast; CONIFERS have almost none. So the
two explanations make opposite predictions ACROSS FOREST TYPE at the same canopy cover:
under (A) the cover effect should be far stronger in deciduous stands; under (B) conifers
should show it as strongly or more, since conifer litter is the more persistent.

FOREST TYPE IS SEPARATED FROM EACH EPOCH'S OWN RETURNS, on an IDENTICAL definition for both
epochs (epoch_canopy_frac.py: fraction of that epoch's returns >2 m above the gen2 bare
earth). Using a 2021-derived covariate to describe the 2008 canopy would build the answer
into the question, and mixing a return fraction with a PyForestScan plant-area product
would confound forest type with method.

    CONIFER   ~ high canopy in BOTH epochs        (f_2008 >= --f1-conifer)
    DECIDUOUS ~ high in 2021, near-none in 2008   (f_2008 <= --f1-decid, f_2021 >= --f2-min)

The joint (f_2008, f_2021) distribution is printed and plotted first, so the reader can
judge whether two populations are actually separable here before reading any offset number.

Errors are spatially cluster-robust (binstats.block_ids, 50 m blocks): cells inside one
woodlot are not independent. Sparse cover bins are reported with honest large error bars,
never dropped.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/conifer_vs_deciduous.py --tile data/derived/elbaext
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from lidar_diff_icp.binstats import block_ids, quantile_edges, binned_stats, nmad

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext")
ap.add_argument("--curv-max", type=float, default=0.015)
ap.add_argument("--inc-max", type=float, default=5.0)
ap.add_argument("--f1-conifer", type=float, default=0.35, help="2008 canopy fraction at/above which a cell is called conifer")
ap.add_argument("--f1-decid", type=float, default=0.10, help="2008 canopy fraction at/below which a cell is called deciduous")
ap.add_argument("--f2-min", type=float, default=0.30, help="2021 canopy fraction a deciduous cell must exceed (it must HAVE a canopy)")
ap.add_argument("--min-ret", type=int, default=10, help="minimum returns per cell per epoch for a valid canopy fraction")
ap.add_argument("--block", type=float, default=50.0, help="spatial block size (m) for cluster-robust errors")
ap.add_argument("--nbins", type=int, default=4)
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))


def grid(tile):
    for fn in ("meta.json", "corrections_geoid.json", "corrections.json"):
        p = f"{tile}/{fn}"
        if os.path.exists(p):
            j = json.load(open(p)); b = j["bounds"]; r = float(j.get("res") or j.get("res_m"))
            ny = int(j.get("ny") or round((b[3] - b[1]) / r))
            nx = int(j.get("nx") or round((b[2] - b[0]) / r))
            return r, ny, nx
    raise SystemExit(f"no grid meta in {tile}")


RES, NY, NX = grid(A.tile)
g1 = np.load(f"{A.tile}/gen1_canopy_frac.npz")
g2 = np.load(f"{A.tile}/gen2_canopy_frac.npz")
f1 = g1["frac"].ravel(); f2 = g2["frac"].ravel()
n1 = g1["n_ret"].ravel(); n2 = g2["n_ret"].ravel()
bld1 = g1["n_bldg"].ravel()                     # gen1 ASPRS class 6: structures, not trees
pfs = np.load(f"{A.tile}/canopy_cover_pfs.npy").ravel()
HAG = float(g1["hag"])

# ---------------------------------------------------------------- reference population
df = pd.read_parquet(f"{A.tile}/beam_offset_table.parquet",
                     columns=["cell", "d_mm_corr", "curv_laplacian", "incidence", "in_grid",
                              "intensity"])
cell = df.cell.to_numpy()
ridge = np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()
keep = (df.in_grid.values
        & (df.curv_laplacian.abs().to_numpy() <= A.curv_max)
        & (df.incidence.to_numpy() < A.inc_max)
        & ridge[cell]
        & np.isfinite(df.d_mm_corr.to_numpy()))
cell = cell[keep]; d = df.d_mm_corr.to_numpy(float)[keep]
inten = df.intensity.to_numpy(float)[keep]
del df

# per-cell median offset -- the unit of analysis
order = np.argsort(cell, kind="stable")
cell, d, inten = cell[order], d[order], inten[order]
ucell, start = np.unique(cell, return_index=True)
dmed = np.array([np.median(d[a:b]) for a, b in zip(start, np.r_[start[1:], cell.size])])
nret = np.diff(np.r_[start, cell.size])

print("=" * 92)
print(f"CONIFER vs DECIDUOUS test of the cover-dependent offset   [{TILE}]")
print(f"reference = divide network, |curv_laplacian| <= {A.curv_max:g}, incidence < {A.inc_max:g} deg")
print(f"offset = d_mm_corr (mm, gen1 - gen2; negative = 2021 surface higher)")
print(f"{keep.sum():,} returns in {ucell.size:,} cells (median {np.median(nret):.0f} returns/cell)")
print("=" * 92)

valid = (n1[ucell] >= A.min_ret) & (n2[ucell] >= A.min_ret) & np.isfinite(f1[ucell]) & np.isfinite(f2[ucell])
nobld = bld1[ucell] == 0
print(f"valid canopy fraction in both epochs: {valid.sum():,} cells "
      f"({100*valid.mean():.1f}%);  gen1 class-6 (building) returns present in "
      f"{int((~nobld & valid).sum()):,} of those -- excluded")
sel = valid & nobld
C1 = f1[ucell]; C2 = f2[ucell]; P = pfs[ucell]

# ---------------------------------------------------------------- joint distribution
E1 = np.array([0, .05, .10, .20, .35, .50, .75, 1.001])
E2 = np.array([0, .15, .30, .45, .60, .75, 1.001])
print("\nJOINT (2008, 2021) CANOPY DISTRIBUTION -- reference cells, count "
      f"(fraction of returns >{HAG:g} m above gen2 bare earth, IDENTICAL definition both epochs)")
hdr = "  f2008 \\ f2021  " + "".join(f"{lo:.2f}-{hi:<5.2f}".rjust(12) for lo, hi in zip(E2[:-1], E2[1:])) + "      row"
print(hdr); print("-" * len(hdr))
for lo, hi in zip(E1[:-1], E1[1:]):
    r = sel & (C1 >= lo) & (C1 < hi)
    row = [int((r & (C2 >= a) & (C2 < b)).sum()) for a, b in zip(E2[:-1], E2[1:])]
    print(f"  {lo:.2f}-{hi:<5.2f}   " + "".join(f"{v:>12,}" for v in row) + f"{int(r.sum()):>9,}")
print("  " + "col".rjust(13) + "".join(f"{int((sel & (C2>=a) & (C2<b)).sum()):>12,}"
                                       for a, b in zip(E2[:-1], E2[1:])) + f"{int(sel.sum()):>9,}")

# threshold sensitivity: how many cells (and independent blocks) does each candidate leave?
print(f"\nCONIFER-THRESHOLD SENSITIVITY (cells / 50 m blocks surviving, f2021 >= {A.f2_min:g})")
for t in (0.20, 0.25, 0.30, 0.35, 0.40, 0.50):
    s = sel & (C1 >= t) & (C2 >= A.f2_min)
    nb = np.unique(block_ids(ucell[s], NX, RES, A.block)).size if s.any() else 0
    print(f"   f2008 >= {t:.2f}: {int(s.sum()):>6,} cells, {nb:>4d} blocks, "
          f"median f2021 {np.median(C2[s]) if s.any() else np.nan:.2f}")

conifer = sel & (C1 >= A.f1_conifer) & (C2 >= A.f2_min)
decid = sel & (C1 <= A.f1_decid) & (C2 >= A.f2_min)
mid = sel & (C2 >= A.f2_min) & ~conifer & ~decid
print(f"\nCLASSES (both require f2021 >= {A.f2_min:g}, i.e. the cell HAS a 2021 canopy):")
for nm, s in (("CONIFER  f2008 >= %.2f" % A.f1_conifer, conifer),
              ("DECIDUOUS f2008 <= %.2f" % A.f1_decid, decid),
              ("intermediate", mid)):
    nb = np.unique(block_ids(ucell[s], NX, RES, A.block)).size if s.any() else 0
    print(f"   {nm:<26s} {int(s.sum()):>7,} cells  {nb:>5d} blocks  "
          f"median f2021 {np.median(C2[s]):.3f}  median PFS cover {np.median(P[s]):.3f}")

# ---------------------------------------------------------------- offset vs cover by class
# Bin edges come from the CONIFER cover distribution (the limiting class, so each bin carries
# comparable conifer weight) and are then applied UNCHANGED to the other classes: a matched
# comparison needs identical bins, not bins convenient to each class.
PAI = np.load(f"{A.tile}/pai_pfs.npy").ravel()[ucell]

# FULL RANGE first, each class over its own whole cover span (no truncation): the matched
# tables below use conifer-derived edges, which necessarily clip the deciduous low end.
print(f"\n{'='*92}\nFULL COVER RANGE, each class binned over ITS OWN span (nothing truncated, "
      "nothing dropped for sparsity)")
for nm, s_ in (("CONIFER", conifer), ("DECIDUOUS", decid), ("intermediate", mid)):
    e = quantile_edges(P[s_], A.nbins + 2)
    b = binned_stats(P[s_], dmed[s_], e, block=block_ids(ucell[s_], NX, RES, A.block))
    print(f"\n  {nm}  (n = {int(s_.sum()):,} cells)")
    print(f"  {'cover bin':>15s} {'mean cover':>11s} {'median d (mm)':>14s} {'SE_cluster':>11s}"
          f" {'cells':>8s} {'blocks':>7s}")
    for i in range(len(b)):
        print(f"  {b.lo[i]:6.3f}-{b.hi[i]:<8.3f} {b.x[i]:>11.3f} {b.y[i]:>14.1f} "
              f"{b.se_block[i]:>11.2f} {b.n[i]:>8,} {b.n_block[i]:>7d}")

AXES = [("2021 canopy cover (PyForestScan)", P),
        ("2021 canopy RETURN FRACTION (same definition as the 2008 axis)", C2),
        ("2021 plant area index (PyForestScan; does NOT saturate like cover)", PAI)]
res = {}
for ax_name, X in AXES:
    edges = quantile_edges(X[conifer], A.nbins)
    print(f"\n{'='*92}\nOFFSET vs {ax_name.upper()}, BY FOREST TYPE")
    print("  bin edges (conifer quantiles, applied to every class): "
          + ", ".join(f"{e:.3f}" for e in edges))
    r = {}
    for nm, s_ in (("CONIFER", conifer), ("DECIDUOUS", decid), ("intermediate", mid)):
        b = binned_stats(X[s_], dmed[s_], edges, block=block_ids(ucell[s_], NX, RES, A.block))
        r[nm] = b
        print(f"\n  {nm}  (n = {int(s_.sum()):,} cells, median {np.median(nret[s_]):.0f} returns/cell)")
        print(f"  {'bin':>15s} {'mean x':>9s} {'median d (mm)':>14s} {'SE_cluster':>11s}"
              f" {'SE_naive':>9s} {'cells':>8s} {'blocks':>7s}")
        for i in range(len(b)):
            print(f"  {b.lo[i]:6.3f}-{b.hi[i]:<8.3f} {b.x[i]:>9.3f} {b.y[i]:>14.1f} "
                  f"{b.se_block[i]:>11.2f} {b.se_return[i]:>9.2f} {b.n[i]:>8,} {b.n_block[i]:>7d}")
    bc, bd = r["CONIFER"], r["DECIDUOUS"]
    print("\n  MATCHED CONTRAST (deciduous - conifer; PHENOLOGY predicts deciduous MORE "
          "negative, i.e. diff < 0)")
    print(f"  {'bin':>15s} {'mean x con':>11s} {'mean x dec':>11s} {'conifer':>9s} {'decid':>9s}"
          f" {'diff':>8s} {'SE':>7s} {'z':>7s}")
    dsum = wsum = 0.0
    for i in range(len(bd)):
        j = np.where(np.isclose(bc.lo, bd.lo[i]))[0]
        if not j.size:
            continue
        j = j[0]
        diff = bd.y[i] - bc.y[j]; se = float(np.hypot(bd.se[i], bc.se[j]))
        dsum += diff / se**2; wsum += 1 / se**2
        print(f"  {bd.lo[i]:6.3f}-{bd.hi[i]:<8.3f} {bc.x[j]:>11.3f} {bd.x[i]:>11.3f} "
              f"{bc.y[j]:>9.1f} {bd.y[i]:>9.1f} {diff:>8.1f} {se:>7.1f} {diff/se:>7.2f}")
    if wsum:
        print(f"  {'POOLED':>15s} {'':>11s} {'':>11s} {'':>9s} {'':>9s} "
              f"{dsum/wsum:>8.1f} {wsum**-0.5:>7.1f} {dsum/wsum*wsum**0.5:>7.2f}")
    res[ax_name] = r
res = res[AXES[0][0]]

# --------------------------------------------------- confounds the contrast has to survive
# (1) RETURN COUNT. Evergreen canopy passes fewer 2008 pulses to the floor, so conifer cells
#     rest on fewer returns than deciduous ones; if a thin per-cell sample biased the median
#     the contrast would be an artefact. Repeat it inside one common return-count band.
# (2) SPATIAL CONCENTRATION. If the conifer class were two woodlots it could carry one local
#     registration residual instead of a forest-type effect.
from scipy.ndimage import label as _label

X = P; lo, hi = np.quantile(X[conifer], [0.0, 0.25])
print(f"\n{'='*92}\nCONFOUND 1: same return-count band, cover-matched ({lo:.3f}-{hi:.3f})")
band = (nret >= 4) & (nret < 12)
print(f"  returns/cell, conifer p10/50/90: "
      f"{np.percentile(nret[conifer], [10,50,90]).astype(int)};  deciduous: "
      f"{np.percentile(nret[decid], [10,50,90]).astype(int)}")
vals = {}
for nm, s_ in (("CONIFER", conifer), ("DECIDUOUS", decid)):
    m_ = s_ & (X >= lo) & (X < hi) & band
    b = binned_stats(np.zeros(int(m_.sum())), dmed[m_], np.array([-1.0, 1.0]),
                     block=block_ids(ucell[m_], NX, RES, A.block))
    vals[nm] = (b.y[0], b.se[0])
    print(f"  {nm:<10s} {int(m_.sum()):>5,} cells, {b.n_block[0]:>4d} blocks, "
          f"mean cover {X[m_].mean():.3f}, median d {b.y[0]:>7.1f} +- {b.se[0]:.1f} mm")
_d = vals["DECIDUOUS"][0] - vals["CONIFER"][0]; _s = float(np.hypot(vals["DECIDUOUS"][1], vals["CONIFER"][1]))
print(f"  deciduous - conifer = {_d:+.1f} +- {_s:.1f} mm  (z = {_d/_s:.2f})")

msk = np.zeros(NY * NX, bool); msk[ucell[conifer]] = True
lab, nl = _label(msk.reshape(NY, NX), structure=np.ones((3, 3)))
sz = np.bincount(lab.ravel())[1:]
iy, ix = np.divmod(ucell[conifer], NX)
print(f"\nCONFOUND 2: the conifer class is {int(conifer.sum()):,} cells in {nl} connected patches "
      f"(largest {np.sort(sz)[::-1][:5]}),\n  spanning easting {ix.min()*RES:.0f}-{ix.max()*RES:.0f} m "
      f"and northing {iy.min()*RES:.0f}-{iy.max()*RES:.0f} m of the tile.")
big = np.argsort(sz)[::-1][:3] + 1
for L in big:
    sub = conifer & np.isin(ucell, np.flatnonzero(lab.ravel() == L))
    print(f"  drop patch of {int(sub.sum()):>3d} cells (median d {np.median(dmed[sub]):>7.1f}): "
          f"remaining conifer median d {np.median(dmed[conifer & ~sub]):>7.1f} mm")

# (3) RETURN INTENSITY. A weak 2008 return can be ranged LATE (timewalk), which reads the
#     ground LOW -- the same sign as duff. Evergreen canopy attenuates every 2008 pulse that
#     reaches the floor, so forest type and return strength are confounded by construction.
#     This is the one confound the test cannot design away, so it is measured per return.
per_cell_class = np.zeros(NY * NX, np.int8)
per_cell_class[ucell[conifer]] = 1; per_cell_class[ucell[decid]] = 2
per_cell_cov = np.full(NY * NX, np.nan); per_cell_cov[ucell] = P
kc = per_cell_class[cell]; band_r = (per_cell_cov[cell] >= lo) & (per_cell_cov[cell] < hi)
print(f"\n{'='*92}\nCONFOUND 3: per-return gen1 intensity, cover-matched ({lo:.3f}-{hi:.3f})")
for nm, kk in (("CONIFER", 1), ("DECIDUOUS", 2)):
    m_ = band_r & (kc == kk)
    print(f"  {nm:<10s} {int(m_.sum()):>7,} returns, intensity p10/50/90 = "
          f"{np.percentile(inten[m_], [10, 50, 90]).astype(int)}")
print("  matched-intensity contrast (deciduous - conifer):")
for a, b_ in ((0, 15), (15, 25), (25, 40), (40, 300)):
    out = []
    for kk in (1, 2):
        m_ = band_r & (kc == kk) & (inten >= a) & (inten < b_)
        if m_.sum() < 50:
            out.append((np.nan, np.nan, int(m_.sum()))); continue
        bb = binned_stats(np.zeros(int(m_.sum())), d[m_], np.array([-1.0, 1.0]),
                          block=block_ids(cell[m_], NX, RES, A.block))
        out.append((bb.y[0], bb.se[0], int(m_.sum())))
    (yc, sc, nc), (yd, sd, nd) = out
    dd = yd - yc; ss = float(np.hypot(sd, sc))
    print(f"   I {a:>3d}-{b_:<3d}: conifer {yc:>7.1f} +- {sc:<5.1f} (n={nc:>6,})  "
          f"deciduous {yd:>7.1f} +- {sd:<5.1f} (n={nd:>7,})  diff {dd:>7.1f} +- {ss:>5.1f}"
          f"  z = {dd/ss:.2f}")

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
h = ax[0].hexbin(C2[sel], C1[sel], gridsize=45, bins="log", cmap="viridis", extent=(0, 1, 0, 1))
ax[0].axhline(A.f1_conifer, color="r", lw=1); ax[0].axhline(A.f1_decid, color="c", lw=1)
ax[0].axvline(A.f2_min, color="k", lw=.8, ls=":")
ax[0].set_xlabel("2021 canopy return fraction"); ax[0].set_ylabel("2008 canopy return fraction")
ax[0].set_title(f"joint canopy, reference cells (n={int(sel.sum()):,})")
plt.colorbar(h, ax=ax[0], label="cells")
for lo, hi, c in ((0.45, 0.60, "tab:blue"), (0.60, 1.01, "tab:orange")):
    s = sel & (C2 >= lo) & (C2 < hi)
    ax[1].hist(C1[s], bins=np.arange(0, 1.02, .04), histtype="step", density=True,
               color=c, label=f"2021 cover {lo:.2f}-{hi:.2f}  (n={int(s.sum()):,})")
ax[1].axvline(A.f1_conifer, color="r", lw=1); ax[1].axvline(A.f1_decid, color="c", lw=1)
ax[1].set_xlabel("2008 canopy return fraction"); ax[1].set_ylabel("density")
ax[1].set_title("is the 2008 canopy bimodal at fixed 2021 cover?"); ax[1].legend(fontsize=7)
for nm, c in (("CONIFER", "tab:green"), ("DECIDUOUS", "tab:red"), ("intermediate", "0.6")):
    b = res[nm]
    ax[2].errorbar(b.x, b.y, yerr=b.se, marker="o", ms=4, capsize=3, color=c,
                   label=f"{nm} (n={int(b.n.sum()):,})")
ax[2].axhline(0, color="k", lw=.7)
ax[2].set_xlabel("2021 canopy cover (PyForestScan)"); ax[2].set_ylabel("median offset d (mm)")
ax[2].set_title("offset vs cover by forest type\n(cluster-robust SE, 50 m blocks)")
ax[2].legend(fontsize=7)
fig.suptitle(f"{TILE}: phenology vs duff -- offset by forest type on non-eroding divides", y=1.0)
fig.tight_layout()
out = f"analysis/ridgelines/conifer_vs_deciduous_{TILE}.png"
fig.savefig(out, dpi=140, bbox_inches="tight")
print(f"\nfigure: {out}")
