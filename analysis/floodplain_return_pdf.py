#!/usr/bin/env python3
"""PRE-TEST: what does the PDF of gen1's floodplain returns look like, and does it change
with pooling width, with the cloth's point selection, or with beam geometry?

Andy, 2026-10-16: "Should we run a pre-test to see the statistical structure of the
floodplain returns? Perhaps over progressively larger areas, seeing what the PDF of returns
looks like? And/or something that uses first/last returns cleverly?"

WHY. After 3x3 pooling with the (fixed) half-sample mode, ~66 mm of DoD structure across
gen1-spread quartiles remains on the flat floodplain. It is not the estimator's count bias
(1.2 mm at 144 returns) and not the absence of returns at the gen2 surface (98% of cells
have one). So gen1 sees the ground but the ground is not its densest cluster. Whether a
LOWEST-cluster estimator could work depends on whether a ground peak exists at all -- that
is a question about the PDF, and it is answered before any estimator is built.

THREE TESTS, one population: flat floodplain cells (floodplain mask AND terrain slope <= 2
deg), stratified by gen1-spread quartile exactly as gen1_coarsened_ground.py does.

  A. PDF OVER GROWING AREAS. Per-cell histograms of the slope-normal residual against the
     gen2 surface (``d_mm_corr``, registration applied) pooled over 1x1, 3x3, 5x5, 9x9 and
     15x15 cells. Reported: the quartile-average PDF per scale (figure), and per cell the
     mode, median, share of returns at or below the gen2 surface, and whether a SECOND peak
     exists below the main one.
  B. CSF-SELECTED vs ALL RETURNS (``--all-table``): the same thing on a table built by the
     same two producers from the raw merged cloud. The CSF subset must reproduce
     ``d_mm_corr`` to the float32 bit before any comparison is read.
  C. BEAM GEOMETRY: the population PDF by |scan angle| bin and by intensity tercile.
     First/last returns are not usable: gen1 is 97.8% single-return.

PARAMETERS THAT ARE MINE, swept rather than chosen: peak-finding smoothing (20 and 40 mm)
and the second-peak prominence (10% and 20% of the main peak's height). Scan-angle bin
edges are 5 deg steps, the instrument's own 1 deg resolution times five.

    ./lidar-icp/bin/python analysis/floodplain_return_pdf.py
    ./lidar-icp/bin/python analysis/floodplain_return_pdf.py --all-table data/derived/elbaext_star_allreturns
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq
from scipy.ndimage import gaussian_filter1d, uniform_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hsm import hsm_from_hist   # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elbaext_delong_star_nosurf")
ap.add_argument("--all-table", default=None,
                help="scratch tile holding the ALL-returns beam_offset_table (test B)")
ap.add_argument("--scales", default="1,3,5,9,15")
ap.add_argument("--figdir", default="figures")
A = ap.parse_args()
S = A.tile
SCALES = [int(v) for v in A.scales.split(",")]
LO, HI, DZ = -1500.0, 1500.0, 20.0
edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
BWS = (20.0, 40.0)          # mm, peak-finding smoothing
PROMS = (0.10, 0.20)        # second peak must reach this fraction of the main peak
SA_EDGES = [0, 5, 10, 15, 90]

print("INPUTS AND ROLES")
print(f"  {S}/beam_offset_table.parquet   gen1 CSF ground returns: cell, d_mm_corr (slope-normal mm vs gen2, registration applied), scan_angle, intensity")
print(f"  {S}/floodplain_mask.npy, slope.npy, gen1_spread_mm.npy   population and its quartiles (same as gen1_coarsened_ground.py)")
if A.all_table:
    print(f"  {A.all_table}/beam_offset_table.parquet   ALL gen1 returns, same producers, raw merged cloud")
print("PARAMETERS")
print(f"  bins {DZ:g} mm over {LO:g}..{HI:g}; scales {SCALES} cells; smoothing {BWS} mm; prominence {PROMS}; |scan| bins {SA_EDGES}")
print("COLUMNS: mode = histogram half-sample mode (tie-symmetric); med = median; below = share of returns with d <= 0;")
print("  2nd peak = share of cells whose smoothed PDF has a local maximum BELOW the main one reaching the prominence;")
print("  sep = main peak minus that lower peak (mm), median over those cells")

sl = np.load(f"{S}/slope.npy"); fp = np.load(f"{S}/floodplain_mask.npy").astype(bool)
sp = np.load(f"{S}/gen1_spread_mm.npy"); ny, nx = sl.shape


def load_table(tile):
    t = pq.read_table(f"{tile}/beam_offset_table.parquet",
                      columns=["cell", "d_mm_corr", "in_grid", "scan_angle", "intensity",
                               "gps_time", "z"])
    g = t["in_grid"].to_numpy().astype(bool)
    out = {c: t[c].to_numpy()[g] for c in t.column_names if c != "in_grid"}
    d = out["d_mm_corr"].astype(float)
    ok = np.isfinite(d) & (d > LO) & (d < HI)
    return {c: v[ok] for c, v in out.items()}


def hist_grid(cell, d):
    H = np.zeros((ny * nx, NB), np.int32)
    np.add.at(H, (cell, np.clip(((d - LO) / DZ).astype(int), 0, NB - 1)), 1)
    return H.reshape(ny, nx, NB)


def pooled(H, k):
    if k == 1:
        return H.astype(np.float32)
    return np.rint(uniform_filter(H.astype(np.float32), size=(k, k, 1), mode="nearest") * k * k)


def peaks(h, bw, prom):
    """(n_peaks_meeting_prom, main_peak_mm, lowest_qualifying_peak_mm_or_nan)."""
    s = gaussian_filter1d(h, bw / DZ)
    if s.max() <= 0:
        return 0, np.nan, np.nan
    loc = np.flatnonzero((s[1:-1] > s[:-2]) & (s[1:-1] >= s[2:])) + 1
    loc = loc[s[loc] >= prom * s.max()]
    main = ctr[int(np.argmax(s))]
    low = ctr[loc.min()] if loc.size and ctr[loc.min()] < main else np.nan
    return int(loc.size), float(main), float(low)


def cell_stats(Hp, idx):
    """Per-cell mode, median, below-share and peak structure for the cells in idx."""
    flat = Hp.reshape(-1, NB)
    out = {"mode": [], "med": [], "below": [], "n": []}
    for bw in BWS:
        for pr in PROMS:
            out[f"two_{bw:g}_{pr:g}"] = []; out[f"sep_{bw:g}_{pr:g}"] = []
    for i in idx:
        h = flat[i]; tot = h.sum()
        if tot < 3:
            continue
        c = np.cumsum(h)
        out["n"].append(tot)
        out["mode"].append(hsm_from_hist(h, ctr))
        out["med"].append(ctr[np.searchsorted(c, 0.5 * tot)])
        out["below"].append(h[ctr <= 0].sum() / tot)
        for bw in BWS:
            for pr in PROMS:
                npk, main, low = peaks(h, bw, pr)
                two = np.isfinite(low)
                out[f"two_{bw:g}_{pr:g}"].append(two)
                out[f"sep_{bw:g}_{pr:g}"].append(main - low if two else np.nan)
    return {k: np.asarray(v, float) for k, v in out.items()}


def avg_pdf(Hp, idx):
    flat = Hp.reshape(-1, NB)[idx]
    tot = flat.sum(1, keepdims=True); ok = tot[:, 0] > 0
    return (flat[ok] / tot[ok]).mean(0) / DZ          # density per mm, each cell weighted equally


# ------------------------------------------------------------------ population
T = load_table(S)
H = hist_grid(T["cell"], T["d_mm_corr"])
n_own = H.sum(2)
flat = fp & (sl <= 2.0) & np.isfinite(sp) & (n_own >= 3)
qs = np.percentile(sp[flat], [25, 50, 75])
qlab = [f"Q1 <{qs[0]:.0f}", f"Q2 {qs[0]:.0f}-{qs[1]:.0f}", f"Q3 {qs[1]:.0f}-{qs[2]:.0f}", f"Q4 >{qs[2]:.0f}"]
qmask = [flat & (sp >= lo) & (sp < hi) for lo, hi in zip([0] + list(qs), list(qs) + [1e9])]
print(f"\nflat floodplain cells: {flat.sum():,}; gen1-spread quartile edges (mm): {qs.round(1)}")
print(f"  CSF returns in these cells: {int(n_own[flat].sum()):,}; per cell median {np.median(n_own[flat]):.0f}")

# ------------------------------------------------------------------ A: growing areas
print("\n" + "=" * 100)
print("A. PDF OVER GROWING AREAS -- CSF-selected returns (the pipeline's population)")
blues = plt.get_cmap("Blues")
fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=True)
hdr = f"  {'quartile':<14}{'scale':>6}{'n/cell':>8}{'mode':>8}{'med':>8}{'below':>7}" + \
      "".join(f"{f'2nd{bw:g}/{pr:g}':>11}{'sep':>6}" for bw in BWS for pr in PROMS)
print(hdr)
for k_i, k in enumerate(SCALES):
    Hp = pooled(H, k)
    for q, (lab, m) in enumerate(zip(qlab, qmask)):
        idx = np.flatnonzero(m.ravel())
        st = cell_stats(Hp, idx)
        row = f"  {lab:<14}{f'{k}x{k}':>6}{np.median(st['n']):>8.0f}{np.median(st['mode']):>8.1f}" \
              f"{np.median(st['med']):>8.1f}{np.median(st['below']):>7.2f}"
        for bw in BWS:
            for pr in PROMS:
                tw = st[f"two_{bw:g}_{pr:g}"]; se = st[f"sep_{bw:g}_{pr:g}"]
                row += f"{100*tw.mean():>10.1f}%{np.nanmedian(se) if tw.any() else np.nan:>6.0f}"
        print(row, flush=True)
        pdf = avg_pdf(Hp, idx)
        axes[q].plot(ctr, pdf, color=blues(0.35 + 0.6 * k_i / (len(SCALES) - 1)), lw=1.5,
                     label=f"{k}x{k} ({k*5} m)")
    del Hp
for q, lab in enumerate(qlab):
    ax = axes[q]; ax.set_title(f"gen1 spread {lab} mm", fontsize=10)
    ax.axvline(0, color="#52514e", lw=0.8, ls=":"); ax.set_xlim(-500, 500)
    ax.set_xlabel("gen1 return height above gen2 surface (mm)")
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
axes[0].set_ylabel("density per mm (cells weighted equally)"); axes[0].legend(frameon=False, fontsize=8, title="pooling")
fig.suptitle("A. Flat floodplain: PDF of gen1 CSF returns vs pooling width, by gen1-spread quartile", fontsize=11)
fig.tight_layout(); fig.savefig(f"{A.figdir}/floodplain_return_pdf_A_scales.png", dpi=130); plt.close(fig)
print(f"  figure: {A.figdir}/floodplain_return_pdf_A_scales.png")

# ------------------------------------------------------------------ C: beam geometry
print("\n" + "=" * 100)
print("C. BEAM GEOMETRY -- population PDF of CSF returns in each quartile, by |scan angle| and by intensity tercile")
cellq = np.zeros(ny * nx, np.int8) - 1
for q, m in enumerate(qmask):
    cellq[m.ravel()] = q
rq = cellq[T["cell"]]
sa = np.abs(T["scan_angle"].astype(float)); it = T["intensity"].astype(float)
inq = rq >= 0
iterc = np.percentile(it[inq], [100 / 3, 200 / 3])
print(f"  intensity terciles over flat-floodplain returns: {iterc.round(0)}")


def pop_stats(d):
    h = np.histogram(d, bins=edges)[0]
    tot = h.sum(); c = np.cumsum(h)
    return hsm_from_hist(h, ctr), ctr[np.searchsorted(c, 0.5 * tot)], h[ctr <= 0].sum() / tot, h / tot / DZ


fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharey=True)
print(f"  {'quartile':<14}{'bin':<16}{'n':>10}{'mode':>8}{'med':>8}{'below':>7}")
for q, lab in enumerate(qlab):
    base = inq & (rq == q)
    bins = [(f"|scan| {a}-{b}", base & (sa >= a) & (sa < b)) for a, b in zip(SA_EDGES[:-1], SA_EDGES[1:])]
    bins += [(f"intens {['low','mid','high'][j]}", base & cond) for j, cond in enumerate(
        (it < iterc[0], (it >= iterc[0]) & (it < iterc[1]), it >= iterc[1]))]
    for j, (bl, bm) in enumerate(bins):
        if bm.sum() < 200:
            continue
        mo, me, be, pdf = pop_stats(T["d_mm_corr"][bm])
        print(f"  {lab:<14}{bl:<16}{int(bm.sum()):>10,}{mo:>8.1f}{me:>8.1f}{be:>7.2f}")
        ax = axes[0 if j < 4 else 1, q]
        ramp = blues(0.35 + 0.6 * (j if j < 4 else j - 4) / 3)
        ax.plot(ctr, pdf, color=ramp, lw=1.5, label=bl)
for q, lab in enumerate(qlab):
    for r in range(2):
        ax = axes[r, q]; ax.axvline(0, color="#52514e", lw=0.8, ls=":"); ax.set_xlim(-500, 500)
        for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    axes[0, q].set_title(f"gen1 spread {lab} mm", fontsize=10)
    axes[1, q].set_xlabel("gen1 return height above gen2 surface (mm)")
axes[0, 0].legend(frameon=False, fontsize=8); axes[1, 0].legend(frameon=False, fontsize=8)
axes[0, 0].set_ylabel("density per mm, by |scan angle|"); axes[1, 0].set_ylabel("density per mm, by intensity tercile")
fig.suptitle("C. Flat floodplain: PDF of gen1 CSF returns by beam geometry", fontsize=11)
fig.tight_layout(); fig.savefig(f"{A.figdir}/floodplain_return_pdf_C_geometry.png", dpi=130); plt.close(fig)
print(f"  figure: {A.figdir}/floodplain_return_pdf_C_geometry.png")

# ------------------------------------------------------------------ B: CSF vs all returns
if A.all_table:
    print("\n" + "=" * 100)
    print("B. CSF-SELECTED vs ALL gen1 RETURNS, same cells, same producers")
    TA = load_table(A.all_table)
    # identity check on the CSF subset: key on (gps_time, z) which the same LAS supplies to both
    key = lambda t: np.round(t["gps_time"] * 1e6).astype(np.int64) * 100003 + np.round(t["z"] * 100).astype(np.int64)
    ka, kc = key(TA), key(T)
    order = np.argsort(ka); pos = np.searchsorted(ka, kc, sorter=order)
    pos = np.clip(pos, 0, ka.size - 1); hit = ka[order[pos]] == kc
    diff = np.abs(TA["d_mm_corr"][order[pos[hit]]].astype(float) - T["d_mm_corr"][hit].astype(float))
    print(f"  CSF returns matched in the all-returns table: {hit.mean()*100:.3f}%  "
          f"max |d_mm_corr difference| {diff.max():.3e} mm  (must be 0)")
    HA = hist_grid(TA["cell"], TA["d_mm_corr"])
    nA = HA.sum(2)
    print(f"  returns per flat cell: CSF {np.median(n_own[flat]):.0f}   all {np.median(nA[flat]):.0f}")
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharey="row")
    print(f"  {'quartile':<14}{'scale':>6}{'table':>6}{'n/cell':>8}{'mode':>8}{'med':>8}{'below':>7}"
          + "".join(f"{f'2nd{bw:g}/{pr:g}':>11}{'sep':>6}" for bw in BWS for pr in PROMS))
    for r, k in enumerate((1, 3)):
        Hp_c, Hp_a = pooled(H, k), pooled(HA, k)
        for q, (lab, m) in enumerate(zip(qlab, qmask)):
            idx = np.flatnonzero(m.ravel())
            for nm, Hp, col in (("CSF", Hp_c, "#2a78d6"), ("all", Hp_a, "#eb6834")):
                st = cell_stats(Hp, idx)
                row = f"  {lab:<14}{f'{k}x{k}':>6}{nm:>6}{np.median(st['n']):>8.0f}{np.median(st['mode']):>8.1f}" \
                      f"{np.median(st['med']):>8.1f}{np.median(st['below']):>7.2f}"
                for bw in BWS:
                    for pr in PROMS:
                        tw = st[f"two_{bw:g}_{pr:g}"]; se = st[f"sep_{bw:g}_{pr:g}"]
                        row += f"{100*tw.mean():>10.1f}%{np.nanmedian(se) if tw.any() else np.nan:>6.0f}"
                print(row, flush=True)
                axes[r, q].plot(ctr, avg_pdf(Hp, idx), color=col, lw=1.5, label=f"{nm} returns")
            ax = axes[r, q]; ax.axvline(0, color="#52514e", lw=0.8, ls=":"); ax.set_xlim(-500, 500)
            for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
            if r == 0: ax.set_title(f"gen1 spread {lab} mm", fontsize=10)
            if r == 1: ax.set_xlabel("gen1 return height above gen2 surface (mm)")
        axes[r, 0].set_ylabel(f"density per mm, {k}x{k} pooling")
        del Hp_c, Hp_a
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("B. Flat floodplain: CSF-selected vs all gen1 returns", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{A.figdir}/floodplain_return_pdf_B_csf_vs_all.png", dpi=130); plt.close(fig)
    print(f"  figure: {A.figdir}/floodplain_return_pdf_B_csf_vs_all.png")
