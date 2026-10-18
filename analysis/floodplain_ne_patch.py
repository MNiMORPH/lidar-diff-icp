#!/usr/bin/env python3
"""What are the erosion patches that SURVIVE the vegetation correction? The NE patch, in detail.

Andy, 2026-10-18: "The correction is in the right direction. There are still significant
'erosional' patches, especially in the NE corner."

Everything here is read from the star-only product, the raw merged gen1 cloud, the decimated
gen2 full-tile cloud, and NAIP 2010. Numbers quoted in FLOODPLAIN_RETURN_PDF.md section J
are this script's stdout. Figures: floodplain_ne_patches.png (six-panel zoom),
floodplain_ne_intensity.png (gen1 intensity beside the corrected DoD).

    ./lidar-icp/bin/python analysis/floodplain_ne_patch.py
"""
import datetime
import json

import laspy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.ndimage import binary_dilation, distance_transform_edt, label

S = "data/derived/elbaext_delong_star_nosurf"
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3
D = np.load(f"{S}/floodplain_vegcorr_delta.npy") * 1e3
fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3
sp = np.load(f"{S}/gen1_spread_mm.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool); z2 = np.load(f"{S}/z_after.npy")
c = json.load(open(f"{S}/corrections.json")); X0, Y0, X1, Y1 = c["bounds"]; res = c["res_m"]; ny, nx = dod.shape
below = fp & np.isfinite(new) & (new < -lod)
lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); ks = np.argsort(sizes[1:])[::-1][:8] + 1
p1 = lab == ks[0]; flat = fp & (sl <= 2)

print("1. LARGEST remaining below-LoD patches in the CORRECTED DoD (whole floodplain mask)")
print(f"  {'rank':>4}{'cells':>7}{'E (m)':>10}{'N (m)':>11}{'flat%':>7}{'slope':>7}{'corrDoD':>8}{'rawDoD':>8}{'delta':>7}{'span':>6}{'penfail%':>9}{'z2':>7}")
for r, k in enumerate(ks):
    m = lab == k; iy, ix = np.nonzero(m)
    print(f"  {r+1:>4}{m.sum():>7,}{X0+ix.mean()*res:>10.0f}{Y0+iy.mean()*res:>11.0f}{100*(m&(sl<=2)).sum()/m.sum():>6.0f}%{np.median(sl[m]):>7.1f}"
          f"{np.median(new[m]):>8.0f}{np.median(dod[m]):>8.0f}{np.median(D[m]):>7.0f}{np.nanmedian(sp[m]):>6.0f}{100*pf[m].mean():>8.0f}%{np.median(z2[m]):>7.1f}")

t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid", "intensity"]); g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); it = t["intensity"].to_numpy()[g].astype(float)
print("\n2. gen1 CSF-ground return heights above the gen2 surface (mm): patch 1 vs the rest of the flat floodplain")
for lab_, m in (("patch 1 (NE)", np.isin(cell, np.flatnonzero(p1.ravel()))), ("rest of flat floodplain", np.isin(cell, np.flatnonzero((flat & ~p1).ravel())))):
    v = d[m]; print(f"  {lab_:<26} n {m.sum():>9,}  p5 {np.percentile(v,5):>5.0f}  p25 {np.percentile(v,25):>5.0f}  p50 {np.percentile(v,50):>5.0f}  p75 {np.percentile(v,75):>5.0f}  p95 {np.percentile(v,95):>5.0f}   share <= 0: {100*np.mean(v<=0):.0f}%   intensity med {np.median(it[m]):.0f}")

cube = np.load("data/derived/elbaext/nearground_cells_sn.npz"); cells = cube["cells"]; e = cube["edges"]; mid = 0.5 * (e[:-1] + e[1:]) * 1e3
g2 = np.load("data/derived/elbaext/nearground_gen2_class_split.npz"); Hall = (g2["Hg"] + g2["Hn"]).astype(float)
inc = np.isin(cells, np.flatnonzero(p1.ravel())); Hc = Hall[inc].sum(0); tot = Hc.sum(); band = (mid > 150) & (mid <= 2000)
print(f"\n3. gen2 2021 near-ground cube in patch 1: {inc.sum():,} of {p1.sum():,} cells covered; returns/cell {tot/inc.sum():.0f}; "
      f"lowveg share (0.15-2 m) {Hc[band].sum()/tot:.3f} (floodplain-wide median 0.18); class-2 share {g2['Hg'][inc].sum()/tot:.2f}")
for lab_, sel in (("patch 1", inc), ("all cube cells on flat floodplain", np.isin(cells, np.flatnonzero(flat.ravel())))):
    Hc = Hall[sel].sum(0); P = Hc / Hc.sum(); C = np.cumsum(P)
    print(f"   gen2 near-ground quantiles vs its own surface (mm), {lab_}: p10 {mid[np.searchsorted(C,.1)]:.0f} p50 {mid[np.searchsorted(C,.5)]:.0f} p90 {mid[np.searchsorted(C,.9)]:.0f}")

f = laspy.read("data/before/elbaext_gen1_merged.laz"); x = np.asarray(f.x); y = np.asarray(f.y); cl = np.asarray(f.classification)
it1 = np.asarray(f.intensity).astype(float); ps = np.asarray(f.point_source_id); gt = np.asarray(f.gps_time); sa = np.asarray(f.scan_angle_rank).astype(float)
ix = ((x - X0) / res).astype(int); iy = ((y - Y0) / res).astype(int); ing = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny); cell1 = np.where(ing, iy * nx + ix, -1)
ring = binary_dilation(p1, iterations=10) & ~p1 & flat
print("\n4. gen1 RAW cloud, all classes: patch 1 vs a ring of flat cells within 50 m. (gps_time here is GPS week seconds; it dates nothing, but equal values = same pass)")
for lab_, m in (("patch 1", ing & np.isin(cell1, np.flatnonzero(p1.ravel()))), ("ring", ing & np.isin(cell1, np.flatnonzero(ring.ravel())))):
    u, cnt = np.unique(ps[m], return_counts=True)
    for p_ in u:
        mm = m & (ps == p_)
        print(f"  {lab_:<8} line {p_}: n {mm.sum():>6,}  gps {gt[mm].min():.0f}-{gt[mm].max():.0f}  |scan| med {np.median(np.abs(sa[mm])):.0f} deg  intensity med {np.median(it1[mm]):.0f}  "
              f"cls2 {100*np.mean(cl[mm]==2):.0f}%  cls9 (water) {100*np.mean(cl[mm]==9):.1f}%  single-return {100*np.mean(np.asarray(f.number_of_returns)[mm]==1):.0f}%")
I = np.full(ny * nx, np.nan); gg = ing & (cl == 2)
med = pd.Series(it1[gg]).groupby(cell1[gg]).median(); I[med.index.to_numpy()] = med.to_numpy(); I = I.reshape(ny, nx)
dist = np.where(p1, -distance_transform_edt(p1), distance_transform_edt(~p1)) * res
print("\n5. PROFILE across patch 1's boundary (m; negative = inside): median raw DoD, corrected DoD, gen1 intensity, gen1 span")
for lo, hi in ((-60, -40), (-40, -20), (-20, -10), (-10, 0), (0, 10), (10, 20), (20, 40), (40, 60)):
    m = flat & (dist >= lo) & (dist < hi) & np.isfinite(dod)
    print(f"  {lo:>4}..{hi:<4} n {m.sum():>5,}  rawDoD {np.median(dod[m]):>6.0f}  corr {np.median(new[m]):>6.0f}  intens {np.nanmedian(I[m]):>4.0f}  span {np.nanmedian(sp[m]):>4.0f}")

nz = np.load("data/naip/naip2010_4m.npz"); ndvi = nz["ndvi"]; T = nz["transform"]
def naip_idx(E, N): return int((N - T[5]) / T[4]), int((E - T[2]) / T[0])
def ndvi_at(mask):
    jy, jx = np.nonzero(mask); out = []
    for e_, n_ in zip(X0 + (jx + 0.5) * res, Y0 + (jy + 0.5) * res):
        r, cc = naip_idx(e_, n_)
        if 0 <= r < ndvi.shape[0] and 0 <= cc < ndvi.shape[1]: out.append(ndvi[r, cc])
    return np.array(out)
a = ndvi_at(p1); b = ndvi_at(flat & ~p1)
print(f"\n6. NAIP 2010 NDVI: patch 1 median {np.nanmedian(a):.2f} (p10 {np.nanpercentile(a,10):.2f}, p90 {np.nanpercentile(a,90):.2f}); rest of flat floodplain median {np.nanmedian(b):.2f} (p10 {np.nanpercentile(b,10):.2f}, p90 {np.nanpercentile(b,90):.2f})")

acc = {k: np.zeros(256, int) for k in ks[:4]}
with laspy.open("data/after/3dep2021_fulltile.laz") as f2:
    for ch in f2.chunk_iterator(5_000_000):
        cx = np.asarray(ch.x); cy = np.asarray(ch.y); ccl = np.asarray(ch.classification)
        jx = ((cx - X0) / res).astype(int); jy = ((cy - Y0) / res).astype(int); ok = (jx >= 0) & (jx < nx) & (jy >= 0) & (jy < ny); cc = np.where(ok, jy * nx + jx, -1)
        for k in acc:
            m = ok & np.isin(cc, np.flatnonzero((lab == k).ravel())); acc[k] += np.bincount(ccl[m], minlength=256)[:256]
print("\n7. gen2 2021 vendor classes (decimated full-tile cloud) inside the four largest patches: 1 unclassified, 2 ground, 9 water")
for r, k in enumerate(acc):
    a_ = acc[k]; tt = a_.sum()
    print(f"  patch {r+1}: " + ("no points" if tt == 0 else "  ".join(f"cls{c_} {100*a_[c_]/tt:.1f}%" for c_ in np.flatnonzero(a_) if a_[c_] / tt > 0.002)))

r0, r1, c0, c1 = 560, 810, 600, 890; ext = [X0 + c0 * res, X0 + c1 * res, Y0 + r0 * res, Y0 + r1 * res]
rgbn = nz["rgbn"]; pr0, pc0 = naip_idx(ext[0], ext[3]); pr1, pc1 = naip_idx(ext[1], ext[2])
pr0, pr1 = max(min(pr0, pr1), 0), min(max(pr0, pr1), rgbn.shape[1]); pc0, pc1 = max(min(pc0, pc1), 0), min(max(pc0, pc1), rgbn.shape[2])
rgb = np.moveaxis(rgbn[:3, pr0:pr1, pc0:pc1], 0, -1).astype(float); rgb = np.clip(rgb / np.percentile(rgb, 99), 0, 1)
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
def show(ax, field, cm, vmin, vmax, title):
    im = ax.imshow(np.where(fp, field, np.nan)[r0:r1, c0:c1], cmap=cm, vmin=vmin, vmax=vmax, origin="lower", extent=ext, interpolation="nearest"); ax.set_title(title, fontsize=10); plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
show(axes[0, 0], dod, "RdBu", -300, 300, "raw DoD (mm)"); show(axes[0, 1], new, "RdBu", -300, 300, "corrected DoD (mm)"); show(axes[0, 2], z2, "viridis", 216, 226, "gen2 ground elevation (m)")
axes[1, 0].imshow(rgb, origin="upper", extent=ext); axes[1, 0].set_title("NAIP 2010 true colour", fontsize=10)
im = axes[1, 1].imshow(ndvi[pr0:pr1, pc0:pc1], cmap="YlGn", vmin=-0.2, vmax=0.9, origin="upper", extent=ext); axes[1, 1].set_title("NAIP 2010 NDVI", fontsize=10); plt.colorbar(im, ax=axes[1, 1], fraction=0.04, pad=0.02)
show(axes[1, 2], sl, "magma", 0, 10, "terrain slope (deg)")
yy, xx = np.mgrid[r0:r1, c0:c1]; Eg = X0 + (xx + 0.5) * res; Ng = Y0 + (yy + 0.5) * res
for ax in axes.ravel():
    for k in ks[:4]: ax.contour(Eg, Ng, (lab[r0:r1, c0:c1] == k).astype(float), levels=[0.5], colors="k", linewidths=0.8)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("NE corner of the Elba floodplain: the four largest remaining below-LoD patches outlined in black", fontsize=11)
fig.tight_layout(); fig.savefig("figures/floodplain_ne_patches.png", dpi=100); plt.close(fig)
fig, axes = plt.subplots(1, 2, figsize=(13, 6))
im = axes[0].imshow(np.where(fp, I, np.nan)[r0:r1, c0:c1], cmap="Greys_r", vmin=5, vmax=45, origin="lower", interpolation="nearest"); axes[0].set_title("gen1 2008 class-2 intensity (cell median), floodplain mask", fontsize=10); plt.colorbar(im, ax=axes[0], fraction=0.04)
im = axes[1].imshow(np.where(fp, new, np.nan)[r0:r1, c0:c1], cmap="RdBu", vmin=-300, vmax=300, origin="lower", interpolation="nearest"); axes[1].set_title("corrected DoD (mm)", fontsize=10); plt.colorbar(im, ax=axes[1], fraction=0.04)
for ax in axes:
    for k in ks[:4]: ax.contour((lab[r0:r1, c0:c1] == k).astype(float), levels=[0.5], colors="k" if ax is axes[1] else "r", linewidths=0.8)
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig("figures/floodplain_ne_intensity.png", dpi=100)
print("\nfigures: figures/floodplain_ne_patches.png, figures/floodplain_ne_intensity.png")
