#!/usr/bin/env python3
"""2025 vs 2008: has the NE patch's October reflectance changed relative to the floodplain?

Andy, 2026-10-19: "Can you check 2025 vs. 2008 reflectance and if the grass changed?"

Like with like: Landsat-5 TM 2008-10-29 against Landsat-8/9 OLI 2025 (2025-10-04 L9, 13%
cloud catalogued; 2025-11-05 L9, 0%; 2025-09-26 L8, 0.3%) at 30 m, and Sentinel-2
2025-10-02 and 2025-10-27 (both 0%) at 10 m beside 2023-10-23 and 2024-10-17. Absolute
levels are NOT comparable across sensors or years (TM vs OLI band passes, atmosphere,
phenology of the date), so the comparison is GAUGE-FREE: each population's median
reflectance is divided by the flat floodplain's median in the same scene and band. A
patch whose relative spectrum is the same in 2008 and 2025 has not changed relative to
its surroundings; one whose relative NIR or SWIR moved has. The known classes (cropland,
grassland, tree) are run the same way as the null -- their relative spectra should be
stable across years if the method is.

Per pixel, for the NE patch and the F2 class: the log-ratio (2025 / 2008) of each band,
minus the floodplain's log-ratio, with NMAD -- the per-pixel version of the same question.

    ./lidar-icp/bin/python analysis/rcg/compare_2008_2025.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import odc.stac
import pandas as pd
import planetary_computer as pc
import pyarrow.parquet as pq
import pyproj
import rioxarray
import xarray as xr
from pystac_client import Client
from scipy.ndimage import label, uniform_filter

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; FIG = "figures"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res5 = 5.0; BUF = 300.0
to_ll = pyproj.Transformer.from_crs("EPSG:26915", "EPSG:4326", always_xy=True)
lo0, la0 = to_ll.transform(X0 - BUF, Y0 - BUF); lo1, la1 = to_ll.transform(X1 + BUF, Y1 + BUF); BBOX = [lo0, la0, lo1, la1]
cat = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
LSB = ["blue", "green", "red", "nir08", "swir16", "swir22"]; S2B = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
for i in ("LC09_L2SP_026029_20251004_02_T1", "LC09_L2SP_026029_20251105_02_T1", "LC08_L2SP_026029_20250926_02_T1"):
    out = f"{R}/ls_{i.split('_')[3]}_bands.tif"
    if os.path.exists(out): continue
    it = list(cat.search(collections=["landsat-c2-l2"], ids=[i]).items())[0]
    ds = odc.stac.load([it], bands=LSB + ["qa_pixel"], bbox=BBOX, crs="EPSG:26915", resolution=30, chunks={})
    qa = ds["qa_pixel"].squeeze(drop=True).astype("uint16"); bad = 0
    for b in (0, 1, 3, 4, 5): bad = bad | ((qa >> b) & 1)
    clear = (bad == 0) & (ds["red"].squeeze(drop=True) > 0)
    arr = xr.concat([(ds[b].squeeze(drop=True).astype("float32") * 0.0000275 - 0.2).where(clear) for b in LSB], dim="band").assign_coords(band=LSB)
    arr.rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(out, compress="deflate"); print(f"  fetched {i}: clear {float(clear.mean()):.3f}")
for d in ("2025-10-02", "2025-10-27"):
    out = f"{R}/s2_{d.replace('-', '')}_bands.tif"
    if os.path.exists(out): continue
    its = sorted([i for i in cat.search(collections=["sentinel-2-l2a"], bbox=BBOX, datetime=f"{d}T00:00:00Z/{d}T23:59:59Z").items() if "T15TWJ" in i.id], key=lambda i: i.id); it = its[-1]
    ds = odc.stac.load([it], bands=S2B + ["SCL"], bbox=BBOX, crs="EPSG:26915", resolution=10, chunks={})
    scl = ds["SCL"].squeeze(drop=True); clear = ~scl.isin([0, 1, 3, 8, 9, 10, 11])
    arr = xr.concat([(ds[b].squeeze(drop=True).astype("float32") / 10000.0).where(clear) for b in S2B], dim="band").assign_coords(band=S2B)
    arr.rio.write_crs("EPSG:26915").rio.write_nodata(np.nan).rio.to_raster(out, compress="deflate"); print(f"  fetched {it.id}: clear {float(clear.mean()):.3f}")

# populations (as band_separability.py)
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; D = np.load(f"{S}/floodplain_vegcorr_delta.npy") * 1e3
fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3
sp = np.load(f"{S}/gen1_spread_mm.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool); ny, nx = dod.shape
LO, HI, DZ = -5000.0, 5000.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
H = np.zeros((ny * nx, NB), np.int32); np.add.at(H, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1); H = H.reshape(ny, nx, NB)
Hp = np.rint(uniform_filter(H.astype(np.float32), size=(3, 3, 1), mode="nearest") * 9); C = np.cumsum(Hp, 2); tot = C[..., -1]
q = lambda p: ctr[np.argmax(C >= p * tot[..., None], axis=2)]; span_pool = np.where(tot > 0, q(0.9) - q(0.1), np.nan); del H, Hp, C
flat = fp & (sl <= 2) & np.isfinite(dod); F1 = flat & (span_pool >= 400) & (dod < -lod) & (D > 0) & ~pf; F2 = flat & pf
below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
wc = rioxarray.open_rasterio(f"{R}/worldcover.tif").squeeze(drop=True); wca = wc.values[::-1]; wcx0, wcy0 = float(wc.x[0]) - 5, float(wc.y[-1]) - 5
def frac_in(mask5, xs, ys, r):
    k = int(round(r / res5)); out = np.full((len(ys), len(xs)), np.nan)
    for j, y in enumerate(ys):
        iy0 = int((y - r / 2 - Y0) / res5); iy1 = iy0 + k
        if iy0 < 0 or iy1 > ny: continue
        for i, x in enumerate(xs):
            ix0 = int((x - r / 2 - X0) / res5); ix1 = ix0 + k
            if ix0 < 0 or ix1 > nx: continue
            out[j, i] = mask5[iy0:iy1, ix0:ix1].mean()
    return out
def at_px(arr, xs_, ys_, x0, y0, r):
    ix = np.clip(((xs_ - x0) / r).astype(int), 0, arr.shape[1] - 1); iy = np.clip(((ys_ - y0) / r).astype(int), 0, arr.shape[0] - 1)
    out = arr[np.ix_(iy, ix)].astype(float); out[(ys_[:, None] < y0) | (ys_[:, None] >= y0 + arr.shape[0] * r) | (xs_[None, :] < x0) | (xs_[None, :] >= x0 + arr.shape[1] * r)] = np.nan
    return out
nmad = lambda v: 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))

def load(tag, bn):
    da = rioxarray.open_rasterio(f"{R}/{tag}_bands.tif"); return da.x.values, da.y.values, {b: da.values[i] for i, b in enumerate(bn)}
popcache = {}
def pops_for(xs, ys, r):
    key = (len(xs), len(ys))
    if key not in popcache:
        wcv = at_px(wca, xs, ys, wcx0, wcy0, 10.0); P = {k: frac_in(v, xs, ys, r) for k, v in dict(F1=F1, F2=F2, NE=NE, flat=flat).items()}
        popcache[key] = [("KNOWN tree", wcv == 10), ("KNOWN cropland", wcv == 40), ("KNOWN grassland", wcv == 30), ("flat floodplain", P["flat"] > 0.5), ("F1 sedge", P["F1"] > 0.5), ("F2 closed canopy", P["F2"] > 0.5), ("NE patch", P["NE"] > 0.5)]
    return popcache[key]

rows = []
print("\nRELATIVE SPECTRA: population median / flat-floodplain median, per band, per scene (gauge-free)")
scenes = [("ls_20081029", LSB, 30.0), ("ls_20250926", LSB, 30.0), ("ls_20251004", LSB, 30.0), ("ls_20251105", LSB, 30.0),
          ("s2_20231023", S2B, 10.0), ("s2_20241017", S2B, 10.0), ("s2_20251002", S2B, 10.0), ("s2_20251027", S2B, 10.0)]
store = {}
for tag, bn, r in scenes:
    xs, ys, B = load(tag, bn); pops = pops_for(xs, ys, r); fin = np.all([np.isfinite(B[b]) for b in bn], axis=0)
    ref = {b: np.nanmedian(B[b][pops[3][1] & fin]) for b in bn}; store[tag] = (xs, ys, B, pops, fin)
    print(f"\n  {tag}  (flat floodplain median: " + ", ".join(f"{b} {ref[b]:.3f}" for b in bn) + ")")
    print(f"  {'population':<18}{'px':>7}" + "".join(f"{b:>8}" for b in bn) + f"{'NDVI':>8}{'NDMI':>8}")
    for name, m in pops:
        m = m & fin; nir, red = (B["nir08"], B["red"]) if "nir08" in B else (B["B08"], B["B04"]); sw = B["swir16"] if "swir16" in B else B["B11"]; nirm = B["nir08"] if "nir08" in B else B["B8A"]
        ndvi = np.nanmedian(((nir - red) / (nir + red))[m]); ndmi = np.nanmedian(((nirm - sw) / (nirm + sw))[m])
        rel = {b: np.nanmedian(B[b][m]) / ref[b] for b in bn}; rows.append(dict(scene=tag, population=name, px=int(m.sum()), NDVI=ndvi, NDMI=ndmi, **{"rel_" + b: rel[b] for b in bn}))
        print(f"  {name:<18}{m.sum():>7,}" + "".join(f"{rel[b]:>8.2f}" for b in bn) + f"{ndvi:>8.2f}{ndmi:>8.2f}")
pd.DataFrame(rows).to_csv(f"{R}/relative_spectra_2008_2025.csv", index=False)

print("\nPER-PIXEL CHANGE, Landsat 30 m: log(2025/2008) per band for the NE patch and F2, MINUS the flat floodplain's median log-ratio (so 0 = changed like its surroundings)")
xs0, ys0, B0, pops0, fin0 = store["ls_20081029"]
for tag in ("ls_20251004", "ls_20251105", "ls_20250926"):
    xs1, ys1, B1, pops1, fin1 = store[tag]
    if not (np.allclose(xs0, xs1) and np.allclose(ys0, ys1)): print(f"  {tag}: grid mismatch, skipped"); continue
    fin = fin0 & fin1; flatm = pops0[3][1] & fin
    print(f"  2008-10-29 -> {tag[3:]}")
    print(f"  {'population':<18}{'px':>6}" + "".join(f"{b:>14}" for b in LSB))
    for name, m in pops0:
        if name not in ("NE patch", "F2 closed canopy", "F1 sedge", "KNOWN grassland", "KNOWN cropland", "KNOWN tree"): continue
        mm = m & fin; line = f"  {name:<18}{mm.sum():>6,}"
        for b in LSB:
            lr = np.log(B1[b] / B0[b]); v = lr[mm] - np.nanmedian(lr[flatm]); line += f"{np.nanmedian(v):>+8.2f}({nmad(v):4.2f})"
        print(line)
fig, axes = plt.subplots(1, 2, figsize=(16, 5.5))
cols = {"KNOWN tree": "#008300", "KNOWN cropland": "#eda100", "KNOWN grassland": "#1baf7a", "F1 sedge": "#4a3aa7", "F2 closed canopy": "#eb6834", "NE patch": "#e34948"}
for ax, tags, bn, wl in ((axes[0], ["ls_20081029", "ls_20251004", "ls_20251105"], LSB, [485, 560, 660, 830, 1650, 2215]), (axes[1], ["s2_20231023", "s2_20241017", "s2_20251002", "s2_20251027"], S2B, [490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190])):
    df = pd.DataFrame(rows)
    for tag, ls in zip(tags, ["-", "--", ":", "-."]):
        for name in ("NE patch", "F1 sedge", "KNOWN grassland"):
            rr = df[(df.scene == tag) & (df.population == name)]
            if len(rr): ax.plot(wl, [rr["rel_" + b].values[0] for b in bn], color=cols[name], ls=ls, marker="o", ms=3, lw=1.6, label=f"{name} {tag[3:]}")
    ax.axhline(1, color="#52514e", lw=0.8); ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("reflectance / flat-floodplain median"); ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.set_title("Landsat 2008 vs 2025" if ax is axes[0] else "Sentinel-2 2023-2025", fontsize=10)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
fig.suptitle("NE patch, sedge cells and grassland: October reflectance relative to the flat floodplain, by year", fontsize=11)
fig.tight_layout(); fig.savefig(f"{FIG}/rcg_relative_spectra_2008_2025.png", dpi=110); print(f"\nfigure: {FIG}/rcg_relative_spectra_2008_2025.png; table: {R}/relative_spectra_2008_2025.csv")
