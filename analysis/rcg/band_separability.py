#!/usr/bin/env python3
"""Redo the reed-canary-grass test the way the Wisconsin DNR did it: ALL reflective bands on
one mid-October scene, SWIR first, unsupervised clustering -- not NDVI.

What the reports actually say (read 2026-10-19; Bernthal & Willis 2004; Hatch & Bernthal
2008): Landsat bands 1,2,3,4,5,7 on a single scene from Oct 9-31; "multiple dates and NDVI
are not needed"; upland masked out first with the wetland inventory; ISODATA with 75
clusters, each cluster assigned by field truth, then maximum likelihood; "Band 5, the midIR
band (1.55-1.75 um), showed the greatest spectral separation"; overall accuracy 61-83%,
user's 72-92% for RCG; confusions: heavily grazed pasture, cranberry bogs; a single shrub
or standing water in a pixel breaks it (NIR-driven).

THIS SCRIPT, per scene (Landsat-5 2008-10-29 at 30 m; Sentinel-2 2024-10-17 and 2023-10-23
at 10 m), with the same KNOWN-first populations as before:
  1. median reflectance per band per population, and SWIR-based indices:
     NDMI = (nir - swir1)/(nir + swir1); NBR = (nir - swir2)/(nir + swir2); NDVI for reference;
     Sentinel-2 adds the red-edge NDRE = (B8A - B05)/(B8A + B05).
  2. SEPARABILITY d' = |median_A - median_B| / sqrt((NMAD_A^2 + NMAD_B^2)/2) for the
     questioned classes (F2 closed canopy, NE patch) against the rest of the flat
     floodplain, grassland, and the F1 sedge cells -- and for KNOWN pairs (cropland vs
     grassland, tree vs grassland) as the null, so a d' can be read against what a real
     class difference looks like in the same scene.
  3. WISCONSIN-STYLE CLUSTERING: k-means (k = 30, MINE; ISODATA's 75 was for a whole
     scene) on all bands over the floodplain-mask pixels only (our "wetland mask"); then
     for each cluster the share of its pixels that are F2 / NE, against the floodplain
     base rate. If the NE patch were a distinct spectral class it would concentrate in a
     few clusters with high purity; if it spreads across clusters at the base rate it is
     spectrally ordinary floodplain.

    ./lidar-icp/bin/python analysis/rcg/band_separability.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import rioxarray
from scipy.ndimage import label, uniform_filter
from sklearn.cluster import KMeans

S = "data/derived/elbaext_delong_star_nosurf"; R = f"{S}/rcg"; FIG = "figures"
X0, Y0, X1, Y1 = 575600.0, 4882200.0, 580050.0, 4886250.0; res5 = 5.0
dod = np.load(f"{S}/dod.npy") * 1e3; new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; D = np.load(f"{S}/floodplain_vegcorr_delta.npy") * 1e3
fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); sl = np.load(f"{S}/slope.npy"); lod = np.load(f"{S}/lod.npy") * 1e3
sp = np.load(f"{S}/gen1_spread_mm.npy"); pf = np.load(f"{S}/penetration_failure.npy").astype(bool); forest = np.load(f"{S}/forest_pfs.npy").astype(bool)
ny, nx = dod.shape
LO, HI, DZ = -5000.0, 5000.0, 20.0; edges = np.arange(LO, HI + DZ, DZ); ctr = 0.5 * (edges[:-1] + edges[1:]); NB = ctr.size
t = pq.read_table(f"{S}/beam_offset_table.parquet", columns=["cell", "d_mm_corr", "in_grid"]); g = t["in_grid"].to_numpy().astype(bool)
cell = t["cell"].to_numpy()[g]; d = t["d_mm_corr"].to_numpy()[g].astype(float); ok = np.isfinite(d) & (d > LO) & (d < HI)
H = np.zeros((ny * nx, NB), np.int32); np.add.at(H, (cell[ok], np.clip(((d[ok] - LO) / DZ).astype(int), 0, NB - 1)), 1); H = H.reshape(ny, nx, NB)
Hp = np.rint(uniform_filter(H.astype(np.float32), size=(3, 3, 1), mode="nearest") * 9); C = np.cumsum(Hp, 2); tot = C[..., -1]
q = lambda p: ctr[np.argmax(C >= p * tot[..., None], axis=2)]; span_pool = np.where(tot > 0, q(0.9) - q(0.1), np.nan); del H, Hp, C
flat = fp & (sl <= 2) & np.isfinite(dod)
F1 = flat & (span_pool >= 400) & (dod < -lod) & (D > 0) & ~pf; F2 = flat & pf; F3 = flat & (sp < 140) & (np.abs(dod) < 30)
below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); NE = lab == (np.argsort(sizes[1:])[::-1][0] + 1)
masks5 = dict(F1=F1, F2=F2, F3=F3, NE=NE, forest_hs=forest & (sl > 10), flat_fp=flat, fp=fp)
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
def dprime(a, b): return abs(np.nanmedian(a) - np.nanmedian(b)) / np.sqrt((nmad(a) ** 2 + nmad(b) ** 2) / 2 + 1e-12)


def run(tag, bands_tif, res, bandnames, idx_defs, k=30):
    da = rioxarray.open_rasterio(bands_tif); xs, ys = da.x.values, da.y.values; B = {b: da.sel(band=i + 1).values if "band" in da.dims and da.band.dtype != object else None for i, b in enumerate(bandnames)}
    B = {b: da.values[i] for i, b in enumerate(bandnames)}
    wcv = at_px(wca, xs, ys, wcx0, wcy0, 10.0)
    P = {kk: frac_in(v, xs, ys, res) for kk, v in masks5.items()}
    thr = 0.5
    pops = [("KNOWN tree (WC10)", wcv == 10), ("KNOWN cropland (WC40)", wcv == 40), ("KNOWN grassland (WC30)", wcv == 30), ("KNOWN forest hillslope", P["forest_hs"] > thr),
            ("KNOWN F3 cropped floodplain", P["F3"] > thr), ("all flat floodplain", P["flat_fp"] > thr), ("F1 sedge", P["F1"] > thr), ("F2 closed canopy", P["F2"] > thr), ("NE patch", P["NE"] > thr)]
    feats = dict(B); feats.update({nm: f(B) for nm, f in idx_defs.items()})
    keys = list(feats); finite = np.all([np.isfinite(feats[kk]) for kk in bandnames], axis=0)
    print(f"\n{'=' * 100}\n{tag}: median reflectance / index per population (NMAD in parentheses)")
    print(f"  {'population':<30}{'px':>7}" + "".join(f"{kk:>12}" for kk in keys))
    popv = {}
    for name, m in pops:
        m = m & finite; popv[name] = {kk: feats[kk][m] for kk in keys}
        print(f"  {name:<30}{m.sum():>7,}" + "".join(f"{np.nanmedian(feats[kk][m]):>7.3f}({nmad(feats[kk][m]):4.3f})" for kk in keys))
    print(f"\n  SEPARABILITY d' per band/index. NULL rows first (what a real class difference looks like in this scene):")
    print(f"  {'pair':<46}" + "".join(f"{kk:>9}" for kk in keys))
    pairs = [("KNOWN cropland (WC40)", "KNOWN grassland (WC30)"), ("KNOWN tree (WC10)", "KNOWN grassland (WC30)"), ("KNOWN F3 cropped floodplain", "all flat floodplain"),
             ("NE patch", "all flat floodplain"), ("NE patch", "KNOWN grassland (WC30)"), ("NE patch", "F1 sedge"), ("F2 closed canopy", "all flat floodplain"), ("F2 closed canopy", "F1 sedge")]
    rows = []
    for a, b in pairs:
        dd = {kk: dprime(popv[a][kk], popv[b][kk]) for kk in keys}; rows.append(dict(scene=tag, a=a, b=b, **dd))
        print(f"  {a[:22] + ' vs ' + b[:20]:<46}" + "".join(f"{dd[kk]:>9.2f}" for kk in keys))
    # Wisconsin-style clustering over the floodplain mask
    fpx = (P["fp"] > 0.2) & finite
    X = np.column_stack([feats[b][fpx] for b in bandnames]); X = (X - X.mean(0)) / X.std(0)
    km = KMeans(n_clusters=k, n_init=4, random_state=0).fit(X); cl = np.full(fpx.shape, -1); cl[fpx] = km.labels_
    base_ne = (P["NE"] > thr)[fpx].mean(); base_f2 = (P["F2"] > thr)[fpx].mean()
    print(f"\n  K-MEANS k={k} on {X.shape[0]:,} floodplain-mask pixels, all {len(bandnames)} bands. Base rates: NE {100*base_ne:.1f}%  F2 {100*base_f2:.1f}%")
    print(f"  clusters holding the NE patch (share of patch pixels, cluster purity = share of cluster that is NE / F2 / flat floodplain, cluster swir1):")
    ne_px = (P["NE"] > thr) & fpx; f2_px = (P["F2"] > thr) & fpx
    tab = []
    for c in range(k):
        mc = cl == c
        tab.append((c, mc.sum(), (ne_px & mc).sum() / max(ne_px.sum(), 1), (ne_px & mc).sum() / mc.sum(), (f2_px & mc).sum() / mc.sum(), ((P["flat_fp"] > thr) & mc).sum() / mc.sum(), np.nanmedian(feats[bandnames[4]][mc]), np.nanmedian(feats["NDVI"][mc])))
    tab.sort(key=lambda r: -r[2]); cum = 0
    print(f"  {'cluster':>8}{'n px':>7}{'NE share':>10}{'purity NE':>11}{'purity F2':>11}{'flat fp':>9}{bandnames[4]:>8}{'NDVI':>7}")
    for r in tab[:8]:
        cum += r[2]; print(f"  {r[0]:>8}{r[1]:>7,}{100*r[2]:>9.1f}%{100*r[3]:>10.1f}%{100*r[4]:>10.1f}%{100*r[5]:>8.0f}%{r[6]:>8.3f}{r[7]:>7.2f}")
    print(f"  -> top 3 clusters hold {100*sum(r[2] for r in tab[:3]):.0f}% of the NE patch; their NE purity {[f'{100*r[3]:.0f}%' for r in tab[:3]]} against base {100*base_ne:.1f}%")
    return pd.DataFrame(rows), popv, keys


L_B = ["blue", "green", "red", "nir08", "swir16", "swir22"]
L_IDX = {"NDVI": lambda B: (B["nir08"] - B["red"]) / (B["nir08"] + B["red"]), "NDMI": lambda B: (B["nir08"] - B["swir16"]) / (B["nir08"] + B["swir16"]), "NBR": lambda B: (B["nir08"] - B["swir22"]) / (B["nir08"] + B["swir22"])}
S_B = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
S_IDX = {"NDVI": lambda B: (B["B08"] - B["B04"]) / (B["B08"] + B["B04"]), "NDMI": lambda B: (B["B8A"] - B["B11"]) / (B["B8A"] + B["B11"]), "NBR": lambda B: (B["B8A"] - B["B12"]) / (B["B8A"] + B["B12"]), "NDRE": lambda B: (B["B8A"] - B["B05"]) / (B["B8A"] + B["B05"])}
out = []
r1, pv1, k1 = run("Landsat-5 2008-10-29 (30 m)", f"{R}/ls_20081029_bands.tif", 30.0, L_B, L_IDX); out.append(r1)
r2, pv2, k2 = run("Sentinel-2 2024-10-17 (10 m)", f"{R}/s2_20241017_bands.tif", 10.0, S_B, S_IDX); out.append(r2)
r3, pv3, k3 = run("Sentinel-2 2023-10-23 (10 m)", f"{R}/s2_20231023_bands.tif", 10.0, S_B, S_IDX); out.append(r3)
pd.concat(out).to_csv(f"{R}/band_separability.csv", index=False)
fig, axes = plt.subplots(1, 3, figsize=(20, 5.5))
cols = {"KNOWN tree (WC10)": "#008300", "KNOWN cropland (WC40)": "#eda100", "KNOWN grassland (WC30)": "#1baf7a", "KNOWN F3 cropped floodplain": "#2a78d6", "all flat floodplain": "#9e9e9e", "F1 sedge": "#4a3aa7", "F2 closed canopy": "#eb6834", "NE patch": "#e34948"}
for ax, (ttl, pv, bn, wl) in zip(axes, (("Landsat-5 2008-10-29", pv1, L_B, [485, 560, 660, 830, 1650, 2215]), ("Sentinel-2 2024-10-17", pv2, S_B, [490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190]), ("Sentinel-2 2023-10-23", pv3, S_B, [490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190]))):
    for name, col in cols.items():
        if name in pv: ax.plot(wl, [np.nanmedian(pv[name][b]) for b in bn], color=col, marker="o", ms=3, lw=1.8, ls="--" if name.startswith("KNOWN") else "-", label=name)
    ax.set_title(f"{ttl}: median reflectance by population", fontsize=10); ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("surface reflectance")
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
axes[0].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{FIG}/rcg_spectra_october.png", dpi=110); print(f"\nfigure: {FIG}/rcg_spectra_october.png; table: {R}/band_separability.csv")
