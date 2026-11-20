#!/usr/bin/env python
"""2021 lidar intensity over the field site at 1 m: mean and standard deviation of ALL returns per cell, plus the count
(Andy, 2026-11-20: "Use all the returns and bin per square meter: mean and standard deviation").
Clouds: data/after/valley/{4342-29-64,4358-29-01}_3dep_full.laz (the files overlap by 35 m; returns east of E 580030
are taken from the second file only). Raw intensity, no range, angle or gain correction; the two swaths (3040, 3041)
are mixed in each cell, so a per-cell swath share is written too. Extent (MINE, argv overrides as E0 E1 N0 N1):
E 579780-580160, N 4885280-4885760, the floodplain field with a margin, EPSG:6344.
Outputs  data/field/2026-10-25/intensity/intensity_2021_1m_{mean,sd,count,share_swath3040}.tif
         data/field/2026-10-25/intensity/intensity_2021_1m_log.txt
         figures/rtk_2026-10-25_floodplain_intensity_1m.png (mean and SD, survey points as small dots)
"""
import os, sys
import numpy as np, pandas as pd, laspy, rasterio
from rasterio.transform import from_origin
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
OUT = "data/field/2026-10-25/intensity"; os.makedirs(OUT, exist_ok=True); LOG = open(f"{OUT}/intensity_2021_1m_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
E0, E1, N0, N1 = [float(v) for v in sys.argv[1:5]] if len(sys.argv) >= 5 else (579780.0, 580160.0, 4885280.0, 4885760.0)
RES = 1.0; nx, ny = int((E1 - E0) / RES), int((N1 - N0) / RES)
acc = {k: np.zeros(ny * nx) for k in ("n", "s", "s2", "n3040")}; n_all = 0
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= E0) & (x < E1) & (y >= N0) & (y < N1)
            if t == "4358-29-01": m &= x > 580030.0
            if not m.any(): continue
            i = np.clip(((N1 - y[m]) / RES).astype(int), 0, ny - 1); j = np.clip(((x[m] - E0) / RES).astype(int), 0, nx - 1); cell = i * nx + j
            inten = np.asarray(ch.intensity)[m].astype(float); ps = np.asarray(ch.point_source_id)[m]; n_all += m.sum()
            acc["n"] += np.bincount(cell, minlength=ny * nx); acc["s"] += np.bincount(cell, weights=inten, minlength=ny * nx); acc["s2"] += np.bincount(cell, weights=inten ** 2, minlength=ny * nx); acc["n3040"] += np.bincount(cell[ps == 3040], minlength=ny * nx)
n = acc["n"].reshape(ny, nx)
with np.errstate(invalid="ignore", divide="ignore"):
    mean = (acc["s"] / acc["n"]).reshape(ny, nx); sd = np.sqrt(np.maximum(acc["s2"] / acc["n"] - (acc["s"] / acc["n"]) ** 2, 0)).reshape(ny, nx); share = (acc["n3040"] / acc["n"]).reshape(ny, nx)
mean[n == 0] = np.nan; sd[n == 0] = np.nan; share[n == 0] = np.nan; sd[n < 2] = np.nan
say(f"extent E {E0:.0f}-{E1:.0f} N {N0:.0f}-{N1:.0f}, {nx} x {ny} cells at {RES:.0f} m; returns {n_all:,}; cells with returns {int((n > 0).sum())} of {nx * ny}; returns per cell median {np.median(n[n > 0]):.0f}, p10 {np.percentile(n[n > 0], 10):.0f}, p90 {np.percentile(n[n > 0], 90):.0f}")
say(f"mean intensity over cells: median {np.nanmedian(mean):.0f}, p02 {np.nanpercentile(mean, 2):.0f}, p98 {np.nanpercentile(mean, 98):.0f}; SD within cell: median {np.nanmedian(sd):.0f}, p02 {np.nanpercentile(sd, 2):.0f}, p98 {np.nanpercentile(sd, 98):.0f}")
tr = from_origin(E0, N1, RES, RES); prof = dict(driver="GTiff", height=ny, width=nx, count=1, crs="EPSG:6344", transform=tr, compress="deflate")
for name, a, dt, nod in (("mean", mean, "float32", -9999.0), ("sd", sd, "float32", -9999.0), ("count", n, "int32", 0), ("share_swath3040", share, "float32", -9999.0)):
    with rasterio.open(f"{OUT}/intensity_2021_1m_{name}.tif", "w", dtype=dt, nodata=nod, **prof) as dst:
        dst.write(np.where(np.isfinite(a), a, nod).astype(dt), 1)
    say(f"wrote {OUT}/intensity_2021_1m_{name}.tif")
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv"); F = P[P.group.str.startswith("floodplain") & (P.status == "FIX")]
fig, axs = plt.subplots(1, 2, figsize=(14, 8)); ext = (E0, E1, N0, N1)
for ax, a, lab, cm, lim in ((axs[0], mean, "mean intensity of all 2021 returns per 1 m cell", "gray", (np.nanpercentile(mean, 2), np.nanpercentile(mean, 98))), (axs[1], sd, "standard deviation of intensity per 1 m cell", "magma", (0, np.nanpercentile(sd, 98)))):
    im = ax.imshow(a, extent=ext, cmap=cm, vmin=lim[0], vmax=lim[1]); plt.colorbar(im, ax=ax, shrink=0.6, label=lab); ax.plot(F.easting_tied, F.northing_tied, ".", color="#d55e00", ms=2); ax.set_aspect("equal"); ax.tick_params(labelsize=7)
fig.tight_layout(); fig.savefig("figures/rtk_2026-10-25_floodplain_intensity_1m.png", dpi=int(os.environ.get("FIGDPI", 170))); say("wrote figures/rtk_2026-10-25_floodplain_intensity_1m.png")
LOG.close()
