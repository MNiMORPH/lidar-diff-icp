#!/usr/bin/env python3
"""Air photos of the NE patch, 2008-2023, from the MnGeo WMS, with the surviving-patch outlines.

Andy, 2026-10-18: "Perhaps use imagery." The MnGeo image server
(https://imageserver.gisdata.mn.gov/cgi-bin/wms) serves statewide FSA/NAIP years as
layers fsa2008, fsa2009, ..., fsa2023 plus cir2008; JPEG only; EPSG:26915. This fetches the
window 579300-580050 E, 4885300-4886250 N at 1 m and overlays the four largest below-LoD
patches that survive the vegetation correction.

Also run here, because it changed a conclusion: scan angle along a flight line is a
function of across-track POSITION (corr(scan angle, easting) = +0.995 on line 138 in this
window), so a patch/ring contrast in return density vs scan angle is geometry, not surface.
An earlier 'specular surface' reading of that contrast is WITHDRAWN.

    ./lidar-icp/bin/python analysis/floodplain_ne_imagery.py   (needs the network)
"""
import json
import os
import subprocess

import laspy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.image import imread
from scipy.ndimage import label

S = "data/derived/elbaext_delong_star_nosurf"; OUT = "data/derived/ne_patch_imagery"; os.makedirs(OUT, exist_ok=True)
E0, E1, N0, N1 = 579300, 580050, 4885300, 4886250
LAYERS = [("cir2008", "2008 CIR"), ("fsa2008", "2008 NAIP"), ("fsa2009", "2009 NAIP"), ("fsa2010", "2010 NAIP"), ("fsa2013", "2013"),
          ("fsa2015", "2015"), ("fsa2017", "2017"), ("fsa2019", "2019"), ("fsa2021", "2021 NAIP"), ("fsa2023", "2023")]
for L, _ in LAYERS:
    p = f"{OUT}/{L}.jpg"
    if not os.path.exists(p):
        subprocess.run(["curl", "-s", "-m", "120", f"https://imageserver.gisdata.mn.gov/cgi-bin/wms?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS={L}&STYLES=&SRS=EPSG:26915"
                        f"&BBOX={E0},{N0},{E1},{N1}&WIDTH={E1-E0}&HEIGHT={N1-N0}&FORMAT=image/jpeg", "-o", p], check=True)
new = np.load(f"{S}/dod_floodplain_vegcorr.npy") * 1e3; fp = np.load(f"{S}/floodplain_mask.npy").astype(bool); lod = np.load(f"{S}/lod.npy") * 1e3
c = json.load(open(f"{S}/corrections.json")); X0, Y0 = c["bounds"][0], c["bounds"][1]; res = c["res_m"]; ny, nx = new.shape
below = fp & np.isfinite(new) & (new < -lod); lab, n = label(below, structure=np.ones((3, 3))); sizes = np.bincount(lab.ravel()); ks = np.argsort(sizes[1:])[::-1][:4] + 1
c0, c1, r0, r1 = int((E0 - X0) / res), int((E1 - X0) / res), int((N0 - Y0) / res), int((N1 - Y0) / res)
yy, xx = np.mgrid[r0:r1, c0:c1]; Eg = X0 + (xx + 0.5) * res; Ng = Y0 + (yy + 0.5) * res
fig, axes = plt.subplots(2, 5, figsize=(22, 11))
for ax, (L, t) in zip(axes.ravel(), LAYERS):
    ax.imshow(imread(f"{OUT}/{L}.jpg"), extent=[E0, E1, N0, N1], origin="upper")
    for k in ks: ax.contour(Eg, Ng, (lab[r0:r1, c0:c1] == k).astype(float), levels=[0.5], colors="yellow", linewidths=1.0)
    ax.set_title(t, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("NE patch (yellow = the four largest surviving below-LoD patches) in MnGeo air photos 2008-2023; gen1 lidar Nov 2008, gen2 May 2021", fontsize=12)
fig.tight_layout(); fig.savefig("figures/floodplain_ne_imagery.png", dpi=90)
f = laspy.read("data/before/elbaext_gen1_merged.laz"); x = np.asarray(f.x); y = np.asarray(f.y); ps = np.asarray(f.point_source_id); sa = np.asarray(f.scan_angle_rank).astype(float)
win = (ps == 138) & (x > E0) & (x < E1) & (y > N0) & (y < N1)
print(f"line 138 in the window: corr(scan angle, easting) = {np.corrcoef(sa[win], x[win])[0,1]:+.3f}  -> scan angle is across-track position")
print("figure: figures/floodplain_ne_imagery.png")
