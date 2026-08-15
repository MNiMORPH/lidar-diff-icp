#!/usr/bin/env python3
"""Quality-check a fetched NAIP mosaic and its k-means cover classification.

Reports band/NDVI distributions and nodata, and writes two figures: RGB/NIR/NDVI
and RGB/clusters. In dissected terrain NAIP carries optical *shadow* on steep
slopes (steep forested walls read as low/negative NDVI), which the figures make
visible -- a caution when using NAIP cover to stratify lidar error.

Example:
    python analysis/naip_qa.py data/naip/naip2010_4m.npz --figdir figures
"""
import argparse, os, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from naip_cover_error import classify   # reuse the committed classifier


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("naip_npz")
    ap.add_argument("--figdir", default="figures")
    args = ap.parse_args()
    stem = Path(args.naip_npz).stem
    Path(args.figdir).mkdir(parents=True, exist_ok=True)

    d = np.load(args.naip_npz, allow_pickle=True)
    rgbn = d["rgbn"]; ndvi = d["ndvi"]
    minx, miny, maxx, maxy = d["bounds"]; res = float(d["res"])
    R, G, B, NIR = (rgbn[i].astype(float) for i in range(4))
    ext = (minx, maxx, miny, maxy)
    nodata = rgbn.sum(0) == 0
    valid = ~nodata
    print(f"{stem}: {ndvi.shape} @ {res} m,  nodata {100*nodata.mean():.1f}%")
    for name, arr in [("R", R), ("NIR", NIR), ("NDVI", ndvi)]:
        p = [round(float(np.percentile(arr[valid], q)), 3) for q in (2, 25, 50, 75, 98)]
        print(f"  {name:5s} pct[2,25,50,75,98]: {p}")

    # figure 1: RGB / NIR / NDVI
    fig, ax = plt.subplots(1, 3, figsize=(15, 6))
    rgb = np.dstack([R, G, B]); rgb = (rgb / max(rgb.max(), 1) * 255).clip(0, 255).astype(np.uint8)
    ax[0].imshow(rgb, extent=ext); ax[0].set_title("NAIP RGB")
    ax[1].imshow(NIR, extent=ext, cmap="gray"); ax[1].set_title("NIR")
    im = ax[2].imshow(ndvi, extent=ext, cmap="RdYlGn", vmin=-0.2, vmax=0.4)
    ax[2].set_title("NDVI (red = shadow/water)"); fig.colorbar(im, ax=ax[2], shrink=0.6)
    f1 = Path(args.figdir) / f"{stem}_qa.png"
    fig.savefig(f1, dpi=90, bbox_inches="tight"); plt.close(fig)

    # figure 2: RGB / k-means clusters
    cl = classify(args.naip_npz)
    lab = cl["lab"]
    print("cluster  frac%  NDVI  NIR  texture  bright")
    for c in range(cl["k"]):
        s = lab == c
        print(f"  {c}  {100*s.mean():5.1f}  {cl['ndvi'][s].mean():5.2f} "
              f"{cl['NIR'][s].mean():4.0f} {cl['tex'][s].mean():7.1f} {cl['bright'][s].mean():6.0f}")
    fig, ax = plt.subplots(1, 2, figsize=(11, 6))
    ax[0].imshow(rgb, extent=ext); ax[0].set_title("NAIP RGB")
    ax[1].imshow(lab, extent=ext, cmap="tab10"); ax[1].set_title("k-means cover clusters")
    f2 = Path(args.figdir) / f"{stem}_clusters.png"
    fig.savefig(f2, dpi=90, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {f1}  and  {f2}")


if __name__ == "__main__":
    main()
