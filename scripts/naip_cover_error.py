#!/usr/bin/env python3
"""Independent (NAIP) land cover vs. lidar inter-swath error.

Motivation: classifying cover from the lidar itself is circular when the goal is
to characterise/correct the lidar's own error - the cover proxy shares the
error's noise. NAIP (an independent optical sensor) breaks that loop. This
script classifies cover from a NAIP mosaic (k-means on NDVI/NIR/texture/
brightness), overlays it on the co-registered inter-swath residual, and:

1. Re-tests the cover->error relationship with the INDEPENDENT stratifier
   (compare the Spearman here with the ~+0.6 obtained from lidar-derived
   veg-height; a large drop means that number was self-reference).
2. Decomposes the residual per cover class into a LOCAL component (variogram
   nugget) and a BETWEEN-swath navigation component (long-range correlated
   sill), testing whether the local part is cover-driven while the navigation
   part is best seen on flat, low-cover ground.

Caveats: in dissected terrain NAIP carries terrain *shadow* on steep slopes, so
the "forest" cluster conflates canopy with slope-aliased error and NAIP's own
shadow artifact. And note precision vs. accuracy: low inter-swath residual over
grass means the (elevated) grass surface is *repeatable* between passes, not
that it is accurate bare earth. Whether that biases the 2008->2021 *difference*
depends on whether the cover changed between epochs, which is NOT established
here; a cover that is unchanged cancels in the difference.

Example:
    python scripts/fetch_naip.py --bounds 577492.8 4882737.6 580035.0 4886238.3 \
        --year 2010 --res 4 --out data/naip/naip2010_4m.npz
    python scripts/naip_cover_error.py data/before/4342-29-64.laz \
        data/naip/naip2010_4m.npz --pair 136 137
"""
import argparse
import numpy as np
from scipy.ndimage import uniform_filter
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
import laspy

from lidar_diff_icp.swathdiff import _median_grid
from lidar_diff_icp import coreg, variogram as vg


def _nmad(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def classify(npz, k=6, seed=0):
    d = np.load(npz, allow_pickle=True)
    rgbn = d["rgbn"]; ndvi = d["ndvi"]
    minx, miny, maxx, maxy = d["bounds"]; res = float(d["res"])
    R, G, B, NIR = (rgbn[i].astype(float) for i in range(4))
    bright = (R + G + B) / 3.0
    m = uniform_filter(NIR, 5); m2 = uniform_filter(NIR * NIR, 5)
    tex = np.sqrt(np.clip(m2 - m * m, 0, None))
    F = np.stack([ndvi, NIR, tex, bright], -1).reshape(-1, 4)
    F = (F - F.mean(0)) / F.std(0)
    lab = KMeans(n_clusters=k, n_init=5, random_state=seed).fit_predict(F).reshape(ndvi.shape)
    return dict(lab=lab, ndvi=ndvi, NIR=NIR, tex=tex, bright=bright,
                bounds=(minx, miny, maxx, maxy), res=res, k=k)


def residual(laz, a, b, res):
    f = laspy.read(laz)
    x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
    psid = np.asarray(f.point_source_id); cls = np.asarray(f.classification)
    rn = np.asarray(f.return_number); nr = np.asarray(f.number_of_returns)
    filt = (nr == 1) & ~np.isin(cls, [5, 6, 9])
    ma = filt & (psid == a); mb = filt & (psid == b)
    x0 = max(x[ma].min(), x[mb].min()); x1 = min(x[ma].max(), x[mb].max())
    y0 = max(y[ma].min(), y[mb].min()); y1 = min(y[ma].max(), y[mb].max())
    nx = int(np.ceil((x1 - x0) / res)); ny = int(np.ceil((y1 - y0) / res))
    zr = _median_grid(x[ma], y[ma], z[ma], res, x0, y0, nx, ny)
    zs = _median_grid(x[mb], y[mb], z[mb], res, x0, y0, nx, ny)
    c = coreg.nuth_kaab(zr, zs, res)
    r = zr - (coreg._shift_grid(zs, c.dx, c.dy, res) + c.dz)
    gy, gx = np.mgrid[0:ny, 0:nx]; mm = np.isfinite(r)
    return x0 + (gx[mm] + 0.5) * res, y0 + (gy[mm] + 0.5) * res, r[mm]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("laz"); ap.add_argument("naip_npz")
    ap.add_argument("--pair", nargs=2, type=int, default=[136, 137])
    args = ap.parse_args()

    cl = classify(args.naip_npz)
    lab = cl["lab"]; ndvi = cl["ndvi"]; NIR = cl["NIR"]
    minx, miny, maxx, maxy = cl["bounds"]; res = cl["res"]; ny, nx = ndvi.shape
    print("cluster  frac%   NDVI  NIR  texture  bright")
    for c in range(cl["k"]):
        s = lab == c
        print(f"  {c}   {100*s.mean():5.1f}  {ndvi[s].mean():5.2f} {NIR[s].mean():4.0f} "
              f"{cl['tex'][s].mean():7.1f} {cl['bright'][s].mean():6.0f}")

    rx, ry, rv = residual(args.laz, args.pair[0], args.pair[1], res)
    col = np.clip(((rx - minx) / res).astype(int), 0, nx - 1)
    row = np.clip(((maxy - ry) / res).astype(int), 0, ny - 1)
    lab_at = lab[row, col]; ndvi_at = ndvi[row, col]; nir_at = NIR[row, col]

    print("\nresidual by cover cluster + variance decomposition:")
    print(f"  {'cl':>3} {'ncells':>7} {'NMAD':>6} {'nugget(local)':>13} {'sill(nav)':>10} {'range':>6}")
    for c in range(cl["k"]):
        s = lab_at == c
        if s.sum() < 800:
            continue
        mod = vg.fit_spherical(*vg.empirical_variogram(rx[s], ry[s], rv[s], 400))
        print(f"  {c:>3} {s.sum():>7,} {_nmad(rv[s]):>6.3f} {mod.nugget:>13.4f} "
              f"{mod.sill:>10.4f} {mod.range_:>6.0f}")

    lit = nir_at > np.percentile(NIR, 50)
    print(f"\nSpearman(NAIP NDVI, |residual|): all={spearmanr(ndvi_at, np.abs(rv))[0]:+.3f}  "
          f"lit-only={spearmanr(ndvi_at[lit], np.abs(rv[lit]))[0]:+.3f}")
    print("Compare with ~+0.60 from lidar-derived veg-height: a large drop = that")
    print("was largely self-reference (the independence problem this test targets).")


if __name__ == "__main__":
    main()
