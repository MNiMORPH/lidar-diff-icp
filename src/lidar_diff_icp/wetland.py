"""Return-statistic water / wetland flag from airborne lidar.

Data-derived (no external hydrography), so it matches the epoch and resolution of
the data -- see [[landcover-attribution-tool]] for why external polygons are the
wrong tool here. It flags, it does not delete: open water is masked from change
detection (its "elevation" is a stage-dependent surface, not ground), while
wetland is labeled low-confidence so a water-level / emergent-vegetation signal is
not counted as terrestrial deposition (the Cook County finding: wetland interiors
rose ~0.19 m gen1->gen2).

Literature signatures (airborne NIR, ~1064 nm):
* OPEN WATER -> laser DROPOUT. Water strongly absorbs the pulse and reflects
  specularly away from the sensor, so returns are sparse or absent -> very low
  point density; ASPRS class 9 where the vendor set it. (Vetter et al. 2009;
  Hooshyar et al. 2015; Smeeckaert et al. 2013.)
* WETLAND / saturated / inundated-below-canopy -> ATTENUATED GROUND returns.
  Standing water and saturated soil beneath emergent/riparian vegetation cut the
  ground-return density, in FLAT, LOW-lying settings. Flatness + topographic-low
  separate this from dense-canopy shadow (also few ground returns, but on slopes /
  uplands). (Lang & McCarty 2009, intensity+returns for below-canopy inundation.)

Intensity -- the other classic water cue (water returns are low-intensity with high
relative variation) -- is NOT used here because the fetched cloud does not carry
it; gridding intensity would sharpen the open-water call.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter, distance_transform_edt as edt

DRY, WETLAND, OPEN_WATER = 0, 1, 2


def return_stats(path, bounds, res, nx, ny, *, chunk=10_000_000):
    """One streaming pass over a cloud -> per-cell (ground_density, total_density,
    class9_frac, dem). ground = ASPRS class 2; dem = mean ground elevation."""
    import laspy
    X0, Y0 = bounds[0], bounds[1]; N = nx * ny
    tot = np.zeros(N, np.int64); ng = np.zeros(N, np.int64)
    nw = np.zeros(N, np.int64); zg = np.zeros(N, float)
    with laspy.open(str(path)) as fh:
        for pts in fh.chunk_iterator(chunk):
            x = np.asarray(pts.x); y = np.asarray(pts.y); z = np.asarray(pts.z)
            cl = np.asarray(pts.classification)
            ix = ((x - X0) / res).astype(np.int64); iy = ((y - Y0) / res).astype(np.int64)
            ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
            f = iy[ok] * nx + ix[ok]; clo = cl[ok]; zo = z[ok]
            tot += np.bincount(f, minlength=N)
            g = clo == 2; ng += np.bincount(f[g], minlength=N)
            zg += np.bincount(f[g], weights=zo[g], minlength=N)
            nw += np.bincount(f[clo == 9], minlength=N)
    tot = tot.reshape(ny, nx); ng = ng.reshape(ny, nx); nw = nw.reshape(ny, nx)
    dem = np.where(ng > 0, zg.reshape(ny, nx) / np.maximum(ng, 1), np.nan)
    return dict(ground_density=ng, total_density=tot,
                class9_frac=nw / np.maximum(tot, 1.0), dem=dem)


def wetland_flag(ground_density, total_density, dem, res, *, class9_frac=None,
                 flat_slope_deg=2.0, tpi_window_m=300.0, dropout_frac=0.15,
                 ground_frac_pctile=20.0):
    """Return-statistic flag: 0 = dry land, 1 = wetland, 2 = open water.

    ``dropout_frac``: open water = total density below this fraction of the tile
    median (the dropout signature), and flat. ``ground_frac_pctile``: wetland =
    ground-return fraction below this percentile (attenuation), and flat, and in a
    topographic low. Both thresholds are data-relative so they travel across tiles.

    Assumes wetland is a MINORITY of the tile (the percentile threshold sits in the
    dry mode). On a wetland-dominated tile, raise ``ground_frac_pctile`` or bring in
    intensity. Without intensity this is a density/topography flag, not a classifier.
    """
    gd = np.asarray(ground_density, float); td = np.asarray(total_density, float)
    valid = td > 0
    Zf = np.asarray(dem, float).copy(); nn = ~np.isfinite(Zf)
    if nn.any():
        Zf = Zf[tuple(edt(nn, return_distances=False, return_indices=True))]
    gy, gx = np.gradient(gaussian_filter(Zf, 1.0), res)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    tpi = Zf - uniform_filter(Zf, size=max(int(2 * tpi_window_m / res), 3), mode="nearest")
    flat = slope < flat_slope_deg; low = tpi < 0

    med_td = np.median(td[valid]) if valid.any() else 1.0
    openw = valid & (td < dropout_frac * med_td) & flat          # dropout signature
    if class9_frac is not None:
        openw = openw | (np.asarray(class9_frac) > 0.5)          # vendor water

    gfrac = gd / np.maximum(td, 1.0)
    cand = valid & ~openw
    thr = np.percentile(gfrac[cand], ground_frac_pctile) if cand.any() else 0.0
    wet = cand & (gfrac < thr) & flat & low                     # attenuated ground, flat, low

    flag = np.full(td.shape, DRY, np.uint8)
    flag[wet] = WETLAND; flag[openw] = OPEN_WATER
    return flag
