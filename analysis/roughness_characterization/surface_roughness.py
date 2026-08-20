#!/usr/bin/env python3
"""Roughness characterization at Battle Creek — how gen1 (2008 MN) and gen2
(3DEP 2021) lidar render within-cell roughness across KNOWN surface types.

Battle Creek is a CHARACTERIZATION site, not a study site: it was chosen because
it carries clean, ground-truthable surface types (sealed pavement, skinned
softball infields, mowed turf, prairie) that natural study sites lack. The goal
is to separate, per epoch, the instrumental ranging noise (measured on
near-zero-roughness built surfaces) from real surface roughness, and to find
where the two epochs agree (hard surfaces -> sensor-limited) versus diverge
(vegetation -> the sparse 2008 survey over-reads canopy).

Roughness metric: NMAD of the residuals to a per-cell least-squares plane
(detrended within-cell roughness; NMAD, not RMS, because RMS is inflated by
blunder returns under/near vegetation). Reported per region as the p20 "floor"
(smoothest cells ~ the surface's clean value) and the median (typical cell).

DATA PROVENANCE (all under data/, git-ignored; reacquire as below)
  frame (UTM 15N, EPSG:26915): 498135 4975136 499365 4976876  (4x the pilot tile)
  gen1: 4 Ramsey metro sub-tiles merged -> data/before_battlecreek/gen1_4tile.laz
        4342-03-32_b_a, 4342-03-32_a_b, 4342-02-32_c_d, 4342-02-32_d_c
        (download via lidar_diff_icp.tiles.download_tile(name, ..., county="ramsey"))
  gen2: data/after_battlecreek/gen2_4tile.laz
        scripts/fetch_3dep_curl.py --auto --bounds <frame> --max-depth 12
  NAIP: data/naip/naip2010_bc4tile.npz
        scripts/fetch_naip.py --bounds <frame> --year 2010 --res 2

Run:  env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/roughness_characterization/surface_roughness.py
"""
import os
import numpy as np
import laspy
import pandas as pd

from lidar_diff_icp import coreg

FRAME = (498135.0, 4975136.0, 499365.0, 4976876.0)      # x0, y0, x1, y1
RES = 5.0
SIGMA = {"gen1": 0.017, "gen2": 0.011}                   # instrumental floor (infield, upper bound)
GEN1 = "data/before_battlecreek/gen1_4tile.laz"
GEN2 = "data/after_battlecreek/gen2_4tile.laz"
NAIP = "data/naip/naip2010_bc4tile.npz"
CACHE = "data/derived/bc_roughness_char.npz"             # rasters, so reruns are instant

# Surface regions: box = (x0, x1, y0, y1) UTM15N; `mat` = material predicate on
# (R, B, ndvi, slope) at 5 m. Boxes located from NAIP + the classified overlay.
REGIONS = {
    "infield":  dict(box=(498980, 499230, 4976180, 4976400),
                     mat=lambda R, B, nd, s: (nd < 0.15) & ((R - B) > 20) & (s < 2),
                     desc="skinned softball diamonds (clay/sand, S~0)"),
    "parking":  dict(box=(498870, 499160, 4976090, 4976145),
                     mat=lambda R, B, nd, s: (nd < 0.15) & (np.abs(R - B) < 15) & (s < 2),
                     desc="single sealed parking lot just south of the ballfields"),
    "grass":    dict(box=(498980, 499230, 4976180, 4976400),
                     mat=lambda R, B, nd, s: (nd > 0.30) & (s < 2),
                     desc="mowed outfield turf between the diamonds"),
    "prairie":  dict(box=(498780, 499180, 4975680, 4976030),
                     mat=lambda R, B, nd, s: (nd > 0.25) & (nd < 0.55) & (s < 3),
                     desc="restored prairie grass (taller/rougher than mowed turf)"),
}


def nmad_roughness(path, class2, x0, y0, nx, ny, res=RES):
    """Per-cell NMAD of residuals to a least-squares plane (ny x nx), plus the
    per-cell point count and the plane-intercept ground elevation (for slope)."""
    f = laspy.read(path)
    x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
    m = (np.asarray(f.classification) == 2) if class2 else \
        (np.asarray(f.return_number) == np.asarray(f.number_of_returns))
    x, y, z = x[m], y[m], z[m]
    inb = (x >= x0) & (x < x0 + nx * res) & (y >= y0) & (y < y0 + ny * res)
    x, y, z = x[inb], y[inb], z[inb]
    ix = np.floor((x - x0) / res).astype(np.int64); iy = np.floor((y - y0) / res).astype(np.int64)
    f2 = iy * nx + ix; N = nx * ny
    u = x - (x0 + (ix + 0.5) * res); v = y - (y0 + (iy + 0.5) * res)
    S = lambda w: np.bincount(f2, w, minlength=N)
    n = S(np.ones_like(z)); Su = S(u); Sv = S(v); Sz = S(z)
    Suu = S(u * u); Suv = S(u * v); Svv = S(v * v); Suz = S(u * z); Svz = S(v * z); Szz = S(z * z)
    M = np.empty((N, 3, 3))
    M[:, 0, 0] = Suu; M[:, 0, 1] = Suv; M[:, 0, 2] = Su
    M[:, 1, 0] = Suv; M[:, 1, 1] = Svv; M[:, 1, 2] = Sv
    M[:, 2, 0] = Su;  M[:, 2, 1] = Sv;  M[:, 2, 2] = n
    r = np.stack([Suz, Svz, Sz], axis=1)
    valid = (n >= 6) & (np.abs(np.linalg.det(M)) > 1e-6)
    beta = np.zeros((N, 3)); beta[valid] = np.linalg.solve(M[valid], r[valid])
    pred = beta[f2, 0] * u + beta[f2, 1] * v + beta[f2, 2]; resid = z - pred
    med = pd.Series(resid).groupby(f2).transform("median")
    nm = pd.Series(np.abs(resid - med.values)).groupby(f2).median()
    rough = np.full(N, np.nan); rough[nm.index.values] = 1.4826 * nm.values; rough[~valid] = np.nan
    grd = np.full(N, np.nan); grd[valid] = beta[valid, 2]
    return rough.reshape(ny, nx), n.reshape(ny, nx).astype(float), grd.reshape(ny, nx)


def naip_grids(x0, y0, nx, ny, res=RES):
    """Resample NAIP red, blue, NDVI to the 5 m grid (nearest cell center)."""
    nz = np.load(NAIP); ndvi = nz["ndvi"]; rgb = nz["rgbn"]; nb = nz["bounds"]; nyN, nxN = ndvi.shape
    R = np.full((ny, nx), np.nan); B = np.full((ny, nx), np.nan); ND = np.full((ny, nx), np.nan)
    for iy in range(ny):
        for ix in range(nx):
            cx = x0 + (ix + 0.5) * res; cy = y0 + (iy + 0.5) * res
            px = int((cx - nb[0]) / (nb[2] - nb[0]) * nxN); py = int((nb[3] - cy) / (nb[3] - nb[1]) * nyN)
            if 0 <= px < nxN and 0 <= py < nyN:
                R[iy, ix] = rgb[0, py, px]; B[iy, ix] = rgb[2, py, px]; ND[iy, ix] = ndvi[py, px]
    return R, B, ND


def load_or_build():
    x0, y0, x1, y1 = FRAME
    nx = int(round((x1 - x0) / RES)); ny = int(round((y1 - y0) / RES))
    if os.path.exists(CACHE):
        d = np.load(CACHE)
        return (x0, y0, nx, ny, d["r1"], d["r2"], d["R"], d["B"], d["ND"], d["slope"])
    print("computing roughness rasters: CSF gen1 (~10 min) + class-2 gen2 ...", flush=True)
    # MODERN processing only: gen1 bare earth = PDAL CSF ground (not last-return,
    # which reads canopy in tall vegetation); gen2 = 3DEP class-2 ground.
    from lidar_diff_icp.ground import classify_ground_csf
    import shutil
    csf = classify_ground_csf(GEN1)
    r1, _, _ = nmad_roughness(csf, True, x0, y0, nx, ny)
    shutil.rmtree(os.path.dirname(csf), ignore_errors=True)
    r2, _, dem = nmad_roughness(GEN2, True, x0, y0, nx, ny)
    R, B, ND = naip_grids(x0, y0, nx, ny)
    Zf = dem.copy(); nanm = ~np.isfinite(Zf)
    if nanm.any():
        from scipy.ndimage import distance_transform_edt as edt
        Zf = Zf[tuple(edt(nanm, return_distances=False, return_indices=True))]
    slope = np.degrees(coreg.slope_aspect(Zf, RES)[0])
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez(CACHE, r1=r1, r2=r2, R=R, B=B, ND=ND, slope=slope)
    return x0, y0, nx, ny, r1, r2, R, B, ND, slope


def main():
    x0, y0, nx, ny, r1, r2, R, B, ND, slope = load_or_build()
    ix = np.arange(nx); iy = np.arange(ny)
    EE, NN = np.meshgrid(x0 + (ix + 0.5) * RES, y0 + (iy + 0.5) * RES)
    print(f"{'surface':9s} {'n':>4} | {'gen1 p20':>8} {'gen1 med':>8} | {'gen2 p20':>8} {'gen2 med':>8} "
          f"| {'S(gen1)':>7} {'S(gen2)':>7} {'ratio':>5}   description")
    for name, spec in REGIONS.items():
        bx = spec["box"]
        inbox = (EE >= bx[0]) & (EE <= bx[1]) & (NN >= bx[2]) & (NN <= bx[3])
        m = inbox & spec["mat"](R, B, ND, slope) & np.isfinite(r1) & np.isfinite(r2)
        a = r1[m]; b = r2[m]
        # surface roughness after removing each epoch's instrumental floor
        s1 = np.sqrt(max(np.median(a) ** 2 - SIGMA["gen1"] ** 2, 0))
        s2 = np.sqrt(max(np.median(b) ** 2 - SIGMA["gen2"] ** 2, 0))
        ratio = s1 / s2 if s2 > 1e-4 else np.nan
        print(f"{name:9s} {m.sum():4d} | {np.percentile(a,20):8.4f} {np.median(a):8.4f} "
              f"| {np.percentile(b,20):8.4f} {np.median(b):8.4f} | {s1:7.4f} {s2:7.4f} {ratio:5.2f}   {spec['desc']}")


if __name__ == "__main__":
    main()
