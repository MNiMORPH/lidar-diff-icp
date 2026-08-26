#!/usr/bin/env python3
"""WORKED EXAMPLE -- analysis/ridgelines/gen2_csf_compare.py, instrumented.

A copy, not a replacement: the original is untouched. The science is byte-identical; the
only additions are the Run() calls. Compare the two outputs to see exactly what the
banner buys.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python trust/example_instrumented.py
"""
import argparse, json, os, sys
import numpy as np, laspy
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lidar_diff_icp.refcells import reference_cells
from trust.provenance import Run

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--csf", default="data/csf_cache/elba_gen2_patch.las")
ap.add_argument("--min-pts", type=int, default=5)
A = ap.parse_args()

R = Run("does classifying BOTH epochs with our CSF close the gen1-gen2 cover-dependent gap?")

j = json.load(open(R.input(f"{A.tile}/corrections.json",
                           role="pipeline registration: bounds, resolution, applied shifts")))
b = j["bounds"]; RES = float(j["res_m"]); X0, Y0 = b[0], b[1]

# THE SUBSTITUTION THAT WENT UNNARRATED: this is the VENDOR class-2 grid, not our CSF.
zf = np.load(R.input(f"{A.tile}/z_after.npy",
                     role="gen2 ground surface as the pipeline has it = per-cell median "
                          "of the VENDOR class-2 points (NOT our CSF)"))
NY, NX = zf.shape
_zf = zf.copy(); _m = ~np.isfinite(_zf)
if _m.any():
    _zf = _zf[tuple(distance_transform_edt(_m, return_distances=False, return_indices=True))]
_gy, _gx = np.gradient(_zf, RES)
gxf = _gx.ravel(); gyf = _gy.ravel(); zflat = _zf.ravel()
nnorm = np.sqrt(gxf**2 + gyf**2 + 1.0)

g = laspy.read(R.input(A.csf, role="gen2 returns classified by OUR CSF, same parameters as gen1"))
x = np.asarray(g.x); y = np.asarray(g.y); z = np.asarray(g.z)
ix = ((x - X0)/RES).astype(np.int64); iy = ((y - Y0)/RES).astype(np.int64)
ing = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
cc = (iy[ing]*NX + ix[ing])
xc = X0 + ((cc % NX) + 0.5)*RES; yc = Y0 + ((cc // NX) + 0.5)*RES
h = (z[ing] - (zflat[cc] + gxf[cc]*(x[ing]-xc) + gyf[cc]*(y[ing]-yc)))/nnorm[cc] * 1000.0

R.param("min_pts", A.min_pts, src="repo")
order = np.argsort(cc, kind="stable")
cs = cc[order]; hs = h[order]
uniq, start = np.unique(cs, return_index=True)
stop = np.r_[start[1:], cs.size]
dz_csf = np.full(NY*NX, np.nan); npts = np.zeros(NY*NX, np.int64)
for u, a, bb in zip(uniq, start, stop):
    if bb - a >= A.min_pts:
        dz_csf[u] = np.median(hs[a:bb])
    npts[u] = bb - a

cover = np.load(R.input(f"{A.tile}/canopy_cover_pfs.npy",
                        role="per-cell canopy cover fraction from PyForestScan plant-area density"))
dod = np.load(R.input(f"{A.tile}/dod.npy",
                      role="DEM of difference, gen2 minus gen1, metres; + = elevation rose")).ravel() * 1000.0

# THE LABEL THAT NAMED NOTHING: slope_max=90 disables the slope cut entirely.
slope_max = R.param("slope_max", 90.0, src="MINE",
                    why="opened the slope cut from the repo default of 12 deg to take the "
                        "whole CSF patch; this admits steep mass-wasting cells that the "
                        "reference population is defined to exclude")
stable, rep = reference_cells(A.tile, slope_max=slope_max)
R.cuts("reference_cells", rep)
R.mask("stable", stable, of=stable.size,
       defn=f"reference_cells(slope_max={slope_max}) -- divide/curvature/building/"
            f"clearcut cuts, slope cut effectively DISABLED")
cover = cover.ravel()
ok = stable & np.isfinite(dz_csf) & np.isfinite(cover) & np.isfinite(dod)
R.mask("ok", ok, of=stable.size, defn="stable AND all three fields finite AND >=min_pts CSF ground points")

R.column("stratum", "canopy cover bin, dimensionless fraction of PyForestScan cover")
R.column("cells", "number of cells in the stratum, count")
R.column("z_csf-z_after", "median over cells of (our CSF ground - vendor class-2 ground), mm, slope-normal")
R.column("DoD_now", "median DoD as the pipeline reports it, mm, gen2-gen1")
R.banner()

rows = [("open   <0.05", -0.01, 0.05), ("light .05-.20", 0.05, 0.20),
        ("mid   .20-.35", 0.20, 0.35), ("dense  >0.35", 0.35, 1.01), ("ALL", -9, 9)]
tab = []
for nm, lo, hi in rows:
    m = ok if nm == "ALL" else (ok & (cover > lo) & (cover <= hi))
    if m.sum() < 20:
        continue
    d = np.median(dz_csf[m]); D = np.median(dod[m])
    tab.append([nm, f"{m.sum():,d}", f"{d:+.1f}", f"{D:+.1f}"])
R.table(["stratum", "cells", "z_csf-z_after", "DoD_now"], tab)
R.done(headline="cover-stratified CSF-minus-vendor gen2 ground offset")
