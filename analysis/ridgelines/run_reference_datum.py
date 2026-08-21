#!/usr/bin/env python3
"""Elba DoD the PRINCIPLED way, no parabola: reuse cached CSF gen1, register gen1->gen2 by
align_swaths -> Nuth&Kaeaeb lateral -> pad const+tilt vertical -> drift (tie='reference',
allow_parabola=False), differenced against the FULL-density gen2 class-2 ground. Emits the
usual two plots. If the pad datum can't run it RAISES (parabola deactivated) rather than
silently degrading.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/run_reference_datum.py
"""
import json, time
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lidar_diff_icp.pipeline import difference_dem
from lidar_diff_icp.detect import detect_change_standard
from lidar_diff_icp.viz import hillshade
from lidar_diff_icp.canopy import leafon_slope_flag, inflate_lod

BEFORE = "data/before/4342-29-64.laz"
AFTER  = "data/after/3dep2021_fd_class2.laz"           # full-density gen2 class-2 ground
BOUNDS = (577492.8, 4882737.6, 580032.8, 4886237.6)
RES = 5.0
OUT = Path("data/derived/elba_refdatum"); OUT.mkdir(parents=True, exist_ok=True)
FIG = "figures/refdatum"; Path(FIG).mkdir(parents=True, exist_ok=True)


def fig_dod_lod(name, Z21, dod, lod, res, X0, Y0, nx, ny, figdir):
    hs = hillshade(Z21, res, X0, Y0, fill_gaps=False)
    ext = (X0, X0 + nx * res, Y0, Y0 + ny * res); v = 0.3
    fig, ax = plt.subplots(1, 2, figsize=(15, 9))
    ax[0].imshow(hs, extent=ext, origin="lower", cmap="gray", alpha=0.6)
    im0 = ax[0].imshow(dod, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v)
    ax[0].set_title(f"{name}: DEM of Difference (gridded ground): gen2 - gen1 (m)\n"
                    "red = erosion, blue = deposition  [PAD DATUM, no parabola]")
    fig.colorbar(im0, ax=ax[0], shrink=0.6, extend="both")
    im1 = ax[1].imshow(lod, extent=ext, origin="lower", cmap="viridis", vmin=0, vmax=0.2)
    ax[1].set_title("level of detection (m)")
    fig.colorbar(im1, ax=ax[1], shrink=0.6, extend="max")
    for a in ax: a.set_xlabel("Easting (m)"); a.set_ylabel("Northing (m)")
    out = f"{figdir}/{name}_dod_lod.png"; fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    return out


def fig_dem_change(name, Z21, dod, change, regions, res, X0, Y0, nx, ny, figdir):
    hs = hillshade(Z21, res, X0, Y0, fill_gaps=True)
    ext = (X0, X0 + nx * res, Y0, Y0 + ny * res); v = 0.3
    over = np.where(change, dod, np.nan)
    fig, ax = plt.subplots(figsize=(10, 11))
    ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
    im = ax.imshow(over, extent=ext, origin="lower", cmap="RdBu", vmin=-v, vmax=v, alpha=0.7)
    ax.set_title(f"{name}: topographic change above detection limits  [PAD DATUM, no parabola]")
    ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
    fig.colorbar(im, ax=ax, shrink=0.6, extend="both", label="detected DoD (m)")
    out = f"{figdir}/{name}_change.png"; fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    return out


t = time.time()
print("difference_dem  tie=reference  allow_parabola=False  stream=False (in-memory)...", flush=True)
r = difference_dem(BEFORE, AFTER, BOUNDS, res=RES, ground="slope_normal",
                   ground_source="csf", after_ground="class2", stream=False,
                   robust_stable=True, csf_cache="data/csf_cache/elba.las",
                   tie="reference", allow_parabola=False)
dod, lod, Z21, stable = r["dod"], r["lod"], r["z_after"], r["stable"]
nx, ny = r["nx"], r["ny"]; X0, Y0 = r["bounds"][0], r["bounds"][1]
ex = np.isfinite(dod)
print(f"done {time.time()-t:.0f}s  tie_method={r['corrections'].get('cross_epoch_tie', r['corrections'].get('method','?'))}  "
      f"sigma={r['stable_sigma']:.3f}  medDoD={np.nanmedian(dod[ex])*1000:+.1f} mm  medLoD={np.nanmedian(lod):.3f}", flush=True)

# leaf-on LoD flag using penetration from the FULL cloud (class-2 file alone has pen==1)
pen = np.load("data/derived/elba_fulldensity/penetration.npy")
Zf = Z21.copy(); _nm = ~np.isfinite(Zf)
if _nm.any():
    from scipy.ndimage import distance_transform_edt as _edt
    Zf = Zf[tuple(_edt(_nm, return_distances=False, return_indices=True))]
sl = np.degrees(np.arctan(np.hypot(*np.gradient(Zf, RES)[::-1])))
leafon = leafon_slope_flag(pen, sl); lod = inflate_lod(lod, leafon)

det = detect_change_standard(dod, lod, stable, RES)
change, regions = det["change"], det["regions"]
print(f"detected {len(regions)} regions  tau_sys={det['tau_sys_m']:.3f} m", flush=True)

np.save(OUT/"dod.npy", dod); np.save(OUT/"lod.npy", lod); np.save(OUT/"z_after.npy", Z21)
np.save(OUT/"change.npy", change); np.save(OUT/"stable.npy", stable)
json.dump(r["corrections"], open(OUT/"corrections.json", "w"), indent=2)
fa = fig_dod_lod("elba", Z21, dod, lod, RES, X0, Y0, nx, ny, FIG)
fb = fig_dem_change("elba", Z21, dod, change, regions, RES, X0, Y0, nx, ny, FIG)
print(f"wrote {fa}\n      {fb}\nsaved -> {OUT}/", flush=True)
