#!/usr/bin/env python3
"""Elba DoD with the vertical datum restricted to REAL OSM hard surfaces (baseball infield,
cemetery, parking) instead of the loose flat_hard default that admits fields. Saves the
reference cells (never lose them again), reports the tight bound + by-hardness consistency.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/run_osm_datum.py
"""
import json, time
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lidar_diff_icp.pipeline import difference_dem
from lidar_diff_icp.detect import detect_change_standard
from lidar_diff_icp.viz import hillshade
from lidar_diff_icp import references
from lidar_diff_icp.canopy import leafon_slope_flag, inflate_lod

BOUNDS = (577492.8, 4882737.6, 580032.8, 4886237.6); RES = 5.0
OUT = Path("data/derived/elba_refdatum"); FIG = "figures/refdatum"; Path(FIG).mkdir(parents=True, exist_ok=True)
z = np.load(str(OUT/"osm_hard_polys.npz"), allow_pickle=True)
polys = [z[k] for k in z.files if k.startswith("poly")]
print(f"{len(polys)} OSM hard-surface polys")

t = time.time()
r = difference_dem("data/before/4342-29-64.laz", "data/after/3dep2021_fd_class2.laz", BOUNDS,
                   res=RES, ground="slope_normal", ground_source="csf", after_ground="class2",
                   stream=False, robust_stable=True, csf_cache="data/csf_cache/elba.las",
                   tie="reference", allow_parabola=False, ref_polys=polys,
                   save_ref_cells=str(OUT/"ref_cells.npz"), datum_tilt=False)
dod, lod, Z21, stable = r["dod"], r["lod"], r["z_after"], r["stable"]
nx, ny = r["nx"], r["ny"]; X0, Y0 = r["bounds"][0], r["bounds"][1]
ce = r["corrections"]["cross_epoch_datum"]
print(f"done {time.time()-t:.0f}s")
print(f"DATUM: n_ref={ce['n_ref']}  const={1000*ce['const_m']:+.1f} mm  "
      f"tilt={1000*ce['tilt_mag_m_per_km']:.1f} mm/km  resid_nmad={1000*ce['resid_nmad_m']:.1f} mm  "
      f"SE_const≈{1000*ce['resid_nmad_m']/np.sqrt(ce['n_ref']):.2f} mm")

# by-hardness consistency (is the offset the same on pavement vs softer ground?)
cells = dict(np.load(str(OUT/"ref_cells.npz")))
bh = references.datum_offset(cells)["by_hardness"]
print("by-hardness offset (median mm, n):", bh)

# leaf-on LoD flag (penetration from the full cloud) + detect
pen = np.load(str(OUT.parent/"elba_fulldensity/penetration.npy"))
Zf = Z21.copy(); m = ~np.isfinite(Zf)
if m.any():
    from scipy.ndimage import distance_transform_edt as edt
    Zf = Zf[tuple(edt(m, return_distances=False, return_indices=True))]
sl = np.degrees(np.arctan(np.hypot(*np.gradient(Zf, RES)[::-1])))
lod = inflate_lod(lod, leafon_slope_flag(pen, sl))
det = detect_change_standard(dod, lod, stable, RES); change = det["change"]

# crest forest/open (does the canopy result survive the clean datum?)
crest = np.load(str(OUT.parent/"elba_fulldensity/crest_mask.npy")); cm = crest & np.isfinite(dod)
fo = cm & (pen < 0.25); op = cm & (pen >= 0.45)
print(f"crest forest={np.median(dod[fo])*1000:+.1f} open={np.median(dod[op])*1000:+.1f} "
      f"contrast={(np.median(dod[fo])-np.median(dod[op]))*1000:+.1f} mm")

np.save(OUT/"dod_osm.npy", dod); np.save(OUT/"lod_osm.npy", lod)
json.dump(r["corrections"], open(OUT/"corrections_osm.json", "w"), indent=2)

# plots: DoD+LoD, change, and the reference-cell locations
ext = (X0, X0+nx*RES, Y0, Y0+ny*RES); hs = hillshade(Zf, RES, X0, Y0, fill_gaps=True)
for tag, arr, kw, title in [("dod_osm", dod, dict(cmap="RdBu", vmin=-.3, vmax=.3),
                             "DoD gen2-gen1 (m) [OSM hard-surface datum]"),
                            ("change_osm", np.where(change, dod, np.nan), dict(cmap="RdBu", vmin=-.3, vmax=.3),
                             "change above detection [OSM hard-surface datum]")]:
    fig, ax = plt.subplots(figsize=(10, 11)); ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
    im = ax.imshow(arr, extent=ext, origin="lower", alpha=0.7, **kw)
    ax.set_title(f"elba: {title}"); fig.colorbar(im, ax=ax, shrink=0.6, extend="both")
    ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
    fig.savefig(f"{FIG}/elba_{tag}.png", dpi=130, bbox_inches="tight"); plt.close(fig)
# reference cells + polys
fig, ax = plt.subplots(figsize=(10, 11)); ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
ax.scatter(cells["x"], cells["y"], s=2, c="red", label=f"datum cells (n={len(cells['x'])})")
for p in polys: ax.plot(*np.vstack([p, p[:1]]).T, "c-", lw=1)
ax.set_title("OSM hard-surface datum reference cells"); ax.legend()
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
fig.savefig(f"{FIG}/elba_datum_cells.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print(f"saved dod_osm/lod_osm/ref_cells + figs -> {OUT} / {FIG}")
