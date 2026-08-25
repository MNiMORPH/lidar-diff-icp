#!/usr/bin/env python3
"""The elba DEM of Difference, with every correction the pipeline knows about.

Supersedes the Aug-21 elba_refdatum/dod_geoid.npy, which passed a HARDCODED geoid triple
-- and the wrong tile's tilt at that (elbaext's +0.00061/-0.00073 rather than elba's
+0.00078/-0.00057; sub-mm across the tile, but a hardcode where the code now derives the
value) -- and predates the geoid-only datum refactor.

Applied here, in the pipeline's order: per-swath internal alignment to the lowest flight
line -> constant lateral Nuth & Kaeaeb tie -> geoid-difference datum AUTO-COMPUTED from the
PROJ grids -> per-swath along-track GNSS drift. Boresight residual is off: the search for
one returned nothing resolvable, so applying a value would be fitting noise.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/run_elba_dod.py
"""
import json, time
import numpy as np
from lidar_diff_icp.pipeline import difference_dem

OUT = "data/derived/elba_fulldensity"
t0 = time.time()
r = difference_dem(
    "data/before/4342-29-64.laz", "data/after/3dep2021_fd_class2.laz",
    (577492.8, 4882737.6, 580032.8, 4886237.6),
    res=5.0, ground="slope_normal", ground_source="csf", after_ground="class2",
    stream=False, robust_stable=True, csf_cache="data/csf_cache/elba.las",
    tie="reference",              # geoid-difference datum after the lateral shift
    geoid_datum=None,             # AUTO-COMPUTE from the PROJ geoid grids (no hardcode)
    along_track_drift=True,
    correct_boresight=False)      # residual searched for, none resolvable
dod = r["dod"]; lod = r["lod"]
np.save(f"{OUT}/dod.npy", dod); np.save(f"{OUT}/lod.npy", lod)
json.dump(r["corrections"], open(f"{OUT}/corrections.json", "w"), indent=2)
ok = np.isfinite(dod)
print(f"elapsed {time.time()-t0:.0f}s   stable_sigma {r['stable_sigma']*1000:.1f} mm")
print(f"DoD: n={ok.sum():,}  median {np.nanmedian(dod[ok])*1000:+.1f} mm  "
      f"NMAD {1.4826*np.nanmedian(np.abs(dod[ok]-np.nanmedian(dod[ok])))*1000:.1f} mm")
print(f"LoD: median {np.nanmedian(lod)*1000:.1f} mm   |DoD|>LoD: "
      f"{100*np.nanmean(np.abs(dod[ok])>lod[ok]):.1f}% of cells")
print(f"datum: {r['corrections']['cross_epoch_datum']}")
print(f"saved {OUT}/dod.npy, lod.npy, corrections.json")
