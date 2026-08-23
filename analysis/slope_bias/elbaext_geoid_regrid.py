#!/usr/bin/env python3
"""Rebuild elbaext on the GEOID datum (GEOID03->GEOID18), matching elba's canonical tie
so elba and elbaext share ONE vertical reference frame. Mirrors run_geoid_datum.py with
elbaext inputs. The prior elbaext product used a reference_plane fit whose 382k flat-hard
cells are dominated by flat rural fields (leaf-on gen2 sits high -> over-drop), so the
deterministic geoid is the trustworthy shared frame. Writes *_geoid outputs alongside the
reference_plane product (non-destructive); z_after/slope are tie-independent (reused).

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/slope_bias/elbaext_geoid_regrid.py
"""
import json, numpy as np
from lidar_diff_icp.pipeline import difference_dem
BOUNDS = (575600.0, 4882200.0, 580050.0, 4886250.0); RES = 5.0
BEFORE = "data/before/elbaext_gen1_merged.laz"
AFTER  = "data/after/elbaext_3dep_fd_class2.laz"     # class2-extracted: loadable non-streaming
G = (0.067, 0.00061, -0.00073)   # GEOID03->GEOID18 const_m, b(E), c(N) m/km -> ADD to gen1 (same as elba)
r = difference_dem(BEFORE, AFTER, BOUNDS, res=RES, ground="slope_normal", ground_source="csf",
    after_ground="class2", stream=False, robust_stable=True, csf_cache="data/csf_cache/elbaext.las",
    tie="reference", geoid_datum=G)
O = "data/derived/elbaext/"
np.save(O+"dod_geoid.npy", r["dod"]); np.save(O+"lod_geoid.npy", r["lod"])
np.save(O+"stable_geoid.npy", r["stable"])
json.dump(r["corrections"], open(O+"corrections_geoid.json", "w"), indent=2)
ce = r["corrections"]["cross_epoch_datum"]
print("cross_epoch_datum method:", ce["method"], " const_mm",
      1000*ce.get("const_m", float("nan")))
d = r["dod"]; ex = np.isfinite(d)
print(f"stable_sigma {r['stable_sigma']:.3f}  medDoD {np.nanmedian(d[ex])*1000:+.1f} mm  ({ex.sum()} finite)")
print("saved -> data/derived/elbaext/  (dod_geoid, lod_geoid, stable_geoid, corrections_geoid)")
