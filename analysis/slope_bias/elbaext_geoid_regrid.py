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
from lidar_diff_icp import references, io
BOUNDS = (575600.0, 4882200.0, 580050.0, 4886250.0); RES = 5.0
BEFORE = "data/before/elbaext_gen1_merged.laz"
AFTER  = "data/after/elbaext_3dep_fd_class2.laz"     # class2-extracted: loadable non-streaming
# ONE shared vertical frame, DERIVED not hardcoded.
#
# geoid_difference() fits a PLANE over the bounds it is given, about that bounds' centroid.
# The true GEOID03-GEOID18 field is curved (4.44 mm ptp over the shared area), so fitting it
# separately on elba's and elbaext's different footprints yields two different linearisations:
# measured, they disagree by 2.70 mm ptp over the shared area, which is larger than the
# tile-to-tile swath-tie agreement we are trying to preserve. "Both auto-computed" is the
# same METHOD but not the same FRAME, and it is the frame that has to match.
#
# So: fit ONCE on ELBA's bounds and re-express that same plane about elbaext's centroid --
# algebraically identical to elba's datum (0.000 mm ptp difference), and the more accurate of
# the two linearisations against the PROJ field (0.44 vs 0.84 mm RMS over the shared area).
ELBA_BOUNDS = (577492.8, 4882737.6, 580032.8, 4886237.6)   # analysis/ridgelines/run_elba_dod.py
_a, _b, _c = references.geoid_difference(ELBA_BOUNDS, io.MN_GEN1_CRS)
_cx_e, _cy_e = 0.5 * (ELBA_BOUNDS[0] + ELBA_BOUNDS[2]), 0.5 * (ELBA_BOUNDS[1] + ELBA_BOUNDS[3])
_cx_x, _cy_x = 0.5 * (BOUNDS[0] + BOUNDS[2]), 0.5 * (BOUNDS[1] + BOUNDS[3])
G = (_a + _b * (_cx_x - _cx_e) / 1000.0 + _c * (_cy_x - _cy_e) / 1000.0, _b, _c)
print(f"shared geoid plane (elba's, re-centred on elbaext): const {G[0]*1000:+.3f} mm, "
      f"tilt ({G[1]*1000:+.3f},{G[2]*1000:+.3f}) mm/km")
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
