#!/usr/bin/env python3
"""Geoid-only datum DoD: gen1 converted GEOID03->GEOID18 (const +67mm + tilt), lateral
Nuth-Kaeaeb shift kept, NO pad const, no parabola. Reusable geoid_datum capability."""
import json, numpy as np
from lidar_diff_icp.pipeline import difference_dem
G = (0.067, 0.00061, -0.00073)   # (N03-N18) const_m, b(E), c(N) m/km -> ADD to gen1
r = difference_dem("data/before/4342-29-64.laz","data/after/3dep2021_fd_class2.laz",
    (577492.8,4882737.6,580032.8,4886237.6),res=5.0,ground="slope_normal",ground_source="csf",
    after_ground="class2",stream=False,robust_stable=True,csf_cache="data/csf_cache/elba.las",
    tie="reference",geoid_datum=G)
dod=r["dod"]; np.save("data/derived/elba_refdatum/dod_geoid.npy",dod)
json.dump(r["corrections"],open("data/derived/elba_refdatum/corrections_geoid.json","w"),indent=2)
crest=np.load("data/derived/elba_fulldensity/crest_mask.npy"); pen=np.load("data/derived/elba_fulldensity/penetration.npy")
osm=np.load("data/derived/elba_refdatum/dod_osm.npy")
for tag,d in [("pad (dod_osm)",osm),("GEOID-only (new)",dod)]:
    cm=crest&np.isfinite(d); fo=cm&(pen<0.25); op=cm&(pen>=0.45)
    print(f"  {tag:18s}: forest {np.median(d[fo])*1000:+6.1f}  open {np.median(d[op])*1000:+6.1f}  "
          f"contrast {(np.median(d[fo])-np.median(d[op]))*1000:+6.1f}  medDoD {np.nanmedian(d[np.isfinite(d)])*1000:+.1f} mm")
print("saved dod_geoid.npy"); 
