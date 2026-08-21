#!/usr/bin/env python3
"""Apply the two post-hoc corrections to the clean OSM-datum DoD:
 (1) GEOID TILT: remove grad(N_GEOID03 - N_GEOID18) = the gen1/gen2 geoid-model difference tilt
     (b=+0.00061, c=-0.00073 m/km; const already absorbed by the OSM datum). Small (~1 mm/km).
 (2) FOREST-FLOOR OFFSET: remove the veg-dependent lidar offset f(veg_frac) (veg_frac = fraction
     of returns >0.5 m; the R^2=0.72 predictor). Open cells (veg~0) unchanged.
Saves data/derived/elba_refdatum/dod_corrected.npy and reports the crest forest/open contrast
at each stage. Toggle with the two flags.
    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python analysis/ridgelines/dod_corrections.py
"""
import numpy as np
APPLY_GEOID_TILT = True
APPLY_VEG_OFFSET = True

D = "data/derived/elba_fulldensity/"
dod = np.load("data/derived/elba_refdatum/dod_osm.npy"); ny, nx = dod.shape
X0, Y0, RES = 577492.8, 4882737.6, 5.0
cx, cy = 578762.8, 4884487.6                       # tile centre (matches corrections norm)
crest = np.load(D + "crest_mask.npy"); pen = np.load(D + "penetration.npy")

def contrast(d, tag):
    cm = crest & np.isfinite(d); fo = cm & (pen < 0.25); op = cm & (pen >= 0.45)
    f, o = np.median(d[fo])*1000, np.median(d[op])*1000
    print(f"  {tag:18s}: forest {f:+6.1f}  open {o:+6.1f}  contrast {f-o:+6.1f} mm  "
          f"(medDoD {np.nanmedian(d[np.isfinite(d)])*1000:+.1f})")

contrast(dod, "raw OSM-datum")
d = dod.copy()

if APPLY_GEOID_TILT:                               # planar geoid-difference tilt, m
    b, c = 0.00061, -0.00073                        # m/km (East, North)
    E = X0 + (np.arange(nx)+0.5)*RES; N = Y0 + (np.arange(ny)+0.5)*RES
    EE, NN = np.meshgrid(E, N)
    tilt = (b*(EE-cx) + c*(NN-cy)) / 1000.0         # m
    d = d - tilt
    contrast(d, "- geoid tilt")

if APPLY_VEG_OFFSET:                                # remove veg-dependent forest-floor offset
    veg = np.load(D + "canopy_struct.npz")["veg_frac"]
    fin = np.isfinite(d) & np.isfinite(veg)
    # robust-ish linear f: dod = m*veg + a; remove only the veg-dependent part (m*veg), keep a as datum
    A = np.c_[veg[fin], np.ones(fin.sum())]
    m_, a_ = np.linalg.lstsq(A, d[fin], rcond=None)[0]
    print(f"  f(veg): slope m={m_*1000:+.1f} mm per unit veg_frac (removing m*veg_frac)")
    d = d - m_*np.nan_to_num(veg)
    contrast(d, "- veg offset")

np.save("data/derived/elba_refdatum/dod_corrected.npy", d)
print("saved data/derived/elba_refdatum/dod_corrected.npy")
