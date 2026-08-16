#!/usr/bin/env python3
"""Compute the level of detection (LoD) by land-cover zone.

Uses the independent NAIP cover classification and the co-registered inter-swath
residual. Because the error correlation length (>=370 m) exceeds any feature we
would map, N_eff ~ 1 and the LoD is ~flat with feature size (spatial averaging
does not help) - so the LoD is set by the point-scale error.

Definitions (z = 1.96 for 95%):
  sigma_2008(zone)     = NMAD_interswath / sqrt(2)     # one-survey vertical error
  LoD_selfconsistency  = z * NMAD_interswath           # detect a 2008 pass-to-pass diff
  LoD_epoch(2008->2021)= z * sqrt(sigma_2008^2 + sigma_2021^2)   # needs sigma_2021

Example:
    python scripts/lod.py data/before/4342-29-64.laz data/naip/naip2010_4m.npz
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from naip_cover_error import classify, residual, _nmad

Z = 1.96
FLAT = "flat-open (grass/crop/ag)"
ROUGH = "forest / steep valley wall"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("laz"); ap.add_argument("naip_npz")
    ap.add_argument("--pair", nargs=2, type=int, default=[136, 137])
    args = ap.parse_args()

    cl = classify(args.naip_npz)
    lab = cl["lab"]; NIR = cl["NIR"]; minx, miny, maxx, maxy = cl["bounds"]
    res = cl["res"]; ny, nx = lab.shape
    rx, ry, rv = residual(args.laz, args.pair[0], args.pair[1], res)
    col = np.clip(((rx - minx) / res).astype(int), 0, nx - 1)
    row = np.clip(((maxy - ry) / res).astype(int), 0, ny - 1)
    lab_at = lab[row, col]

    # group clusters into flat-open (smooth, lit) vs rough (forest/shadow/steep)
    # by their NIR-texture signature computed in classify()
    tex_mean = {c: cl["tex"][lab == c].mean() for c in range(cl["k"])}
    nir_mean = {c: NIR[lab == c].mean() for c in range(cl["k"])}
    flat_clusters = [c for c in range(cl["k"]) if tex_mean[c] < 20 and nir_mean[c] > 90]

    groups = {FLAT: np.isin(lab_at, flat_clusters),
              ROUGH: ~np.isin(lab_at, flat_clusters)}
    print(f"pair {args.pair[0]}-{args.pair[1]};  flat clusters = {flat_clusters}\n")
    print(f"  {'zone':>28} {'ncells':>8} {'NMAD':>6} {'sig2008':>7} {'LoD_self95':>10}")
    sig = {}
    for name, mask in groups.items():
        nm = _nmad(rv[mask]); s08 = nm / np.sqrt(2); sig[name] = s08
        print(f"  {name:>28} {mask.sum():>8,} {nm:>6.3f} {s08:>7.3f} {Z*nm:>10.3f}")

    print("\n  projected 2008->2021 change LoD95 = z*sqrt(sig2008^2 + sig2021^2)")
    print("  (sigma_2021 provisional until the 3DEP is measured):")
    print(f"  {'sigma_2021 ->':>28} " + "".join(f"{s:>9.2f}m" for s in [0.05, 0.10, 0.15]))
    for name in groups:
        row_ = [Z * np.sqrt(sig[name]**2 + s21**2) for s21 in (0.05, 0.10, 0.15)]
        print(f"  {name:>28} " + "".join(f"{v:>9.3f}m" for v in row_))
    print("\n  Note: LoD is ~flat with feature size (correlation length >=370 m,")
    print("  N_eff~1). A cover unchanged between epochs cancels in the difference;")
    print("  only a *change* in cover would add a term (not established here).")


if __name__ == "__main__":
    main()
