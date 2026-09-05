#!/usr/bin/env python3
"""Compare the two swath-network weightings on a tile's real cloud.

    w = n           the shipped behaviour until 2026-09-05: a tie is trusted in proportion
                    to how many cells it saw
    w = 1/dz_var    information weighting: a tie is trusted in proportion to how well it is
                    DETERMINED, which for the intercept tie is dominated by how far it
                    extrapolates to dtan = 0, not by n

Both solved on the SAME edge set, so the only difference is the weighting.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/swath_network_weighting.py --csf data/csf_cache/battlecreek.las
"""
import argparse

import laspy
import numpy as np

from lidar_diff_icp import coreg, io as lio

ap = argparse.ArgumentParser()
ap.add_argument("--csf", required=True, help="the tile's CSF-classified gen1 cloud")
ap.add_argument("--tie", default="intercept", choices=("intercept", "overlap_median"))
A = ap.parse_args()

f = laspy.read(A.csf)
ps = np.asarray(f.point_source_id)
pc = lio.PointCloud(np.asarray(f.x), np.asarray(f.y), np.asarray(f.z), ps,
                    np.asarray(f.classification), np.zeros(len(ps)),
                    np.asarray(f.scan_angle).astype(float) * 0.006, "EPSG:26915")
ref = int(ps.min())
corr, edges, _ = coreg.align_swaths(pc, ref=ref, tie=A.tie)

idx = {s: i for i, s in enumerate(sorted(set(ps.tolist())))}
E, N = len(edges), len(idx)
print(f"{A.csf}: {N} swaths, {E} edges -> "
      f"{'LOOPS, so the weighting matters' if E > N - 1 else 'a TREE: any positive weights give the SAME solution'}")


def solve(weights):
    A_ = np.zeros((E, N)); O = np.zeros(E)
    for e, ed in enumerate(edges):
        A_[e, idx[ed[0]]] = -1.0; A_[e, idx[ed[1]]] = 1.0; O[e] = ed[4]
    sw = np.sqrt(np.asarray(weights, float))
    c, *_ = np.linalg.lstsq(A_ * sw[:, None], O * sw, rcond=None)
    c = c - c[idx[ref]]
    return c, A_ @ c - O


old, mis_old = solve([e[5] for e in edges])          # w = n
new, mis_new = solve([e[6] for e in edges])          # w = 1/dz_var

tot = sum(e[6] for e in edges) or 1.0
print(f"\n{'pair':14s} {'dz (mm)':>12s} {'n':>9s} {'w = 1/dz_var':>14s} {'w share':>10s}")
for e in sorted(edges, key=lambda e: -e[6]):
    print(f"  {e[0]}-{e[1]:<8d} {1000*e[4]:>12.3f} {int(e[5]):>9d} {e[6]:>14.4e} "
          f"{100*e[6]/tot:>9.4f}%")

print(f"\n{'swath':10s} {'w = n (mm)':>14s} {'w = 1/var (mm)':>16s} {'change (mm)':>13s}")
for s, i in sorted(idx.items()):
    print(f"  {s:<8d} {1000*old[i]:>14.4f} {1000*new[i]:>16.4f} {1000*(new[i]-old[i]):>+13.4f}")

j = int(np.argmax(np.abs(mis_new)))
good = [i for i in range(E) if i != j]
print(f"\nlargest |misclosure|            w=n {1000*np.max(np.abs(mis_old)):9.3f} mm   "
      f"w=1/var {1000*np.max(np.abs(mis_new)):9.3f} mm")
print(f"  the w=1/var maximum sits on {edges[j][0]}-{edges[j][1]}, "
      f"weight share {100*edges[j][6]/tot:.4f}% -- a downweighted outlier SHOWS UP as a "
      f"large residual instead of moving the solution")
if good:
    print(f"largest EXCLUDING that edge     w=n {1000*np.max(np.abs(mis_old[good])):9.3f} mm   "
          f"w=1/var {1000*np.max(np.abs(mis_new[good])):9.3f} mm   <- the real consistency")
