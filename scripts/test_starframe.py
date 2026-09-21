#!/usr/bin/env python3
"""Synthetic recovery test for starframe.fit_star_shifts -- it LOCKS the sign convention.

A sign error here is invisible in real data: the fit still reports small standard errors
and the product still looks plausible, it is just registered the wrong way. So the
convention is established by displacing a known cloud and reading back what the estimator
says, never by reasoning about it in a docstring.

Construction: a smooth synthetic surface with real relief is the "gen2" reference. A
"gen1" cloud is sampled from the SAME surface and then displaced by a known
(sx, sy, sz). Whatever the estimator returns is compared with that displacement, and the
correction that ACTUALLY restores the cloud is reported.

    ./lidar-icp/bin/python scripts/test_starframe.py
"""
import numpy as np

from lidar_diff_icp.starframe import fit_star_shifts

RNG = np.random.default_rng(0)
NX = NY = 160
RES = 5.0
X0 = Y0 = 0.0
grid = (X0, Y0, RES, NX, NY)

ex = (np.arange(NX) + 0.5) * RES
ny_ = (np.arange(NY) + 0.5) * RES
EX, NYY = np.meshgrid(ex, ny_)


def surface(xx, yy):
    """Relief with varied aspect, so both gradient components are identifiable."""
    return (18.0 * np.sin(xx / 190.0) + 13.0 * np.cos(yy / 145.0)
            + 6.0 * np.sin((xx + yy) / 95.0))


Zref = surface(EX, NYY)


def ground_of(xx, yy, zz):
    """Per-cell median -- the same rule applied to both epochs, which is what matters."""
    ix = np.clip(((xx - X0) / RES).astype(int), 0, NX - 1)
    iy = np.clip(((yy - Y0) / RES).astype(int), 0, NY - 1)
    out = np.full(NY * NX, np.nan)
    k = iy * NX + ix
    order = np.argsort(k, kind="stable")
    ks, zs = k[order], zz[order]
    edges = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    for a, b in zip(edges, np.r_[edges[1:], ks.size]):
        out[ks[a]] = np.median(zs[a:b])
    return out.reshape(NY, NX)


n = 400_000
px = RNG.uniform(X0, X0 + NX * RES, n)
py = RNG.uniform(Y0, Y0 + NY * RES, n)
pz = surface(px, py)
sid = RNG.integers(0, 2, n)

SHIFT = {0: (0.60, -0.25, 0.030), 1: (-0.40, 0.15, -0.020)}   # metres, KNOWN
qx, qy, qz = px.copy(), py.copy(), pz.copy()
for s, (sx, sy, sz) in SHIFT.items():
    m = sid == s
    qx[m] += sx; qy[m] += sy; qz[m] += sz

corr, rows = fit_star_shifts(qx, qy, qz, sid, np.ones(n, bool), Zref, ground_of, grid,
                             np.ones((NY, NX), bool), min_cells=200, verbose=False)

print("SYNTHETIC RECOVERY -- the cloud was displaced by a KNOWN amount\n")
print(f"  {'swath':>6}{'applied shift (m)':>28}{'estimator returned (m)':>30}")
for s in sorted(SHIFT):
    a = SHIFT[s]; b = corr[s]
    print(f"  {s:>6}   ({a[0]:+.3f}, {a[1]:+.3f}, {a[2]:+.3f})        "
          f"({b[0]:+.3f}, {b[1]:+.3f}, {b[2]:+.3f})")

print("\n  ratio returned/applied, per component:")
for s in sorted(SHIFT):
    a = np.array(SHIFT[s]); b = np.array(corr[s])
    print(f"    swath {s}:  x {b[0]/a[0]:+.3f}   y {b[1]/a[1]:+.3f}   z {b[2]/a[2]:+.3f}")

print("\n  VERDICT -- every component must come back as MINUS the applied shift, so that")
print("  adding the returned value restores the cloud. Any '+' below is a sign bug:")
bad = []
for comp, i in (("x", 0), ("y", 1), ("z", 2)):
    r = np.mean([corr[s][i] / SHIFT[s][i] for s in SHIFT])
    good = r < 0
    bad.append(not good)
    print(f"    {comp}: returned/applied = {r:+.3f}  "
          f"{'OK, add to correct' if good else 'SIGN BUG -- adding it would DOUBLE the error'}")

print("\n  PROOF: apply the returned correction and re-fit; it must come back ~0")
rx, ry, rz = qx.copy(), qy.copy(), qz.copy()
for s, (cx, cy, cz) in corr.items():
    m = sid == s
    rx[m] += cx; ry[m] += cy; rz[m] += cz
corr2, _ = fit_star_shifts(rx, ry, rz, sid, np.ones(n, bool), Zref, ground_of, grid,
                           np.ones((NY, NX), bool), min_cells=200, verbose=False)
for s in sorted(corr2):
    b = corr2[s]
    print(f"    swath {s} residual after correction: "
          f"({1000*b[0]:+.1f}, {1000*b[1]:+.1f}, {1000*b[2]:+.1f}) mm")
resid_ok = all(abs(v) < 0.05 for s in corr2 for v in corr2[s])

mag_ok = all(abs(abs(corr[s][i] / SHIFT[s][i]) - 1.0) < 0.12 for s in SHIFT for i in range(3))
ok = mag_ok and not any(bad) and resid_ok
print(f"\n  magnitudes within 12%: {mag_ok}   signs all correct: {not any(bad)}   "
      f"residual < 50 mm: {resid_ok}")
print(f"  {'PASS' if ok else 'FAIL'}")
raise SystemExit(0 if ok else 1)
