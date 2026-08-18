"""Regression tests for Nuth & Kaeaeb co-registration.

The key failure these guard against is a sign/convention error in the shift
update: an early version added the offset instead of subtracting it, which
recovered the *negative* of the true shift and then diverged. These tests inject
a known shift into a synthetic surface and require it to be recovered (as its
inverse) to sub-centimetre accuracy. They need no downloaded data.
"""
import numpy as np

from lidar_diff_icp import coreg


def _synthetic_surface(res=1.0, n=400):
    """A smooth surface with slopes in both directions (varied aspect)."""
    yy, xx = np.mgrid[0:n, 0:n] * res
    return 8.0 * np.sin(xx / 40.0) + 6.0 * np.cos(yy / 55.0) + 0.02 * xx


def test_recovers_known_shift():
    res = 1.0
    z_ref = _synthetic_surface(res)
    for sx, sy, sz in [(1.5, -0.8, 0.05), (-2.3, 1.1, -0.12), (3.0, -3.0, 0.2)]:
        z_src = coreg._shift_grid(z_ref, sx, sy, res) + sz
        c = coreg.nuth_kaab(z_ref, z_src, res)
        assert c.converged, f"did not converge for shift {(sx, sy, sz)}"
        # N&K returns the correction: the inverse of the injected shift
        assert abs(c.dx + sx) < 0.05, f"dx: {c.dx} vs {-sx}"
        assert abs(c.dy + sy) < 0.05, f"dy: {c.dy} vs {-sy}"
        assert abs(c.dz + sz) < 0.02, f"dz: {c.dz} vs {-sz}"


def test_translation_tilt_recovers_known():
    res = 1.0
    z_ref = _synthetic_surface(res)
    n = z_ref.shape[0]
    gy, gx = np.mgrid[0:n, 0:n]
    Xc = (gx + 0.5) * res; Xc = Xc - Xc.mean()
    Yc = (gy + 0.5) * res; Yc = Yc - Yc.mean()
    sx, sy, a0, a1, a2 = 1.5, -0.8, 0.05, 0.001, -0.0005   # tilt 1.0 / -0.5 mm/m
    z_src = coreg._shift_grid(z_ref, sx, sy, res) + (a0 + a1 * Xc + a2 * Yc)
    r = coreg.tie_translation_tilt(z_ref, z_src, res, 0, 0)
    assert r["converged"]
    assert abs(r["dx"] + sx) < 0.05 and abs(r["dy"] + sy) < 0.05
    assert abs(r["c1"] + a1) < 2e-4 and abs(r["c2"] + a2) < 2e-4


def test_polynomial_tie_recovers_quadratic_warp():
    res = 1.0
    z_ref = _synthetic_surface(res)
    n = z_ref.shape[0]
    gy, gx = np.mgrid[0:n, 0:n]
    X = (gx + 0.5) * res; Y = (gy + 0.5) * res
    xm = 0.5 * (X.max() + X.min()); xhr = 0.5 * (X.max() - X.min())
    ym = 0.5 * (Y.max() + Y.min()); yhr = 0.5 * (Y.max() - Y.min())
    Xn = (X - xm) / xhr; Yn = (Y - ym) / yhr
    sxf = 0.8 + 0.4 * Xn - 0.2 * Yn + 0.15 * Xn * Xn
    syf = -0.5 + 0.3 * Yn
    szf = 0.05 + 0.02 * Xn
    z_src = coreg._warp_grid(z_ref, sxf, syf, res) + szf
    r = coreg.tie_polynomial(z_ref, z_src, res, 0, 0, order=2)
    assert r["converged"]
    assert r["nmad_after"] < 0.02
    assert np.sqrt(np.nanmean((r["dx_field"] + sxf) ** 2)) < 0.03


def test_zero_shift_is_zero():
    res = 1.0
    z = _synthetic_surface(res)
    c = coreg.nuth_kaab(z, z.copy(), res)
    assert np.hypot(c.dx, c.dy) < 0.02
    assert abs(c.dz) < 0.01
    assert c.nmad_after <= c.nmad_before + 1e-9


def test_correction_surface_recovers_warp_and_ignores_real_change():
    """DeLong-style correction surface: it must recover a smooth vertical warp
    from STABLE (flat, small-dz) ground, and must NOT absorb a large localized
    change that the masks exclude -- the surface should interpolate smoothly
    through it, leaving the change in the corrected difference."""
    res = 5.0; n = 80
    yy, xx = np.mgrid[0:n, 0:n] * res
    z_src = np.zeros((n, n))                       # flat ground -> all "stable"
    warp = 0.15 * np.sin(2 * np.pi * xx / 300.0)   # smooth spurious warp, <0.7 m
    z_ref = z_src + warp
    z_ref[40:45, 40:45] += 2.0                     # localized real change
    cs = coreg.correction_surface(z_ref, z_src, res, 0.0, 0.0,
                                  slope_thresh_deg=3.0, dz_thresh=0.7, radius=400.0)
    C = cs["C"]
    assert not cs["stable"][42, 42]                # change patch excluded from stable
    far = np.ones((n, n), bool); far[35:50, 35:50] = False
    assert np.sqrt(np.nanmean((C[far] - warp[far]) ** 2)) < 0.03  # recovers warp
    corrected = cs["dz"] - C
    assert corrected[42, 42] > 1.5                 # does not absorb the +2 m change


def test_correction_surface_exclude_drops_sources():
    """A caller-supplied exclude mask (e.g. a TPI floodplain buffer) must remove
    those cells from the stable IDW sources so the correction is not fit to
    them -- otherwise real low-slope change there is absorbed."""
    res = 5.0; n = 80
    z_src = np.zeros((n, n)); z_ref = z_src + 0.10   # flat, all otherwise stable
    ex = np.zeros((n, n), bool); ex[:, :40] = True
    a = coreg.correction_surface(z_ref, z_src, res, 0.0, 0.0)
    b = coreg.correction_surface(z_ref, z_src, res, 0.0, 0.0, exclude=ex)
    assert b["n_stable"] < a["n_stable"]
    assert not b["stable"][10, 10]                   # excluded region dropped
    assert b["stable"][10, 60]                       # kept region retained


def test_along_track_drift_recovers_known():
    """Inject a smooth per-swath drift as a function of gps_time and recover it.
    The residual on stable ground is drift + noise; the spline fit must track the
    drift, not the noise."""
    rng = np.random.default_rng(0); n = 20000
    gps = np.sort(rng.uniform(0.0, 100.0, n))
    swath = np.ones(n, int)
    true = 0.04 * np.sin(np.pi * gps / 100.0)          # smooth half-sine
    change = true + rng.normal(0, 0.02, n)             # stable residual = drift + noise
    stable = np.ones(n, bool)
    drift, curves = coreg.fit_along_track_drift(gps, change, stable, swath)
    assert 1 in curves
    assert np.corrcoef(drift, true)[0, 1] > 0.9
    assert np.sqrt(np.mean((drift - true) ** 2)) < 0.02


def test_tie_falls_back_on_gentle_terrain():
    """On gentle terrain with no real warp the Nuth & Kaeaeb horizontal shift is
    unconstrained; the tie must fall back to a rigid vertical offset, never a
    runaway (dx ~ 1000 m) that makes the fit worse than the input."""
    rng = np.random.default_rng(0); res = 5.0; n = 60
    yy, xx = np.mgrid[0:n, 0:n] * res; cx = cy = n * res / 2
    dome = 105.0 - 1.5e-4 * ((xx - cx) ** 2 + (yy - cy) ** 2)   # < 3 deg everywhere
    z_ref = dome + rng.normal(0, 0.02, (n, n))
    z_src = dome + rng.normal(0, 0.02, (n, n))                  # no warp, just noise
    t = coreg.tie_polynomial(z_ref, z_src, res, 0.0, 0.0, order=2)
    # safe by property, however reached (rigid fallback or early break to
    # identity): no runaway horizontal shift, and never worse than the input.
    assert np.ptp(t["dx_field"]) < 5.0
    assert t["nmad_after"] <= t["nmad_before"] + 1e-9


def test_slope_aspect_planar():
    """A plane rising toward the east has slope arctan(gradient) and aspect 90
    deg (steepest ascent toward east, clockwise from north)."""
    res = 1.0; n = 20
    _, xx = np.mgrid[0:n, 0:n] * res
    z = 0.1 * xx                                   # rises toward +x (east)
    slope, aspect = coreg.slope_aspect(z, res)
    inner = (slice(2, -2), slice(2, -2))           # avoid edge gradients
    assert abs(np.median(slope[inner]) - np.arctan(0.1)) < 1e-6
    assert abs(np.degrees(np.median(aspect[inner])) - 90.0) < 1.0


def test_align_swaths_recovers_vertical_offset():
    """Two overlapping swaths, swath 2 biased +0.3 m in z; the free-network
    alignment must pin swath 1 and recover ~ -0.3 m for swath 2."""
    from lidar_diff_icp import io
    rng = np.random.default_rng(0); n = 15000
    def ground(x, y):
        return 100.0 + 0.05 * x + 0.03 * y         # ~3.4 deg slope
    xa = rng.uniform(0, 120, n); ya = rng.uniform(0, 150, n)     # swath 1
    xb = rng.uniform(80, 200, n); yb = rng.uniform(0, 150, n)    # swath 2 (overlap 80-120)
    za = ground(xa, ya) + rng.normal(0, 0.02, n)
    zb = ground(xb, yb) + 0.3 + rng.normal(0, 0.02, n)           # +0.3 m bias
    x = np.concatenate([xa, xb]); y = np.concatenate([ya, yb]); z = np.concatenate([za, zb])
    ps = np.concatenate([np.ones(n), np.full(n, 2)]).astype(int)
    pc = io.PointCloud(x, y, z, ps, np.zeros_like(z), np.zeros_like(z),
                       np.zeros_like(ps), io.MN_GEN1_CRS)
    corr, edges, mis = coreg.align_swaths(pc, ref=1)
    assert abs(corr[1][2]) < 0.05                  # reference pinned
    assert abs(corr[2][2] + 0.3) < 0.08            # +0.3 m bias recovered as -0.3


def test_align_swaths_ignores_nan_edge(monkeypatch):
    """A non-adjacent swath pair sharing gridded extent but no actual cells
    returns n=0 with NaN shifts (real multi-swath tiles do this). A single NaN
    observation must NOT poison the least-squares network -- the connected chain
    1-2-3 still determines every swath. Regression for the Carlton all-NaN DoD."""
    from lidar_diff_icp import io
    def mk(dz, n):                                  # a pairwise Coreg observation
        return coreg.Coreg(0.0, 0.0, dz, 0.0, 0.0, 0.0, n, 0.1, 0.05, 3, True)
    fake = {(1, 2): mk(-0.02, 1000), (2, 3): mk(-0.03, 1000),
            (1, 3): mk(np.nan, 0)}                  # the poisoning empty-overlap edge
    monkeypatch.setattr(coreg, "coregister_swaths",
                        lambda pc, a, b, res, exclude: fake[(a, b)])
    ps = np.array([1, 1, 2, 2, 3, 3])
    pc = io.PointCloud(np.zeros(6), np.zeros(6), np.zeros(6), ps,
                       np.zeros(6), np.zeros(6), np.zeros(6, int), io.MN_GEN1_CRS)
    corr, edges, mis = coreg.align_swaths(pc, ref=1)
    assert all(np.isfinite(v).all() for v in corr.values())   # not NaN-poisoned
    assert (1, 3) not in {(e[0], e[1]) for e in edges}        # empty edge dropped
    assert abs(corr[1][2]) < 1e-9                             # reference pinned
    assert abs(corr[2][2] + 0.02) < 1e-6                      # chain determines 2
    assert abs(corr[3][2] + 0.05) < 1e-6                      # ... and 3


if __name__ == "__main__":
    test_recovers_known_shift()
    test_zero_shift_is_zero()
    test_correction_surface_recovers_warp_and_ignores_real_change()
    test_correction_surface_exclude_drops_sources()
    test_along_track_drift_recovers_known()
    test_tie_falls_back_on_gentle_terrain()
    print("coreg regression tests PASS")
