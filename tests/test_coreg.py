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


if __name__ == "__main__":
    test_recovers_known_shift()
    test_zero_shift_is_zero()
    print("coreg regression tests PASS")
