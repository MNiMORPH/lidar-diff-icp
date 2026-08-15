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
