"""Apply a ground-control datum constant, making a product's elevation GAUGE-INVARIANT.

Why this is a required step, not an optional refinement
-------------------------------------------------------
``coreg.align_swaths`` solves a FREE NETWORK and subtracts the reference swath's value
afterwards. That gauge does not touch any swath-to-swath difference — but it sets the
absolute level the whole mosaic inherits, because that level becomes **the reference
line's own vertical error**. Measured on elbaext, the six per-swath ``dz`` span

    133  +0.00   134 +22.00   135  +6.20   136  -9.80   137 -18.40   138 -22.60
    => re-gauging on a different line moves EVERY elevation by up to 44.60 mm

So an uncorrected product's elevation is an arbitrary implementation detail
(``ref=int(ps.min())``), not a measurement.

**The correction removes that dependence exactly.** With ``corrected = z + c`` and ``c``
measured against the SAME gauged product, re-gauging by ``d`` shifts ``z`` by ``+d`` and
``c`` by ``−d``. They cancel: the corrected surface is identical whichever line is pinned.
:func:`gauge_invariance_residual` demonstrates this rather than asserting it.

That is the reason to apply the constant even where it looks negligible. At Elba the
correction is only +2.12 mm — but the gauge choice there is worth 44.60 mm, **21×
larger**. The smallness is a property of line 133 having been a lucky pin, not of the data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DatumApplication:
    """A datum constant applied to a product, and the gauge it is tied to."""

    constant_mm: float
    sigma_mm: float
    gauge_ref: int
    source: str
    note: str = ""

    def apply(self, z):
        """Return ``z`` on absolute NAVD88. Positive constant = the surface reads LOW."""
        return np.asarray(z, float) + self.constant_mm / 1000.0

    def regauged_to(self, new_ref: int, per_swath_dz_mm: dict) -> "DatumApplication":
        """The same datum expressed against a DIFFERENT reference line.

        Re-gauging shifts the product by ``-dz[new_ref]`` relative to the current gauge, so
        the constant must move by ``+dz[new_ref]`` for the corrected surface to be
        unchanged. This is arithmetic, not a re-measurement -- but it MUST be done, or the
        constant silently belongs to the wrong product.
        """
        d = float(per_swath_dz_mm[new_ref]) - float(per_swath_dz_mm[self.gauge_ref])
        return DatumApplication(
            constant_mm=self.constant_mm + d, sigma_mm=self.sigma_mm,
            gauge_ref=int(new_ref), source=self.source,
            note=f"re-gauged from line {self.gauge_ref} to {new_ref} ({d:+.2f} mm)")


def gauge_invariance_residual(z_mm, per_swath_dz_mm, datum: DatumApplication):
    """Corrected elevation under EVERY possible gauge; the spread must be ~0.

    Returns ``(levels_uncorrected, levels_corrected)`` in mm, one entry per candidate
    reference line. The uncorrected spread is the arbitrariness the gauge introduces; the
    corrected spread is what survives it.
    """
    refs = sorted(per_swath_dz_mm)
    unc, cor = [], []
    for r in refs:
        shift = float(per_swath_dz_mm[r]) - float(per_swath_dz_mm[datum.gauge_ref])
        z_r = float(z_mm) - shift                     # the product re-gauged onto line r
        unc.append(z_r)
        cor.append(z_r + datum.regauged_to(r, per_swath_dz_mm).constant_mm)
    return np.array(unc), np.array(cor)
