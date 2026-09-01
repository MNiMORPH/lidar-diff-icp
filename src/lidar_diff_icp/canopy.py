"""Canopy-closure / ground-penetration flag for leaf-on acquisitions.

Where a gen2 (newer) survey was flown under leaf-on / green-up canopy, its ground
returns are STARVED (canopy blocks the laser), and on slopes the sparse ground reads
high -- a bare-earth bias that grows with slope and lives under CLOSED canopy (verified
at Elba: gen2 flown 2021-05-01 green-up; DoD bias +tan(slope) concentrated where ground
penetration is low, absent on bare/open ground). This is a data/acquisition limit, not
a griddable artifact, so we FLAG it (widen the level-of-detection) rather than correct.

The covariate is per-cell GROUND PENETRATION = ground returns / total returns. Low
penetration on a slope = leaf-on-under-forest = unreliable bare-earth. See
analysis/slope_bias/ and mn-3dep-audit (which flags the acquisition itself from its date).

STATUS OF THAT MECHANISM -- read before citing it (2026-09-01)
--------------------------------------------------------------
The flag itself stands: it is a conservative widening, and the OBSERVATION behind it -- a
slope-dependent DoD bias concentrated in forest -- has held. The MECHANISM in the first
paragraph, gen2 ground starvation under leaf-on canopy, has not:

* `analysis/slope_bias/README.md` has read OPEN since 2026-08-21: the best-motivated
  mechanism "is not proven; a competing vegetation-penetration mechanism is now in play."
* gen2 is not ground-starved in absolute terms. `analysis/gen1_grass_lift.py` records
  5.78 gen2 against 0.49 gen1 ground returns per m^2 -- about 12x MORE ground, not less.
* `analysis/ridgelines/FRAME_2026-08-26.md` finds the offset "tracks the 2008 canopy, not
  2021 cover", with phenology refuted as the driver, which points the mechanism at gen1
  rather than at gen2 leaf-on.
* `analysis/ridgelines/AUDIT_findings.md` puts `penetration.npy` in Tier 1 as a
  gen2-derived variable and warns against binning any gen1 quantity against it;
  `gen1_own_penetration.py` is the gen1-internal counterpart.

So: use this to WIDEN a detection limit, which is what it does. Do not cite it as the
explanation of the slope bias, and do not use `penetration` as a canopy measure --
`canopy_cover_pfs` is the cover measure, computed identically on every tile.
"""
from __future__ import annotations

import numpy as np


def ground_penetration(after_laz, bounds, res, nx, ny, *, ground_class=2, noise_class=7,
                       chunk_points=20_000_000):
    """Per-cell fraction of returns that reached the ground (class ``ground_class`` /
    total non-noise). Low = closed/leaf-on canopy.

    A cell with NO returns is ``nan``, NOT ``0.0``: zero means "measured, nothing reached
    the ground", which every downstream cut of the form ``pen < 0.25`` reads as maximally
    closed canopy. ``tests/test_canopy.py`` pins the distinction.

    Read in chunks rather than whole. ``chunk_points`` is a MEMORY setting and cannot
    change the answer: the per-cell counts are accumulated with ``bincount``, which is
    additive, and the single division happens once at the end -- so any chunk size gives
    bit-identical output. It matters because a full-density tile is ~1.8e8 points.
    """
    import laspy
    X0, Y0, X1, Y1 = bounds
    tot = np.zeros(nx * ny, float)
    gnd = np.zeros(nx * ny, float)
    with laspy.open(str(after_laz)) as fh:
        for p in fh.chunk_iterator(chunk_points):
            x = np.asarray(p.x); y = np.asarray(p.y); cl = np.asarray(p.classification)
            m = (x >= X0) & (x < X1) & (y >= Y0) & (y < Y1) & (cl != noise_class)
            cid = ((y[m] - Y0) / res).astype(int) * nx + ((x[m] - X0) / res).astype(int)
            tot += np.bincount(cid, minlength=nx * ny)
            gnd += np.bincount(cid[cl[m] == ground_class], minlength=nx * ny)
    frac = np.where(tot > 0, gnd / np.maximum(tot, 1), np.nan)
    return frac.reshape(ny, nx)


def leafon_slope_flag(penetration, slope_deg, *, min_penetration=0.25, min_slope=12.0):
    """Boolean flag: cells where ground penetration is poor AND the slope is steep
    enough that sparse ground biases the surface -- the leaf-on-under-forest zone whose
    bare-earth DoD should be treated as low-confidence.

    An unmeasured cell (``penetration`` NaN) is NOT flagged: absent data is not evidence
    of closed canopy.

    BOTH thresholds are conventions, not derivations. ``min_slope=12.0`` is the repo-wide
    gentle-ground cut -- the same 12 deg used by ``refcells.reference_cells`` and by the
    cover and near-ground drivers -- so it is at least consistent, though no source in
    this repo derives it. ``min_penetration=0.25`` is the same 0.25 as the forest side of
    the strata cut (``pen < 0.25`` forest, ``>= 0.45`` open); whether that reuse is
    deliberate is not recorded. Pass them explicitly if a site needs a different bar.
    """
    return (np.isfinite(penetration) & (penetration < min_penetration)
            & (slope_deg > min_slope))


def inflate_lod(lod, flag, *, factor=2.0):
    """Widen the level-of-detection on flagged (leaf-on-forest-slope) cells so change
    there is held to a higher bar. Returns a copy; flagged cells' LoD *= factor.

    ``factor=2.0`` HAS NO DERIVATION anywhere in this repo -- searched. It is a round
    number, and it sets what counts as detected change on every flagged cell at every
    site ``scripts/run_all_sites.py`` runs, since that caller does not override it. Doubling
    the detection limit roughly quarters the area that can register significant change
    there, so this is a load-bearing choice presented as a default. It should either be
    calibrated against measured stable-ground scatter under closed canopy, or passed
    explicitly per site with the value stated in the run's output.
    """
    out = lod.copy()
    out[flag] = out[flag] * factor
    return out
