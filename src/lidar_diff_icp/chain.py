"""The cross-epoch correction chain, as independent steps a driver runs in order.

Andy, 2026-09-11: "The pipeline should have each step be a Python module, function,
and/or class. A driver should then run them in order. This way, 'our' driver can be built
rapidly, as can a 'DeLong' driver, and from the same code base."

Every step here takes the same :class:`Ctx`, mutates the gen1 point arrays IN PLACE, and
writes what it did into ``ctx.record``. A driver is then just a list:

    OUR_CHAIN    = [LateralShift(), GeoidConversion(...), AlongTrackDrift()]
    DELONG_CHAIN = [LateralShift(), CorrectionSurface()]

and ``pipeline.ROUTES`` names the two. Steps compose because none of them reads another's
internals -- only ``Ctx`` and their own arguments.

ORDER IS PART OF THE METHOD, not a detail the driver may shuffle freely:

* :class:`LateralShift` must run FIRST. Get x, y right before z is touched, or terrain
  slope leaks into the elevation difference. Its own magnitude shows why: the lateral
  term is 10.25 mm mean on slope <= 3 deg and 212.60 mm on the steepest ground.
* :class:`CorrectionSurface` must precede :class:`AlongTrackDrift` if both run. A drift
  fit on a residual that still holds a spatial field absorbs that field PER LINE, and a
  position-only surface cannot undo a per-line offset once made. Measured: drift -> surface
  leaves the mean per-line step at 19.91 mm, surface -> drift at 5.49 mm, against 5.07 mm
  from ``align_swaths`` alone. :func:`run_chain` REFUSES the wrong order rather than
  silently producing the worse product.

MUTATION IS DELIBERATE. At statewide scale these are tens of millions of points and a copy
per step is real memory on a shared machine. The arrays a caller passes in ARE the arrays
that come out; copy first if the originals are needed.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import coreg

__all__ = ["Ctx", "Correction", "LateralShift", "GeoidConversion", "CorrectionSurface",
           "AlongTrackDrift", "run_chain"]


@dataclass
class Ctx:
    """Everything a correction step may read, and the record it writes into.

    ``x``, ``y``, ``z`` are the gen1 points and are MUTATED. ``Zref`` is the gen2 grid --
    every step except :class:`GeoidConversion` reads it, which is precisely why none of
    this is gen1-internal, whatever a step is named.
    """

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    ground: np.ndarray                 # bare-earth mask into x/y/z
    Zref: np.ndarray                   # gen2 gridded ground
    ground_of: object                  # the SAME estimator that produced Zref
    grid: tuple                        # (X0, Y0, res, nx, ny)
    bounds: tuple
    gps_time: np.ndarray | None = None
    source_id: np.ndarray | None = None
    stable: np.ndarray | None = None   # geometric stable mask, grid-shaped
    floodplain: np.ndarray | None = None
    verbose: bool = True
    record: dict = field(default_factory=dict)
    #: Grid-shaped arrays a step APPLIED, kept so the product can be audited. ``record``
    #: goes to corrections.json and cannot hold an array, so a correction that is applied
    #: and then discarded is invisible: on the DeLong route the correction surface is the
    #: term that replaces BOTH the geoid and the drift, and until 2026-09-17 nothing
    #: downstream could see it. Reconstructing it from the beam table does not work -- that
    #: table carries only the dz_* columns, so the surface is exactly what is missing.
    grids: dict = field(default_factory=dict)

    def cell_index(self):
        """Point -> cell indices, clipped to the grid. Several steps need this and it
        must be computed the same way in each, or they disagree by a cell at the edge."""
        X0, Y0, res, nx, ny = self.grid
        ix = np.clip(((self.x - X0) / res).astype(int), 0, nx - 1)
        iy = np.clip(((self.y - Y0) / res).astype(int), 0, ny - 1)
        return iy, ix

    def gen1_grid(self):
        """The gen1 ground gridded by the caller's estimator, as the steps see it now."""
        b = self.ground
        return self.ground_of(self.x[b], self.y[b], self.z[b])


class Correction:
    """One step. Subclasses set ``name`` and implement ``apply``."""

    name = "correction"

    def apply(self, ctx: Ctx) -> None:      # pragma: no cover - interface
        raise NotImplementedError

    def __repr__(self):
        return f"{type(self).__name__}()"


class LateralShift(Correction):
    """Nuth & Kaeaeb order-0: ONE constant (dx, dy) from the full topography.

    Order 0 is a single constant shift, NOT the removed order-2 parabola; only the
    horizontal part of the fit is used. Registering the whole DEM recovers it because
    drainage divides do not move to first order, and the aspect dipole it fits is
    erosion-robust since diffuse erosion has no aspect dependence.
    """

    name = "lateral_shift"

    def apply(self, ctx):
        X0, Y0, res, nx, ny = ctx.grid
        hs = coreg.tie_polynomial(ctx.Zref, ctx.gen1_grid(), res, X0, Y0, order=0)
        # y is shifted using the ALREADY-UPDATED x, as the inline code did. Kept exactly:
        # changing it would move every product by a sub-cell amount for no reason.
        ctx.x += coreg.eval_poly_field(hs["a"], ctx.x, ctx.y, hs["norm"], 0)
        ctx.y += coreg.eval_poly_field(hs["b"], ctx.x, ctx.y, hs["norm"], 0)
        ctx.record[self.name] = {
            "horizontal_shift_m": [round(float(hs["a"][0]), 4),
                                   round(float(hs["b"][0]), 4)]}


class GeoidConversion(Correction):
    """The geoid-model datum shift N_gen1 - N_gen2, e.g. GEOID03 -> GEOID18.

    A required geodetic offset whenever gen1's surface is built on its own terms, and
    auto-computed from the PROJ grids when not supplied -- no hard-coded constant.

    NOT in the DeLong chain. The correction surface registers gen1 onto gen2 directly and
    absorbs any smooth epoch offset: a 54.87 mm geoid ERROR moves the DoD by -54.870 mm
    without the surface and by 0.000 mm with it.
    """

    name = "geoid"

    def __init__(self, gen1_geoid=None, geoid_datum=None):
        self.gen1_geoid = gen1_geoid
        self.geoid_datum = geoid_datum

    def apply(self, ctx):
        from . import references
        # FRAME FIRST. This refusal used to be absent, and its absence shipped three sites
        # in the wrong vertical frame. It is the first statement so a wrong call costs
        # nothing -- the caller may have streamed a multi-gigabyte tile to get here.
        if self.geoid_datum is None and self.gen1_geoid is None:
            raise ValueError(
                "gen1_geoid is required: the PROJ geoid grid gen1 was reduced in. It used "
                "to default to GEOID03, and because nothing overrode it every site was "
                "differenced as if gen1 were the 2008 SE-Minnesota survey -- adding "
                "+54.87 mm to Battle Creek's gen1, +27.75 to Carlton's and +26.39 to "
                "Cook's, unnoticed until 2026-09-07. Resolve it from the Site: "
                "acquisitions.for_project(site.gen1_project).geoid_grid.")
        gd = self.geoid_datum
        if gd is None:
            gd = references.geoid_difference(ctx.bounds, 26915,
                                             before_geoid=self.gen1_geoid)
        gc, gb, gcc = gd     # (const_m, b East, c North), m and m/km of N_gen1 - N_gen2
        cx = 0.5 * (ctx.bounds[0] + ctx.bounds[2])
        cy = 0.5 * (ctx.bounds[1] + ctx.bounds[3])
        ctx.z += gc + gb * (ctx.x - cx) / 1000.0 + gcc * (ctx.y - cy) / 1000.0
        ctx.record[self.name] = {"const_m": gc, "tilt_b_m_per_km": gb,
                                 "tilt_c_m_per_km": gcc, "centroid": [cx, cy],
                                 "grid": self.gen1_geoid}
        if ctx.verbose:
            print(f"  geoid-difference datum: const {1000*gc:+.1f} mm, tilt "
                  f"({1000*gb:+.3f},{1000*gcc:+.3f}) mm/km", flush=True)


class CorrectionSurface(Correction):
    """DeLong et al. (2022) masked-stable IDW correction surface.

    Masks the cells that may hold real change -- slope > 3 deg, |dz| > 0.7 m, and whatever
    ``ctx.floodplain`` excludes -- then interpolates the remaining stable residual and
    subtracts it. It is a function of GRID POSITION ONLY, so two flight lines sharing a
    cell get the identical correction: it provably cannot change a per-line step, and
    ``align_swaths`` remains the only thing that can.
    """

    name = "correction_surface"

    def __init__(self, radius_m=400.0):
        self.radius_m = radius_m

    def apply(self, ctx):
        X0, Y0, res, nx, ny = ctx.grid
        C = coreg.correction_surface(ctx.Zref, ctx.gen1_grid(), res, X0, Y0,
                                     radius=self.radius_m, exclude=ctx.floodplain)["C"]
        iy, ix = ctx.cell_index()
        Cpt = C[iy, ix]
        good = np.isfinite(Cpt)
        ctx.z[good] += Cpt[good]
        ctx.grids[self.name] = C
        ctx.record[self.name] = {"radius_m": self.radius_m,
                                 "points_corrected": int(good.sum()),
                                 "points_total": int(good.size)}


class AlongTrackDrift(Correction):
    """Per-swath vertical drift as a smooth function of ``gps_time``.

    Motivated as GNSS trajectory drift, but the fit is empirical: it absorbs ANYTHING that
    varies smoothly along-track on stable ground, including a datum gradient (see task
    #63). Measured to make per-line agreement WORSE -- align_swaths takes the mean per-line
    step to 5.07 mm and this returns it to ~13-21 mm -- so it is not in the default chain.
    """

    name = "along_track_drift"

    def apply(self, ctx):
        if ctx.gps_time is None or ctx.source_id is None or ctx.stable is None:
            raise ValueError("AlongTrackDrift needs gps_time, source_id and stable; one "
                             "is None. It is a per-swath, time-domain fit on stable "
                             "ground and cannot be run without all three.")
        iy, ix = ctx.cell_index()
        resid = ctx.Zref - ctx.gen1_grid()
        chg = resid[iy, ix]
        stab_pt = (ctx.ground & ctx.stable[iy, ix] & np.isfinite(chg)
                   & (np.abs(chg) < 0.15))
        drift, curves = coreg.fit_along_track_drift(ctx.gps_time, chg, stab_pt,
                                                    ctx.source_id)
        ctx.z += drift
        ctx.record[self.name] = {"curves": curves,
                                 "n_stable_points": int(stab_pt.sum())}


def run_chain(chain, ctx: Ctx) -> Ctx:
    """Run steps in order, refusing an order the measurements rule out.

    The one ordering constraint that is not taste: a drift fit on a residual that still
    holds a spatial field converts that field into per-line offsets, and a position-only
    surface cannot undo them afterwards. Refused rather than warned, because the result is
    a quietly worse product rather than an error.
    """
    names = [s.name for s in chain]
    if "along_track_drift" in names and "correction_surface" in names:
        if names.index("along_track_drift") < names.index("correction_surface"):
            raise ValueError(
                "AlongTrackDrift is placed BEFORE CorrectionSurface. A drift fit on a "
                "residual that still holds a spatial field absorbs that field per line, "
                "and a spatial surface cannot undo a per-line offset once made: measured "
                "mean per-line step 19.91 mm that way against 5.49 mm the other, on a "
                "5.07 mm baseline from align_swaths. Put CorrectionSurface first.")
    if names and names[0] != "lateral_shift":
        raise ValueError(
            f"the chain starts with {names[0]!r}, not 'lateral_shift'. Get x, y right "
            f"before z is touched, or terrain slope leaks into the elevation difference "
            f"(the lateral term is 212.60 mm mean on the steepest ground).")
    for step in chain:
        step.apply(ctx)
    ctx.record["chain"] = names
    return ctx


# ============================ THE GEN2 CHAIN =================================
# Andy, 2026-09-11: "we alter gen2 via the vegetation correction ... before fitting gen1
# to it. The DeLong step is the same. It is just matching a different gen2."
#
# That sentence is the whole design. gen1's chain above is unchanged; what changes is the
# TARGET it is registered onto. So gen2 needs a chain of its own, run BEFORE the gen1 one.
#
# WHY THIS IS A NEW SEAM AND NOT A FLAG. Until now gen2's only correction ran AFTER
# apply_datum (pipeline.correct_reference, called at the end), so gen1 was always fitted to
# the UNCORRECTED gen2 and the correction was applied to the finished difference. Reversing
# that is an ORDER change, which no argument could express. The rule recorded in
# pipeline's architecture note was: split when you can name the second implementation that
# would use the seam. Here it is named -- median gen2 against cover-percentile gen2 -- so
# the seam is earned rather than speculative.

__all__ += ["Gen2Ctx", "Gen2Correction", "SpreadPercentile", "run_gen2_chain"]


@dataclass
class Gen2Ctx:
    """The gen2 grid a step may replace, and what it needs to do so.

    ``Zref`` is the gen2 ground the gen1 chain will be registered onto. A step REPLACES it
    (rebinds, not in-place) because a percentile change is a new surface, not a nudge.
    """

    Zref: np.ndarray
    Z21: np.ndarray                    # the q=0.50 grid; the reference the columns hang off
    after_laz: str
    grid: tuple                        # (X0, Y0, res, nx, ny)
    tile_dir: str | None = None
    verbose: bool = True
    record: dict = field(default_factory=dict)
    grids: dict = field(default_factory=dict)   # per-cell diagnostics, for corrections.json


class Gen2Correction:
    """One gen2-side step. Subclasses set ``name`` and implement ``apply``."""

    name = "gen2_correction"

    def apply(self, ctx: Gen2Ctx) -> None:       # pragma: no cover - interface
        raise NotImplementedError

    def __repr__(self):
        return f"{type(self).__name__}()"


class SpreadPercentile(Gen2Correction):
    """Read gen2's ground at a percentile set by the cell's own CLASS-2 SPREAD.

    THIS IS THE ADOPTED gen2 MOVER. ``difference_dem`` reaches it as
    ``ground_q="calibrated"`` with ``gen2_curve=<curve .npz>``; the curve is an isotonic
    regression fitted from SURVEYED CONTROL MARKS (:mod:`groundq`), so unlike the per-tile
    cover fit it transfers between tiles.

    NO DEFAULT CURVE, deliberately. Which point types were fitted decides what the curve
    means, and groundq records the hazard: pooled over all 519 marks the curve looked like
    a 16% RMS gain, but the falling limb is entirely the 162 VVA marks, which are sited
    UNDER VEGETATION by design. Applied to ordinary ground it widened the DoD's scatter at
    both sites tested (Elba NMAD 74.8 -> 79.1 mm, Whitewater 85.0 -> 92.4). On open ground
    the NVA-only curve measured WORSE than the plain median, RMS 52.5 against 49.1 held
    out. So naming the curve is a scientific choice and this step will not make it.
    """

    name = "spread_percentile"

    def __init__(self, curve, *, chunk=8_000_000):
        self.curve = curve
        self.chunk = chunk

    def apply(self, ctx):
        from . import groundq
        X0, Y0, res, nx, ny = ctx.grid
        surf = groundq.surface_from_grid(ctx.Z21, X0, Y0, res)
        gq = groundq.correct_gen2(ctx.after_laz, self.curve, surf=surf, chunk=self.chunk)
        nn = surf["nnorm"].reshape(ny, nx)
        dv = (gq["h2_mm"] - gq["h2_median_mm"]).reshape(ny, nx) / 1000.0 * nn
        declined = int((np.isfinite(ctx.Zref) & ~np.isfinite(dv)).sum())
        ctx.Zref = ctx.Zref + dv
        ctx.grids["gen2_ground_q"] = gq["q"].reshape(ny, nx)
        ctx.grids["gen2_class2_sd_mm"] = gq["sd_mm"].reshape(ny, nx)
        ctx.grids["gen2_correction_m"] = dv
        ctx.record[self.name] = {
            "curve": str(self.curve),
            "median_correction_mm": float(np.nanmedian(dv) * 1000.0),
            "cells_lowered": int(np.sum(dv < -0.001)), "cells_declined": declined}
        if ctx.verbose:
            print(f"  gen2 spread percentile: median correction "
                  f"{np.nanmedian(dv)*1000:+.1f} mm; {int(np.sum(dv < -0.001)):,} cells "
                  f"lowered; {declined:,} DECLINED", flush=True)


def run_gen2_chain(chain, ctx: Gen2Ctx) -> Gen2Ctx:
    """Run gen2 steps in order, then hand ``ctx.Zref`` to the gen1 chain.

    ONLY ``SpreadPercentile`` belongs here. Andy, 2026-09-11: "The ground_q method is the
    only one that should be there." A CoverPercentile step was briefly added and is gone:
    FRAME says the cover route's outputs are "dod_cover_q2.npy / lod_cover_q2.npy, NEVER
    dod.npy / lod.npy", and anything in this chain moves Zref and therefore dod.npy. The
    cover route keeps its proper home in the alongside group
    (analysis/modules/vegetation_correction/dod_cover_corrected.py), which applies it to a
    finished DoD rather than to the surface gen1 is registered onto.

    The chain exists for ORDER, not for variety: gen2 corrected BEFORE gen1 is registered
    onto it. That was unreachable while the only gen2 correction ran after apply_datum, and
    the order is declared here instead of being implied by where a call sits among three
    hundred lines -- which is how the old post-hoc order went unnoticed.
    """
    for step in chain:
        step.apply(ctx)
    ctx.record["gen2_chain"] = [s.name for s in chain]
    return ctx
