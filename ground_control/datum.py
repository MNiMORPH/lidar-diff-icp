"""The datum constant at a site, for one epoch, from that epoch's OWN control.

The rule this module exists to enforce, and the reason it was written:

    **2008 control for gen1.  2021 control for gen2.  No cross-epoch control.**

The artifact this replaces, ``data/derived/elba_fulldensity/z_before_absolute.{npy,json}``,
ties gen1 at two *2021* marks chained through five and six gen1 flight-line links.  That
tie therefore carries thirteen years of real ground change at a road shoulder -- the very
quantity the DoD exists to measure -- plus the chain error that imported swath constants
incur.  It was not a design decision: git shows the product was written at 15:53 on
2026-08-26 and gen1's own 2008 control entered the repo at 17:17, eighty-three minutes
later.  It used the only control that existed at the time.

Sign convention throughout: ``tie = surveyed - z_lidar``.  **POSITIVE = the surface reads
LOW**, so the constant is what you ADD.  Verified per row by
:func:`control.verify_sign_convention` rather than inherited.

WHAT THE UNCERTAINTY IS THE UNCERTAINTY OF
------------------------------------------
Never "the SE".  Each estimate names its own quantity:

``sd_field_mm``
    the error in predicting *the spatially correlated component of the delivered
    surface's offset at this coordinate*, nugget filtered out as uncorrelated noise.
    This is the systematic offset of the surface at the site.  It is the right number
    for a datum constant.
``sd_new_mark_mm``
    the error in predicting *what a single new control mark placed here would read*.
    Includes the nugget, so it carries micro-scale, mark-siting and survey variance.
    Larger, and NOT the datum's uncertainty.

Neither is an SE of a sample mean.  The sample mean is the wrong estimator here: its SE
falls with n while the mean itself moves further, because the quantity varies spatially.

WHAT IS DELIBERATELY NOT DECIDED HERE
-------------------------------------
* **Cover treatment.**  It is the largest single lever in the problem -- adjusting gen1
  on open ground rather than on all marks pooled moves the answer by 57.20 mm, and gen2
  by 38.40 mm.  :func:`datum_at_site` takes ONE treatment and reports it;
  :func:`sweep_treatments` runs several so the spread is visible.  Nothing here picks one.
* **Which gen2 surface.**  gen2 publishes four (delivered DEM and delivered cloud, at two
  quality-level blocks) and they are four different answers; ``surface=`` is required.
* **The variogram range.**  The control does not determine it -- fits pin to the largest
  lag centre.  Every estimate is reported over a SWEEP of ``max_lag`` and estimator, never
  a single fit.  A fitted range with no empirical variogram under it is unfalsifiable.

WHAT THIS DOES NOT INCLUDE
--------------------------
The **bridge** from the delivered surface the control measures to our CSF-reprocessed,
swath-aligned, geoid-shifted reconstruction.  It is not zero and must be added
separately; ``DatumEstimate.bridge_mm`` is where it goes, and it is None until supplied.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, asdict, field
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "src"))

import control  # noqa: E402
from lidar_diff_icp.groundtruth import residual_field as RF  # noqa: E402


#: Cover treatments, per epoch.  ``marks`` selects which classes enter the field;
#: ``predict_at`` names the class the answer is FOR.  ``drift`` puts cover in the
#: universal-kriging drift instead of filtering on it.
TREATMENTS = {
    "gen1": {
        "open":            dict(marks=("L1O",), predict_at="L1O", drift=False),
        "open+urban":      dict(marks=("L1O", "L5U"), predict_at="L1O", drift=False),
        "cover_covariate": dict(marks=("L1O", "L2T", "L3B", "L4F", "L5U"),
                                predict_at="L1O", drift=True),
    },
    "gen2": {
        "open":            dict(marks=("NVA",), predict_at="NVA", drift=False),
        "cover_covariate": dict(marks=("NVA", "VVA"), predict_at="NVA", drift=True),
    },
}


@dataclass(frozen=True)
class VariogramRow:
    """One row of the sweep: a fit, and the value it produced."""

    max_lag_m: float
    n_lags: int
    estimator: str
    nugget: float
    sill: float
    range_m: float
    range_pinned_to_largest_lag: bool
    value_mm: float
    sd_field_mm: float
    sd_new_mark_mm: float


@dataclass
class DatumEstimate:
    """The constant to ADD to this epoch's surface at this site, and what it rests on."""

    epoch: str
    route: str
    site_easting: float
    site_northing: float
    crs: str
    surface: str | None
    cover_treatment: str
    predict_at_cover: str
    n_marks: int

    constant_mm: float                 # median over the variogram sweep
    constant_min_mm: float
    constant_max_mm: float
    sd_field_mm: float                 # median over the sweep
    sd_field_min_mm: float             # -- reported because the nugget/sill split that
    sd_field_max_mm: float             #    sets it is NOT identified by this control
    sd_new_mark_mm: float
    sd_is_of: str

    sweep: list = field(default_factory=list)
    bridge_mm: float | None = None
    bridge_sd_mm: float | None = None
    notes: list = field(default_factory=list)

    @property
    def constant_with_bridge_mm(self):
        """The constant carried onto OUR reconstructed surface, or None."""
        return None if self.bridge_mm is None else self.constant_mm + self.bridge_mm

    def to_dict(self):
        d = asdict(self)
        d["constant_with_bridge_mm"] = self.constant_with_bridge_mm
        d["sign_convention"] = (
            "tie = surveyed - z_lidar; POSITIVE = the surface reads LOW, "
            "so this constant is what you ADD to the surface")
        return d


def datum_at_site(epoch: str, *, easting: float, northing: float,
                  treatment: str, surface: str | None = None,
                  max_lags_m, n_lags: int, n_pairs: int, estimators, seed: int,
                  bridge_mm: float | None = None,
                  bridge_sd_mm: float | None = None) -> DatumEstimate:
    """Predict the datum constant AT a site from this epoch's own published residuals.

    Every sweep parameter is required.  A single ``max_lag`` would privilege one fitted
    range, and the control does not determine the range.
    """
    if epoch not in TREATMENTS:
        raise ValueError(f"epoch={epoch!r} must be one of {tuple(TREATMENTS)}")
    if treatment not in TREATMENTS[epoch]:
        raise ValueError(
            f"treatment={treatment!r} unknown for {epoch}; "
            f"known: {tuple(TREATMENTS[epoch])}")
    spec = TREATMENTS[epoch][treatment]

    load = control.load_control(epoch, surface=surface)
    cr = load.residuals
    keep = np.isin(cr.cover, np.array(spec["marks"]))
    x, y, v = cr.easting[keep], cr.northing[keep], cr.resid_mm[keep]
    cov = cr.cover[keep]
    if x.size < 3:
        raise ValueError(f"{epoch}/{treatment}: only {x.size} marks; nothing to fit")

    if spec["drift"]:
        X, labels, ev = RF.cover_design(cov, spec["marks"])
        x0_drift = ev(spec["predict_at"])
    else:
        X, labels, x0_drift = None, ("const",), None

    rows = []
    for ml in max_lags_m:
        for est in estimators:
            model, centers, gamma, counts = RF.fit_field(
                x, y, v, max_lag_m=ml, n_lags=n_lags, n_pairs=n_pairs,
                estimator=est, seed=seed)
            kr = RF.krige(x, y, v, model, easting, northing,
                          X=X, x0_drift=x0_drift, drift_labels=tuple(labels))
            finite = centers[np.isfinite(gamma) & (counts > 0)]
            pinned = bool(finite.size and abs(model.range_ - finite.max()) < 1.0)
            rows.append(VariogramRow(
                max_lag_m=float(ml), n_lags=int(n_lags), estimator=str(est),
                nugget=float(model.nugget), sill=float(model.sill),
                range_m=float(model.range_), range_pinned_to_largest_lag=pinned,
                value_mm=float(kr.value_mm), sd_field_mm=float(kr.sd_field_mm),
                sd_new_mark_mm=float(kr.sd_new_mark_mm)))

    vals = np.array([r.value_mm for r in rows])
    sdf = np.array([r.sd_field_mm for r in rows])
    notes = []
    # The kriging variance depends ENTIRELY on the nugget/sill partition, and that
    # partition is set by the variogram's behaviour at short lags -- where a control
    # set spread over ~200 km has very few pairs.  Report the sweep range for it for
    # the same reason the value's range is reported: a single number would assert an
    # identifiability the data does not have.
    notes.append(
        f"sd_field over the sweep: {sdf.min():.2f} .. {sdf.max():.2f} mm "
        f"(median {np.median(sdf):.2f}); it is set by the fitted nugget/sill split, "
        f"which ranges nugget {min(r.nugget for r in rows):.0f}..{max(r.nugget for r in rows):.0f} "
        f"and sill {min(r.sill for r in rows):.0f}..{max(r.sill for r in rows):.0f} "
        f"on these same marks")
    n_pin = sum(r.range_pinned_to_largest_lag for r in rows)
    if n_pin:
        notes.append(
            f"the fitted range pinned to the largest lag centre in {n_pin} of "
            f"{len(rows)} sweep rows: the control does not determine a range, so the "
            f"sweep spread -- not any single fit -- is the honest statement")
    return DatumEstimate(
        epoch=epoch, route="published_residual_field",
        site_easting=float(easting), site_northing=float(northing), crs=load.crs,
        surface=surface, cover_treatment=treatment,
        predict_at_cover=spec["predict_at"], n_marks=int(x.size),
        constant_mm=float(np.median(vals)),
        constant_min_mm=float(vals.min()), constant_max_mm=float(vals.max()),
        sd_field_mm=float(np.median(sdf)),
        sd_field_min_mm=float(sdf.min()), sd_field_max_mm=float(sdf.max()),
        sd_new_mark_mm=float(np.median([r.sd_new_mark_mm for r in rows])),
        sd_is_of=("the error in predicting the spatially correlated component of the "
                  "DELIVERED surface's offset at this coordinate, nugget filtered out "
                  "as uncorrelated noise; NOT an SE of a mean over marks, and NOT what "
                  "a single new mark here would read (that is sd_new_mark_mm)"),
        sweep=[asdict(r) for r in rows],
        bridge_mm=bridge_mm, bridge_sd_mm=bridge_sd_mm, notes=notes)


def sweep_treatments(epoch: str, *, surface=None, treatments=None, **kw):
    """Every cover treatment for this epoch, so the choice's cost is visible."""
    names = tuple(TREATMENTS[epoch]) if treatments is None else tuple(treatments)
    return {t: datum_at_site(epoch, treatment=t, surface=surface, **kw) for t in names}
