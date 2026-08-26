"""Combine independent ties into one datum constant, with a budget that keeps
common-mode error at full size.

Why this is not ``np.mean``
---------------------------
Two ties at Elba -- 2210 on the west chain and 2036 on the east -- agree to 7.5 mm.
That number is smaller than either mark's own radius uncertainty (12.4 and 27.0 mm), so
it is a **consistency check, not an accuracy**. Averaging the two and quoting
``sigma/sqrt(2)`` would turn a check into a claim.

The reason is that the two ties do not have independent errors. They share:

* one lateral (Nuth & Kaeaeb) shift, measured at Elba and extrapolated 7-16 km;
* one alignment estimator, whose extent-dependent repeatability is measured at
  12.4 mm RMS on dz (``analysis/MISSION_TIME_DRIFT.md`` section 4);
* one ground source, one quantile, one surface order;
* one un-applied along-track drift term.

A term shared by every tie does not average down with *n*. So each budget term carries a
``kind``:

``"random"``
    Independent between marks -- averages down as the inverse-variance weighting says.
    A mark's radius spread and its own chain's link errors are random.
``"common"``
    Shared by every tie in the combination -- enters the total at **full size**, no
    matter how many marks are added. The lateral extrapolation and the alignment
    estimator's repeatability are common.
``"unmodelled"``
    A known gap with no measured distribution -- a bound, not a sigma. Reported on its
    own line and **never** folded into the quadrature total, because doing so would
    dress a knowledge gap up as a measurement.

:func:`combine_ties` therefore returns a :class:`DatumConstant` whose headline is one
number and whose uncertainty is a **table**. ``sigma_total_mm`` exists, but it is the
quadrature of the random and common terms only, and ``unmodelled_mm`` sits beside it
unabsorbed.

Sign convention, unchanged from :mod:`~lidar_diff_icp.groundtruth.tie`: the constant is
what you **ADD to gen1** to place it on the surveyed datum.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

KINDS = ("random", "common", "unmodelled")


@dataclass(frozen=True)
class BudgetTerm:
    """One line of the uncertainty budget.

    ``value_mm`` is a 1-sigma-equivalent magnitude in millimetres; ``kind`` decides how
    (or whether) it combines; ``source`` says where the number was measured, and is
    printed with it so a term cannot appear without its provenance.
    """

    name: str
    value_mm: float
    kind: str
    source: str
    applies_to: str = "all ties"

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(f"kind must be one of {KINDS}, not {self.kind!r}")
        if not self.source or len(self.source.split()) < 2:
            raise ValueError(
                f"budget term {self.name!r} needs a source saying where the number was "
                "measured; an unattributed error term is not a budget")
        if not np.isfinite(self.value_mm) or self.value_mm < 0:
            raise ValueError(f"budget term {self.name!r}: value_mm must be finite and >= 0")


@dataclass
class DatumConstant:
    """The combined datum constant and the budget behind it.

    ``value_mm`` is the constant to ADD to gen1. ``sigma_total_mm`` combines the random
    and common terms in quadrature; ``unmodelled_mm`` is listed separately and is NOT in
    it.
    """

    value_mm: float
    ties: list                       # [(point_id, tie_mm, sigma_mm), ...]
    terms: list                      # [BudgetTerm, ...]
    weighting: str
    convention: str = ("constant to ADD to gen1 (already in the reference-swath frame "
                       "and geoid-shifted to the checkpoints' geoid model) to place it "
                       "on surveyed NAVD88(GEOID18)")
    notes: list = field(default_factory=list)

    # ------------------------------------------------------------------ the three sums
    @property
    def random_mm(self) -> float:
        """Quadrature of the random terms AFTER the inverse-variance weighting.

        The weighted mean of *k* ties with weights ``w_i`` has variance
        ``sum(w_i^2 * s_i^2) / (sum w_i)^2``; with ``w_i = 1/s_i^2`` that is
        ``1/sum(1/s_i^2)``. Random terms that are not per-mark sigmas (a chain's link
        error, say) are divided by ``sqrt(k)`` -- the same averaging, stated explicitly.
        """
        k = max(len(self.ties), 1)
        per_mark = [s for _, _, s in self.ties if np.isfinite(s) and s > 0]
        v = 1.0 / np.sum([1.0 / s ** 2 for s in per_mark]) if per_mark else 0.0
        extra = [t.value_mm ** 2 / k for t in self.terms if t.kind == "random"]
        return float(np.sqrt(v + np.sum(extra)))

    @property
    def common_mm(self) -> float:
        """Quadrature of the common-mode terms. **Not** divided by anything.

        This is the property the regression test pins: adding more marks must never move
        this number.
        """
        v = [t.value_mm ** 2 for t in self.terms if t.kind == "common"]
        return float(np.sqrt(np.sum(v))) if v else 0.0

    @property
    def unmodelled_mm(self) -> float:
        """Largest single unmodelled bound. Reported, never added into the total."""
        v = [t.value_mm for t in self.terms if t.kind == "unmodelled"]
        return float(max(v)) if v else 0.0

    @property
    def sigma_total_mm(self) -> float:
        """Random and common in quadrature. Excludes ``unmodelled_mm`` by design."""
        return float(np.hypot(self.random_mm, self.common_mm))

    @property
    def spread_mm(self) -> float:
        """Max-min of the input ties -- the consistency check, not the accuracy."""
        v = [t for _, t, _ in self.ties]
        return float(max(v) - min(v)) if len(v) > 1 else float("nan")

    # ------------------------------------------------------------------------- reports
    @staticmethod
    def table_columns() -> dict:
        return {
            "term": "what the uncertainty term is",
            "kind": ("random = independent between marks (averages down); common = "
                     "shared by every tie (does NOT average down); unmodelled = a bound "
                     "with no measured distribution (reported, never added in)"),
            "mm": "1-sigma-equivalent magnitude, mm",
            "applies_to": "which ties the term acts on",
            "source": "where the number was measured",
        }

    def table_rows(self) -> list:
        rows = []
        for kind in KINDS:
            for t in self.terms:
                if t.kind == kind:
                    rows.append([t.name, t.kind, f"{t.value_mm:.1f}", t.applies_to, t.source])
        rows.append(["-- random subtotal (weighted)", "random", f"{self.random_mm:.1f}",
                     f"{len(self.ties)} ties", "quadrature of the random terms after weighting"])
        rows.append(["-- common subtotal", "common", f"{self.common_mm:.1f}", "all ties",
                     "quadrature; independent of the number of ties"])
        rows.append(["== sigma_total (random + common)", "", f"{self.sigma_total_mm:.1f}",
                     "the datum constant",
                     "quadrature of the two subtotals; EXCLUDES the unmodelled rows"])
        if self.unmodelled_mm:
            rows.append(["== largest unmodelled bound", "unmodelled",
                         f"{self.unmodelled_mm:.1f}", "the datum constant",
                         "NOT in sigma_total -- a gap, not a measurement"])
        return rows

    def to_dict(self) -> dict:
        """JSON-serialisable record, for the product sidecar."""
        return dict(
            datum_constant_mm=self.value_mm,
            sign_convention=self.convention,
            weighting=self.weighting,
            ties=[dict(point_id=p, tie_mm=t, sigma_mm=s) for p, t, s in self.ties],
            tie_spread_mm=self.spread_mm,
            uncertainty_budget=[dict(name=t.name, value_mm=t.value_mm, kind=t.kind,
                                     source=t.source, applies_to=t.applies_to)
                                for t in self.terms],
            sigma_random_mm=self.random_mm,
            sigma_common_mm=self.common_mm,
            sigma_total_mm=self.sigma_total_mm,
            unmodelled_bound_mm=self.unmodelled_mm,
            notes=list(self.notes),
        )

    def to_json(self, path) -> str:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        return str(path)


def combine_ties(ties, terms, *, weighting="inverse-variance on each mark's radius sigma",
                 notes=()) -> DatumConstant:
    """Combine per-mark ties into one datum constant.

    ``ties``   iterable of ``(point_id, tie_mm, sigma_mm)``. ``sigma_mm`` is the mark's
               own uncertainty (in this project, half its radius spread), used for the
               inverse-variance weight. A non-finite or non-positive sigma gets equal
               weight and says so in the notes rather than being dropped -- a control
               point is never silently discarded.
    ``terms``  iterable of :class:`BudgetTerm`. Terms are NOT invented here; the caller
               supplies them with their sources.

    Nothing is filtered. If a tie should not be an anchor, the caller excludes it
    explicitly and records why.
    """
    ties = [(str(p), float(t), float(s)) for p, t, s in ties]
    if not ties:
        raise ValueError("combine_ties needs at least one tie")
    terms = list(terms)
    for t in terms:
        if not isinstance(t, BudgetTerm):
            raise TypeError(f"terms must be BudgetTerm, got {type(t).__name__}")
    notes = list(notes)
    sig = np.array([s for _, _, s in ties], float)
    val = np.array([t for _, t, _ in ties], float)
    good = np.isfinite(sig) & (sig > 0)
    if not good.all():
        notes.append(
            "equal weight given to " +
            ", ".join(p for (p, _, _), g in zip(ties, good) if not g) +
            ": no positive sigma. Kept, not dropped.")
        sig = np.where(good, sig, np.nanmedian(sig[good]) if good.any() else 1.0)
    w = 1.0 / sig ** 2
    value = float(np.sum(w * val) / np.sum(w))
    return DatumConstant(value_mm=value, ties=ties, terms=terms, weighting=weighting,
                         notes=notes)
