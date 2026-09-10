"""The 2008 control residual as a SPATIAL FIELD, and its value at a named site.

Why this module exists
----------------------
Every estimate of gen1's vertical offset made so far reports ``sd / sqrt(n)`` of a
sample mean over whatever control marks happened to be readable. That answers *"how
well do we know the average over those marks"*. It does not answer *"what is the value
at Elba"*, and the two differ whenever the residual varies in space -- which it does:
the sample mean of the vendor's own published residual moves further than its own SE as
the search radius grows.

So model the field. The quantity modelled is the vendor's ``dnr_error_m``
(``Control Z - Surface Z``, positive = the delivered 2008 surface reads LOW), which
involves no lidar processing of ours at all. Fit a variogram to it with
:mod:`lidar_diff_icp.variogram`, krige to the site, and report a **prediction variance
at that location** rather than an SE of a mean.

What this module refuses to do
------------------------------
1. **It invents no cut.** There is no default radius, bin width, lag count, pair count,
   block size, minimum ``n`` or cover selection anywhere below. Every such quantity is a
   required argument of the function that uses it. The driver
   ``analysis/control_residual_field.py`` sweeps them and prints every point of the
   sweep.
2. **It does not choose a cover stratum.** :func:`stratify` returns the named strata
   side by side; picking one is the caller's act.
3. **It does not hide which uncertainty it is reporting.** :class:`KrigeResult` carries
   *two* standard deviations with different meanings, both named in the dataclass
   docstring, because the difference between them is the nugget and that is not a
   detail.

Sign convention, matching :mod:`lidar_diff_icp.groundtruth.tie`: positive means the
surface reads LOW, i.e. the value is the constant to ADD to the surface.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..variogram import VariogramModel, empirical_variogram, fit_spherical

_DATA = Path(__file__).with_name("data")

#: The bundled 2008 control set, as delivered by ``parse_mndnr_2008_control.py``.
DEFAULT_CONTROL_CSV = _DATA / "mn_dnr_2008_control_semn.csv"

#: The MnDNR land-cover classes, as they appear in the ``point_type`` column.
COVER_CLASSES = ("L1O", "L2T", "L3B", "L4F", "L5U", "other")

COVER_MEANING = {
    "L1O": "open terrain",
    "L2T": "tall weeds and crops",
    "L3B": "brush and low trees",
    "L4F": "forested",
    "L5U": "urban",
    "other": "unclassed (FSA targets)",
}


# --------------------------------------------------------------------- loading

@dataclass(frozen=True)
class ControlResiduals:
    """De-duplicated control marks and the vendor's own residual at each.

    Attributes
    ----------
    point_id, county, cover
        Per-mark identity. ``cover`` is the CSV's ``point_type`` column verbatim.
    easting, northing
        EPSG:26915 metres.
    resid_mm
        ``dnr_error_m * 1000`` = ``Control Z - Surface Z`` in millimetres. Positive =
        the delivered 2008 surface reads BELOW the surveyed mark.
    n_rows_in, n_dup_rows, n_dup_groups
        Rows read, rows that were part of a duplicate group, and the number of groups.
    """

    point_id: np.ndarray
    county: np.ndarray
    cover: np.ndarray
    easting: np.ndarray
    northing: np.ndarray
    resid_mm: np.ndarray
    n_rows_in: int
    n_dup_rows: int
    n_dup_groups: int

    def __len__(self) -> int:
        return int(self.resid_mm.size)


def load_residuals(csv_path=DEFAULT_CONTROL_CSV) -> ControlResiduals:
    """Read the bundled control CSV and de-duplicate on exact ``(E, N, Z)``.

    A mark on a county line is printed in both counties' validation reports; those rows
    are the same physical mark and must not enter a variogram twice. De-duplication is
    on the exact triple, keeping the first occurrence -- no tolerance, because a
    tolerance would be an invented parameter.
    """
    import csv as _csv

    rows = []
    with open(csv_path, newline="") as f:
        for r in _csv.DictReader(f):
            rows.append(r)
    n_rows_in = len(rows)

    seen: dict[tuple, int] = {}
    keep: list[dict] = []
    dup_rows = 0
    for r in rows:
        key = (r["easting"], r["northing"], r["elevation"])
        if key in seen:
            seen[key] += 1
            dup_rows += 1
            continue
        seen[key] = 1
        keep.append(r)
    n_dup_groups = sum(1 for v in seen.values() if v > 1)

    return ControlResiduals(
        point_id=np.array([r["point_id"] for r in keep]),
        county=np.array([r["county"] for r in keep]),
        cover=np.array([r["point_type"] for r in keep]),
        easting=np.array([float(r["easting"]) for r in keep]),
        northing=np.array([float(r["northing"]) for r in keep]),
        resid_mm=np.array([float(r["dnr_error_m"]) * 1000.0 for r in keep]),
        n_rows_in=n_rows_in,
        n_dup_rows=dup_rows,
        n_dup_groups=n_dup_groups,
    )


def check_sign_convention(csv_path=DEFAULT_CONTROL_CSV, *, tol_m: float):
    """Re-derive, on every row, which subtraction ``dnr_error_m`` actually is.

    ``tol_m`` has no default: it is the arithmetic tolerance the caller is willing to
    call "exact", and it belongs to the caller. Returns
    ``(n_rows, n_control_minus_surface, max_resid_cms, n_surface_minus_control,
    max_resid_smc)``.
    """
    import csv as _csv

    a, b = [], []
    with open(csv_path, newline="") as f:
        for r in _csv.DictReader(f):
            z = float(r["elevation"])
            s = float(r["dnr_surface_z_m"])
            e = float(r["dnr_error_m"])
            a.append((z - s) - e)
            b.append((s - z) - e)
    a = np.abs(np.array(a))
    b = np.abs(np.array(b))
    return len(a), int((a <= tol_m).sum()), float(a.max()), int((b <= tol_m).sum()), float(b.max())


def stratify(cr: ControlResiduals, covers) -> np.ndarray:
    """Boolean mask of marks whose cover class is in ``covers`` (an explicit sequence)."""
    covers = tuple(covers)
    unknown = set(covers) - set(COVER_CLASSES)
    if unknown:
        raise ValueError(f"unknown cover class(es) {sorted(unknown)}; known: {COVER_CLASSES}")
    return np.isin(cr.cover, covers)


# ----------------------------------------------------------------- the variogram

def fit_field(x, y, v, *, max_lag_m, n_lags, n_pairs, estimator, seed):
    """Empirical variogram + weighted spherical fit. Every argument is required.

    Returns ``(model, centers, gamma, counts)`` so the empirical points can be shown
    beside the fitted parameters -- a fitted range with no empirical variogram under it
    is unfalsifiable.
    """
    centers, gamma, counts = empirical_variogram(
        np.asarray(x, float), np.asarray(y, float), np.asarray(v, float),
        max_lag=max_lag_m, n_lags=n_lags, n_pairs=n_pairs,
        estimator=estimator, seed=seed,
    )
    model = fit_spherical(centers, gamma, counts)
    return model, centers, gamma, counts


def _gamma(h, model: VariogramModel):
    """Spherical semivariogram of the fitted model; gamma(0) = 0 exactly."""
    h = np.asarray(h, float)
    hr = np.clip(h / model.range_, 0.0, 1.0)
    g = model.nugget + model.sill * (1.5 * hr - 0.5 * hr ** 3)
    return np.where(h == 0.0, 0.0, g)


# ------------------------------------------------------------------- kriging

@dataclass(frozen=True)
class KrigeResult:
    """A kriged value at one location, with the two uncertainties named apart.

    Attributes
    ----------
    value_mm
        The kriged residual at the target, in mm. Sign as in
        :class:`ControlResiduals`.
    sd_new_mark_mm
        Standard deviation of the error in predicting **the residual a new control mark
        placed at this location would show**. Includes the nugget, i.e. it carries the
        micro-scale + mark-siting + survey variance that one mark realises.
    sd_field_mm
        Standard deviation of the error in predicting **the spatially correlated
        component of the field at this location**, with the nugget treated as
        uncorrelated noise and filtered out. This is the smaller of the two and is the
        uncertainty of the systematic offset of the delivered surface at the target, as
        distinct from what any single mark there would read. At a target that coincides
        exactly with a datum, kriging is an exact interpolator, ``sd_new_mark_mm`` goes to
        zero and this quantity is clamped there too; that case does not arise when
        predicting at an unsampled site or under leave-one-out.
    n_marks, sum_weights_pos
        Marks entering the system, and the sum of the positive kriging weights (a
        readable measure of how local the prediction is).
    drift_labels, drift_values
        The universal-kriging drift terms evaluated at the target, named.
    """

    value_mm: float
    sd_new_mark_mm: float
    sd_field_mm: float
    n_marks: int
    sum_weights_pos: float
    drift_labels: tuple
    drift_values: tuple


def _uk_system(x, y, model, X):
    """Left-hand side of the universal-kriging system in SEMIVARIOGRAM form.

    ``X`` is (n, p) of drift basis columns; ``X[:, 0]`` must be the constant 1 column
    (ordinary kriging is the p = 1 case). Returns the (n+p, n+p) matrix.
    """
    n = x.size
    p = X.shape[1]
    H = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
    A = np.zeros((n + p, n + p))
    A[:n, :n] = _gamma(H, model)
    A[:n, n:] = X
    A[n:, :n] = X.T
    return A


def krige(x, y, v, model: VariogramModel, x0, y0, *, X=None, x0_drift=None,
          drift_labels=("const",)):
    """Universal (or, with ``X = None``, ordinary) kriging of ``v`` at ``(x0, y0)``.

    No search neighbourhood, no maximum number of neighbours, no distance cut: the
    whole set enters every system. Those would all be invented parameters, and at
    n < 1000 they buy nothing.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    n = x.size
    if X is None:
        X = np.ones((n, 1))
        x0_drift = np.ones(1)
    X = np.asarray(X, float)
    x0_drift = np.asarray(x0_drift, float)
    p = X.shape[1]

    A = _uk_system(x, y, model, X)
    h0 = np.hypot(x - x0, y - y0)
    b = np.zeros(n + p)
    b[:n] = _gamma(h0, model)
    b[n:] = x0_drift
    sol = np.linalg.solve(A, b)
    lam = sol[:n]
    mu = sol[n:]

    value = float(lam @ v)
    # Semivariogram-form OK/UK variance: sigma^2 = lam.g0 + mu.f0
    var_new = float(lam @ b[:n] + mu @ x0_drift)
    # The nugget is uncorrelated between the target and every datum, so the kriging
    # weights for the signal are identical and the signal variance is var_new - nugget.
    var_field = var_new - model.nugget
    return KrigeResult(
        value_mm=value,
        sd_new_mark_mm=float(np.sqrt(max(var_new, 0.0))),
        sd_field_mm=float(np.sqrt(max(var_field, 0.0))),
        n_marks=n,
        sum_weights_pos=float(lam[lam > 0].sum()),
        drift_labels=tuple(drift_labels),
        drift_values=tuple(float(t) for t in x0_drift),
    )


def krige_many(x, y, v, model: VariogramModel, x0, y0, *, X=None, X0=None):
    """Krige at many targets from ONE factorisation of the training system.

    Same estimator as :func:`krige`, batched over the right-hand sides. Returns
    ``(value_mm, var_new_mm2)`` arrays. Used by the cross-validators, where refactorising
    per target would dominate the cost and change nothing.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    x0 = np.atleast_1d(np.asarray(x0, float))
    y0 = np.atleast_1d(np.asarray(y0, float))
    n = x.size
    if X is None:
        X = np.ones((n, 1))
        X0 = np.ones((x0.size, 1))
    X = np.asarray(X, float)
    X0 = np.asarray(X0, float).reshape(x0.size, -1)
    p = X.shape[1]

    A = _uk_system(x, y, model, X)
    H0 = np.hypot(x[None, :] - x0[:, None], y[None, :] - y0[:, None])   # (m, n)
    B = np.zeros((n + p, x0.size))
    B[:n, :] = _gamma(H0, model).T
    B[n:, :] = X0.T
    sol = np.linalg.solve(A, B)
    lam = sol[:n, :]
    mu = sol[n:, :]
    value = lam.T @ v
    var_new = np.einsum("ij,ij->j", lam, B[:n, :]) + np.einsum("ij,ij->j", mu, X0.T)
    return value, var_new


def loo_errors(x, y, v, model: VariogramModel, *, X=None):
    """Exact leave-one-out kriging errors and variances, from one matrix inverse.

    Block-inverse (Dubrule) identity for the augmented system: with
    ``B = A^-1``, the LOO error at ``i`` is ``-(B z)_i / B_ii`` and the LOO kriging
    variance is ``-1 / B_ii``, where ``z`` is the data vector padded with zeros in the
    drift rows. Verified against brute-force refitting by
    :func:`verify_loo_shortcut`.

    Returns ``(err_mm, var_mm2)``; ``err = prediction - observation``.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    n = x.size
    if X is None:
        X = np.ones((n, 1))
    X = np.asarray(X, float)
    p = X.shape[1]
    A = _uk_system(x, y, model, X)
    B = np.linalg.inv(A)
    dB = np.diag(B)[:n]
    z = np.zeros(n + p)
    z[:n] = v
    err = -(B @ z)[:n] / dB
    var = -1.0 / dB
    return err, var


def verify_loo_shortcut(x, y, v, model: VariogramModel, idx, *, X=None):
    """Refit kriging from scratch with each ``idx`` held out; return max |difference|.

    A shortcut that is not checked is a claim. ``idx`` is caller-supplied: there is no
    default sample of indices to verify on.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    n = x.size
    if X is None:
        X = np.ones((n, 1))
    X = np.asarray(X, float)
    fast_err, fast_var = loo_errors(x, y, v, model, X=X)
    de, dv = 0.0, 0.0
    for i in idx:
        k = np.ones(n, bool)
        k[i] = False
        r = krige(x[k], y[k], v[k], model, x[i], y[i], X=X[k], x0_drift=X[i])
        de = max(de, abs((r.value_mm - v[i]) - fast_err[i]))
        dv = max(dv, abs(r.sd_new_mark_mm ** 2 - fast_var[i]))
    return de, dv


def block_cv(x, y, v, *, block_m, max_lag_m, n_lags, n_pairs, estimator, seed, X=None,
             refit_variogram, variogram_on):
    """Spatially blocked cross-validation on a square grid of side ``block_m``.

    Each block of the grid is one held-out fold; the model is fitted on the rest. With
    ``refit_variogram=True`` the variogram is re-estimated inside every training fold,
    which is the honest version; with ``False`` a single model fitted on all the data is
    reused, which is optimistic and is reported as such.

    ``variogram_on`` is ``"raw"`` (fit the variogram to ``v`` itself) or
    ``"ols_residual"`` (fit it to ``v`` minus the least-squares fit of the drift basis
    ``X``, which is what a drift term requires and what the universal-kriging system then
    assumes). It has no default because the right choice depends on whether ``X`` carries
    anything but the constant.

    Returns ``(err_mm, block_id, n_blocks)`` with one error per mark.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    n = x.size
    if X is None:
        X = np.ones((n, 1))
    X = np.asarray(X, float)

    bi = np.floor((x - x.min()) / block_m).astype(int)
    bj = np.floor((y - y.min()) / block_m).astype(int)
    bid = bi * (bj.max() + 1) + bj
    blocks = np.unique(bid)

    if variogram_on not in ("raw", "ols_residual"):
        raise ValueError("variogram_on must be 'raw' or 'ols_residual'")

    def _vgm_target(sel):
        if variogram_on == "raw":
            return v[sel]
        beta, *_ = np.linalg.lstsq(X[sel], v[sel], rcond=None)
        return v[sel] - X[sel] @ beta

    full_model = None
    if not refit_variogram:
        allsel = np.ones(n, bool)
        full_model, *_ = fit_field(x, y, _vgm_target(allsel), max_lag_m=max_lag_m,
                                   n_lags=n_lags, n_pairs=n_pairs, estimator=estimator,
                                   seed=seed)

    err = np.full(n, np.nan)
    for b in blocks:
        test = bid == b
        train = ~test
        if train.sum() <= X.shape[1]:
            continue
        if refit_variogram:
            model, *_ = fit_field(x[train], y[train], _vgm_target(train),
                                  max_lag_m=max_lag_m, n_lags=n_lags, n_pairs=n_pairs,
                                  estimator=estimator, seed=seed)
        else:
            model = full_model
        pred, _ = krige_many(x[train], y[train], v[train], model, x[test], y[test],
                             X=X[train], X0=X[test])
        err[test] = pred - v[test]
    return err, bid, int(blocks.size)


def constant_null_errors(v, fold_id):
    """Errors of the null model "one global constant", fold by fold.

    For each fold the constant is the mean (and, second column, the median) of the
    training data -- the same folds the kriging saw, so the comparison is like for like.
    Returns ``(err_mean_mm, err_median_mm)``.
    """
    v = np.asarray(v, float)
    fold_id = np.asarray(fold_id)
    em = np.full(v.size, np.nan)
    ed = np.full(v.size, np.nan)
    for f in np.unique(fold_id):
        test = fold_id == f
        train = ~test
        if train.sum() == 0:
            continue
        em[test] = np.mean(v[train]) - v[test]
        ed[test] = np.median(v[train]) - v[test]
    return em, ed


def cover_design(cover, classes):
    """Full-rank drift basis: constant + one indicator per class after the first.

    ``classes`` is caller-supplied and ordered; the FIRST entry is the reference class
    absorbed into the constant. Returns ``(X, labels, evaluator)`` where ``evaluator(c)``
    gives the drift row at which to predict for cover class ``c``.
    """
    classes = tuple(classes)
    cover = np.asarray(cover)
    cols = [np.ones(cover.size)]
    labels = [f"const(={classes[0]})"]
    for c in classes[1:]:
        cols.append((cover == c).astype(float))
        labels.append(f"is_{c}")
    X = np.column_stack(cols)

    def evaluator(c):
        row = np.zeros(len(classes))
        row[0] = 1.0
        if c != classes[0]:
            row[classes.index(c)] = 1.0
        return row

    return X, tuple(labels), evaluator


# ---------------------------------------------------------------------------
# EPOCH-AWARE CONTROL LOADING
#
# Folded in from ground_control/control.py on 2026-09-10, per INTEGRATION.md. That module
# existed only because its estimators had to stay untouched: "Making residual_field itself
# epoch-agnostic would mean editing src/, which this session may not do." That constraint
# is gone, so the adapter moves in and the module is deleted. No estimator was ever
# reimplemented there and none is reimplemented here.
#
# NAMED load_control_residuals, not load_control: gen1_datum.load_control already exists in
# this subpackage and returns a ControlSet, a different object. Two functions with one name
# and different return types is the trap this project has paid for twice (two datum.py, two
# `estimate`). The name says what it returns.
#
# The machinery below (fit_field, krige, loo_errors, block_cv, cover_design) takes plain
# (x, y, v) and was always epoch-agnostic; only these EDGES were gen1-schema-bound.
# ---------------------------------------------------------------------------

#: The control table each epoch is read from.
#:
#: ⚠ GEN1_CSV IS ONE SURVEY. It is the SE-Minnesota 2008 table, and gen1 is FOUR
#: acquisitions: `mn_dnr_control_swmn2010_metro2011.csv` (mnrv, battlecreek) and
#: `mn_dnr_control_arrowhead2011_duluth2012.csv` (cook, carlton) sit beside it in this
#: directory and are NOT reachable through here. So `load_control_residuals("gen1")`
#: means "the SE-MN 2008 control", not "gen1's control" -- correct for elba and
#: whitewater, silently empty of the other four sites' marks.
#:
#: This is the same shape as the geoid default that cost +54.87 mm at Ramsey: a name that
#: says "gen1" while meaning one survey. Making `project` explicit here is a BEHAVIOUR
#: change with ~8 call sites and is deliberately not bundled into the move that brought
#: this code here. Recorded rather than quietly carried forward.
GEN1_CSV = _DATA / "mn_dnr_2008_control_semn.csv"
GEN2_CSV = _DATA / "mn_se_driftless_2021_control.csv"

EPOCHS = ("gen1", "gen2")

#: gen2's four delivered surfaces.  ``<quality-level>_<product>``.
GEN2_SURFACES = ("ql1_dem", "ql1_laz", "ql0_dem", "ql0_laz")

#: The horizontal CRS each table's coordinates are published in.  These are carried and
#: asserted, never silently mixed: PROJ returns a null transform between them here, but
#: the NAD83(1986) -> NAD83(2011) realization difference is not modelled by that null
#: transform and is therefore NOT accounted for by treating them as equal.
EPOCH_CRS = {"gen1": "EPSG:26915", "gen2": "EPSG:6344"}

#: gen1's cover taxonomy is MnDNR's five classes; gen2's is USGS's binary accuracy class.
#: They are DIFFERENT TAXONOMIES.  Nothing in this module maps one onto the other; any
#: such mapping is a scientific choice and belongs to the caller, stated as one.
GEN1_COVER_CLASSES = ("L1O", "L2T", "L3B", "L4F", "L5U", "other")
GEN2_COVER_CLASSES = ("NVA", "VVA", "LCP")


def _f(v):
    v = (v or "").strip()
    return float(v) if v else None


def _rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# ------------------------------------------------------------------ sign convention

@dataclass(frozen=True)
class SignCheck:
    """Which subtraction the table's residual column actually is, re-derived per row."""

    epoch: str
    column: str
    n_rows_checked: int
    n_exact_surveyed_minus_lidar: int
    n_exact_lidar_minus_surveyed: int
    worst_miss_surveyed_minus_lidar_m: float
    worst_miss_lidar_minus_surveyed_m: float

    @property
    def is_surveyed_minus_lidar(self) -> bool:
        return (self.n_rows_checked > 0
                and self.n_exact_surveyed_minus_lidar == self.n_rows_checked)


def _sign_check(epoch, rows, z_col, err_col, column_label, tol_m):
    okf = okr = n = 0
    wf = wr = 0.0
    for r in rows:
        z, s, e = _f(r["elevation"]), _f(r[z_col]), _f(r[err_col])
        if None in (z, s, e):
            continue
        n += 1
        df, dr = abs((z - s) - e), abs((s - z) - e)
        okf += df < tol_m
        okr += dr < tol_m
        wf, wr = max(wf, df), max(wr, dr)
    return SignCheck(epoch, column_label, n, okf, okr, wf, wr)


def verify_sign_convention(epoch: str, *, tol_m: float, surface: str | None = None):
    """Re-derive the sign convention on every row.  ``tol_m`` is required, not assumed."""
    if epoch == "gen1":
        return [_sign_check("gen1", _rows(GEN1_CSV), "dnr_surface_z_m", "dnr_error_m",
                            "dnr_error_m", tol_m)]
    if epoch != "gen2":
        raise ValueError(f"epoch={epoch!r} must be one of {EPOCHS}")
    rows = _rows(GEN2_CSV)
    surfaces = GEN2_SURFACES if surface is None else (surface,)
    return [_sign_check("gen2", rows, f"usgs_{s}_z_m", f"usgs_{s}_error_m",
                        f"usgs_{s}_error_m", tol_m) for s in surfaces]


# ------------------------------------------------------------------------- loading

@dataclass(frozen=True)
class ControlLoad:
    """A loaded control set, plus everything about HOW it was selected.

    ``residuals`` is a :class:`residual_field.ControlResiduals`, so every estimator in
    that module accepts it unchanged.  The rest of the fields exist so a run cannot
    report a number without also being able to report which rows produced it.
    """

    residuals: object          # residual_field.ControlResiduals
    epoch: str
    crs: str
    surface: str | None
    n_rows_in: int
    n_with_residual: int
    n_marks_out: int
    n_dup_rows: int
    n_dup_groups: int
    roles_present: dict
    covers_present: dict
    dropped_no_residual: int


def load_control_residuals(epoch: str, *, surface: str | None = None,
                 roles: tuple | None = None) -> ControlLoad:
    """Load one epoch's control into ``residual_field.ControlResiduals``.

    Parameters
    ----------
    epoch : {"gen1", "gen2"}
    surface : str, required for gen2
        Which delivered surface the residual describes; one of :data:`GEN2_SURFACES`.
        There is no default -- gen2 publishes four and they are four different answers.
    roles : tuple, optional
        Restrict to these ``role`` values (gen2 only; gen1's table has no role column).
        Left as ``None`` this does NOT silently filter: it keeps every row that carries a
        residual, and reports in ``roles_present`` which roles survived.  gen2's 143 LCPs
        carry no residual in any of the four columns, so they drop out as a FACT of the
        table rather than as a choice made here -- see REPORT.md.

    De-duplication is on the exact ``(easting, northing, elevation)`` triple, no
    tolerance, keeping the first occurrence -- a tolerance would be an invented parameter.
    A mark on a county line is printed in both counties' validation reports and must not
    enter a variogram twice.
    """

    if epoch not in EPOCHS:
        raise ValueError(f"epoch={epoch!r} must be one of {EPOCHS}")
    if epoch == "gen2" and surface is None:
        raise ValueError(
            "epoch='gen2' requires surface=; it publishes four residual columns "
            f"({', '.join(GEN2_SURFACES)}) and they are four different answers. "
            "Picking one here would hide a choice that changes the number.")
    if epoch == "gen1" and surface is not None:
        raise ValueError("epoch='gen1' publishes ONE surface; surface= is meaningless")
    if epoch == "gen2" and surface not in GEN2_SURFACES:
        raise ValueError(f"surface={surface!r} must be one of {GEN2_SURFACES}")

    path = GEN1_CSV if epoch == "gen1" else GEN2_CSV
    rows = _rows(path)
    n_rows_in = len(rows)
    err_col = "dnr_error_m" if epoch == "gen1" else f"usgs_{surface}_error_m"

    keep_rows = [r for r in rows if _f(r[err_col]) is not None]
    n_with_residual = len(keep_rows)
    if roles is not None:
        if epoch == "gen1":
            raise ValueError("gen1's table has no 'role' column; roles= is meaningless")
        keep_rows = [r for r in keep_rows if r["role"] in roles]

    seen, keep, dup = {}, [], 0
    for r in keep_rows:
        k = (r["easting"], r["northing"], r["elevation"])
        if k in seen:
            seen[k] += 1
            dup += 1
            continue
        seen[k] = 1
        keep.append(r)
    n_dup_groups = sum(1 for v in seen.values() if v > 1)

    county = ([r.get("county", "") for r in keep] if epoch == "gen1"
              else [r.get("va_blocks", "") for r in keep])

    cr = ControlResiduals(
        point_id=np.array([r["point_id"] for r in keep]),
        county=np.array(county),
        cover=np.array([r["point_type"] for r in keep]),
        easting=np.array([float(r["easting"]) for r in keep]),
        northing=np.array([float(r["northing"]) for r in keep]),
        resid_mm=np.array([float(r[err_col]) * 1000.0 for r in keep]),
        n_rows_in=n_rows_in,
        n_dup_rows=dup,
        n_dup_groups=n_dup_groups,
    )

    def _count(key):
        out = {}
        for r in keep:
            out[r.get(key, "")] = out.get(r.get(key, ""), 0) + 1
        return dict(sorted(out.items()))

    return ControlLoad(
        residuals=cr, epoch=epoch, crs=EPOCH_CRS[epoch], surface=surface,
        n_rows_in=n_rows_in, n_with_residual=n_with_residual, n_marks_out=len(keep),
        n_dup_rows=dup, n_dup_groups=n_dup_groups,
        roles_present=_count("role") if epoch == "gen2" else {},
        covers_present=_count("point_type"),
        dropped_no_residual=n_rows_in - n_with_residual,
    )


def cover_from_point_id_prefix(point_id):
    """The BROKEN way to read gen1's cover, kept so its breakage can be demonstrated.

    21 gen1 marks carry a ``point_id`` beginning ``L10`` with a DIGIT ZERO rather than
    ``L1O`` with the letter O.  Parsing cover from the prefix silently drops them, giving
    209 open marks instead of the CSV's own 230.  Never use this to select marks; it
    exists so ``tests/`` can assert the difference and keep the trap documented.
    """
    return str(point_id)[:3]
