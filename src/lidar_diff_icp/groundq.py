"""The ground percentile as a function of the cell's own class-2 spread.

WHAT THIS IS FOR, AND WHAT IT IS NOT. The pipeline takes the per-cell MEDIAN of the
ground-class returns, ``ground_q = 0.50``, and on OPEN GROUND that is the right answer:
against the 227 NVA control marks it lands -3.5 mm from surveyed truth, and the curve fitted
here is measurably WORSE -- held out, RMS 52.5 mm against 49.1 for the plain median. So this
module does NOT supply a default. ``difference_dem`` uses 0.50 unless a curve is named, and
:func:`load_curve` refuses every shortcut that would let one be applied unexamined.

The relation itself is real. It is a statement about VEGETATED ground, because that is where
the marks carrying it were sited: pooled over all 519 control marks the curve looked like a
16% RMS improvement, but the control set is three populations -- NVA (open, n=227, class-2
median -3.5 mm from truth), VVA (sited UNDER vegetation by design, n=162, +103.3 mm) and LCP
(the acquisition's own calibration points, n=130, -23.1 mm). The pooled curve's entire falling
limb is the VVA marks. Applying it to ordinary ground widened the DoD's scatter at both sites
tested (Elba NMAD 74.8 -> 79.1 mm, Whitewater 85.0 -> 92.4 on common cells).

The full account, with the measurements that overturned the earlier one, is in
``analysis/GROUND_Q_FROM_CLASS2_SPREAD.md``.

WHAT A CURVE IS. An isotonic (monotone non-increasing) regression of

    rank of the surveyed ground within a mark's class-2 returns

on ``log(class-2 standard deviation in mm)``. Monotone because more contamination cannot mean
a HIGHER ground rank; isotonic rather than a fitted form because the shape should come from
the data rather than from a functional family or a threshold.

TWO REFUSALS, NOT CONVENTIONS. A curve is valid only for the epoch AND the control point
types it was fitted on, and :func:`load_curve` raises rather than assume either. 2008 was
flown leaf-off in November and 2021 at green-up in May, with different classifiers; and NVA,
VVA and LCP are different populations that the project's datum work has required be stated
since ``ground_control/run_bridge_gen2.py`` -- "--point-types NVA is the consequence, not a
preference". Pooling them is the specific mistake this module was built around and then made.
"""
from __future__ import annotations

import os

import numpy as np

#: Where calibrate_ground_q.py writes its curves.
CURVE_DIR = os.path.join("data", "derived")


def curve_path(epoch, point_types):
    """Canonical path for a curve, e.g. ('gen2_2021_control', ['NVA']).

    The point types are IN THE FILENAME because they decide what the curve means. NVA, VVA
    and LCP marks are three different populations -- measured on this control, the class-2
    median sits -3.5 mm from truth at NVA, +103.3 mm at VVA and -23.1 mm at LCP -- so a
    curve is only interpretable if you can see which went into it.
    """
    tag = "-".join(sorted(t.upper() for t in point_types))
    return os.path.join(CURVE_DIR, f"ground_q_vs_class2sd_{epoch}_{tag}.npz")


def load_curve(path_or_epoch, *, expect_epoch=None):
    """Load a calibration curve, refusing a mismatch rather than silently applying it.

    Accepts a path or an epoch name. ``expect_epoch`` is checked against the epoch recorded
    IN the file, not against the filename, so a renamed file cannot smuggle the wrong
    calibration through.
    """
    p = path_or_epoch if str(path_or_epoch).endswith(".npz") else curve_path(path_or_epoch)
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"no ground-q curve at {p}. It is produced by\n"
            f"    ./lidar-icp/bin/python analysis/calibrate_ground_q.py --set <epoch>\n"
            f"from the surveyed control marks. There is no default curve and no fallback: "
            f"applying an uncalibrated percentile would put an unmeasured bias into every "
            f"elevation on the tile.")
    z = np.load(p, allow_pickle=True)
    if "point_types" not in z.files:
        raise ValueError(
            f"{p} records no point_types. Curves written before 2026-09-04 pooled NVA, VVA "
            f"and LCP marks, which are three different populations: the class-2 median sits "
            f"-3.5 mm from truth at NVA, +103.3 mm at VVA (sited UNDER vegetation by design) "
            f"and -23.1 mm at LCP. A pooled curve's falling limb is the VVA marks, so it "
            f"cannot be read as a vegetation correction for ordinary ground. Refit stating "
            f"the types:\n"
            f"    analysis/calibrate_ground_q.py --set <epoch> --point-types NVA")
    got = str(z["set"]) if "set" in z else None
    if expect_epoch is not None and got != expect_epoch:
        raise ValueError(
            f"{p} was calibrated on {got!r}, but {expect_epoch!r} was asked for. A curve is "
            f"valid only for its own epoch -- 2008 is leaf-off November, 2021 is green-up "
            f"May, and the classifiers differ. Calibrate the epoch you mean.")
    return {"log_sd_mm": np.asarray(z["log_sd_mm"], float),
            "q": np.asarray(z["q"], float),
            "epoch": got,
            "point_types": str(z["point_types"]),
            "n_marks": int(z["n_marks"]) if "n_marks" in z else None,
            "path": p,
            "provenance": {k: str(z[k]) for k in
                           ("fitted_on", "response", "covariate", "shape", "cv",
                            "known_limits") if k in z}}


def q_from_spread(sd_mm, curve, *, min_count=None, count=None):
    """Ground percentile per cell, from that cell's class-2 spread in MILLIMETRES.

    Outside the calibrated range the curve is held at its end values rather than
    extrapolated: an isotonic fit has no meaningful behaviour beyond its data, and a linear
    continuation would run q out of [0, 1] -- the failure that dogged every cover-relation
    version of this correction.

    Returns NaN where the spread is unusable, and where ``count`` is below ``min_count`` if
    both are given. NaN means "this method declines to estimate here", which the caller must
    handle; it is deliberately not the 0.50 default, because falling back silently is how an
    uncorrected cell would end up looking corrected.
    """
    sd = np.asarray(sd_mm, float)
    q = np.full(sd.shape, np.nan)
    ok = np.isfinite(sd) & (sd > 0)
    if count is not None and min_count is not None:
        ok &= np.asarray(count, float) >= min_count
    if ok.any():
        q[ok] = np.interp(np.log(sd[ok]), curve["log_sd_mm"], curve["q"])
    return q


def describe(curve, q=None):
    """One block of text naming the curve and what it did -- for a run's own output."""
    out = [f"ground_q from the class-2 spread: {curve['path']}",
           f"  epoch {curve['epoch']}, calibrated on {curve['n_marks']} surveyed marks"]
    out.insert(2, f"  point types {curve['point_types']} -- NVA/VVA/LCP are different "
                 f"populations and are never pooled")
    for k in ("fitted_on", "shape", "cv", "known_limits"):
        if k in curve["provenance"]:
            out.append(f"  {k}: {curve['provenance'][k]}")
    if q is not None:
        qq = np.asarray(q, float)
        fin = np.isfinite(qq)
        if fin.any():
            out.append(f"  q applied: median {np.nanmedian(qq):.3f}  "
                       f"p10 {np.nanpercentile(qq, 10):.3f}  min {np.nanmin(qq):.3f}; "
                       f"{int(np.sum(qq < 0.45)):,} cells below 0.45, "
                       f"{int((~fin).sum()):,} cells declined")
    return "\n".join(out)


# ---------------------------------------------------------------------------------------
# APPLYING THE CURVE TO A TILE
#
# Lifted from analysis/ridgelines/dod_cover_corrected.py, which is the code that actually
# produced Elba's corrected DoD. It is moved here rather than rewritten, because the one
# time it was rewritten -- inside pipeline._stream_ground -- the residual frame silently
# changed (no plane on one call, the REGIONAL plane on the other, an anchor-relative window
# instead of a fixed one) and the correction made stable ground worse: Whitewater's
# stable_sigma went 86.3 -> 93.0 mm against 86.3 uncorrected.
#
# THE FRAME IS THE METHOD. At the control marks the spread and the rank are measured on the
# slope-normal residual to an order-2 surface fitted through the mark's own box
# (analysis/control_mode_shift.py:90). On a tile the analogue is the residual to the gen2
# grid plus its local gradient. A curve indexed by a spread measured in any other frame is
# indexed by a different quantity.
#
# THIS IS INTRINSICALLY TWO-PASS. The reference surface IS the q = 0.50 gen2 grid
# (`z_after.npy`, saved from pipeline's own Z21), so the correction cannot be folded into
# the pass that builds that grid: grid first, then re-read the cloud and correct against it.


def surface_from_grid(z, x0, y0, res):
    """A gridded ground surface and its per-cell local plane -- the frame the correction
    lives in.

    NaN cells are filled from their nearest finite neighbour before the gradient is taken,
    so a hole does not propagate a NaN plane into every cell around it. The fill affects the
    PLANE only; ``valid`` records which cells the grid actually had, and a cell with no
    returns still ends with a NaN ground.
    """
    from scipy.ndimage import distance_transform_edt

    ny, nx = z.shape
    valid = np.isfinite(z)
    zf = z.copy()
    if not valid.all():
        zf = zf[tuple(distance_transform_edt(~valid, return_distances=False,
                                             return_indices=True))]
    gy, gx = np.gradient(zf, res)
    return {"x0": x0, "y0": y0, "res": float(res), "ny": ny, "nx": nx,
            "z": zf.ravel(), "valid": valid.ravel(),
            "dzde": gx.ravel(), "dzdn": gy.ravel(),
            "nnorm": np.sqrt(gx.ravel() ** 2 + gy.ravel() ** 2 + 1.0)}


def reference_surface(tile_dir, *, grid_file="z_after.npy"):
    """:func:`surface_from_grid` for a tile on disk, reading its own grid geometry."""
    from . import registration as reg

    j = reg.read_corrections(tile_dir)
    b = j["bounds"]
    surf = surface_from_grid(np.load(os.path.join(tile_dir, grid_file)),
                             b[0], b[1], float(j["res_m"]))
    surf.update(tile=tile_dir, grid_file=grid_file)
    return surf


def column_histogram(cloud, surf, *, zlo=-1.0, zhi=2.0, dz=0.02, classes=(2,),
                     chunk=3_000_000):
    """Per-cell histogram of the slope-normal residual of ``cloud`` to ``surf``.

    ``classes=None`` takes every return; ``(2,)`` takes the ASPRS ground class. One
    streaming pass, peak RAM O(cells x bins). Returns ``(H, n_in)``.
    """
    import laspy

    x0, y0, res, ny, nx = surf["x0"], surf["y0"], surf["res"], surf["ny"], surf["nx"]
    zflat, gxf, gyf, nn = surf["z"], surf["dzde"], surf["dzdn"], surf["nnorm"]
    nz = int(round((zhi - zlo) / dz))
    H = np.zeros((ny * nx, nz), np.int32)
    n_in = 0
    with laspy.open(str(cloud)) as f:
        for pts in f.chunk_iterator(chunk):
            cl = np.asarray(pts.classification)
            keep = np.ones(cl.shape, bool) if classes is None else np.isin(cl, classes)
            x = np.asarray(pts.x)[keep]; y = np.asarray(pts.y)[keep]; z = np.asarray(pts.z)[keep]
            ix = ((x - x0) / res).astype(np.int64); iy = ((y - y0) / res).astype(np.int64)
            ing = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
            cc = iy[ing] * nx + ix[ing]
            xc = x0 + ((cc % nx) + 0.5) * res
            yc = y0 + ((cc // nx) + 0.5) * res
            h = (z[ing] - (zflat[cc] + gxf[cc] * (x[ing] - xc)
                           + gyf[cc] * (y[ing] - yc))) / nn[cc]
            zi = np.floor((h - zlo) / dz).astype(np.int64)
            m = (zi >= 0) & (zi < nz)
            np.add.at(H, (cc[m], zi[m]), 1)
            n_in += int(m.sum())
    return H, n_in


def spread_from_histogram(H, zlo, dz, *, min_count=20):
    """The per-cell class-2 spread in MM -- the statistic the curve is indexed by.

    The plain standard deviation about the column's own mean, which is what ``np.std(hg)``
    measures at the marks. NOT a robust spread: a robust one is a different statistic, and
    indexing the curve with it is what pushed Whitewater's ground down by ~1 m in the worst
    tenth of cells. Cells with fewer than ``min_count`` returns get NaN -- the same minimum
    the calibration required of a mark. Returns ``(sd_mm, count)``.
    """
    nz = H.shape[1]
    ntot = np.cumsum(H, 1).astype(float)[:, -1]
    ctr = zlo + (np.arange(nz) + 0.5) * dz
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_h = (H * ctr).sum(1) / np.maximum(ntot, 1)
        var_h = (H * (ctr - mean_h[:, None]) ** 2).sum(1) / np.maximum(ntot, 1)
    return np.where(ntot >= min_count, np.sqrt(np.maximum(var_h, 0)) * 1000.0, np.nan), ntot


def ground_at_q(H, zlo, dz, q):
    """Slope-normal ground height in MM at the per-cell quantile ``q``.

    Read off the per-cell CDF with in-bin linear interpolation. A NaN ``q`` propagates to a
    NaN height: a cell the curve declined to estimate is left out, never quietly given the
    default.
    """
    nc, nz = H.shape
    C = np.cumsum(H, 1).astype(float)
    ntot = C[:, -1]
    idx = np.arange(nc)
    qq = np.asarray(q, float)
    r = qq * ntot
    k = (C >= r[:, None]).argmax(1)
    below = np.where(k > 0, C[idx, np.maximum(k - 1, 0)], 0.0)
    inbin = C[idx, k] - below
    frac = np.where(inbin > 0, (r - below) / np.maximum(inbin, 1e-9), 0.0)
    g = np.where(ntot > 0, (zlo + (k + np.clip(frac, 0, 1)) * dz) * 1000.0, np.nan)
    # A NaN q must NOT fall through to bin 0. `C >= nan` is False everywhere, so argmax
    # returns 0 and the cell silently takes the FLOOR of the column: measured on Elba,
    # 1,131 declined cells carried a median DoD of -1.091 m, a metre of fake erosion,
    # where their plain gen2 median said +0.003 m. Decline explicitly instead.
    return np.where(np.isfinite(qq), g, np.nan)


def ground_at_median(H, zlo, dz):
    """The plain per-cell median height in MM -- the uncorrected estimate, for comparison.

    Bin CENTRE, not interpolated: this is the control the correction is judged against, and
    it is reproduced exactly as the Elba producer computed it.
    """
    C = np.cumsum(H, 1).astype(float)
    ntot = C[:, -1]
    k = (C >= 0.5 * ntot[:, None]).argmax(1)
    return np.where(ntot > 0, (zlo + (k + 0.5) * dz) * 1000.0, np.nan)


def correct_gen2(gen2_cloud, curve, *, surf=None, tile_dir=None, zlo=-1.0, zhi=2.0,
                 dz=0.02, min_count=20, chunk=3_000_000, verbose=True):
    """The whole correction: cloud in, corrected gen2 ground out.

    Give either ``surf`` (:func:`surface_from_grid`, for a caller that already holds the
    grid) or ``tile_dir`` (to read the tile's ``z_after.npy``). ``curve`` is a loaded curve
    (:func:`load_curve`) or a path/epoch. Returns a dict with ``h2_mm`` (corrected,
    slope-normal, MM above the plane), ``h2_median_mm`` (the uncorrected control), ``q``,
    ``sd_mm``, ``count``, the histogram ``H`` and the reference surface ``surf``.
    """
    if not isinstance(curve, dict):
        curve = load_curve(curve)
    if surf is None:
        if tile_dir is None:
            raise ValueError("correct_gen2 needs either surf= (an in-memory grid) or "
                             "tile_dir= (to read that tile's z_after.npy)")
        surf = reference_surface(tile_dir)
    H, n_in = column_histogram(gen2_cloud, surf, zlo=zlo, zhi=zhi, dz=dz, chunk=chunk)
    if verbose:
        print(f"gen2: {n_in:,} class2 returns in {zlo:+.2f}..{zhi:+.2f} m over "
              f"{int((H.sum(1) > 0).sum()):,} cells", flush=True)
    sd_mm, count = spread_from_histogram(H, zlo, dz, min_count=min_count)
    q = q_from_spread(sd_mm, curve, min_count=min_count, count=count)
    n_oob = int(np.sum((q < 0) | (q > 1)))
    if n_oob and verbose:
        print(f"q outside [0,1] on {n_oob:,} cells; clipped to the column's ends, which is a "
              f"FLOOR/CEILING not a fit -- those cells take the extreme return, not a "
              f"percentile")
    q = np.clip(q, 0.0, 1.0)
    return {"h2_mm": ground_at_q(H, zlo, dz, q), "h2_median_mm": ground_at_median(H, zlo, dz),
            "q": q, "sd_mm": sd_mm, "count": count, "H": H, "surf": surf, "n_in": n_in,
            "curve": curve, "column": {"classes": "class2", "window_lo_m": zlo,
                                       "window_hi_m": zhi, "bin_m": dz,
                                       "min_count": min_count}}

# ---------------------------------------------------------------------------------------
# CALIBRATING THE CURVE -- the other half of the same object
#
# The application above and the calibration below MUST measure the same statistic, or the
# curve is indexed by a quantity it was not fitted on. That is not a hypothetical: indexing
# it with a robust spread (1.4826*(p75-p25)/1.349) of a residual to a different plane pushed
# Whitewater's ground down by ~1 m in the worst tenth of cells. So the definition lives here,
# once, and tests/test_groundq.py pins that spread_from_histogram converges to it.
#
# The mark-reading itself stays in analysis/calibrate_ground_q.py: it depends on the control
# CSVs and the per-mark boxes, which are survey bookkeeping rather than method.


def mark_statistics(hg, mu_true):
    """The covariate and the response for ONE control mark.

    ``hg``: that mark's ground-class returns as slope-normal heights above its own fitted
    order-2 surface, in METRES. ``mu_true``: the surveyed ground in that same frame, metres.

    ``sd``, ``nmad`` and ``iqr`` come back in METRES and ``med_minus_truth`` in MM, which is
    how the calibration has always carried them. ``sd`` is the PLAIN standard deviation --
    the quantity :func:`spread_from_histogram` reproduces on a tile.
    """
    hg = np.asarray(hg, float)
    return {"sd": float(np.std(hg)),
            "rank": float(np.mean(hg < mu_true)),
            "nmad": float(1.4826 * np.median(np.abs(hg - np.median(hg)))),
            "iqr": float(np.percentile(hg, 75) - np.percentile(hg, 25)),
            "skew": float(np.mean(((hg - hg.mean()) / max(np.std(hg), 1e-9)) ** 3)),
            "med_minus_truth": (float(np.median(hg)) - mu_true) * 1000}


def spatial_folds(easting, northing, *, n_folds=5, block_m=10_000, seed=0):
    """Cross-validation folds of whole ``block_m`` blocks, so a fold never sees its own area.

    Marks near each other share a flight line, a phenology and a crew. Splitting them at
    random would let the curve be scored on a mark it effectively already saw, which is how
    a spatial calibration flatters itself. Returns ``(fold_index_per_mark, blocks)``.
    """
    e = (np.asarray(easting, float) // block_m).astype(int).astype(str)
    n = (np.asarray(northing, float) // block_m).astype(int).astype(str)
    blk = np.char.add(np.char.add(e, "_"), n)
    _, first = np.unique(blk, return_index=True)
    ub = blk[np.sort(first)]                 # order of first appearance, as pandas .unique()
    np.random.default_rng(seed).shuffle(ub)
    fold = {b: i % n_folds for i, b in enumerate(ub)}
    return np.array([fold[b] for b in blk], int), ub


def fit_curve(sd_mm, rank):
    """Isotonic, monotone NON-INCREASING, rank on log(spread).

    Monotone because more contamination cannot mean a HIGHER ground rank -- the only
    constraint imposed, and it is physical. Isotonic rather than a functional family so the
    flat-then-falling shape comes from the data instead of from a chosen form or a threshold.
    """
    from sklearn.isotonic import IsotonicRegression
    return IsotonicRegression(increasing=False, out_of_bounds="clip").fit(
        np.log(np.asarray(sd_mm, float)), np.asarray(rank, float))


def save_curve(path, iso, *, n_marks, epoch, fitted_on, response, covariate, shape, cv,
               known_limits, **extra):
    """Write a curve WITH its provenance, so a consumer can see what it may be applied to.

    ``extra`` carries any further provenance the caller must record -- notably
    ``point_types``, since the control types are different populations and which ones were
    fitted on decides what the curve means.
    """
    np.savez(path, log_sd_mm=iso.f_.x, q=iso.f_.y, n_marks=n_marks, set=epoch,
             fitted_on=fitted_on, response=response, covariate=covariate, shape=shape,
             cv=cv, known_limits=known_limits, **extra)
    return path
