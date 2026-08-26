"""Flight-line vertical ties measured WHERE THEY ARE USED, and chained line to line.

The problem this exists for
---------------------------
``coreg.align_swaths`` solves one constant per flight line from the overlaps inside
whatever tile it is handed. At Elba that tile is ~4.5 x 4 km, and the constants land in
``data/derived/elbaext/corrections_geoid.json``. Those constants have then been applied to
surveyed marks up to **62.9 km** away along the same lines. That is an extrapolation of a
locally-fitted quantity, and this project has already measured the analogous quantity
going bad with distance: ``analysis/CROSS_LINE_FIT.md`` §6 shows the across-track
coefficient ``c_s`` running +23 to +429 mm per unit tangent across 3.5 km of ONE pair --
five to ten times its own fitted standard error -- and concludes that ``c_s`` "is a
property of *the piece of overlap it was fitted on*, not of the flight line".

Nothing had tested whether the *constant* ``k`` behaves the same way, because there was no
machinery to measure a tie **at a stated place**. This module is that machinery. It does
not decide anything about ``k``; it makes ``k`` measurable as a function of position so
the question can be answered with data.

What it does NOT do
-------------------
It does not reimplement the tie. Every vertical offset in this module comes out of
:func:`lidar_diff_icp.coreg.coregister_swaths` (via :func:`lidar_diff_icp.coreg.align_swaths`
for the network case), the same estimator the pipeline uses, in one of its two tie modes.
The only thing added is **where the points come from**: a window about a stated location
instead of a whole tile.

Two ties, and the module records which
--------------------------------------
``tie="overlap_median"`` is the shipped Nuth & Kaeaeb vertical reduction (median of the
overlap difference). ``tie="intercept"`` is :func:`lidar_diff_icp.coreg.across_track_tie`,
the LAD intercept at across-track position zero, which does not depend on which part of
the sidelap the extent happens to cover; it is what
``data/derived/elbaext/corrections_geoid.json`` records as ``swath_tie`` and is the
pipeline default. Both are supported and every result carries ``tie_mode``.

**The intercept tie needs real scan angles.** ``groundtruth.chain`` builds its pair clouds
with ``scan_angle`` set to zeros, which silently degrades ``tie="intercept"`` to the
median tie (``across_track_tie`` returns NaN on a zero-range predictor and
``coregister_swaths`` then falls back). This module therefore carries its own tile reader,
which keeps the scan angle -- read through
:func:`lidar_diff_icp.groundtruth.tie.scan_angle_deg`, which RAISES rather than returning
zeros when the dimension is absent. That is the one deliberate fork from ``chain.py``, and
it is why: the same cache cannot serve both, because ``chain.py``'s cache does not contain
the field.

No defaults
-----------
Window size, grid resolution, tie mode, window shape and the excluded classes are all
**required keyword arguments**. There is no default window: the size of the window is the
quantity under investigation, so choosing one inside the library would hide the answer.
:func:`window_ladder` exists to make the dependence visible, the way
``groundtruth.tie.radius_ladder`` does for a mark tie.

Sign conventions
----------------
``coregister_swaths(pc, swath_ref=a, swath_src=b)`` returns the shift that moves ``b``
onto ``a``, so ``LocalTie.dz_m`` is what to **add to line ``line_src``'s z** to bring it
into ``line_ref``'s frame at that window. Walking a chain outward from the reference,
``offset[b] = offset[a] + dz``; :attr:`LocalChain.dz_total_m` is what to add to the
source line's z to reach the target line's frame **at the mark**.

The question this module exists to answer, in code
--------------------------------------------------
"For a mark at ``(E, N)`` on line ``L``, what constant puts it in line 137's frame **at
this location**, and how well is that constant known?"::

    cache = TileCache(cache_dir="data/derived/localtie_cache")
    ch = chain_local(tiles, easting=E, northing=N, source_line=L, target_line=137,
                     half_width_m=400.0, shape="square", res_m=2.0, tie="intercept",
                     exclude=(5, 6, 9), cache=cache,
                     ladder_half_widths_m=[100.0, 200.0, 400.0, 800.0, 1200.0])
    ch.dz_total_mm            # add this to the mark's gen1 elevation
    ch.dz_sigma_window_m      # the error to quote: per-link window-ladder spreads
    ch.dz_sigma_formal_m      # coreg's own sigma -- measured to be ~30x too small
    ch.max_solve_distance_m   # how far from the mark the farthest link had to be solved
    compare_to_constants(ch, imported)   # against a constant fitted somewhere else

``ladder_half_widths_m`` is optional and has no default; without it
``dz_sigma_window_m`` is NaN and only the formal sigma is available, which
``analysis/LOCAL_TIE_CHAINING.md`` measures to be optimistic by one to two orders of
magnitude.

A degenerate window is reported, not hidden
-------------------------------------------
``coreg.nuth_kaab`` abandons its fit when fewer than 100 grid cells exceed its 3-degree
slope floor, and in that branch it returns ``dz = 0.0`` exactly, with ``n = 0`` -- not the
overlap median, and not NaN. On a whole tile that never happens; on a small window over
flat ground it happens easily, and a silent 0.0 mm tie is the worst possible failure for
this module's purpose. Every :class:`LocalTie` therefore carries ``degenerate`` (the fit
used zero cells) and ``dz_overlap_median_m``, an independent read of the same overlap by
:func:`lidar_diff_icp.swathdiff.swath_difference`. Nothing is dropped or corrected on that
basis -- the caller is told.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

from . import coreg, swathdiff
from .io import MN_GEN1_CRS, PointCloud
from .groundtruth.provenance import Param
from .groundtruth.tie import scan_angle_deg


# --------------------------------------------------------------------- tile reading

@dataclass
class TileCache:
    """Per-tile point arrays, kept in memory and optionally on disk.

    Unlike ``groundtruth.chain``'s cache this keeps **every class** and the **scan
    angle**, so the class filter and the tie mode stay the caller's choice at solve time
    rather than being baked into the cache. Coordinates are stored as float32 offsets
    from a per-tile integer origin (~1 mm resolution over a 3.5 km tile), which is what
    ``chain.py`` does and is well inside the millimetre level these ties are read at.

    ``release()`` drops the in-memory arrays; this is a shared laptop and a tile is
    ~7-10 million points.
    """

    cache_dir: str | None = None
    _arrays: dict = field(default_factory=dict, repr=False)
    _headers: dict = field(default_factory=dict, repr=False)

    def header_bbox(self, path: str):
        """(x0, y0, x1, y1) from the LAS header -- no decompression."""
        path = str(path)
        if path not in self._headers:
            import laspy
            with laspy.open(path) as f:
                h = f.header
                self._headers[path] = (float(h.mins[0]), float(h.mins[1]),
                                       float(h.maxs[0]), float(h.maxs[1]))
        return self._headers[path]

    def _npz(self, path: str):
        if not self.cache_dir:
            return None
        st = os.stat(path)
        return os.path.join(
            self.cache_dir,
            f"{os.path.basename(path)}.{st.st_size}.{int(st.st_mtime)}.localtie.npz")

    def tile(self, path: str) -> dict:
        """``dict`` of ``x, y, z, psid, cls, ang`` (scan angle in DEGREES) for one tile."""
        path = str(path)
        if path in self._arrays:
            return self._arrays[path]
        npz = self._npz(path)
        if npz and os.path.exists(npz):
            d = np.load(npz)
            o = d["origin"]
            arr = {"x": d["x"].astype(np.float64) + o[0],
                   "y": d["y"].astype(np.float64) + o[1],
                   "z": d["z"].astype(np.float64) + o[2],
                   "psid": d["psid"], "cls": d["cls"], "ang": d["ang"].astype(np.float64)}
        else:
            import laspy
            f = laspy.read(path)
            x = np.asarray(f.x); y = np.asarray(f.y); z = np.asarray(f.z)
            arr = {"x": x, "y": y, "z": z,
                   "psid": np.asarray(f.point_source_id).astype(np.int32),
                   "cls": np.asarray(f.classification).astype(np.uint8),
                   "ang": scan_angle_deg(f)}
            self._headers[path] = (float(x.min()), float(y.min()),
                                   float(x.max()), float(y.max()))
            if npz:
                os.makedirs(self.cache_dir, exist_ok=True)
                o = np.array([np.floor(x.min()), np.floor(y.min()), np.floor(z.min())])
                np.savez(npz, origin=o,
                         x=(x - o[0]).astype(np.float32), y=(y - o[1]).astype(np.float32),
                         z=(z - o[2]).astype(np.float32), psid=arr["psid"],
                         cls=arr["cls"], ang=arr["ang"].astype(np.float32))
            del f
        self._arrays[path] = arr
        return arr

    def release(self, path: str | None = None):
        if path is None:
            self._arrays.clear()
        else:
            self._arrays.pop(str(path), None)


def _window_mask(x, y, easting, northing, half_width_m, shape):
    """Points inside the window. ``shape`` is required and has no default."""
    m = (np.abs(x - easting) <= half_width_m) & (np.abs(y - northing) <= half_width_m)
    if shape == "square":
        return m
    if shape == "disk":
        m[m] = ((x[m] - easting) ** 2 + (y[m] - northing) ** 2) <= half_width_m ** 2
        return m
    raise ValueError(f"shape={shape!r} must be 'square' or 'disk'")


@dataclass
class Window:
    """gen1 returns within ``half_width_m`` of a location, with their scan angles."""

    pc: PointCloud
    easting: float
    northing: float
    half_width_m: float
    shape: str
    tiles: list
    line_counts: dict            # point_source_id -> points in the window (all classes)

    @property
    def lines(self) -> list:
        return sorted(self.line_counts)


def window_cloud(tile_paths, *, easting, northing, half_width_m, shape,
                 cache: TileCache | None = None, lines=None) -> Window:
    """Crop the tiles to a window about ``(easting, northing)``.

    Only tiles whose header bbox reaches the window are decompressed. ``lines``, when
    given, keeps only those ``point_source_id`` values. Classification is NOT filtered
    here -- that is the ``exclude`` argument of the tie, applied by ``coreg`` itself.
    """
    cache = TileCache() if cache is None else cache
    half_width_m = float(half_width_m)
    xs, ys, zs, ps, cl, an, used = [], [], [], [], [], [], []
    for path in [str(p) for p in tile_paths]:
        x0, y0, x1, y1 = cache.header_bbox(path)
        if (x1 < easting - half_width_m or x0 > easting + half_width_m
                or y1 < northing - half_width_m or y0 > northing + half_width_m):
            continue
        a = cache.tile(path)
        m = _window_mask(a["x"], a["y"], easting, northing, half_width_m, shape)
        if lines is not None:
            m &= np.isin(a["psid"], np.asarray(list(lines), a["psid"].dtype))
        if not m.any():
            continue
        xs.append(a["x"][m]); ys.append(a["y"][m]); zs.append(a["z"][m])
        ps.append(a["psid"][m]); cl.append(a["cls"][m]); an.append(a["ang"][m])
        used.append(path)
    if not xs:
        raise ValueError(
            f"no points within {half_width_m} m of ({easting}, {northing}) in "
            f"{len(list(tile_paths))} tile(s)")
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
    psid = np.concatenate(ps); cls = np.concatenate(cl); ang = np.concatenate(an)
    pc = PointCloud(x=x, y=y, z=z, point_source_id=psid, classification=cls,
                    gps_time=np.zeros(x.size), scan_angle=ang, crs=MN_GEN1_CRS)
    u, c = np.unique(psid, return_counts=True)
    return Window(pc=pc, easting=float(easting), northing=float(northing),
                  half_width_m=half_width_m, shape=shape, tiles=used,
                  line_counts={int(a): int(b) for a, b in zip(u, c)})


# ------------------------------------------------------------------- the pairwise tie

@dataclass
class LocalTie:
    """One flight-line pair's vertical tie, measured on one window.

    ``dz_m`` is what to ADD to ``line_src``'s z to bring it into ``line_ref``'s frame at
    this window. Everything needed to know how much to believe it is beside it.
    """

    line_ref: int
    line_src: int
    easting: float
    northing: float
    half_width_m: float
    shape: str
    tie_mode: str
    res_m: float
    exclude: tuple
    dz_m: float
    dz_sigma_m: float             # coreg's nmad_after / sqrt(n): formal, and optimistic
    dx_m: float
    dy_m: float
    n_nk_cells: int               # cells the Nuth & Kaeaeb fit used (0 = degenerate)
    dz_overlap_median_m: float    # independent read: median of the unshifted overlap
    n_overlap_cells: int
    nmad_overlap_m: float
    overlap_area_km2: float
    nmad_before_m: float
    nmad_after_m: float
    converged: bool
    degenerate: bool              # the N&K fit used ZERO cells; see the module docstring
    dtan_min: float               # across-track coordinate actually sampled in this window
    dtan_max: float
    dtan_median: float
    c_mm_per_tan: float           # the across-track slope fitted with the intercept
    k_check_m: float              # across_track_tie's intercept: EQUALS dz_m when tie='intercept'
    extrapolated: bool            # dtan = 0, where the intercept is read, is OUTSIDE the sample
    n_points_ref: int
    n_points_src: int
    tiles: list

    @property
    def dz_mm(self) -> float:
        return self.dz_m * 1000.0

    def row(self) -> list:
        return [f"{self.line_ref}-{self.line_src}", f"{self.half_width_m:.0f}",
                f"{self.dz_m * 1000:+.1f}", f"{self.dz_sigma_m * 1000:.1f}",
                f"{self.dz_overlap_median_m * 1000:+.1f}", self.n_nk_cells,
                self.n_overlap_cells, f"{self.overlap_area_km2:.4f}",
                f"{self.nmad_after_m * 1000:.0f}",
                f"{self.dtan_min:+.3f}", f"{self.dtan_max:+.3f}",
                f"{self.c_mm_per_tan:+.0f}",
                "YES" if self.extrapolated else "",
                "YES" if self.degenerate else ""]

    @staticmethod
    def columns() -> dict:
        return {
            "pair": "flight-line pair, ref-src; dz is added to src to reach ref's frame",
            "half_width_m": "window half-width about the stated location, m",
            "dz_mm": "the tie, from coreg.coregister_swaths in the stated tie mode, mm",
            "dz_sig_mm": "coreg's formal 1-sigma, nmad_after/sqrt(n_cells), mm -- optimistic",
            "dz_median_mm": "independent check: median of the UNSHIFTED overlap difference, mm",
            "nk_cells": "cells the Nuth & Kaeaeb fit used; 0 means the fit was abandoned",
            "ovl_cells": "cells where both lines have terrain returns in this window",
            "ovl_km2": "that overlap's area, km^2",
            "nmad1_mm": "robust scatter left in the overlap after the tie, mm",
            "dtan_lo": "smallest across-track coordinate tan(scan_ref)-tan(scan_src) in "
                       "the window, dimensionless; the intercept tie is read at dtan = 0",
            "dtan_hi": "largest across-track coordinate in the window, dimensionless "
                       "(a difference of tangents)",
            "c_mm_tan": "across-track slope fitted alongside the intercept, mm per unit "
                        "tangent (coreg.across_track_tie's second return)",
            "extrap": "YES/blank flag: YES when dtan = 0 lies OUTSIDE [dtan_lo, "
                      "dtan_hi], i.e. the intercept tie is an extrapolation here",
            "degenerate": "YES/blank flag: YES when the Nuth & Kaeaeb fit used zero "
                          "cells, in which case dz is coreg's abandoned-fit 0.0 m",
        }


def _across_track_diagnostics(pc, line_ref, line_src, res_m, exclude, dx, dy):
    """The across-track coordinate this window actually samples, and the slope on it.

    This MIRRORS the geometry of ``coreg.coregister_swaths``'s intercept branch -- same
    overlap bbox, same grids, same shift, same estimator -- so that ``k_check`` reproduces
    the tie that function returns in ``tie="intercept"`` mode. ``tests/test_localtie.py``
    asserts that equality, which is what keeps the mirroring honest rather than assumed.

    Its purpose is diagnostic: the intercept is read at ``dtan = 0``, and on a small
    window near the edge of a sidelap the sampled ``dtan`` may not reach zero at all, so
    the tie is an extrapolation. Nothing is corrected or dropped on that basis.
    """
    terr = ~np.isin(pc.classification, exclude)
    ma = terr & (pc.point_source_id == line_ref)
    mb = terr & (pc.point_source_id == line_src)
    x, y = pc.x, pc.y
    x0 = max(x[ma].min(), x[mb].min()); x1 = min(x[ma].max(), x[mb].max())
    y0 = max(y[ma].min(), y[mb].min()); y1 = min(y[ma].max(), y[mb].max())
    nx = int(np.ceil((x1 - x0) / res_m)); ny = int(np.ceil((y1 - y0) / res_m))
    sa = np.asarray(pc.scan_angle, float)
    t_ref = swathdiff._median_grid(x[ma], y[ma], np.tan(np.radians(sa[ma])),
                                   res_m, x0, y0, nx, ny)
    t_src = swathdiff._median_grid(x[mb], y[mb], np.tan(np.radians(sa[mb])),
                                   res_m, x0, y0, nx, ny)
    z_ref = swathdiff._median_grid(x[ma], y[ma], pc.z[ma], res_m, x0, y0, nx, ny)
    z_src = swathdiff._median_grid(x[mb], y[mb], pc.z[mb], res_m, x0, y0, nx, ny)
    dtan = t_ref - coreg._shift_grid(t_src, dx, dy, res_m)
    dh = z_ref - coreg._shift_grid(z_src, dx, dy, res_m)
    k, c, _n = coreg.across_track_tie(dh, dtan)
    m = np.isfinite(dtan) & np.isfinite(dh)
    if not m.any():
        return float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), False
    lo = float(np.nanmin(dtan[m])); hi = float(np.nanmax(dtan[m]))
    return lo, hi, float(np.nanmedian(dtan[m])), c * 1000.0, k, not (lo <= 0.0 <= hi)


def local_pair_tie(window: Window, line_ref: int, line_src: int, *,
                   res_m: float, tie: str, exclude) -> LocalTie:
    """The vertical tie of one pair, on one window, by ``coreg.coregister_swaths``.

    ``res_m``, ``tie`` and ``exclude`` are required: the repo's own values are 2.0,
    ``"intercept"`` (``corrections_geoid.json`` ``swath_tie``) and ``(5, 6, 9)``, but
    naming them here would hide them from the caller's provenance record.

    Raises ``ValueError`` -- from ``coreg`` -- when the two lines do not overlap inside
    the window. That is a property of the window, and it is a reportable result.
    """
    pc = window.pc
    exclude = tuple(exclude)
    c = coreg.coregister_swaths(pc, line_ref, line_src, res_m, exclude, tie=tie)
    sd = swathdiff.swath_difference(pc, line_ref, line_src, res_m, exclude)
    lo, hi, mid, c_mm, k_chk, extrap = _across_track_diagnostics(
        pc, line_ref, line_src, res_m, exclude, c.dx, c.dy)
    terr = ~np.isin(pc.classification, exclude)
    return LocalTie(
        line_ref=int(line_ref), line_src=int(line_src),
        easting=window.easting, northing=window.northing,
        half_width_m=window.half_width_m, shape=window.shape,
        tie_mode=tie, res_m=float(res_m), exclude=exclude,
        dz_m=float(c.dz), dz_sigma_m=float(c.dz_sigma),
        dx_m=float(c.dx), dy_m=float(c.dy), n_nk_cells=int(c.n),
        dz_overlap_median_m=float(sd.median_offset), n_overlap_cells=int(sd.n_cells),
        nmad_overlap_m=float(sd.robust_std),
        overlap_area_km2=sd.n_cells * res_m * res_m / 1e6,
        nmad_before_m=float(c.nmad_before), nmad_after_m=float(c.nmad_after),
        converged=bool(c.converged), degenerate=(int(c.n) == 0),
        dtan_min=lo, dtan_max=hi, dtan_median=mid, c_mm_per_tan=c_mm,
        k_check_m=k_chk, extrapolated=bool(extrap),
        n_points_ref=int((terr & (pc.point_source_id == line_ref)).sum()),
        n_points_src=int((terr & (pc.point_source_id == line_src)).sum()),
        tiles=list(window.tiles))


def pair_tie_at(tile_paths, line_ref: int, line_src: int, *, easting, northing,
                half_width_m, shape, res_m, tie, exclude,
                cache: TileCache | None = None) -> LocalTie:
    """Crop and tie in one call. See :func:`window_cloud` and :func:`local_pair_tie`."""
    w = window_cloud(tile_paths, easting=easting, northing=northing,
                     half_width_m=half_width_m, shape=shape, cache=cache,
                     lines=(line_ref, line_src))
    return local_pair_tie(w, line_ref, line_src, res_m=res_m, tie=tie, exclude=exclude)


# -------------------------------------------------------------------- window ladder

@dataclass
class WindowLadder:
    """How the tie moves as the window grows -- the analogue of the radius ladder.

    ``groundtruth.tie`` reports a mark tie over a ladder of radii and takes the *spread
    over the ladder* as the headline uncertainty, because on that data the radius term
    dominates the fit's own standard error. The same question is open here, and the same
    answer shape: a tie that moves by 40 mm between a 300 m and a 900 m window is not a
    40 mm-accurate number whatever its formal sigma says.
    """

    ties: list
    line_ref: int
    line_src: int
    easting: float
    northing: float
    half_widths_m: list

    @property
    def dz_mm(self) -> np.ndarray:
        return np.array([t.dz_m * 1000.0 for t in self.ties], float)

    @property
    def spread_mm(self) -> float:
        """max - min over the ladder, NaN-safe. NaN if fewer than two finite ties."""
        v = self.dz_mm[np.isfinite(self.dz_mm)]
        return float(v.max() - v.min()) if v.size > 1 else float("nan")

    @property
    def sd_mm(self) -> float:
        v = self.dz_mm[np.isfinite(self.dz_mm)]
        return float(v.std(ddof=1)) if v.size > 1 else float("nan")

    def rows(self):
        return [t.row() for t in self.ties]


def window_ladder(tile_paths, line_ref: int, line_src: int, *, easting, northing,
                  half_widths_m, shape, res_m, tie, exclude,
                  cache: TileCache | None = None) -> WindowLadder:
    """Tie the same pair at the same place over a caller-supplied list of window sizes.

    ``half_widths_m`` has no default. Windows that cannot be tied (the pair does not
    overlap inside them, or no points fall in them) are still reported -- as a
    :class:`LocalTie` with NaN ``dz_m`` and ``n_nk_cells = 0`` -- rather than skipped,
    because "too small to tie here" is the answer to the question being asked.
    """
    cache = TileCache() if cache is None else cache
    out = []
    for h in half_widths_m:
        try:
            out.append(pair_tie_at(tile_paths, line_ref, line_src, easting=easting,
                                   northing=northing, half_width_m=h, shape=shape,
                                   res_m=res_m, tie=tie, exclude=exclude, cache=cache))
        except ValueError as e:
            out.append(LocalTie(
                line_ref=int(line_ref), line_src=int(line_src),
                easting=float(easting), northing=float(northing), half_width_m=float(h),
                shape=shape, tie_mode=tie, res_m=float(res_m), exclude=tuple(exclude),
                dz_m=float("nan"), dz_sigma_m=float("nan"), dx_m=float("nan"),
                dy_m=float("nan"), n_nk_cells=0, dz_overlap_median_m=float("nan"),
                n_overlap_cells=0, nmad_overlap_m=float("nan"), overlap_area_km2=0.0,
                nmad_before_m=float("nan"), nmad_after_m=float("nan"), converged=False,
                degenerate=True, dtan_min=float("nan"), dtan_max=float("nan"),
                dtan_median=float("nan"), c_mm_per_tan=float("nan"),
                k_check_m=float("nan"), extrapolated=False,
                n_points_ref=0, n_points_src=0, tiles=[str(e)]))
    return WindowLadder(ties=out, line_ref=int(line_ref), line_src=int(line_src),
                        easting=float(easting), northing=float(northing),
                        half_widths_m=[float(h) for h in half_widths_m])


# ------------------------------------------------------------------- local network

@dataclass
class LocalNetwork:
    """Every line in one window, solved into one line's frame at that window.

    This is ``coreg.align_swaths`` run on a window instead of a tile: the same free-network
    least squares over the same pairwise observations, gauged on ``ref_line``.
    ``constants_m[s]`` is what to ADD to line ``s``'s z to reach ``ref_line``'s frame here.
    """

    ref_line: int
    easting: float
    northing: float
    half_width_m: float
    shape: str
    tie_mode: str
    res_m: float
    constants_m: dict
    edges: list
    misclosure_mm: np.ndarray
    line_counts: dict
    tiles: list
    warnings: list = field(default_factory=list)

    def constants_mm(self) -> dict:
        return {s: v[2] * 1000.0 for s, v in self.constants_m.items()}


def local_network(tile_paths, lines, *, easting, northing, half_width_m, shape,
                  res_m, tie, exclude, ref_line,
                  cache: TileCache | None = None) -> LocalNetwork:
    """Solve all of ``lines`` into ``ref_line``'s frame using only this window's overlaps.

    Reuses :func:`lidar_diff_icp.coreg.align_swaths` unchanged -- the gauge, the weighting
    and the misclosure are its own. On a set of adjacent parallel lines the overlap graph
    is a chain, so the misclosure is identically zero and carries no information; it is
    returned anyway because it *is* informative wherever a loop exists (a cross line, or a
    pair that overlaps non-adjacently).
    """
    import warnings as _w
    w = window_cloud(tile_paths, easting=easting, northing=northing,
                     half_width_m=half_width_m, shape=shape, cache=cache, lines=lines)
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        corr, edges, mis = coreg.align_swaths(w.pc, res_m, tuple(exclude),
                                              ref=int(ref_line), tie=tie)
    return LocalNetwork(ref_line=int(ref_line), easting=w.easting, northing=w.northing,
                        half_width_m=w.half_width_m, shape=shape, tie_mode=tie,
                        res_m=float(res_m),
                        constants_m={int(k): tuple(float(q) for q in v)
                                     for k, v in corr.items()},
                        edges=edges, misclosure_mm=np.asarray(mis) * 1000.0,
                        line_counts=w.line_counts, tiles=w.tiles,
                        warnings=[str(c.message) for c in caught])


# ------------------------------------------------------------------------- chaining

def _cell_ids(x, y, res):
    """Cell ids on a GLOBAL grid anchored at (0, 0), so ids from different tiles agree.

    Written out rather than imported because ``groundtruth.chain._cells`` is private and
    that module belongs to another author; the encoding is deliberately identical so the
    two modules' overlap cells are the same cells.
    """
    return (np.floor(y / res).astype(np.int64) << 32) + np.floor(x / res).astype(np.int64)


@dataclass
class OverlapPoint:
    """Where a pair can be tied closest to a target location."""

    line_a: int
    line_b: int
    easting: float
    northing: float
    distance_m: float        # from the target location to that overlap cell centre
    n_overlap_cells: int
    median_easting: float    # centre of the pair's overlap in these tiles ...
    median_northing: float   # ... which is NOT where the nearest cell is: on a N-S strip
                             # the nearest cell to a point west of it sits on its west EDGE,
                             # and an edge window samples the across-track term one-sidedly
    tiles: list


def nearest_overlap_point(tile_paths, line_a: int, line_b: int, *, easting, northing,
                          res_m, exclude, cache: TileCache | None = None) -> OverlapPoint:
    """The point of the pair's overlap closest to ``(easting, northing)``.

    A chain from a mark's line to a target line moves sideways across the flight
    direction, so the far links cannot be solved *at* the mark; they can only be solved at
    the nearest place the pair actually overlaps. Making that displacement an explicit
    output is the point: it is the distance over which the link's tie is being assumed
    constant, and it is exactly the quantity that was left unstated when Elba's constants
    were applied 62.9 km away.
    """
    cache = TileCache() if cache is None else cache
    ids_a, ids_b = [], []
    used = []
    for path in [str(p) for p in tile_paths]:
        a = cache.tile(path)
        terr = ~np.isin(a["cls"], tuple(exclude))
        ma = terr & (a["psid"] == line_a)
        mb = terr & (a["psid"] == line_b)
        if not ma.any() or not mb.any():
            continue
        ids_a.append(np.unique(_cell_ids(a["x"][ma], a["y"][ma], res_m)))
        ids_b.append(np.unique(_cell_ids(a["x"][mb], a["y"][mb], res_m)))
        used.append(path)
    if not ids_a or not ids_b:
        raise ValueError(f"lines {line_a} and {line_b} share no tile with terrain returns")
    inter = np.intersect1d(np.unique(np.concatenate(ids_a)),
                           np.unique(np.concatenate(ids_b)))
    if inter.size == 0:
        raise ValueError(f"lines {line_a} and {line_b} do not overlap in these tiles")
    cx = (inter & 0xFFFFFFFF).astype(np.int64)
    cx = np.where(cx >= 1 << 31, cx - (1 << 32), cx)          # sign-extend the low word
    cy = inter >> 32
    ex = (cx + 0.5) * res_m
    ny_ = (cy + 0.5) * res_m
    d = np.hypot(ex - easting, ny_ - northing)
    k = int(np.argmin(d))
    return OverlapPoint(line_a=int(line_a), line_b=int(line_b), easting=float(ex[k]),
                        northing=float(ny_[k]), distance_m=float(d[k]),
                        n_overlap_cells=int(inter.size),
                        median_easting=float(np.median(ex)),
                        median_northing=float(np.median(ny_)), tiles=used)


@dataclass
class ChainLink:
    tie: LocalTie
    solve_easting: float
    solve_northing: float
    solve_distance_m: float       # mark -> where this link had to be solved
    ladder: WindowLadder | None = None


@dataclass
class LocalChain:
    """A mark's line carried into a target line's frame, link by link, near the mark."""

    easting: float
    northing: float
    source_line: int
    target_line: int
    nodes: list
    links: list
    dz_total_m: float
    dz_sigma_formal_m: float          # per-link coreg sigmas in quadrature: optimistic
    dz_sigma_window_m: float          # from the per-link window ladders, when requested
    half_width_m: float
    tie_mode: str
    res_m: float
    max_solve_distance_m: float
    degenerate_links: list
    params: list = field(default_factory=list)

    @property
    def dz_total_mm(self) -> float:
        return self.dz_total_m * 1000.0

    def rows(self):
        return [[f"{l.tie.line_ref}-{l.tie.line_src}",
                 f"{l.solve_easting:.0f}", f"{l.solve_northing:.0f}",
                 f"{l.solve_distance_m:.0f}", f"{l.tie.dz_m * 1000:+.1f}",
                 f"{l.tie.dz_sigma_m * 1000:.1f}",
                 f"{l.tie.dz_overlap_median_m * 1000:+.1f}", l.tie.n_nk_cells,
                 ("" if l.ladder is None else f"{l.ladder.spread_mm:.1f}")]
                for l in self.links]

    @staticmethod
    def columns() -> dict:
        return {
            "link": "flight-line pair, ref-src",
            "solve_E": "easting of the window this link was solved in, m",
            "solve_N": "northing of that window, m",
            "dist_m": "distance from the mark to that window centre, m -- the "
                      "displacement over which this link's tie is assumed constant",
            "dz_mm": "the link's tie, mm",
            "dz_sig_mm": "coreg's formal 1-sigma for it, mm",
            "dz_median_mm": "independent overlap-median read of the same link, mm",
            "nk_cells": "cells the Nuth & Kaeaeb fit used; 0 means abandoned",
            "ladder_mm": "spread of this link's tie over the requested window ladder, mm",
        }


def chain_local(tile_paths, *, easting, northing, source_line, target_line,
                half_width_m, shape, res_m, tie, exclude, path=None,
                ladder_half_widths_m=None, cache: TileCache | None = None) -> LocalChain:
    """Chain ``source_line`` into ``target_line``'s frame, each link solved near the mark.

    ``path`` is the ordered list of lines to walk, ``[source_line, ..., target_line]``.
    When it is ``None`` the route is planned by ``groundtruth.chain`` -- its
    :func:`~lidar_diff_icp.groundtruth.chain.plan_path` breadth-first search over the
    measured overlap graph, which minimises **link count** because each link adds error --
    on an inventory built from this module's cache (see :func:`inventory_from_cache`).
    That module is imported and reused, not modified.

    Each link is solved in a window of ``half_width_m`` centred on
    :func:`nearest_overlap_point` for that pair: the closest place to the mark where the
    pair can be tied at all. ``ladder_half_widths_m``, when given, runs a
    :class:`WindowLadder` per link as well and returns their spreads summed in quadrature
    as ``dz_sigma_window_m`` -- an empirical alternative to the formal sigma, which
    ``analysis/ABSOLUTE_BASIS_ELBA.md`` §2 shows to be optimistic by roughly an order of
    magnitude for this quantity.
    """
    cache = TileCache() if cache is None else cache
    tile_paths = [str(p) for p in tile_paths]
    source_line = int(source_line); target_line = int(target_line)
    if path is None:
        nodes = plan_path_local(tile_paths, source_line, target_line,
                                res_m=res_m, exclude=exclude, cache=cache)
    else:
        nodes = [int(n) for n in path]
        if nodes[0] != source_line or nodes[-1] != target_line:
            raise ValueError(f"path {nodes} does not run {source_line} -> {target_line}")
    links = []
    dz = 0.0
    var_formal = 0.0
    var_window = 0.0
    have_ladder = ladder_half_widths_m is not None
    # Walk from the TARGET inward, so each step's `dz` adds the next line into the frame
    # already accumulated -- the sign convention of coreg.coregister_swaths(ref=a, src=b).
    walk = nodes[::-1]
    for a, b in zip(walk[:-1], walk[1:]):
        op = nearest_overlap_point(tile_paths, a, b, easting=easting, northing=northing,
                                   res_m=res_m, exclude=exclude, cache=cache)
        t = pair_tie_at(tile_paths, a, b, easting=op.easting, northing=op.northing,
                        half_width_m=half_width_m, shape=shape, res_m=res_m, tie=tie,
                        exclude=exclude, cache=cache)
        lad = None
        if have_ladder:
            lad = window_ladder(tile_paths, a, b, easting=op.easting,
                                northing=op.northing, half_widths_m=ladder_half_widths_m,
                                shape=shape, res_m=res_m, tie=tie, exclude=exclude,
                                cache=cache)
            if np.isfinite(lad.spread_mm):
                var_window += (lad.spread_mm / 1000.0) ** 2
        links.append(ChainLink(tie=t, solve_easting=op.easting, solve_northing=op.northing,
                               solve_distance_m=op.distance_m, ladder=lad))
        dz += t.dz_m
        var_formal += t.dz_sigma_m ** 2
    params = [
        Param("half_width_m", float(half_width_m), "caller",
              "window half-width for every link; no library default exists"),
        Param("res_m", float(res_m), "caller",
              "grid resolution passed to coreg.coregister_swaths"),
        Param("tie", tie, "caller",
              "coreg tie mode: 'overlap_median' or 'intercept'"),
        Param("exclude", tuple(exclude), "caller",
              "classes excluded from the terrain proxy"),
        Param("shape", shape, "caller", "window shape: 'square' or 'disk'"),
        Param("path", list(nodes), "caller" if path is not None else "repo",
              "walked route source -> target"
              + ("" if path is not None else "; fewest-link route from "
                 "groundtruth.chain.plan_path over the measured overlap graph")),
    ]
    return LocalChain(
        easting=float(easting), northing=float(northing), source_line=source_line,
        target_line=target_line, nodes=nodes, links=links, dz_total_m=float(dz),
        dz_sigma_formal_m=float(np.sqrt(var_formal)),
        dz_sigma_window_m=float(np.sqrt(var_window)) if have_ladder else float("nan"),
        half_width_m=float(half_width_m), tie_mode=tie, res_m=float(res_m),
        max_solve_distance_m=max((l.solve_distance_m for l in links), default=0.0),
        degenerate_links=[(l.tie.line_ref, l.tie.line_src) for l in links
                          if l.tie.degenerate],
        params=params)


def inventory_from_cache(tile_paths, *, exclude, cache: TileCache):
    """A ``groundtruth.chain.SwathInventory`` backed by this module's already-read tiles.

    ``chain.build_inventory`` would decompress every tile a second time into its own cache
    format. This fills the same structure from arrays already in memory, so
    ``chain.overlap_graph`` and ``chain.plan_path`` can be reused as they are.
    """
    from .groundtruth import chain as gt_chain

    inv = gt_chain.SwathInventory(cache_dir=None, exclude=tuple(exclude))
    for p in [str(q) for q in tile_paths]:
        a = cache.tile(p)
        keep = ~np.isin(a["cls"], tuple(exclude))
        x = a["x"][keep]; y = a["y"][keep]; z = a["z"][keep]; ps = a["psid"][keep]
        lines = {}
        for s in np.unique(ps):
            m = ps == s
            lines[int(s)] = {"n": int(m.sum()),
                             "bbox": [float(x[m].min()), float(y[m].min()),
                                      float(x[m].max()), float(y[m].max())]}
        inv.tiles[p] = gt_chain.TileLines(p, lines, int(a["cls"].size), int(keep.sum()))
        inv._points[p] = (x, y, z, ps)
    return inv


def plan_path_local(tile_paths, source_line: int, target_line: int, *, res_m, exclude,
                    cache: TileCache) -> list:
    """Fewest-link route from ``source_line`` to ``target_line``, via ``groundtruth.chain``."""
    from .groundtruth import chain as gt_chain

    inv = inventory_from_cache(tile_paths, exclude=exclude, cache=cache)
    graph = gt_chain.overlap_graph(inv, res=res_m)
    paths = gt_chain.plan_path(graph, inv, [int(source_line)], [int(target_line)])
    if not paths:
        raise ValueError(
            f"no overlap route from line {source_line} to line {target_line} in "
            f"{len(list(tile_paths))} tile(s)")
    return list(paths[0].nodes)


# --------------------------------------------------- comparing against fitted constants

@dataclass
class TieComparison:
    """A local chain against a constant fitted somewhere else -- e.g. Elba's."""

    local_mm: float
    imported_mm: float
    difference_mm: float
    source_line: int
    target_line: int
    easting: float
    northing: float
    max_solve_distance_m: float
    note: str = ""


def compare_to_constants(chain: LocalChain, constants_dzm: dict, *, note="") -> TieComparison:
    """Difference between a local chain and an imported per-line constant set.

    ``constants_dzm`` maps line id -> the vertical correction ADDED to that line by
    whatever fit produced it (the ``dz`` component of
    ``per_swath_internal_alignment_dxdydz_m``). The imported prediction for putting the
    source line into the target line's frame is ``dz[source] - dz[target]``; the gauge the
    imported set was solved on cancels in that difference, so two products gauged on
    different lines are compared correctly.
    """
    s = float(constants_dzm[chain.source_line])
    t = float(constants_dzm[chain.target_line])
    imported = (s - t) * 1000.0
    return TieComparison(local_mm=chain.dz_total_mm, imported_mm=imported,
                         difference_mm=chain.dz_total_mm - imported,
                         source_line=chain.source_line, target_line=chain.target_line,
                         easting=chain.easting, northing=chain.northing,
                         max_solve_distance_m=chain.max_solve_distance_m, note=note)
