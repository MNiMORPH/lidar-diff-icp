"""Measure BOTH epochs' datum constants at any MN site, ready for the pipeline.

This is the reusable method the brief asked for. It runs, for one site, exactly what was
run at Elba:

  gen1  every 2008 control mark inside a tile on disk, assigned to a flight line by its
        OWN ground returns, restricted to the site's lines, reused point_source_ids
        disambiguated by track collinearity, combined with the LINE as the unit of
        replication -- then carried onto our surface with the bridge.
  gen2  its own 2021 held-out checkpoints, kriged to the site on the chosen delivered
        surface.
  geoid the GEOID03->GEOID18 term the pipeline adds to gen1, from the PROJ grids.

and returns the ``absolute_datum`` dict ``pipeline.difference_dem`` accepts.

Nothing here has a default that could silently set the answer. ``covers``,
``collinear_sigma``, ``gen2_surface``, ``bridge_mm`` and ``gauge_ref`` are all required:
each of them moved the Elba answer by more than the correction itself.

**The bridge is the caller's.** Measuring it needs CSF at every mark (~9 min at Elba), so
it is an argument, not something this function invents. Pass a value measured AT THIS SITE,
or pass Elba's with the assumption stated -- ``bridge_source`` records which.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, asdict, field
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "src"))

import apply_datum as AD  # noqa: E402
import datum as DAT  # noqa: E402
import lines as L  # noqa: E402
import same_line as S  # noqa: E402
from lidar_diff_icp import references  # noqa: E402
from lidar_diff_icp.groundtruth import gen1_datum as G  # noqa: E402


@dataclass
class SiteDatum:
    """Both epochs' constants at one site, with every input that set them."""

    site: tuple
    gen1_delivered_mm: float
    gen1_delivered_se_mm: float
    gen1_n_marks: int
    gen1_n_lines: int
    bridge_mm: float
    bridge_source: str
    gen1_our_surface_mm: float
    gen2_delivered_mm: float
    gen2_surface: str
    geoid_mm: float
    gauge_ref: int
    dod_shift_mm: float
    params: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def to_pipeline(self, source: str) -> dict:
        return AD.datum_for_pipeline(
            self.gen1_our_surface_mm, self.geoid_mm, self.gen2_delivered_mm,
            self.gauge_ref, source, gen1_sigma_mm=self.gen1_delivered_se_mm)

    def to_dict(self):
        return asdict(self)


def measure_site_datum(*, easting, northing, bounds, crs, psids, tile_dirs, tracks_path,
                       covers, collinear_sigma, res, gen2_surface, bridge_mm,
                       bridge_source, gauge_ref, max_lags_m, n_lags, n_pairs, estimators,
                       seed) -> SiteDatum:
    """Measure gen1's and gen2's constants at a site. Every argument is required."""
    ts = L.load_tracks(tracks_path)
    control = G.load_control()

    meas, kept, rejected, est = S.estimate_by_returns(
        ts, psids=psids, easting=easting, northing=northing, covers=covers,
        tile_dirs=tile_dirs, res=res, collinear_sigma=collinear_sigma, control=control)

    warn = []
    if est.n_lines < 3:
        warn.append(f"only {est.n_lines} flight line(s) carry a mark: the SE over lines "
                    f"is barely defined")
    if est.n_marks < 5:
        warn.append(f"only {est.n_marks} marks survived; at Elba 8 marks gave +/-23 mm")
    missing = sorted(set(int(p) for p in psids) - {int(m.line_id) for m in kept})
    if missing:
        warn.append(f"lines with NO mark: {missing} -- their level is unconstrained and "
                    f"they may carry substantial weight in the tile's ground")

    g2 = DAT.datum_at_site("gen2", easting=easting, northing=northing, treatment="open",
                           surface=gen2_surface, max_lags_m=max_lags_m, n_lags=n_lags,
                           n_pairs=n_pairs, estimators=estimators, seed=seed)
    a0, bx, cy = references.geoid_difference(list(bounds), crs)
    cx, cyy = (bounds[0] + bounds[2]) / 2.0, (bounds[1] + bounds[3]) / 2.0
    geoid_mm = (a0 + bx * (easting - cx) / 1000.0
                + cy * (northing - cyy) / 1000.0) * 1000.0

    c1_ours = est.value_mm + float(bridge_mm)
    return SiteDatum(
        site=(float(easting), float(northing)),
        gen1_delivered_mm=float(est.value_mm), gen1_delivered_se_mm=float(est.se_mm),
        gen1_n_marks=int(est.n_marks), gen1_n_lines=int(est.n_lines),
        bridge_mm=float(bridge_mm), bridge_source=str(bridge_source),
        gen1_our_surface_mm=float(c1_ours),
        gen2_delivered_mm=float(g2.constant_mm), gen2_surface=str(gen2_surface),
        geoid_mm=float(geoid_mm), gauge_ref=int(gauge_ref),
        dod_shift_mm=float(g2.constant_mm - (c1_ours - geoid_mm)),
        params=dict(covers=list(covers), collinear_sigma=collinear_sigma, res_m=res,
                    psids=list(psids), gen2_surface=gen2_surface,
                    max_lags_m=list(max_lags_m), n_lags=n_lags, n_pairs=n_pairs,
                    estimators=list(estimators), seed=seed),
        warnings=warn)
