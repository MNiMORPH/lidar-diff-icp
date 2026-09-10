"""How many of a site's own survey control marks are REACHABLE in tiles on disk.

The gen1 side had no equivalent of :mod:`completeness`. That gap is not theoretical: on
2026-09-10 mnrv built completely, passed all six pipeline steps, and could not have a datum
at all -- **zero** control marks, of any cover, fall inside its single gen1 tile -- and
nothing in the graph said so. It took a failed run to find out. elba and whitewater reach
29 open marks because 47 tiles were fetched by hand in an earlier session; mnrv reaches 0
because they were not.

**Reported, never judged**, exactly as ``completeness`` is: this module counts and names,
and applies no threshold. How many marks a datum needs is a scientific decision, and the
radius follows from it rather than the other way round.

WHY TILES UNDER MARKS, NOT A RADIUS OF COVERAGE
-----------------------------------------------
A mark's tie needs a small window: ``measure_site`` derives ``crop_half_width_m`` as
``5*res``, a 50 m box at the 5 m grid, and the bridge's CSF reconstruction needs 300 m.
Against a 2410 m tile pitch, that is a few hundred metres of a ~22 MB tile. But MnGeo
serves plain LAZ -- HTTP range reads are used for the 512-byte header only, and there is no
COPC or .lax index -- so a whole tile is the unit of download.

That makes covering a DISC wasteful and covering the MARKS cheap. At mnrv, 11 open marks
within 25 km sit in 11 distinct tiles, ~242 MB; the 25 km disc would be ~85 tiles and
~1.9 GB. So this reports the tiles the marks are IN, by name, and a fetch pulls exactly
those. It is what was done by hand for SE-MN.

The survey is the Site's own (:mod:`lidar_diff_icp.acquisitions`), never a default: gen1 is
four acquisitions on two geoid models, and a mark from the wrong survey is not a weaker
tie, it is a different datum.
"""
from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field
from pathlib import Path

__all__ = ["ControlReach", "reach_for_site", "summary_line"]


@dataclass(frozen=True)
class ControlReach:
    """What a site can and cannot reach, with the tiles named."""

    site: str
    project_id: str
    covers: tuple[str, ...]
    n_marks_in_survey: int
    n_marks_considered: int          # after the cover restriction, before any radius
    radius_m: float | None
    n_marks_in_radius: int
    on_disk: tuple[str, ...] = ()    # tiles holding a considered mark, present locally
    missing: tuple[str, ...] = ()    # tiles holding a considered mark, NOT present
    marks_reachable: int = 0         # marks whose tile is on disk
    tile_dirs: tuple[str, ...] = ()
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def approx_missing_mb(self) -> int:
        """Rough download size. 22 MB is the observed size of gen1 tiles on disk
        (20-24 MB across data/before), not a specification."""
        return 22 * len(self.missing)


def reach_for_site(site, *, covers=("L1O",), radius_m=None, tile_dirs=None,
                   cache=None) -> ControlReach:
    """Which of ``site``'s survey's marks are in tiles on disk, and which tiles are not.

    ``covers``    restrict to these cover classes. Defaults to open ground only, which is
                  the standing rule for a datum -- pooling vegetated marks bakes canopy
                  response into the level. Pass ``None`` to count every cover.
    ``radius_m``  optional distance cut from the site centre. There is NO default: how far
                  to reach is a scientific choice about how many marks the datum needs, and
                  the count at each radius is what this exists to show.
    ``tile_dirs`` where gen1 tiles live. Defaults to the directory of the Site's own gen1.
    """
    from . import acquisitions, tiles as T
    from .groundtruth import gen1_datum as G

    acq = acquisitions.for_project(site.gen1_project)     # raises on an unknown survey
    cset = G.control_for_survey(site.gen1_project)
    marks = list(cset.marks)
    n_all = len(marks)
    if covers is not None:
        marks = [m for m in marks if m.cover_class in set(covers)]
    n_considered = len(marks)

    b = site.bounds
    if b is None:
        import json
        b = json.load(open(os.path.join(site.tile_dir, "corrections.json")))["bounds"]
    cx, cy = 0.5 * (b[0] + b[2]), 0.5 * (b[1] + b[3])
    if radius_m is not None:
        marks = [m for m in marks
                 if ((m.checkpoint.easting - cx) ** 2
                     + (m.checkpoint.northing - cy) ** 2) ** 0.5 <= float(radius_m)]

    dirs = tuple(tile_dirs) if tile_dirs else (os.path.dirname(site.gen1),)
    have = set()
    for d in dirs:
        for p in glob.glob(os.path.join(d, "*.laz")):
            have.add(Path(p).stem)

    kw = {} if cache is None else {"cache": cache}
    on_disk, missing, reachable = set(), set(), 0
    for m in marks:
        name = T.find_tile(m.checkpoint.easting, m.checkpoint.northing, **kw)
        # A metro tile is stored under a SUFFIXED filename (find_tile's own docstring),
        # so match by prefix as well as exact stem or a present tile reads as missing.
        hit = next((h for h in have if h == name or h.startswith(name)), None)
        if hit:
            on_disk.add(hit); reachable += 1
        else:
            missing.add(name)

    notes = []
    if n_considered and not marks:
        notes.append(f"no mark of cover {covers} within radius_m={radius_m}")
    if covers is not None and n_considered < n_all:
        notes.append(f"{n_all - n_considered} of {n_all} marks excluded by cover "
                     f"{covers}; they are counted, not discarded")
    return ControlReach(
        site=site.name, project_id=acq.project_id,
        covers=tuple(covers) if covers else (), n_marks_in_survey=n_all,
        n_marks_considered=n_considered, radius_m=radius_m,
        n_marks_in_radius=len(marks), on_disk=tuple(sorted(on_disk)),
        missing=tuple(sorted(missing)), marks_reachable=reachable,
        tile_dirs=dirs, notes=tuple(notes))


def summary_line(r: ControlReach) -> str:
    """One line, in the shape completeness.summary_line prints."""
    rad = "no radius cut" if r.radius_m is None else f"within {r.radius_m/1000:.0f} km"
    return (f"[{r.site}] gen1 control reach: {r.marks_reachable} of {r.n_marks_in_radius} "
            f"{'/'.join(r.covers) or 'all-cover'} marks {rad} are in tiles on disk "
            f"({r.project_id}); {len(r.missing)} tile(s) missing "
            f"~{r.approx_missing_mb} MB")
