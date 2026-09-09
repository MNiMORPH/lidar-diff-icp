"""Which SURVEY a gen1 tile came from, and the vertical frame that survey was reduced in.

gen1 is not one dataset. The six pilot sites fall in **four acquisitions on two geoid
models**, and until 2026-09-07 the pipeline differenced every one of them GEOID03 ->
GEOID18 because :func:`references.geoid_difference` defaulted ``before_geoid`` to GEOID03
and ``pipeline.difference_dem`` never overrode it. Measured cost of that single default,
as the constant added to gen1:

    battlecreek  lidar_metro2011      applied +71.85  correct +16.98  ERROR +54.87 mm
    carlton      lidar_duluth2012     applied +78.14  correct +50.39  ERROR +27.75 mm
    cook         lidar_arrowhead2011  applied +51.68  correct +25.29  ERROR +26.39 mm

The datum is ADDED to gen1, so adding too much makes gen1 read high and the DoD
(``gen2 - gen1``) read LOW by the error. At Battle Creek that is 61% of the site's own
90 mm LoD, and larger than its entire zero-line lever (37.50 mm).

**The geoid is never defaulted here.** Every lookup raises on an acquisition it does not
know, because the failure mode is silent: a wrong geoid does not look like an error from
inside the tile, it looks like tens of millimetres of erosion.

**Provenance.** The validation reports themselves state no datum -- this is a DATASET-level
assertion, and each record therefore carries the metadata page it was read from and the
sentence that asserts it, so the claim can be checked without re-deriving it.

``counties`` lists only the counties whose metadata we have READ. It is not the project's
full extent: MnGeo organises by county, and we have confirmed one or eight per project as
the sites required. Adding a site means reading its county's page and extending the tuple.
"""
from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Acquisition", "ACQUISITIONS", "GEOID_GRID", "for_project", "for_county",
           "GEN2_3DEP_GEOID", "GEN2_3DEP_GRID"]

#: PROJ grid names for the NAVD88 realizations gen1 was reduced in.
GEOID_GRID = {"GEOID03": "us_noaa_geoid03_conus.tif",
              "GEOID09": "us_noaa_geoid09_conus.tif",
              "GEOID12B": "us_noaa_g2012bu0.tif",
              "GEOID18": "us_noaa_g2018u0.tif"}

#: gen2 is USGS 3DEP 2021 at every site we hold, and 3DEP publishes on GEOID18. Named
#: here so the after-epoch frame is a stated fact rather than a default buried in a
#: signature. If a gen2 ever arrives on another realization this is where it goes.
GEN2_3DEP_GEOID = "GEOID18"
GEN2_3DEP_GRID = GEOID_GRID[GEN2_3DEP_GEOID]


@dataclass(frozen=True)
class Acquisition:
    """One gen1 lidar survey, and the evidence for the frame it is on."""

    project_id: str
    geoid_model: str
    collected: str
    metadata_page: str
    datum_quote: str
    counties: tuple[str, ...]
    #: The MnGeo validation reports that carry this survey's checkpoint tables, by the
    #: stem the parser matches. NOT always per-county: the Arrowhead and Duluth surveys
    #: file theirs under `projects/<name>/`, by BLOCK, which is why the county
    #: directories look as though no report exists. Empty means none has been located.
    report_regions: tuple[str, ...] = ()

    @property
    def geoid_grid(self) -> str:
        """The PROJ grid name for this survey's geoid, for ``geoid_difference``."""
        try:
            return GEOID_GRID[self.geoid_model]
        except KeyError:
            raise KeyError(
                f"{self.project_id} is recorded on {self.geoid_model}, which has no PROJ "
                f"grid in GEOID_GRID. Add it rather than substituting a near neighbour: "
                f"GEOID03 and GEOID09 differ by 54.87 mm at Ramsey.") from None


ACQUISITIONS: dict[str, Acquisition] = {
    a.project_id: a for a in (
        Acquisition(
            "lidar_semn2008", "GEOID03", "2008", "lidar_semn2008.html",
            "Vertical datum: NAVD88 (Geoid03)",
            ("dodge", "fillmore", "houston", "mower", "olmsted", "steele",
             "wabasha", "winona"),
            report_regions=("dodge", "fillmore", "houston", "mower", "olmsted",
                            "steele", "wabasha", "winona")),
        Acquisition(
            "lidar_swmn2010", "GEOID03", "2010", "lidar_swmn2010.html",
            "The NAVD88, Geoid03 vertical datum was used.",
            ("lesueur",), report_regions=("lesueur",)),
        Acquisition(
            "lidar_metro2011", "GEOID09", "2011", "lidar_metro2011.html",
            "The NAVD88 (Geoid09) vertical datum was used.",
            ("ramsey",), report_regions=("ramsey",)),
        Acquisition(
            "lidar_arrowhead2011", "GEOID09", "2011", "lidar_arrowhead2011.html",
            "The geoid used to reduce satellite derived elevations to orthometric "
            "heights was Geoid09.",
            ("cook",),
            # projects/arrowhead/block_3/Arrowhead_block_3_validation_report.pdf --
            # nearest mark 5.08 km from the cook tile. Blocks 1,2,4,5 are 88-161 km away.
            report_regions=("arrowhead_block3",)),
        Acquisition(
            "lidar_duluth2012", "GEOID09", "2012", "lidar_duluth2012.html",
            "Lidar data are in the UTM Zone 15 coordinate system, NAD83 96, NAVD88 "
            "Geoid09 meters tiled by USGS 1/16, 1:24,000 quadrangles.",
            ("carlton",),
            # projects/duluth_fall_2012/duluth_2012_vertical_validation_report.pdf --
            # 508 checkpoints, nearest 3.10 km from the carlton tile.
            report_regions=("duluth2012",)),
    )
}

_BY_COUNTY = {c: a for a in ACQUISITIONS.values() for c in a.counties}
#: report-stem -> Acquisition, for the checkpoint parser. A "region" is a county for the
#: surveys MnGeo files that way and a project block for the ones it does not.
BY_REPORT_REGION = {r: a for a in ACQUISITIONS.values() for r in a.report_regions}


def for_project(project_id: str) -> Acquisition:
    """The acquisition record, or raise. Never falls back to a default geoid."""
    try:
        return ACQUISITIONS[project_id]
    except KeyError:
        raise KeyError(
            f"no acquisition recorded for {project_id!r}. Add it to ACQUISITIONS with the "
            f"geoid model and the sentence from its MnGeo metadata page that asserts the "
            f"datum. Do NOT default: applying GEOID03 to a Geoid09 survey moved gen1 by "
            f"+54.87 mm at Ramsey, unnoticed, until 2026-09-07. "
            f"Known: {sorted(ACQUISITIONS)}") from None


def for_county(county: str) -> Acquisition:
    """The acquisition covering a MnGeo county directory name, or raise."""
    try:
        return _BY_COUNTY[county]
    except KeyError:
        raise KeyError(
            f"no acquisition recorded for county {county!r}. Read that county's MnGeo "
            f"metadata page and extend the covering Acquisition's `counties`. "
            f"Known: {sorted(_BY_COUNTY)}") from None
