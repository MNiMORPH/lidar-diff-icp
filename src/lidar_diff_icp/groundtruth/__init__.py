"""Pin a floating lidar vertical datum to surveyed ground control.

The 2008 MN DNR lidar (gen1) carries no absolute vertical tie at Elba: the pipeline
registers gen1 to gen2 laterally, applies the GEOID03 -> GEOID18 datum shift, and
aligns the flight lines to one another, but the whole *group* is still free to float
in z (see :func:`lidar_diff_icp.coreg.align_swaths` -- "the group's absolute offset
from another epoch must be tied separately").

This package supplies that tie from surveyed checkpoints, in three separable parts:

``checkpoints``
    Load surveyed control with its datum metadata explicit -- horizontal CRS, vertical
    datum, geoid model, units, point type (NVA on open ground / VVA under vegetation).
    Refuses to hand back a point whose geoid model is unknown.

``tie``
    Given one checkpoint and a lidar cloud, return the lidar ground elevation there with
    an uncertainty, a **radius curve**, and diagnostics. The estimator is the project's
    slope-normal ground read (:func:`lidar_diff_icp.pipeline.difference_dem`,
    ``ground="slope_normal"``) generalised from a cell centre to an arbitrary point.

``datum``
    Combine independent ties into ONE constant with an uncertainty budget that keeps
    common-mode error at full size: a term shared by every tie (the extrapolated lateral
    shift, the alignment estimator's repeatability) does not average down with the number
    of marks, and an unmodelled gap is reported beside the total rather than inside it.

``gen1_datum``
    gen1's absolute vertical datum at ANY Minnesota site, from the 2008 MnGeo control
    the acquisition was itself validated on -- same datum, same geoid, no conversion and
    no cross-epoch term. Discovers marks near a site or a set of flight lines, names the
    tiles (never fetches them), assigns each mark to a flight line FROM THE RETURNS, and
    combines with the flight line as the unit of replication.

``chain``
    Flight lines that cover a checkpoint are usually not the lines that cover the study
    area. Overlapping swaths chain, so a tie measured on one line propagates. The solver
    tries the zero-link along-swath case first and only then searches for the shortest
    chain, because every link adds error and a chain has no internal redundancy.

Method, sign conventions and the radius pathology: ``docs/groundtruth.md``.
"""
from __future__ import annotations

from .checkpoints import (
    Checkpoint,
    CheckpointSet,
    UnknownDatumError,
    load_bundled,
    read_3dep_va_shapefile,
    list_bundled,
)
from .datum import BudgetTerm, DatumConstant, combine_ties
from .provenance import Param, declare
from .tie import (
    GroundReturns,
    RadiusEstimate,
    TieEstimate,
    csf_ground_near,
    estimate_tie,
    geoid_shift_for,
    ground_elevation_at,
    radius_ladder,
    scan_angle_deg,
    vendor_ground_near,
)
from .gen1_datum import (
    ControlMark,
    ControlSet,
    DatumMismatchError,
    Gen1DatumEstimate,
    LineAssignment,
    LineGroup,
    MarkMeasurement,
    MarkSite,
    SitingScreen,
    TileNeed,
    TileResolution,
    assert_no_geoid_conversion,
    assign_line_from_returns,
    combine_datum,
    discover_near_lines,
    discover_near_point,
    load_control,
    measure_site,
    measure_sites,
    resolve_tiles,
    swath_constants_from_corrections,
)
from .chain import (
    ChainPath,
    ChainSolution,
    Link,
    SwathInventory,
    build_inventory,
    compare_paths,
    covering_lines,
    overlap_graph,
    plan_path,
    solve_chain,
    solve_link,
)
__all__ = [
    "Checkpoint", "CheckpointSet", "UnknownDatumError", "load_bundled",
    "read_3dep_va_shapefile", "list_bundled",
    "BudgetTerm", "DatumConstant", "combine_ties",
    "Param", "declare",
    "GroundReturns", "RadiusEstimate", "TieEstimate", "csf_ground_near", "estimate_tie",
    "geoid_shift_for", "ground_elevation_at", "radius_ladder", "scan_angle_deg",
    "vendor_ground_near",
    "ControlMark", "ControlSet", "DatumMismatchError", "Gen1DatumEstimate",
    "LineAssignment", "LineGroup", "MarkMeasurement", "MarkSite", "SitingScreen",
    "TileNeed", "TileResolution", "assert_no_geoid_conversion",
    "assign_line_from_returns", "combine_datum", "discover_near_lines",
    "discover_near_point", "load_control", "measure_site", "measure_sites",
    "resolve_tiles", "swath_constants_from_corrections",
    "ChainPath", "ChainSolution", "Link", "SwathInventory", "build_inventory",
    "compare_paths", "covering_lines", "overlap_graph", "plan_path", "solve_chain",
    "solve_link",
]
