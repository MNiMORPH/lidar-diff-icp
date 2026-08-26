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
from .provenance import Param, declare
__all__ = [
    "Checkpoint", "CheckpointSet", "UnknownDatumError", "load_bundled",
    "read_3dep_va_shapefile", "list_bundled",
    "Param", "declare",
]
