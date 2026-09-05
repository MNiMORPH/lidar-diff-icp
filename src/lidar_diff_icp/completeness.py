"""Is a site's point cloud COMPLETE, and does the record say so?

WHY THIS EXISTS. A truncated fetch does not look like an error. Whitewater's gen2 file
averaged 11.39 returns per m2 over a tile that was 15.45 west and 5.52 east of a seam --
it would have passed any tile-average density test, and it produced a DoD with a
flight-line-shaped artifact that read as terrain. Elba was worse: its gen2 was 10.90x
thinner than every other site's, 19.27x on the ground class, so one site was differenced
from a different instrument than its peers. Neither was visible from inside the tile.

WHAT THIS MODULE DOES, AND WHAT IT DELIBERATELY DOES NOT. It records a completeness
measurement per site and REFUSES when there is no record. It does not decide what ratio is
good enough: no threshold is set here, because a cut-off nobody stated is exactly the kind
of invented parameter that turns a judgement into a silent one. `check` reports the number
and says whether a record exists; what to do about a low ratio is the caller's -- and,
where it matters, Andy's -- decision.

The measurement itself comes from ``analysis/ept_coverage_check.py``, which asks the SOURCE
rather than the file: the EPT hierarchy JSON carries a point count per node, so the points
available over a bbox can be summed without downloading any, and compared with what our
file holds. That script writes the record; this module reads it and gates on it.
"""
from __future__ import annotations

import json
import os

__all__ = ["RECORD", "record_path", "read", "write", "check", "CompletenessUnknown"]

#: One JSON file per tile directory, so the record travels with the products it describes.
RECORD = "data_completeness.json"


class CompletenessUnknown(LookupError):
    """No completeness record for this site. UNKNOWN IS NOT A PASS."""


def record_path(tile_dir):
    return os.path.join(str(tile_dir), RECORD)


def read(tile_dir):
    """The record, or None if there is none. Callers that must not proceed use `check`."""
    p = record_path(tile_dir)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def write(tile_dir, *, epoch, cloud, points_in_file, points_available,
          measured_by, note=""):
    """Record one epoch's completeness for a tile.

    `points_available` is the AREA-WEIGHTED share of the source's node set inside the
    bounding box. Nodes clip the box, so a raw node-count sum overstates what a complete
    fetch would hold -- that overstatement once read two complete sites as short.

    The ratio is stored alongside its two inputs, never alone: a bare 0.67 cannot be
    rechecked, and which of the two numbers moved is the whole diagnosis.
    """
    if points_available <= 0:
        raise ValueError(
            f"points_available={points_available!r} for {epoch} at {tile_dir}: a ratio "
            f"against zero is not a measurement. Record the failure instead of a number.")
    rec = read(tile_dir) or {"tile": os.path.basename(str(tile_dir).rstrip("/")),
                             "epochs": {}}
    rec["epochs"][epoch] = {
        "cloud": cloud,
        "points_in_file": int(points_in_file),
        "points_available_in_bbox": float(points_available),
        "ratio": float(points_in_file) / float(points_available),
        "measured_by": measured_by,
        "note": note,
    }
    os.makedirs(str(tile_dir), exist_ok=True)
    with open(record_path(tile_dir), "w") as f:
        json.dump(rec, f, indent=2, sort_keys=True)
    return rec


def check(tile_dir, *, epochs=("gen1", "gen2"), require=True):
    """The completeness of each epoch, and whether it is known at all.

    Returns ``{epoch: {...record..., "known": bool}}``. With ``require=True`` (the default)
    a MISSING record raises: a site whose completeness has never been measured must not
    pass silently, because that is indistinguishable from one measured and found whole.

    NO THRESHOLD IS APPLIED. The ratio is reported; whether 0.67 is acceptable for a
    particular question is a judgement, and this module does not make it. Recording the
    number is what lets someone else make it.
    """
    rec = read(tile_dir)
    out = {}
    missing = []
    for e in epochs:
        r = (rec or {}).get("epochs", {}).get(e)
        if r is None:
            missing.append(e)
            out[e] = {"known": False}
        else:
            out[e] = dict(r, known=True)
    if require and missing:
        raise CompletenessUnknown(
            f"no completeness record for {', '.join(missing)} at {tile_dir}. UNKNOWN IS "
            f"NOT A PASS: a truncated fetch does not look like an error from inside the "
            f"tile -- whitewater's gen2 averaged 11.39 returns/m2 over a tile that was "
            f"15.45 west and 5.52 east of a seam. Run "
            f"analysis/ept_coverage_check.py --write to measure it, or pass require=False "
            f"to proceed while stating that you know it is unmeasured.")
    return out


def summary_line(tile_dir, epoch, c):
    """One line for a report. States the ratio AND its two inputs, never the ratio alone."""
    r = c.get(epoch, {})
    if not r.get("known"):
        return f"  {epoch}: completeness UNKNOWN (no record; unknown is not a pass)"
    return (f"  {epoch}: {r['points_in_file']:,} of an estimated "
            f"{r['points_available_in_bbox']:,.0f} available in the bbox "
            f"= {r['ratio']:.3f}  [{r['measured_by']}]")
