"""The derived-product chain, declared: what each step needs, makes, and runs.

The optional half of this pipeline is a dependency graph that lived only in whoever
remembered it. Getting the order wrong is not loud: running the convexity producer after
the curvature one silently drops three columns from ``ridgecrest_pixels.npz``, and a q2 fit
read from a ``beam_offset_table`` older than ``corrections.json`` carries registration terms
that are no longer in force. Both happened.

So the graph is data here, not lore. Each :class:`Step` declares what it REQUIRES and what
it PRODUCES, both as filenames inside a tile directory, and the order is derived from those
declarations rather than written down. That buys three things:

* ``plan``   -- the commands to run, in an order that satisfies the dependencies
* ``state``  -- per step: MISSING, STALE (an output older than one of its inputs), or OK
* a test that the graph is well formed, which fails if a step is added with a requirement
  nothing produces

This module does NOT run anything. It prints what to run and what is out of date, because
several steps take tens of minutes and read multi-gigabyte clouds, and because deciding to
regenerate an adopted product is a judgement, not a scheduling detail.

    lidar-diff-workflow --tile data/derived/elba_fulldensity --check
    lidar-diff-workflow --tile data/derived/elbaext --plan --gen1 <laz> --gen2 <laz>
"""
from __future__ import annotations

import argparse
import datetime as _dt
import os
from dataclasses import dataclass, field

PY = "env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python"

#: Products of ``pipeline.difference_dem`` (via ``analysis/ridgelines/run_elba_dod.py`` or
#: ``scripts/run_all_sites.py``). Everything below builds on these; the workflow treats them
#: as given rather than pretending to schedule the DoD itself, whose vertical frame is a
#: per-region decision (see ``analysis/slope_bias/elbaext_geoid_regrid.py``).
BASE_INPUTS = ("corrections.json", "z_after.npy", "dod.npy", "lod.npy")


@dataclass(frozen=True)
class Step:
    name: str
    produces: tuple[str, ...]
    requires: tuple[str, ...]
    command: str                      # {tile}, {gen1}, {gen2} are substituted
    optional: bool = False
    note: str = ""
    needs: tuple[str, ...] = field(default=())     # external args the command needs
    #: Files this step MODIFIES IN PLACE rather than creating -- it adds columns to a file
    #: another step owns. They are not `produces` (the other step owns the file) and must
    #: not be `requires` either, because running the step makes them newer than its own
    #: outputs and it would then report itself STALE forever. Declared so the augmentation
    #: is visible rather than implicit, and so the ordering it demands is documented.
    mutates: tuple[str, ...] = field(default=())


STEPS: tuple[Step, ...] = (
    Step("slope",
         produces=("slope.npy",),
         requires=("z_after.npy",),
         command=f"{PY} scripts/make_slope.py --tile {{tile_name}}",
         note="Surface slope from the gap-filled gen2 grid. It was NOT a base input: "
              "difference_dem does not write it, and its only producer was hardcoded to "
              "elba, which is why carlton and cook have every other base product but no "
              "slope."),
    Step("ridge_mask",
         produces=("ridge_mask.npy",),
         requires=("z_after.npy",),
         command=f"{PY} analysis/ridgelines/trace_ridgelines.py {{tile}} --out ridge_mask.npy",
         note="Scherler & Schwanghart divide network; the candidate ridge cells."),
    Step("convexity",
         produces=("floodplain_mask.npy", "crest_mask.npy", "kappa_L10.npy", "kappa_L20.npy",
                   "kappa_L30.npy", "ridgecrest_pixels.npz", "ridgecrest_pixels.csv"),
         requires=("z_after.npy", "slope.npy", "ridge_mask.npy"),
         command=f"{PY} analysis/ridgelines/convexity_dod_landcover.py --tile {{tile_name}} "
                 f"--dod {{dod}} --without penetration",
         note="floodplain_mask gates reference_cells, so this comes before any q2 work. "
              "Drop --without penetration once that layer exists to get the forest/open "
              "crest split."),
    Step("curvature",
         produces=("curv_xx.npy", "curv_yy.npy", "curv_laplacian.npy"),
         requires=("z_after.npy", "crest_mask.npy"),
         mutates=("ridgecrest_pixels.npz",),
         command=f"{PY} analysis/ridgelines/curvature_diffusion.py --tile {{tile_name}}",
         note="MUST follow convexity, which crest_mask.npy already enforces: this ADDS "
              "curv_xx/curv_yy/curv_laplacian columns to ridgecrest_pixels.npz in place, "
              "and the convexity producer rewrites that file from scratch, so the reverse "
              "order drops the three columns silently."),
    Step("pfs_cover",
         produces=("canopy_cover_pfs.npy", "forest_pfs.npy", "open_pfs.npy", "pai_pfs.npy"),
         requires=("z_after.npy",),
         command=f"{PY} analysis/forest_metrics_pfs.py {{tile}} {{gen2}}",
         needs=("gen2",),
         note="PyForestScan cover -- the cover measure computed identically on every tile."),
    Step("penetration",
         produces=("penetration.npy",),
         requires=("z_after.npy",),
         command=f"{PY} scripts/make_penetration.py --tile {{tile_name}} --after {{gen2}}",
         optional=True, needs=("gen2",),
         note="gen2 ground-return fraction. AUDIT_findings.md warns it is gen2-derived and "
              "should not drive gen1-internal conclusions; only strata_core and the "
              "forest/open crest split need it."),
    Step("gen1_angles",
         produces=("gen1_csf_angles.npz",),
         requires=("z_after.npy", "corrections.json"),
         command=f"{PY} analysis/ridgelines/gen1_save_angles_slope.py {{tile_name}} {{gen1}} "
                 f"--without penetration,core_forest,core_open,floodplain_mask",
         needs=("gen1",),
         note="Per-return gen1 CSF ground offsets with beam geometry. Drop the --without "
              "list for whichever strata the tile actually has."),
    Step("beam_table",
         produces=("beam_offset_table.parquet", "beam_offset_table.head.csv"),
         requires=("gen1_csf_angles.npz", "corrections.json", "canopy_cover_pfs.npy",
                   "curv_laplacian.npy"),
         command=f"{PY} analysis/ridgelines/beam_offset_table.py {{tile}} {{gen1}}",
         needs=("gen1",),
         note="Applies the four registration terms. If corrections.json is newer than this "
              "table, every q2 number downstream is on superseded registration."),
    Step("nearground",
         produces=("nearground_cells_sn.npz",),
         requires=("z_after.npy", "curv_laplacian.npy"),
         command=f"{PY} analysis/ridgelines/nearground_cells.py --tile {{tile}} "
                 f"--gen1 {{gen1}} --gen2 {{gen2}} --out nearground_cells_sn.npz",
         needs=("gen1", "gen2"),
         note="Slope-normal near-ground column, both epochs."),
    Step("nearground_split",
         produces=("nearground_gen2_class_split.npz",),
         requires=("z_after.npy", "nearground_cells_sn.npz"),
         command=f"{PY} analysis/ridgelines/nearground_class_split.py --tile {{tile}} "
                 f"--gen2 {{gen2}}",
         needs=("gen2",),
         note="gen2 class-2 near-ground histogram; the column q2 indexes into."),
    Step("q2_fit",
         produces=("q2_cover_fit.json",),
         requires=("nearground_cells_sn.npz", "z_after.npy", "canopy_cover_pfs.npy",
                   "beam_offset_table.parquet", "nearground_gen2_class_split.npz",
                   "floodplain_mask.npy", "curv_laplacian.npy"),
         command=f"{PY} analysis/ridgelines/q2_cover_fit.py --tile {{tile}}",
         note="PER SITE -- the relation depends on each pair's phenology and undergrowth, "
              "so there is no slope to carry between tiles."),
    Step("dod_cover",
         produces=("dod_cover_q2.npy", "dod_cover_q2.json", "dod_gen2_median.npy",
                   "gen2_q2_used.npy"),
         requires=("q2_cover_fit.json", "canopy_cover_pfs.npy", "z_after.npy"),
         command=f"{PY} analysis/ridgelines/dod_cover_corrected.py --tile {{tile}} "
                 f"--gen2 {{gen2}}",
         needs=("gen2",),
         note="The cover-corrected DoD. Reads its slope from q2_cover_fit.json."),
    Step("lod_cover",
         produces=("lod_cover_q2.npy",),
         requires=("dod_cover_q2.npy", "slope.npy", "curv_laplacian.npy"),
         command=f"{PY} analysis/ridgelines/lod_cover_q2.py --tile {{tile}}",
         note="LoD refitted on the corrected DoD."),
    Step("cover_calibration",
         produces=("cover_offset_calibration.json",),
         requires=("beam_offset_table.parquet", "canopy_cover_pfs.npy", "slope.npy",
                   "curv_laplacian.npy"),
         command=f"{PY} analysis/ridgelines/cover_offset_reference.py --tile {{tile}}",
         optional=True,
         note="offset-vs-cover on non-eroding ground; dod_cover_attribution.py reads it. "
              "Writes cover_offset_calibration_<tile>.json for tiles other than elba."),
    Step("canopy_struct",
         produces=("canopy_struct.npz",),
         requires=("z_after.npy",),
         command=f"{PY} analysis/ridgelines/canopy_struct.py --tile {{tile_name}} "
                 f"--after {{gen2}}",
         optional=True, needs=("gen2",),
         note="Per-cell canopy structure from the full unclassified gen2 cloud: veg_frac, "
              "understory/midstory fractions, canopy_height_p95, low_gap. Streams ~1.8e8 "
              "points."),
    Step("strata_core",
         produces=("core_forest.npy", "core_open.npy"),
         requires=("penetration.npy", "floodplain_mask.npy", "canopy_struct.npz",
                   "z_after.npy"),
         command=f"{PY} analysis/ridgelines/strata_core.py --tile {{tile_name}}",
         optional=True,
         note="BLOCKED: canopy_struct.npz has no producer in this repo."),
)


def mutation_order_ok(steps=STEPS):
    """Every mutated file must be produced by a step that runs EARLIER."""
    seq = [s.name for s in order(steps)]
    made = {f: s.name for s in steps for f in s.produces}
    bad = []
    for s in steps:
        for f in s.mutates:
            owner = made.get(f)
            if owner is None or seq.index(owner) >= seq.index(s.name):
                bad.append((s.name, f, owner))
    return bad


def order(steps=STEPS):
    """Steps in an order that satisfies every declared requirement.

    Raises if a requirement is neither a base input nor produced by some step, and if the
    remaining requirements cannot be satisfied (a cycle, or a missing producer).
    """
    made = {f: s.name for s in steps for f in s.produces}
    unknown = {r for s in steps for r in s.requires
               if r not in made and r not in BASE_INPUTS}
    if unknown:
        raise ValueError(
            f"these requirements are neither base inputs nor produced by any step: "
            f"{sorted(unknown)}. Either add the producing step or list the file in "
            f"BASE_INPUTS if it comes from difference_dem.")
    have, out, left = set(BASE_INPUTS), [], list(steps)
    while left:
        ready = [s for s in left if all(r in have for r in s.requires)]
        if not ready:
            stuck = {s.name: sorted(set(s.requires) - have) for s in left}
            raise ValueError(f"cannot order the graph; unmet requirements: {stuck}")
        ready.sort(key=lambda s: steps.index(s))
        s = ready[0]
        out.append(s); have.update(s.produces); left.remove(s)
    return out


def _mtime(path):
    return os.path.getmtime(path) if os.path.exists(path) else None


def state(tile_dir, steps=STEPS):
    """Per step: ('MISSING'|'STALE'|'OK', detail). STALE = an output older than an input."""
    res = {}
    for s in steps:
        outs = {f: _mtime(os.path.join(tile_dir, f)) for f in s.produces}
        ins = {f: _mtime(os.path.join(tile_dir, f)) for f in s.requires}
        absent_in = [f for f, m in ins.items() if m is None]
        absent_out = [f for f, m in outs.items() if m is None]
        if absent_out:
            res[s.name] = ("MISSING", f"no {', '.join(sorted(absent_out)[:3])}"
                           + (f" (+{len(absent_out)-3} more)" if len(absent_out) > 3 else "")
                           + (f"; inputs absent: {', '.join(sorted(absent_in))}" if absent_in else ""))
            continue
        oldest_out = min(outs.values())
        newer = sorted(f for f, m in ins.items() if m is not None and m > oldest_out)
        if newer:
            res[s.name] = ("STALE", f"input(s) newer than output: {', '.join(newer)}")
        else:
            res[s.name] = ("OK", "")
    return res


def plan(tile_dir, *, gen1=None, gen2=None, dod=None, steps=STEPS, include_optional=True):
    """The commands, in dependency order, with the tile substituted."""
    name = os.path.basename(str(tile_dir).rstrip("/"))
    supplied = {"gen1": gen1, "gen2": gen2}
    cmds = []
    for s in order(steps):
        if s.optional and not include_optional:
            continue
        missing = [n for n in s.needs if supplied.get(n) is None]
        c = (s.command.replace("{tile_name}", name).replace("{tile}", str(tile_dir))
             .replace("{gen1}", gen1 or "<--gen1 NOT GIVEN>")
             .replace("{gen2}", gen2 or "<--gen2 NOT GIVEN>")
             .replace("{dod}", dod or f"{tile_dir}/dod.npy"))
        cmds.append((s, c, missing))
    return cmds


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tile", required=True, help="tile directory under data/derived/")
    ap.add_argument("--gen1", help="gen1 CSF-cached cloud, for the steps that read it")
    ap.add_argument("--gen2", help="gen2 full-return cloud, for the steps that read it")
    ap.add_argument("--dod", help="DoD grid for the convexity step (default <tile>/dod.npy)")
    ap.add_argument("--check", action="store_true", help="report each step's state")
    ap.add_argument("--plan", action="store_true", help="print the commands in order")
    ap.add_argument("--skip-optional", action="store_true")
    a = ap.parse_args(argv)
    if not (a.check or a.plan):
        a.check = True

    if a.check:
        st = state(a.tile)
        base = [f for f in BASE_INPUTS if not os.path.exists(os.path.join(a.tile, f))]
        print(f"{a.tile}")
        if base:
            print(f"  BASE INPUTS ABSENT: {', '.join(base)} -- run difference_dem first")
        print(f"  {'step':<20} {'state':<8} detail")
        for s in order():
            k, d = st[s.name]
            tag = s.name + (" *" if s.optional else "")
            print(f"  {tag:<20} {k:<8} {d}")
        print("  * = optional. STALE means an output is older than one of its own inputs.")

    if a.plan:
        print(f"\n# {a.tile} -- in dependency order")
        for s, c, missing in plan(a.tile, gen1=a.gen1, gen2=a.gen2, dod=a.dod,
                                  include_optional=not a.skip_optional):
            print(f"\n# {s.name}{' (optional)' if s.optional else ''}: {s.note}")
            if missing:
                print(f"#   NEEDS --{' --'.join(missing)}")
            print(c)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
