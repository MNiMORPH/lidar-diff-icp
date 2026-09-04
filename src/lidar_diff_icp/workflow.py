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
import hashlib
import json
import os
import re
import subprocess
from dataclasses import dataclass, field

PY = "env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python"

#: PyForestScan is NOT installed in the venv that runs everything else, so `pfs_cover` --
#: and only `pfs_cover` -- needs the conda environment. Until 2026-09-02 this step was
#: printed with `PY` like every other, which meant `--plan` emitted a command that dies on
#: `ModuleNotFoundError: No module named 'pyforestscan'` at every tile: a plan that cannot
#: be run is worse than no plan, because it looks authoritative.
#:
#: PROJ_DATA/GDAL_DATA are SET here rather than unset (the inverse of `PY`): the conda
#: env's own grids are needed, because the base proj.db it would otherwise find is stale.
#: Override with LIDAR_DIFF_PFS_PYTHON if the environment lives elsewhere.
_PFS_ENV = os.environ.get("LIDAR_DIFF_PFS_ENV", os.path.expanduser("~/anaconda3/envs/lidar-icp"))
PY_PFS = os.environ.get(
    "LIDAR_DIFF_PFS_PYTHON",
    f"PROJ_DATA={_PFS_ENV}/share/proj GDAL_DATA={_PFS_ENV}/share/gdal {_PFS_ENV}/bin/python")

#: Products of ``pipeline.difference_dem`` (via ``analysis/ridgelines/run_elba_dod.py`` or
#: ``scripts/run_all_sites.py``). Everything below builds on these; the workflow treats them
#: as given rather than pretending to schedule the DoD itself, whose vertical frame is a
#: per-region decision (see ``analysis/slope_bias/elbaext_geoid_regrid.py``).
BASE_INPUTS = ("corrections.json", "z_after.npy", "dod.npy", "lod.npy")


#: What produces the base inputs, for the code-vs-product check. difference_dem writes
#: them; run_all_sites.py is the driver. A change to either invalidates every tile's DoD and
#: LoD at once -- the failure this check exists for, since nothing inside a tile directory
#: shows it.
BASE_CODE = ("src/lidar_diff_icp/pipeline.py", "scripts/run_all_sites.py",
             "src/lidar_diff_icp/groundq.py")

#: The ground-q calibration is an INPUT to the base products, not a tile product: since
#: 2026-09-04 difference_dem takes each cell's percentile from it, so re-calibrating changes
#: every tile's z_after and DoD. It lives outside the tile directory, which is exactly why it
#: needs listing here -- nothing inside a tile reveals that its ground came from an older
#: curve. Produced by analysis/calibrate_ground_q.py.
BASE_GLOBAL_INPUTS = ("data/derived/ground_q_vs_class2sd_gen2_2021_control.npz",)


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
    #: Extra source files whose change invalidates this step's outputs, beyond the script
    #: named in `command` (which is parsed out automatically). Library modules go here.
    code: tuple[str, ...] = field(default=())


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
         command=f"{PY_PFS} analysis/forest_metrics_pfs.py {{tile}} {{gen2}}",
         optional=True, needs=("gen2",),
         note="PyForestScan cover. OPT-IN (Andy, 2026-09-02): run it only when a COVER "
              "CORRECTION is actually wanted at this site -- do not build it as a matter "
              "of course. Nothing else needs it: slope, ridge_mask, convexity and curvature "
              "never read it, and gen1_angles and beam_table only CARRY it as a column, "
              "omitting that column cleanly when it is absent. Its real consumers are "
              "q2_fit / dod_cover / cover_calibration, where it is definitional. Building "
              "it anyway is not free -- a large tile needs an untwine COPC first, ~14 GB of "
              "scratch for a 1 GB cloud."),
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
         command=f"{PY} analysis/ridgelines/gen1_save_angles_slope.py {{tile_name}} {{gen1}}",
         needs=("gen1",),
         note="Per-return gen1 CSF ground offsets with beam geometry. No flags needed: a "
              "stratum the tile lacks is simply omitted from the archive and said to be "
              "omitted. --without is for EXCLUDING a layer that is present."),
    Step("beam_table",
         produces=("beam_offset_table.parquet", "beam_offset_table.head.csv"),
         requires=("gen1_csf_angles.npz", "corrections.json", "curv_laplacian.npy"),
         command=f"{PY} analysis/ridgelines/beam_offset_table.py {{tile}} {{gen1}}",
         needs=("gen1",),
         note="Applies the four registration terms. Does NOT require canopy cover: it only "
              "carries it as a column, omitted when the layer is absent. If corrections.json "
              "is newer than this table, every q2 number downstream is on superseded "
              "registration."),
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


#: Where a tile's content hashes are remembered between checks.
MANIFEST = ".workflow_hashes.json"


def _mtime(path):
    return os.path.getmtime(path) if os.path.exists(path) else None


def _sha256(path, _chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(_chunk), b""):
            h.update(blk)
    return h.hexdigest()


def _load_manifest(tile_dir):
    p = os.path.join(tile_dir, MANIFEST)
    try:
        return json.load(open(p))
    except Exception:
        return {}


def _save_manifest(tile_dir, man):
    try:
        json.dump(man, open(os.path.join(tile_dir, MANIFEST), "w"), indent=1, sort_keys=True)
    except OSError:
        pass


def effective_mtime(tile_dir, name, man, *, update=True):
    """The time this file's CONTENT last changed, not the time it was last written.

    A producer that re-runs and writes a byte-identical file bumps its mtime, and a
    plain mtime comparison then cascades false staleness through everything downstream --
    at Carlton that meant a six-minute PyForestScan rebuild that would have changed
    nothing. So: mtime is the cheap trigger, the hash is the arbiter. When a file's
    (size, mtime) has moved but its sha256 has not, the recorded content-time is kept.

    LIMIT, stated because it is easy to over-read: this cannot tell that two files written
    seconds apart in the SAME run belong together. It forgives an identical rewrite; it
    does not reconstruct provenance the producers never recorded.
    """
    p = os.path.join(tile_dir, name)
    if not os.path.exists(p):
        return None
    st = os.stat(p)
    rec = man.get(name)
    if rec and rec.get("size") == st.st_size and rec.get("mtime") == st.st_mtime:
        return rec.get("content_time", st.st_mtime)      # unchanged since last seen
    digest = _sha256(p)
    if rec and rec.get("sha256") == digest:              # rewritten, byte-identical
        content_time = rec.get("content_time", st.st_mtime)
    else:
        content_time = st.st_mtime                       # genuinely new content
    if update:
        man[name] = {"size": st.st_size, "mtime": st.st_mtime,
                     "sha256": digest, "content_time": content_time}
    return content_time


_SCRIPT_RE = re.compile(r"(?:^|\s)((?:analysis|scripts|src)/[\w/]+\.py)")


def script_of(step):
    """Source files whose change should invalidate this step's outputs."""
    found = tuple(_SCRIPT_RE.findall(step.command))
    return tuple(dict.fromkeys(found + tuple(step.code)))


def _git_commit_time(path):
    try:
        out = subprocess.run(["git", "log", "-1", "--format=%ct", "--", path],
                             capture_output=True, text=True, timeout=10)
        return float(out.stdout.strip()) if out.stdout.strip() else None
    except Exception:
        return None


def _is_dirty(path):
    try:
        out = subprocess.run(["git", "status", "--porcelain", "--", path],
                             capture_output=True, text=True, timeout=10)
        return bool(out.stdout.strip())
    except Exception:
        return False


def code_time(path):
    """When this source last CHANGED, in a way that survives a fresh checkout.

    The last commit that touched it, because a clone or a branch switch rewrites every
    working-tree mtime and would otherwise mark the whole world stale. If the file is
    modified against HEAD, its working-tree mtime is used instead -- an uncommitted edit is
    a real change.
    """
    if not os.path.exists(path):
        return None
    if _is_dirty(path):
        return os.path.getmtime(path)
    return _git_commit_time(path) or os.path.getmtime(path)


def code_state(tile_dir, steps=STEPS):
    """Per step: source files that changed AFTER the step's outputs were written.

    The blind spot this closes: comparing product mtimes only against each other cannot see
    that the CODE moved on. A tile reports every step OK while its products are several
    commits behind, and a code change invalidates every tile at once, so nothing inside one
    tile directory reveals it.
    """
    out = {}
    for s in steps:
        outs = [_mtime(os.path.join(tile_dir, f)) for f in s.produces]
        if any(m is None for m in outs):
            continue                                   # MISSING is state()'s to report
        oldest = min(outs)
        newer = [(c, code_time(c)) for c in script_of(s)]
        hits = sorted(c for c, t in newer if t is not None and t > oldest)
        if hits:
            out[s.name] = hits
    return out


def base_code_state(tile_dir):
    """The same check for the base inputs, which no Step produces."""
    outs = [_mtime(os.path.join(tile_dir, f)) for f in BASE_INPUTS]
    present = [m for m in outs if m is not None]
    if not present:
        return []
    oldest = min(present)
    stale = [c for c in BASE_CODE if (code_time(c) or 0) > oldest]
    # The calibration curve goes through the SAME content-time path as source, deliberately:
    # re-running the calibration and getting the same curve back must NOT invalidate every
    # tile in the project. Only a curve whose CONTENT changed does.
    stale += [c for c in BASE_GLOBAL_INPUTS if (code_time(c) or 0) > oldest]
    return sorted(stale)


def state(tile_dir, steps=STEPS):
    """Per step: ('MISSING'|'STALE'|'OK', detail). STALE = an output older than an input."""
    man = _load_manifest(tile_dir)
    res = {}
    for s in steps:
        # ASYMMETRIC, and it has to be. An OUTPUT is judged by when it was last PRODUCED
        # (its mtime): re-running a step against changed inputs clears its staleness even
        # if the bytes come out identical -- which is the normal case for a step the change
        # does not reach. An INPUT is judged by when its CONTENT last changed, so a
        # producer that rewrites a file identically does not cascade.
        #
        # Using content_time on both sides was the bug: gen1_csf_angles.npz was rebuilt
        # after corrections.json changed, came out byte-identical, kept its old
        # content_time, and so could never stop reporting stale however often it was run.
        outs = {f: _mtime(os.path.join(tile_dir, f)) for f in s.produces}
        ins = {f: effective_mtime(tile_dir, f, man) for f in s.requires}
        for f in s.produces:                       # keep the manifest current for consumers
            effective_mtime(tile_dir, f, man)
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
    _save_manifest(tile_dir, man)
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
        cs = code_state(a.tile)
        base = [f for f in BASE_INPUTS if not os.path.exists(os.path.join(a.tile, f))]
        print(f"{a.tile}")
        if base:
            print(f"  BASE INPUTS ABSENT: {', '.join(base)} -- run difference_dem first")
        print(f"  {'step':<20} {'state':<12} detail")
        for s in order():
            k, d = st[s.name]
            if s.name in cs:                       # code moved on, whatever the files say
                k = "CODE-STALE" if k == "OK" else k + "+CODE"
                d = (d + "; " if d else "") + "changed since: " + ", ".join(
                    os.path.basename(c) for c in cs[s.name])
            tag = s.name + (" *" if s.optional else "")
            print(f"  {tag:<20} {k:<12} {d}")
        bc = base_code_state(a.tile)
        if bc:
            print(f"  {'(base inputs)':<20} {'CODE-STALE':<12} "
                  f"dod/lod/z_after predate: {', '.join(os.path.basename(c) for c in bc)}")
        print("  * = optional. STALE = an output older than one of its own inputs.")
        print("  CODE-STALE = the SOURCE that produced it changed since. A code change "
              "invalidates every")
        print("  tile at once, so nothing inside one tile directory reveals it.")

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
