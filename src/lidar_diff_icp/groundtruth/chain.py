"""Propagating a tie from the line that covers a checkpoint to the lines that cover the site.

The problem
-----------
Surveyed checkpoints and study sites are not on the same flight lines. At Elba the gen1
lines run north-south with a measured swath half-width of ~710-730 m, and every nearby
3DEP checkpoint is displaced 3.1-15.6 km **east-west**. No amount of extending the study
area moves a line sideways, so the four lines over Elba (135-138) cover no checkpoint and
never can (``analysis/ELBAEXT2_SCOPE.md`` §2).

But adjacent lines **overlap** -- 1.31-1.48 km² per pair at ~960 m spacing -- and
:func:`lidar_diff_icp.coreg.coregister_swaths` measures the vertical offset across an
overlap. So a tie measured on a line that does cover a checkpoint propagates, link by
link, to a line that covers the site.

Search order: cheapest first
----------------------------
1. **Along-swath, zero links.** If a target line itself covers the checkpoint, the tie
   applies directly and no cross-swath transfer happens at all. :func:`plan_path` tests
   this first, every time -- at Elba the answer is no, but the module is meant to be
   reused where the answer is yes and the chain is then free.
2. **Otherwise the shortest chain**, minimising **link count** -- not distance, not tile
   count -- because each link contributes its own alignment error and a chain has no
   internal redundancy to absorb it (see :func:`lidar_diff_icp.coreg.align_swaths`:
   misclosure is identically zero on a tree and carries no information).

Each link is solved on **its overlap only**: the pair's shared cells, not the tiles that
contain them. Only the tiles named by the chosen path are ever read.

Error
-----
Along a chain, per-link uncertainties add in quadrature and there is nothing to check
them against. The only real check is a **second, independent path** to the same target --
at Elba, the western chain to line 128 versus the eastern chain to line 144. Their
disagreement is the error bar, and :func:`compare_paths` reports it as such.

Sign convention
---------------
``coregister_swaths(pc, swath_ref=a, swath_src=b)`` returns the shift that moves ``b``
onto ``a``. So walking a path outward from a reference line ``r`` (whose offset is 0 by
definition of the frame), the offset of the next line is ``offset[b] = offset[a] + dz``.
``ChainSolution.dz_total_m`` is what must be ADDED to the far line's z to bring it into
the reference line's frame.
"""
from __future__ import annotations

import json
import os
from collections import deque
from dataclasses import dataclass, field

import numpy as np

from .provenance import Param

#: Classes excluded from the terrain proxy, exactly as
#: :func:`lidar_diff_icp.coreg.coregister_swaths` does it: high vegetation, building,
#: water. Everything else -- including unclassified (0), ground (2), model key point (8)
#: and overlap (12) -- is kept, so no CSF run is needed for the chain.
DEFAULT_EXCLUDE = (5, 6, 9)

#: Grid resolution for the Nuth & Kaeaeb link solve and the overlap-area measurement.
#: The default of :func:`lidar_diff_icp.coreg.coregister_swaths` and
#: :func:`~lidar_diff_icp.coreg.align_swaths`.
DEFAULT_RES = 2.0


# --------------------------------------------------------------------- tile inventory

@dataclass
class TileLines:
    path: str
    lines: dict          # point_source_id -> {"n": int, "bbox": [x0, y0, x1, y1]}
    n_points: int
    n_terrain: int


@dataclass
class SwathInventory:
    """Which flight lines are in which tiles, and where their points are.

    Built by one pass over each tile. ``cache_dir`` stores the terrain points per tile as
    float32 offsets so a re-run costs a load instead of a LAZ decompress; the inventory
    itself is JSON next to them.
    """

    tiles: dict = field(default_factory=dict)         # path -> TileLines
    cache_dir: str | None = None
    exclude: tuple = DEFAULT_EXCLUDE
    _points: dict = field(default_factory=dict, repr=False)

    @property
    def lines(self) -> list:
        out = set()
        for t in self.tiles.values():
            out |= set(t.lines)
        return sorted(out)

    def tiles_for(self, line: int) -> list:
        return [p for p, t in self.tiles.items() if line in t.lines]

    def tiles_for_pair(self, a: int, b: int) -> list:
        """Tiles holding both lines -- the only tiles a link's overlap can live in."""
        return [p for p, t in self.tiles.items() if a in t.lines and b in t.lines]

    def bbox(self, line: int):
        bb = [t.lines[line]["bbox"] for t in self.tiles.values() if line in t.lines]
        if not bb:
            raise KeyError(line)
        a = np.asarray(bb, float)
        return (a[:, 0].min(), a[:, 1].min(), a[:, 2].max(), a[:, 3].max())

    def points(self, path: str):
        """Terrain points of one tile as ``(x, y, z, psid)``, cached in memory."""
        if path in self._points:
            return self._points[path]
        arr = _load_tile_points(path, self.cache_dir, self.exclude)
        self._points[path] = arr
        return arr

    def release(self, path: str | None = None):
        """Drop cached point arrays (all of them, or one tile's)."""
        if path is None:
            self._points.clear()
        else:
            self._points.pop(path, None)


def _cache_name(path, cache_dir):
    if not cache_dir:
        return None
    st = os.stat(path)
    return os.path.join(cache_dir,
                        f"{os.path.basename(path)}.{st.st_size}.{int(st.st_mtime)}.terrain.npz")


def _load_tile_points(path, cache_dir, exclude):
    npz = _cache_name(path, cache_dir)
    if npz and os.path.exists(npz):
        d = np.load(npz)
        o = d["origin"]
        return (d["x"].astype(np.float64) + o[0], d["y"].astype(np.float64) + o[1],
                d["z"].astype(np.float64) + o[2], d["psid"])
    import laspy
    f = laspy.read(str(path))
    keep = ~np.isin(np.asarray(f.classification), exclude)
    x = np.asarray(f.x)[keep]; y = np.asarray(f.y)[keep]; z = np.asarray(f.z)[keep]
    ps = np.asarray(f.point_source_id)[keep].astype(np.int32)
    if npz:
        os.makedirs(cache_dir, exist_ok=True)
        o = np.array([np.floor(x.min()), np.floor(y.min()), np.floor(z.min())])
        np.savez(npz, origin=o, x=(x - o[0]).astype(np.float32),
                 y=(y - o[1]).astype(np.float32), z=(z - o[2]).astype(np.float32), psid=ps)
    return x, y, z, ps


def build_inventory(tile_paths, *, cache_dir=None, exclude=DEFAULT_EXCLUDE,
                    inventory_json=None) -> SwathInventory:
    """One pass per tile: which lines it holds, their bboxes and counts.

    ``inventory_json``, if given, is read when it is current for every tile (size and
    mtime match) and written otherwise, so repeated runs do not re-decompress.
    """
    tile_paths = [str(p) for p in tile_paths]
    inv = SwathInventory(cache_dir=cache_dir, exclude=tuple(exclude))
    stamp = {p: [os.path.getsize(p), int(os.path.getmtime(p))] for p in tile_paths}
    if inventory_json and os.path.exists(inventory_json):
        try:
            j = json.load(open(inventory_json))
            if j.get("stamp") == {k: v for k, v in stamp.items()} and \
                    j.get("exclude") == list(exclude):
                for p, t in j["tiles"].items():
                    inv.tiles[p] = TileLines(p, {int(k): v for k, v in t["lines"].items()},
                                             t["n_points"], t["n_terrain"])
                return inv
        except Exception:
            pass
    import laspy
    for p in tile_paths:
        f = laspy.read(p)
        cl = np.asarray(f.classification)
        keep = ~np.isin(cl, exclude)
        x = np.asarray(f.x)[keep]; y = np.asarray(f.y)[keep]
        ps = np.asarray(f.point_source_id)[keep]
        lines = {}
        for s in np.unique(ps):
            m = ps == s
            lines[int(s)] = {"n": int(m.sum()),
                             "bbox": [float(x[m].min()), float(y[m].min()),
                                      float(x[m].max()), float(y[m].max())]}
        inv.tiles[p] = TileLines(p, lines, int(len(cl)), int(keep.sum()))
        del f, cl, keep, x, y, ps
    if inventory_json:
        os.makedirs(os.path.dirname(inventory_json) or ".", exist_ok=True)
        json.dump({"stamp": stamp, "exclude": list(exclude),
                   "tiles": {p: {"lines": {str(k): v for k, v in t.lines.items()},
                                 "n_points": t.n_points, "n_terrain": t.n_terrain}
                             for p, t in inv.tiles.items()}},
                  open(inventory_json, "w"), indent=1)
    return inv


# ------------------------------------------------------------------------- the graph

def _cells(x, y, res):
    """Cell ids on a GLOBAL grid anchored at (0, 0), so ids from different tiles agree
    and a cell counted twice in a tile overlap collapses to one."""
    return (np.floor(y / res).astype(np.int64) << 32) + np.floor(x / res).astype(np.int64)


def overlap_graph(inv: SwathInventory, *, res=DEFAULT_RES, lines=None):
    """Overlap area (km²) for every pair of lines that share cells.

    Only pairs that co-occur in some tile are tested (a pair in no common tile has no
    measurable overlap here), and the cell ids are global, so a pair whose overlap
    straddles two tiles is measured once, not twice.
    """
    lines = set(inv.lines) if lines is None else set(lines)
    pairs = {}
    for path, t in inv.tiles.items():
        present = sorted(set(t.lines) & lines)
        if len(present) < 2:
            continue
        x, y, _, ps = inv.points(path)
        cells = {s: np.unique(_cells(x[ps == s], y[ps == s], res)) for s in present}
        for i, a in enumerate(present):
            for b in present[i + 1:]:
                inter = np.intersect1d(cells[a], cells[b], assume_unique=True)
                if inter.size:
                    key = (a, b)
                    pairs[key] = (np.union1d(pairs[key], inter) if key in pairs else inter)
        del cells
    return {k: {"cells": int(v.size), "area_km2": v.size * res * res / 1e6}
            for k, v in pairs.items()}


def covering_lines(inv: SwathInventory, easting, northing, radius, *, res=DEFAULT_RES):
    """Lines with terrain returns within ``radius`` of a point, and how many.

    This is the direct test -- returns on the ground near the mark -- rather than a
    cross-track distance against a fitted nadir track, so it needs no track model and
    cannot be fooled by a line that ends before it gets there.
    """
    out = {}
    for path, t in inv.tiles.items():
        x0, y0, x1, y1 = (min(v["bbox"][0] for v in t.lines.values()),
                          min(v["bbox"][1] for v in t.lines.values()),
                          max(v["bbox"][2] for v in t.lines.values()),
                          max(v["bbox"][3] for v in t.lines.values()))
        if not (x0 - radius <= easting <= x1 + radius and y0 - radius <= northing <= y1 + radius):
            continue
        x, y, _, ps = inv.points(path)
        near = (np.abs(x - easting) <= radius) & (np.abs(y - northing) <= radius)
        near &= np.hypot(x - easting, y - northing) <= radius
        for s, c in zip(*np.unique(ps[near], return_counts=True)):
            out[int(s)] = out.get(int(s), 0) + int(c)
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


# ---------------------------------------------------------------------------- paths

@dataclass
class ChainPath:
    """A route from a checkpoint's line to a target line, and why it was chosen."""

    nodes: list                    # [source_line, ..., target_line]
    edges: list                    # [(a, b), ...]
    along_swath: bool              # True when nodes has length 1: zero links
    n_links: int
    overlap_km2: list              # per edge, from the graph
    tiles: list                    # tiles that must be read to solve this path
    note: str = ""

    def __str__(self):             # pragma: no cover - display only
        return " - ".join(str(n) for n in self.nodes) + f"  ({self.n_links} link(s))"


def _bfs_paths(adj, sources, targets, *, max_paths=8):
    """All fewest-link paths from any source to any target, breadth-first.

    Level-by-level, so the first level at which a target is reached yields every path of
    that length and the search stops there. Link count is the only cost: distance and
    tile count do not enter, because what accumulates along a chain is alignment error
    per link.
    """
    sources = list(dict.fromkeys(int(s) for s in sources))
    targets = {int(t) for t in targets}
    level = [[s] for s in sources]
    depth = {s: 0 for s in sources}
    found = []
    d = 0
    while level and not found and d <= len(adj) + 1:
        nxt = []
        for p in level:
            for nb in sorted(adj.get(p[-1], ())):
                if nb in p or depth.get(nb, 1 << 30) < d + 1:
                    continue
                depth[nb] = d + 1
                q = p + [nb]
                (found if nb in targets else nxt).append(q)
        level = nxt
        d += 1
    return found[:max_paths]


def plan_path(graph, inv: SwathInventory, source_lines, target_lines, *, max_paths=8):
    """Plan routes from ``source_lines`` (cover the checkpoint) to ``target_lines``.

    **Step 1, always: the along-swath test.** If a source line is itself a target line,
    the tie needs no cross-swath transfer -- returned as a zero-link
    :class:`ChainPath` and nothing else is searched. **Step 2:** breadth-first over the
    overlap graph, so the first solutions found are the ones with the fewest links.

    Returns a list of :class:`ChainPath`, shortest first. Empty when the graph does not
    connect the two sets, which is a reportable result, not an error.
    """
    src = list(dict.fromkeys(int(s) for s in source_lines))
    tgt = list(dict.fromkeys(int(t) for t in target_lines))
    direct = [s for s in src if s in tgt]
    if direct:
        return [ChainPath([d], [], True, 0, [], inv.tiles_for(d),
                          note="along-swath: the checkpoint is under a target line, "
                               "so no cross-swath transfer is needed")
                for d in direct]
    adj = {}
    for (a, b) in graph:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    out = []
    for nodes in _bfs_paths(adj, src, tgt, max_paths=max_paths):
        edges = list(zip(nodes[:-1], nodes[1:]))
        areas = [graph.get(e, graph.get(e[::-1], {})).get("area_km2", float("nan"))
                 for e in edges]
        tiles = []
        for a, b in edges:
            for p in inv.tiles_for_pair(a, b):
                if p not in tiles:
                    tiles.append(p)
        out.append(ChainPath(nodes, edges, False, len(edges), areas, tiles))
    return out


# ---------------------------------------------------------------------- link solving

@dataclass
class Link:
    a: int
    b: int
    dz_m: float          # add to b's z to bring it into a's frame
    dx_m: float
    dy_m: float
    dz_sigma_m: float
    n_cells: int
    overlap_km2: float
    nmad_before_m: float
    nmad_after_m: float   # the residual left in the overlap after the link is applied
    converged: bool
    tiles: list


@dataclass
class ChainSolution:
    path: ChainPath
    links: list
    dz_total_m: float          # add to the far line's z to reach the reference frame
    dx_total_m: float
    dy_total_m: float
    dz_sigma_m: float          # per-link sigmas in quadrature
    reference_line: int
    far_line: int
    params: list = field(default_factory=list)

    def table_rows(self):
        return [[f"{l.a}-{l.b}", f"{l.overlap_km2:.3f}", l.n_cells, f"{l.dz_m*1000:+.1f}",
                 f"{l.dz_sigma_m*1000:.1f}", f"{l.nmad_before_m*1000:.0f}",
                 f"{l.nmad_after_m*1000:.0f}", f"{l.dx_m:+.3f}", f"{l.dy_m:+.3f}"]
                for l in self.links]

    @staticmethod
    def table_columns() -> dict:
        return {
            "link": "the flight-line pair, ref-src",
            "area_km2": "measured overlap area of the pair on the solve grid, km^2",
            "cells": "grid cells used in the Nuth & Kaeaeb fit, count",
            "dz_mm": "vertical shift to add to src to align it onto ref, mm",
            "dz_sig_mm": "1-sigma of that shift (nmad_after / sqrt(n)), mm",
            "nmad0_mm": "robust scatter of the overlap difference before the link, mm",
            "nmad1_mm": "and after -- the residual the link could not remove, mm",
            "dx_m": "eastward shift to add to src, m",
            "dy_m": "northward shift to add to src, m",
        }


def _pair_cloud(inv: SwathInventory, a: int, b: int, *, res, pad):
    """Points of lines ``a`` and ``b`` **inside their overlap only**, across the tiles
    that hold both. ``pad`` metres of margin are kept because Nuth & Kaeaeb needs slope
    from the neighbouring cells of the overlap edge."""
    from ..io import PointCloud

    tiles = inv.tiles_for_pair(a, b)
    if not tiles:
        raise ValueError(f"lines {a} and {b} share no tile in this inventory")
    xs, ys, zs, ps = [], [], [], []
    box = None
    for path in tiles:
        x, y, z, s = inv.points(path)
        ma = s == a; mb = s == b
        if not ma.any() or not mb.any():
            continue
        bb = (max(x[ma].min(), x[mb].min()) - pad, max(y[ma].min(), y[mb].min()) - pad,
              min(x[ma].max(), x[mb].max()) + pad, min(y[ma].max(), y[mb].max()) + pad)
        if bb[2] <= bb[0] or bb[3] <= bb[1]:
            continue
        sel = ((s == a) | (s == b)) & (x >= bb[0]) & (x <= bb[2]) & (y >= bb[1]) & (y <= bb[3])
        xs.append(x[sel]); ys.append(y[sel]); zs.append(z[sel]); ps.append(s[sel])
        box = bb if box is None else (min(box[0], bb[0]), min(box[1], bb[1]),
                                      max(box[2], bb[2]), max(box[3], bb[3]))
    if not xs:
        raise ValueError(f"lines {a} and {b} do not overlap in any shared tile")
    x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
    s = np.concatenate(ps)
    pc = PointCloud(x=x, y=y, z=z, point_source_id=s,
                    classification=np.full(x.size, 2, np.uint8),   # pre-filtered terrain
                    gps_time=np.zeros(x.size), scan_angle=np.zeros(x.size, np.int8),
                    crs="EPSG:26915")
    return pc, tiles, box


def solve_link(inv: SwathInventory, a: int, b: int, *, res=DEFAULT_RES, graph=None,
               pad=None) -> Link:
    """Solve one link with the repo's own :func:`lidar_diff_icp.coreg.coregister_swaths`.

    Only the pair's overlap region is loaded and gridded. ``pad`` defaults to ``2*res``
    -- two cells of margin, the stencil Nuth & Kaeaeb's slope/aspect gradients need at
    the overlap edge.
    """
    from .. import coreg

    pad = 2.0 * res if pad is None else float(pad)
    pc, tiles, _ = _pair_cloud(inv, a, b, res=res, pad=pad)
    c = coreg.coregister_swaths(pc, a, b, res, exclude=())     # already terrain-only
    area = float("nan")
    if graph is not None:
        g = graph.get((a, b), graph.get((b, a)))
        if g:
            area = g["area_km2"]
    return Link(a=a, b=b, dz_m=float(c.dz), dx_m=float(c.dx), dy_m=float(c.dy),
                dz_sigma_m=float(c.dz_sigma), n_cells=int(c.n), overlap_km2=area,
                nmad_before_m=float(c.nmad_before), nmad_after_m=float(c.nmad_after),
                converged=bool(c.converged), tiles=tiles)


def solve_chain(inv: SwathInventory, path: ChainPath, *, res=DEFAULT_RES, graph=None,
                reference="last") -> ChainSolution:
    """Solve every link of ``path`` and accumulate the shift into one frame.

    ``reference`` names which END of the path is the frame everything is brought into:
    ``"last"`` (the default) makes ``path.nodes[-1]`` the reference, which is what you
    want when :func:`plan_path` was called with the STUDY lines as targets --
    ``dz_total_m`` is then what to ADD to the checkpoint's line to reach the study
    frame. ``"first"`` reverses it.
    """
    if reference not in ("first", "last"):
        raise ValueError("reference must be 'first' or 'last'")
    nodes = list(path.nodes) if reference == "first" else list(path.nodes)[::-1]
    links = []
    dz = dx = dy = 0.0
    var = 0.0
    for a, b in zip(nodes[:-1], nodes[1:]):
        L = solve_link(inv, a, b, res=res, graph=graph)
        links.append(L)
        dz += L.dz_m; dx += L.dx_m; dy += L.dy_m
        var += L.dz_sigma_m ** 2
    params = [
        Param("link_res_m", res, "repo",
              "coreg.coregister_swaths / align_swaths default grid resolution"),
        Param("link_exclude_classes", tuple(inv.exclude), "repo",
              "coreg.coregister_swaths terrain proxy ~isin(classification,(5,6,9)) -- the "
              "VENDOR classification, so the chain needs no CSF run"),
        Param("path", list(nodes), "repo",
              "fewest-link route found by breadth-first search over the measured overlap "
              "graph; link count is the cost because each link adds error"),
    ]
    return ChainSolution(path=path, links=links, dz_total_m=dz, dx_total_m=dx,
                         dy_total_m=dy, dz_sigma_m=float(np.sqrt(var)),
                         reference_line=nodes[0], far_line=nodes[-1], params=params)


def compare_paths(solutions) -> dict:
    """Disagreement between independent chains to the same frame -- the real error bar.

    A chain's misclosure is identically zero, so its internal residuals cannot detect an
    accumulated error. Two independent routes can. Returns the per-solution totals, their
    spread, and the quadrature-summed formal sigma for comparison -- if the spread is much
    larger than the formal sigma, the formal sigma is not to be believed.
    """
    tot = {f"{s.reference_line}->{s.far_line} ({s.path.n_links} links)": s.dz_total_m * 1000.0
           for s in solutions}
    vals = list(tot.values())
    formal = [s.dz_sigma_m * 1000.0 for s in solutions]
    return {"dz_total_mm": tot,
            "spread_mm": (max(vals) - min(vals)) if len(vals) > 1 else float("nan"),
            "formal_sigma_mm": formal,
            "n_paths": len(solutions)}
