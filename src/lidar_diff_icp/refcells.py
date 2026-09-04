"""The stable reference-cell population used to calibrate gen1 against gen2.

Everything that estimates *how* gen1's returns map onto gen2's -- the near-ground rank
gap, the matching quantile, the cover curve -- is calibrated on cells that should not
have changed between epochs. Getting that population right is the whole experiment, so it
is defined ONCE, here, instead of being re-derived per script.

**The governing constraint: the stability test must not look at the vertical offset.**
The quantity being calibrated (a few mm in the open, up to ~130 mm under dense canopy) is
the same size as the LoD (median ~105 mm on these cells). So excluding cells where
``|DoD| > LoD`` -- the standard change detector -- would strip out precisely the
high-cover cells carrying the signal, and the calibration would then "prove" that
vegetation has no effect. That is circular, and it is available here only as an explicit
sensitivity (``gross_change_mm`` set low), never as the default.

Every default criterion is therefore either geometric or vegetation-structural:

* **divide cells, low concavity** (``ridge_mask``, ``|curv_laplacian| <= curv_max``) --
  hilltops shed rather than collect, so they are the geomorphic no-change population. This
  is THE population. A slope + TPI "low-gradient upland" proxy is not a substitute: it keeps
  ground that both receives and sheds, and it was removed from the tree on 2026-09-04.
* **valley floor cut by ELEVATION**, at the antimode of this tile's own elevation histogram
  -- not by the TPI floodplain mask, whose extent depends on the window width and which
  keeps flat terrace ground at valley level.
* **gentle slope** (``slope_max``) -- excludes cells where mass wasting is plausible and
  where intra-cell relief dominates the return column. On these tiles it does not move
  any answer; it trims the noise.
* **no building returns in either epoch** (``n_bldg``) -- structures come and go, and
  their returns are not ground.
* **not in a known retreat zone** (``blufftop_margin_mask``, where the tile has one).
* **no gross change** (``|DoD| > gross_change_mm``, default 500 mm) -- quarries, road
  construction, fills. This DOES look at the offset, but at 25x its scale and ~4x the
  LoD, so it cannot remove a vegetation-driven cell. It bites 0.65-0.67% at Elba.
* **not clear-cut between epochs** (``frac_gen1 - frac_gen2 > clearcut_drop``) -- canopy
  present in 2008 and gone in 2021 means the ground itself was worked over. The test is
  ONE-SIDED on purpose: gen1 is leaf-off November and gen2 leaf-on May, so a deciduous
  stand legitimately gains canopy fraction between epochs. Excluding on ``|delta frac|``
  would strip deciduous stands for their phenology and bias the cover strata.
"""
from __future__ import annotations

import os

import numpy as np

__all__ = ["reference_cells", "floodplain_by_elevation", "valley_top_from_histogram",
           "landscape_of", "valley_top_for_landscape", "VALLEY_TOP_M"]


def landscape_of(tile_dir, root="data/derived"):
    """Every tile that shares ground with this one, transitively. The LANDSCAPE.

    A valley top is a property of the landscape, not of the rectangle we happened to cut
    from it, so it must be computed on ground the tiles share. Computed per tile it is not
    even self-consistent: refcells' own record has elba at 228.9 m and elbaext at 237.1 m
    ON OVERLAPPING GROUND, a 31,242-cell disagreement about which cells are floodplain in
    the same valley.

    Membership is MEASURED from the tiles' own bounds -- overlap, then transitive closure --
    never declared. Returns a sorted list of tile directory names including this one.
    """
    import json
    B = {}
    for d in sorted(os.listdir(root)):
        p = os.path.join(root, d, "corrections.json")
        if os.path.exists(p):
            try:
                B[d] = json.load(open(p))["bounds"]
            except Exception:
                pass
    name = os.path.basename(str(tile_dir).rstrip("/"))
    if name not in B:
        return [name]

    def hit(a, b):
        return (min(B[a][2], B[b][2]) > max(B[a][0], B[b][0])
                and min(B[a][3], B[b][3]) > max(B[a][1], B[b][1]))

    group, frontier = {name}, [name]
    while frontier:
        a = frontier.pop()
        for b in B:
            if b not in group and hit(a, b):
                group.add(b); frontier.append(b)
    return sorted(group)


def valley_top_for_landscape(tile_dir, root="data/derived", bins=60):
    """The valley top from the POOLED elevations of the whole landscape.

    Every tile sharing ground gets the same cut by construction, which per-tile computation
    cannot promise. Returns ``(valley_top_m, members)``; the elevation is None where the
    pooled histogram has no minimum above its dominant mode.
    """
    members = landscape_of(tile_dir, root=root)
    zs = []
    for d in members:
        p = os.path.join(root, d, "z_after.npy")
        if os.path.exists(p):
            z = np.load(p).ravel()
            zs.append(z[np.isfinite(z)])
    if not zs:
        return None, members
    return valley_top_from_histogram(np.concatenate(zs), bins=bins), members


def valley_top_from_histogram(z, *, bins=60):
    """Valley top from the elevation histogram: the FIRST LOCAL MINIMUM ABOVE THE DOMINANT
    MODE. Returns the elevation, or None if the tile has no minimum above its mode.

    The dominant mode is the valley floor wherever the valley is the flattest, most
    repeated elevation in the tile -- which is what a floodplain is. Walking up from it to
    the first dip finds where that population stops.

    It replaces "the trough between the tallest bin in each half of the elevation RANGE",
    which required a SECOND mode and had no test that one existed. Measured 2026-09-04:

        elba        dominant 222.2 m -> 226.9 m (17.2% of tile);
                    old rule 231.5 m (19.4%);  established VALLEY_TOP 230.0 m (18.8%)
        whitewater  dominant 202.3 m -> 207.7 m (35.8%);
                    old rule 297.9 m (76.5%) -- it found no second mode and cut the tile in
                    three quarters, taking the whole valley SIDE with the floor

    KNOWN LIMITS, to revisit (Andy, 2026-09-04):
      * THE PREMISE IS THAT THE COMMONEST ELEVATION IS THE VALLEY FLOOR, and it fails where
        it is not. Measured across the six sites 2026-09-04:

            elba        mode 222.2 m of 214-353 -> cut 226.9 m, removes 17.2%   ok
            whitewater  mode 202.3 m of 201-360 -> cut 207.7 m, removes 35.8%   ok
            carlton     mode 226.4 m of 209-317 -> cut 233.5 m, removes 22.2%   ok
            mnrv        mode 229.4 m of 226-304 -> cut 232.0 m, removes 28.5%   ok
            cook        mode 588.4 m of 447-606 -> NO minimum above it          refuses
            battlecreek mode 281.6 m of 244-304 -> cut 283.6 m, removes 75.2%   WRONG

        Cook is a lake-studded plateau: its commonest elevation is the upland, so there is
        no dip above the mode and this returns None, which makes reference_cells refuse.
        BATTLE CREEK IS THE BUILT ENVIRONMENT (Andy, 2026-09-04): graded lots and fill set
        the modal elevation, the cut lands 2 m above it, and 75.2% of the tile goes. Revisit
        for urban sites and for landscapes that are not valley-dominated.
      * A tile with NO valley still returns a cut, just above whatever its commonest
        elevation is. Check the fraction it removes before trusting it on flat ground.
      * It is a per-tile threshold, so it is NOT comparable between tiles. Use VALLEY_TOP_M
        where a landscape has an established value.
    """
    zc = np.asarray(z, float).ravel()
    zc = zc[np.isfinite(zc)]
    if zc.size <= 100:
        return None
    h, e = np.histogram(zc, bins=bins)
    c = 0.5 * (e[:-1] + e[1:])
    pk = int(np.argmax(h))
    i = pk + 1
    while i + 1 < len(h) and not (h[i] <= h[i - 1] and h[i] <= h[i + 1]):
        i += 1
    return float(c[i]) if i + 1 < len(h) else None

#: The valley-top ELEVATION per site, in metres. STATED, never inferred.
#:
#: Only values with a source are here. A site absent from this table makes the valley cut
#: RAISE rather than guess: a wrong valley top silently changes which ground counts as
#: "no-change", and every calibration downstream rests on that population.
VALLEY_TOP_M = {
    # Elba: established before this session -- analysis/steady_state/
    # run_steady_state_strata.py VALLEY_TOP, and ALLFOREST_BLUFFLAND.md.
    "elba": 230.0,
    "elba_fulldensity": 230.0,
    "elbaext": 230.0,
}


def floodplain_by_elevation(z, valley_top_m):
    """Valley-floor mask: cells below a STATED elevation. Returns ``(mask, valley_top_m)``.

    The elevation is given, never inferred. An earlier version read it off the antimode of
    the tile's own elevation histogram; that is withdrawn (Andy, 2026-09-04). It worked at
    elba only because that tile is genuinely bimodal -- a 25,601-cell valley peak at 222.2 m
    and a 20,570-cell upland peak at 338.3 m with a sharp 2,919-cell trough at 231.5 m. At
    whitewater the valley peak is 98,783 cells at 202.3 m and the "upland peak" is 6,545 at
    353.7 m, the tail end rather than a mode, so the histogram simply decays and the
    procedure returned an arbitrary point in it, cutting 76% of the tile.

    And TPI is not an alternative: measured at whitewater, the TPI mask misses 150,873 cells
    of valley floor while removing 16,218 upland hollows at a median 305.5 m. Never use TPI
    for floodplain (Andy, 2026-09-04).
    """
    z = np.asarray(z, float)
    t = float(valley_top_m)
    return np.isfinite(z) & (z < t), t


def _opt(d, name):
    p = os.path.join(d, name)
    return np.load(p) if os.path.exists(p) else None


def reference_cells(tile_dir, *, cells=None, curv_max=0.015, slope_max=12.0,
                    gross_change_mm=500.0, clearcut_drop=0.30, require_ridge=True,
                    exclude_valley=True, valley_top_m=None,
                    use_floodplain_mask=False):
    """Boolean mask of stable reference cells, plus a report of what each cut removed.

    ``cells`` is an optional array of flat cell indices (e.g. the near-ground cube's
    ``cells``); the mask is returned over those, otherwise over the whole raveled grid.
    Any criterion whose input the tile does not carry is skipped and recorded as such.

    Returns ``(mask, report)`` where ``report`` maps criterion name -> number of cells it
    removed from the running mask, in application order, with ``"start"`` and ``"kept"``.
    """
    sl = np.load(os.path.join(tile_dir, "slope.npy")).ravel()
    idx = np.arange(sl.size) if cells is None else np.asarray(cells)
    m = np.ones(idx.size, bool)
    rep = {"start": int(m.sum())}

    def cut(name, ok):
        nonlocal m
        before = int(m.sum())
        m = m & ok
        rep[name] = before - int(m.sum())

    ridge = _opt(tile_dir, "ridge_mask.npy")
    if require_ridge and ridge is not None:
        cut("not a divide cell", ridge.astype(bool).ravel()[idx])
    curv = _opt(tile_dir, "curv_laplacian.npy")
    if curv is not None:
        c = np.abs(curv.ravel()[idx])
        cut(f"|curv| > {curv_max:g}", np.isfinite(c) & (c <= curv_max))
    cut(f"slope >= {slope_max:g} deg", np.isfinite(sl[idx]) & (sl[idx] < slope_max))

    g1 = _opt(tile_dir, "gen1_canopy_frac.npz")
    g2 = _opt(tile_dir, "gen2_canopy_frac.npz")
    if g1 is not None and g2 is not None:
        cut("building returns", (g1["n_bldg"].ravel()[idx] == 0)
            & (g2["n_bldg"].ravel()[idx] == 0))
        cut(f"clear-cut (frac drop > {clearcut_drop:g})",
            (g1["frac"].ravel()[idx] - g2["frac"].ravel()[idx]) <= clearcut_drop)

    bluff = _opt(tile_dir, "blufftop_margin_mask.npy")
    if bluff is not None:
        cut("blufftop retreat margin", ~bluff.astype(bool).ravel()[idx])

    # Valley floor. Ridge-ness, curvature and slope do NOT exclude flat valley-bottom
    # ground, so terraces and floodplain enter the "divide" population and dominate it:
    # a 19% valley limb produced an easting gradient of -84.8 mm/km against +3.5 on the
    # upland (analysis/STABLE_POINT_TILT_AUDIT.md), and a 27%-floodplain flat-slope bin
    # produced a spurious +7.9 mm rise. Excluded by default for anything divide-based.
    # Cut by ELEVATION (Andy, 2026-09-04), not by the TPI floodplain mask. The mask is a
    # topographic-position heuristic -- TPI over an 800 m window < -2 m -- so what it removes
    # depends on the window and on how wide the valley is, and it keeps flat terrace ground
    # sitting at valley level. The ANTIMODE of this tile's own elevation histogram separates
    # the upland plateau from the valley floor on the quantity that actually defines a
    # floodplain. use_floodplain_mask now defaults to FALSE; pass True for an older
    # population.
    if exclude_valley:
        if use_floodplain_mask:
            fld = _opt(tile_dir, "floodplain_mask.npy")
            if fld is None:
                # It does not run rather than running differently. Skipping this in silence
                # made two tiles' populations differ by 39,038 cells with nothing to notice,
                # and every comparison between them was invalid without saying so.
                raise FileNotFoundError(
                    f"{tile_dir}/floodplain_mask.npy is missing, so the floodplain cut "
                    f"cannot be applied. A population WITHOUT it is not comparable with one "
                    f"that has it. Either produce the mask for this tile, or pass "
                    f"use_floodplain_mask=False to state that you are deliberately working "
                    f"without it.")
            cut("floodplain mask", ~fld.astype(bool).ravel()[idx])
        zf = _opt(tile_dir, "z_after.npy")
        if zf is None:
            raise FileNotFoundError(
                f"{tile_dir}/z_after.npy is missing, so the valley cut cannot be applied.")
        # THE CALLER ALWAYS SAYS WHICH (Andy, 2026-09-04). No silent fallback: a stated
        # value and a computed one are different claims, and a chain that quietly
        # substitutes one for the other hides which was used. Elba rebuilt on a computed
        # 226.9 m instead of its cited 230.0 m exactly that way, and nothing said so.
        #
        #   valley_top_m = <float>       an elevation in metres, stated
        #                = "registry"    VALLEY_TOP_M for this tile; raises if absent
        #                = "histogram"   computed from the LANDSCAPE's pooled elevations
        if valley_top_m is None:
            raise ValueError(
                "exclude_valley needs valley_top_m and will not choose for you. Pass an "
                "elevation in metres, or 'registry' to use the established value for this "
                "tile, or 'histogram' to compute it from the landscape's pooled "
                "elevations. Known registry values: "
                + ", ".join(f"{k}={v:g} m" for k, v in sorted(VALLEY_TOP_M.items()))
                + ". Or exclude_valley=False to work without the cut deliberately.")
        if valley_top_m == "registry":
            name = os.path.basename(str(tile_dir).rstrip("/"))
            valley_top_m = VALLEY_TOP_M.get(name)
            if valley_top_m is None:
                raise ValueError(
                    f"no established valley top for {name!r}. Add one to "
                    f"refcells.VALLEY_TOP_M with its source, or ask for 'histogram'.")
            src = "registry"
        elif valley_top_m == "histogram":
            valley_top_m, members = valley_top_for_landscape(tile_dir)
            if valley_top_m is None:
                raise ValueError(
                    f"the pooled elevation histogram for this landscape "
                    f"({', '.join(members)}) has no minimum above its dominant mode, so no "
                    f"valley top can be computed. State one, or exclude_valley=False.")
            src = f"histogram over {len(members)} tile(s)"
        else:
            src = "stated"
        z = zf.ravel()[idx]
        cut(f"below valley top {float(valley_top_m):.1f} m ({src})",
            ~(np.isfinite(z) & (z < float(valley_top_m))))

    dod = _opt(tile_dir, "dod.npy")
    if dod is not None and gross_change_mm is not None:
        d = dod.ravel()[idx] * 1000.0
        cut(f"|DoD| > {gross_change_mm:g} mm", np.isfinite(d)
            & (np.abs(d) <= gross_change_mm))

    rep["kept"] = int(m.sum())
    return m, rep
