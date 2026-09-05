"""One record per site: the clouds, the extent, and how its valley cut is decided.

WHY THIS EXISTS. The per-site configuration was split across scripts/run_all_sites.py --
SITES held the clouds and bounds, VALLEY_TOP_SOURCE held the cut, tile directories and CSF
cache paths were rebuilt from the name at each use -- and difference_dem itself had no
notion of WHICH TILE it was running on. When terrain.terrain_masks needed tile identity to
resolve a valley top, that gap showed up as a parameter called `tile_dir_for_landscape`: a
patch over a missing concept. A Site carries the identity, so nothing has to reconstruct it.
"""
from __future__ import annotations

import os
from dataclasses import dataclass

__all__ = ["Site", "SITES", "site"]


@dataclass(frozen=True)
class Site:
    """Everything a run needs to know about one tile, in one place.

    ``valley_top`` is how the floodplain cut is DECIDED, never the cut itself unless it is a
    number: an elevation in metres, ``"registry"`` (an established, cited value), or
    ``"histogram"`` (computed from the landscape's pooled elevations). It is stated here per
    site so the choice is visible in one table rather than implied by an argument default.
    """
    name: str
    gen1: str
    gen2: str
    valley_top: object
    bounds: tuple | None = None
    stream: bool = True
    derived_root: str = "data/derived"
    csf_root: str = "data/csf_cache"

    @property
    def tile_dir(self) -> str:
        return os.path.join(self.derived_root, self.name)

    @property
    def csf_cache(self) -> str:
        return os.path.join(self.csf_root, f"{self.name}.las")


SITES: dict[str, Site] = {
    # gen2 was 3dep2021_fulltile.laz until 2026-09-02: 16,785,971 points, 1.89 pts/m2,
    # against 182,923,322 (20.55 pts/m2) sitting beside it in the same extent -- a 10.90x
    # thinning, and 19.27x on the ground class (0.30 vs 5.78 pts/m2). Every other site
    # already read its full cloud, so Elba alone was differenced from something 6.03-27.53x
    # thinner than its peers: a cross-site method artifact before it is an uncertainty
    # problem. It also put 109,894 of 339,950 DoD cells (32.33%) below cell_plane_roughness
    # min_n=6, so they lost their LoD entirely.
    # valley_top "registry": its landscape has an established 230.0 m
    # (run_steady_state_strata.py VALLEY_TOP, ALLFOREST_BLUFFLAND.md).
    "elba": Site("elba", "data/before/4342-29-64.laz",
                 "data/after/3dep2021_fulldensity.laz", "registry",
                 (577492.8, 4882737.6, 580032.8, 4886237.6)),
    # 3dep_4358_fulltile.laz was a TRUNCATED fetch: 5.52 returns/m2 east of easting 586362
    # against 15.45 west, same five flight lines both sides. Replaced 2026-09-04 by an
    # uncapped re-fetch, 148,050,625 points, west/east density ratio 1.06.
    "whitewater": Site("whitewater", "data/before/4358-26-03.laz",
                       "data/after/3dep_4358_fulldensity.laz", "histogram"),
    "mnrv": Site("mnrv", "data/before_mnrv/4342-23-01.laz",
                 "data/after_mnrv/mnrv_3dep2021.laz", "histogram"),
    # cook's pooled histogram has NO minimum above its dominant mode -- a lake-studded
    # plateau, mode 588.4 m of a 447-606 m range -- so "histogram" RAISES here. It needs a
    # stated elevation before it can build. Left as-is rather than guessed.
    "cook": Site("cook", "data/before_ne/1158-31-59.laz",
                 "data/after_ne/ne_3dep_fulldensity.laz", "histogram",
                 (709531.0, 5323589.0, 711986.0, 5327144.0)),
    "carlton": Site("carlton", "data/before_carlton/2742-12-53.laz",
                    "data/after_carlton/carlton_3dep.laz", "histogram",
                    (547805.0, 5163676.0, 550225.0, 5167166.0)),
    # battlecreek's histogram cut removes 75.2% of the tile: a BUILT ENVIRONMENT, where
    # graded lots set the modal elevation rather than a valley floor. It needs a stated
    # elevation or a fraction guard; the value below is not trusted.
    #
    # stream=True since 2026-09-05 (Andy: "unblock"). It was False because the tile is the
    # smallest here -- 615 x 870 m, 8.4 M points -- and fits in memory. But the in-memory
    # path grids with pandas groupby.quantile, which takes ONE quantile for all cells, so
    # ground_q="calibrated" refuses on it: battlecreek alone could not run the vegetation
    # correction, and a six-site run raised there. Streaming costs nothing at this size.
    # EXPECT ITS PRODUCTS TO MOVE SLIGHTLY on the next rebuild: the streaming route bins the
    # column, so `spread` and the ground grid differ from the exact in-memory computation by
    # under one bin. That is a route change, not an error, and it makes battlecreek
    # comparable with the other five, which have always streamed.
    "battlecreek": Site("battlecreek", "data/before_battlecreek/4342-03-32_b_a.laz",
                        "data/after_battlecreek/battlecreek_3dep.laz", "histogram",
                        (498750.0, 4975136.0, 499365.0, 4976006.0)),
}


def site(name: str) -> Site:
    """The Site record for ``name``, refusing rather than inventing one."""
    try:
        return SITES[name]
    except KeyError:
        raise KeyError(f"no site {name!r}; known: {', '.join(sorted(SITES))}") from None
