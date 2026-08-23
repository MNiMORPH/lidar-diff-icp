"""Driver: zero-curvature "steady-state" DoD check on TWO additional strata.

Reuses steady_state_cells.py (module) exactly as run_steady_state.py does, but
applies the same geomorphic steady-state selection to:

  (1) ALL FOREST  -- every forest cell (penetration < 0.25 & NOT floodplain),
      NOT just the eroded/core-forest interior.  Tests whether the +19.8 mm
      core-forest steady-state offset holds over the full forest population,
      which includes canopy edges and marginal cells.

  (2) BLUFFLAND FARMLAND -- open/farmland (penetration >= 0.45 & NOT floodplain)
      located in the DISSECTED BLUFF terrain, as distinct from the flat UPLAND
      plateau.  Definition (stated & defended below).

DATUM: GEOID tie ONLY  (data/derived/elba_refdatum/dod_geoid.npy, gen2 - gen1).
The parabola tie is a flexible warp that absorbs real hillslope change and is
deliberately NOT used here.

For context the figure also overlays the known reference curves:
  core forest         (+19.8 mm)
  upland farmland      (-4.9 mm, prior 66th-pct elevation split of core_open;
                        here reproduced as core_open on the plateau, z >= 330 m)

CONFOUNDS (flagged, not corrected):
  * DoD, curvature, and slope are all derived from the gen2 DEM -> a shared-frame
    confound: the "planar" cells are planar in gen2, not independently in gen1.
  * Any per-cell scan-incidence / point-density artifact between epochs is not
    removed and can bias a canopy-vs-open contrast.

Run:
    ./lidar-icp/bin/python analysis/steady_state/run_steady_state_strata.py
"""

from __future__ import annotations

import os
import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from steady_state_cells import (
    steady_state_mask, extract_diff, diff_stats, eps_curv_from_quantile,
)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FD = os.path.join(ROOT, "data", "derived", "elba_fulldensity")
RD = os.path.join(ROOT, "data", "derived", "elba_refdatum")
FIGDIR = os.path.join(ROOT, "figures", "refdatum")
os.makedirs(FIGDIR, exist_ok=True)

MAX_SLOPE = 15.0       # deg, conservative mass-wasting exclusion
CENTRAL_FRAC = 0.30    # central 30% |curvature| band -> "genuinely planar"

# bluffland-definition thresholds (see block comment in the bluffland section)
PLATEAU_BASE = 330.0   # m; open plateau sits ~330-348 m (dense mode)
VALLEY_TOP = 230.0     # m; valley-floor / near-floodplain terrace sits <~230 m
BLUFF_LO, BLUFF_HI = VALLEY_TOP, PLATEAU_BASE   # 230 <= z < 330 -> dissected bluff

# ---------------------------------------------------------------------------
# load  (GEOID datum only)
# ---------------------------------------------------------------------------
curv = np.load(os.path.join(FD, "curv_laplacian.npy"))
slope = np.load(os.path.join(FD, "slope.npy"))
z_after = np.load(os.path.join(FD, "z_after.npy"))           # gen2 DEM
pen = np.load(os.path.join(FD, "penetration.npy"))
floodplain = np.load(os.path.join(FD, "floodplain_mask.npy")).astype(bool)
core_forest = np.load(os.path.join(FD, "core_forest.npy")).astype(bool)
core_open = np.load(os.path.join(FD, "core_open.npy")).astype(bool)
dod_geoid = np.load(os.path.join(RD, "dod_geoid.npy"))       # PRIMARY & ONLY datum

finite_pen = np.isfinite(pen)
finite_z = np.isfinite(z_after)

# ---------------------------------------------------------------------------
# cover masks
# ---------------------------------------------------------------------------
# ALL FOREST: every forest cell, not just the eroded interior.
all_forest = (pen < 0.25) & (~floodplain) & finite_pen

# open / farmland base per spec
open_farm = (pen >= 0.45) & (~floodplain) & finite_pen

# BLUFFLAND FARMLAND definition
# -----------------------------
# The open-farmland elevation distribution is strongly bimodal: a dense flat
# PLATEAU at ~330-348 m and a low VALLEY FLOOR / near-floodplain terrace at
# ~215-225 m, with the DISSECTED BLUFF terrain spanning the sparse middle
# (~230-325 m) where the steep (>20 deg) bluff faces are.  We therefore define
# BLUFFLAND farmland as open cells at MID elevation, 230 <= z < 330 m:
#   - below the plateau base (excludes the flat upland plateau), and
#   - above the valley-floor terrace (excludes the high-scatter alluvial /
#     floodplain-margin band, which carries a hydrologic, not hillslope, signal).
# Justification that this band really is bluff terrain (not an arbitrary slice):
#   52.5% of these cells lie within 30 m of a steep (>20 deg) bluff face, vs only
#   3% of the plateau cells -- printed below.
bluffland = open_farm & finite_z & (z_after >= BLUFF_LO) & (z_after < BLUFF_HI)

# UPLAND plateau farmland (context / known -4.9 mm reference).  The prior
# reference used a 66th-pct elevation split of core_open; here reproduced as
# core_open restricted to the plateau (z >= 330 m).
upland_open = core_open & finite_z & (z_after >= PLATEAU_BASE)

# adjacency diagnostic for the bluffland definition
steep = (slope > 20.0) & np.isfinite(slope)
near_bluff = ndimage.binary_dilation(steep, iterations=6)   # ~30 m of a steep face
plateau_open = open_farm & finite_z & (z_after >= PLATEAU_BASE)
frac_bluff_adj = (bluffland & near_bluff).sum() / max(bluffland.sum(), 1)
frac_plat_adj = (plateau_open & near_bluff).sum() / max(plateau_open.sum(), 1)

# ---------------------------------------------------------------------------
# per-stratum steady-state selection (eps_curv from that stratum's own band)
# ---------------------------------------------------------------------------
def run_stratum(name, cover):
    base = (cover & (slope < MAX_SLOPE)
            & np.isfinite(curv) & np.isfinite(slope) & np.isfinite(dod_geoid))
    eps = eps_curv_from_quantile(curv, base, central_frac=CENTRAL_FRAC)
    mask = steady_state_mask(curv, slope, dod_geoid, cover, eps, MAX_SLOPE)
    diff = extract_diff(dod_geoid, mask)
    s = diff_stats(diff)
    return dict(name=name, cover_n=int(cover.sum()), base_n=int(base.sum()),
                eps=eps, mask=mask, diff=diff, stats=s)

r_forest = run_stratum("ALL FOREST", all_forest)
r_bluff = run_stratum("BLUFFLAND FARMLAND", bluffland)
r_core = run_stratum("core forest (context)", core_forest)
r_upl = run_stratum("upland farmland (context)", upland_open)

# ---------------------------------------------------------------------------
# print
# ---------------------------------------------------------------------------
def _row(r):
    s = r["stats"]
    return (f"{r['name']:26s} cover_n={r['cover_n']:6d}  ss_n={s['n']:5d}  "
            f"median={s['median_mm']:+7.1f}  NMAD={s['nmad_mm']:6.1f}  "
            f"mean={s['mean_mm']:+7.1f}  eps={r['eps']:.5f} 1/m  (mm)")

print("=" * 90)
print("Steady-state (zero-curvature, low-slope) DoD check -- GEOID datum only "
      "(gen2 - gen1)")
print("=" * 90)
print(f"max_slope = {MAX_SLOPE:.0f} deg ; eps_curv = central {int(CENTRAL_FRAC*100)}% "
      f"|kappa| band, computed per stratum")
print()
print("BLUFFLAND definition: open farmland (penetration>=0.45 & NOT floodplain) at "
      f"mid elevation")
print(f"  {BLUFF_LO:.0f} <= z < {BLUFF_HI:.0f} m  (below the ~{PLATEAU_BASE:.0f} m "
      f"plateau base, above the ~{VALLEY_TOP:.0f} m valley-floor terrace).")
print(f"  Adjacency check: {frac_bluff_adj*100:.1f}% of bluffland cells within 30 m "
      f"of a steep (>20 deg) bluff face,")
print(f"                   vs {frac_plat_adj*100:.1f}% of plateau farmland cells "
      f"-> band is genuinely dissected bluff terrain.")
print()
print(_row(r_forest))
print(_row(r_bluff))
print(_row(r_core))
print(_row(r_upl))
print()
print("CONFOUNDS: DoD/curv/slope all from the gen2 DEM (shared-frame; 'planar' is")
print("  planar in gen2, not independently in gen1). Per-cell scan-incidence /")
print("  point-density differences between epochs are not removed.")

# ---------------------------------------------------------------------------
# figure: overlaid PDFs
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=110)
bins = np.linspace(-250, 250, 101)  # mm


def _hist(diff, **kw):
    ax.hist(diff * 1e3, bins=bins, density=True, **kw)


# context (thin, faded) first
_hist(r_upl["diff"], color="C2", histtype="step", lw=1.3, ls=":",
      alpha=0.85,
      label=f"upland farmland [ctx] (n={r_upl['stats']['n']}, "
            f"med {r_upl['stats']['median_mm']:+.1f})")
_hist(r_core["diff"], color="0.45", histtype="step", lw=1.3, ls=":",
      alpha=0.9,
      label=f"core forest [ctx] (n={r_core['stats']['n']}, "
            f"med {r_core['stats']['median_mm']:+.1f})")
# primary strata (bold)
_hist(r_forest["diff"], color="C0", histtype="step", lw=2.2,
      label=f"ALL forest (n={r_forest['stats']['n']}, "
            f"med {r_forest['stats']['median_mm']:+.1f})")
_hist(r_bluff["diff"], color="C3", histtype="step", lw=2.2,
      label=f"bluffland farmland (n={r_bluff['stats']['n']}, "
            f"med {r_bluff['stats']['median_mm']:+.1f})")

ax.axvline(0, color="k", lw=0.9, alpha=0.7)
ax.axvline(r_forest["stats"]["median_mm"], color="C0", lw=1.1, ls="--", alpha=0.8)
ax.axvline(r_bluff["stats"]["median_mm"], color="C3", lw=1.1, ls="--", alpha=0.8)

ann = (f"ALL forest:        median {r_forest['stats']['median_mm']:+6.1f} mm  "
       f"NMAD {r_forest['stats']['nmad_mm']:5.1f}  n={r_forest['stats']['n']}\n"
       f"bluffland farmland: median {r_bluff['stats']['median_mm']:+6.1f} mm  "
       f"NMAD {r_bluff['stats']['nmad_mm']:5.1f}  n={r_bluff['stats']['n']}\n"
       f"core forest [ctx]:  median {r_core['stats']['median_mm']:+6.1f} mm  "
       f"n={r_core['stats']['n']}\n"
       f"upland farm [ctx]:  median {r_upl['stats']['median_mm']:+6.1f} mm  "
       f"n={r_upl['stats']['n']}")
ax.text(0.015, 0.975, ann, transform=ax.transAxes, va="top", ha="left",
        fontsize=7.6, family="monospace",
        bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))

ax.set_xlabel("gen2 - gen1 elevation difference (mm)   [GEOID datum]")
ax.set_ylabel("probability density (1/mm)")
ax.set_title("Steady-state (planar, slope<15 deg) DoD PDFs: all forest vs bluffland "
             "farmland\n"
             f"eps_curv = central {int(CENTRAL_FRAC*100)}% |kappa| band per stratum "
             "; GEOID datum only")
ax.legend(fontsize=7.0, loc="upper right")
ax.set_xlim(-250, 250)
fig.tight_layout()

out = os.path.join(FIGDIR, "steady_state_allforest_bluffland_pdf.png")
fig.savefig(out, dpi=110)
w, h = fig.canvas.get_width_height()
print()
print(f"Wrote {out}  ({w}x{h} px)")

# emit machine-readable numbers for the markdown writer
print()
print("MDNUMS " + repr(dict(
    forest=r_forest["stats"], bluff=r_bluff["stats"],
    core=r_core["stats"], upland=r_upl["stats"],
    eps_forest=r_forest["eps"], eps_bluff=r_bluff["eps"],
    frac_bluff_adj=frac_bluff_adj, frac_plat_adj=frac_plat_adj,
    cover_forest=r_forest["cover_n"], cover_bluff=r_bluff["cover_n"],
)))
