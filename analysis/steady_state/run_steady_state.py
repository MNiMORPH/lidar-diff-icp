"""Driver: steady-state (planar, low-slope) forest cells as a datum / DoD check.

Selects core-forest cells that are locally planar (|grad^2 z| < eps_curv) and
low-slope (< 15 deg), then reports the PDF of the gen2 - gen1 elevation
difference over them for BOTH datum products (geoid tie [primary] and parabola
tie).  See steady_state_cells.py for the geomorphic rationale and caveats.

Run:
    ./lidar-icp/bin/python analysis/steady_state/run_steady_state.py
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from steady_state_cells import (
    steady_state_mask, extract_diff, diff_stats,
    eps_curv_from_quantile, eps_curv_from_diffusion,
)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FD = os.path.join(ROOT, "data", "derived", "elba_fulldensity")
RD = os.path.join(ROOT, "data", "derived", "elba_refdatum")
FIGDIR = os.path.join(ROOT, "figures", "refdatum")
os.makedirs(FIGDIR, exist_ok=True)

MAX_SLOPE = 15.0          # deg, conservative mass-wasting exclusion
CENTRAL_FRAC = 0.30       # central 30% |curvature| band -> "genuinely planar"
DT_YR = 13.0              # 2008 -> 2021 interval

# ---------------------------------------------------------------------------
# load
# ---------------------------------------------------------------------------
curv = np.load(os.path.join(FD, "curv_laplacian.npy"))
slope = np.load(os.path.join(FD, "slope.npy"))
core_forest = np.load(os.path.join(FD, "core_forest.npy"))
z_after = np.load(os.path.join(FD, "z_after.npy"))
dod_geoid = np.load(os.path.join(RD, "dod_geoid.npy"))      # PRIMARY datum
dod_parab = np.load(os.path.join(FD, "dod.npy"))            # parabola-tie datum

# base population for choosing eps_curv: core forest, low slope, finite data
base = (core_forest & (slope < MAX_SLOPE)
        & np.isfinite(curv) & np.isfinite(slope) & np.isfinite(dod_geoid))

# ---------------------------------------------------------------------------
# eps_curv: principled choice + diffusive-signal sanity budget
# ---------------------------------------------------------------------------
eps_quant = eps_curv_from_quantile(curv, base, central_frac=CENTRAL_FRAC)

print("=" * 78)
print("eps_curv selection")
print("=" * 78)
print(f"Base population (core forest, slope<{MAX_SLOPE:.0f} deg, finite): {base.sum()} cells")
print(f"Curvature over base: median={np.median(curv[base]):+.5f}  "
      f"NMAD={1.4826*np.median(np.abs(curv[base]-np.median(curv[base]))):.5f} 1/m")
print()
print(f"PRIMARY eps_curv = central {int(CENTRAL_FRAC*100)}% |kappa| band = "
      f"{eps_quant:.5f} 1/m")
print()
print("Diffusive-signal budget check  (dz/dt = K*kappa*dt, dt=13 yr):")
print("  eps such that a planar-band cell's diffusion signal stays < 5 mm:")
for K in (0.002, 0.010, 0.050):
    eps_diff = eps_curv_from_diffusion(K, DT_YR, 0.005)
    sig_at_eps = K * eps_quant * DT_YR * 1e3  # signal at our chosen band edge (mm)
    print(f"    K={K:5.3f} m2/yr -> eps_diff={eps_diff:.4f} 1/m ; "
          f"signal at chosen band edge = {sig_at_eps:5.2f} mm")
print("  (Even for stiff K, the diffusion signal at our planar-band edge is far")
print("   below the ~50 mm forest DoD NMAD, so the steady-state read is not")
print("   limited by residual diffusion but by datum + scan-geometry noise.)")

# report cell counts for several eps choices
print()
print("Cells selected (core forest, slope<15) vs eps_curv choice:")
for frac in (0.20, 0.30, 0.40):
    e = eps_curv_from_quantile(curv, base, central_frac=frac)
    m = steady_state_mask(curv, slope, dod_geoid, core_forest, e, MAX_SLOPE)
    print(f"    central {int(frac*100)}%: eps={e:.5f} 1/m -> n={m.sum()}")

EPS = eps_quant  # use central-30% band

# ---------------------------------------------------------------------------
# steady-state mask + stats for both datums
# ---------------------------------------------------------------------------
mask_geoid = steady_state_mask(curv, slope, dod_geoid, core_forest, EPS, MAX_SLOPE)
mask_parab = steady_state_mask(curv, slope, dod_parab, core_forest, EPS, MAX_SLOPE)

# all-core-forest contrast set (slope<15, finite dod), no curvature restriction
allcf_geoid = (core_forest & (slope < MAX_SLOPE) & np.isfinite(dod_geoid))
allcf_parab = (core_forest & (slope < MAX_SLOPE) & np.isfinite(dod_parab))

diff_ss_g = extract_diff(dod_geoid, mask_geoid)
diff_ss_p = extract_diff(dod_parab, mask_parab)
diff_all_g = extract_diff(dod_geoid, allcf_geoid)
diff_all_p = extract_diff(dod_parab, allcf_parab)

s_ss_g = diff_stats(diff_ss_g)
s_ss_p = diff_stats(diff_ss_p)
s_all_g = diff_stats(diff_all_g)
s_all_p = diff_stats(diff_all_p)

def _row(label, s):
    return (f"{label:34s} n={s['n']:5d}  median={s['median_mm']:+7.1f}  "
            f"NMAD={s['nmad_mm']:6.1f}  mean={s['mean_mm']:+7.1f}  "
            f"IQR={s['iqr_mm']:6.1f}  [Q25 {s['q25_mm']:+.1f}, Q75 {s['q75_mm']:+.1f}]  (mm)")

print()
print("=" * 78)
print(f"gen2 - gen1 elevation difference over cells   (eps_curv={EPS:.5f} 1/m, "
      f"slope<{MAX_SLOPE:.0f} deg)")
print("=" * 78)
print(_row("STEADY-STATE, geoid datum", s_ss_g))
print(_row("STEADY-STATE, parabola datum", s_ss_p))
print(_row("ALL core forest, geoid datum", s_all_g))
print(_row("ALL core forest, parabola datum", s_all_p))

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=120)
bins = np.linspace(-250, 250, 101)  # mm

def _hist(ax, diff, **kw):
    ax.hist(diff * 1e3, bins=bins, density=True, **kw)

# all-core-forest contrast (geoid) as a faint filled backdrop
_hist(ax, diff_all_g, color="0.75", alpha=0.6, histtype="stepfilled",
      label=f"all core forest, geoid (n={s_all_g['n']})")
# steady-state, geoid (primary)
_hist(ax, diff_ss_g, color="C0", histtype="step", lw=2.0,
      label=f"steady-state, geoid (n={s_ss_g['n']})")
# steady-state, parabola (datum comparison)
_hist(ax, diff_ss_p, color="C3", histtype="step", lw=1.6, ls="--",
      label=f"steady-state, parabola (n={s_ss_p['n']})")

ax.axvline(0, color="k", lw=0.8, alpha=0.6)
ax.axvline(s_ss_g["median_mm"], color="C0", lw=1.2, ls=":")
ax.axvline(s_ss_p["median_mm"], color="C3", lw=1.2, ls=":")

ann = (f"steady-state (geoid): median {s_ss_g['median_mm']:+.1f} mm, "
       f"NMAD {s_ss_g['nmad_mm']:.1f} mm, n={s_ss_g['n']}\n"
       f"steady-state (parabola): median {s_ss_p['median_mm']:+.1f} mm, "
       f"NMAD {s_ss_p['nmad_mm']:.1f} mm, n={s_ss_p['n']}\n"
       f"all core forest (geoid): median {s_all_g['median_mm']:+.1f} mm, "
       f"NMAD {s_all_g['nmad_mm']:.1f} mm")
ax.text(0.015, 0.975, ann, transform=ax.transAxes, va="top", ha="left",
        fontsize=8.2, family="monospace",
        bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))

ax.set_xlabel("gen2 - gen1 elevation difference (mm)")
ax.set_ylabel("probability density (1/mm)")
ax.set_title("Steady-state planar low-slope forest cells: elevation-difference PDF\n"
             f"|grad^2 z| < {EPS:.4f} 1/m (central {int(CENTRAL_FRAC*100)}% band), "
             f"slope < {MAX_SLOPE:.0f} deg, core forest")
ax.legend(fontsize=7.6, loc="upper right")
ax.set_xlim(-250, 250)
fig.tight_layout()

out = os.path.join(FIGDIR, "steady_state_diff_pdf.png")
fig.savefig(out, dpi=120)
w, h = fig.canvas.get_width_height()
print()
print(f"Wrote {out}  ({w}x{h} px)")
