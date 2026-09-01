#!/usr/bin/env python3
"""The cover calibration, as measured: what the binned medians do, and which candidate
form follows them.

Reads ONLY the cover_offset_calibration*.json files that cover_offset_reference.py writes,
so nothing here re-fits anything -- every point, error bar and curve is that file's own
content. Two things it makes visible that the per-tile figures do not:

  1. The two tiles side by side. Their linear coefficients now agree to 0.11 mm per unit
     cover; the constants previously TYPED into dod_cover_attribution.py (49.6 and 48.4)
     implied a tile-to-tile difference of 1.2. Those two are annotated as numbers, NOT
     drawn as curves -- only their slopes were ever recorded, so a line through them would
     need an intercept I would have to invent.
  2. Why the linear form is the wrong one. The residual panels show the misfit that the
     AIC comparison reports: linear leaves a systematic arc, the selected form does not.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/plot_cover_calibration.py
"""
import argparse, json, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tiles", nargs="+",
                default=["data/derived/elba_fulldensity", "data/derived/elbaext"])
ap.add_argument("--out", default="figures/refdatum/cover_calibration_forms.png")
# The constants this replaced, for annotation only. They are historical literals from
# dod_cover_attribution.py before commit 97798e4, given as DoD-gain k = -b.
PRIOR_K = {"elba_fulldensity": 49.6, "elbaext": 48.4}
A = ap.parse_args()

DESIGN = {                       # the same three forms cover_offset_reference.py compares
    "linear":        lambda c: np.c_[np.ones_like(c), c],
    "quadratic":     lambda c: np.c_[np.ones_like(c), c, c**2],
    "optical depth": lambda c: np.c_[np.ones_like(c), -np.log(1 - np.clip(c, 0, 0.98))],
}
COL = {"elba_fulldensity": "C0", "elbaext": "C1"}
FORM_STYLE = {"linear": ("--", "C3"), "quadratic": ("-", "C2"),
              "optical depth": ("-", "C4")}


def load(tile_dir):
    name = os.path.basename(tile_dir.rstrip("/"))
    for fn in (f"cover_offset_calibration_{name}.json", "cover_offset_calibration.json"):
        p = os.path.join(tile_dir, fn)
        if os.path.exists(p):
            return name, json.load(open(p)), p
    raise SystemExit(f"no cover_offset_calibration*.json under {tile_dir}; produce it with "
                     f"cover_offset_reference.py --tile {tile_dir}")


def short(form_key):
    return form_key.split()[0] if not form_key.startswith("optical") else "optical depth"


cals = [load(t) for t in A.tiles]

fig = plt.figure(figsize=(13.5, 6.2), dpi=130)
gs = fig.add_gridspec(len(cals), 2, width_ratios=[1.55, 1.0], hspace=0.32, wspace=0.24)
axA = fig.add_subplot(gs[:, 0])

for i, (name, cal, path) in enumerate(cals):
    b = cal["binned"]
    cc = np.asarray(b["cover"]); mm = np.asarray(b["median_mm"])
    se = np.asarray(b["robust_se_mm"]); nn = np.asarray(b["n"])
    sel = cal["selected_form"]
    coefs = {short(k): np.asarray(v) for k, v in cal["forms"].items()}
    xs = np.linspace(0, cc.max(), 300)

    axA.errorbar(cc, mm, yerr=se, fmt="o", ms=5, lw=1.3, capsize=3, color=COL[name],
                 label=f"{name}  ({cal['population']['n_returns']:,} returns)", zorder=3)
    for fk in ("linear", short(sel)):
        ls, _ = FORM_STYLE[fk]
        axA.plot(xs, DESIGN[fk](xs) @ coefs[fk], ls, lw=2.0 if fk != "linear" else 1.4,
                 color=COL[name], alpha=1.0 if fk != "linear" else 0.55, zorder=2)

    # ---- residual panel: every candidate form, in units of the bin's own SE ----
    ax = fig.add_subplot(gs[i, 1])
    for fk, coef in coefs.items():
        ls, c = FORM_STYLE[fk]
        r = (mm - DESIGN[fk](cc) @ coef) / se
        ax.plot(cc, r, "o" + ls, ms=4, lw=1.2, color=c,
                label=f"{fk}{'  (selected)' if fk == short(sel) else ''}")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_ylabel("residual / SE")
    ax.set_title(f"{name} — misfit by form", fontsize=9)
    ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=1)
    if i == len(cals) - 1:
        ax.set_xlabel("PyForestScan canopy cover (fraction)")

axA.axhline(0, color="k", lw=0.7)
axA.set_xlabel("PyForestScan canopy cover (fraction)")
axA.set_ylabel("median offset d_mm_corr (mm)   [gen1 − gen2]")
axA.set_title("offset vs canopy cover on non-eroding reference ground\n"
              "dashed = linear; solid = the form each tile's own AIC selects", fontsize=10)
axA.grid(alpha=0.3)

note = ["linear b (mm per unit cover), and the constant it replaced:"]
for name, cal, _ in cals:
    lin = [v for k, v in cal["forms"].items() if k.startswith("linear")][0]
    note.append(f"   {name:<18s} read {-lin[1]:6.2f}    previously typed {PRIOR_K[name]:.1f}")
note.append("selected form, by AIC on the weighted binned medians:")
for name, cal, _ in cals:
    note.append(f"   {name:<18s} {short(cal['selected_form'])}")
axA.text(0.02, 0.03, "\n".join(note), transform=axA.transAxes, fontsize=7.5,
         family="monospace", va="bottom",
         bbox=dict(fc="white", ec="0.7", alpha=0.9, boxstyle="round,pad=0.45"))
axA.legend(fontsize=8, loc="upper right")

os.makedirs(os.path.dirname(A.out), exist_ok=True)
fig.savefig(A.out, bbox_inches="tight"); plt.close(fig)
print(f"wrote {A.out}")
for name, cal, path in cals:
    lin = [v for k, v in cal["forms"].items() if k.startswith("linear")][0]
    print(f"  {name:<18s} linear b {lin[1]:+.4f}  ->  k {-lin[1]:.2f}   "
          f"selected {short(cal['selected_form'])}   (from {path})")
