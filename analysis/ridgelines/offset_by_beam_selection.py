#!/usr/bin/env python3
"""Which BEAM GEOMETRY reads the ground truest as terrain steepens?

For every 5 m cell, pick the one return closest to each of four end-member geometries and
follow its offset against slope. The four are two scanner-frame selections and two
surface-frame ones -- the distinction that matters on slopes, where a near-vertical beam
is NOT the one most perpendicular to the ground:

    near-nadir       min |scan_angle|   beam closest to vertical
    most-horizontal  max |scan_angle|   beam closest to horizontal (swath edge)
    most-perpendicular  min incidence   beam closest to the surface NORMAL
    most-parallel       max incidence   beam closest to grazing the surface

Selecting WITHIN a cell is the point: all four look at the same ground with the same real
change and the same datum, so the four curves differ only by beam geometry. A cell needs
several returns for the picks to be distinct, hence --min-n.

Defaults to the REGISTRATION-CORRECTED offset (d_mm_corr): the raw d_mm carries per-swath
misalignment, which is itself a between-flight-line effect and therefore contaminates
exactly this comparison -- the four picks often come from different flight lines.

    env -u PROJ_DATA -u GDAL_DATA ./lidar-icp/bin/python \
        analysis/ridgelines/offset_by_beam_selection.py --tile data/derived/elbaext
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tile", default="data/derived/elba_fulldensity")
ap.add_argument("--offset", default="corr", choices=("raw", "corr"))
ap.add_argument("--curv-max", type=float, default=None)
ap.add_argument("--ridge", action="store_true")
ap.add_argument("--min-n", type=int, default=3, help="minimum returns per cell to pick from")
ap.add_argument("--min-scan-spread", type=float, default=0.0,
                help="require this within-cell |scan angle| RANGE (deg). Half the cells are "
                     "single-flight-line with ~1 deg of spread, so the four selections there "
                     "are nearly the same beam and cannot discriminate; raising this keeps "
                     "the flight-line-overlap cells where the comparison has leverage")
ap.add_argument("--min-inc-spread", type=float, default=0.0,
                help="require this within-cell incidence RANGE (deg)")
ap.add_argument("--bin", type=float, default=2.0, help="slope bin width (deg)")
A = ap.parse_args()
TILE = os.path.basename(A.tile.rstrip("/"))
TAG = ("" if TILE == "elba_fulldensity" else f"_{TILE}")
TAG += "" if A.offset == "corr" else "_raw"
if A.curv_max is not None: TAG += f"_curv{A.curv_max:g}"
if A.min_scan_spread > 0: TAG += f"_sp{A.min_scan_spread:g}"
if A.ridge: TAG += "_ridge"

DCOL = "d_mm_corr" if A.offset == "corr" else "d_mm"
SLOPE_EDGES = np.arange(0, 46, A.bin)
MIN_BIN = 100

df = pd.read_parquet(f"{A.tile}/beam_offset_table.parquet",
                     columns=["cell", DCOL, "slope", "canopy_cover", "scan_angle",
                              "incidence", "curv_laplacian", "in_grid"])
df = df[df.in_grid.values].copy()
lab = f"{TILE}; {'registration-corrected' if A.offset == 'corr' else 'RAW (pre-registration)'}"
if A.curv_max is not None:
    df = df[(df.curv_laplacian.abs() <= A.curv_max).to_numpy()].copy()
    lab += f"; |Laplacian|<={A.curv_max:g}"
if A.ridge:
    rm = np.load(f"{A.tile}/ridge_mask.npy").astype(bool).ravel()
    df = df[rm[df.cell.to_numpy()]].copy(); lab += "; RIDGELINES"
n_per = df.groupby("cell")["slope"].transform("size").to_numpy()
df = df[n_per >= A.min_n].copy()
lab += f"; cells with >={A.min_n} returns"
if A.min_scan_spread > 0:
    g = df.groupby("cell")["scan_angle"]
    rng = (g.transform("max") - g.transform("min")).to_numpy()
    df = df[rng >= A.min_scan_spread].copy()
    lab += f"; scan spread >={A.min_scan_spread:g} deg"
if A.min_inc_spread > 0:
    g = df.groupby("cell")["incidence"]
    rng = (g.transform("max") - g.transform("min")).to_numpy()
    df = df[rng >= A.min_inc_spread].copy()
    lab += f"; incidence spread >={A.min_inc_spread:g} deg"

cell = df.cell.to_numpy(); d = df[DCOL].to_numpy(float)
slope = df.slope.to_numpy(float); cover = df.canopy_cover.to_numpy(float)
absa = np.abs(df.scan_angle.to_numpy(float)); inc = df.incidence.to_numpy(float)
print("=" * 88)
print(f"OFFSET BY BEAM SELECTION  [{lab}]")
print(f"{len(df):,} returns in {len(np.unique(cell)):,} cells")
print("=" * 88)


def pick(key, largest=False):
    """Index of the return that minimises (or maximises) `key` within each cell."""
    k = -np.asarray(key, float) if largest else np.asarray(key, float)
    order = np.lexsort((k, cell))
    cs = cell[order]
    first = np.ones(cs.size, bool); first[1:] = cs[1:] != cs[:-1]
    return order[first]


SEL = {"near-nadir": pick(absa), "most-horizontal": pick(absa, largest=True),
       "most-perpendicular": pick(inc), "most-parallel": pick(inc, largest=True)}
COLORS = {"near-nadir": "C0", "most-horizontal": "C3",
          "most-perpendicular": "C2", "most-parallel": "C1"}
for k, i in SEL.items():
    print(f"  {k:19s} n={i.size:,}  |scan| median {np.median(absa[i]):5.1f} deg   "
          f"incidence median {np.median(inc[i]):5.1f} deg")


def curve(idx, sub=None):
    """(centres, medians, counts) of offset vs slope for a selection, optional cover mask."""
    s, v = slope[idx], d[idx]
    if sub is not None:
        m = sub[idx]; s, v = s[m], v[m]
    c, med, nn = [], [], []
    for lo, hi in zip(SLOPE_EDGES[:-1], SLOPE_EDGES[1:]):
        m = (s >= lo) & (s < hi)
        if m.sum() >= MIN_BIN:
            c.append(0.5*(lo+hi)); med.append(np.median(v[m])); nn.append(int(m.sum()))
    return np.array(c), np.array(med), np.array(nn)


STRATA = [("all returns", None), ("open (cc<0.10)", np.isfinite(cover) & (cover < 0.10)),
          ("forest (cc>0.50)", np.isfinite(cover) & (cover > 0.50))]
for sname, smask in STRATA:
    print(f"\nmedian offset (mm) vs slope -- {sname}")
    rows = {k: curve(i, smask) for k, i in SEL.items()}
    cs = sorted({c for k in rows for c in rows[k][0]})
    print("   slope  " + "".join(f"{k:>20s}" for k in SEL))
    for c in cs:
        line = f"   {c:5.1f}  "
        for k in SEL:
            cc_, mm_, nn_ = rows[k]
            j = np.where(cc_ == c)[0]
            line += (f"{mm_[j[0]]:+9.1f}({nn_[j[0]]:>6,})" if j.size else f"{'--':>20s}")
        print(line)

fig, ax = plt.subplots(2, 2, figsize=(15, 10), dpi=130)
a0 = ax[0, 0]
hb = a0.hexbin(slope, d, gridsize=(60, 60), bins="log", cmap="viridis",
               extent=(0, 45, -300, 300), mincnt=1)
fig.colorbar(hb, ax=a0, label="log10 count")
bc, bm = [], []
for lo, hi in zip(np.arange(0, 45, 1.0), np.arange(1, 46, 1.0)):
    m = (slope >= lo) & (slope < hi)
    if m.sum() >= 50: bc.append(0.5*(lo+hi)); bm.append(np.median(d[m]))
a0.plot(bc, bm, "o-", color="crimson", ms=3.5, lw=1.3, label="binned median (all returns)")
a0.axhline(0, color="k", lw=.6); a0.set_xlim(0, 45); a0.set_ylim(-300, 300)
a0.set_xlabel("surface slope (deg)"); a0.set_ylabel("offset d (mm), gen1 vs gen2")
a0.set_title("ALL returns: per-beam offset vs slope"); a0.legend(fontsize=8)

for axi, (sname, smask) in zip((ax[0, 1], ax[1, 0], ax[1, 1]), STRATA):
    for k, i in SEL.items():
        c, m, nn = curve(i, smask)
        if c.size:
            axi.plot(c, m, "o-", ms=4, lw=1.5, color=COLORS[k], label=f"{k} (n={nn.sum():,})")
    axi.axhline(0, color="k", lw=.6); axi.set_xlim(0, 45)
    axi.set_xlabel("surface slope (deg)"); axi.set_ylabel("median offset d (mm)")
    axi.set_title(f"best beam per cell, by geometry — {sname}")
    axi.legend(fontsize=7.5); axi.grid(alpha=.3)
fig.suptitle(f"gen1 offset by BEAM SELECTION vs slope — {lab}", y=1.0)
out = f"figures/refdatum/offset_by_beam_selection{TAG}.png"
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"\nwrote {out}")
