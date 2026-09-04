"""Per-mark curve inputs PLUS the flight line each mark sits on, and the lines in a tile."""
import os, sys, numpy as np, pandas as pd, laspy
sys.path.insert(0, "analysis")
from control_mode_shift import STRUCT, BOX, marks
from lidar_diff_icp import groundq
from lidar_diff_icp.groundtruth.tie import _design

SET = "gen2_2021_control"

# --- 1. the flight lines present in Whitewater's gen2 tile ---------------------------
WW = "data/after/3dep_4358_fulltile.laz"
cnt = {}
with laspy.open(WW) as f:
    for pts in f.chunk_iterator(5_000_000):
        cl = np.asarray(pts.classification); ps = np.asarray(pts.point_source_id)
        u, c = np.unique(ps[cl == 2], return_counts=True)
        for k, v in zip(u.tolist(), c.tolist()):
            cnt[k] = cnt.get(k, 0) + v
tot = sum(cnt.values())
print(f"WHITEWATER gen2 class-2 returns {tot:,} over {len(cnt)} point_source_id(s)")
for k in sorted(cnt, key=lambda k: -cnt[k]):
    print(f"   psid {k:6d}   {cnt[k]:12,}   {100*cnt[k]/tot:5.1f}%")
np.save("data/derived/whitewater/gen2_psid_counts.npy",
        np.array(sorted(cnt.items()), dtype=np.int64))

# --- 2. every mark: the curve's two numbers, plus its own line -----------------------
rows = []
for t in marks(SET).itertuples():
    sp, bp = f"{STRUCT}/{SET}__{t.point_id}.npz", f"{BOX}/{SET}__{t.point_id}.laz"
    if not (os.path.exists(sp) and os.path.exists(bp)):
        continue
    z = np.load(sp); coef = z["surface_coef"]
    E, N, R = float(z["easting"]), float(z["northing"]), float(z["struct_radius"])
    f = laspy.read(bp)
    x, y, zz = np.asarray(f.x), np.asarray(f.y), np.asarray(f.z)
    cl, ps = np.asarray(f.classification), np.asarray(f.point_source_id)
    g = (np.hypot(x - E, y - N) <= R) & (cl == 2)
    if g.sum() < 20:
        continue
    nn = np.sqrt(1 + coef[1]**2 + coef[2]**2)
    hg = np.sort((zz[g] - (_design(x[g]-E, y[g]-N, 2) @ coef)) / nn)
    mu = (float(t.elevation) - float(coef[0])) / nn
    u, c = np.unique(ps[g], return_counts=True)
    rows.append(dict(point_id=t.point_id, easting=t.easting, northing=t.northing, mu=mu,
                     psid_dominant=int(u[np.argmax(c)]),
                     psid_dominant_frac=float(c.max() / c.sum()),
                     psid_all="|".join(str(int(v)) for v in u),
                     n_class2=int(g.sum()), **groundq.mark_statistics(hg, mu)))
F = pd.DataFrame(rows)
F.to_csv("data/derived/control_marks_by_line_gen2_2021_control.csv", index=False)
print(f"\nwrote data/derived/control_marks_by_line_gen2_2021_control.csv  ({len(F)} marks)")
