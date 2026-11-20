#!/usr/bin/env python
"""How often does one 2021 pulse give more than one return over the floodplain, and how close can two returns of one
pulse be? (Andy, 2026-11-20.) Pulses are identified from the RAW clouds (full-precision gps_time; the per-return CSVs
round it to 0.1 ms, which merges hundreds of pulses at a 2 MHz pulse rate) by (point_source_id, gps_time) inside the
bounding box of the 147 floodplain FIX points, tiles 4342-29-64 and 4358-29-01 (overlap kept from the first tile).
For each point, the returns within 3 m; a pulse counts as fully seen when the returns in the window equal its
number_of_returns. Points are split into sward (canopy share below 5 %, from the shape table) and tree edge; 5 % is MINE.
Outputs  data/field/2026-10-25/returns/floodplain_multiplicity_log.txt
"""
import os
import numpy as np, pandas as pd, laspy
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."; os.chdir(ROOT)
D = "data/field/2026-10-25/returns"; LOG = open(f"{D}/floodplain_multiplicity_log.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")
S = pd.read_csv(f"{D}/floodplain_return_shape_points.csv")[["group", "point", "E", "N", "canopy_share"]]
P = pd.read_csv("data/field/2026-10-25/rtk_points_all.csv")[["group", "point", "navd88_tied_m"]]; S = S.merge(P, on=["group", "point"]); S["kind"] = np.where(S.canopy_share < 0.05, "sward", "tree edge")
bbox = (S.E.min() - 4, S.E.max() + 4, S.N.min() - 4, S.N.max() + 4); parts = []
for t in ("4342-29-64", "4358-29-01"):
    with laspy.open(f"data/after/valley/{t}_3dep_full.laz") as r:
        for ch in r.chunk_iterator(5_000_000):
            x = np.asarray(ch.x); y = np.asarray(ch.y); m = (x >= bbox[0]) & (x <= bbox[1]) & (y >= bbox[2]) & (y <= bbox[3])
            if t == "4358-29-01": m &= x > 580030.0
            if m.any(): parts.append(pd.DataFrame(dict(x=x[m], y=y[m], z=np.asarray(ch.z)[m], t=np.asarray(ch.gps_time)[m], ps=np.asarray(ch.point_source_id)[m], rn=np.asarray(ch.return_number)[m], nr=np.asarray(ch.number_of_returns)[m], c=np.asarray(ch.classification)[m])))
G = pd.concat(parts, ignore_index=True); say(f"{len(G):,} returns in the points' box; gps_time precision in the file: smallest nonzero step between sorted times {np.diff(np.unique(G.t))[:5].min():.2e} s")
pulses = G.groupby(["ps", "t"]); G["n_seen"] = pulses.rn.transform("size")
say(f"pulses in the box {pulses.ngroups:,}; returns per pulse from the record: " + ", ".join(f"{int(k)}: {100 * v:.1f} %" for k, v in G.nr.value_counts(normalize=True).sort_index().items()) + f"; pulses whose returns all appear in the box {100 * (G.groupby(['ps', 't']).apply(lambda q: len(q) == q.nr.iloc[0]).mean()):.1f} %")
rows = []
for _, p in S.iterrows():
    d = np.hypot(G.x - p.E, G.y - p.N); w = G[d <= 3].copy(); w["dz"] = w.z - p.navd88_tied_m; w["kind"] = p.kind; w["pt"] = f"{p.group}:{p.point}"; rows.append(w)
W = pd.concat(rows, ignore_index=True)
for kind, g in W.groupby("kind"):
    say(f"\n{kind} ({S[S.kind == kind].shape[0]} points): returns {len(g)}; number_of_returns share " + ", ".join(f"{int(k)}: {100 * v:.1f} %" for k, v in g.nr.value_counts(normalize=True).sort_index().items()))
    pul = g.groupby(["pt", "ps", "t"]); say(f"   pulses {pul.ngroups}; pulses with >1 return on record {100 * (pul.nr.first() > 1).mean():.1f} %")
    full = g[g.n_seen == g.nr]; fp = full[full.nr > 1].groupby(["pt", "ps", "t"]); gaps = []
    for _, q in fp:
        q = q.sort_values("rn"); gaps += list(-np.diff(q.z.to_numpy()))
    gaps = np.array(gaps)
    if len(gaps): say(f"   multi-return pulses fully seen: {fp.ngroups}; gaps between consecutive returns: min {gaps.min():.2f} m, p01 {np.percentile(gaps, 1):.2f}, p05 {np.percentile(gaps, 5):.2f}, p10 {np.percentile(gaps, 10):.2f}, median {np.median(gaps):.2f}, p90 {np.percentile(gaps, 90):.2f} m; gaps under 0.5 m {int((gaps < 0.5).sum())}, under 1.0 m {int((gaps < 1.0).sum())}, under 2.0 m {int((gaps < 2.0).sum())}")
    gr = g[g.dz < 0.05]; say(f"   returns within 5 cm of the ground: {len(gr)} ({100 * len(gr) / len(g):.2f} % of returns): single {100 * (gr.nr == 1).mean():.0f} %, last-of-several {100 * ((gr.nr > 1) & (gr.rn == gr.nr)).mean():.0f} %, earlier-of-several {100 * ((gr.nr > 1) & (gr.rn < gr.nr)).mean():.0f} %")
    last = g[(g.nr > 1) & (g.rn == g.nr)]
    if len(last): say(f"   last-of-several returns: {len(last)}, height above ground median {last.dz.median():+.3f} m, p10 {last.dz.quantile(.1):+.3f}; within 5 cm of ground {100 * (last.dz < 0.05).mean():.1f} %")
    say(f"   pulse density: median over points of pulses within 3 m / 28.27 m^2 = {(g.groupby('pt').apply(lambda q: q.groupby(['ps', 't']).ngroups).median() / (np.pi * 9)):.1f} pulses/m^2; returns/m^2 {(g.groupby('pt').size().median() / (np.pi * 9)):.1f}")
LOG.close()
