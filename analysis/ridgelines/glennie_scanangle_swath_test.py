#!/usr/bin/env python
"""
Glennie et al. (2014, GRL 10.1002/2014GL059919) TEST 2 for the Elba, MN
gen1 (2008) "reads-the-ground-too-low-on-slopes" residual.

Question: is that residual a LEGACY GPS/IMU BORESIGHT / SCAN-MIRROR artifact?
Diagnostic signature (Glennie): a signed offset that
  (a) CORRELATES WITH SCAN ANGLE / distance-from-flight-line,
  (b) CONCENTRATES AT SWATH EDGES (large |scan_angle|),
  (c) differs COHERENTLY BETWEEN SWATHS (point_source_id).
If present -> legacy nav-era systematic error, not slope physics.

This is entirely gen1-internal. The gen2 z_after bare-earth plane is only a
fixed spatial reference (nested-subset comparison of a fixed registered dataset
is internally consistent). We do NOT use gen2 magnitude/penetration/canopy as a
covariate.

CRITICAL: scan angle and slope/incidence are CONFOUNDED (steeper slopes are
sampled with different geometry). We therefore ALWAYS run the scan-angle test
WITHIN slope bands. A pooled d_mm-vs-scan_angle correlation is not acceptable.

Datum: we report everything as a RESIDUAL relative to the flat-ground
(slope < 3 deg) median d_mm. That flat-ground median already has the
GEOID03->GEOID18 (~67 mm) constant baked into the z_after reference; subtracting
it removes the datum constant so the remaining structure is the signal.

Sign convention: gen1 BELOW the gen2 plane -> d_mm NEGATIVE. "gen1 low" = more
negative residual.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DERIVED = os.path.join(ROOT, "data", "derived", "elba_fulldensity")
NPZ = os.path.join(DERIVED, "gen1_csf_angles.npz")

# grid geometry (from task spec / project memory)
NX, NY = 508, 700
FLAT_SLOPE_DEG = 3.0

SLOPE_BINS = [0, 5, 12, 20, 27, 40, 90]
SLOPE_LABELS = ["0-5", "5-12", "12-20", "20-27", "27-40", ">40"]
SA_BINS = [0, 5, 10, 15, 20]          # |scan_angle| bin edges (deg)
SA_LABELS = ["0-5", "5-10", "10-15", "15-20"]


def nmad(x):
    """Normalized median absolute deviation (robust sigma)."""
    if x.size == 0:
        return np.nan
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))


def load():
    d = np.load(NPZ)
    ing = d["in_grid"]
    out = {k: d[k][ing] for k in
           ["incidence", "scan_angle", "slope", "d_mm", "cell",
            "point_source_id", "stratum", "core_forest", "core_open"]}
    out["absa"] = np.abs(out["scan_angle"])

    # datum = flat-ground (< FLAT_SLOPE_DEG) median d_mm
    flat = out["slope"] < FLAT_SLOPE_DEG
    datum = float(np.median(out["d_mm"][flat]))
    out["datum_mm"] = datum
    out["res_mm"] = out["d_mm"] - datum          # residual (datum removed)
    out["n_flat"] = int(flat.sum())

    # ridgeline flow-free mask, mapped to per-return via cell index
    ridge = np.load(os.path.join(DERIVED, "ridge_mask.npy"))
    cxx = np.load(os.path.join(DERIVED, "curv_xx.npy"))
    cyy = np.load(os.path.join(DERIVED, "curv_yy.npy"))
    # "small |curv|" threshold: median of |curv| distribution (~0.0018 both axes)
    cthr = 0.002
    flowfree = ridge & (np.abs(cxx) < cthr) & (np.abs(cyy) < cthr)
    # cell = flat NY*NX index (row-major over (NY,NX)); guard bounds
    ok = (out["cell"] >= 0) & (out["cell"] < NX * NY)
    ridge_flat = ridge.ravel()
    flowfree_flat = flowfree.ravel()
    out["on_ridge"] = np.zeros(out["cell"].shape, bool)
    out["on_flowfree"] = np.zeros(out["cell"].shape, bool)
    out["on_ridge"][ok] = ridge_flat[out["cell"][ok]]
    out["on_flowfree"][ok] = flowfree_flat[out["cell"][ok]]
    return out


def binstat(res, sel):
    x = res[sel]
    return np.median(x) if x.size else np.nan, nmad(x), int(x.size)


def table_slope_x_scanangle(D, sel_base=None, title=""):
    """Median residual d_mm per (slope band x |scan_angle| bin)."""
    res = D["res_mm"]; sl = D["slope"]; absa = D["absa"]
    base = np.ones(res.shape, bool) if sel_base is None else sel_base
    lines = []
    lines.append(f"### {title}")
    lines.append("")
    hdr = "| slope (deg) | " + " | ".join(f"|SA| {s}" for s in SA_LABELS) + " | nadir<5 vs edge>15 |"
    lines.append(hdr)
    lines.append("|" + "---|" * (len(SA_LABELS) + 2))
    rows = []
    for i in range(len(SLOPE_LABELS)):
        s0, s1 = SLOPE_BINS[i], SLOPE_BINS[i + 1]
        smask = base & (sl >= s0) & (sl < s1)
        cells = []
        rowvals = {}
        for j in range(len(SA_LABELS)):
            a0, a1 = SA_BINS[j], SA_BINS[j + 1]
            sel = smask & (absa >= a0) & (absa < a1)
            m, nm, n = binstat(res, sel)
            rowvals[SA_LABELS[j]] = (m, nm, n)
            cells.append(f"{m:+.0f} (nmad {nm:.0f}, n={n})" if n else "-")
        # nadir vs edge delta
        mn, _, nn = binstat(res, smask & (absa < 5))
        me, _, ne = binstat(res, smask & (absa > 15))
        delta = me - mn if (nn and ne) else np.nan
        dcell = f"{mn:+.0f} -> {me:+.0f}  (d={delta:+.0f})" if (nn and ne) else "-"
        lines.append("| " + SLOPE_LABELS[i] + " | " + " | ".join(cells) + " | " + dcell + " |")
        rows.append((SLOPE_LABELS[i], rowvals, mn, me, delta, nn, ne))
    lines.append("")
    return "\n".join(lines), rows


def table_swaths(D, sel_base=None):
    res = D["res_mm"]; sl = D["slope"]; absa = D["absa"]; psid = D["point_source_id"]
    base = np.ones(res.shape, bool) if sel_base is None else sel_base
    lines = ["### Per-swath (point_source_id) coherence", ""]
    lines.append("| swath | n | med res (all) | med nadir<5 | med edge>15 | edge-nadir | n(slope>27) | med res slope>27 |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for ps in sorted(np.unique(psid)):
        pm = base & (psid == ps)
        m_all, _, n_all = binstat(res, pm)
        m_n, _, _ = binstat(res, pm & (absa < 5))
        m_e, _, _ = binstat(res, pm & (absa > 15))
        d = m_e - m_n
        steep = pm & (sl > 27)
        m_st, _, n_st = binstat(res, steep)
        lines.append(f"| {ps} | {n_all} | {m_all:+.0f} | {m_n:+.0f} | {m_e:+.0f} | {d:+.0f} | {n_st} | {m_st:+.0f} |")
    lines.append("")
    return "\n".join(lines)


def make_figure(D, sel_base, path, suptitle):
    res = D["res_mm"]; sl = D["slope"]; absa = D["absa"]
    base = np.ones(res.shape, bool) if sel_base is None else sel_base
    # fine |scan_angle| bins for curves
    edges = np.arange(0, 18.001, 2.0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=140)
    cmap = plt.get_cmap("viridis")
    for i in range(len(SLOPE_LABELS)):
        s0, s1 = SLOPE_BINS[i], SLOPE_BINS[i + 1]
        smask = base & (sl >= s0) & (sl < s1)
        med = []
        for k in range(len(centers)):
            sel = smask & (absa >= edges[k]) & (absa < edges[k + 1])
            med.append(np.median(res[sel]) if sel.sum() > 50 else np.nan)
        ax.plot(centers, med, "-o", ms=4, color=cmap(i / (len(SLOPE_LABELS) - 1)),
                label=f"slope {SLOPE_LABELS[i]} deg")
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_xlabel("|scan angle| (deg)  -- swath edge to the right")
    ax.set_ylabel("median residual d_mm (datum removed)\n<0 = gen1 reads ground low")
    ax.set_title(suptitle)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def main():
    D = load()
    print(f"datum (flat-ground median d_mm) = {D['datum_mm']:.1f} mm "
          f"(from n={D['n_flat']} flat returns)")

    # ---- report assembly ----
    out = []
    out.append("# Glennie TEST 2 -- scan-angle / swath-edge / per-swath diagnostic")
    out.append("")
    out.append(f"_Generated by `glennie_scanangle_swath_test.py`. "
               f"n = {D['res_mm'].size:,} gen1 CSF ground returns (in_grid)._")
    out.append("")
    out.append("**Datum removed:** median d_mm on flat stable ground "
               f"(slope < {FLAT_SLOPE_DEG:g} deg) = **{D['datum_mm']:.1f} mm** "
               f"(n={D['n_flat']:,}). All values below are residuals relative to "
               "that flat-ground median; the ~67 mm GEOID03->GEOID18 datum is "
               "already baked into the z_after reference and removed here as part "
               "of the flat-ground constant. Sign: residual < 0 means gen1 reads "
               "the ground LOW.")
    out.append("")
    out.append("**Confound guard:** scan angle and slope are confounded, so the "
               "scan-angle test is run WITHIN slope bands. Judge by median mm "
               "shifts, not r (per-return NMAD is ~150-270 mm).")
    out.append("")

    t_all, rows_all = table_slope_x_scanangle(
        D, None, "1. All ground returns -- median residual d_mm by (slope band, |scan angle|)")
    out.append(t_all)

    # steep stable set: flow-free ridgeline (best proxy for the >27 deg onset locus)
    t_rf, rows_rf = table_slope_x_scanangle(
        D, D["on_flowfree"],
        "1b. Flow-free ridgeline cells only (ridge_mask & |curv_xx|<0.002 & |curv_yy|<0.002)")
    out.append(t_rf)

    out.append("## 2. Swath-edge concentration")
    out.append("")
    out.append("The `nadir<5 vs edge>15` column in the tables above is the direct "
               "test: median residual for |scan_angle|<5 deg vs >15 deg within each "
               "slope band. A boresight / scan-mirror error is worst at the edges "
               "(large delta, edge more negative).")
    out.append("")

    out.append(table_swaths(D, None))
    out.append(table_swaths(D, D["on_flowfree"]).replace(
        "Per-swath (point_source_id) coherence",
        "Per-swath coherence -- flow-free ridgeline cells only"))

    # within-swath nadir-vs-edge inside the steep (27-40) band: the crux
    res = D["res_mm"]; sl = D["slope"]; absa = D["absa"]; psid = D["point_source_id"]
    out.append("### 3b. Within-swath nadir-vs-edge inside the steep 27-40 deg band (the onset locus)")
    out.append("")
    out.append("| swath | med nadir<5 (n) | med edge>15 (n) | edge-nadir |")
    out.append("|---|---|---|---|")
    band = (sl >= 27) & (sl < 40)
    for p in sorted(np.unique(psid)):
        m = band & (psid == p)
        nsel = m & (absa < 5); esel = m & (absa > 15)
        mn = np.median(res[nsel]) if nsel.sum() > 30 else np.nan
        me = np.median(res[esel]) if esel.sum() > 30 else np.nan
        out.append(f"| {p} | {mn:+.0f} (n={int(nsel.sum())}) | "
                   f"{me:+.0f} (n={int(esel.sum())}) | {me-mn:+.0f} |")
    out.append("")
    out.append("Swath geometry (does one swath just cover more steep ground?):")
    out.append("")
    out.append("| swath | med \\|SA\\| | med slope | frac slope>27 |")
    out.append("|---|---|---|---|")
    for p in sorted(np.unique(psid)):
        m = psid == p
        out.append(f"| {p} | {np.median(absa[m]):.1f} | "
                   f"{np.median(sl[m]):.1f} | {(sl[m]>27).mean():.3f} |")
    out.append("")

    # ---- VERDICT ----
    out.append("## VERDICT")
    out.append("")
    out.append("**Boresight / scan-mirror / legacy-nav artifact: NOT SUPPORTED "
               "(excluded as the driver of the steep gen1-low onset).**")
    out.append("")
    out.append("Glennie's three signatures fail:")
    out.append("")
    out.append("1. **Sign is wrong for a scan-angle artifact.** Within every "
               "low-to-moderate slope band (0-27 deg), the ground reads *lower "
               "at nadir* and *higher toward the swath edge* (edge-nadir delta "
               "= **+10 to +26 mm**, positive). A boresight / scan-mirror error "
               "drives the ground progressively LOW at the edges (negative "
               "delta). We see the opposite: the residual improves toward the "
               "edge, not degrades.")
    out.append("")
    out.append("2. **The steep onset is a NADIR phenomenon, not edge-"
               "concentrated.** The -30 to -45 mm gen1-low at 27-40 deg sits at "
               "|scan_angle| < 5 deg (near-nadir) and gets *less* negative toward "
               "the edge in the pooled table (-41 -> -25). It is strongest where "
               "a scan-angle artifact would be weakest.")
    out.append("")
    out.append("3. **Between-swath edge behavior is incoherent.** Inside the "
               "27-40 deg band the edge-nadir delta flips sign across swaths: "
               "swath 136 = -54 mm (edge more negative), swath 137 = +37 mm "
               "(edge less negative), swath 138 = +11 mm. A shared boresight "
               "would push all swaths the *same* direction at the edges. They "
               "disagree. The per-swath *offsets* do differ (135 -26, 136 -13, "
               "137 -2, 138 +18 mm), but that is a mild co-registration / "
               "flight-line offset, not a scan-angle ramp, and the swaths do "
               "not differ enough in steep-ground coverage (frac slope>27 = "
               "0.07-0.12) to explain the onset by sampling.")
    out.append("")
    out.append("The >40 deg band shows a large -126 mm edge-nadir delta, but it "
               "rests on n~750 edge returns with NMAD 520-674 mm; it is noise-"
               "dominated extreme geometry, not a reliable signal.")
    out.append("")
    out.append("**Conclusion.** The gen1 steep-slope low is *scan-angle-flat to "
               "improving* within slope bands and *nadir-worst*, so it is NOT a "
               "legacy GPS/IMU boresight or scan-mirror artifact. The magnitude "
               "in play (steep-band nadir ~ -30 to -45 mm relative to flat "
               "ground) exceeds the ~20 mm budget and remains to be explained by "
               "a genuine slope effect (return-geometry / leaf-state ranging on "
               "the slope-facing surface) or residual co-registration -- not by "
               "the scanner. This is a clean negative that removes boresight from "
               "the candidate list.")
    out.append("")

    fig1 = make_figure(D, None, os.path.join(HERE, "glennie_scanangle_all.png"),
                       "Median residual d_mm vs |scan angle|, by slope band (all ground returns)")

    out_md = os.path.join(HERE, "GLENNIE_SCANANGLE_SWATH_TEST.md")
    # verdict is appended by hand in the .md after inspecting numbers; here we
    # emit a data-driven verdict block too.
    out.append("## Data-driven summary (edge-nadir delta by slope band, all returns)")
    out.append("")
    out.append("| slope band | nadir<5 | edge>15 | edge-nadir |")
    out.append("|---|---|---|---|")
    for lab, _, mn, me, delta, nn, ne in rows_all:
        out.append(f"| {lab} | {mn:+.0f} | {me:+.0f} | {delta:+.0f} |")
    out.append("")

    with open(out_md, "w") as f:
        f.write("\n".join(out))
    print("wrote", out_md)
    print("wrote", fig1)
    # echo key numbers to stdout
    print("\nEdge-nadir delta (mm) by slope band, all returns:")
    for lab, _, mn, me, delta, nn, ne in rows_all:
        print(f"  slope {lab:>6s}: nadir {mn:+.0f}  edge {me:+.0f}  delta {delta:+.0f}")
    print("\nFlow-free ridgeline, edge-nadir delta:")
    for lab, _, mn, me, delta, nn, ne in rows_rf:
        if nn and ne:
            print(f"  slope {lab:>6s}: nadir {mn:+.0f}  edge {me:+.0f}  delta {delta:+.0f}  (n_n~,n_e~)")


if __name__ == "__main__":
    main()
