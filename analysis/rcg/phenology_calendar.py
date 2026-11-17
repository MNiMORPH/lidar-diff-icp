"""Phenology calendar for the plants of the Whitewater floodplain and bluff, from the USA-NPN Nature's Notebook
observations the literature agent pulled (analysis/rcg/usanpn/usanpn_earliest_first_yes_per_site_year.csv: the
earliest first-yes day of year per species, phenophase, site and year, MN/WI/IA, 2009-2026), with the two lidar
flights, the RTK survey and the Sentinel-2 scenes used in data/field/2026-10-25 drawn as dates. One point per
site-year; the box is the interquartile range, the line the median; n beside each row. No numbers typed by hand
except the dates of our own acquisitions. Output analysis/rcg/phenology_calendar_2026-11-18.png and a CSV of the
medians."""
import pandas as pd, numpy as np, datetime as dt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = pd.read_csv("analysis/rcg/usanpn/usanpn_earliest_first_yes_per_site_year.csv")
d = d[(d.latitude >= 42.5) & (d.latitude <= 47.5)]            # the band of latitudes around the site (MINE: +-2.5 deg of 44.1 N)
species = [("reed canarygrass", "reed canary grass"), ("smooth brome", "smooth brome"), ("Kentucky bluegrass", "Kentucky bluegrass"), ("big bluestem", "big bluestem"), ("common buckthorn", "common buckthorn")]
phases = [("Initial growth (grasses/sedges)", "first growth"), ("Breaking leaf buds", "first growth"), ("Open flowers (grasses/sedges)", "open flowers"), ("Open flowers", "open flowers"),
          ("Ripe fruits", "ripe fruits"), ("Colored leaves", "colored leaves"), ("Falling leaves", "falling leaves"), ("All leaves withered (herbaceous and grasses)", "leaves withered")]
order = ["first growth", "open flowers", "ripe fruits", "colored leaves", "falling leaves", "leaves withered"]
colors = {"first growth": "#2a9d8f", "open flowers": "#e9c46a", "ripe fruits": "#f4a261", "colored leaves": "#e76f51", "falling leaves": "#9c6644", "leaves withered": "#6c757d"}
rows = []; fig, ax = plt.subplots(figsize=(11, 5.2)); y = 0; yt = []; yl = []
for key, name in species:
    s = d[d.common_name == key]
    for ph, lab in phases:
        v = s[s.phenophase == ph].earliest_first_yes_doy.values
        if len(v) < 3: continue
        q1, med, q3 = np.percentile(v, [25, 50, 75]); rows.append(dict(species=name, phenophase=lab, n=len(v), doy_p25=q1, doy_median=med, doy_p75=q3, doy_min=v.min(), doy_max=v.max()))
        ax.scatter(v, np.full(len(v), y) + np.random.default_rng(0).uniform(-0.18, 0.18, len(v)), s=9, color=colors[lab], alpha=0.55, linewidths=0)
        ax.plot([q1, q3], [y, y], color=colors[lab], lw=6, solid_capstyle="butt", alpha=0.9); ax.plot([med, med], [y - 0.3, y + 0.3], color="k", lw=1.2)
        ax.text(372, y, f"n {len(v)}", va="center", fontsize=7, color="#444"); yt.append(y); yl.append(f"{name}: {lab}"); y += 1
    y += 0.6
for doy, lab, c in ((dt.date(2021, 5, 1).timetuple().tm_yday, "2021 lidar (1 May)", "#1f77b4"), (dt.date(2026, 10, 25).timetuple().tm_yday, "RTK survey (25 Oct)", "k")):
    ax.axvline(doy, color=c, lw=1.2, ls="--"); ax.text(doy, -0.8, lab, rotation=90, va="bottom", ha="center", fontsize=7.5, color=c)
ax.axvspan(dt.date(2008, 11, 1).timetuple().tm_yday, dt.date(2008, 11, 30).timetuple().tm_yday, color="#1f77b4", alpha=0.10); ax.text(dt.date(2008, 11, 15).timetuple().tm_yday, -0.8, "2008 lidar (November)", rotation=90, va="bottom", ha="center", fontsize=7.5, color="#1f77b4")
for dd in ("2026-04-10", "2026-07-14", "2026-08-08", "2026-10-17", "2026-11-08"):
    ax.axvline(dt.date.fromisoformat(dd).timetuple().tm_yday, color="#aaa", lw=0.6)
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=8); ax.set_ylim(y - 0.4, -0.6); ax.set_xlim(60, 400)
months = [dt.date(2021, m, 1).timetuple().tm_yday for m in range(3, 13)]; ax.set_xticks(months); ax.set_xticklabels(["Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
ax.set_xlabel("day of year (USA-NPN first observation per site and year, 42.5-47.5 N; thin gray lines: Sentinel-2 scenes used)")
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.savefig("analysis/rcg/phenology_calendar_2026-11-18.png", dpi=150, bbox_inches="tight"); pd.DataFrame(rows).to_csv("analysis/rcg/phenology_calendar_medians.csv", index=False)
print(pd.DataFrame(rows).round(0).to_string(index=False))
