"""
build_timeline.py
Renders a clean experiment timeline from the Heat_stress Timeline sheet:
temperature history per treatment, phase bands, measurement timepoints and
sampling. Dates transcribed from the sheet.

Changes from previous version:
  - incidents lane and annotations removed
  - 17C photo dates corrected (were 13/19/26/28 May; sheet says 14/19/22/26/28)
  - imaging and IPAM separated into their own lanes, since they do not
    coincide: IPAM was also run on each elevated treatment's ramp-start day
    (8 May TP3, 12 May TP2), which has no matching photo
  - pre-ramp imaging days added (23 Apr - 11 May), previously omitted
"""
import datetime as dt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

D = lambda s: dt.date.fromisoformat(s)
C = {"14C": "#2c7fb8", "17C": "#f0a500", "20C": "#d7301f"}

# --- daily temperature history (from the Timeline sheet) ---
days, base = [], []
d = D("2026-03-23")
ramp = {D("2026-03-23"): 10, D("2026-03-26"): 12, D("2026-03-30"): 12.5,
        D("2026-03-31"): 13, D("2026-04-01"): 13.5, D("2026-04-02"): 14}
cur = 10
while d <= D("2026-05-30"):
    if d in ramp:
        cur = ramp[d]
    days.append(d); base.append(cur); d += dt.timedelta(days=1)


def series(diverge):
    """base until divergence, then follow `diverge` {date: temp}, else hold."""
    out, hold = [], 14
    keys = sorted(diverge)
    for day, b in zip(days, base):
        if keys and day >= keys[0]:
            for k in keys:
                if day >= k:
                    hold = diverge[k]
            out.append(hold)
        else:
            out.append(b)
    return out


t20 = series({D("2026-05-08"): 15, D("2026-05-09"): 16, D("2026-05-10"): 17,
              D("2026-05-11"): 18, D("2026-05-12"): 19, D("2026-05-13"): 20})
t17 = series({D("2026-05-12"): 15, D("2026-05-13"): 16, D("2026-05-14"): 17})
t14 = base

# --- events ---
phases = [("Prep & crossing", "2026-03-23", "2026-03-29", "#dfe6ec"),
          ("Ramp-up & acclimation at 14C", "2026-03-30", "2026-05-07", "#eef3ee"),
          ("Heat stress & measurement", "2026-05-08", "2026-05-30", "#f6e7e4")]

# imaging, all treatments together before the ramps diverge
photo_pre = ["2026-04-23", "2026-04-24", "2026-04-27", "2026-04-30",
             "2026-05-01", "2026-05-04", "2026-05-05", "2026-05-06",
             "2026-05-07", "2026-05-08", "2026-05-11"]
photo = {"20C": ["2026-05-13", "2026-05-18", "2026-05-21", "2026-05-25",
                 "2026-05-27"],
         "17C": ["2026-05-14", "2026-05-19", "2026-05-22", "2026-05-26",
                 "2026-05-28"],
         "14C": ["2026-05-15", "2026-05-20", "2026-05-29"]}
ipam = {"20C": ["2026-05-08", "2026-05-13", "2026-05-18", "2026-05-21",
                "2026-05-25", "2026-05-27"],
        "17C": ["2026-05-12", "2026-05-14", "2026-05-19", "2026-05-22",
                "2026-05-26", "2026-05-28"],
        "14C": ["2026-05-15", "2026-05-20", "2026-05-29"]}
samples = [("DNA", "2026-05-04"), ("DNA", "2026-05-05"),
           ("DNA", "2026-05-06"), ("DNA", "2026-05-07"),
           ("RNA", "2026-05-27"), ("RNA", "2026-05-28"),
           ("RNA", "2026-05-29")]

fig, (ax, ev) = plt.subplots(2, 1, figsize=(13, 7.4), sharex=True,
                             gridspec_kw=dict(height_ratios=[3, 1.6]))

for name, s, e, col in phases:
    ax.axvspan(D(s), D(e), color=col, zorder=0)
    ax.text(D(s) + (D(e) - D(s)) / 2, 20.8, name, ha="center", va="bottom",
            fontsize=9, color="#555")

ax.plot(days, t14, color=C["14C"], lw=2, label="14C")
ax.plot(days, t17, color=C["17C"], lw=2, label="17C")
ax.plot(days, t20, color=C["20C"], lw=2, label="20C")
ax.set_ylabel("temperature (C)")
ax.set_ylim(9, 22)
ax.legend(loc="lower right", fontsize=9)
ax.set_title("Giant kelp heat-stress experiment timeline", fontsize=13,
             fontweight="bold")

# event lanes
lanes = {"Imaging": 2, "IPAM": 1, "Sampling": 0}
for lab, y in lanes.items():
    ev.text(D("2026-03-22"), y, lab, ha="right", va="center", fontsize=8)

ev.scatter([D(x) for x in photo_pre], [2] * len(photo_pre),
           facecolor="white", edgecolor="#777", s=32, lw=0.9, zorder=3)
for t, dates in photo.items():
    ev.scatter([D(x) for x in dates], [2] * len(dates), color=C[t], s=40,
               zorder=3)
for t, dates in ipam.items():
    ev.scatter([D(x) for x in dates], [1] * len(dates), color=C[t], s=40,
               marker="D", zorder=3)

seen = set()
for kind, x in samples:
    m = "s" if kind == "DNA" else "^"
    ev.scatter(D(x), 0, marker=m, color="#444", s=45, zorder=3)
    if kind not in seen:
        ev.text(D(x), 0.32, kind, fontsize=6.5, ha="left", color="#444")
        seen.add(kind)

ev.text(D("2026-04-23"), 2.42, "pre-ramp imaging (all treatments)",
        fontsize=6.5, color="#777", ha="left")

ev.set_ylim(-0.6, 2.8)
ev.set_yticks([])
for s in ("left", "right", "top"):
    ev.spines[s].set_visible(False)
ev.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MO))
ev.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
plt.setp(ev.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ev.set_xlim(D("2026-03-20"), D("2026-06-01"))

fig.tight_layout()
fig.savefig("experiment_timeline.png", dpi=150, bbox_inches="tight")
print("wrote experiment_timeline.png")