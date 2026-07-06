"""
build_timeline.py
Renders a clean experiment timeline from the Heat_stress Timeline sheet:
temperature history per treatment, phase bands, measurement timepoints,
sampling, and logged incidents. Dates transcribed from the sheet.
"""
import datetime as dt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Patch

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
          ("Heat challenge & measurement", "2026-05-08", "2026-05-30", "#f6e7e4")]
measure = {"20C": ["2026-05-13", "2026-05-18", "2026-05-21", "2026-05-25", "2026-05-27"],
           "17C": ["2026-05-13", "2026-05-19", "2026-05-26", "2026-05-28"],
           "14C": ["2026-05-15", "2026-05-20", "2026-05-29"]}
samples = [("DNA", "2026-05-04"), ("DNA", "2026-05-05"),
           ("DNA", "2026-05-06"), ("DNA", "2026-05-07"),
           ("RNA", "2026-05-27"), ("RNA", "2026-05-28"),
           ("RNA", "2026-05-29")]
incidents = [("2026-05-07", "331B/348B\ncross-contam."),
             ("2026-05-13", "WC skipped\nRack 3B"),
             ("2026-05-17", "52C died\n(waterbath 6 leak)"),
             ("2026-05-26", "50% WC\ninstead of 90%")]

fig, (ax, ev) = plt.subplots(2, 1, figsize=(13, 7), sharex=True,
                             gridspec_kw=dict(height_ratios=[3, 1.4]))

for name, s, e, col in phases:
    ax.axvspan(D(s), D(e), color=col, zorder=0)
    ax.text(D(s) + (D(e) - D(s)) / 2, 20.8, name, ha="center", va="bottom",
            fontsize=9, color="#555")

ax.plot(days, t14, color=C["14C"], lw=2, label="14C (TP1)")
ax.plot(days, t17, color=C["17C"], lw=2, label="17C (TP2)")
ax.plot(days, t20, color=C["20C"], lw=2, label="20C (TP3)")
ax.set_ylabel("temperature (C)")
ax.set_ylim(9, 22)
ax.legend(loc="lower right", fontsize=9)
ax.set_title("Giant kelp heat-stress experiment timeline", fontsize=13,
             fontweight="bold")

for x, lab in incidents:
    ax.axvline(D(x), color="#b3001b", ls=":", lw=1, zorder=1)
    ax.annotate(lab, (D(x), 9.4), fontsize=7, color="#b3001b", ha="center",
                va="bottom")

# event lanes
lanes = {"IPAM + photo": 2, "Sampling": 1, "Incidents": 0}
for lab, y in lanes.items():
    ev.text(D("2026-03-22"), y, lab, ha="right", va="center", fontsize=8)
for t, dates in measure.items():
    ev.scatter([D(x) for x in dates], [2] * len(dates), color=C[t], s=40,
               zorder=3)
seen = set()
for kind, x in samples:
    m = "s" if kind == "DNA" else "^"
    ev.scatter(D(x), 1, marker=m, color="#444", s=45, zorder=3)
    if kind not in seen:
        ev.text(D(x), 1.32, kind, fontsize=6.5, ha="left", color="#444")
        seen.add(kind)
for x, _ in incidents:
    ev.scatter(D(x), 0, marker="x", color="#b3001b", s=55, zorder=3)

ev.set_ylim(-0.6, 2.7)
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