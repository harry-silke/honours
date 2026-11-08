"""
build_timeline_simple.py
Experiment timeline, sized to sit as a panel beneath the BioRender workflow
and the crossing design figure.

Changes from the first version:
  - the three treatments are identical until 8 May, so the shared period is
    drawn as one grey line rather than three overlapping coloured ones;
    colours diverge only where the treatments actually diverge
  - measurement days marked during the acclimation and growth period as well
    as during heat stress
  - larger fonts for print
  - panel label (set PANEL = None to omit)

Measurement days are the imaging sessions that have data in
master_photos_with_parents.csv, not the eleven listed in the Timeline sheet:
a dot should mean a measurement that can be analysed. 30 April, 4, 5 and 6 May
fall while every treatment was still at 14 C, so they sit on the shared grey
line. 8 and 11 May fall after the 20 C ramp began, so they are drawn on each
treatment's own line at that treatment's temperature.
"""
import datetime as dt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

PANEL = "c"           # set to None to omit the panel label
SHARED = "#6b6a65"    # grey for the period before treatments diverge

D = lambda s: dt.date.fromisoformat(s)
C = {"14C": "#2c7fb8", "17C": "#f0a500", "20C": "#d7301f"}
SPLIT = D("2026-05-08")   # first day any treatment leaves 14 C

# --- daily temperature history ---
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


temps = {
    "20C": series({D("2026-05-08"): 15, D("2026-05-09"): 16,
                   D("2026-05-10"): 17, D("2026-05-11"): 18,
                   D("2026-05-12"): 19, D("2026-05-13"): 20}),
    "17C": series({D("2026-05-12"): 15, D("2026-05-13"): 16,
                   D("2026-05-14"): 17}),
    "14C": base,
}

# heat-stress measurement days, aligned to days since each treatment reached
# its target temperature
measure = {"20C": ["2026-05-13", "2026-05-18", "2026-05-21", "2026-05-25",
                   "2026-05-27"],
           "17C": ["2026-05-14", "2026-05-19", "2026-05-22", "2026-05-26",
                   "2026-05-28"],
           "14C": ["2026-05-15", "2026-05-20", "2026-05-29"]}

# acclimation and growth period. pre_shared sits on the grey line; pre_split
# falls after the 20 C ramp started, so it goes on the treatment lines.
pre_shared = ["2026-04-30", "2026-05-04", "2026-05-05", "2026-05-06"]
pre_split  = {"14C": ["2026-05-08"],
              "17C": ["2026-05-08", "2026-05-11"],
              "20C": ["2026-05-08", "2026-05-11"]}
for t, extra in pre_split.items():
    measure[t] = extra + measure[t]

temp_on   = {t: dict(zip(days, v)) for t, v in temps.items()}
shared_on = dict(zip(days, base))
i_split   = days.index(SPLIT)

fig, ax = plt.subplots(figsize=(11, 3.6))

ax.axvspan(D("2026-03-23"), SPLIT, color="#f0f3f5", zorder=0)
ax.axvspan(SPLIT, D("2026-05-30"), color="#faeeeb", zorder=0)
ax.text(D("2026-04-15"), 21.2, "Acclimation and growth", ha="center",
        fontsize=11, color="#555")
ax.text(D("2026-05-19"), 21.2, "Heat stress", ha="center", fontsize=11,
        color="#555")

# shared history: one grey line, all treatments identical here
ax.plot(days[:i_split + 1], base[:i_split + 1], color=SHARED, lw=2.2,
        zorder=2, label="all treatments")

# measurement days during acclimation and growth, on the shared line
pts = [D(x) for x in pre_shared]
ax.scatter(pts, [shared_on[p] for p in pts], color=SHARED, s=34,
           edgecolor="white", lw=1.0, zorder=3)

# diverging period only
for t in ("14C", "17C", "20C"):
    ax.plot(days[i_split - 1:], temps[t][i_split - 1:], color=C[t], lw=2.2,
            zorder=2, label=t)
    pts = [D(x) for x in measure[t]]
    ax.scatter(pts, [temp_on[t][p] for p in pts], color=C[t], s=34,
               edgecolor="white", lw=1.0, zorder=3)

ax.set_ylabel("Temperature (\u00b0C)", fontsize=12)
ax.set_ylim(9, 22)
ax.legend(loc="center left", fontsize=10.5, frameon=False)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
ax.xaxis.set_minor_locator(mdates.WeekdayLocator(byweekday=mdates.MO))
ax.set_xlim(D("2026-03-20"), D("2026-06-01"))
ax.tick_params(labelsize=11)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

if PANEL:
    ax.text(-0.085, 1.10, f"({PANEL})", transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="top", ha="left")

fig.tight_layout()
fig.savefig("experiment_timeline_simple.png", dpi=300, bbox_inches="tight",
            facecolor="white")
print("wrote experiment_timeline_simple.png")