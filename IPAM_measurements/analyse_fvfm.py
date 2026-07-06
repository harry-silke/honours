"""
analyse_fvfm.py
---------------
Fv/Fm heat-resistance analysis from master_fvfm.csv.

Days-in-treatment = measurement date minus ramp-end (20C 13 May, 17C 14 May,
14C 15 May). Pre-ramp baseline reads (negative days) are dropped. Known
mechanical deaths in EXCLUDE_IDS are removed. Low-tissue and FvFm<=0 rows
are dropped.

Outputs (into --out, default figures_fvfm/):
  fvfm_timecourse_by_population.png  one panel per population, a line per
                                     treatment over days, SD error bars
  fvfm_resistance_ranking.png        populations ranked by endpoint resistance
  resistance_ranking.csv             the same numbers

Resistance = each cross's endpoint (max day) 20C Fv/Fm minus its own 14C
sibling, aggregated by female population. Higher = holds up better vs its
own control.
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAMP_END = {"14C": "2026-05-15", "17C": "2026-05-14", "20C": "2026-05-13"}
EXCLUDE_IDS = {"52C"}          # mechanical death (waterbath 6 freshwater leak)
TREAT_ORDER = ["14C", "17C", "20C"]
TREAT_COLOUR = {"14C": "#2c7fb8", "17C": "#f0a500", "20C": "#d7301f"}


def load(path):
    d = pd.read_csv(path)
    d = d[(~d["low_tissue"].astype(bool)) & (d["FvFm"] > 0)
          & (~d["individual"].isin(EXCLUDE_IDS))].copy()
    d["dt"] = pd.to_datetime(d["date_iso"])
    d["day"] = (d["dt"] - pd.to_datetime(d["treatment"].map(RAMP_END))).dt.days
    d = d[d["day"] >= 0]
    return d


def timecourse(d, outdir):
    pops = sorted(d["female_pop"].dropna().unique())
    ncol = 3
    nrow = int(np.ceil(len(pops) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3 * nrow),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, pop in zip(axes, pops):
        sub = d[d["female_pop"] == pop]
        for t in TREAT_ORDER:
            g = (sub[sub["treatment"] == t].groupby("day")["FvFm"]
                 .agg(["mean", "std"]).reset_index())
            if g.empty:
                continue
            ax.errorbar(g["day"], g["mean"], yerr=g["std"], marker="o",
                        ms=3, capsize=2, lw=1.2, color=TREAT_COLOUR[t], label=t)
        ax.set_title(pop, fontsize=9)
        ax.set_ylim(0, 0.8)
    for ax in axes[len(pops):]:
        ax.set_visible(False)
    axes[0].legend(title="treatment", fontsize=8)
    fig.supxlabel("days in treatment")
    fig.supylabel("Fv/Fm (mean +/- SD)")
    fig.tight_layout()
    p = os.path.join(outdir, "fvfm_timecourse_by_population.png")
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


def resistance(d, outdir):
    end = (d.sort_values("day")
           .groupby(["cross_id", "female_pop", "treatment"]).tail(1))
    piv = end.pivot_table(index=["cross_id", "female_pop"],
                          columns="treatment", values="FvFm").reset_index()

    # baseline A: each cross's own day-0 14C reading (fallback: pop mean)
    base14 = d[(d.treatment == "14C") & (d.day == d.day.min())]
    base_by_cross = base14.groupby("cross_id")["FvFm"].mean()
    base_by_pop = base14.groupby("female_pop")["FvFm"].mean()
    piv["day0_baseline"] = piv["cross_id"].map(base_by_cross)
    missing = piv["day0_baseline"].isna()
    if missing.any():
        piv.loc[missing, "day0_baseline"] = piv.loc[
            missing, "female_pop"].map(base_by_pop)
        print(f"{int(missing.sum())} cross(es) had no day-0 14C reading; "
              f"used their population's day-0 mean instead.")

    # baseline B: that same cross's own endpoint (day-14) 14C reading
    for t in ("17C", "20C"):
        if t in piv:
            piv[f"vsDay0_{t}"] = piv[t] - piv["day0_baseline"]
            if "14C" in piv:
                piv[f"vsEndpoint_{t}"] = piv[t] - piv["14C"]

    def agg(t, kind, with_n):
        col = f"{kind}_{t}"
        d = {f"{kind}_{t}": (col, "median"),
             f"{kind}_{t}_se": (col, lambda s: s.std() / np.sqrt(len(s)))}
        if with_n:
            d[f"n_{t}"] = (col, "size")
        return piv.dropna(subset=[col]).groupby("female_pop").agg(**d).round(3)

    rank = (agg("20C", "vsDay0", True).join(agg("17C", "vsDay0", True))
            .join(agg("20C", "vsEndpoint", False))
            .join(agg("17C", "vsEndpoint", False)))
    rank = rank.sort_values("vsDay0_20C", ascending=False)
    rank.to_csv(os.path.join(outdir, "resistance_ranking.csv"))

    order = rank.index[::-1]
    y = np.arange(len(order))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)

    def panel(ax, kind, title, xlabel):
        b17 = ax.barh(y + 0.2, rank.loc[order, f"{kind}_17C"], height=0.4,
                      xerr=rank.loc[order, f"{kind}_17C_se"], capsize=3,
                      color=TREAT_COLOUR["17C"], label="17C")
        b20 = ax.barh(y - 0.2, rank.loc[order, f"{kind}_20C"], height=0.4,
                      xerr=rank.loc[order, f"{kind}_20C_se"], capsize=3,
                      color=TREAT_COLOUR["20C"], label="20C")
        for yi, pop in zip(y, order):
            n20 = int(rank.loc[pop, "n_20C"])
            n17 = int(rank.loc[pop, "n_17C"])
            ax.text(ax.get_xlim()[1], yi, f"  n={n17}/{n20}",
                    fontsize=7, va="center", ha="left", color="#555")
        ax.axvline(0, color="k", lw=.8)
        ax.set_yticks(y)
        ax.set_yticklabels(order)
        ax.set_xlabel(xlabel)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8, loc="lower right" if kind == "vsDay0" else
                  "upper left")

    panel(ax1, "vsDay0",
          "Comparison A: last day vs FIRST day\n"
          "(total change over the whole experiment, own 14C day-0 baseline)",
          "Fv/Fm: last day minus first-day 14C baseline")
    panel(ax2, "vsEndpoint",
          "Comparison B: last day vs LAST day\n"
          "(heat effect only, own 14C control measured on the same "
          "final day)",
          "Fv/Fm: last day minus same-day 14C control")

    fig.suptitle("Heat resistance by female population: two ways to define "
                 "the control", fontsize=12)
    fig.tight_layout()
    p = os.path.join(outdir, "fvfm_resistance_ranking.png")
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return rank, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("master", nargs="?", default="master_fvfm.csv")
    ap.add_argument("--out", default="figures_fvfm")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = load(args.master)
    tp = timecourse(d, args.out)
    rank, rp = resistance(d, args.out)
    print("Heat resistance by population (endpoint 20C - own 14C):")
    print(rank.to_string())
    print(f"\nWrote {tp}\n      {rp}\n      "
          f"{os.path.join(args.out, 'resistance_ranking.csv')}")


if __name__ == "__main__":
    main()