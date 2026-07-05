"""
analyse_fvfm_vs_photos.py
-------------------------
Cross-compare PAM Fv/Fm against the photo metrics (coverage, area,
reflectance) per individual at the endpoint (~day 14 in treatment).

The two datasets were shot on different days and different date formats, so
they are matched on individual, each reduced to its last day-in-treatment
reading, then joined. Days-in-treatment uses ramp-end 20C 13 May, 17C 14 May,
14C 15 May. Reflectance is the mean of the three calibrated channels.

Usage:
    python analyse_fvfm_vs_photos.py master_fvfm.csv master_results.csv \
        --out figures_cross/
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAMP_END = {"14C": "2026-05-15", "17C": "2026-05-14", "20C": "2026-05-13"}
SUF_TREAT = {"": "14C", "A": "14C", "B": "17C", "C": "20C"}
TREAT_COLOUR = {"14C": "#2c7fb8", "17C": "#f0a500", "20C": "#d7301f"}


def days(df, date_col, fmt):
    dt = pd.to_datetime(df[date_col].astype(str), format=fmt, errors="coerce")
    ramp = pd.to_datetime(df["treatment"].map(RAMP_END))
    return (dt - ramp).dt.days


def photo_endpoint(path):
    p = pd.read_csv(path, dtype={"id1": str, "id2": str, "date": str})
    chans = ["_kelp_mean_reflectance_red", "_kelp_mean_reflectance_green",
             "_kelp_mean_reflectance_blue"]

    def side(pre, idc):
        return pd.DataFrame(dict(
            individual=p[idc].astype(str), date=p["date"],
            coverage=p[pre + "_kelp_coverage"], area=p[pre + "_kelp_area_cm2"],
            refl=p[[pre + c for c in chans]].mean(axis=1)))
    ph = pd.concat([side("left", "id1"), side("right", "id2")],
                   ignore_index=True)
    ph = ph[ph.individual.notna() & (ph.individual != "1000")].reset_index(
        drop=True)
    suf = ph.individual.str.extract(r"([A-Za-z])$")[0].fillna("")
    ph["treatment"] = suf.map(SUF_TREAT)
    ph = ph[ph.treatment.notna()].reset_index(drop=True)
    ph["cross_id"] = pd.to_numeric(
        ph.individual.str.replace(r"[^0-9]", "", regex=True),
        errors="coerce").astype("Int64")
    ph["day"] = days(ph, "date", "%Y%m%d")
    ph = ph[ph.day >= 0]
    return (ph.sort_values("day").groupby("individual").tail(1)
            [["cross_id", "treatment", "coverage", "area", "refl"]])


def fvfm_endpoint(path):
    f = pd.read_csv(path)
    f = f[(~f["low_tissue"].astype(bool)) & (f["FvFm"] > 0)].copy()
    f["day"] = days(f, "date_iso", "%Y-%m-%d")
    f = f[f.day >= 0]
    return (f.sort_values("day").groupby("individual").tail(1)
            [["cross_id", "treatment", "FvFm", "necrosis_pct"]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fvfm")
    ap.add_argument("photos")
    ap.add_argument("--out", default="figures_cross")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    m = photo_endpoint(args.photos).merge(
        fvfm_endpoint(args.fvfm), on=["cross_id", "treatment"], how="inner")
    print(f"Matched individuals: {len(m)}")
    print("matched by treatment:", m.treatment.value_counts().to_dict())

    panels = [("coverage", "photo coverage"),
              ("area", "photo area (cm2)"),
              ("refl", "photo reflectance (calibrated)"),
              ("necrosis_pct", "PAM necrosis %")]
    draw_order = ["20C", "17C", "14C"]      # 14C last so it sits on top
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, (col, lab) in zip(axes.ravel(), panels):
        for t in draw_order:
            s = m[m.treatment == t]
            ax.scatter(s["FvFm"], s[col], s=12, alpha=.5,
                       color=TREAT_COLOUR[t], edgecolor="none", label=t)
        r = m[["FvFm", col]].corr().iloc[0, 1]
        ax.set_xlabel("Fv/Fm (endpoint)")
        ax.set_ylabel(lab)
        ax.set_title(f"{lab}   r = {r:.2f}", fontsize=10)
    axes[0, 0].legend(title="treatment", fontsize=8)
    fig.suptitle("Fv/Fm vs photo metrics per individual at ~day 14")
    fig.tight_layout()
    p = os.path.join(args.out, "fvfm_vs_photos.png")
    fig.savefig(p, dpi=140)
    plt.close(fig)

    corr = {c: round(m[["FvFm", c]].corr().iloc[0, 1], 3)
            for c, _ in panels}
    print("Correlations with Fv/Fm:", corr)
    print(f"Wrote {p}")


if __name__ == "__main__":
    main()