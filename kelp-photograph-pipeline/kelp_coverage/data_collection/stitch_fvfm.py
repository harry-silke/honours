"""
stitch_fvfm.py
--------------
Combines every per-date fvfm_results.csv (from measure_fvfm.py) into one
self-contained master_fvfm.csv, joining treatment, female population, and
parent cultures from crosses.csv.

PAM dates are DDMMYY (e.g. 080526 = 2026-05-08); a date_iso column is added
so downstream tools can sort and rebase unambiguously. Treatment from the id
suffix: none = 14C (A), B = 17C (B), C = 20C (C). Ids outside those (mislabels)
are dropped. Rows with tissue_px < MIN_TISSUE are kept but flagged low_tissue,
since Fv/Fm off very few pixels is unreliable.

Usage
-----
    python stitch_fvfm.py <root> --crosses crosses.csv --out master_fvfm.csv
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd

TREAT = {"": "14C", "A": "14C", "B": "17C", "C": "20C"}
ABC = {"": "A", "A": "A", "B": "B", "C": "C"}
MIN_TISSUE = 200

FINAL = ["photo_filename", "pair_key", "date", "date_iso", "side",
         "individual", "cross_id", "treatment", "treat_abc", "female_pop",
         "roi_px", "tissue_px", "low_tissue", "Fo", "Fm", "FvFm",
         "necrosis_pct", "px_below_0.10", "px_below_0.15", "px_below_0.20",
         "fo_floor", "necrosis_cut", "parent1_culture", "parent2_culture"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--crosses", default="crosses.csv")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    files = sorted(glob.glob(
        os.path.join(args.root, "**", "fvfm_results.csv"), recursive=True))
    files = [f for f in files if "master" not in os.path.basename(f).lower()]
    if not files:
        sys.exit(f"No fvfm_results.csv found under {args.root}")
    print(f"Found {len(files)} per-date file(s).")

    frames = []
    for f in files:
        d = pd.read_csv(f, dtype={"individual": str, "date": str})
        if "roi_px" not in d.columns and "rect_px" in d.columns:
            d = d.rename(columns={"rect_px": "roi_px"})  # old rectangle output
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)

    # DDMMYY -> real date (restore any dropped leading zero first)
    df["date_iso"] = pd.to_datetime(
        df["date"].str.strip().str.zfill(6), format="%d%m%y",
        errors="coerce").dt.strftime("%Y-%m-%d")
    if df["date_iso"].isna().any():
        bad = sorted(df.loc[df["date_iso"].isna(), "date"].unique())
        print(f"NOTE: {int(df['date_iso'].isna().sum())} rows had unparseable "
              f"DDMMYY dates: {bad}")

    suffix = df["individual"].str.extract(r"([A-Za-z])$")[0].fillna("")
    df["treatment"] = suffix.map(TREAT)
    df["treat_abc"] = suffix.map(ABC)
    dropped = df["treatment"].isna()
    if dropped.any():
        print(f"Dropped {int(dropped.sum())} rows with non 14/17/20C ids: "
              f"{sorted(df.loc[dropped, 'individual'].dropna().unique())}")
    df = df[~dropped].copy()

    df["cross_id"] = pd.to_numeric(
        df["individual"].str.replace(r"[^0-9]", "", regex=True),
        errors="coerce").astype("Int64")
    df["low_tissue"] = df["tissue_px"].fillna(0) < MIN_TISSUE

    if not os.path.isfile(args.crosses):
        sys.exit(f"{args.crosses} not found; make it with "
                 "build_metadata_from_master.py")
    cr = pd.read_csv(args.crosses).rename(columns={"parent2_pop": "female_pop"})
    cr = cr[["cross_id", "female_pop", "parent1_culture", "parent2_culture"]]
    cr["cross_id"] = cr["cross_id"].astype("Int64")
    df = df.merge(cr, on="cross_id", how="left")

    for c in FINAL:
        if c not in df.columns:
            df[c] = pd.NA
    df = df[FINAL]

    out = args.out or os.path.join(args.root, "master_fvfm.csv")
    df.to_csv(out, index=False)
    print(f"\nWrote {len(df)} rows -> {out}")
    print(f"Dates: {df['date_iso'].nunique()}  "
          f"Individuals: {df['individual'].nunique()}  "
          f"Populations: {df['female_pop'].dropna().nunique()}  "
          f"low_tissue flagged: {int(df['low_tissue'].sum())}")


if __name__ == "__main__":
    main()