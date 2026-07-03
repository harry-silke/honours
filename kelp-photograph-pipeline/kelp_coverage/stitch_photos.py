"""
stitch_photos.py
----------------
Combines every per-folder kelp_photo_results.csv (one row per photo per side,
written by measure_kelp_coverage.py) into one self-contained master_photos.csv,
joining treatment, female population, and parent cultures from crosses.csv.

crosses.csv comes from build_metadata_from_master.py (cross_id, parent1_culture,
parent2_culture, parent2_pop = female population).

Treatment from the id suffix: none = 14C (A), B = 17C (B), C = 20C (C). Rows
whose id is not one of these (mislabels, CALIBRATION) are dropped.

Usage
-----
    python stitch_photos.py <root> --crosses crosses.csv --out master_photos.csv
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd

TREAT = {"": "14C", "B": "17C", "C": "20C"}
ABC = {"": "A", "B": "B", "C": "C"}

FINAL_COLS = ["photo_filename", "photo_n", "pair_key", "date", "side",
              "roi_name", "individual", "cross_id", "treatment", "treat_abc",
              "female_pop", "coverage", "area_cm2", "refl_R", "refl_G",
              "refl_B", "refl_bright", "refl_rg", "flag", "sentinel",
              "parent1_culture", "parent2_culture"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--crosses", default="crosses.csv")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    files = sorted(glob.glob(
        os.path.join(args.root, "**", "kelp_photo_results.csv"), recursive=True))
    if not files:
        sys.exit(f"No kelp_photo_results.csv found under {args.root}")
    print(f"Found {len(files)} per-photo file(s).")

    df = pd.concat([pd.read_csv(f, dtype={"individual": str, "date": str})
                    for f in files], ignore_index=True)

    suffix = df["individual"].str.extract(r"([A-Za-z])$")[0].fillna("")
    df["treatment"] = suffix.map(TREAT)
    df["treat_abc"] = suffix.map(ABC)
    dropped = df["treatment"].isna()
    if dropped.any():
        bad = sorted(df.loc[dropped, "individual"].dropna().unique())
        print(f"Dropped {int(dropped.sum())} rows with non 14/17/20C ids: {bad}")
    df = df[~dropped].copy()
    df["cross_id"] = pd.to_numeric(
        df["individual"].str.replace(r"[^0-9]", "", regex=True),
        errors="coerce").astype("Int64")

    if not os.path.isfile(args.crosses):
        sys.exit(f"{args.crosses} not found; make it with "
                 "build_metadata_from_master.py")
    cr = pd.read_csv(args.crosses)
    cr = cr.rename(columns={"parent2_pop": "female_pop"})[
        ["cross_id", "female_pop", "parent1_culture", "parent2_culture"]]
    cr["cross_id"] = cr["cross_id"].astype("Int64")
    df = df.merge(cr, on="cross_id", how="left")

    n_unmatched = int(df["female_pop"].isna().sum())
    if n_unmatched:
        print(f"NOTE: {n_unmatched} rows had a cross_id not in {args.crosses}")

    df = df[FINAL_COLS]
    out = args.out or os.path.join(args.root, "master_photos.csv")
    df.to_csv(out, index=False)
    print(f"\nWrote {len(df)} photo-rows -> {out}")
    print(f"Individuals: {df['individual'].nunique()}  "
          f"Dates: {df['date'].nunique()}  "
          f"Populations: {df['female_pop'].dropna().nunique()}")


if __name__ == "__main__":
    main()