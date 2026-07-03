"""
stitch_results.py
-----------------
Combines every per-folder kelp_coverage_results.csv (written by
measure_kelp_coverage.py, one per date folder) into a single
master_results.csv, renaming columns to the schema analyse_kelp.R reads.

measure_kelp_coverage.py column  ->  master (R) column
    coverage_<side>              ->  <side>_kelp_coverage
    area_cm2_<side>              ->  <side>_kelp_area_cm2
    refl_R/G/B_<side>            ->  <side>_kelp_mean_reflectance_red/green/blue
Sides are lower-cased (Left -> left) to match the R script.

Usage
-----
    python stitch_results.py /Volumes/Expansion/Heat_stress_experiment_photos
    python stitch_results.py <root> --out <root>/master_results.csv
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd

CHAN = {"R": "red", "G": "green", "B": "blue"}


def rename_map(cols):
    out = {}
    sides = set()
    for c in cols:
        m = re.match(r"coverage_(.+)$", c)
        if m:
            sides.add(m.group(1))
    for s in sides:
        lo = s.lower()
        if lo not in ("left", "right"):
            print(f"NOTE: side {s!r} is not left/right; kept as {lo!r}")
        out[f"coverage_{s}"] = f"{lo}_kelp_coverage"
        out[f"area_cm2_{s}"] = f"{lo}_kelp_area_cm2"
        for k, name in CHAN.items():
            out[f"refl_{k}_{s}"] = f"{lo}_kelp_mean_reflectance_{name}"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    files = sorted(glob.glob(
        os.path.join(args.root, "**", "kelp_coverage_results.csv"),
        recursive=True))
    if not files:
        sys.exit(f"No kelp_coverage_results.csv found under {args.root}")

    print(f"Found {len(files)} result file(s):")
    for f in files:
        print(f"  {os.path.relpath(f, args.root)}")

    frames = [pd.read_csv(f, dtype={"id1": str, "id2": str, "date": str})
              for f in files]
    combined = pd.concat(frames, ignore_index=True)
    combined = combined.rename(columns=rename_map(combined.columns))

    out = args.out or os.path.join(args.root, "master_results.csv")
    combined.to_csv(out, index=False)
    print(f"\nStitched {len(combined)} rows -> {out}")

    dates = sorted(combined["date"].dropna().unique())
    print(f"Dates: {', '.join(map(str, dates))}")
    if dates and not re.match(r"^\d{1,2}\.\d{1,2}\.\d{4}$", str(dates[0])):
        print("WARNING: dates are not DD.MM.YYYY. analyse_kelp.R uses ymd(); "
              "if these are YYYYMMDD, change ymd() to dmy() in the R script.")


if __name__ == "__main__":
    main()