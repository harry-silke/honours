#!/usr/bin/env python3
"""
Stage 4: Results Aggregator
============================
Reads all raw ImageJ CSVs from results/raw/, parses filenames to extract
sample IDs and metadata, maps left/right ROI labels to actual sample IDs,
and outputs one clean analysis-ready CSV.

Filename convention: sampleID1_sampleID2_DDMMYY_channel.csv
Example: C1_C2_100226_red.csv
  → left ROI  = C1
  → right ROI = C2
  → date      = 100226
  → channel   = red

Usage:
    python scripts/aggregate_results.py
    python scripts/aggregate_results.py --output results/final.csv
    python scripts/aggregate_results.py --wide
"""

import argparse
import re
import sys
import warnings
from pathlib import Path
from datetime import datetime

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Run: pip install pyyaml")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas not installed. Run: pip install pandas")
    sys.exit(1)


# ── Filename parser ────────────────────────────────────────────────────────────

FILENAME_PATTERN = re.compile(
    r"^(?P<sample_id1>[A-Za-z0-9]+)"
    r"_(?P<sample_id2>[A-Za-z0-9]+)"
    r"_(?P<date>\d{6})"
    r"_(?P<channel>red|green|blue)"
    r"\.csv$",
    re.IGNORECASE,
)


def parse_filename(filename: str) -> dict | None:
    match = FILENAME_PATTERN.match(filename)
    if not match:
        return None
    return {
        "sample_id1": match.group("sample_id1").upper(),
        "sample_id2": match.group("sample_id2").upper(),
        "date": match.group("date"),
        "channel": match.group("channel").lower(),
    }


def parse_date(date_str: str) -> str:
    """Convert DDMMYY to DD.MM.20YY for readability."""
    try:
        return datetime.strptime(date_str, "%d%m%y").strftime("%d.%m.%Y")
    except ValueError:
        return date_str


# ── CSV reader ─────────────────────────────────────────────────────────────────

def read_channel_csv(csv_path: Path, metadata: dict) -> pd.DataFrame:
    """
    Read one ImageJ results CSV, map left/right ROI labels to sample IDs,
    and attach metadata columns.
    """
    df = pd.read_csv(csv_path)

    # Normalise column names
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    # Rename ImageJ column variants to standard names
    rename_map = {
        "mean": "mean_reflectance",
        "mean_gray_value": "mean_reflectance",
        "stddev": "std_dev",
        "min": "min_gray",
        "max": "max_gray",
        "min_gray_value": "min_gray",
        "max_gray_value": "max_gray",
        "mode": "modal_gray",
        "modal_gray_value": "modal_gray",
        "label": "roi_label",
    }
    df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns}, inplace=True)

    # Map left/right ROI labels to actual sample IDs
    roi_to_sample = {
        "left":  metadata["sample_id1"],
        "right": metadata["sample_id2"],
    }

    if "roi_label" in df.columns:
        df["sample_id"] = df["roi_label"].str.lower().str.strip().map(roi_to_sample)
        unmapped = df["sample_id"].isna().sum()
        if unmapped > 0:
            warnings.warn(
                f"  {csv_path.name}: {unmapped} ROI(s) could not be mapped to sample IDs. "
                f"Expected 'left' and 'right' — check your ROI map labels."
            )

    # Attach metadata
    df["sample_id1"]   = metadata["sample_id1"]
    df["sample_id2"]   = metadata["sample_id2"]
    df["date"]         = metadata["date"]
    df["date_readable"] = parse_date(metadata["date"])
    df["channel"]      = metadata["channel"]
    df["source_file"]  = csv_path.name

    return df


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate ImageJ results CSVs into one analysis-ready dataframe"
    )
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--output", default=None,
                        help="Output CSV path (default: results/aggregated_TIMESTAMP.csv)")
    parser.add_argument("--wide", action="store_true",
                        help="Pivot to wide format (one row per sample per timepoint, columns per channel)")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: Config not found: {config_path}")
        sys.exit(1)

    config      = yaml.safe_load(config_path.read_text())
    results_dir = Path(config["paths"]["results_dir"])

    if not results_dir.exists():
        print(f"ERROR: Results directory not found: {results_dir}")
        sys.exit(1)

    csv_files = sorted(results_dir.glob("*.csv"))
    if not csv_files:
        print(f"No CSV files found in {results_dir}")
        sys.exit(0)

    print(f"Found {len(csv_files)} CSV files\n")

    all_dfs, skipped = [], []

    for csv_path in csv_files:
        metadata = parse_filename(csv_path.name)
        if metadata is None:
            warnings.warn(
                f"Skipping '{csv_path.name}' — expected format: "
                f"sampleID1_sampleID2_DDMMYY_channel.csv"
            )
            skipped.append(csv_path.name)
            continue

        try:
            df = read_channel_csv(csv_path, metadata)
            all_dfs.append(df)
            print(f"  [OK]  {csv_path.name}  ({len(df)} ROIs)")
        except Exception as e:
            warnings.warn(f"Error reading '{csv_path.name}': {e}")
            skipped.append(csv_path.name)

    if not all_dfs:
        print("\nNo valid CSVs to aggregate.")
        sys.exit(1)

    combined = pd.concat(all_dfs, ignore_index=True)

    # Column order: identifiers first, then measurements
    id_cols   = ["sample_id", "date", "date_readable", "channel", "roi_label", "sample_id1", "sample_id2"]
    meas_cols = ["mean_reflectance", "std_dev", "min_gray", "max_gray", "modal_gray", "area", "kelp_coverage", "kelp_area", "kelp_mean_reflectance"]    
    other     = [c for c in combined.columns if c not in id_cols + meas_cols + ["source_file"]]
    col_order = (
        [c for c in id_cols if c in combined.columns]
        + [c for c in meas_cols if c in combined.columns]
        + other
        + ["source_file"]
    )
    combined = combined[col_order]

    # Sort by date, sample, channel
    sort_cols = [c for c in ["date", "sample_id", "channel"] if c in combined.columns]
    combined.sort_values(sort_cols, inplace=True)

    # Optional wide pivot
    if args.wide:
        print("\nPivoting to wide format...")
        combined = combined.pivot_table(
            index=["sample_id", "date", "date_readable"],
            columns="channel",
            values=[c for c in meas_cols if c in combined.columns],
            aggfunc="first",
        )
        combined.columns = [f"{stat}_{chan}" for stat, chan in combined.columns]
        combined.reset_index(inplace=True)

    # Output path
    if args.output:
        output_path = Path(args.output)
    else:
        timestamp   = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = results_dir.parent / f"aggregated_{timestamp}.csv"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, index=False)

    # Summary
    print(f"\n{'='*50}")
    print(f"Aggregation complete")
    print(f"  Rows         : {len(combined)}")
    print(f"  Unique dates : {combined['date'].nunique() if 'date' in combined.columns else 'N/A'}")
    print(f"  Unique samples: {combined['sample_id'].nunique() if 'sample_id' in combined.columns else 'N/A'}")
    print(f"  Saved to     : {output_path}")

    if skipped:
        print(f"\n  Skipped {len(skipped)} file(s):")
        for f in skipped:
            print(f"    {f}")

    print(f"\nLoad in R:      df <- read.csv('{output_path}')")
    print(f"Load in Python: df = pd.read_csv('{output_path}')")


if __name__ == "__main__":
    main()
