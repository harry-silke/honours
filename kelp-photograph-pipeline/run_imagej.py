#!/usr/bin/env python3
"""
Stage 3: Direct TIF Measurement (Pure Python, no Fiji)
=======================================================
Reads TIF files directly from each date folder, applies the per-channel
Rodbard calibration, measures reflectance and kelp coverage within each
ImageJ ROI, groups photos of the same sample pair, and writes
analysis-ready CSVs.

Registration is intentionally skipped — the camera position is fixed.

Filename convention for input TIFs:  ID1_ID2_DATE_N.tif
  - N (1-4) is the photo number within the pair (missing _N is treated as N=1).
  - DATE is treated as an opaque string (no parsing/format check here).
  - Photos sharing the same ID1_ID2_DATE form a "pair group"; per-photo stats
    are averaged within the group, giving one row per group.

Outputs (one row per pair group):
  - Per-channel calibrated stats per side (left/right): mean, stddev, min, modal
  - Per-side kelp_coverage           (fraction of ROI classified as kelp)
  - Per-side per-channel kelp_mean_reflectance (calibrated, kelp-only pixels)
  - Per-side kelp_area_cm2           (kelp pixel count → cm² via CM2_PER_PIXEL)
  - n_photos and a `flag` column noting unexpected group sizes or sentinel IDs

Sentinel "1000" (unreadable lid) in either ID drops that side's data for the
whole group; the other side is still computed.

Output files:
  - <results_dir>/<label>_results.csv         per timepoint
  - <results_dir.parent>/master_results.csv   across all timepoints

A full run (no --timepoint) rebuilds master_results.csv from scratch from
that run's per-timepoint CSVs. A --timepoint run updates/replaces only that
date's rows in the existing master, preserving the rest.

Prerequisites:
  - Calibration .txt files in calibrations_dir for each time point
  - Per-timepoint ROI zip files in roi_maps_dir (specified in config.yaml)
  - TIF files in images_dir/<label>/ named ID1_ID2_DATE[_N].tif

Dependencies:
  pip install tifffile roifile scipy scikit-image pyyaml numpy pandas

Usage:
    python run_imagej.py
    python run_imagej.py --timepoint 10.03.2026
    python run_imagej.py --dry-run
"""

import argparse
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

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

try:
    import tifffile
except ImportError:
    print("ERROR: tifffile not installed. Run: pip install tifffile")
    sys.exit(1)

try:
    import roifile
except ImportError:
    print("ERROR: roifile not installed. Run: pip install roifile")
    sys.exit(1)

try:
    from scipy.optimize import curve_fit
except ImportError:
    print("ERROR: scipy not installed. Run: pip install scipy")
    sys.exit(1)

try:
    from skimage.draw import polygon
except ImportError:
    print("ERROR: scikit-image not installed. Run: pip install scikit-image")
    sys.exit(1)


# ── Constants ────────────────────────────────────────────────────────────────────

CHANNELS = ("red", "green", "blue")
SIDES    = ("left", "right")                # ROI labels in zip + column prefixes
STATS    = ("mean", "stddev", "min", "modal")

# Grouping / filename convention
EXPECTED_COUNTS = (2, 3, 4)
SENTINEL_ID     = "1000"

# Spatial calibration — ported from measure_kelp_coverage.py.
# Calibrite ColorChecker ruler: 30 pixels per mm.
MM_PER_PIXEL  = 1 / 30
CM2_PER_PIXEL = (MM_PER_PIXEL ** 2) / 100   # mm² -> cm²


def build_column_order() -> list[str]:
    """Stable column order for both per-timepoint and master CSVs."""
    cols = ["pair_key", "id1", "id2", "date", "n_photos", "flag"]
    for side in SIDES:
        for ch in CHANNELS:
            for stat in STATS:
                cols.append(f"{side}_{ch}_{stat}")
    for side in SIDES:
        cols.append(f"{side}_kelp_coverage")
    for side in SIDES:
        for ch in CHANNELS:
            cols.append(f"{side}_kelp_mean_reflectance_{ch}")
    for side in SIDES:
        cols.append(f"{side}_kelp_area_cm2")
    return cols


COLUMN_ORDER = build_column_order()


# ── Config ───────────────────────────────────────────────────────────────────────

def load_config(config_path: Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


# ── Filename parsing ─────────────────────────────────────────────────────────────

def parse_stem(stem: str):
    """Parse ID1_ID2_DATE[_N] filename stem.

    Returns dict with id1, id2, date, n, pair_key — or None if it doesn't match.
    A missing _N suffix is treated as n=1 (older single-photo convention).
    Mirrors parse_stem() in measure_kelp_coverage.py."""
    parts = stem.split("_")
    if len(parts) == 3:
        id1, id2, date = parts
        n = 1
    elif len(parts) == 4 and parts[3].isdigit():
        id1, id2, date, n_str = parts
        n = int(n_str)
    else:
        return None
    return {"id1": id1, "id2": id2, "date": date, "n": n,
            "pair_key": f"{id1}_{id2}_{date}"}


# ── TIF discovery ────────────────────────────────────────────────────────────────

def list_tifs(folder: Path, exclude_prefix: str = "") -> list[Path]:
    """All TIFs in folder (any case), with AppleDouble (._) sidecars and any
    file whose stem starts with exclude_prefix removed."""
    files = (sorted(folder.glob("*.tif"))  + sorted(folder.glob("*.TIF"))
           + sorted(folder.glob("*.tiff")) + sorted(folder.glob("*.TIFF")))
    out = []
    pfx = exclude_prefix.lower()
    for f in files:
        if f.name.startswith("._"):
            continue
        if pfx and f.stem.lower().startswith(pfx):
            continue
        out.append(f)
    return out


# ── TIF loading ──────────────────────────────────────────────────────────────────

def open_tif(path: Path, retries: int = 5, delay: float = 60.0) -> np.ndarray:
    """
    Load a TIF and return (H, W, 3) uint8 array.
    Handles multi-page, (H,W,3) and (3,H,W) layouts.
    Converts 16-bit to 8-bit to match calibration range.
    """
    for attempt in range(retries):
        try:
            img = tifffile.imread(str(path))
            if img.ndim == 4:
                img = img[0]
            if img.ndim == 3 and img.shape[2] == 3:
                if img.dtype != np.uint8:
                    img = (img / 256).astype(np.uint8)
                return img
            elif img.ndim == 3 and img.shape[0] == 3:
                img = np.moveaxis(img, 0, -1)
                if img.dtype != np.uint8:
                    img = (img / 256).astype(np.uint8)
                return img
            raise ValueError(f"Unexpected image shape {img.shape}: {path}")
        except OSError as e:
            if attempt < retries - 1:
                print(f"  [RETRY {attempt+1}/{retries}] {path.name}: {e}")
                time.sleep(delay)
            else:
                raise


# ── ROI masks ────────────────────────────────────────────────────────────────────

def load_roi_masks(roi_zip: Path, image_shape: tuple) -> dict:
    """
    Parse ROI zip and return {roi_name_lower: boolean mask (H, W)}.
    Swaps coords[:, 1] / coords[:, 0] for ellipse ROIs where roifile
    returns (col, row) instead of (row, col).
    """
    h, w = image_shape[:2]
    rois = roifile.roiread(str(roi_zip))
    if not isinstance(rois, list):
        rois = [rois]

    masks = {}
    for roi in rois:
        if roi.roitype == roifile.ROI_TYPE.RECT:
            mask = np.zeros((h, w), dtype=bool)
            top    = max(0, int(roi.top))
            bottom = min(h, int(roi.bottom))
            left   = max(0, int(roi.left))
            right  = min(w, int(roi.right))
            mask[top:bottom, left:right] = True
        else:
            coords = roi.coordinates()
            if coords is None or len(coords) < 3:
                print(f"  [WARN] ROI '{roi.name}' has no usable coordinates — skipped")
                continue
            # roifile returns (col, row) for ellipse ROIs — swap to (row, col)
            rows = coords[:, 1].astype(float)
            cols = coords[:, 0].astype(float)
            rr, cc = polygon(rows, cols, shape=(h, w))
            mask = np.zeros((h, w), dtype=bool)
            mask[rr, cc] = True

        masks[roi.name.lower()] = mask
        print(f"  ROI '{roi.name}': {mask.sum()} pixels")

    return masks


# ── Calibration ──────────────────────────────────────────────────────────────────

def find_calibration_txt(calibrations_dir: Path, date_short: str, channel: str) -> Path | None:
    expected = calibrations_dir / f"{date_short}_calibration_{channel}.txt"
    return expected if expected.exists() else None


def load_calibration_values(cal_file: Path) -> list[float]:
    lines = cal_file.read_text().strip().splitlines()
    return [float(l.split()[0]) for l in lines if l.strip()]


def rodbard(x, a, b, c, d):
    return d + (a - d) / (1.0 + (x / c) ** b)


def fit_rodbard(raw_values: list[float], reflectance_values: list[float]):
    x  = np.array(raw_values)
    y  = np.array(reflectance_values)
    p0 = [y.max(), 1.0, np.median(x), y.min()]
    popt, _ = curve_fit(rodbard, x, y, p0=p0, maxfev=10000)
    return lambda vals: rodbard(np.asarray(vals, dtype=float), *popt)


# ── Kelp mask ────────────────────────────────────────────────────────────────────

def get_kelp_mask(rgb: np.ndarray, roi_mask: np.ndarray, threshold: float) -> np.ndarray:
    """
    Return a boolean mask of pixels within roi_mask where:
      - green channel > 30 (excludes dark rim/shadow pixels)
      - R:G ratio > threshold (classifies as kelp)
    """
    r = rgb[:, :, 0].astype(float)
    g = rgb[:, :, 1].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rg = np.where(g > 30, r / g, 0.0)
    return roi_mask & (rg > threshold)


# ── Per-ROI measurement ──────────────────────────────────────────────────────────

def measure_roi(channel: np.ndarray, mask: np.ndarray, cal_fn) -> dict:
    """Calibrated reflectance statistics for one ROI on one channel."""
    if mask.sum() == 0:
        return {
            "area":   0,
            "mean":   float("nan"),
            "stddev": float("nan"),
            "min":    float("nan"),
            "modal":  float("nan"),
        }

    raw_pixels = channel[mask].astype(float)
    cal_pixels = cal_fn(raw_pixels)

    raw_int   = channel[mask].astype(np.uint16)
    counts    = np.bincount(raw_int)
    raw_mode  = int(np.argmax(counts))
    modal_cal = float(cal_fn(np.array([raw_mode]))[0])

    return {
        "area":   int(mask.sum()),
        "mean":   float(np.mean(cal_pixels)),
        "stddev": float(np.std(cal_pixels)),
        "min":    float(np.min(cal_pixels)),
        "modal":  modal_cal,
    }


# ── Per-photo measurement ────────────────────────────────────────────────────────

def measure_photo(img: np.ndarray, side_masks: dict, cal_fns: dict,
                  coverage_threshold: float) -> dict:
    """One photo → per-side per-channel calibrated stats + kelp metrics.

    Returns {side: {channel: {mean, stddev, min, modal, kelp_mean_reflectance},
                    kelp_coverage, kelp_area_px}}
    for each side present in side_masks ('left' / 'right')."""
    out = {}
    for side, roi_mask in side_masks.items():
        kelp_mask = get_kelp_mask(img, roi_mask, coverage_threshold)
        n_roi   = int(roi_mask.sum())
        n_kelp  = int(kelp_mask.sum())
        coverage = float(n_kelp / n_roi) if n_roi > 0 else float("nan")

        side_data = {"kelp_coverage": coverage, "kelp_area_px": n_kelp}
        for ch_idx, channel in enumerate(CHANNELS):
            ch_data    = img[:, :, ch_idx]
            stats      = measure_roi(ch_data, roi_mask, cal_fns[channel])
            kelp_stats = measure_roi(ch_data, kelp_mask, cal_fns[channel])
            side_data[channel] = {
                "mean":                  stats["mean"],
                "stddev":                stats["stddev"],
                "min":                   stats["min"],
                "modal":                 stats["modal"],
                "kelp_mean_reflectance": kelp_stats["mean"],
            }
        out[side] = side_data
    return out


# ── Group aggregation ────────────────────────────────────────────────────────────

def _nanmean(values) -> float:
    """np.nanmean without the empty-slice warning. Returns NaN on empty."""
    arr = np.array([v for v in values if v is not None], dtype=float)
    arr = arr[~np.isnan(arr)]
    return float(arr.mean()) if arr.size else float("nan")


def aggregate_group(pair_key: str, photos: list, sides_present: set) -> dict:
    """Average per-photo measurements within a pair_key group. Returns one
    flat dict whose keys match COLUMN_ORDER."""
    n_photos = len(photos)
    parsed0  = photos[0]["parsed"]
    id1, id2, date = parsed0["id1"], parsed0["id2"], parsed0["date"]

    drop = set()
    if id1 == SENTINEL_ID:
        drop.add("left")
    if id2 == SENTINEL_ID:
        drop.add("right")

    flags = []
    if n_photos not in EXPECTED_COUNTS:
        flags.append(f"unexpected_count:{n_photos}")
    if drop:
        flags.append("sentinel_dropped:" + ",".join(sorted(drop)))
    flag = ";".join(flags) if flags else "ok"

    row: dict = {
        "pair_key": pair_key, "id1": id1, "id2": id2, "date": date,
        "n_photos": n_photos, "flag": flag,
    }

    nan = float("nan")
    for side in SIDES:
        side_dropped = (side in drop) or (side not in sides_present)
        for ch in CHANNELS:
            for stat in STATS:
                if side_dropped:
                    row[f"{side}_{ch}_{stat}"] = nan
                else:
                    row[f"{side}_{ch}_{stat}"] = _nanmean(
                        p["measurements"][side][ch][stat] for p in photos
                    )
        if side_dropped:
            row[f"{side}_kelp_coverage"] = nan
        else:
            row[f"{side}_kelp_coverage"] = _nanmean(
                p["measurements"][side]["kelp_coverage"] for p in photos
            )
        for ch in CHANNELS:
            if side_dropped:
                row[f"{side}_kelp_mean_reflectance_{ch}"] = nan
            else:
                row[f"{side}_kelp_mean_reflectance_{ch}"] = _nanmean(
                    p["measurements"][side][ch]["kelp_mean_reflectance"]
                    for p in photos
                )
        if side_dropped:
            row[f"{side}_kelp_area_cm2"] = nan
        else:
            row[f"{side}_kelp_area_cm2"] = _nanmean(
                p["measurements"][side]["kelp_area_px"] * CM2_PER_PIXEL
                for p in photos
            )
    return row


# ── Master CSV ───────────────────────────────────────────────────────────────────

def write_master_csv(master_path: Path, run_dfs: dict, single_timepoint: bool) -> None:
    """Write the master CSV across all timepoints.

    - Full run (single_timepoint=False): rebuild from this run's per-timepoint
      DataFrames, replacing any existing master.
    - --timepoint run (single_timepoint=True): if a master exists, keep all
      rows whose `date` is NOT in this run's outputs and concatenate this
      run's rows on top. If no master exists, create one from this run.

    `run_dfs` is {label: DataFrame}. NB: empty run_dfs is a no-op so a fully
    skipped --timepoint doesn't clobber the existing master.
    """
    if not run_dfs:
        return

    run_concat = pd.concat(list(run_dfs.values()), ignore_index=True)

    if single_timepoint and master_path.exists():
        existing = pd.read_csv(master_path, dtype={"date": str})
        run_dates = set(run_concat["date"].astype(str).unique())
        keep = existing[~existing["date"].astype(str).isin(run_dates)]
        combined = pd.concat([keep, run_concat], ignore_index=True)
    else:
        combined = run_concat

    # Stable column order, then sort for predictable diffs
    cols = [c for c in COLUMN_ORDER if c in combined.columns]
    extras = [c for c in combined.columns if c not in cols]
    combined = combined[cols + extras]
    if {"date", "pair_key"}.issubset(combined.columns):
        combined = combined.sort_values(["date", "pair_key"]).reset_index(drop=True)

    master_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(master_path, index=False)


# ── Preflight ────────────────────────────────────────────────────────────────────

def preflight(tp: dict, config: dict, roi_zip: Path) -> list[str]:
    errors  = []
    paths   = config["paths"]
    label   = tp["label"]
    cal_prefix = config.get("calibration_filename_prefix", "calibration_")

    date_dir = Path(paths["images_dir"]) / label
    if not date_dir.exists():
        errors.append(f"Date folder not found: {date_dir}")
    else:
        if not list_tifs(date_dir, exclude_prefix=cal_prefix):
            errors.append(f"No TIF files found in {date_dir}")

    if not roi_zip.exists():
        errors.append(f"ROI map not found: {roi_zip}")

    cal_dir = Path(paths["calibrations_dir"])
    for channel in CHANNELS:
        if find_calibration_txt(cal_dir, tp["date_short"], channel) is None:
            errors.append(
                f"Calibration file missing: {tp['date_short']}_calibration_{channel}.txt "
                f"(expected in {cal_dir})"
            )

    return errors


# ── Per-timepoint processing ─────────────────────────────────────────────────────

def process_timepoint(tp: dict, config: dict, roi_zip: Path,
                      cal_fns: dict, coverage_threshold: float) -> pd.DataFrame | None:
    """Process one timepoint folder and return the per-group DataFrame.
    Returns None on hard failure (no usable images / no ROIs)."""
    label      = tp["label"]
    cal_prefix = config.get("calibration_filename_prefix", "calibration_")
    date_dir   = Path(config["paths"]["images_dir"]) / label

    tif_files = list_tifs(date_dir, exclude_prefix=cal_prefix)
    print(f"  {len(tif_files)} TIF file(s) found")
    if not tif_files:
        print("  [FAIL] No TIF files to process")
        return None

    # Open the first openable file to determine ROI dimensions
    first_img = None
    first_idx = 0
    for i, tif in enumerate(tif_files):
        try:
            first_img = open_tif(tif)
            first_idx = i
            break
        except Exception as e:
            print(f"  [WARN] Could not open {tif.name}: {e}")
    if first_img is None:
        print("  [FAIL] Could not open any TIF in this folder")
        return None

    roi_masks = load_roi_masks(roi_zip, first_img.shape)
    side_masks = {name: mask for name, mask in roi_masks.items() if name in SIDES}
    if not side_masks:
        print(f"  [FAIL] No 'left' or 'right' ROI in {roi_zip.name} "
              f"(found: {sorted(roi_masks)})")
        return None
    if set(side_masks) != set(SIDES):
        print(f"  [WARN] ROI map only has {sorted(side_masks)} — missing side(s) "
              "will be NaN in output")
    print(f"  Image shape: {first_img.shape[:2]}  |  coverage threshold: {coverage_threshold}")

    sides_present = set(side_masks.keys())
    groups: dict = defaultdict(list)
    unparseable, errors = [], []
    images_ok = 0

    for i, tif in enumerate(tif_files):
        parsed = parse_stem(tif.stem)
        if parsed is None:
            unparseable.append(tif.name)
            continue
        try:
            img = first_img if i == first_idx else open_tif(tif)
        except Exception as e:
            errors.append((tif.name, str(e)))
            continue

        try:
            measurements = measure_photo(img, side_masks, cal_fns, coverage_threshold)
        except Exception as e:
            errors.append((tif.name, str(e)))
            continue

        groups[parsed["pair_key"]].append({
            "parsed": parsed, "measurements": measurements, "filename": tif.name,
        })
        images_ok += 1

    if unparseable:
        print(f"  [WARN] {len(unparseable)} file(s) did not match ID1_ID2_DATE[_N]:")
        for name in unparseable:
            print(f"    {name}")
    if errors:
        print(f"  [WARN] {len(errors)} file(s) errored during processing:")
        for name, err in errors:
            print(f"    {name}: {err}")

    if not groups:
        print("  [FAIL] No usable images processed")
        return None

    rows = []
    for pair_key in sorted(groups):
        photos = sorted(groups[pair_key], key=lambda p: p["parsed"]["n"])
        rows.append(aggregate_group(pair_key, photos, sides_present))

    df = pd.DataFrame(rows, columns=COLUMN_ORDER)
    flagged = (df["flag"] != "ok").sum()
    print(f"  Processed: {images_ok} image(s) → {len(df)} group(s) "
          f"(flagged: {flagged})")
    return df


# ── Main ──────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Stage 3: measure TIF images directly (no Fiji, no registration). "
            "Photos sharing an ID1_ID2_DATE stem are grouped and averaged into "
            "one row each. Writes <label>_results.csv per timepoint and "
            "master_results.csv across all timepoints in this run."
        )
    )
    parser.add_argument("--config",    default="config.yaml")
    parser.add_argument("--timepoint", default=None,
                        help="Process only this time point (e.g. 10.03.2026). "
                             "Updates that date's rows in master_results.csv "
                             "instead of rebuilding it.")
    parser.add_argument("--dry-run",   action="store_true",
                        help="Preflight check only — do not process images or "
                             "write CSVs.")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: Config not found: {config_path}")
        sys.exit(1)

    config             = load_config(config_path)
    cal_dir            = Path(config["paths"]["calibrations_dir"])
    results_dir        = Path(config["paths"]["results_dir"])
    standards          = config["reflectance_standards"]
    meas_cfg           = config.get("measurements", {})
    coverage_threshold = meas_cfg.get("kelp_coverage_threshold", 1.05)
    cal_prefix         = config.get("calibration_filename_prefix", "calibration_")

    results_dir.mkdir(parents=True, exist_ok=True)

    time_points = config.get("time_points", [])
    if args.timepoint:
        time_points = [tp for tp in time_points if tp["label"] == args.timepoint]
        if not time_points:
            print(f"ERROR: Time point '{args.timepoint}' not found in config.")
            sys.exit(1)

    total_ok = total_skip = total_fail = 0
    run_dfs: dict = {}

    for tp in time_points:
        label      = tp["label"]
        date_short = tp["date_short"]
        print(f"── {label}")

        # Resolve ROI zip — per-timepoint takes priority over global fallback
        roi_map_name = tp.get("roi_map", config.get("roi_map"))
        if roi_map_name is None:
            print(f"  [SKIP] No roi_map defined for {label} in config.yaml")
            total_skip += 1
            print()
            continue
        roi_zip = Path(config["paths"]["roi_maps_dir"]) / roi_map_name

        errors = preflight(tp, config, roi_zip)
        if errors:
            for e in errors:
                print(f"  [SKIP] {e}")
            total_skip += 1
            print()
            continue

        if args.dry_run:
            date_dir  = Path(config["paths"]["images_dir"]) / label
            tif_files = list_tifs(date_dir, exclude_prefix=cal_prefix)
            print(f"  [DRY RUN] Preflight passed — {len(tif_files)} TIF(s) "
                  "would be processed")
            print(f"  ROI map: {roi_zip.name}")
            print()
            continue

        # Fit calibration curves
        cal_fns = {}
        for channel in CHANNELS:
            cal_file = find_calibration_txt(cal_dir, date_short, channel)
            raw_vals = load_calibration_values(cal_file)
            cal_fns[channel] = fit_rodbard(raw_vals, standards)
        print("  Calibration curves fitted")

        df = process_timepoint(tp, config, roi_zip, cal_fns, coverage_threshold)
        if df is None:
            total_fail += 1
            print()
            continue

        out_path = results_dir / f"{label}_results.csv"
        df.to_csv(out_path, index=False)
        print(f"  Wrote {out_path}")
        run_dfs[label] = df
        total_ok += 1
        print()

    # Master CSV
    if not args.dry_run and run_dfs:
        master_path = results_dir.parent / "master_results.csv"
        write_master_csv(master_path, run_dfs, single_timepoint=bool(args.timepoint))
        print(f"Master CSV: {master_path} "
              f"({'updated for ' + args.timepoint if args.timepoint else 'rebuilt'})")

    print("=" * 50)
    print("Stage 3 complete.")
    print(f"  OK     : {total_ok}")
    print(f"  Skipped: {total_skip}")
    print(f"  Failed : {total_fail}")

    if total_ok > 0:
        print(f"\nPer-timepoint CSVs: {results_dir}")
        print(f"Master CSV       : {results_dir.parent / 'master_results.csv'}")

    if total_fail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
