#!/usr/bin/env python3
"""
Stage 3: Direct TIF Measurement (Pure Python, no Fiji)
=======================================================
Reads TIF files directly from each date folder, applies the existing
per-channel Rodbard calibration, measures reflectance and kelp coverage
within each ImageJ ROI, and writes per-channel CSVs.

Registration is intentionally skipped — the camera position is fixed.

Outputs per ROI per channel:
  - mean: calibrated reflectance across all ROI pixels
  - kelp_coverage: fraction of ROI pixels classified as kelp (R:G > threshold)
  - kelp_mean_reflectance: calibrated reflectance of kelp-only pixels

Prerequisites:
  - Calibration .txt files in calibrations_dir for each time point
  - Per-timepoint ROI zip files in roi_maps_dir (specified in config.yaml)
  - TIF files in images_dir/<label>/

Dependencies:
  pip install tifffile roifile scipy scikit-image pyyaml numpy

Usage:
    python run_imagej.py
    python run_imagej.py --timepoint 10.03.2026
    python run_imagej.py --dry-run
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Run: pip install pyyaml")
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


# ── Config ───────────────────────────────────────────────────────────────────────

def load_config(config_path: Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


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
    Parse ROI zip and return {roi_name: boolean mask (H, W)}.
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

        masks[roi.name] = mask
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


# ── CSV writer ───────────────────────────────────────────────────────────────────

FIELDNAMES = [
    "label", "area", "mean", "stddev", "min", "modal",
    "kelp_coverage", "kelp_area", "kelp_mean_reflectance"
]


def write_channel_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


# ── Preflight ────────────────────────────────────────────────────────────────────

def preflight(tp: dict, config: dict, roi_zip: Path) -> list[str]:
    errors  = []
    paths   = config["paths"]
    label   = tp["label"]

    date_dir = Path(paths["images_dir"]) / label
    if not date_dir.exists():
        errors.append(f"Date folder not found: {date_dir}")
    else:
        tif_files = sorted(date_dir.glob("*.tif")) + sorted(date_dir.glob("*.tiff"))
        if not tif_files:
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


# ── Main ──────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Stage 3: measure TIF images directly (no Fiji, no registration)"
    )
    parser.add_argument("--config",    default="config.yaml")
    parser.add_argument("--timepoint", default=None,
                        help="Process only this time point (e.g. 10.03.2026)")
    parser.add_argument("--dry-run",   action="store_true",
                        help="Preflight check only — do not process images")
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
            tif_files = sorted(date_dir.glob("*.tif")) + sorted(date_dir.glob("*.tiff"))
            print(f"  [DRY RUN] Preflight passed — {len(tif_files)} TIF(s) would be processed")
            print(f"  ROI map: {roi_zip.name}")
            print()
            continue

        # Fit calibration curves
        cal_fns = {}
        for channel in CHANNELS:
            cal_file = find_calibration_txt(cal_dir, date_short, channel)
            raw_vals = load_calibration_values(cal_file)
            cal_fns[channel] = fit_rodbard(raw_vals, standards)
        print(f"  Calibration curves fitted")

        # Collect TIF files (exclude calibration images)
        date_dir  = Path(config["paths"]["images_dir"]) / label
        tif_files = sorted(
            f for f in (
                list(date_dir.glob("*.tif")) + list(date_dir.glob("*.tiff"))
            )
            if not f.stem.lower().startswith(cal_prefix.lower())
        )
        print(f"  {len(tif_files)} TIF file(s) found")

        if not tif_files:
            print(f"  [FAIL] No TIF files to process")
            total_fail += 1
            print()
            continue

        # Load ROI masks once per timepoint
        try:
            first_img = open_tif(tif_files[0])
        except Exception as e:
            print(f"  [FAIL] Could not open first TIF to read dimensions: {e}")
            total_fail += 1
            print()
            continue

        roi_masks = load_roi_masks(roi_zip, first_img.shape)
        if not roi_masks:
            print(f"  [FAIL] No usable ROIs in {roi_zip}")
            total_fail += 1
            print()
            continue

        print(f"  Image shape: {first_img.shape[:2]}  |  coverage threshold: {coverage_threshold}")

        images_ok = images_fail = 0
        pending = [(tif_files[0], first_img)] + [(p, None) for p in tif_files[1:]]

        for tif_path, cached_img in pending:
            try:
                img = cached_img if cached_img is not None else open_tif(tif_path)
            except Exception as e:
                print(f"  [WARN] Could not open {tif_path.name}: {e}")
                images_fail += 1
                continue

            channel_rows = {ch: [] for ch in CHANNELS}

            for roi_name, roi_mask in roi_masks.items():
                kelp_mask = get_kelp_mask(img, roi_mask, coverage_threshold)
                n_roi     = int(roi_mask.sum())
                n_kelp    = int(kelp_mask.sum())
                coverage  = float(n_kelp / n_roi) if n_roi > 0 else 0.0

                for ch_idx, channel in enumerate(CHANNELS):
                    ch_data    = img[:, :, ch_idx]
                    stats      = measure_roi(ch_data, roi_mask, cal_fns[channel])
                    kelp_stats = measure_roi(ch_data, kelp_mask, cal_fns[channel])

                    channel_rows[channel].append({
                        "label":                 roi_name,
                        "area":                  stats["area"],
                        "mean":                  round(stats["mean"],   3),
                        "stddev":                round(stats["stddev"], 3),
                        "min":                   round(stats["min"],    3),
                        "modal":                 round(stats["modal"],  3),
                        "kelp_coverage":         round(coverage,        4),
                        "kelp_area":             n_kelp,
                        "kelp_mean_reflectance": round(kelp_stats["mean"], 3)
                            if not np.isnan(kelp_stats["mean"]) else "",
                    })

            for channel in CHANNELS:
                out_path = results_dir / f"{tif_path.stem}_{channel}.csv"
                write_channel_csv(out_path, channel_rows[channel])

            print(f"  [OK]   {tif_path.name}")
            images_ok += 1

        print(f"  Processed: {images_ok} OK, {images_fail} failed")
        total_ok += 1
        print()

    print("=" * 50)
    print("Stage 3 complete.")
    print(f"  OK     : {total_ok}")
    print(f"  Skipped: {total_skip}")
    print(f"  Failed : {total_fail}")

    if total_ok > 0:
        print(f"\nCSVs written to: {results_dir}")
        print("Next: python aggregate_results.py")

    if total_fail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()