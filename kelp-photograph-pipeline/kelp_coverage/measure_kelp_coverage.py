"""
measure_kelp_coverage.py
------------------------
Measures green-reflectance-based kelp coverage and physical area (cm²)
from a folder of TIFF images using ROI masks from a single zip file in
the same folder.

Kelp classifier: R:G > 1.05, excluding pixels where G < 30 (dark/shadow).
This matches the main pipeline's kelp_coverage metric.

Usage
-----
    python measure_kelp_coverage.py <folder>

    e.g.  python measure_kelp_coverage.py ./images/17.03.2026/

Requirements
------------
    pip install tifffile numpy pandas roifile tqdm Pillow

Output
------
    Prints per-image results to stdout.
    Writes kelp_coverage_results.csv into the image folder.
"""

import sys
import zipfile
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import roifile
from tqdm import tqdm

warnings.filterwarnings("ignore")

# --- Classifier parameters (keep in sync with main pipeline) ---
RG_THRESHOLD  = 1.05   # R:G ratio threshold for kelp classification
DARK_CUTOFF   = 30     # exclude pixels where green channel < this value

# --- Spatial calibration ---
# Measured from Calibrite ColorChecker ruler: 30 pixels per mm
MM_PER_PIXEL  = 1 / 30
CM2_PER_PIXEL = (MM_PER_PIXEL ** 2) / 100   # mm² -> cm²


def find_roi_zip(folder: Path) -> Path:
    zips = list(folder.glob("*.zip"))
    if len(zips) == 0:
        raise FileNotFoundError(f"No zip file found in {folder}")
    if len(zips) > 1:
        raise ValueError(
            f"Multiple zip files found in {folder}: {[z.name for z in zips]}\n"
            "Place exactly one ROI zip in the folder."
        )
    return zips[0]


def load_roi_masks(roi_zip: Path, image_shape: tuple) -> dict:
    from PIL import Image, ImageDraw

    masks = {}
    with zipfile.ZipFile(roi_zip) as zf:
        roi_names = [n for n in zf.namelist() if n.endswith(".roi")]
        with tempfile.TemporaryDirectory() as tmp:
            zf.extractall(tmp)
            for name in roi_names:
                roi_path = Path(tmp) / name
                roi = roifile.roiread(str(roi_path))
                coords = roi.coordinates()
                if coords is None or len(coords) == 0:
                    continue
                # roifile returns (col, row) for ellipse ROIs — swap to (row, col)
                rows = coords[:, 1].astype(float)
                cols = coords[:, 0].astype(float)
                row_idx = np.clip(rows.astype(int), 0, image_shape[0] - 1)
                col_idx = np.clip(cols.astype(int), 0, image_shape[1] - 1)
                img = Image.new("L", (image_shape[1], image_shape[0]), 0)
                poly = list(zip(col_idx.tolist(), row_idx.tolist()))
                if len(poly) >= 3:
                    ImageDraw.Draw(img).polygon(poly, fill=1)
                masks[Path(name).stem] = np.array(img, dtype=bool)
    return masks


def classify_kelp(r: np.ndarray, g: np.ndarray, mask: np.ndarray):
    roi_r = r[mask].astype(float)
    roi_g = g[mask].astype(float)

    valid = roi_g >= DARK_CUTOFF
    roi_r = roi_r[valid]
    roi_g = roi_g[valid]

    if len(roi_g) == 0:
        return np.nan, np.nan

    kelp = (roi_r / roi_g) > RG_THRESHOLD
    return float(kelp.sum() / len(kelp)), float(kelp.sum())


def process_image(tif_path: Path, masks: dict) -> dict:
    rgb = tifffile.imread(str(tif_path))[0]
    r = rgb[:, :, 0]
    g = rgb[:, :, 1]
    image_shape = r.shape

    coverages = []
    areas_cm2 = []
    roi_results = {}

    for roi_name, mask in masks.items():
        if mask.shape != image_shape:
            roi_results[f"coverage_{roi_name}"] = np.nan
            roi_results[f"area_cm2_{roi_name}"] = np.nan
            continue
        cov, kelp_px = classify_kelp(r, g, mask)
        area = kelp_px * CM2_PER_PIXEL if not np.isnan(kelp_px) else np.nan
        roi_results[f"coverage_{roi_name}"] = cov
        roi_results[f"area_cm2_{roi_name}"] = area
        if not np.isnan(cov):
            coverages.append(cov)
        if not np.isnan(area):
            areas_cm2.append(area)

    row = {
        "filename":      tif_path.name,
        "mean_coverage": float(np.nanmean(coverages)) if coverages else np.nan,
        "mean_area_cm2": float(np.nanmean(areas_cm2)) if areas_cm2 else np.nan,
    }
    row.update(roi_results)
    return row


def main():
    if len(sys.argv) != 2:
        print("Usage: python measure_kelp_coverage.py <folder>")
        sys.exit(1)

    folder = Path(sys.argv[1]).resolve()
    if not folder.is_dir():
        print(f"Not a directory: {folder}")
        sys.exit(1)

    tif_files = (sorted(folder.glob("*.tif"))  + sorted(folder.glob("*.TIF"))
               + sorted(folder.glob("*.tiff")) + sorted(folder.glob("*.TIFF")))
    if not tif_files:
        print(f"No TIFF files found in {folder}")
        sys.exit(1)

    roi_zip = find_roi_zip(folder)
    print(f"Folder:      {folder}")
    print(f"ROI zip:     {roi_zip.name}")
    print(f"Images:      {len(tif_files)}")
    print(f"Calibration: {1/MM_PER_PIXEL:.0f} px/mm  ({CM2_PER_PIXEL*1e6:.4f} mm^2/px)")
    print()

    print("Loading ROI masks ...")
    img0 = tifffile.imread(str(tif_files[0]))[0]
    image_shape = img0.shape[:2]
    del img0

    masks = load_roi_masks(roi_zip, image_shape)
    print(f"ROIs loaded: {list(masks.keys())}")
    print()

    rows = []
    for tif in tqdm(tif_files, desc="Processing"):
        try:
            row = process_image(tif, masks)
        except Exception as e:
            row = {"filename": tif.name, "mean_coverage": np.nan,
                   "mean_area_cm2": np.nan, "error": str(e)}
        rows.append(row)

    df = pd.DataFrame(rows)

    print("\nPer-image results:")
    print(df[["filename", "mean_coverage", "mean_area_cm2"]].to_string(index=False))

    folder_mean_cov  = df["mean_coverage"].mean()
    folder_mean_area = df["mean_area_cm2"].mean()
    folder_std_area  = df["mean_area_cm2"].std()
    print(f"\nFolder average kelp coverage: {folder_mean_cov*100:.2f}%")
    print(f"Folder average kelp area:     {folder_mean_area:.4f} cm2  (SD {folder_std_area:.4f})")

    out_csv = folder / "kelp_coverage_results.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nResults saved to: {out_csv}")


if __name__ == "__main__":
    main()