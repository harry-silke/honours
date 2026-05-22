"""
measure_kelp_coverage.py
------------------------
Measures green-reflectance-based kelp coverage and physical area (cm²)
from a folder of TIFF images using ROI masks from a single zip file in
the same folder.

Photos are grouped by their ID1_ID2_DATE stem (the renaming pipeline names
them ID1_ID2_DDMMYY_N.tif with N in 2..4) and coverage is averaged across
the 2-4 photos of each pair. Files without a _N suffix are treated as a
single-photo group. Groups with an unexpected count (not 2-4) are still
processed and marked in the flag column. A sample ID of "1000" means the
lid was unreadable upstream: that side's ROI is dropped for the whole
group, but the other side is still processed.

Output is one row per pair+date group.

Usage
-----
    python measure_kelp_coverage.py <folder>

    e.g.  python measure_kelp_coverage.py ./images/17.03.2026/

Requirements
------------
    pip install tifffile numpy pandas roifile tqdm Pillow

Output
------
    Prints per-group results to stdout.
    Writes kelp_coverage_results.csv into the image folder.
"""

import sys
import zipfile
import tempfile
import warnings
from collections import defaultdict
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

# --- Grouping / filename convention ---
EXPECTED_COUNTS   = (2, 3, 4)   # photos per pair+date; anything else is flagged
SENTINEL_ID       = "1000"      # unreadable lid; drop that side
ROI_FOR_LEFT_ID   = "left"      # ROI name for sampleID1 (first id in filename)
ROI_FOR_RIGHT_ID  = "right"     # ROI name for sampleID2 (second id in filename)


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


def parse_stem(stem: str):
    """Parse ID1_ID2_DATE[_N] filename stem.

    Returns dict with id1, id2, date, n, pair_key — or None if it doesn't match.
    A missing _N (old single-photo format) is treated as n=1."""
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


def process_image(tif_path: Path, masks: dict) -> dict:
    """Per-photo result: {roi_name: (coverage, area_cm2)}. NaN if the ROI mask
    shape doesn't match or the ROI is fully dark."""
    rgb = tifffile.imread(str(tif_path))[0]
    r = rgb[:, :, 0]
    g = rgb[:, :, 1]
    image_shape = r.shape

    out = {}
    for roi_name, mask in masks.items():
        if mask.shape != image_shape:
            out[roi_name] = (np.nan, np.nan)
            continue
        cov, kelp_px = classify_kelp(r, g, mask)
        area = kelp_px * CM2_PER_PIXEL if not np.isnan(kelp_px) else np.nan
        out[roi_name] = (cov, area)
    return out


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
    roi_names = list(masks.keys())
    print(f"ROIs loaded: {roi_names}")
    print()

    # Stage 1: process every photo, bucket results by pair_key
    groups = defaultdict(list)   # pair_key -> list of dicts {n, filename, rois}
    unparseable = []
    errors = []
    for tif in tqdm(tif_files, desc="Processing"):
        parsed = parse_stem(tif.stem)
        if parsed is None:
            unparseable.append(tif.name)
            continue
        try:
            rois = process_image(tif, masks)
        except Exception as e:
            errors.append((tif.name, str(e)))
            rois = {name: (np.nan, np.nan) for name in roi_names}
        groups[parsed["pair_key"]].append({
            "n": parsed["n"], "filename": tif.name, "rois": rois,
            "id1": parsed["id1"], "id2": parsed["id2"], "date": parsed["date"],
        })

    if unparseable:
        print(f"\nWARNING: {len(unparseable)} file(s) did not match "
              f"ID1_ID2_DATE[_N]:")
        for name in unparseable:
            print(f"  {name}")
    if errors:
        print(f"\nWARNING: {len(errors)} file(s) errored during processing:")
        for name, err in errors:
            print(f"  {name}: {err}")

    # Stage 2: average across photos within each group
    rows = []
    for pair_key in sorted(groups.keys()):
        photos = sorted(groups[pair_key], key=lambda d: d["n"])
        n_photos = len(photos)
        id1, id2, date = photos[0]["id1"], photos[0]["id2"], photos[0]["date"]

        drop = set()
        if id1 == SENTINEL_ID:
            drop.add(ROI_FOR_LEFT_ID)
        if id2 == SENTINEL_ID:
            drop.add(ROI_FOR_RIGHT_ID)

        flags = []
        if n_photos not in EXPECTED_COUNTS:
            flags.append(f"unexpected_count:{n_photos}")
        if drop:
            flags.append("sentinel_dropped:" + ",".join(sorted(drop)))
        flag = ";".join(flags) if flags else "ok"

        row = {
            "pair_key": pair_key,
            "id1": id1, "id2": id2, "date": date,
            "n_photos": n_photos,
            "flag": flag,
        }

        roi_cov, roi_area = {}, {}
        for name in roi_names:
            if name in drop:
                roi_cov[name]  = np.nan
                roi_area[name] = np.nan
            else:
                covs  = [p["rois"][name][0] for p in photos
                         if not np.isnan(p["rois"][name][0])]
                areas = [p["rois"][name][1] for p in photos
                         if not np.isnan(p["rois"][name][1])]
                roi_cov[name]  = float(np.mean(covs))  if covs  else np.nan
                roi_area[name] = float(np.mean(areas)) if areas else np.nan
            row[f"coverage_{name}"] = roi_cov[name]
            row[f"area_cm2_{name}"] = roi_area[name]

        # Mean across the (non-dropped) ROIs — i.e. across the two samples
        # in the pair. NaNs are skipped.
        cov_vals  = [v for v in roi_cov.values()  if not np.isnan(v)]
        area_vals = [v for v in roi_area.values() if not np.isnan(v)]
        row["mean_coverage"] = float(np.mean(cov_vals))  if cov_vals  else np.nan
        row["mean_area_cm2"] = float(np.mean(area_vals)) if area_vals else np.nan

        rows.append(row)

    df = pd.DataFrame(rows)
    # Tidy column order
    front = ["pair_key", "id1", "id2", "date", "n_photos", "flag",
             "mean_coverage", "mean_area_cm2"]
    rest = [c for c in df.columns if c not in front]
    df = df[front + rest]

    print("\nPer-group results:")
    print(df[front].to_string(index=False))

    folder_mean_cov  = df["mean_coverage"].mean()
    folder_mean_area = df["mean_area_cm2"].mean()
    folder_std_area  = df["mean_area_cm2"].std()
    print(f"\nGroups: {len(df)}  "
          f"(flagged: {(df['flag'] != 'ok').sum()})")
    print(f"Folder average kelp coverage: {folder_mean_cov*100:.2f}%")
    print(f"Folder average kelp area:     {folder_mean_area:.4f} cm2  "
          f"(SD {folder_std_area:.4f})")

    out_csv = folder / "kelp_coverage_results.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nResults saved to: {out_csv}")


if __name__ == "__main__":
    main()