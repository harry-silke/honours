#!/usr/bin/env python3
"""
ROI placement verification — checks that each timepoint's ROI zip
lands correctly on its images.

Saves one QC image per timepoint to Desktop for visual inspection.

Usage:
    python test_roi_placement.py
"""

import numpy as np
import matplotlib.pyplot as plt
import tifffile
import roifile
from skimage.draw import polygon
from pathlib import Path

# ── One entry per timepoint ───────────────────────────────────────────────────
TIMEPOINTS = [
    {
        "label":    "10.03.2026",
        "image":    "/Users/harrysilke/Desktop/test_images/C66_C63_100326.tif",
        "roi_zip":  "roi_maps/ROI_Map_100326.zip",
    },
    {
        "label":    "17.03.2026",
        "image":    "/Users/harrysilke/Desktop/test_images/C89_C86_170326.tif",
        "roi_zip":  "roi_maps/ROI_Map_170326.zip",
    },
    {
        "label":    "27.03.2026",
        "image":    "/Users/harrysilke/Desktop/test_images/C16_C13_270326.tif",
        "roi_zip":  "roi_maps/ROI_Map_270326.zip",
    },
]

THRESHOLD       = 1.05
MIN_GREEN_VALUE = 30
# ─────────────────────────────────────────────────────────────────────────────


def load_image(path: str) -> np.ndarray:
    img = tifffile.imread(path)
    if img.ndim == 4:
        img = img[0]
    if img.dtype != np.uint8:
        img = (img / 256).astype(np.uint8)
    if img.ndim == 3 and img.shape[0] == 3:
        img = np.moveaxis(img, 0, -1)
    return img


def load_masks(roi_zip: str, shape: tuple) -> dict:
    h, w = shape[:2]
    rois = roifile.roiread(roi_zip)
    if not isinstance(rois, list):
        rois = [rois]

    masks = {}
    for roi in rois:
        if roi.roitype == roifile.ROI_TYPE.RECT:
            mask = np.zeros((h, w), dtype=bool)
            mask[max(0,int(roi.top)):min(h,int(roi.bottom)),
                 max(0,int(roi.left)):min(w,int(roi.right))] = True
        else:
            coords = roi.coordinates()
            if coords is None or len(coords) < 3:
                print(f"  [WARN] ROI '{roi.name}' skipped — no coordinates")
                continue
            # roifile returns (col, row) for ellipse ROIs — swap to (row, col)
            rows = coords[:, 1].astype(float)
            cols = coords[:, 0].astype(float)
            rr, cc = polygon(rows, cols, shape=(h, w))
            mask = np.zeros((h, w), dtype=bool)
            mask[rr, cc] = True

        masks[roi.name] = mask
        print(f"    ROI '{roi.name}': {mask.sum()} pixels")

    return masks


for tp in TIMEPOINTS:
    print(f"\n── {tp['label']}")

    if not Path(tp["image"]).exists():
        print(f"  [SKIP] Image not found: {tp['image']}")
        continue
    if not Path(tp["roi_zip"]).exists():
        print(f"  [SKIP] ROI zip not found: {tp['roi_zip']}")
        continue

    img   = load_image(tp["image"])
    h, w  = img.shape[:2]
    print(f"  Image shape: {img.shape}")

    masks = load_masks(tp["roi_zip"], img.shape)
    if not masks:
        print("  [SKIP] No usable ROIs")
        continue

    # R:G ratio
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rg = np.where(g > MIN_GREEN_VALUE, r / g, np.nan)

    combined_mask = np.zeros((h, w), dtype=bool)
    for mask in masks.values():
        combined_mask |= mask

    # Kelp classification overlay
    kelp_mask = combined_mask & (rg > THRESHOLD)
    overlay   = img.copy()
    overlay[kelp_mask] = [255, 80, 0]

    # Per-ROI coverage
    for name, mask in masks.items():
        n_roi  = int(mask.sum())
        n_kelp = int((mask & (rg > THRESHOLD)).sum())
        print(f"  {name}: {n_kelp/n_roi:.1%} kelp ({n_kelp}/{n_roi} px)")

    # Plot: original with ROI outline | kelp classification
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f"{tp['label']}  —  {Path(tp['image']).name}", fontsize=11)

    # Panel 1: original + ROI outlines
    axes[0].imshow(img)
    axes[0].set_title("Original + ROI outlines")
    for name, mask in masks.items():
        rows = np.where(mask.any(axis=1))[0]
        cols = np.where(mask.any(axis=0))[0]
        if len(rows) and len(cols):
            axes[0].add_patch(plt.Circle(
                ((cols[0] + cols[-1]) / 2, (rows[0] + rows[-1]) / 2),
                (cols[-1] - cols[0]) / 2,
                color="cyan", fill=False, linewidth=2
            ))
            axes[0].text(
                (cols[0] + cols[-1]) / 2,
                rows[0] - 20,
                name, color="cyan", fontsize=10, ha="center",
                bbox=dict(facecolor="black", alpha=0.4, pad=2, edgecolor="none")
            )
    axes[0].axis("off")

    # Panel 2: kelp classification
    n_kelp_total = int(kelp_mask.sum())
    n_roi_total  = int(combined_mask.sum())
    axes[1].imshow(overlay)
    axes[1].set_title(
        f"Kelp classification (orange)\n"
        f"{n_kelp_total} px, {n_kelp_total/n_roi_total:.1%} of ROI"
    )
    axes[1].axis("off")

    plt.tight_layout()
    out = Path(f"~/Desktop/roi_check_{tp['label'].replace('.', '')}.png").expanduser()
    plt.savefig(out, dpi=100, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")