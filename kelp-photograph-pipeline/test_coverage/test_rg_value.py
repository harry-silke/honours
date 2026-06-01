#!/usr/bin/env python3
"""
Auto ROI detection — finds petri dishes automatically using Hough circle transform.
Run this on several images to validate before integrating into the pipeline.

Usage:
    python test_auto_roi.py
"""

import numpy as np
import matplotlib.pyplot as plt
import tifffile
import cv2

# ── Configuration ─────────────────────────────────────────────────────────────
TIF_PATH  = "/Users/harrysilke/Desktop/test_images/C16_C13_270326.tif"
THRESHOLD = 1.05
MIN_GREEN_VALUE = 30

EXPECTED_CIRCLES = 4    # total circles including lids
KEEP_SMALLEST    = 2    # keep the 2 smallest (the actual kelp dishes)
MIN_RADIUS       = 550
MAX_RADIUS       = 900  # back to 900 to catch everything
MIN_DIST         = 800
PARAM2           = 25
SHRINK_FACTOR = 0.90   # shrink detected circle to avoid rim
SCALE         = 0.25    # downsample factor for Hough (speeds up detection)
# ─────────────────────────────────────────────────────────────────────────────

def detect_circles(img: np.ndarray, scale: float = SCALE) -> np.ndarray | None:
    """
    Detect circles in img. Returns array of (cx, cy, r) in full-res coordinates,
    or None if nothing found.
    """
    h, w = img.shape[:2]
    grey    = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
    blurred = cv2.GaussianBlur(grey, (31, 31), 0)  # was (15, 15)    
    small   = cv2.resize(blurred, (int(w * scale), int(h * scale)))
    circles = cv2.HoughCircles(
        small,
        cv2.HOUGH_GRADIENT,
        dp=1,
        minDist=int(MIN_DIST * scale),
        param1=50,
        param2=PARAM2,
        minRadius=int(MIN_RADIUS * scale),
        maxRadius=int(MAX_RADIUS * scale),
    )

    if circles is None:
        return None

    # Scale back to full resolution
    return (circles[0] / scale).astype(int)


def build_circle_mask(cx: int, cy: int, r: int, shape: tuple) -> np.ndarray:
    h, w = shape[:2]
    yy, xx = np.ogrid[:h, :w]
    return ((xx - cx)**2 + (yy - cy)**2) <= r**2


# Load image
img = tifffile.imread(TIF_PATH)
if img.ndim == 4:
    img = img[0]
if img.dtype != np.uint8:
    img = (img / 256).astype(np.uint8)
if img.ndim == 3 and img.shape[0] == 3:
    img = np.moveaxis(img, 0, -1)

h, w = img.shape[:2]
print(f"Image shape: {img.shape}")

circles = detect_circles(img)

# Sort by radius ascending, keep the 2 smallest (kelp dishes, not lids)
circles = circles[circles[:, 2].argsort()]
circles = circles[:KEEP_SMALLEST]

# Then sort left to right for labelling
circles = circles[circles[:, 0].argsort()]

if circles is None:
    print("ERROR: No circles detected. Try lowering PARAM2 or adjusting MIN/MAX_RADIUS.")
    exit(1)

print(f"Detected {len(circles)} circle(s)")

if len(circles) != EXPECTED_CIRCLES:
    print(f"WARNING: Expected {EXPECTED_CIRCLES}, got {len(circles)}.")
    print("Adjust PARAM2, MIN_DIST, or MIN/MAX_RADIUS.")

# Sort left to right for consistent Left/Right labelling
circles = circles[circles[:, 0].argsort()]
labels  = ["Left", "Right"] if len(circles) == 2 else [f"Dish_{i+1}" for i in range(len(circles))]

roi_masks = {}
for label, (cx, cy, r) in zip(labels, circles):
    r_shrunk = int(r * SHRINK_FACTOR)
    mask = build_circle_mask(cx, cy, r_shrunk, img.shape)
    roi_masks[label] = mask
    print(f"  {label}: centre=({cx}, {cy})  r={r}  shrunk r={r_shrunk}  pixels={mask.sum()}")

# R:G ratio
r_ch = img[:, :, 0].astype(float)
g_ch = img[:, :, 1].astype(float)
with np.errstate(divide="ignore", invalid="ignore"):
    rg = np.where(g_ch > MIN_GREEN_VALUE, r_ch / g_ch, np.nan)

combined_mask = np.zeros((h, w), dtype=bool)
for mask in roi_masks.values():
    combined_mask |= mask

rg_roi   = rg[combined_mask]
rg_valid = rg_roi[~np.isnan(rg_roi)]

print(f"\nValid ROI pixels : {len(rg_valid)}")
print(f"R:G median       : {np.median(rg_valid):.3f}")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle(TIF_PATH.split("/")[-1], fontsize=11)

# Panel 1: detected circles
axes[0].imshow(img)
axes[0].set_title("Detected circles\n(grey=full, cyan=shrunk ROI)")
for label, (cx, cy, r) in zip(labels, circles):
    r_shrunk = int(r * SHRINK_FACTOR)
    axes[0].add_patch(plt.Circle((cx, cy), r,        color="grey", fill=False, linewidth=1.5, linestyle="--"))
    axes[0].add_patch(plt.Circle((cx, cy), r_shrunk, color="cyan", fill=False, linewidth=2))
    axes[0].text(cx, cy - r_shrunk - 20, label, color="cyan", fontsize=10, ha="center",
                 bbox=dict(facecolor="black", alpha=0.4, pad=2, edgecolor="none"))
axes[0].axis("off")

# Panel 2: R:G heatmap
rg_display = np.full((h, w), np.nan)
rg_display[combined_mask] = rg[combined_mask]
im = axes[1].imshow(rg_display, cmap="RdYlGn_r", vmin=0.8, vmax=1.3)
axes[1].set_title("R:G ratio within detected ROIs")
plt.colorbar(im, ax=axes[1], fraction=0.03)
axes[1].axis("off")

# Panel 3: kelp classification
kelp_mask = combined_mask & (rg > THRESHOLD)
overlay   = img.copy()
overlay[kelp_mask] = [255, 80, 0]
n_kelp = int(kelp_mask.sum())
n_roi  = int(combined_mask.sum())
axes[2].imshow(overlay)
axes[2].set_title(f"Kelp at threshold={THRESHOLD}\n{n_kelp} px, {n_kelp/n_roi:.1%} of ROI")
axes[2].axis("off")

plt.tight_layout()
plt.savefig("auto_roi_check.png", dpi=150)
print("\nSaved: auto_roi_check.png")

# Per-ROI breakdown
print("\nPer-ROI coverage:")
for label, mask in roi_masks.items():
    k = int((mask & (rg > THRESHOLD)).sum())
    n = int(mask.sum())
    print(f"  {label}: {k/n:.1%} kelp ({k} px)")

# Threshold sweep
print("\nThreshold sweep:")
for t in [0.95, 1.0, 1.05, 1.1, 1.15, 1.2]:
    frac = float((rg_valid > t).mean())
    print(f"  {t:.2f} -> {frac:.1%}")

plt.show()