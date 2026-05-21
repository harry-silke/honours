import numpy as np
import tifffile
import roifile
import zipfile
import tempfile
from PIL import Image, ImageDraw
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IMG_PATH = "/Users/harrysilke/Desktop/test_images/tiff/DSC05480.tif"
ZIP_PATH = "/Users/harrysilke/Desktop/test_images/tiff/RoiSet.zip"
OUT_PATH = "/Users/harrysilke/Desktop/test_images/qc_check.png"

img = tifffile.imread(IMG_PATH)[0]
print("Image shape:", img.shape)
r = img[:, :, 0].astype(float)
g = img[:, :, 1].astype(float)
image_shape = img.shape[:2]

overlay = img.copy()

with zipfile.ZipFile(ZIP_PATH) as zf:
    with tempfile.TemporaryDirectory() as tmp:
        zf.extractall(tmp)
        for roi_name in ["Left", "Right"]:
            roi = roifile.roiread(f"{tmp}/{roi_name}.roi")
            coords = roi.coordinates()
            rows = coords[:, 1].astype(float)
            cols = coords[:, 0].astype(float)
            row_idx = np.clip(rows.astype(int), 0, image_shape[0] - 1)
            col_idx = np.clip(cols.astype(int), 0, image_shape[1] - 1)
            mask_img = Image.new("L", (image_shape[1], image_shape[0]), 0)
            ImageDraw.Draw(mask_img).polygon(list(zip(col_idx.tolist(), row_idx.tolist())), fill=1)
            mask = np.array(mask_img, dtype=bool)

            roi_r = r[mask]
            roi_g = g[mask]
            valid = roi_g >= 30
            kelp = np.zeros(mask.sum(), dtype=bool)
            kelp[valid] = (roi_r[valid] / roi_g[valid]) > 1.05
            coverage = kelp.sum() / valid.sum()
            area_cm2 = kelp.sum() * ((1/30)**2 / 100)
            print(f"{roi_name} — Coverage: {coverage*100:.2f}%  Area: {area_cm2:.4f} cm2")

            valid_mask = np.zeros(image_shape, dtype=bool)
            valid_mask[mask] = valid
            kelp_mask = np.zeros(image_shape, dtype=bool)
            kelp_mask[valid_mask] = kelp[valid]
            overlay[kelp_mask] = [255, 140, 0]

plt.figure(figsize=(14, 9))
plt.imshow(overlay)
plt.title("Kelp classification — Left + Right ROI (orange = kelp)")
plt.axis("off")
plt.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
print(f"Saved to {OUT_PATH}")

thresholds = [1.00, 1.01, 1.02, 1.03, 1.04, 1.05]

fig, axes = plt.subplots(1, len(thresholds), figsize=(30, 8))

for ax, thresh in zip(axes, thresholds):
    overlay = img.copy()
    with zipfile.ZipFile(ZIP_PATH) as zf:
        with tempfile.TemporaryDirectory() as tmp:
            zf.extractall(tmp)
            for roi_name in ["Left", "Right"]:
                roi = roifile.roiread(f"{tmp}/{roi_name}.roi")
                coords = roi.coordinates()
                rows = coords[:, 1].astype(float)
                cols = coords[:, 0].astype(float)
                row_idx = np.clip(rows.astype(int), 0, image_shape[0] - 1)
                col_idx = np.clip(cols.astype(int), 0, image_shape[1] - 1)
                mask_img = Image.new("L", (image_shape[1], image_shape[0]), 0)
                ImageDraw.Draw(mask_img).polygon(list(zip(col_idx.tolist(), row_idx.tolist())), fill=1)
                mask = np.array(mask_img, dtype=bool)

                roi_r = r[mask]
                roi_g = g[mask]
                valid = roi_g >= 30
                kelp = np.zeros(mask.sum(), dtype=bool)
                kelp[valid] = (roi_r[valid] / roi_g[valid]) > thresh

                valid_mask = np.zeros(image_shape, dtype=bool)
                valid_mask[mask] = valid
                kelp_mask = np.zeros(image_shape, dtype=bool)
                kelp_mask[valid_mask] = kelp[valid]
                overlay[kelp_mask] = [255, 140, 0]

                coverage = kelp.sum() / valid.sum()
                area_cm2 = kelp.sum() * ((1/30)**2 / 100)
                print(f"thresh {thresh} {roi_name} — {coverage*100:.2f}%  {area_cm2:.4f} cm2")

    ax.imshow(overlay)
    ax.set_title(f"threshold {thresh}", fontsize=14)
    ax.axis("off")

plt.tight_layout()
plt.savefig("/Users/harrysilke/Desktop/test_images/threshold_comparison.png", dpi=120, bbox_inches="tight")
print("Saved to ~/Desktop/test_images/threshold_comparison.png")