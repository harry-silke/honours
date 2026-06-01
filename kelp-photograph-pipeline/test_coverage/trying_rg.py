import tifffile
import numpy as np
import matplotlib.pyplot as plt
import roifile
from skimage.draw import polygon
from pathlib import Path

img_path = "/Users/harrysilke/Library/CloudStorage/OneDrive-AustralianNationalUniversity/Honours/Kelp_growth_photos/17.03.2026/registered_uncompressed/C6_C3_170326.tif"
img = tifffile.imread(img_path)

roi_path = "/Users/harrysilke/Desktop/UNI/Honours/imas_hons/kelp-photograph-pipeline/roi_maps/ROI_Map.zip"

img = tifffile.imread(img_path)
if img.ndim == 3 and img.shape[2] == 3:
    red, green, blue = img[:,:,0], img[:,:,1], img[:,:,2]
else:
    red, green, blue = img[0], img[1], img[2]

rois = roifile.roiread(roi_path)
if not isinstance(rois, list):
    rois = [rois]

fig, axes = plt.subplots(1, len(rois), figsize=(12, 4))
if len(rois) == 1:
    axes = [axes]

for ax, roi in zip(axes, rois):
    coords = roi.coordinates()
    rows, cols = coords[:,0], coords[:,1]
    rr, cc = polygon(rows, cols, shape=red.shape)
    
    r_pixels = red[rr, cc].astype(float)
    g_pixels = green[rr, cc].astype(float)
    rg = np.where(g_pixels > 0, r_pixels / g_pixels, 0)
    
    ax.hist(rg, bins=100, range=(0.7, 1.5))
    ax.set_title(f"ROI: {roi.name}")
    ax.set_xlabel("R:G ratio")
    ax.set_ylabel("Pixel count")
    ax.axvline(x=1.05, color='red', linestyle='--', label='threshold=1.05')
    ax.legend()

plt.tight_layout()
plt.savefig("rg_roi_histogram.png", dpi=150)
print("Saved rg_roi_histogram.png")