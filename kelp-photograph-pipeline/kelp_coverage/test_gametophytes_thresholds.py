from pathlib import Path

import numpy as np
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IMG_DIR = Path("/Users/harrysilke/Desktop/gametophyte_test")
IMG_PATHS = [
    IMG_DIR / "test_1.tiff",
    IMG_DIR / "test_2.tiff",
    IMG_DIR / "test_3.tiff",
]
OUT_DIR = IMG_DIR

N_THRESH = 8          # grid below is 2 x 4, so keep this at 8
RATIO_MIN = 1.00      # first threshold in the sweep
RATIO_MAX = 1.07      # last threshold in the sweep (0.01 steps at 8 values)
QC_THRESH = 1.05      # single threshold used for the qc image
MIN_G = 30            # pixels with green below this are ignored (too dark to trust the ratio)
PX_PER_MM = 30        # 30 was for the kelp photos; set from the microscope calibration for slides

px_area_cm2 = (1 / PX_PER_MM) ** 2 / 100
thresholds = np.linspace(RATIO_MIN, RATIO_MAX, N_THRESH)


def process_image(img_path):
    stem = img_path.stem
    print(f"\n===== {img_path.name} =====")

    img = tifffile.imread(img_path)
    if img.ndim == 4:     # multi page / stack: take the first frame
        img = img[0]
    print("Image shape:", img.shape)

    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)

    valid = g >= MIN_G
    ratio = np.zeros(g.shape)
    ratio[valid] = r[valid] / g[valid]

    print(f"Valid pixels: {valid.sum()} of {valid.size} ({valid.mean()*100:.1f}%)")
    print("R/G ratio percentiles over valid pixels:")
    for p in [1, 5, 25, 50, 75, 95, 99]:
        print(f"  p{p:<3} {np.percentile(ratio[valid], p):.3f}")

    def classify(thresh):
        return valid & (ratio > thresh)

    def report(label, mask):
        coverage = mask.sum() / valid.sum()
        print(f"{label}: coverage {coverage*100:.2f}%  area {mask.sum()*px_area_cm2:.4f} cm2")

    def overlay_for(mask):
        ov = img.copy()
        ov[mask] = [255, 140, 0]
        return ov

    # Single threshold QC image
    qc = classify(QC_THRESH)
    report(f"thresh {QC_THRESH:.2f}", qc)
    plt.figure(figsize=(14, 9))
    plt.imshow(overlay_for(qc))
    plt.title(f"{stem}: threshold {QC_THRESH:.2f} (orange = gametophyte)")
    plt.axis("off")
    plt.savefig(OUT_DIR / f"{stem}_qc_check.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Eight threshold sweep
    fig, axes = plt.subplots(2, 4, figsize=(28, 12))
    for ax, thresh in zip(axes.ravel(), thresholds):
        m = classify(thresh)
        report(f"thresh {thresh:.3f}", m)
        ax.imshow(overlay_for(m))
        ax.set_title(f"threshold {thresh:.3f}", fontsize=14)
        ax.axis("off")
    fig.suptitle(stem, fontsize=18)
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"{stem}_threshold_comparison.png", dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Saved {stem}_qc_check.png and {stem}_threshold_comparison.png to {OUT_DIR}")


for path in IMG_PATHS:
    process_image(path)